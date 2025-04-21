library(cfda)

# Parametres du processus de Markov a 3 etats

K = 3
labels = c("A", "B", "C") # Noms des etats


# Matrice de transition (probabilites egales vers les autres etats)
P = (1 - diag(K)) / (K - 1)    # off-diagonales = 0.5, diagonales = 0
# Taux de sortie (lambda) pour chaque etat
lambda_out = c(0.1, 0.1, 0.1)  # groupe hors cluster: changements moderees
                              #la duree moyenne passee dans chacun des 
                              #trois etats est la meme 

lambda_clu = c(0.01, 0.2, 0.2) # groupe cluster: etat "A" tres persistant, autres moderes

# Distribution initiale des etats
pi0_out = c(1/3, 1/3, 1/3)  # initial uniformement aleatoire
                            #chaque individu demarre dans l’un des 3 etats 
                            #avec probabilite egale

pi0_clu = c(1, 0, 0)        # initial toujours dans l'etat "A" (cluster)

# Simulation du groupe hors cluster (70 individus)
set.seed(123)  # pour reproductibilite
data_out = generate_Markov(n = 70, K = K, P = P, lambda = lambda_out, 
                         pi0 = pi0_out, Tmax = 50, labels = labels)
# Simulation du groupe cluster (30 individus)
data_clu = generate_Markov(n = 30, K = K, P = P, lambda = lambda_clu, 
                         pi0 = pi0_clu, Tmax = 50, labels = labels)

# Reindexer les id du groupe cluster (pour eviter les doublons)
data_clu$id = data_clu$id + 70

# Combiner les deux groupes 
data_all = rbind(data_out, data_clu)

summary_cfd(data_all)

head(data_all,20)
#visualiser un echantillon des trajectoires
#plotData(data_all[data_all$id <= 80, ])

# couper a Tmax = 50
data_cut = cut_data(data_all, Tmax = 50) # Toutes les trajectoires se terminent a Tmax=50
# pour que l’estimation de l’encodage et le calcul des scores 
#soient comparables d’un individu a l’autre.

summary_cfd(data_cut)




# Creer une base de fonctions B-splines sur [0,50]
Tmax = 50
m = 10  # nombre de fonctions de base B-splines
basis = create.bspline.basis(rangeval = c(0, Tmax), nbasis = m, norder = 4)

# Calculer l'encodage optimal (FMCA) sur les donnees categorielles

fmca_res = compute_optimal_encoding(data_cut, basis,nCores = 4)# nCores = 4 indique le nombre de cœurs CPU
summary(fmca_res)  # apercu des variances expliquees, etc.

par(mfrow=c(1,2))
#graphique valeurs propres cumulees normalisees
#Pour savoir combien de composantes sont necessaires pour capter l’essentiel 
#de la variabilite des trajectoires
plotEigenvalues(fmca_res, cumulative = TRUE, normalize = TRUE)

#graphique de la projection des trajectoires (ou individus) dans le plan factoriel .
plotComponent(fmca_res, comp = c(1, 2), addNames = FALSE)

# graphique de la fonction d'encodage optimale sur la PC1.
plot(fmca_res)
plot(fmca_res, addCI = TRUE)

# Calculer les scores des individus sur les composantes principales
pc_scores = predict(fmca_res)  
pc_scores
# Extraire le score de la premiere composante pour chaque individu
score_CP1 = pc_scores[, 1]

# Generer des coordonnees spatiales aleatoires pour chaque individu
n_total = 100
coords = matrix(NA, nrow = n_total, ncol = 2)
colnames(coords) = c("x", "y")

# Positions pour les 70 individus hors cluster 
#repartis uniformement dans le carre [0,0.6]x[0.0.6]

coords[1:70, "x"] = runif(70, min = 0.0, max = 0.6)
coords[1:70, "y"] = runif(70, min = 0.0, max = 0.6)

# Positions pour les 30 individus du cluster 
#repartis uniformement dans le carre [0.8,1.0]x[0.8,1.0]
coords[71:100, "x"] = runif(30, min = 0.8, max = 1.0)
coords[71:100, "y"] = runif(30, min = 0.8, max = 1.0)

#Visualisation rapide des positions. pch=1: cercle vide, pch=17:triangle plein
# asp=1 : meme ecartement sur les axes x et y
plot(coords[,1], coords[,2], col = c(rep("blue",70), rep("red",30)),
     pch = c(rep(1,70), rep(17,30)), asp = 1,
     main = "Positions spatiales des individus")

library(HDSpatialScan)

# Preparer les donnees pour le scan
score_vector = score_CP1         # vecteur des scores (numerique, long. 100)
coords_matrix = coords           # matrice des coordonnees (100 x 2)
# Executer le scan spatial non parametrique (WMW) avec 999 permutations Monte Carlo
# "UNP"= univariate nonparametric
#type_minimaxi = "sites/indiv" precise que ces valeurs mini/maxi 
#sont comptees en nombre de sites (individus).
res = SpatialScan(method = "UNP", data = score_vector, sites_coord = coords_matrix,
                   system = "Euclidean", 
                   mini = 1, maxi = n_total/2, type_minimaxi = "sites/indiv", 
                   MC = 999)

# Resultats du scan pour la méthode UNP
res_unp = res$UNP  # extraire l'objet de resultats pour la methode univariee non-param
res_unp
str(res_unp)


#UNP scan procedure 
#################### 
#2 significant clusters have been detected by the scan procedure with p-values of 0.001, 0.001 


summary(res_unp, type_summ = "nparam", only.MLC = TRUE)


#Cluster 1
#p-value     0.001
#Radius      0.168

#$complete_summary     
#Overall                           Inside cluster 1 Outside cluster 1
#Number of sites 100.000           30.000            70.000
#Q25              -5.318           -6.822             0.593
#Median            0.732           -6.822             2.344
#Q75               3.822           -6.201             5.058


library(ggplot2)
library(ggforce)

# indice du MLC (le plus probable, premier element)
members = res_unp$sites_clusters[[1]] # les membres du cluster 1

# centre (x0, y0) et rayon du MLC
center  = res_unp$centres_clusters[1, ]  # le centre du cluster 1
radius  = res_unp$radius_clusters[1] # le rayon du cluster 1    

# p‑valeur
pval_mlc = res_unp$pval_clusters[1] # la p-valu du cluster 1


# preparer le data.frame
df = data.frame(
  x         = res_unp$sites_coord[,1],
  y         = res_unp$sites_coord[,2],
  inMLC     = factor(seq_len(nrow(res_unp$sites_coord)) %in% members,
                     levels = c(FALSE, TRUE),
                     labels = c("hors MLC", "dans MLC"))
)

# tracer
# tous les points qui verifient : sqrt((x-0.8457784)^2+(y-0.8967244)^2)<=0.168
# seront dans le cluster
ggplot(df, aes(x, y, color = inMLC)) +
  geom_point(size = 2) +
  geom_circle(aes(x0 = center[1], y0 = center[2], r = radius),
              inherit.aes = FALSE,
              color = "black", linetype = "dashed") +
  coord_equal() +
  scale_color_manual(values = c("grey60", "red")) +
  labs(
    title    = "Most Likely Cluster (WMW Scan)",
    subtitle = paste0("p-value = ", round(pval_mlc, 3),
                      " | taille = ", length(members), " indiv.")
  ) +
  theme_minimal()





##### Implementation a la main #############

N = length(score_CP1)

#### 1. calculez les rangs une fois
ranks = rank(score_CP1, ties.method = "average")

# 1.2. matrice des distances euclidiennes
D = as.matrix(dist(coords))  # N×N


##### 2. Generation des fenetres candidates
# stocker les fenetres : liste de listes (centre, rayon, indices)

#Pour chaque site i (centre), on prend en candidat tous les cercles 
#de rayon egal a chacune des distances a un autre site j, triees par ordre 
#croissant. On s’arrete des que la fenetre contient N/2 individus.
windows = vector("list", N)# liste vide comportant N sites 

for(i in 1:N){
  # distances de i aux autres points
  d_i = D[i, ]
  # ordonne distances croissantes (y compris à soi = 0)
  ord = order(d_i)
  # cumul du nombre de points inclus
  wlist = list()
  for(j in ord){
    radius = d_i[j]
    inside = which(d_i <= radius)
    if(length(inside) > N/2) break
    # on ajoute a wlist la fenetre centree en i de rayon ‘radius’ contenant les indices ‘members’
    wlist[[length(wlist)+1]] = list(center=i, 
                                     radius=radius,
                                     members=inside)
  }
  # Stocke la liste des fenetres candidates generees pour le site i  
  windows[[i]] = wlist
}
# liste de toutes les fenetres candidates pour tous les sites
candidates = do.call(c, windows)  #100 centres (un par individu).Chacun 
#on stocke 50 fenetres. Ce qui donne 5000 fenetres 
candidates


##### 3.Calcul de la statistique WMW pour chaque fenetre
compute_pval = function(members, ranks, N){
  n_z = length(members)
  Wz = sum(ranks[members])
  Ez = n_z * (N+1)/2
  Vz = n_z*(N - n_z)*(N+1)/12
  Tz = (Wz - Ez) / sqrt(Vz)
  p_right = 1 - pnorm(Tz)  # unilateral a droite
  p_left  = pnorm(Tz)      # unilateral a gauche
  # ici on recherche cluster de scores faibles ⇒ unilat a gauche
  #return(p_right)
  return(p_left)
}

# calculer p_min et MLC en brut
best = list(pval=1, window=NULL) # conserver au fil de la boucle 
#la plus petite p‑valeur rencontree et la fenetre associee
for(w in candidates){
  p = compute_pval(w$members, ranks, N)
  if(p < best$pval){
    best$pval = p
    best$window = w
  }
}
# best$pval = p_min_obs=1.746323e-09, best$window = MLC
# Le most likely cluster est la fenetre centree sur le site n°15, 
#de rayon ~ 0,398, contenant 45 individus


#### 4. Monte Carlo (B = 999) ####
# On permute les rangs (ou directement les scores) B fois et pour chaque perm :

B = 999
pmin_sim = numeric(B)# vecteur de longueur B  qui contient les B 
# statistiques de scan calculees sur les donnees permutees
for(b in 1:B){
  perm_ranks = sample(ranks)   # permutation aleatoire
  pm = Inf # Initialiser la variable pm a l’infini,
  for(w in candidates){
    p   = compute_pval(w$members, perm_ranks, N)
    pm  = min(pm, p)
  }
  pmin_sim[b] = pm
}
# p‑valeur globale
R = sum(pmin_sim <= best$pval) + 1
p_global = R / (B + 1)



#### 5. Visualisation du MLC ##########


# extraire centre et rayon
i0     = best$window$center
radius = best$window$radius
members = best$window$members
center = coords[i0, ]

# base R
par(mfrow=c(1,2))
plot(coords, col = ifelse(seq_len(N) %in% members, "red","grey"),
     pch=16, xlab="x", ylab="y", main=paste0("MLC WMW p=",round(p_global,3)))
symbols(center[1], center[2], circles=radius,
        add=TRUE, inches=FALSE, lwd=2, lty=2)

# ou ggplot2 + ggforce

df = data.frame(coords, inMLC = seq_len(N) %in% members)
ggplot(df, aes(x,y,color=inMLC)) +
  geom_point() +
  geom_circle(aes(x0=center[1], y0=center[2], r=radius),
              inherit.aes=FALSE, linetype="dashed") +
  coord_equal() +
  scale_color_manual(values=c("grey","red"))+
  labs(
    title    = "Most Likely Cluster (WMW Scan)",
    subtitle = paste0("p-value = ", round(pval_mlc, 3),
                      " | taille = ", length(members), " indiv.")
  ) +
  theme_minimal()



