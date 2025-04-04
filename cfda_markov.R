library(cfda)
library(ggplot2)

n_total = 400       # Nombre total d'individus
n_noncluster = 300  # Nombre d'individus hors cluster
n_cluster = 100     # Nombre d'individus dans le cluster
K = 3               # Nombre d'etats, par exemple "A", "B", "C"
Tmax = 12           # Duree maximale (par exemple, 12 mois)
labels = c("A", "B", "C") # Noms des etats

# Parametres pour le groupe non-cluster

#matrice de transition.
#Chaque ligne de cette matrice represente les probabilites de passer 
#de l’etat courant a chacun des etats possibles, et la somme des 
#probabilites de chaque ligne est egale a 1.
P1 = matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
lambda1 = rep(1, K) #time spent in each state
pi0_1 = c(1, 0, 0) # Tous les individus commencent en etat "A"

# Generer les donnees pour le groupe non-cluster
set.seed(42)
data_noncluster = generate_Markov(n = n_noncluster, K = K, P = P1, lambda = lambda1,
                                  pi0 = pi0_1, Tmax = Tmax, labels = labels)

# Parametres pour le groupe cluster : 
# Par exemple, on modifie la matrice de transition pour induire un decalage dans les trajectoires
P2 =  matrix(1 / 5, nrow = K, ncol = K) - diag(rep(1 / 5, K))
lambda2 = rep(1, K)   # meme duree
pi0_2 = c(1, 0, 0)

# Generer les donnees pour le groupe cluster
data_cluster = generate_Markov(n = n_cluster, K = K, P = P2, lambda = lambda2,
                               pi0 = pi0_2, Tmax = Tmax, labels = labels)


# Je decale les ids de data_cluster pour qu'ils ne se chevauchent pas avec ceux de data_noncluster.
data_cluster$id = data_cluster$id + max(data_noncluster$id)

# Combiner les donnees
data_total = rbind(data_noncluster, data_cluster)
summary_cfd(data_total)

head(data_total,20)


#visualiser un echantillon des trajectoires
plotData(data_total[data_total$id <= 50, ])


#Calcul de la Duree de Suivi pour chaque individu.
duration = compute_duration(data_total)

#couper les trajectoires pour qu'elles se terminent a Tmax=10
data_cut = cut_data(data_total, 10)

# Creation de la base B-spline pour representer les trajectoires
basis = create.bspline.basis(c(0, 10), nbasis = 10, norder = 4)

# Calcul de l'encodage optimal
fmca = compute_optimal_encoding(data_cut, basis, nCores = 1)
summary(fmca)

# Visualisation des valeurs propres
plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE)
plot(fmca)

# Selection des composantes principales pour expliquer 90 % de la variance
nPc90 = which(cumsum(prop.table(fmca$eigenvalues)) > 0.9)[1]
cat("Nombre de composantes sélectionnées pour 90% de la variance :", nPc90, "\n")

############### Clustering

# Extraction des scores correspondants

pc_scores = fmca$pc[, 1:nPc90]

# Clustering hirarchique sur les scores
hc = hclust(dist(pc_scores), method = "average")
plot(hc, labels = FALSE, main = "Dendrogramme des clusters", ylab = "Hauteur")

# Decouper le dendrogramme pour obtenir 2 clusters 
clusters = cutree(hc, k = 2)



# Visualiser les trajectoires selon les clusters obtenus
plotData(data_cut, group = clusters, addId = FALSE, addBorder = FALSE, sort = TRUE)+
  ggtitle("Trajectoires par clusters (méthode average)")



