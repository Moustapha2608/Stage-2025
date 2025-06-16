library(cfda)
library(fda)
set.seed(42)
# Simulate the Jukes-Cantor model of nucleotide replacement
K <- 3
Tmax <- 1
d_JK <- generate_Markov(n = 100, K = K, Tmax = Tmax)
d_JK2 <- cut_data(d_JK, Tmax)
# create basis object
m <- 10
base <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
# compute encoding
encoding <- compute_optimal_encoding(d_JK2, base, computeCI = FALSE, nCores = 8)

phi    <- eval.basis(encoding$pt$t,base) # matrice   96 x m
coeffs <- encoding$alpha # length(coeffs) = 30
probas <- encoding$pt$pt # probas marginales : matrice K x 96
temps  <- encoding$pt$t  # time values des donnees
scores      <- encoding$pc    # composantes principales 100 x 30

nComp <- length(encoding$eigenvalues)   # Nombres de composantes principales

Tn  <- nrow(phi)             # 96

# probas doit être de dimension T × K :
#if (!all(dim(probas) == c(Tn,K))) probas <- t(probas)

############ Calcul des fonctions propres
A <- array(0, c(Tn, K, nComp), dimnames = list(NULL, paste0("state",1:K), paste0("PC",1:nComp)))
for(i in seq_len(nComp)) {
  # coeffs[[i]] est une matrix m×K, chaque colonne x donne α(·,x) pour la iᵉ composante
  A[,,i] <- phi %*% coeffs[[i]]  # donne T × K
}

# a_i^x : encoding = \sum_j \alpha_{(x, j)} * phi_j(t)
#a_i_t_x <- list()
#for (j in seq_len(nComp)) {
 # a <- fd(encoding$alpha[[j]], encoding$basisobj)
 # a_i_t_x[[j]] <- eval.fd(temps, a)
#}

# Compute p^x(t) =  1 / (sum_{i \geq 1} a_i^x(t)^2 * lambda_i)
denominator <- 1
for (j in seq_len(nComp)) {
  denominator <- denominator + A[,,j]^2 * encoding$eigenvalues[j]
}
pt <- 1 / denominator
##########   ∑_{i=1..M} Z_i * a^x_i(t)/ p^x(t)
n_ind <- nrow(scores)
reconstruction1  <- array(NA, c(n_ind, Tn, K))

##############  ∑_{i=1..M} Z_i * a^x_i(t) * p^x(t)
reconstruction2 <- array(NA, c(n_ind, Tn, K))

##############  ∑_{i=1..M} Z_i * a^x_i(t) * p^x(t) + p^x(t)
reconstruction3 <- array(NA, c(n_ind, Tn, K))

for(i in 1:n_ind) {
  Zi <- scores[i, 1:nComp]          # vecteur M
  for(x in 1:K) {
    num1 <- A[, x, ] %*% Zi   # vecteur longitud T
    P <- pt[, x]       # p^x(t)
    reconstruction1[i, , x] <- num1 / P #  la reconstruction de l’indicatrice 1_X_i(t)=x
    reconstruction2[i, , x] <- P  * num1
    reconstruction3 [i, , x] <- P  * num1 + P
  }
}



############# (∑_{i=1..M} Z_i * a^x_i(t) + 1 )* p^x(t)
ri <- reconstructIndicators(encoding, nComp)
ri# ri est un data.frame long avec colonnes time, id, state1, state2, state3, state
# on transforme en array n_ind × T × K
reconstruction4 <- array(NA, c(n_ind, Tn, K))
for(x in 1:K) {
  mat <- matrix(ri[[paste0("state",x)]], n_ind, Tn, byrow = TRUE)
  reconstruction4[,,x] <- mat
}

###########
pred1 <- apply(reconstruction1, c(1,2), which.max) 
pred2 <- apply(reconstruction2, c(1,2), which.max)
pred3 <- apply(reconstruction3, c(1,2), which.max)
pred4 <- apply(reconstruction3, c(1,2), which.max)

# construire la vérité : matrix n_ind × Tn  
# à partir des données originales d_JK2
true_mat <- matrix(NA, n_ind, Tn)
for(i in 1:n_ind) {
  # extraire la trajectoire i et remplir true_mat[i, ] par l’état à chaque temps
  sub   <- d_JK2[d_JK2$id == i, ]        # les événements de l'individu i
  times <- sub$time
  states<- sub$state
  
  # Pour chaque t dans `temps`, findInterval renvoie l'indice du plus grand
  # times[j] <= t. Si findInterval = 0, on prend l'état initial (states[1]).
  idx   <- findInterval(temps, times)
  
  # Si jamais idx[k] == 0 (temps[k] < times[1]), on reste dans l'état states[1]
  # tant qu'il n'y a pas de saut, on reste dans l'état initial
  idx[idx == 0] <- 1
  
  # On remplit la i-ème ligne avec l'état correspondant :
  true_mat[i, ] <- states[idx]
}


acc1 <- mean(pred1 == true_mat)
acc2 <- mean(pred2 == true_mat)
acc3 <- mean(pred3 == true_mat)
acc4 <- mean(pred4 == true_mat)

df <- data.frame(
  formule  = c("1/p * sum(Z a)", "sum(a p Z)", "reconstructIndicators", "sum(a p Z) + p" ),
  accuracy = c(acc1, acc2, acc3,acc4)
)



#pred <- apply(indicators[,paste0("state",1:3)], 1, which.max)
#mean(pred == indicators$state)  # taux de bonnes reconstitutions


#)


library(ggplot2)

# Tracer le graphique
ggplot(df, aes(x = formule, y = accuracy * 100, fill = formule)) +
  geom_col(width = 0.6) +
  labs(
    title = "Comparaison des méthodes de reconstruction de l’indicatrice",
    x = "Méthode",
    y = "Taux de précision (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0, 100)

