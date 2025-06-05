#### PACKAGE #####
library(sp)
library(cfda)
library(HDSpatialScan)

### Analyse de l'objet " departements" ###

summary(departements)
str(departements)
class(departements)

### Code des departements ###
code_dept <- as.numeric(as.character(departements@data$CODE_DEPT)) # 94

### Recuperation des departs d'Ile de France ###
idx <- which(departements@data$NOM_REG == "ILE-DE-FRANCE")
code_dept_clus <- code_dept[idx]
## c(91, 93, 92, 95, 94, 77, 75, 78)

## Code des departs hors cluster
code_dept_out <- setdiff(code_dept, code_dept_clus)

### Nombre d'individus
n_ind_per_dept <- 10
n_clus <- length(code_dept_clus)* n_ind_per_dept # n_clus = 80
n_out <- length(code_dept_out)* n_ind_per_dept # n_out = 860

### PARAMETRES ###

k <- 3
Tmax <- 18

### Matrices de transition ###
## SIMU 1 ##

### Hors cluster ###
P1_hors_clus <- matrix(0, nrow = k, ncol = k)
P1_hors_clus[1, 2] <- 0.5; P1_hors_clus[1, 3] <- 0.5
P1_hors_clus[2, 1] <- 0.5; P1_hors_clus[2, 3] <- 0.5
P1_hors_clus[3, 1] <- 0.5; P1_hors_clus[3, 2] <- 0.5

### Weak ###
P1_cluster_weak <- P1_hors_clus

### Strong ###
P1_cluster_strong <- P1_hors_clus

### SIMU 2 ###

### Hors cluster ###
P2_hors_clus <- matrix(0, nrow = k, ncol = k)
P2_hors_clus[1, 2] <- 1; P2_hors_clus[2, 1] <- 0.5
P2_hors_clus[2, 3] <- 0.5; P2_hors_clus[3, 1] <- 0.5
P2_hors_clus[3, 2] <- 0.5

### Weak ###
P2_cluster_weak <- P2_hors_clus

### Strong ###
P2_cluster_strong <- P2_hors_clus


## taux de sejour Lambda ##
## Simu 1 ##

### hors cluster ###
lambda1_hors_clus <- c(0.2, 0.2, 0.2)

### Weak ###
lambda1_cluster_weak <- c(1, 1, 1)

### Strong ###
lambda1_cluster_strong <- c(2, 2, 2)

## Simu 2 ##

### hors cluster ###
lambda2_hors_clus <- c(0.2, 0.2, 1e-7)

### Weak ###
lambda2_cluster_weak <- c(0.3, 0.3, 1e-7)

### Strong ###
lambda2_cluster_strong <- c(0.4, 0.4, 1e-7 )



### definition de la fonction ###
detect_cluster <- function(scenario, intensity, simu_id)
{
  set.seed(20 + simu_id)  # Graine différente à chaque simulation
  
    # Choix des matrices de transition et des lambdas appropriés pour chaque scénario et chaque intensité
     if(scenario == "Simu1") {
       P_out <- P1_hors_clus
       P_in  <- P1_hors_clus 
       lambda_out  <- lambda1_hors_clus
       lambda_in <- if (intensity == "weak") lambda1_cluster_weak else lambda1_cluster_strong

    } else if (scenario == "Simu2") {
        P_out <- P2_hors_clus
        P_in <-  P2_hors_clus
        lambda_out  <- lambda2_hors_clus
        lambda_in <- if (intensity == "weak") lambda2_cluster_weak else lambda2_cluster_strong
    }


       ## 1. Simulation des trajectoires
       data_cluster <- generate_Markov(n = n_clus, K = k, P = P_in, 
                                lambda = lambda_in, Tmax = Tmax)
       data_out <- generate_Markov(n = n_out, K = k, P = P_out, 
                            lambda = lambda_out, Tmax = Tmax)


       # Reindexer les id du groupe hors_cluster (pour eviter les doublons)
       data_out$id <- data_out$id + n_clus

      # Combinez les données
       data_all <- rbind(data_cluster, data_out)
       
       ## selection des individus pour le codage ##
       data_cut <- cut_data(data_all, Tmax)
       
       ### creation de la base B-spline ###
       basis <- create.bspline.basis(c(0, Tmax), nbasis = 10, norder = 4)
       
       ### Encodage ###
       fmca <- compute_optimal_encoding(data_cut, basis, nCores = 20)
       
       ### Recuperation des valeurs propres ###
       VP <- fmca$eigenvalues
       
       ### Scores des individus ###
       scores <- fmca$pc
       
       ### fixer l'inertie a 90 ###
       q_90 <- which(cumsum(prop.table(fmca$eigenvalues)) > 0.9)[1]
       
       ### Preparation des data de la function SpatialScan ###
       scores_q <- scores[, 1 : q_90, drop = FALSE]
       
       ## Coordonnees des departs ###
       coords <- coordinates(departements)
       
       ### coordonnees des individus ###
       all_dept <- c(rep(code_dept_clus, each = n_ind_per_dept), rep(code_dept_out, each = n_ind_per_dept))
       coord_ind <- coords[match(all_dept, code_dept),]
       
       
       #### Application du scan spatial MNP ###
       res_scan <- SpatialScan(method = "MNP", data = scores_q, sites_coord = coord_ind, system = "Euclidean", mini = 1, maxi = nrow(coord_ind)/2, MC = 99, nbCPU = 20)
       
       ### Extraction du MLC ###
       if(length(res_scan$MNP$sites_clusters) > 0)
       {
         mlc_sites <- res_scan$MNP$sites_clusters[[1]]
         p_val <- res_scan$MNP$pval_clusters[1]
       }
       else
       {
         mlc_sites <- integer(0)
         p_val <- NA
       }
       mlc_depts <- unique(all_dept[mlc_sites])
       
       # IDs des individus détectés
       id_detectes <- mlc_sites
       
       # Département d'appartenance pour chaque individu détecté
       dept_individus_detectes <- all_dept[mlc_sites]
       
       
       
       ### Enregistrement ###
       results <- list(
         scenario = scenario,
         intensity = lambda_in,
         simulation = simu_id,
         dept_simules = code_dept_clus,
         dept_detectes = mlc_depts,
         individus_detectes = length(mlc_depts)* n_ind_per_dept,
         ids_detectes = id_detectes,
         dept_par_individu = table(dept_individus_detectes),
         pvalue = p_val
         
       )
       file_path <- file.path("results_simulation/", paste0("simu_", scenario, "_", intensity, "_", simu_id, ".rds"))
       
       ## enregister results a l'endroit indique par file_path ##
       saveRDS(results, file = file_path)
       
       return(results)
}


       
    
 
    
   
