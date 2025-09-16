#------------------------------------------#
#-- Fit spaMM to simulated data --#
#------------------------------------------#
dir.out = "~/Simulation/"
library(tidyverse)
library(Seurat)
library(smiDE)
library(spaMM)
library(peakRAM)

# Function to get results -------------------------------------------------------------
get_spaMM_results <- function(gene, kprop, dist, method, data,Ycounts, Ygaus,scen){
  benchmark_and_return <- function(expr) {
    result <- NULL
    mem <- peakRAM(result <- eval(expr))
    list(result = result, memory = mem)
  }
  n <- nrow(data)
  p <- ncol(Ygaus)
  
  #---- Gets clusters ----
  k <- round(nrow(data) * kprop)
  if (kprop == 1) {
    data <- data %>%
      mutate(cluster = 1:nrow(data),
             sdimxclust = sdimx,
             sdimyclust = sdimy)
  } else {
    set.seed(4371)
    info_clust <- kmeans(data[, c('sdimx', 'sdimy')], k, iter.max = 50)
    
    cluster_assignments <- data.frame(
      cell_ID = data$cell_ID,
      cluster = as.integer(as.factor(info_clust$cluster))
    )
    
    centers_df <- data.frame(
      cluster = as.integer(as.factor(1:k)),
      sdimxclust = info_clust$centers[,1],
      sdimyclust = info_clust$centers[,2]
    )
    
    clust_info <- left_join(cluster_assignments, centers_df, by = "cluster")
    data <- left_join(data, clust_info, by = "cell_ID")
  }
  rownames(data) <- data$cell_ID
  
  #----- Run Models ----
  if(dist == "negbin"){
    dataclust <- cbind(Ycounts, data)
    dataclust$geneinterest <- dataclust %>% pull(gene)
    
    if(method == "IndepCluster"){
      mod_spa_ind <- benchmark_and_return(quote(
        fitme(geneinterest ~ offset(log(totalcounts)) + niche + (1 | cluster ),
              data = dataclust, family = negbin(),verbose=c(TRACE=F))
      ))
      coefs <- get_info_spamm(scen, gene, kprop, mod_spa_ind$result,"Spatial Indep", dist,
                              mod_spa_ind$memory$Elapsed_Time_sec,
                              mod_spa_ind$memory$Total_RAM_Used_MiB,
                              mod_spa_ind$memory$Peak_RAM_Used_MiB)
      
    }else if(method == "GP"){
      mod_spa_gp <- benchmark_and_return(quote(
        fitme(geneinterest ~ offset(log(totalcounts)) + niche + Matern(1 | sdimxclust + sdimyclust),
              data = dataclust, family = negbin(),verbose=c(TRACE=F),
              fixed=list(nu=0.5)
        )
      ))
      
      coefs <- get_info_spamm(scen, gene, kprop,mod_spa_gp$result ,"Spatial GP", dist,
                              mod_spa_gp$memory$Elapsed_Time_sec,
                              mod_spa_gp$memory$Total_RAM_Used_MiB,
                              mod_spa_gp$memory$Peak_RAM_Used_MiB)
    }
    
  }else if(dist == "gaussian"){
    dataclust <- cbind(Ycounts, data)
    dataclust$geneinterest <- dataclust %>% pull(gene)
    
    if(method == "IndepCluster"){
      mod_spa_ind <- benchmark_and_return(quote(
        fitme(geneinterest ~ niche + (1 | cluster ),
              data = dataclust, verbose=c(TRACE=F))
      ))
      coefs <- get_info_spamm(scen, gene, kprop, mod_spa_ind$result,"Spatial Indep", dist,
                              mod_spa_ind$memory$Elapsed_Time_sec,
                              mod_spa_ind$memory$Total_RAM_Used_MiB,
                              mod_spa_ind$memory$Peak_RAM_Used_MiB)
    }else if(method == "GP"){
      mod_spa_gp <- benchmark_and_return(quote(
        fitme(geneinterest ~ niche + Matern(1 | sdimxclust + sdimyclust),
              data = dataclust, verbose=c(TRACE=F),
              fixed=list(nu=0.5)
        )
      ))
      
      coefs <- get_info_spamm(scen, gene, kprop, mod_spa_gp$result ,"Spatial GP", dist,
                              mod_spa_gp$memory$Elapsed_Time_sec,
                              mod_spa_gp$memory$Total_RAM_Used_MiB,
                              mod_spa_gp$memory$Peak_RAM_Used_MiB)
    }
  }
  
  coefs
}

# Function to extract spaMM results
get_info_spamm <- function(scen, gene, kprop, model, name, dist,time,totalRAM,peakRAM){
  coefs <- summary(model)$beta_table[,1:2]
  nu <- model$ranFix$corrPars$`1`$nu
  rho <- model$ranFix$corrPars$`1`$rho
  phi <- model$phi
  lambda <- model$lambda
  coefs <- as.data.frame(cbind(scen, gene,kprop,name,dist, rownames(coefs),
                               coefs,NA,NA,
                               ifelse(length(nu) == 0, NA, nu),
                               ifelse(length(rho) == 0, NA, rho),
                               ifelse(length(phi) == 0, NA, phi),
                               ifelse(length(lambda) == 0, NA, lambda),
                               time,totalRAM,peakRAM, NA
  ))
  colnames(coefs) <- c("Scenario","Gene","kprop","Model","Distribution","Variable","Estimate","SE",
                       "LB","UB","Nu","Rho","Phi","Lambda","Time","TotalRAM","PeakRAM","pvalue")
  coefs <- coefs %>% mutate(LB = as.numeric(Estimate) - qnorm(0.975) * as.numeric(SE),
                            UB = as.numeric(Estimate) + qnorm(0.975) * as.numeric(SE),
                            pvalue = 2*pnorm(-abs(as.numeric(Estimate)/as.numeric(SE))))
  return(coefs)
}

# Read data ---------------------------------------------------------------
load(paste0(dir.out,paste0("SimulatedY.Rdata")))
n <- nrow(data)
data$cell_ID <- as.character(data$cell_ID)

# Run ---------------------------------------------------------------------
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)
simparams <- expand_grid(gene = colnames(Ycounts),
                         method = c("IndepCluster","GP"),
                         kprop = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                                   0.08, 0.09,  0.1, 0.15,  0.2,
                                   0.25, 0.3,   0.4,   0.5,  1),
                         dist = c("negbin","gaussian")
) 

results <- get_spaMM_results(gene=simparams$gene[args],
                             kprop=simparams$kprop[args],
                             dist=simparams$dist[args], 
                             method = simparams$method[args],
                             data,Ycounts, Ygaus,scen)

write.csv(results, file = paste0(dir.out,paste("Results/PartialSimulation/spaMM/spaMM",
                                               paste(sapply(simparams[args, ], as.character), collapse = "_") ,".csv",sep="_")))

