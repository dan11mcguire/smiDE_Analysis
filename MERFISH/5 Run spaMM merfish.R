#--------------------#
#- Run spaMM Models -#
#--------------------#

data_dir <- "~MERFISH/"
ncores<- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')) -1

library(tidyverse)
library(Seurat)
library(smiDE)
library(peakRAM)
library(future)
library(future.apply)
library(emmeans)
library(broom)
library(parallel)
library(spaMM)
library(furrr)
benchmark_and_return <- function(expr) {
  result <- NULL
  mem <- peakRAM(result <- eval(expr))
  list(result = result, memory = mem)
}
RankNorm <-function (u, k = 0.375){
  n <- length(u)
  r <- rank(u)
  out <- qnorm((r - k)/(n - 2 * k + 1))
  return(out)
}


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
params <- expand_grid(kprop = c(0.05,0.25,0.5,1),
                      methodrun = c("indep_full_nb","gp_full_nb","indep_full_gaus","gp_full_gaus"), 
                      var = c("kernelMicroglia1000"))
methodrun <- params$methodrun[as.numeric(slurm_arrayid)]
kprop <- params$kprop[as.numeric(slurm_arrayid)]
var <- params$var[as.numeric(slurm_arrayid)]

load(paste0(data_dir, "Allen/WholeBrain_seurat_distance.Rdata"))
load(file = paste0(data_dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
filteredgenes <- overlap_metrics_allen$result %>% filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1, ratio < 1) %>%pull(target)

# Set functions -----------------------------------------------------------

get_spaMM_results <- function(var, data,kprop,type, dist, method,neighbor_expr_list,
                              targets=NULL,nCores){
  data@meta.data$DEvar <- data@meta.data[,var]
  
  p <- nrow(data@assays$RNA)
  n <- ncol(data@assays$RNA)
  #---- Gets clusters ----
  k <- round(nrow(data) * kprop)
  if(kprop == 1){
    data@meta.data <- data@meta.data %>% dplyr::mutate(cluster = 1:nrow(data@meta.data),
                                                       sdimx = x,
                                                       sdimy = y)
  }else{
    set.seed(4371)
    info_clust <- kmeans(data@meta.data[,c('x','y')], k, iter.max = 50)
    colnames(info_clust$centers) <- c("sdimx","sdimy")
    
    data@meta.data <- data@meta.data %>% dplyr::mutate(cluster = info_clust$cluster) %>%
      left_join(data.frame(info_clust$centers, cluster = c(1:k)),
                by = "cluster")
  }
  rownames(data@meta.data) <- data$cell_ID
  
  if(is.null(targets))targets <- rownames(data@assays$RNA$counts)
  
  #----- Define function ----
  getspamm <- function(var, target, data_target,assay, neighbor_expr_list,
                       dist, type, kprop, method){
    cat("Starting target:", target, "\n")
    data_target$otherct_expr <- RankNorm(as.numeric(neighbor_expr_list[target,data_target$cell_ID]))
    
    if(dist == "negbin"){
      if(type == "onlyspatial"){
        form <- "y ~ DEvar + offset(log(nCount_RNA)) "
      }else if(type == "full"){
        form <- "y ~ DEvar + otherct_expr + offset(log(nCount_RNA))"
      }
      
      data_target$y <- assay$counts[target,]
      
      if(method == "IndepCluster"){
        mod_spamm <- tryCatch({
          spaMM::fitme(formula(paste(form, "(1 | cluster )", sep = "+")),
                       data = data_target, family = "negbin",verbose=c(TRACE=F))
        }, error = function(e) {
          return(NULL)
        })
        
      }else if(method == "GP"){
        mod_spamm <- tryCatch({
          spaMM::fitme(formula(paste(form, "Matern(1 | sdimx + sdimy)", sep = "+")),
                       data = data_target, family = "negbin",verbose=c(TRACE=F),
                       fixed=list(nu=0.5))
        }, error = function(e) {
          return(NULL)
        })
      }
      
    }else if(dist == "gaussian"){
      if(type == "onlyspatial"){
        form <- "y ~ DEvar"
      }else if(type == "full"){
        form <- "y ~ DEvar + otherct_expr"
      }
      
      data_target$y <- assay$data[target,]
      if(method == "IndepCluster"){
        mod_spamm <- tryCatch({
          spaMM::fitme(formula(paste(form, "(1 | cluster )", sep = "+")),
                       data = data_target, verbose=c(TRACE=F))
        }, error = function(e) {
          return(NULL)
        })
      }else if(method == "GP"){
        mod_spamm <- tryCatch({
          spaMM::fitme(formula(paste(form, "Matern(1 | sdimx + sdimy)", sep = "+")),
                       data = data_target,verbose=c(TRACE=F),
                       fixed=list(nu=0.5))
        }, error = function(e) {
          return(NULL)
        })
      }
    }
    if (!is.null(mod_spamm)) {
      coefs <- summary(mod_spamm)$beta_table
      
      nu <- mod_spamm$ranFix$corrPars$`1`$nu
      rho <- mod_spamm$ranFix$corrPars$`1`$rho
      phi <- mod_spamm$phi
      lambda <- mod_spamm$lambda
      coefs <- as.data.frame(cbind(var,target,dist, type, kprop,method,
                                   rownames(coefs), coefs,NA,NA,NA,
                                   ifelse(length(nu) == 0, NA, nu),
                                   ifelse(length(rho) == 0, NA, rho),
                                   ifelse(length(phi) == 0, NA, phi),
                                   ifelse(length(lambda) == 0, NA, lambda),
                                   FALSE
      ))
      colnames(coefs) <- c("DEvar","Gene","Distribution","Type","kprop","Method","Variable","Estimate","SE","tvalue",
                           "LB","UB","pvalue","nu","rho","phi","lambda","failed")
      coefs <- coefs %>% mutate(LB = as.numeric(Estimate) - qnorm(0.975) * as.numeric(SE),
                                UB = as.numeric(Estimate) + qnorm(0.975) * as.numeric(SE),
                                pvalue = 2*pnorm(-abs(as.numeric(Estimate)/as.numeric(SE))))
    }else{
      coefs <- data.frame(DEvar=var,
                          Gene = target,
                          Distribution = dist,
                          Type = type,
                          kprop = kprop,
                          Method = method,
                          Variable = NA,
                          Estimate = NA,
                          SE = NA,
                          tvalue = NA,
                          LB = NA,
                          UB = NA,
                          pvalue = NA,
                          nu = NA,
                          rho = NA,
                          phi = NA,
                          lambda = NA,
                          failed = T
      )
    }
    return(coefs)
  }
  
  #----- Run Models ----
  if(nCores > 1){
    neighbor_expr_list <- as.matrix(neighbor_expr_list)
    
    plan(multisession, workers = nCores)
    results <- future_map(targets, function(target) {
      getspamm(var=var,
               target = target,
               data_target = data@meta.data,
               assay = data@assays$RNA,
               neighbor_expr_list = neighbor_expr_list,
               method = method,
               dist = dist,
               type = type,
               kprop = kprop
      )
    }, .options = furrr_options(seed = TRUE))
    plan("sequential")
  }else{
    results <- lapply( targets, function(target) {
      getspamm(var=var,
               target = target,
               data_target = data@meta.data,
               assay = data@assays$RNA,
               neighbor_expr_list = neighbor_expr_list,
               method = method,
               dist = dist,
               type = type,
               kprop = kprop
      )
    })
  }
  results <- do.call(rbind, results)
  results
}



# Get results -------------------------------------------------------------
# Indep Clusters ----------------------------------------------------------
# Only spatial
# Negative Binomials
if(methodrun == "indep_spatial_nb"){
  indep_spatial_nb <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "onlyspatial", dist="negbin",
                      method="IndepCluster",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      targets = filteredgenes,
                      nCores = ncores)
  ))
  indep_spatial_nb <- indep_spatial_nb$result %>%
    mutate(Time =  indep_spatial_nb$memory$Elapsed_Time_sec,
           TotalRAM= indep_spatial_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  indep_spatial_nb$memory$Peak_RAM_Used_MiB)
  write.csv(indep_spatial_nb, file = paste0("ResultsDistance/spaMM/spaMM_indep_spatial_nb",var,"_kprop",kprop,".csv"))
}

# Gaussian
if(methodrun == "indep_spatial_gaus"){
  indep_spatial_gaus <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "onlyspatial", dist="gaussian",
                      method="IndepCluster",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      # targets = rownames(dataDE@assays$RNA$counts)[1:2],
                      nCores = ncores)
  ))
  indep_spatial_gaus <- indep_spatial_gaus$result %>%
    mutate(Time =  indep_spatial_gaus$memory$Elapsed_Time_sec,
           TotalRAM= indep_spatial_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  indep_spatial_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(indep_spatial_gaus, file = paste0("ResultsDistance/spaMM/spaMM_indep_spatial_gaus",var,"_kprop",kprop,".csv"))
}
# Full
# Negative Binomial
if(methodrun == "indep_full_nb"){
  indep_full_nb <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "full", dist="negbin",
                      method="IndepCluster",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      targets = filteredgenes,
                      nCores = ncores)
  ))
  indep_full_nb <- indep_full_nb$result %>%
    mutate(Time =  indep_full_nb$memory$Elapsed_Time_sec,
           TotalRAM= indep_full_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  indep_full_nb$memory$Peak_RAM_Used_MiB)
  write.csv(indep_full_nb, file = paste0("ResultsDistance/spaMM/spaMM_indep_full_nb",var,"_kprop",kprop,".csv"))
}

# Gaussian
if(methodrun == "indep_full_gaus"){
  indep_full_gaus <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "full", dist="gaussian",
                      method="IndepCluster",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      targets = filteredgenes,
                      nCores = ncores)
  ))
  indep_full_gaus <- indep_full_gaus$result %>%
    mutate(Time =  indep_full_gaus$memory$Elapsed_Time_sec,
           TotalRAM= indep_full_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  indep_full_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(indep_full_gaus, file = paste0("ResultsDistance/spaMM/spaMM_indep_full_gaus",var,"_kprop",kprop,".csv"))
}

# GP Clusters ----------------------------------------------------------
# Only spatial
# Negative Binomials
if(methodrun == "gp_spatial_nb"){
  gp_spatial_nb <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "onlyspatial", dist="negbin",
                      method="GP",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      targets = filteredgenes,
                      nCores = ncores)
  ))
  gp_spatial_nb <- gp_spatial_nb$result %>%
    mutate(Time =  gp_spatial_nb$memory$Elapsed_Time_sec,
           TotalRAM= gp_spatial_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  gp_spatial_nb$memory$Peak_RAM_Used_MiB)
  write.csv(gp_spatial_nb, file = paste0("ResultsDistance/spaMM/spaMM_gp_spatial_nb",var,"_kprop",kprop,".csv"))
}

# Gaussian
if(methodrun == "gp_spatial_gaus"){
  gp_spatial_gaus <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "onlyspatial", dist="gaussian",
                      method="GP",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      targets = filteredgenes,
                      nCores = ncores)
  ))
  gp_spatial_gaus <- gp_spatial_gaus$result %>%
    mutate(Time =  gp_spatial_gaus$memory$Elapsed_Time_sec,
           TotalRAM= gp_spatial_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  gp_spatial_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(gp_spatial_gaus, file = paste0("ResultsDistance/spaMM/spaMM_gp_spatial_gaus",var,"_kprop",kprop,".csv"))
}

# Full
# Negative Binomial
if(methodrun == "gp_full_nb"){
  gp_full_nb <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "full", dist="negbin",
                      method="GP",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      targets = filteredgenes,
                      nCores = ncores)
  ))
  gp_full_nb <- gp_full_nb$result %>%
    mutate(Time =  gp_full_nb$memory$Elapsed_Time_sec,
           TotalRAM= gp_full_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  gp_full_nb$memory$Peak_RAM_Used_MiB)
  write.csv(gp_full_nb, file = paste0("ResultsDistance/spaMM/spaMM_gp_full_nb",var,"_kprop",kprop,".csv"))
}

# Gaussian
if(methodrun == "gp_full_gaus"){
  gp_full_gaus <- benchmark_and_return(quote(
    get_spaMM_results(var=var,data=dataDE,kprop = kprop,type = "full", dist="gaussian",
                      method="GP",
                      neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                      targets = filteredgenes,
                      nCores = ncores)
  ))
  gp_full_gaus <- gp_full_gaus$result %>%
    mutate(Time =  gp_full_gaus$memory$Elapsed_Time_sec,
           TotalRAM= gp_full_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  gp_full_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(gp_full_gaus, file = paste0("ResultsDistance/spaMM/spaMM_gp_full_gaus",var,"_kprop",kprop,".csv"))
}