#---------------------#
#- Run BYM2 Models -#
#---------------------#
# BYM2
data_dir <- "~MERFISH/"
ncores<- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')) -1

library(tidyverse)
library(peakRAM)
library(future)
library(future.apply)
library(furrr)
library(spData)
library(spdep)
library(INLA)
library(fmesher)
library(broom)
library(parallel)

benchmark_and_return <- function(expr) {
  result <- NULL
  mem <- peakRAM(result <- eval(expr))
  list(result = result, memory = mem)
}

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
params <- expand_grid(kprop = c(0.05,0.25,0.5,1),
                      methodrun = c("full_gaus","full_nb"),
                      prior_phi = c("iid","spatial","balanced"),
                      prior_prec = c("conservative","moderate","standard","permissive","flat"),
                      adjmethod = "delauneytrig",
                      kadj = c(NA),
                      var = c("kernelMicroglia1000")
)
methodrun <- params$methodrun[as.numeric(slurm_arrayid)]
kprop <- params$kprop[as.numeric(slurm_arrayid)]
prior_phi <- params$prior_phi[as.numeric(slurm_arrayid)]
prior_prec <- params$prior_prec[as.numeric(slurm_arrayid)]
adjmethod <- params$adjmethod[as.numeric(slurm_arrayid)]
kadj <- params$kadj[as.numeric(slurm_arrayid)]
var <- params$var[as.numeric(slurm_arrayid)]


load(file = paste0(data_dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
load(paste0(data_dir, "Allen/WholeBrain_df_distance.Rdata"))
filteredgenes <- overlap_metrics_allen$result %>% filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1, ratio < 1) %>%pull(target)


# Set functions -----------------------------------------------------------
get_bym2_results <- function(var,metadata, counts, norm,kprop,type, dist, meshscenario,neighbor_expr_list,
                             prior_phi, prior_prec, adjmethod, kadj,targets=NULL,nCores){
  metadata$DEvar <- metadata[,var]
  # load inside to avoid issues when in parallel
  require(INLA)
  require(spdep)
  library(Matrix)
  RankNorm <-
    function (u, k = 0.375)
    {
      n <- length(u)
      r <- rank(u)
      out <- stats::qnorm((r - k)/(n - 2 * k + 1))
      return(out)
    }
  n <- nrow(metadata)
  #---- Gets clusters ----
  k <- round(nrow(metadata) * kprop)
  if (kprop == 1) {
    metadata <- metadata %>%
      mutate(clusterloc = as.integer(as.factor(1:nrow(metadata))),
             sdimx = x,
             sdimy = y)
  } else {
    set.seed(4371)
    info_clust <- kmeans(metadata[, c("x", "y")], centers = k, iter.max = 50)
    
    # Assign consistent cluster IDs
    cluster_assignments <- data.frame(
      cell_ID = metadata$cell_ID,
      clusterloc = as.integer(as.factor(info_clust$cluster))
    )
    
    centers_df <- data.frame(
      clusterloc = 1:k,
      sdimx = info_clust$centers[,1],
      sdimy = info_clust$centers[,2]
    )
    
    clust_info <- left_join(cluster_assignments, centers_df, by = "clusterloc")
    metadata <- left_join(metadata, clust_info, by = "cell_ID")
  }
  
  
  #---- Setting parameters ----
  # Create adjacency graph between centroids
  # Delaunay triangulation
  if(adjmethod == "delauneytrig"){
    nb <- metadata %>% dplyr::select(clusterloc, sdimx, sdimy) %>%
      distinct %>% arrange(clusterloc) %>% dplyr::select(sdimx,sdimy) %>% distinct %>% spdep::tri2nb()
  }else if(adjmethod == "knn"){
    auxmeta <- metadata %>% dplyr::select(clusterloc, sdimx, sdimy) %>%
      distinct %>% arrange(clusterloc)
    kn <- knearneigh(auxmeta[,c("sdimx","sdimy")], k = kadj)
    nb <- knn2nb(kn, row.names = auxmeta$clusterloc)
  }
  adjmat <- nb2mat(nb, style = "B", zero.policy = TRUE)
  graph_path <- Matrix(adjmat, sparse = TRUE)
  
  if(is.null(targets))targets <- rownames(counts)
  #----- Define function ----
  getbym2 <- function(var, target, data_target,counts,norm, neighbor_expr_list,
                      mesh, dist, type, kprop, prior_phi, prior_prec, adjmethod, kadj,graph_path){
    RankNorm <-
      function (u, k = 0.375)
      {
        n <- length(u)
        r <- rank(u)
        out <- stats::qnorm((r - k)/(n - 2 * k + 1))
        return(out)
      }
    
    INLA::inla.setOption(num.threads = "1")
    cat(paste0("Fitting model to target ",target, "."))
    data_target$otherct_expr <- RankNorm(neighbor_expr_list[target,data_target$cell_ID])
    
    if(dist == "negbin"){
      data_target$y <- counts[target,]
      
      prior_phi_values <- dplyr::case_when(
        prior_phi == "iid" ~ c(0.3, 0.5),       # Shrinks toward i.i.d. (often better for NegBin as it interacts with theta)
        prior_phi == "spatial" ~ c(0.7, 0.5),   # Shrinks toward ICAR
        prior_phi == "balanced" ~ c(0.5, 0.5)   # 50% chance it is greater than 0.5
      )
      prior_prec_values <- dplyr::case_when(
        prior_prec == "conservative" ~ c(0.3, 0.01),   # shrink to σ < 0.3 → strong control
        prior_prec == "moderate"     ~ c(0.5, 0.01),   # shrink to σ < 0.5 → safe for NB
        prior_prec == "standard"     ~ c(0.75, 0.01),  # balanced, allows spatial variation
        prior_prec == "permissive"   ~ c(1.0, 0.01),   # wider variation
        prior_prec == "flat"         ~ c(2.0, 0.01)    # very permissive
      )
      
      
      if(type == "onlyspatial"){
        form <- y ~ DEvar + offset(log(nCount_RNA)) +
          f(clusterloc, model = "bym2", graph = graph_path,
            scale.model = TRUE,constr = TRUE,
            hyper = list(
              prec = list(prior = "pc.prec", param = prior_prec_values),
              phi  = list(prior = "pc", param = prior_phi_values)
            ))
      }else if(type == "full"){
        form <- y ~ DEvar + otherct_expr + offset(log(nCount_RNA))+
          f(clusterloc, model = "bym2", graph = graph_path,
            scale.model = TRUE,constr = TRUE,
            hyper = list(
              prec = list(prior = "pc.prec", param = prior_prec_values),
              phi  = list(prior = "pc", param = prior_phi_values)
            ))
      }
      
      
      modBYM2 <- tryCatch({
        INLA::inla(form,
                   data = data_target,
                   control.predictor = list(compute = TRUE),
                   family = "nbinomial",
                   num.threads = "1:1" ,
                   control.compute = list(dic = TRUE, waic = TRUE),
                   quantiles = c(0.025/426,0.025/258, 0.025, 0.5, 1-0.025,1-0.025/258, 1-0.025/426),
                   control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
        )
      }, error = function(e) {
        return(NULL)
      })
      
    }else if(dist == "gaussian"){
      data_target$y <- norm[target,]
      prior_phi_values <- dplyr::case_when(
        prior_phi == "iid" ~ c(0.3, 0.7),       # Shrinks toward i.i.d. (often better for NegBin as it interacts with theta)
        prior_phi == "spatial" ~ c(0.7, 0.7),   # Shrinks toward ICAR
        prior_phi == "balanced" ~ c(0.5, 0.5)   # 50% chance it is greater than 0.5
      )
      prior_prec_values <- dplyr::case_when(
        prior_prec == "conservative" ~ c(0.5, 0.01),   # shrink to σ < 0.5 (relative to outcome scale)
        prior_prec == "moderate"     ~ c(1.0, 0.01),   # σ < 1.0
        prior_prec == "standard"     ~ c(1.5, 0.01),   # looser than NB
        prior_prec == "permissive"   ~ c(2.0, 0.01),   # quite loose
        prior_prec == "flat"         ~ c(3.0, 0.01)    # nearly flat
      )
      
      if(type == "onlyspatial"){
        form <- y ~ DEvar  +
          f(clusterloc, model = "bym2", graph = graph_path,
            scale.model = TRUE,
            hyper = list(
              prec = list(prior = "pc.prec", param = prior_prec_values),
              phi  = list(prior = "pc", param = prior_phi_values)
            ))
      }else if(type == "full"){
        form <- y ~ DEvar + otherct_expr +
          f(clusterloc, model = "bym2", graph = graph_path,
            scale.model = TRUE,
            hyper = list(
              prec = list(prior = "pc.prec", param = prior_prec_values),
              phi  = list(prior = "pc", param = prior_phi_values)
            ))
      }
      
      
      modBYM2 <- tryCatch({
        INLA::inla(form,
                   data = data_target,
                   control.predictor = list(compute = TRUE),
                   family = "gaussian",
                   num.threads = "1:1" ,
                   control.compute = list(dic = TRUE, waic = TRUE),
                   verbose = F,
                   quantiles = c(0.025/426,0.025/258, 0.025, 0.5, 1-0.025,1-0.025/258, 1-0.025/426)
        )
      }, error = function(e) {
        return(NULL)
      })
    }
    
    
    if (!is.null(modBYM2)) {
      coefs <- modBYM2$summary.fixed[,c(1:9)]
      coefs <- as.data.frame(cbind(var, target,dist, type, kprop,
                                   prior_phi, prior_prec, adjmethod, kadj,
                                   rownames(coefs), coefs,
                                   overdis = 1/modBYM2$summary.hyperpar[1,4],
                                   rangeest = modBYM2$summary.hyperpar[2,4],
                                   stdev = modBYM2$summary.hyperpar[3,4],
                                   failed = F
      ))
      colnames(coefs) <- c("DEvar","Gene","Distribution","Type","kprop","prior_phi", "prior_prec", "adjmethod", "kadj",
                           "Variable","Estimate","SE","LBbonf_all","LBbonf_cont",
                           "LB","Median","UB","UBbonf_cont","UBbonf_all","overdis","rangeest","stdev","failed")
      pseudo_p <- function(marg) {
        p_le_0 <- INLA::inla.pmarginal(0, marg)   # P(beta <= 0)
        2 * min(p_le_0, 1-p_le_0)
      }
      coefs$pseudopvalue <- pseudo_p(modBYM2$marginals.fixed[["DEvar"]] )
      return(coefs)
    }else{
      coefs <- data.frame(DEvar=var,
                          Gene=target,Distribution=dist,Type=type,kprop=kprop,
                          prior_phi=prior_phi,prior_prec=prior_prec,
                          Variable=NA,Estimate=NA,SE=NA,LBbonf_all=NA,LBbonf_cont=NA,
                          LB=NA,Median=NA,UB=NA,UBbonf_cont=NA,UBbonf_all=NA,overdis=NA,rangeest=NA,stdev=NA,failed = T,
                          pseudopvalue = NA
      )
    }
  }
  
  #----- Run Models ----
  if(nCores > 1){
    options(future.globals.maxSize = 5 * 1024^3)  # 5 GiB
    neighbor_expr_list <- as.matrix(neighbor_expr_list)
    if (future::supportsMulticore()) {
      future::plan(future::multicore, workers = ncores)
    } else {
      future::plan(future::multisession, workers = ncores)
    }
    results <- future_map(targets, function(target) {
      getbym2(var=var,
              target = target,
              data_target = metadata,
              counts=counts, norm = norm,
              neighbor_expr_list = neighbor_expr_list,
              mesh=mesh,
              dist = dist,
              type = type,
              kprop = kprop,
              prior_phi=prior_phi, prior_prec=prior_prec,
              adjmethod=adjmethod, kadj=kadj,graph_path = graph_path
      )
    }, .options = furrr_options(seed = TRUE, globals = TRUE,
                                packages = c("INLA", "stats", "spdep") ))
    plan("sequential")
  }else{
    results <- lapply(targets,function(target) {
      getbym2(var=var,
              target = target,
              data_target = metadata,
              counts=counts, norm = norm,
              neighbor_expr_list = neighbor_expr_list,
              mesh=mesh,
              dist = dist,
              type = type,
              kprop = kprop,
              prior_phi=prior_phi, prior_prec=prior_prec,
              adjmethod=adjmethod, kadj=kadj,graph_path = graph_path
      )}
    )
  }
  results <- do.call(rbind, results)
  results
}



# Get results -------------------------------------------------------------
# Only spatial
# Negative Binomial
if(methodrun == "spatial_nb"){
  print("Only spatial")
  print("Negative binomial")
  spatial_nb <- benchmark_and_return(quote(
    get_bym2_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = kprop,type = "onlyspatial", dist="negbin",
                     meshscenario=1,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_phi=prior_phi, prior_prec=prior_prec,
                     adjmethod=adjmethod, kadj=kadj,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  spatial_nb <- spatial_nb$result %>%
    mutate(Time =  spatial_nb$memory$Elapsed_Time_sec,
           TotalRAM= spatial_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  spatial_nb$memory$Peak_RAM_Used_MiB)
  write.csv(spatial_nb, file = paste0("ResultsDistance/BYM2/BYM2_spatial_nb",var,"_kprop",kprop,prior_phi,prior_prec,adjmethod,kadj,".csv"))
}

# Gaussian
if(methodrun == "spatial_gaus"){
  print("Gaussian")
  spatial_gaus <- benchmark_and_return(quote(
    get_bym2_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = kprop,type = "onlyspatial", dist="gaussian",
                     meshscenario=1,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_phi=prior_phi, prior_prec=prior_prec,
                     adjmethod=adjmethod, kadj=kadj,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  spatial_gaus <- spatial_gaus$result %>%
    mutate(Time =  spatial_gaus$memory$Elapsed_Time_sec,
           TotalRAM= spatial_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  spatial_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(spatial_gaus, file = paste0("ResultsDistance/BYM2/BYM2_spatial_gaus",var,"_kprop",kprop,prior_phi,prior_prec,adjmethod,kadj,".csv"))
}

# Full
# Negative Binomial
if(methodrun == "full_nb"){
  print("Negative binomial")
  full_nb <- benchmark_and_return(quote(
    get_bym2_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = kprop,type = "full", dist="negbin",
                     meshscenario=1,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_phi=prior_phi, prior_prec=prior_prec,
                     adjmethod=adjmethod, kadj=kadj,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  full_nb <- full_nb$result %>%
    mutate(Time =  full_nb$memory$Elapsed_Time_sec,
           TotalRAM= full_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  full_nb$memory$Peak_RAM_Used_MiB)
  write.csv(full_nb, file = paste0("ResultsDistance/BYM2/BYM2_full_nb",var,"_kprop",kprop,prior_phi,prior_prec,adjmethod,kadj,".csv"))
}
# Gaussian
if(methodrun == "full_gaus"){
  print("Gaussian")
  full_gaus <- benchmark_and_return(quote(
    get_bym2_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = kprop,type = "full", dist="gaussian",
                     meshscenario=1,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_phi=prior_phi, prior_prec=prior_prec,
                     adjmethod=adjmethod, kadj=kadj,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  full_gaus <- full_gaus$result %>%
    mutate(Time =  full_gaus$memory$Elapsed_Time_sec,
           TotalRAM= full_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  full_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(full_gaus, file = paste0("ResultsDistance/BYM2/BYM2_full_gaus",var,"_kprop",kprop,prior_phi,prior_prec,adjmethod,kadj,".csv"))
}