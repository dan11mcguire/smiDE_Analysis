#--------------------------------#
#-- Fit BYM2 to simulated data --#
#--------------------------------#
dir.out = "~/Simulation/"
library(tidyverse)
library(INLA)
library(fmesher)
library(spdep)
library(peakRAM, lib.loc)


# Function to get results -------------------------------------------------------------
get_BYM2_results <- function(gene, kprop, dist, adjmethod, kadj,prior_phi, prior_prec,
                             data, Ycounts, Ygaus,scen){
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
  
  #---- Setting parameters ----
  # Create adjacency graph between centroids
  # Delaunay triangulation
  if(adjmethod == "delauneytrig"){
    nb <- data %>% dplyr::select(cluster, sdimxclust, sdimyclust) %>%
      distinct %>% arrange(cluster) %>% dplyr::select(sdimxclust,sdimyclust) %>% tri2nb
  }else if(adjmethod == "knn"){
    kn <- knearneigh(data[, c("sdimxclust", "sdimyclust")], k = kadj)
    nb <- knn2nb(kn, row.names = dataclust$cluster)
  }
  adjmat <- nb2mat(nb, style = "B", zero.policy = TRUE)
  library(Matrix)
  graph <- Matrix(adjmat, sparse = TRUE)
  
  # distance based
  # nb <- dnearneigh(coords, d1 = 0, d2 = 100, row.names = data$cluster)
  
  #----- Run Models ----
  if(dist == "negbin"){
    dataclust <- cbind(Ycounts, data)
    dataclust$geneinterest <- dataclust %>% pull(gene)
  }else if(dist == "gaussian"){
    dataclust <- cbind(Ycounts, data)
    dataclust$geneinterest <- dataclust %>% pull(gene)
  }
  
  # Set priors
  # mixing between spatial and iif random effects, phi = 0 fully iid, phi = 1 fully spatial (ICAR)
  prior_phi_values <- case_when(
    prior_phi == "iid" ~ c(0.3, 0.5),       # Shrinks toward i.i.d.
    prior_phi == "spatial" ~ c(0.7, 0.5),   # Shrinks toward ICAR
    prior_phi == "balanced" ~ c(0.5, 0.5)   # 50% chance it is greater than 0.5
  )
  
  if(dist == "negbin"){
    # total precision of the spatial effects
    prior_prec_values <- case_when(
      prior_prec == "conservative" ~ c(0.3, 0.01),   # shrink to σ < 0.3 → strong control
      prior_prec == "moderate"     ~ c(0.5, 0.01),   # shrink to σ < 0.5 → safe for NB
      prior_prec == "standard"     ~ c(0.75, 0.01),  # balanced, allows spatial variation
      prior_prec == "permissive"   ~ c(1.0, 0.01),   # wider variation
      prior_prec == "flat"         ~ c(2.0, 0.01)    # very permissive
    )
  }else{
    prior_prec_values <- case_when(
      prior_prec == "conservative" ~ c(0.5, 0.01),   # shrink to σ < 0.5 (relative to outcome scale)
      prior_prec == "moderate"     ~ c(1.0, 0.01),   # σ < 1.0
      prior_prec == "standard"     ~ c(1.5, 0.01),   # looser than NB
      prior_prec == "permissive"   ~ c(2.0, 0.01),   # quite loose
      prior_prec == "flat"         ~ c(3.0, 0.01)    # nearly flat
    )
  }
  
  
  # Run
  formula <- geneinterest ~ niche +
    f(cluster, model = "bym2", graph = graph,
      scale.model = TRUE,
      hyper = list(
        prec = list(prior = "pc.prec", param = prior_prec_values),
        phi  = list(prior = "pc", param = prior_phi_values)
      ))
  if(dist == "negbin"){
    modBYM2 <- benchmark_and_return(quote(
      inla(formula,
           E = totalcounts,
           data = dataclust,
           control.predictor = list(compute = TRUE),
           family = "nbinomial",
           control.compute = list(dic = TRUE, waic = TRUE),
           verbose = F,
           control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
      )
    ))
  }else{
    modBYM2 <- benchmark_and_return(quote(
      inla(formula,
           data = dataclust,
           control.predictor = list(compute = TRUE),
           family = "gaussian",
           control.compute = list(dic = TRUE, waic = TRUE),
           verbose = F
      )
    ))
  }
  
  coefs <- modBYM2$result$summary.fixed[,c(4,3,5)]
  print(coefs)
  coefs <- as.data.frame(cbind(scen, gene, kprop, "BYM2",dist,adjmethod, kadj,
                               prior_phi, prior_prec,
                               rownames(coefs), coefs[,1],NA,coefs[,2:3],
                               1/modBYM2$result$summary.hyperpar[1,4],
                               modBYM2$result$summary.hyperpar[2,4],
                               modBYM2$result$summary.hyperpar[3,4],
                               modBYM2$memory$Elapsed_Time_sec,
                               modBYM2$memory$Total_RAM_Used_MiB,
                               modBYM2$memory$Peak_RAM_Used_MiB
  ))
  colnames(coefs) <- c("Scenario","Gene","kprop","Model","Distribution","AdjacencyMethod","Kadj","prior_phi","prior_prec",
                       "Variable","Estimate","SE","LB","UB","overdis","rangeest","stdev","Time","TotalRAM","PeakRAM")
  pseudo_p <- function(marg) {
    p_le_0 <- INLA::inla.pmarginal(0, marg)
    p_ge_0 <- 1 - p_le_0  
    2 * min(p_le_0, p_ge_0)
  }
  coefs$pseudopvalue <- pseudo_p(modBYM2$result$marginals.fixed[["niche"]] )
  
  coefs
}

# Read data ---------------------------------------------------------------
load(paste0(dir.out,paste0("SimulatedY.Rdata")))
n <- nrow(data)
data$cell_ID <- as.character(data$cell_ID)
data$niche <- data$nichenum

# Run ---------------------------------------------------------------------
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)
simparams <- expand_grid(gene = colnames(Ycounts),
                         dist = c("negbin","gaussian"),
                         prior_phi = c("iid","spatial","balanced"),
                         prior_prec = c("conservative","moderate", "standard", "permissive", "flat"),
                         kprop = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                                   0.08, 0.09,  0.1, 0.15,  0.2,0.3,0.4,0.5,1)
) 

for(adjmethod in c("delauneytrig")){
  results<- with(simparams[args,],
         get_BYM2_results(gene=gene, kprop=kprop, dist=dist, adjmethod = adjmethod,
                          kadj = ifelse(adjmethod == "delauneytrig",NA,4),
                          prior_phi=prior_phi,prior_prec=prior_prec,  data=data, Ycounts=Ycounts, Ygaus=Ygaus,scen=scen))
  
write.csv(results, file = paste0(dir.out,paste("Results/PartialSimulation/BYM2/BYM2",adjmethod,
                                                   paste(sapply(simparams[args, ], as.character), collapse = "_") ,".csv",sep="_")))
}


