#--------------------------------#
#-- Fit SPDE to simulated data --#
#--------------------------------#
dir.out = "~/Simulation/"

library(tidyverse)
library(INLA)
library(fmesher)
library(peakRAM)

# Function to get results -------------------------------------------------------------
get_INLA_results <- function(gene, dist, meshscenario, prior_range, prior_sigma, data, Ycounts, Ygaus,scen){
  benchmark_and_return <- function(expr) {
    result <- NULL
    mem <- peakRAM(result <- eval(expr))
    list(result = result, memory = mem)
  }
  n <- nrow(data)
  p <- ncol(Ygaus)
  rownames(data) <- data$cell_ID
  
  #---- Setting parameters ----
  library(FNN)
  mnn <- {
    nn <- get.knn(as.matrix(cbind(data$sdimx, data$sdimy)), k = 1)$nn.dist[,1]
    median(nn, na.rm = TRUE)
  }
  
  # Mesh ---
  range <- max(max(data$sdimx) - min(data$sdimx), max(data$sdimy) - min(data$sdimy))
  cutoff <- range/sqrt(nrow(data))
  cutoff <- cutoff * meshscenario
  mesh <- fm_mesh_2d(loc = data[, c("sdimx","sdimy")],
                     offset = c(0.05 * range,  0.1 * range),
                     max.edge = c(cutoff * 5, cutoff * 15),
                     cutoff = cutoff
  )
  
  
  #----- Run Models ----
  if(dist == "negbin"){
    dataclust <- cbind(Ycounts, data)
    dataclust$geneinterest <- dataclust %>% pull(gene)
  }else if(dist == "gaussian"){
    dataclust <- cbind(Ycounts, data)
    dataclust$geneinterest <- dataclust %>% pull(gene)
  }
  
  # Set A
  A <- fm_basis(mesh, loc = data.matrix(dataclust[, c("sdimx","sdimy")]))
  stack <- inla.stack(
    tag = "est",
    data = list(geneinterest = dataclust$geneinterest),
    effects = list(
      s = 1:ncol(A),
      niche = dataclust$niche,
      totalcounts = dataclust$totalcounts
    ),
    A = list(A, 1, 1)
  )
  
  # Set priors
  prior_range_values <- case_when(
    prior_range == "short0.1"      ~ c( 0.1, 0.5),   
    prior_range == "standard1"     ~ c(1, 0.5),  
    prior_range == "long2"         ~ c(2, 0.5),   
    prior_range == "weak5"         ~ c(5, 0.5) 
  )
  
  if(dist == "negbin"){
    prior_sigma_values <- case_when(
      prior_sigma == "conservative" ~ c(0.5, 0.01),   # 99% chance that sigma is < 0.5 (prevents spatial field dominating the model, shrinks to low variance)
      prior_sigma == "moderate"     ~ c(0.75, 0.01),  # mild shrinkage
      prior_sigma == "standard"     ~ c(1, 0.01),     # 2.7-fold variation (standard)
      prior_sigma == "permissive"   ~ c(1.5, 0.01),   # allows moderate spatial variation
      prior_sigma == "weak"         ~ c(3, 0.01)      # almost not informative
    )
  }else{
    sd_gene <- sd(dataclust$geneinterest)
    prior_sigma_values <- case_when(
      prior_sigma == "conservative" ~ c(0.3 * sd_gene, 0.01),   # 99% chance that sigma is < 30% of total variation
      prior_sigma == "moderate" ~ c(0.5 * sd_gene, 0.01),       #
      prior_sigma == "standard" ~ c(sd_gene, 0.01),             # balanced
      prior_sigma == "permissive"   ~ c(1.5 * sd_gene, 0.01),   # allows large spatial SD
      prior_sigma == "weak" ~ c(3 * sd_gene, 0.01)              # almost not informative
    )
  }
  
  spde <- inla.spde2.pcmatern(mesh, alpha = 2,
                              prior.range = prior_range_values,
                              prior.sigma = prior_sigma_values,
                              constr = T)
  
  # Run
  formula <- geneinterest ~ offset(log(totalcounts)) + niche + f(s, model = spde)
  
  if(dist == "negbin"){
    modSPDE <- benchmark_and_return(quote(
      inla(formula,
           data = inla.stack.data(stack, spde = spde),
           control.predictor = list(A = inla.stack.A(stack), compute = T),
           family = "nbinomial",
           control.inla = list(int.strategy = "eb"),
           verbose = F,
           control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
      )
    ))
  }else{
    modSPDE <- benchmark_and_return(quote(
      inla(formula,
           data = inla.stack.data(stack, spde = spde),
           control.predictor = list(A = inla.stack.A(stack), compute = T),
           family = "gaussian",
           control.inla = list(int.strategy = "eb"),
           verbose = F
      )
    ))
  }
  
  coefs <- modSPDE$result$summary.fixed[,c(4,3,5)]
  print(coefs)
  coefs <- as.data.frame(cbind(scen,gene,  "SPDE",dist,meshscenario,
                               prior_range, prior_sigma,
                               rownames(coefs), coefs[,1],NA,coefs[,2:3],
                               1/modSPDE$result$summary.hyperpar[1,4],
                               modSPDE$result$summary.hyperpar[2,4],
                               modSPDE$result$summary.hyperpar[3,4],
                               modSPDE$memory$Elapsed_Time_sec,
                               modSPDE$memory$Total_RAM_Used_MiB,
                               modSPDE$memory$Peak_RAM_Used_MiB
  ))
  colnames(coefs) <- c("Scenario","Gene","Model","Distribution","mesh","prior_range","prior_sigma",
                       "Variable","Estimate","SE","LB","UB","overdis","rangeest","stdev","Time","TotalRAM","PeakRAM")
  
  pseudo_p <- function(marg) {
    p_le_0 <- INLA::inla.pmarginal(0, marg)   # P(beta <= 0)
    p_ge_0 <- 1 - p_le_0                      # P(beta > 0)
    2 * min(p_le_0, p_ge_0)
  }
  coefs$pseudopvalue <- pseudo_p(modSPDE$result$marginals.fixed[["niche"]] )
  
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
                         prior_range = c("short0.1","standard1","long2","weak5"),
                         prior_sigma = c("conservative", "moderate","standard", "permissive", "weak")
) 

results<-do.call(rbind, lapply(c(0.5,2), function(meshscenario){
  with(simparams[args,],
       get_INLA_results(gene=gene, dist=dist, meshscenario = meshscenario,
                        prior_range=prior_range,prior_sigma=prior_sigma, data=data,
                        Ycounts=Ycounts, Ygaus=Ygaus,scen=scen))
}))
write.csv(results, file = paste0(dir.out,paste("Results/PartialSimulation/INLA/INLA",
                                 paste(sapply(simparams[args, ], as.character), collapse = "_") ,".csv",sep="_")))



