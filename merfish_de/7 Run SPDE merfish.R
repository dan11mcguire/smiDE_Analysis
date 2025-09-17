#-------------------#
#- Run SPDE Models -#
#-------------------#
data_dir <- "~MERFISH/"
ncores<- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')) -1

library(tidyverse)
library(peakRAM)
library(future)
library(future.apply)
library(furrr)
library(INLA)
library(fmesher)
library(broom)
benchmark_and_return <- function(expr) {
  result <- NULL
  mem <- peakRAM(result <- eval(expr))
  list(result = result, memory = mem)
}
RankNorm <-
  function (u, k = 0.375)
  {
    n <- length(u)
    r <- rank(u)
    out <- qnorm((r - k)/(n - 2 * k + 1))
    return(out)
  }

slurm_arrayid <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
params <- expand_grid(meshscenario = c(3,0.5),
                      prior_range = c("short","standard","long","weak"),
                      prior_sigma = c("conservative","moderate","standard","permissive","flat"),
                      methodrun = c("full_nb","full_gaus"),
                      var = c("kernelMicroglia1000")
)

methodrun <- params$methodrun[slurm_arrayid]
meshscenario <- params$meshscenario[slurm_arrayid]
prior_range <- params$prior_range[slurm_arrayid]
prior_sigma <- params$prior_sigma[slurm_arrayid]
var <- params$var[slurm_arrayid]


load(file = paste0(data_dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
load(paste0(data_dir, "Allen/WholeBrain_df_distance.Rdata"))
filteredgenes <- overlap_metrics_allen$result %>% filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1, ratio < 1) %>%pull(target)

# Set functions -----------------------------------------------------------

get_INLA_results <- function(var, metadata, counts, norm,kprop,type, dist, meshscenario,neighbor_expr_list,
                             prior_range, prior_sigma,targets=NULL,nCores){
  metadata$DEvar <- metadata[,var]
  n <- nrow(metadata)
  metadata <- metadata %>% dplyr::mutate(sdimx = x,sdimy = y)
  rownames(metadata) <- metadata$cell_ID
  
  #---- Setting parameters ----
  library(FNN)
  mnn <- {
    nn <- get.knn(as.matrix(cbind(metadata$x, metadata$y)), k = 1)$nn.dist[,1]
    median(nn, na.rm = TRUE)
  }
  # Mesh ---
  range <- max(max(metadata$x) - min(metadata$x),
               max(metadata$y) - min(metadata$y))
  cutoff <- range/sqrt(n)
  cutoff <- cutoff * meshscenario
  mesh <- fm_mesh_2d(loc = metadata[, c("sdimx","sdimy")],
                     offset = c(0.05 * range,  0.1 * range),
                     cutoff = cutoff
  )
  if(is.null(targets))targets <- rownames(counts)
  
  #----- Define function ----
  getinla <- function(var,target, data_target,counts,norm, neighbor_expr_list,
                      mesh,meshscenario, dist, type, kprop, prior_range, prior_sigma){
    INLA::inla.setOption(num.threads = "1:1")
    message("Threads in worker: ", INLA::inla.getOption("num.threads"))
    cat(paste0("Fitting model to target ",target, "."))
    data_target$otherct_expr <- RankNorm(neighbor_expr_list[target,data_target$cell_ID])
    
    if(dist != "gaussian"){
      data_target$y <- counts[target,]
      
      # Set A
      A <- fm_basis(mesh, loc = data.matrix(data_target[, c("sdimx","sdimy")]))
      stack <- inla.stack(
        tag = "est",
        data = list(y = data_target$y),
        effects = list(
          s = 1:ncol(A),
          DEvar = data_target$DEvar,
          nCount_RNA = data_target$nCount_RNA,
          otherct_expr = data_target$otherct_expr
        ),
        A = list(A, 1,1,1)
      )
      
      prior_range_values <- case_when(
        prior_range == "short"    ~ c( 5 * mnn, 0.5),   # P(ρ < 5·mnn)  = 0.5  (more local)
        prior_range == "standard" ~ c(10 * mnn, 0.5),   # P(ρ < 10·mnn) = 0.5  (default)
        prior_range == "long"     ~ c(20 * mnn, 0.5),   # P(ρ < 20·mnn) = 0.5  (smoother)
        prior_range == "weak"     ~ c(range/3, 0.5)     # very broad, tied to slide span
      )
      prior_sigma_values <- case_when(
        prior_sigma == "conservative" ~ c(0.5, 0.01),   # 99% chance that sigma is < 0.5 (prevents spatial field dominating the model, shrinks to low variance)
        prior_sigma == "moderate"     ~ c(0.75, 0.01),  # mild shrinkage
        prior_sigma == "standard" ~ c(1, 0.01),         # 2.7-fold variation (standard)
        prior_sigma == "permissive" ~ c(1.5, 0.01),     # allows moderate spatial variation
        prior_sigma == "flat" ~ c(3, 0.01)              # almost not informative
      )
      spde <- inla.spde2.pcmatern(mesh, alpha = 2,
                                  prior.range = prior_range_values,
                                  prior.sigma = prior_sigma_values,
                                  constr = T)
      
      if(type == "onlyspatial"){
        form <- y ~ DEvar + offset(log(nCount_RNA)) + f(s, model = spde)
        formnull <-  y ~ offset(log(nCount_RNA)) + f(s, model = spde)
      }else if(type == "full"){
        form <- y ~ DEvar + otherct_expr + offset(log(nCount_RNA))+ f(s, model = spde)
        formnull <- y ~ otherct_expr + offset(log(nCount_RNA))+ f(s, model = spde)
      }
      
      modSPDE <- tryCatch({
          inla(form,
               data = inla.stack.data(stack, spde = spde),
               control.predictor = list(A = inla.stack.A(stack), compute = T),
               family = "nbinomial",
               control.inla = list(strategy = "gaussian"),
               verbose = F,
               num.threads = "1:1" ,
               quantiles = c(0.025/426,0.025/258, 0.025, 0.5, 1-0.025,1-0.025/258, 1-0.025/426),
               control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
          )
        }, error = function(e) {
          return(NULL)
        })
    }else if(dist == "gaussian"){
      data_target$y <- norm[target,]
      
      # Set A
      A <- fm_basis(mesh, loc = data.matrix(data_target[, c("sdimx","sdimy")]))
      stack <- inla.stack(
        tag = "est",
        data = list(y = data_target$y),
        effects = list(
          s = 1:ncol(A),
          DEvar = data_target$DEvar,
          nCount_RNA = data_target$nCount_RNA,
          otherct_expr = data_target$otherct_expr
        ),
        A = list(A, 1,1,1)
      )
      range <- max(max(data_target$sdimx) - min(data_target$sdimx),
                   max(data_target$sdimy) - min(data_target$sdimy))
      prior_range_values <- case_when(
        prior_range == "short"    ~ c( 5 * mnn, 0.5),   # P(ρ < 5·mnn)  = 0.5  (more local)
        prior_range == "standard" ~ c(10 * mnn, 0.5),   # P(ρ < 10·mnn) = 0.5  (default)
        prior_range == "long"     ~ c(20 * mnn, 0.5),   # P(ρ < 20·mnn) = 0.5  (smoother)
        prior_range == "weak"     ~ c(range/3, 0.5)     # very broad, tied to slide span
      )
      sd_gene <- sd(data_target$y)
      prior_sigma_values <- case_when(
        prior_sigma == "conservative" ~ c(0.3 * sd_gene, 0.01),   # 99% chance that sigma is < 30% of total variation
        prior_sigma == "moderate" ~ c(0.5 * sd_gene, 0.01),       #
        prior_sigma == "standard" ~ c(sd_gene, 0.01),             # balanced
        prior_sigma == "permissive"   ~ c(1.5 * sd_gene, 0.01),   # allows large spatial SD
        prior_sigma == "flat" ~ c(3 * sd_gene, 0.01)              # almost not informative
      )
      spde <- inla.spde2.pcmatern(mesh, alpha = 2,
                                  prior.range = prior_range_values,
                                  prior.sigma = prior_sigma_values,
                                  constr = T)
      
      if(type == "onlyspatial"){
        form <- y ~ DEvar  + f(s, model = spde)
        formnull <- y ~  f(s, model = spde)
      }else if(type == "full"){
        form <- y ~ DEvar + otherct_expr + f(s, model = spde)
        formnull <- y ~ otherct_expr + f(s, model = spde)
      }
      
      
      modSPDE <- tryCatch({
        inla(form,
             data = inla.stack.data(stack, spde = spde),
             control.predictor = list(A = inla.stack.A(stack), compute = T),
             family = "gaussian",
             control.inla = list(strategy = "gaussian"),
             verbose = F,
             num.threads = "1:1" ,
             quantiles = c(0.025/426,0.025/258, 0.025, 0.5, 1-0.025,1-0.025/258, 1-0.025/426)
        )
      }, error = function(e) {
        return(NULL)
      })
    }
    
    if (!is.null(modSPDE)) {
      print("modSPDE ")
      coefs <- modSPDE$summary.fixed[,c(1:9)]
      coefs <- as.data.frame(cbind(var,target,dist, type, kprop,meshscenario,
                                   prior_range, prior_sigma,
                                   rownames(coefs), coefs,
                                   1/modSPDE$summary.hyperpar[1,4],
                                   modSPDE$summary.hyperpar[2,4],
                                   modSPDE$summary.hyperpar[3,4],
                                   FALSE
      ))
      colnames(coefs)<- c("DEvar","Gene","Distribution","Type","kprop","meshscenario","prior_range","prior_sigma",
                          "Variable","Estimate","SE","LBbonf_all","LBbonf_cont",
                          "LB","Median","UB","UBbonf_cont","UBbonf_all","overdis","rangeest","stdev","failed")
      pseudo_p <- function(marg) {
        p_le_0 <- INLA::inla.pmarginal(0, marg)   
        p_ge_0 <- 1 - p_le_0                     
        2 * min(p_le_0, p_ge_0)
      }
      coefs$pseudopvalue <- pseudo_p(modSPDE$marginals.fixed[["DEvar"]] )
      return(coefs)
    }else{
      coefs <- data.frame(
        DEvar=var,Gene=target,Distribution=dist,Type=type,kprop=kprop,meshscenario=meshscenario,
        prior_range=prior_range,prior_sigma=prior_sigma,
        Variable=NA,Estimate=NA,SE=NA,LBbonf_all=NA,LBbonf_cont=NA,
        LB=NA,Median=NA,UB=NA,UBbonf_cont=NA,UBbonf_all=NA,overdis=NA,rangeest=NA,stdev=NA,
        failed=T,pseudopvalue = NA
      )
    }
  }
  
  #----- Run Models ----
  if(nCores > 1){
    neighbor_expr_list <- as.matrix(neighbor_expr_list)
    plan(multisession, workers = nCores)
    results <- future_map(targets, function(target) {
      getinla(var=var,
              target = target,
              data_target = metadata,
              counts = counts, norm = norm,
              neighbor_expr_list = neighbor_expr_list,
              mesh=mesh,meshscenario=meshscenario,
              dist = dist,
              type = type,
              kprop = kprop,
              prior_range=prior_range,
              prior_sigma=prior_sigma
      )
    }, .options = furrr_options(seed = TRUE))
    plan("sequential")
  }else{
    results <- lapply(targets,function(target) {
      getinla(var=var,
              target = target,
              data_target = metadata,
              counts = counts, norm = norm,
              neighbor_expr_list = neighbor_expr_list,
              mesh=mesh,meshscenario=meshscenario,
              dist = dist,
              type = type,
              kprop = kprop,
              prior_range=prior_range,
              prior_sigma=prior_sigma
      )})
  }
  results <- do.call(rbind, results)
  results
}



# Get results -------------------------------------------------------------
#for(kprop in c(0.05, 0.1,1)){
print(paste("meshscenario",meshscenario))
# Only spatial
print("Only spatial")
# Negative Binomials
if(methodrun == "spatial_nb"){
  print("Negative binomial")
  spatial_nb <- benchmark_and_return(quote(
    get_INLA_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = 1,
                     type = "onlyspatial", dist="negbin",
                     meshscenario=meshscenario,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_range = prior_range,
                     prior_sigma = prior_sigma,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  spatial_nb <- spatial_nb$result %>%
    mutate(Time =  spatial_nb$memory$Elapsed_Time_sec,
           TotalRAM= spatial_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  spatial_nb$memory$Peak_RAM_Used_MiB)
  write.csv(spatial_nb, file = paste0("ResultsDistance/INLA/INLA_spatial_nb",var,"_mesh",meshscenario,prior_range,prior_sigma,".csv"))
}
# Gaussian
if(methodrun == "spatial_gaus"){
  print("Gaussian")
  spatial_gaus <- benchmark_and_return(quote(
    get_INLA_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = 1,type = "onlyspatial", dist="gaussian",
                     meshscenario=meshscenario,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_range = prior_range,
                     prior_sigma = prior_sigma,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  spatial_gaus <- spatial_gaus$result %>%
    mutate(Time =  spatial_gaus$memory$Elapsed_Time_sec,
           TotalRAM= spatial_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  spatial_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(spatial_gaus, file = paste0("ResultsDistance/INLA/INLA_spatial_gaus",var,"_mesh",meshscenario,prior_range,prior_sigma,".csv"))
}

# Full
print("Full")
# Negative Binomial
if(methodrun == "full_nb"){
  print("Negative binomial")
  full_nb <- benchmark_and_return(quote(
    get_INLA_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = 1,type = "full", dist="negbin",
                     meshscenario=meshscenario,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_range = prior_range,
                     prior_sigma = prior_sigma,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  full_nb <- full_nb$result %>%
    mutate(Time =  full_nb$memory$Elapsed_Time_sec,
           TotalRAM= full_nb$memory$Total_RAM_Used_MiB,
           PeakRAM =  full_nb$memory$Peak_RAM_Used_MiB)
  write.csv(full_nb, file = paste0("ResultsDistance/INLA/INLA_full_nb",var,"_mesh",meshscenario,prior_range,prior_sigma,".csv"))
}
# Gaussian
if(methodrun == "full_gaus"){
  print("Gaussian")
  full_gaus <- benchmark_and_return(quote(
    get_INLA_results(var=var,metadata=dataDEmetadata, counts=dataDEcounts, norm=dataDEnorm,kprop = 1,type = "full", dist="gaussian",
                     meshscenario=meshscenario,
                     neighbor_expr_list=pre_de_allen$result$nblist$neighbor_expr_byct$otherct,
                     prior_range = prior_range,
                     prior_sigma = prior_sigma,
                     targets = filteredgenes,
                     nCores = ncores)
  ))
  full_gaus <- full_gaus$result %>%
    mutate(Time =  full_gaus$memory$Elapsed_Time_sec,
           TotalRAM= full_gaus$memory$Total_RAM_Used_MiB,
           PeakRAM =  full_gaus$memory$Peak_RAM_Used_MiB)
  write.csv(full_gaus, file = paste0("ResultsDistance/INLA/INLA_full_gaus",var,"_mesh",meshscenario,prior_range,prior_sigma,".csv"))
}


