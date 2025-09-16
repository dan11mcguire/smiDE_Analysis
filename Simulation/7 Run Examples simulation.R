#------------------------------#
#-- Getting examples --#
#------------------------------#
dir.out = "/home/avasconc/RA-Ali/DifferentialExpression/Simulation/"

library(tidyverse)
library(MASS)
library(spaMM)

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)
simparams <- expand.grid(sim = c(1,10, 99))
sim <-  simparams$sim[args]

# Read in data
# load(paste0(dir.out,"Lung5-3data.Rdata"))
# data <-  fulldata <- fulldata %>% mutate(auxID = 1:nrow(fulldata)) %>%
#   filter(cell_type == "macrophage",
#          niche %in% c("myeloid-enriched stroma","stroma"))
# data$niche <- ifelse(data$niche == "myeloid-enriched stroma",1,0)
# n <- nrow(data)
# kprops <- c(0.05, 0.25, 0.50)
# dataclust <- data
# for (kp in kprops) {
#   k <- round(nrow(data) * kp)
#
#   info_clust <- kmeans(data[, c("sdimx", "sdimy")], k, iter.max = 50)
#
#   # cluster assignment
#   cluster_assignments <- data.frame(
#     cell_ID = data$cell_ID,
#     cluster = as.integer(as.factor(info_clust$cluster))
#   )
#
#   # centers
#   centers_df <- data.frame(
#     cluster = (as.factor(1:k)),
#     sdimx = info_clust$centers[,1],
#     sdimy = info_clust$centers[,2]
#   )
#
#   # join cluster + centers
#   clust_info <- left_join(cluster_assignments, centers_df, by = "cluster")
#
#   # rename columns with suffix (5, 25, 50, etc.)
#   suffix <- as.character(kp*100)
#   names(clust_info)[names(clust_info) == "cluster"] <- paste0("cluster", suffix)
#   names(clust_info)[names(clust_info) == "sdimx"]   <- paste0("sdimx", suffix)
#   names(clust_info)[names(clust_info) == "sdimy"]   <- paste0("sdimy", suffix)
#
#   # merge back into dataclust
#   dataclust <- left_join(dataclust, clust_info, by = "cell_ID")
# }
# rownames(dataclust) <- dataclust$cell_ID
# save(dataclust, file=paste0(dir.out,"dataclust.Rdata"))
load(paste0(dir.out,"dataclust.Rdata"))
load(paste0(dir.out,"intensities.Rdata"))
n <- nrow(dataclust)

#-- Scenario --
if(sim == 99){
  set.seed(67983*sim)
  lambda1 <- dataclust$totalcounts * 0.01 * exp(rnorm(n,0,0.2))
  lambda2 <-  5 * as.numeric(((dataclust$sdimy +17 )^2 < 0.5 & (dataclust$sdimx - 17 )^2 < 0.5 ))
  set.seed(7695*sim)
  dataclust <- dataclust %>% dplyr::mutate(Ycounts = rpois(n,lambda1) + rpois(n,lambda2))
  
}else{
  load(paste0(dir.out,paste0("SimulatedY.Rdata")))
  dataclust$Ycounts <- Ycounts[,sim]
}

ggplot(dataclust, aes(x = sdimx, y = sdimy, col = Ycounts)) +
  geom_point() + scale_color_viridis_c()

dataclust <- dplyr::ungroup(dataclust)
# Naive -------------------------------------------------------------------
## No spatial ------
mod_ind <- MASS::glm.nb(Ycounts ~ niche,data = dataclust)

# Cluster based ------------------------------------------------------------------
dataclust$cluster5 <- as.factor(dataclust$cluster5)
dataclust$cluster25 <- as.factor(dataclust$cluster25)
dataclust$cluster50 <- as.factor(dataclust$cluster50)


#----Indep clusters --
mod_spa_ind5 <- fitme(Ycounts ~ niche + (1 | cluster5 ),data = dataclust, family = negbin())
mod_spa_ind25 <- fitme(Ycounts ~ niche + (1 | cluster25 ),data = dataclust, family = negbin())
mod_spa_ind50 <- fitme(Ycounts ~ niche + (1 | cluster50 ),data = dataclust, family = negbin())

#---- GP --
mod_spa_gp5 <- fitme(Ycounts ~ niche + Matern(1 | sdimx5 + sdimy5),
                     data = dataclust, family = negbin(),verbose=c(TRACE=F),
                     fixed=list(nu=0.5))
mod_spa_gp25 <- fitme(Ycounts ~ niche + Matern(1 | sdimx25 + sdimy25),
                      data = dataclust, family = negbin(),verbose=c(TRACE=F),
                      fixed=list(nu=0.5))
mod_spa_gp50 <- fitme(Ycounts ~ niche + Matern(1 | sdimx50 + sdimy50),
                      data = dataclust, family = negbin(),verbose=c(TRACE=F),
                      fixed=list(nu=0.5))

# BYM2 --------------------------------------------------------------------
library(Matrix)

# Set priors
prior_phi_values <- c(0.5, 0.5)
prior_prec_values <- c(0.75, 0.01)

# 5%
nb5 <- dataclust[!duplicated(dataclust$cluster5), c("sdimx5", "sdimy5")]
nb5 <- nb5[order(dataclust$cluster5[!duplicated(dataclust$cluster5)]), ]
nb5 <- tri2nb(nb5)
adjmat5 <- nb2mat(nb5, style = "B", zero.policy = TRUE)
graph5 <- Matrix(adjmat5, sparse = TRUE)

formula <- Ycounts ~ niche +
  f(cluster5, model = "bym2", graph = graph5,
    scale.model = TRUE,
    hyper = list(
      prec = list(prior = "pc.prec", param = prior_prec_values),
      phi  = list(prior = "pc", param = prior_phi_values)
    ))
modBYM2_5 <- inla(formula,
                  data = dataclust,
                  control.predictor = list(compute = TRUE),
                  family = "nbinomial",
                  control.compute = list(dic = TRUE, waic = TRUE),
                  verbose = F,
                  control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
)

# 25%
nb25 <- dataclust[!duplicated(dataclust$cluster25),c("sdimx25", "sdimy25")]
nb25 <- nb25[order(dataclust$cluster25[!duplicated(dataclust$cluster25)]), ]
nb25 <- tri2nb(nb25)
adjmat25 <- nb2mat(nb25, style = "B", zero.policy = TRUE)
graph25 <- Matrix(adjmat25, sparse = TRUE)

formula <- Ycounts ~ niche +
  f(cluster25, model = "bym2", graph = graph25,
    scale.model = TRUE,
    hyper = list(
      prec = list(prior = "pc.prec", param = prior_prec_values),
      phi  = list(prior = "pc", param = prior_phi_values)
    ))
modBYM2_25 <- inla(formula,
                   data = dataclust,
                   control.predictor = list(compute = TRUE),
                   family = "nbinomial",
                   control.compute = list(dic = TRUE, waic = TRUE),
                   verbose = F,
                   control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
)

# 50%
nb50 <- dataclust[!duplicated(dataclust$cluster50),c("sdimx50", "sdimy50")]
nb50 <- nb50[order(dataclust$cluster50[!duplicated(dataclust$cluster50)]), ]
nb50 <- tri2nb(nb50)
adjmat50 <- nb2mat(nb50, style = "B", zero.policy = TRUE)
graph50 <- Matrix(adjmat50, sparse = TRUE)

formula <- Ycounts ~ niche +
  f(cluster50, model = "bym2", graph = graph50,
    scale.model = TRUE,
    hyper = list(
      prec = list(prior = "pc.prec", param = prior_prec_values),
      phi  = list(prior = "pc", param = prior_phi_values)
    ))
modBYM2_50 <- inla(formula,
                   data = dataclust,
                   control.predictor = list(compute = TRUE),
                   family = "nbinomial",
                   control.compute = list(dic = TRUE, waic = TRUE),
                   verbose = F,
                   control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
)

# INLA -------------------------------------------------------------------
library(FNN)
mnn <- {
  nn <- get.knn(as.matrix(cbind(dataclust$sdimx, dataclust$sdimy)), k = 1)$nn.dist[,1]
  median(nn, na.rm = TRUE)
}

# Dense mesh ---
range <- max(max(dataclust$sdimx) - min(dataclust$sdimx), max(dataclust$sdimy) - min(dataclust$sdimy))
cutoff <- range/sqrt(nrow(dataclust))
mesh <- fm_mesh_2d(loc = dataclust[, c("sdimx","sdimy")],
                   offset = c(0.05 * range,  0.1 * range),
                   cutoff =  0.05 * cutoff
)

# Set A
A <- fm_basis(mesh, loc = data.matrix(dataclust[, c("sdimx","sdimy")]))
stack <- inla.stack(
  tag = "est",
  data = list(Ycounts = dataclust$Ycounts),
  effects = list(
    s = 1:ncol(A),
    niche = dataclust$niche
  ),
  A = list(A, 1)
)

# Set priors
prior_range_values <- c(10 * mnn, 0.5)
prior_sigma_values <- c(1, 0.01)

spde <- inla.spde2.pcmatern(mesh, alpha = 2,
                            prior.range = prior_range_values,
                            prior.sigma = prior_sigma_values,
                            constr = T)

# Run
formula <- Ycounts ~ niche + f(s, model = spde)
modSPDEdense <- inla(formula,
                     data = inla.stack.data(stack, spde = spde),
                     control.predictor = list(A = inla.stack.A(stack), compute = T),
                     family = "nbinomial",
                     control.inla = list(int.strategy = "eb"),
                     verbose = F,
                     control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
)

# Dense mesh ---
range <- max(max(dataclust$sdimx) - min(dataclust$sdimx), max(dataclust$sdimy) - min(dataclust$sdimy))
cutoff <- range/sqrt(nrow(dataclust))
meshdense <- fm_mesh_2d(loc = dataclust[, c("sdimx","sdimy")],
                        offset = c(0.05 * range,  0.1 * range),
                        cutoff =  0.05 * cutoff
)

# Set A
A <- fm_basis(meshdense, loc = data.matrix(dataclust[, c("sdimx","sdimy")]))
stack <- inla.stack(
  tag = "est",
  data = list(Ycounts = dataclust$Ycounts),
  effects = list(
    s = 1:ncol(A),
    niche = dataclust$niche
  ),
  A = list(A, 1)
)

# Set priors
spde <- inla.spde2.pcmatern(meshdense, alpha = 2,
                            prior.range = prior_range_values,
                            prior.sigma = prior_sigma_values,
                            constr = T)

# Run
formula <- Ycounts ~ niche + f(s, model = spde)
modSPDEdense <- inla(formula,
                     data = inla.stack.data(stack, spde = spde),
                     control.predictor = list(A = inla.stack.A(stack), compute = T),
                     family = "nbinomial",
                     control.inla = list(int.strategy = "eb"),
                     verbose = F,
                     control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01))))
)

# Coarse mesh ---
range <- max(max(dataclust$sdimx) - min(dataclust$sdimx), max(dataclust$sdimy) - min(dataclust$sdimy))
cutoff <- range/sqrt(nrow(dataclust))
meshcoarse <- fm_mesh_2d(loc = dataclust[, c("sdimx","sdimy")],
                         offset = c(0.05 * range,  0.1 * range),
                         cutoff =  3 * cutoff
)

# Set A
A <- fm_basis(meshcoarse, loc = data.matrix(dataclust[, c("sdimx","sdimy")]))
stack <- inla.stack(
  tag = "est",
  data = list(Ycounts = dataclust$Ycounts),
  effects = list(
    s = 1:ncol(A),
    niche = dataclust$niche
  ),
  A = list(A, 1)
)

# Set priors
spde <- inla.spde2.pcmatern(meshcoarse, alpha = 2,
                            prior.range = prior_range_values,
                            prior.sigma = prior_sigma_values,
                            constr = T)

# Run
formula <- Ycounts ~ niche + f(s, model = spde)
modSPDEcoarse <- inla(formula,
                      data = inla.stack.data(stack, spde = spde),
                      control.predictor = list(A = inla.stack.A(stack), compute = T),
                      family = "nbinomial",
                      control.inla = list(int.strategy = "eb"),
                      verbose = F,
                      control.family = list(hyper = list( theta = list(prior = "loggamma",param = c(1, 0.01)))))


# Coefficients ------------------------------------------------------------
coefs1 <- data.frame(Estimate = c(summary(mod_ind)$coefficients[2,1],
                                 summary(mod_spa_ind5)$beta_table[2,1],
                                 summary(mod_spa_ind25)$beta_table[2,1],
                                 summary(mod_spa_ind50)$beta_table[2,1],
                                 summary(mod_spa_gp5)$beta_table[2,1],
                                 summary(mod_spa_gp25)$beta_table[2,1],
                                 summary(mod_spa_gp50)$beta_table[2,1]
),
SE = c(summary(mod_ind)$coefficients[2,2],
       summary(mod_spa_ind5)$beta_table[2,2],
       summary(mod_spa_ind25)$beta_table[2,2],
       summary(mod_spa_ind50)$beta_table[2,2],
       summary(mod_spa_gp5)$beta_table[2,2],
       summary(mod_spa_gp25)$beta_table[2,2],
       summary(mod_spa_gp50)$beta_table[2,2]
)) %>%
  mutate(LB = Estimate - qnorm(0.975) * SE,
         UB = Estimate + qnorm(0.975) * SE,
         Model = c("Naive","IndepClust5","IndepClust25",#"IndepClust50",
                   "GP5","GP25","GP50"
         )) %>%
  dplyr::select(!SE)

coefs2 <- data.frame(Estimate = c(modSPDEdense$summary.fixed$`0.5quant`[2],
                                  modSPDEcoarse$summary.fixed$`0.5quant`[2],
                                  modBYM2_5$summary.fixed$`0.5quant`[2],
                                  modBYM2_25$summary.fixed$`0.5quant`[2],
                                  modBYM2_50$summary.fixed$`0.5quant`[2]
),
LB = c(modSPDEdense$summary.fixed$`0.025quant`[2],
       modSPDEcoarse$summary.fixed$`0.025quant`[2],
       modBYM2_5$summary.fixed$`0.025quant`[2],
       modBYM2_25$summary.fixed$`0.025quant`[2],
       modBYM2_50$summary.fixed$`0.025quant`[2]
),
UB = c(modSPDEdense$summary.fixed$`0.975quant`[2],
       modSPDEcoarse$summary.fixed$`0.975quant`[2],
       modBYM2_5$summary.fixed$`0.975quant`[2],
       modBYM2_25$summary.fixed$`0.975quant`[2],
       modBYM2_50$summary.fixed$`0.975quant`[2]
),
Model = c("SPDEdense","SPDEcoarse","BYM2_5","BYM2_25","BYM2_50"))

coefs <- rbind(coefs1,coefs2)


# Spatial pattern ---------------------------------------------------------
# INLA
projdense <- fmesher::fm_evaluator(meshdense, loc = as.matrix(dataclust[,c("sdimx","sdimy")]))
field.projdense <- fmesher::fm_evaluate(projdense, modSPDEdense$summary.random$s$mean)

projcoarse <- fmesher::fm_evaluator(meshcoarse, loc = as.matrix(dataclust[,c("sdimx","sdimy")]))
field.projcoarse <- fmesher::fm_evaluate(projcoarse, modSPDEcoarse$summary.random$s$mean)

dataranef_inla <- dataclust %>%
  select(cell_ID, sdimx, sdimy, cluster5, cluster25, cluster50) %>%
  mutate(spdedense = field.projdense,
         spdecoarse = field.projcoarse) %>%
  left_join(data.frame(cluster5 = modBYM2_5$summary.random[[1]]$ID,
                       bym25 = modBYM2_5$summary.random[[1]]$mean)) %>%
  left_join(data.frame(cluster25 = modBYM2_25$summary.random[[1]]$ID,
                       bym225 = modBYM2_25$summary.random[[1]]$mean)) %>%
  left_join(data.frame(cluster50 = modBYM2_50$summary.random[[1]]$ID,
                       bym250 = modBYM2_50$summary.random[[1]]$mean)) 

# Frequentist
dataranef <- dataclust %>%
  dplyr::select(cell_ID, sdimx, sdimy, cluster5, cluster25, cluster50,
                sdimx5,sdimy5,sdimx25,sdimy25,sdimx50,sdimy50) %>%
  mutate(glmres = residuals(mod_ind)) %>%
  left_join(data.frame(cluster5 = names(ranef(mod_spa_ind5, type = "correlated")[[1]]),
                       indep5 = ranef(mod_spa_ind5, type = "correlated")[[1]])) %>%
  left_join(data.frame(cluster25 = names(ranef(mod_spa_ind25, type = "correlated")[[1]]),
                       indep25 = ranef(mod_spa_ind25, type = "correlated")[[1]])) #%>%
  left_join(data.frame(cluster50 = names(ranef(mod_spa_ind50, type = "correlated")[[1]]),
                       indep50 = ranef(mod_spa_ind50, type = "correlated")[[1]])) %>%
  left_join(dataranef_inla)

# Gaussian process
#5
re5 <- ranef(mod_spa_gp5, type = "correlated")[[1]]
coords5 <- do.call(rbind, strsplit(names(re5), ":"))
df_gp5 <- data.frame(
  sdimx5 = round(as.numeric(coords5[,1]),8),
  sdimy5 = round(as.numeric(coords5[,2]),8),
  gp5 = as.numeric(re5)
)

#25
re25 <- ranef(mod_spa_gp25, type = "correlated")[[1]]
coords25 <- do.call(rbind, strsplit(names(re25), ":"))
df_gp25 <- data.frame(
  sdimx25 =  round(as.numeric(coords25[,1]),8),
  sdimy25 =  round(as.numeric(coords25[,2]),8),
  gp25 = as.numeric(re25)
)
#50
re50 <- ranef(mod_spa_gp50, type = "correlated")[[1]]
coords50 <- do.call(rbind, strsplit(names(re50), ":"))
df_gp50 <- data.frame(
  sdimx50 = round(as.numeric(coords50[,1]),5),
  sdimy50 = round(as.numeric(coords50[,2]),5),
  gp50 = as.numeric(re50)
)

dataranef <- dataranef %>%
  mutate(sdimx5 = round(sdimx5,8),
         sdimy5 = round(sdimy5,8),
         sdimx25 = round(sdimx25,8),
         sdimy25 = round(sdimy25,8)) %>%
  left_join(df_gp5) %>%
  left_join(df_gp25)

save(mod_ind,mod_spa_ind5,mod_spa_ind25,mod_spa_ind50,
     mod_spa_gp5,mod_spa_gp25,mod_spa_gp50,
     coefs, dataranef, file = paste0(dir.out, "Examples/Example","_",sim,".Rdata"))
