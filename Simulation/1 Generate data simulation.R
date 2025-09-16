#-------------------#
#-- Simulate data --#
#-------------------#
#https://nanostring.com/resources/smi-ffpe-dataset-giotto-object/
dir.out = "~/Simulation/"

library(Giotto)
library(tidyverse)
# Data --------------------------------------------------------------------
load(paste0(dir.out,"Processed Data Giotto Object/SMI_Giotto_Object.RData"))
all.equal(rownames(t(gem@expression$rna$raw)), gem@cell_metadata$rna$cell_ID)
lung5 <- gem@cell_metadata$rna$cell_ID[gem@cell_metadata$rna$Run_Tissue_name == "Lung5_Rep3"]
counts <- as.matrix(t(gem@expression$rna$raw)[lung5,])
fulldata <- gem@cell_metadata$rna %>%
  left_join(gem@spatial_locs$raw) %>% 
  filter(cell_ID %in% lung5) %>%
  cbind(.,counts)
save(fulldata,file = paste0(dir.out,"Lung5-3data.Rdata") )

data <-  fulldata %>% mutate(auxID = 1:nrow(fulldata)) %>%
  filter(cell_type == "macrophage",
         niche %in% c("myeloid-enriched stroma","stroma"))
data <- data %>% mutate(niche = as.factor(niche),
                        nichenum = as.numeric(niche ==  "stroma" ))

# Spatial random fields
MLdistMat <- telefit::maternCov(dist(fulldata[,c("sdimx", "sdimy")]),
                               smoothness = 0.5, scale = 1)
set.seed(4315)
intensities <- t(MASS::mvrnorm(n=5000, mu = rep(0,nrow(fulldata)), Sigma=MLdistMat))
save(intensities, paste0(dir.out,"intensities.Rdata"))

# Simulation parameters

params <- expand.grid(beta = c(rep(0,5), seq(0,1.5, by = 0.1)), # For power c(seq(0,0.5, by = 0.03))
                      sim = 1:100) %>% group_by(sim) %>% 
  dplyr::mutate(id = 1:n(),
                aux = as.numeric(paste0(sim,0,id)))

params <- params %>%
  left_join(data.frame(aux = unique(params$aux),
                       id2 = 1:length(unique(params$aux))),
            by = "aux") %>%
  mutate(gene_name = paste("Gene",beta,sim,id,sep="_"))

# Simulate counts
Ycounts <- do.call(cbind, lapply(1:nrow(params), function(args){
  sim <- params$sim[args]
  beta <- params$beta[args]
  id <- params$id2[args]
  genename <- sim
  
  lambda <- data$totalcounts * 0.01 * exp(beta * data$nichenum + intensities[,id])
  set.seed(67983*id)
  Ycount <- rpois(nrow(fulldata),lambda)
  return(Ycount)
}))
colnames(Ycounts) <- params$gene_name
p <- ncol(Ycounts)

# Normalized counts
Ygaus <- do.call(cbind, lapply(1:p, function(gene){
  Ycounts[,gene]/data$totalcounts * mean(data$totalcounts)
}))
colnames(Ygaus) <- colnames(Ycounts)
save(Ycounts, Ygaus,params,data, file=paste(dir.out,paste0("SimulatedY.Rdata")))

