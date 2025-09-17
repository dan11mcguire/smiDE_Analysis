#------------------------#
#-- Results Simulation --#
#------------------------#

dir.out = "~/Simulation/"
library(tidyverse)
library(cowplot)
load(file = paste(dir.out,"SimulatedY.Rdata"))

# Combining simulation results --------------------------------------------
## Cside ---
load(paste0(dir.out, "Results/PartialSimulation/ResultsCside.Rdata"))
df_cside <- results_cside$result@de_results$all_gene_list$main %>%
  rownames_to_column(var = "gene_name") %>%
  left_join(params) %>%
  mutate(Sign = ifelse(p_val < 0.05, 1, 0),
         kprop = 1,
         Model = "C-SIDE",
         Estimate = log_fc,
         pvalue = p_val,
         Gene = sim,
         Distribution = "Negative Binomial",
         Time = results_cside$memory$Elapsed_Time_sec/2100,
         TotalRAM = results_cside$memory$Total_RAM_Used_MiB,
         PeakRAM = results_cside$memory$Peak_RAM_Used_MiB,
  ) %>%
  dplyr::select(Gene,gene_id=gene_name, beta, sim,id, Model, Distribution,Estimate,pvalue, Sign, Time, TotalRAM, PeakRAM, kprop)

## Naive  ---
load(paste0(dir.out, "Results/PartialSimulation/ResultsNaiveModels.Rdata"))
df_naive <- results_naive %>%
  left_join(params, by = c("gene_id"="gene_name")) %>%
  mutate(Sign = ifelse(p_val < 0.05, 1, 0),
         kprop = 1,
         Time = Time/2100,
         model = paste("Naive",model),
         Estimate = avg_log2FC,
         pvalue = p_val,
         Gene = sim,
         dist = factor(dist, levels = c("negbin","gaussian"),
                       labels = c("Negative Binomial","Gaussian"))) %>%
  filter(model %in% c("Naive wilcox","Naive MAST","Naive deseq2","Naive negbinom"))%>%
  dplyr::select(Gene,gene_id, beta, sim,id,Model = model, Distribution=dist,Estimate,pvalue, Sign, Time, TotalRAM, PeakRAM, kprop)


## INLA  ---
files_inla <- list.files( paste0(dir.out,"Results/PartialSimulation/INLA/"),pattern = '\\.csv')
df_inla0 <- do.call(rbind, lapply(files_inla, function(file)data.table::fread(paste0( paste0(dir.out,"Results/PartialSimulation/INLA/"),file))))[,-1]

df_inla <- df_inla0 %>% left_join(params, by = c("Gene"="gene_name")) %>%
  filter(Variable %in% c("nichestroma","niche")) %>%
  mutate(Sign = ifelse(UB < 0 | LB > 0, 1, 0), Gene = sim,
         pvalue = pseudopvalue,
         Distribution = factor(Distribution, levels = c("negbin","gaussian"),
                               labels = c("Negative Binomial","Gaussian"))) %>%
  dplyr::select(Gene, beta, sim,id, Model, Distribution,
                mesh, prior_range, prior_sigma,Estimate,pvalue, Sign, Time, TotalRAM, PeakRAM)


# BYM2  ---
files_bym2 <- list.files(paste0(dir.out,"Results/PartialSimulation/BYM2/"),pattern = '\\.csv')
df_bym20 <- do.call(rbind, lapply(files_bym2, function(file)data.table::fread(paste0( paste0(dir.out,"Results/PartialSimulation/BYM2/"),file))))[,-1]

df_bym2 <- df_bym20 %>% left_join(params2, by = c("Gene"="gene_name")) %>%
  filter(Variable %in% c("niche")) %>%
  mutate(Sign = ifelse(UB < 0 | LB > 0, 1, 0),
         pvalue = pseudopvalue,
         Gene = sim,
         Distribution = factor(Distribution, levels = c("negbin","gaussian"),
                               labels = c("Negative Binomial","Gaussian"))) %>%
  select(Gene, beta, sim,id, Model, Distribution,prop=kprop,AdjacencyMethod, prior_phi, prior_prec,Estimate,pvalue, Sign, Time, TotalRAM, PeakRAM)


# spaMM  ---
files_spamm <- list.files(paste0(dir.out,"Results/PartialSimulation/spaMM/"),pattern = '\\.csv')
df_spamm0 <- bind_rows(lapply(files_spamm, function(file) {data.table::fread(paste0(paste0(dir.out,"Results/PartialSimulation/spaMM/"),file))}))

df_spamm <- df_spamm0 %>% left_join(params, by = c("Gene"="gene_name")) %>%
  filter(Variable %in% c("nichestroma","niche")) %>%
  mutate(Sign = ifelse(pvalue < 0.05, 1, 0),
         Gene = sim,
         Distribution = factor(Distribution, levels = c("negbin","gaussian"),
                               labels = c("Negative Binomial","Gaussian"))) %>%
  dplyr::select(Gene, beta, sim,id, Model, Distribution,prop=kprop,Estimate,pvalue, Sign, Time, TotalRAM, PeakRAM)


# Combine data  ---
dataplot <- bind_rows(df_bym2,df_inla,df_spamm,df_naive,df_cside ) %>%
  dplyr::group_by(Model, prop, Distribution, Gene, prior_phi, prior_prec,
                  mesh, prior_range, prior_sigma) %>%
  dplyr::mutate(Rank = cor(beta, Estimate, method = "spearman",use="pairwise.complete.obs"),
                Type1Error = ifelse(beta == 0,Sign,NA )
  ) %>%
  dplyr::ungroup(Gene) %>%
  dplyr::summarise(Type1Error = mean(Type1Error,na.rm = T),
                   Rank = mean(Rank, na.rm = T),
                   Time = mean(Time, na.rm = T)) %>%
  filter(Distribution %in% c("Negative Binomial", "Gaussian")) %>%
  mutate(label = ifelse(Model)) %>%
  filter((prior_range == "standard" | is.na(prior_range)),
         (prior_sigma == "standard" | is.na(prior_sigma)),
         (prior_phi == "balanced" | is.na(prior_phi)),
         (prior_prec == "standard" | is.na(prior_prec))) %>%
  mutate(Model=factor(Model, levels = c("Naive deseq2","Naive negbinom",
                                        "Naive wilcox","Naive MAST",
                                        "C-SIDE","Spatial GP","Spatial Indep","INLA-GP","BYM2"),
                      labels = c("DESeq2","NB GLM","Wilcox","MAST",
                                 "C-SIDE","GP","Independent\nclusters","SPDE-INLA","BYM2"))
  ) 
write.csv(dataplot, paste0(dir.out, "ResultsSimulation.Rdata"))
