#-------------------------#
#- Combine results -#
#-------------------------#
dir <-  "~MERFISH/"

library(tidyverse)
library(smiDE)
library(ggrepel)
library(nicheDE)

load(file = paste0(dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
filtered_genes <- overlap_metrics_allen$result %>%
  filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1, ratio < 1) %>% pull(target)

# smiDE -------------------------------------------------------------------
## Naive ----
# Negative Binomial
load(paste0(dir,"ResultsDistance/NaivesmiDE/Naive_smiDE_nb_kernelMicroglia1000.Rdata"))

smide_naive_nb_wb <- do.call(rbind, lapply(naivesmiDE_nb$result,
                                           function(x)smiDE::results(x, comparisons = "model_summary", variable = "DEvar")[[1]])) %>%
  filter(term == "DEvar") %>%
  mutate(Distribution =  "negbin",
         Type = "Naive",
         Method = "smiDE",
         DEvar="kernelMicroglia1000",
         Time = naivesmiDE_nb$memory$Elapsed_Time_sec,
         PeakRAM = naivesmiDE_nb$memory$Peak_RAM_Used_MiB,
         TotalRAM = naivesmiDE_nb$memory$Total_RAM_Used_MiB,
         Estimate = est,
         kprop = 1) %>%
  select(target, Type, Distribution,kprop, Method, DEvar, Estimate, pval,Time, PeakRAM, TotalRAM) %>%
  dplyr::rename(pvalue =pval, Gene = target)

# Gaussian
load(paste0(dir,"ResultsDistance/NaivesmiDE/Naive_smiDE_gaus_kernelMicroglia1000.Rdata"))
smide_naive_gaus_wb <- do.call(rbind,lapply(naivesmiDE_gaus$result,
                                            function(x)smiDE::results(x, comparisons = "model_summary", variable = "DEvar")[[1]])) %>%
  filter(term == "DEvar") %>%
  mutate(Distribution =  "gaussian",
         Type = "Naive",
         Method = "smiDE",
         DEvar="kernelMicroglia1000",
         Time = naivesmiDE_gaus$memory$Elapsed_Time_sec,
         PeakRAM = naivesmiDE_gaus$memory$Peak_RAM_Used_MiB,
         TotalRAM = naivesmiDE_gaus$memory$Total_RAM_Used_MiB,
         Estimate = est,
         kprop = 1) %>%
  select(target, Type, Distribution,kprop, Method, DEvar, Estimate, pval,Time, PeakRAM, TotalRAM) %>%
  dplyr::rename(pvalue = pval, Gene = target)


# Semi Naive ----
# Negative Binomial
load(paste0(dir,"ResultsDistance/NaivesmiDE/semiNaive_smiDE_nb_kernelMicroglia1000.Rdata"))

smide_seminaive_nb_wb <- do.call(rbind,lapply(seminaivesmiDE_nb$result,
                                              function(x)smiDE::results(x, comparisons = "model_summary", variable = "DEvar")[[1]])) %>%
  filter(term == "DEvar") %>%
  mutate(Distribution =  "negbin",
         Type = "Semi naive",
         Method = "smiDE",
         DEvar="kernelMicroglia1000",
         Time = seminaivesmiDE_nb$memory$Elapsed_Time_sec,
         PeakRAM = seminaivesmiDE_nb$memory$Peak_RAM_Used_MiB,
         TotalRAM = seminaivesmiDE_nb$memory$Total_RAM_Used_MiB,
         Estimate = est,
         kprop = 1) %>%
  select(target, Type, Distribution,kprop, Method, DEvar, Estimate, pval,Time, PeakRAM, TotalRAM) %>%
  dplyr::rename(pvalue = pval, Gene = target)

# Gaussian
load(paste0(dir,"ResultsDistance/NaivesmiDE/semiNaive_smiDE_gaus_kernelMicroglia1000.Rdata"))

smide_seminaive_gaus_wb <- do.call(rbind,lapply(seminaivesmiDE_gaus$result,
                                                function(x)smiDE::results(x, comparisons = "model_summary", variable = "DEvar")[[1]])) %>%
  filter(term == "DEvar") %>%
  mutate(Distribution =  "gaussian",
         Type = "Semi Naive",
         Method = "smiDE",
         DEvar="kernelMicroglia1000",
         Time = seminaivesmiDE_gaus$memory$Elapsed_Time_sec,
         PeakRAM = seminaivesmiDE_gaus$memory$Peak_RAM_Used_MiB,
         TotalRAM = seminaivesmiDE_gaus$memory$Total_RAM_Used_MiB,
         Estimate = est,
         kprop = 1) %>%
  select(target,Type, Distribution,kprop, Method, DEvar, Estimate, pval,Time, PeakRAM, TotalRAM) %>%
  dplyr::rename(pvalue = pval, Gene = target)

results_smiDE_wb <- rbind(smide_naive_nb_wb,
                          smide_naive_gaus_wb,
                          smide_seminaive_nb_wb,
                          smide_seminaive_gaus_wb)

# spaMM -------------------------------------------------------------------
filepath_spamm <- paste0(dir,"ResultsDistance/spaMM/")
files_spamm <- list.files(filepath_spamm,pattern = '\\.csv')
results_spamm_wb <- bind_rows(lapply(files_spamm, 
                                     function(f) {data.table::fread(file.path(filepath_spamm, f), 
                                                                    colClasses = list(character = "V1"))}))  %>%
  filter(Variable == "DEvar" | is.na(Variable)) %>%
  select(Gene,Type, Distribution,kprop, Method, DEvar, Estimate, pvalue,Time, PeakRAM, TotalRAM,failed)

# INLA --------------------------------------------------------------------
filepath_inla <- paste0(dir,"ResultsDistance/INLA/")
files_inla <- list.files(filepath_inla,pattern = '\\.csv')
results_inla_wb <- bind_rows(lapply(files_inla, function(f) {
  data.table::fread(file.path(filepath_inla, f),
                    colClasses = list(character = "rangeest"),
                    na.strings = c("NA", "", ".")) %>% mutate(overdis = as.numeric(overdis))}))  %>%
  filter(Variable == "DEvar"|is.na(Variable)) %>%
  mutate(Method = "INLA",
         pvalue = ifelse(LB > 0 | UB < 0, 0, 1),
         padj =  ifelse(LBbonf_cont > 0 | UBbonf_cont < 0, 0, 1)) %>%
  select(Gene,Type, Distribution,kprop, Method, DEvar , Estimate, pvalue,padj,pseudopvalue,Time, PeakRAM, TotalRAM,failed,meshscenario,prior_range, prior_sigma )

# BYM2 --------------------------------------------------------------------
filepath_bym2 <- paste0(dir,"ResultsDistance/BYM2/")
files_bym2 <- list.files(filepath_bym2,pattern = '\\.csv')
results_bym2_wb <- bind_rows(lapply(files_bym2, function(f) {data.table::fread(file.path(filepath_bym2, f), colClasses = list(character = "V1"))}))  %>%
  filter(Variable == "DEvar" | is.na(Variable)) %>%
  mutate(Method = "BYM2",
         pvalue = ifelse(LB > 0 | UB < 0, 0, 1),
         padj =  ifelse(LBbonf_cont > 0 | UBbonf_cont < 0, 0, 1)) %>%
  select(Gene,Type, Distribution,kprop, Method, DEvar, Estimate, pvalue,padj,
         pseudopvalue,Time, PeakRAM, TotalRAM,failed, prior_phi, prior_prec,failed)

# NicheDE --------------------------------------------------------------------
library(nicheDE)
results_nicheDE_wb <- do.call(rbind, lapply(c("negbin","gaus"), function(dist){
  load(file = paste0(dir,"ResultsDistance/NicheDE/nicheDE_Microglia_1000_",dist,".Rdata"))
  get_niche_DE_genes(NDE_obj$result,test.level = "I",
                     index='Neurons', niche = 'Microglia',positive = T,alpha = 1) %>%
    mutate(Estimate = 1) %>%
    rbind(get_niche_DE_genes(NDE_obj$result,
                             test.level = "I",
                             index='Neurons',
                             niche = 'Microglia',positive = F,alpha = 1)%>%
            mutate(Estimate = -1)) %>%
    group_by(Genes) %>%
    summarise(Pvalues.Interaction = min(Pvalues.Interaction),
              Estimate = mean(Estimate)) %>%
    mutate(Distribution = dist,
           Type = "nicheDE",
           Method = "nicheDE",
           DEvar = "kernelMicroglia1000",
           Time = NDE_obj$memory$Elapsed_Time_sec,
           PeakRAM = NDE_obj$memory$Peak_RAM_Used_MiB,
           TotalRAM = NDE_obj$memory$Total_RAM_Used_MiB,
           kprop = 1) %>%
    select(Gene = Genes,Type, Distribution,kprop, Method, DEvar, Estimate,
           padj = Pvalues.Interaction,Time, PeakRAM, TotalRAM)
}))


# Cside -------------------------------------------------------------------
load(paste0(dir, "ResultsDistance/ResultsCside.Rdata"))
df_cside <- results_cside$result@de_results$all_gene_list$`01 IT-ET Glut` %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Type = "Cside",
         Distribution = "negbin",
         kprop = 1,
         Method = "C-SIDE",
         DEvar = "kernelMicroglia1000",
         Estimate = log_fc,
         padj = p_val,
         Time = results_cside$memory$Elapsed_Time_sec,
         TotalRAM = results_cside$memory$Total_RAM_Used_MiB,
         PeakRAM = results_cside$memory$Peak_RAM_Used_MiB,
  ) %>%
  dplyr::select(Gene, Type,Distribution,kprop, Method,DEvar, Estimate,padj, Time, TotalRAM, PeakRAM, kprop)

# Combining results -------------------------------------------------------
genes <- read_csv(paste0(dir, "Allen/gene.csv"))
genes$gene_identifier  <- as.character(genes$gene_identifier )

results_all<- bind_rows(results_smiDE_wb,
                        results_spamm_wb) %>%
  filter(!((Type != "Naive") & !(Gene %in% filtered_genes))) %>%
  group_by(Method, Type, kprop, Distribution,DEvar) %>%
  mutate(padj = p.adjust(pvalue, method = "bonferroni"),
         pvalue = p.adjust(pvalue, method = "BH")) %>%
  bind_rows(results_nicheDE_wb, results_inla_wb,results_bym2_wb,df_cside) %>%
  filter(failed == F | is.na(failed)) %>%
  mutate(DE = as.numeric(padj < 0.05),
         DEfdr = as.numeric(pvalue < 0.05)
  ) %>%
  filter(Distribution %in% c("gaus","gaussian","negbin") ) %>%
  mutate(Distribution = factor(Distribution, levels = c("gaus","gaussian","negbin"),
                               labels = c("Gaussian","Gaussian","Negative Binomial"))) %>%
  mutate(Type = factor(Type, levels = c("Naive","Semi naive", "Semi Naive", "full","nichdeDE" ),
                       labels = c("Naive","Seminaive", "Seminaive", "SpatiallyAware","nichdeDE" )))%>%
  mutate(Method = ifelse(Method == "smiDE" & Type == "Naive", "smiDE Naive",
                         ifelse(Method == "smiDE" & Type == "Seminaive", "smiDE Seminaive",Method)),
         Method = factor(Method, levels=c("smiDE Naive","smiDE Seminaive","IndepCluster","BYM2","GP","INLA","sdmTMB","nicheDE","C-SIDE"),
                         labels = c("smiDE Naive","smiDE Seminaive","Independent\nclusters","BYM2","GP","SPDE","GP-sdmTMB","nicheDE","C-SIDE")))%>%
  mutate(Analysis = ifelse((Method %in% c("smiDE Naive","smiDE Seminaive","nicheDE","C-SIDE")) |
                             (kprop == 0.5 & Method == "Independent\nclusters") |
                             (kprop == 0.5 & Method == "GP") | 
                             (Method == "SPDE" & meshscenario %in% c(0.5) &
                                prior_range == "standard" & prior_sigma == "standard") |
                             (Method == "BYM2" & kprop %in% c(1) &
                                prior_phi == "balanced" & prior_prec == "standard"),
                           "Main","Secondary")
  )%>% filter(DEvar == "kernelMicroglia1000") %>%
  filter(!(Method == "Independent\nclusters" & kprop == 1),
         meshscenario %in% c(0.5,3)| is.na(meshscenario))%>%
  left_join(genes[,1:2], by = c("Gene" = "gene_identifier"))
write.csv(results_all, file = paste0(dir,"Results_realdata_Microglia100_all.csv"))
