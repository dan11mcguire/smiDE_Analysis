#-------------------------#
#- Get real data results -#
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


# Figures -----------------------------------------------------------------
library(tidyverse)
library(smiDE)
library(ggrepel)
library(nicheDE)


# Neurons and Microglia ---------------------------------------------------
load( paste0(dir, "Allen/WholeBrain_seurat.Rdata"))
colnames(seurat_allen@meta.data)[4] <- "cell_ID"
dataDE2 <- subset(seurat_allen, subset = (class %in% c("01 IT-ET Glut")|subclass == "334 Microglia NN"))
dataDE2 <- subset(dataDE2, subset = division == "Isocortex")

center_id <- "1017155956100820114-1"
# kernel density
load( paste0(dir, "Allen/WholeBrain_seurat_distance.Rdata"))
load(file = paste0(dir,"ResultsDistance/NicheDE/nicheDE_Microglia_1000_negbin",".Rdata"))
df <- dataDE2@meta.data
cx <- df$x[df$cell_ID == center_id]; cy <- df$y[df$cell_ID ==center_id]
circle_df <- function(cx, cy, r, n = 360) {
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = cx + r*cos(t), y = cy + r*sin(t))
}

ring <- circle_df(cx, cy, 1000/NDE_obj$result@scale)
plot_circle <- dataDE@meta.data %>%filter(x > 3.5, x<5.5, y>3, y<4)%>%
  ggplot(aes( x= x, y = y)) +
  geom_point(cex = .5,aes(col = "Excitatory neuron", pch = "Excitatory neuron")) +
  geom_path(data = ring, aes(x, y), inherit.aes = FALSE, lty = 2) +
  geom_point(data = data.frame(x = cx, y = cy), aes(x, y),
             shape = 19, fill = "white", size = 1, stroke = 0.7,
             inherit.aes = FALSE) +
  geom_point(data = dataDE2@meta.data %>% filter(class == "34 Immune",
                                                 x > 3.5, x<5.5, y>3, y<4),
             cex=0.7,
             aes(col = "Microglia", pch = "Microglia")) +
  theme_bw() +coord_fixed()+
  scale_color_manual(values = c("darkgrey","red")) +
  scale_shape_manual(values = c(19,4))+
  labs(col = "Cell type", pch = "Cell type") +
  theme(strip.background = element_blank(),
        text = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size = 2)));plot_circle
ggsave(paste0(dir, "Images/Final/kernelVascular1000_circle.jpeg"), height = 90, width = 180, units = "mm")

plot_density <- dataDE@meta.data %>%
  ggplot(aes( x= x, y = y, col = kernelMicroglia1000)) +
  geom_point(cex = .2) +
  scale_color_viridis_c()+
  labs(col = "Density")+
  theme_bw() +coord_fixed()+
  theme(strip.background = element_blank())+
  theme(text = element_text(size = 11)) ;plot_density
ggsave(paste0(dir, "Images/Final/kernelMicroglia1000.jpeg"), height = 90, width = 180, units = "mm")



# Results -----------------------------------------------------------------
my_colors <- c( "#724690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05",
                "#CC503E", "#94346E", "#6F4070", "#994E95","#666666")

## nDE by time --------------------------------------------------------------------
plot_nDE <- (results_all %>%
               group_by(Method, Type, kprop, Distribution,DEvar, meshscenario,
                        prior_range, prior_sigma, prior_phi, prior_prec) %>%
               summarise(nDE = sum(DE, na.rm = T),
                         time = mean(Time, na.rm = T)) %>%
               filter(prior_range == "standard" |is.na(prior_range) ,
                 prior_sigma == "standard"|is.na(prior_sigma),
                 prior_phi == "balanced" |is.na(prior_phi) | 
                   (prior_phi == "iid" & Method == "BYM2" & kprop == 1 & Distribution == "Gaussian") ,
                 prior_prec == "standard" |is.na(prior_prec)) %>%
               mutate(label = ifelse(Method %in% c("GP-sdmTMB","SPDE"),
                                     ifelse(meshscenario == 0.5, "Dense","Coarse"),
                                     ifelse(Method %in% c("Independent\nclusters","GP","BYM2"),
                                            paste(kprop*100,"%"),"")))%>%
               ggplot(aes(x = time, y = nDE, col = Method, group= paste(Method, Distribution),
                          label = paste(meshscenario,kprop))) +
               geom_point(cex=2) + 
               theme_bw()+
               facet_wrap(Distribution~.)+
               geom_label_repel(aes(label = label),show.legend = FALSE, size = 3,max.overlaps = Inf,
                                box.padding = 0.5, label.padding = 0.1, min.segment.length = 0.1) +
               scale_x_continuous(trans=scales::log_trans(),
                                  breaks = c(0, 60,60*2,60*5,60*10,60*30,60*60,2*60*60,4*60*60,6*60*60),
                                  labels = c("0s","1m","2m","5m","10m","30m","1h","2h","4h","6h")) +
               scale_color_manual(limits = c("C-SIDE","nicheDE","smiDE Naive","smiDE Seminaive",
                                             "Independent\nclusters",
                                             "BYM2","GP","SPDE"),
                                  values = my_colors[c(1,12,2,11,3,5,7,8)]) +
               scale_shape_manual(values = c(19, 3)) +
               labs(x = "Time", y= "Number DE genes"))+
  theme(text = element_text(size = 11),
        strip.background = element_blank());plot_nDE
ggsave(paste0(dir, "Images/Final/nDEbytime.jpeg"), height = 120, width = 220, units = "mm")

## Upset --------------------------------------------------------------------
library(ComplexUpset)

sig_sets_nb <- results_all %>%
  filter(Analysis == "Main",
         Distribution == "Negative Binomial") %>%
  group_by(Method) %>%
  summarise(
    n_total = n_distinct(gene_symbol),
    genes_sig = list(unique(gene_symbol[padj <= 0.05])),
    n_sig = lengths(genes_sig),
    .groups = "drop"
  )

x_nb = list( `C-SIDE` = sig_sets_nb$genes_sig[[8]],
             `nicheDE` = sig_sets_nb$genes_sig[[7]],
             `smiDE Naive` = sig_sets_nb$genes_sig[[1]],
             `smiDE SemiNaive` = sig_sets_nb$genes_sig[[2]],
             `Independent\nclusters (50%)` = sig_sets_nb$genes_sig[[3]],
             `BYM2 (100%)` = sig_sets_nb$genes_sig[[4]],
             `GP (50%)` = sig_sets_nb$genes_sig[[5]],
             `SPDE (Dense)` = sig_sets_nb$genes_sig[[6]]
             
)

# upset
incidence <- UpSetR::fromList(x_nb) %>%           
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(across(-gene, as.logical))
plot_upset <- ComplexUpset::upset(
  incidence,
  intersect = names(x_nb),
  height_ratio=1,
  width_ratio=0.2,
  set_sizes=F,
  encode_sets = FALSE, 
  queries = list(
    ComplexUpset::upset_query(
      intersect = c("nicheDE","smiDE Naive","smiDE SemiNaive","Independent\nclusters (50%)",
                    "BYM2 (100%)","C-SIDE",
                    "GP (50%)","SPDE (Dense)"),   
      color = "red",
      fill  = "red",
      only_components = c("intersections_matrix")
    ),
    ComplexUpset::upset_query(
      intersect = c("nicheDE","smiDE Naive","smiDE SemiNaive","Independent\nclusters (50%)",
                    "C-SIDE",
                    "GP (50%)","SPDE (Dense)"),   
      color = "red",
      fill  = "red",
      only_components = c("intersections_matrix")
    )
  )
)+
  labs(x = "Gene sets", y = "Intersection size") +
  theme(
    axis.text.y = element_text(size = 12),   
    axis.title.x = element_text(size = 12)
  );plot_upset
ggsave(paste0(dir, "Images/Final/UpsetNB.jpeg"), height = 160, width = 230, units = "mm")

# Gaussian
sig_sets_g <- results_all %>%
  filter(Analysis == "Main",
         Distribution != "Negative Binomial") %>%
  group_by(Method) %>%
  summarise(
    n_total = n_distinct(gene_symbol),
    genes_sig = list(unique(gene_symbol[padj <= 0.05])),
    n_sig = lengths(genes_sig),
    .groups = "drop"
  )

x_gau = list( `nicheDE` = sig_sets_g $genes_sig[[7]],
              `smiDE Naive` = sig_sets_g $genes_sig[[1]],
              `smiDE SemiNaive` = sig_sets_g $genes_sig[[2]],
              `Independent\nclusters (50%)` = sig_sets_g $genes_sig[[3]],
              `BYM2 (100%)` = sig_sets_g $genes_sig[[4]],
              `GP (50%)` = sig_sets_g $genes_sig[[5]],
              `SPDE (Dense)` = sig_sets_g $genes_sig[[6]]
)

incidence_gau <- UpSetR::fromList(x_gau) %>%         
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(across(-gene, as.logical))
ComplexUpset::upset(
  incidence_gau,
  intersect = names(x_gau),
  height_ratio=1,
  width_ratio=0.2,
  set_sizes=F,
  encode_sets = FALSE, 
  queries = list(
    ComplexUpset::upset_query(
      intersect = c("nicheDE","smiDE Naive","smiDE SemiNaive","Independent\nclusters (50%)",
                    "BYM2 (100%)",
                    "GP (50%)","SPDE (Dense)"), 
      color = "red",
      fill  = "red",
      only_components = c("intersections_matrix")
    ),
    ComplexUpset::upset_query(
      intersect = c("nicheDE","smiDE Naive","smiDE SemiNaive",
                    "BYM2 (100%)","SPDE (Dense)"),   
      color = "red",
      fill  = "red",
      only_components = c("intersections_matrix")
    ),
    ComplexUpset::upset_query(
      intersect = c("smiDE Naive","smiDE SemiNaive","Independent\nclusters (50%)",
                    "BYM2 (100%)","GP (50%)","SPDE (Dense)"),   
      color = "red",
      fill  = "red",
      only_components = c("intersections_matrix")
    ),
    ComplexUpset::upset_query(
      intersect = c("smiDE Naive","smiDE SemiNaive","nicheDE","SPDE (Dense)"), 
      color = "red",
      fill  = "red",
      only_components = c("intersections_matrix")
    )
    
    
  )
)+
  labs(x = "Gene sets", y = "Intersection size")
ggsave(paste0(dir, "Images/Final/UpsetGaus.jpeg"), height = 120, width = 200, units = "mm")


# Rank --------------------------------------------------------------------
df_rank <- results_all %>% ungroup%>%
  filter(Distribution == "Negative Binomial",
         DEvar == "kernelMicroglia1000",
         prior_range == "standard" |is.na(prior_range) ,
         prior_sigma == "standard"|is.na(prior_sigma),
         prior_phi == "balanced" |is.na(prior_phi) ,
         prior_prec == "standard"|is.na(prior_prec),
         meshscenario %in% c(0.5,1,2,3)|is.na(meshscenario)
  ) %>%
  mutate(pseudopvalue = ifelse(is.na(pseudopvalue),padj, pseudopvalue)) %>%
  mutate(group_id = ifelse(is.na(meshscenario),
                           ifelse(is.na(kprop), Method,
                                  paste0(Method,paste0(" (",kprop*100,"%)"))),
                           paste0(Method,ifelse(meshscenario == 0.5, " (Dense)"," (Coarse)")))) %>%
  filter(!is.na(pseudopvalue)) %>%
  mutate(score = -log10(pseudopvalue)) %>%
  select(Gene, group_id, score)
levels(as.factor(df_rank$group_id))

mat_wide <- df_rank %>% tidyr::pivot_wider(names_from = group_id, values_from = score)
gene_ids <- mat_wide$Gene
mat <- mat_wide %>% select(-Gene) %>% as.data.frame()
apply(mat,2,function(x)order)
cor_mat <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
cor_df <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("group_1", "group_2", "corr")

cor_df$group_1 <- factor(cor_df$group_1,
                         levels = c("C-SIDE (100%)","nicheDE (100%)","smiDE Naive (100%)","smiDE Seminaive (100%)",
                                    "Independent\nclusters (5%)",
                                    "Independent\nclusters (25%)",
                                    "Independent\nclusters (50%)",
                                    "BYM2 (5%)",
                                    "BYM2 (25%)",
                                    "BYM2 (50%)",
                                    "BYM2 (100%)",
                                    "GP (5%)" ,
                                    "GP (25%)" ,
                                    "GP (50%)" ,
                                    "SPDE (Coarse)",
                                    "SPDE (Dense)"
                         ),
                         labels = c("C-SIDE","nicheDE","smiDE Naive","smiDE Seminaive",
                                    "Independent\nclusters (5%)",
                                    "Independent\nclusters (25%)",
                                    "Independent\nclusters (50%)",
                                    "BYM2 (5%)",
                                    "BYM2 (25%)",
                                    "BYM2 (50%)","BYM2 (100%)",
                                    "GP (5%)" ,
                                    "GP (25%)" ,
                                    "GP (50%)" ,
                                    "SPDE (Coarse)",
                                    "SPDE (Dense)"
                         ))
cor_df$group_2 <- factor(cor_df$group_2,
                         levels = c("C-SIDE (100%)","nicheDE (100%)","smiDE Naive (100%)","smiDE Seminaive (100%)",
                                    "Independent\nclusters (5%)",
                                    "Independent\nclusters (25%)",
                                    "Independent\nclusters (50%)",
                                    "BYM2 (5%)",
                                    "BYM2 (25%)",
                                    "BYM2 (50%)","BYM2 (100%)",
                                    "GP (5%)" ,
                                    "GP (25%)" ,
                                    "GP (50%)" ,
                                    "SPDE (Coarse)",
                                    "SPDE (Dense)"
                         ),
                         labels = c("C-SIDE","nicheDE","smiDE Naive","smiDE Seminaive",
                                    "Independent\nClusters (5%)",
                                    "Independent\nClusters (25%)",
                                    "Independent\nClusters (50%)",
                                    "BYM2 (5%)",
                                    "BYM2 (25%)",
                                    "BYM2 (50%)","BYM2 (100%)",
                                    "GP (5%)" ,
                                    "GP (25%)" ,
                                    "GP (50%)" ,
                                    "SPDE (Coarse)",
                                    "SPDE (Dense)"
                         ))

plot_rank <- ggplot(cor_df, aes(group_1, group_2, fill = corr)) +
  geom_tile() +
  geom_label(aes(label = sprintf("%.2f", corr)), size = 4, fill = "white") +
  scale_fill_gradient2(limits = c(-1, 1), na.value = "grey90") +
  coord_equal() +
  scale_y_discrete(limits=rev)+
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        text = element_text(size = 15));plot_rank
ggsave(paste0(dir, "Images/Final/RanksCorrelation_NB.jpeg"), height = 200, width = 200, units = "mm")



# Volcano naive -----------------------------------------------------------
load(paste0(dir,"ResultsDistance/NaivesmiDE/Naive_smiDE_nb_kernelMicroglia1000.Rdata"))

smide_naive_nb <- do.call(rbind, lapply(naivesmiDE_nb$result,
                                        function(x)smiDE::results(x,  variable = "DEvar")[[1]])) %>%
  left_join(genes[,1:2], by = c("target" = "gene_identifier")) %>%
  mutate(padj = p.adjust(p.value, method = "bonferroni"))

signNaive <- results_all %>% filter(Method == "smiDE Naive",
                                    Distribution == "Negative Binomial",
                                    DE ==1) %>% pull(Gene)
signSemiNaive <-  results_all %>% filter(Method == "smiDE Seminaive",
                                         Distribution == "Negative Binomial",
                                         DE ==1) %>% pull(Gene)
signSpatial <-  results_all %>% filter(Method == "SPDE",
                                       Distribution == "Negative Binomial",
                                       meshscenario == 0.5,
                                       DE ==1) %>% pull(Gene)
removed_filtering <- signNaive[!signNaive %in% filtered_genes]
removed_contamination <- signNaive[!signNaive %in% signSemiNaive & ! signNaive %in% removed_filtering]
removed_spatial <- signNaive[!(signNaive %in% signSpatial) & !(signNaive %in% removed_contamination) &
                               !(signNaive %in% removed_filtering)]

plot_volcano <- smide_naive_nb %>%
  mutate(label = ifelse(target %in% signNaive, "Significant", "Not signficant"),
         label = ifelse(target %in% removed_contamination, "Removed: contamination (n=19)", label),
         label = ifelse(target %in% removed_filtering, "Removed: filtering (n=56)", label),
         label = ifelse(target %in% removed_spatial, "Removed: spatial (n=100)", label),
         label = ifelse(target %in% signSpatial, "SPDE DE genes (n=3)",label)
  ) %>%
  mutate(label = factor(label,levels = c("Not signficant",
                                         "Removed: filtering (n=56)","Removed: contamination (n=19)",
                                         "Removed: spatial (n=100)","SPDE DE genes (n=3)"),
                        labels = c("Not signficant",
                                   "Removed: filtering\n(n=56)","Removed: contamination\n(n=19)",
                                   "Removed: spatial\n(n=100)","SPDE DE genes\n(n=3)"))) %>%
  arrange(label == "SPDE DE genes\n(n=3)") %>%
  ggplot(aes(x = fold_change, y = -log10(padj), col = label, pch = (label ==  "SPDE DE genes (n=3)"))) +
  geom_point(cex=1) +
  geom_hline(yintercept = -log10(0.05), lty = 2, cex = 0.3) +
  geom_label_repel(data=smide_naive_nb %>% filter(target %in% signSpatial) %>% mutate(label =  "SPDE DE genes\n(n=3)"),
                   aes(label = gene_symbol),show.legend = FALSE,size=3,
                   box.padding = 0.5, label.padding = 0.1, min.segment.length = 0.1,
                   force = 10)+
  theme_bw() +
  ggtitle("Naive DE genes (n = 178)")+
  labs(x = "Fold change", y = "-log10(pvalue)",col="")+
  scale_shape_manual(values = c(4,19),guide = "none") +
  scale_color_manual(values = c("lightgrey", my_colors[12],my_colors[2],my_colors[7],my_colors[8])) +
  theme(text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.key.height = unit(1.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3))) ;plot_volcano
ggsave(paste0(dir, "Images/Final/Volcano_Naive.jpeg"), height = 110, width = 200, units = "mm")


## Scatterplot ------------------------------------------------------------
library(cowplot)

# top DE in each method
results_all %>% filter(#Analysis == "Main",
  Distribution == "Negative Binomial",
  DE == 1) %>%
  mutate(pseudopvalue = ifelse(is.na(pseudopvalue),pvalue, pseudopvalue)) %>%
  ungroup %>%
  group_by(Method,meshscenario,kprop) %>%
  slice_min(order_by = pseudopvalue, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  select(Method,meshscenario,kprop, gene_symbol, Gene) %>%
  arrange(Method)%>%print(n=50)

results_compare <-  results_all %>% filter(Analysis == "Main",
                                           Distribution == "Negative Binomial") %>%
  filter(DEvar %in% c("kernelMicroglia1000")) %>%
  mutate(pseudopvalue = ifelse(is.na(pseudopvalue),ifelse(is.na(pvalue),padj,pvalue), pseudopvalue)) %>%
  ungroup %>%
  select(Gene,gene_symbol, Method,Estimate, pseudopvalue,DEvar ) %>%
  pivot_wider(names_from = Method, values_from = c(Estimate,pseudopvalue))

# SPDE vs naive
top_genes3 <- c("Whrn","Fosl2","Caln1","Zbtb16","Rasgrp1","Ndst4")
plot_grid(
  results_compare %>%
    ggplot(aes(x = `Estimate_SPDE` , y = `Estimate_smiDE Seminaive`
    )) +
    geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0,lty = 2)+
    geom_point(aes(col = gene_symbol %in% c("Fosl2","Zbtb16","Rasgrp1"))) + geom_abline(lty = 2) +
    geom_label_repel(
      data = subset(results_compare, gene_symbol %in% c(top_genes3)),
      aes(label = gene_symbol),size = 3,max.overlaps = Inf,
      box.padding = 0.5, label.padding = 0.1, min.segment.length = 0.1) +
    scale_color_manual(values = c("black","red")) +
    
    #facet_wrap(DEvar~.)+
    labs(x = "Estimate SPDE", y = "Estimate Seminaive") +
    theme_bw()+
    theme(legend.position = "none")
  ,
  results_compare %>%
    ggplot(aes(x = -log10(`pseudopvalue_SPDE`) , y = -log10(`pseudopvalue_smiDE Seminaive`))) +
    geom_hline(yintercept = -log10(0.05/258), lty=2) + geom_vline(xintercept = -log10(0.05/258), lty=2) +
    geom_point(aes(col = gene_symbol %in% c("Fosl2","Zbtb16","Rasgrp1"))) + geom_abline(lty=2) +
    geom_label_repel(
      data = subset(results_compare, gene_symbol %in% c(top_genes3)),
      aes(label = gene_symbol),size = 3,max.overlaps = Inf,
      box.padding = 0.5, label.padding = 0.1, min.segment.length = 0.1
    ) +
    scale_color_manual(values = c("black","red")) +
    labs(x = "-log10(pseudopvalue SPDE)", y = "-log10(pvalue smiDE Seminaive)")+
    theme_bw() +
    theme(legend.position = "none")
)
ggsave(paste0(dir, "Images/Final/Comp_inla_naive_distance.jpeg"), height = 110, width = 170, units = "mm")

# CPM plots spatial -------------------------------------------------------
load( paste0(dir, "Allen/WholeBrain_seurat_distance.Rdata"))
load( file = paste0(dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
colnames(seurat_allen@meta.data)[4] <- "cell_ID"

dataDEmetadata <- dataDE@meta.data
dataDEcounts <- as.matrix(dataDE@assays$RNA$counts)
dataDEnorm <- dataDE@assays$RNA$data
aux_selected<- genes %>% filter(gene_symbol %in% top_genes3)

for(id in 1:nrow(aux_selected)){
  selected_genes <- aux_selected$gene_identifier[id]
  label <- aux_selected$gene_symbol[id]
  
  
  library(ggpubr)
  dataDEmetadata %>% cbind(t(dataDEcounts)) %>%
    select(all_of(c("x","y","nCount_RNA",selected_genes))) %>%
    pivot_longer(all_of(selected_genes), names_to= "gene", values_to = "expr") %>%
    arrange(expr/nCount_RNA) %>%
    group_by(gene) %>% #mutate(expr = scale(expr/nCount_RNA )) %>%
    ggplot(aes( x= x, y = y, col = expr/nCount_RNA*1e6)) +
    geom_point(cex = 1) +
    scale_color_viridis_c()+
    facet_wrap(paste(label)~.) +
    labs(col = "CPM")+
    theme_bw() + coord_fixed()+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=17))
  ggsave(paste0(dir, "Images/Final/Expr_",label,".jpeg"), height = 90, width = 180, units = "mm")
}

# Benchmark ---------------------------------------------------------------
(results_all %>% 
   group_by(Method, Type, kprop, Distribution,DEvar, meshscenario,
            prior_range, prior_sigma, prior_phi, prior_prec) %>%
   summarise(nDE = sum(DE, na.rm = T),
             time = mean(Time, na.rm = T),
             PeakRAM = mean(PeakRAM, na.rm = T),
             TotalRAM = mean(TotalRAM)
   ) %>%
   filter(prior_range == "standard" |is.na(prior_range) ,
     prior_sigma == "standard"|is.na(prior_sigma),
     prior_phi == "balanced" |is.na(prior_phi) ,
     prior_prec == "standard" |is.na(prior_prec)) %>%
   mutate(label = ifelse(Method %in% c("SPDE"),paste("m =",meshscenario),
                         ifelse(Method %in% c("Independent\nclusters","GP","BYM2"),
                                paste("k =",kprop),"")))%>%
   ggplot(aes(x = time, y = PeakRAM, col = Method, group= paste(Method, Distribution),
              label = paste(meshscenario,kprop),
              pch = Distribution, lty = Distribution
   )) +
   geom_hline(yintercept=c(50,100,150,1500), col = "lightgrey")+
   geom_point(cex=2) +
   theme_bw()+
   geom_label_repel(aes(label = label),show.legend = FALSE, size = 3,max.overlaps = Inf,
                    box.padding = 0.5, label.padding = 0.1, min.segment.length = 0.1) +
   scale_x_continuous(trans=scales::log_trans(),
                      breaks = c(0, 60,60*2,60*5,60*10,60*30,60*60,2*60*60,4*60*60,16*60*60),
                      labels = c("0s","1m","2m","5m","10m","30m","1h","2h","4h","16h")) +
   scale_y_continuous(breaks=c(0,50,100,150,500,1000,1500))+
   scale_color_manual(limits = c("C-SIDE","nicheDE","smiDE Naive","smiDE Seminaive",
                                 "Independent\nclusters",
                                 "BYM2","GP","SPDE"),
                      values = my_colors[c(1,12,2,11,3,5,7,8)]) +
   labs(x = "Time", y= "PeakRAM (MiB)")+
   theme(text = element_text(size = 15)))
ggsave(paste0(dir, "Images/Final/Benchmark.jpeg"), height = 150, width = 300, units = "mm")

results_all%>%
  filter(Distribution != "Negative Binomial",
         meshscenario %in% c(0.5,2) | is.na(meshscenario),
         prior_range == "standard" |is.na(prior_range) , prior_sigma == "standard"|is.na(prior_sigma)) %>%
  dplyr::group_by(Method, Type, kprop,,meshscenario ) %>%
  summarise(nDEbonf = sum(DEbonf, na.rm = T),
            nDEfdr = sum(DEfdr, na.rm = T))



# Priors ------------------------------------------------------------------
# SPDE ------------------------------------------------------------------
avg_neurons <- overlap_metrics_allen$result %>%
  filter(class == "01 IT-ET Glut") %>%
  select(target, avg_cluster)

# run results_inla_wb on 4_Results_Distance.R
results_inla_wb %>% filter(failed == T)
results_inla_wb %>% filter(Distribution %in% c("negbin","gaussian"), meshscenario == 0.5,
                           !is.na(padj)) %>%
  dplyr::group_by(Distribution, prior_range, prior_sigma ) %>%
  summarise(nDE = sum(padj == 0, na.rm = T),
            n = n()) %>%
  mutate(label = paste(nDE,n,sep = "/")) %>%
  select(!c(nDE,n)) %>%
  pivot_wider(names_from = prior_sigma, values_from = label)

results_inla_wb %>% filter(Distribution == "negbin") %>%
  filter(padj == 0) %>%
  dplyr::group_by(meshscenario) %>%
  dplyr::count(Gene) %>%
  print(n=50)

results_inla_priors <- results_inla_wb %>% left_join(avg_neurons, by = c("Gene"="target")) %>%
  filter(Distribution %in% c("negbin","gaussian"),
         meshscenario %in% c(0.5)) %>%
  mutate(meshscenario = factor(meshscenario, levels = c(0.5,3),
                               labels = c("Dense mesh","Coarse mesh")),
         prior_range = factor(prior_range, levels = c("short","standard","long","weak"),
                              labels = c("Short","Standard","Long","Weak")),
         prior_sigma = factor(prior_sigma, levels = c("conservative","moderate","standard","permissive","flat"),
                              labels =c("Conservative","Moderate","Standard","Permissive","Flat")),
         Distribution = factor(Distribution, levels = c("negbin","gaussian"),
                               labels = c("Negative binomial","Gaussian")))


## Prior range -------------------------------------------------------------
(results_inla_priors %>%
   select(Gene,Distribution, meshscenario, prior_range, prior_sigma, pseudopvalue) %>%
   pivot_wider(names_from =prior_range, values_from = pseudopvalue  ) %>%
   pivot_longer(c(Long,Short,Weak), names_to = "prior_range", values_to = "pseudopvalue") %>%
   filter(prior_sigma == "Standard") %>%
   mutate(prior_range = factor(prior_range, levels = c("Short","Long","Weak"))) %>%
   ggplot(aes(x = -log10(Standard), y = -log10(pseudopvalue),label = Gene)) +
   geom_abline(slope = 1, intercept = 0, lty = 2, col = "darkgrey")+
   geom_point(cex=0.5) +
   facet_grid(Distribution ~ prior_range)) +
  theme_bw() +
  labs(x = "-log10(pseudopvalue Standard)",y = "-log10(pseudopvalue other)" )+
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/Final/SuppSPDE_PriorRange_densemesh.jpeg"),height = 120, width = 250,units = "mm")

# varying prior_range
results_inla_priors %>%
  select(Gene,Distribution, meshscenario, prior_range, prior_sigma, pseudopvalue,avg_cluster) %>%
  pivot_wider(names_from =prior_range, values_from = pseudopvalue  ) %>%
  pivot_longer(c(Long,Short,Weak), names_to = "prior_range", values_to = "pseudopvalue") %>%
  mutate(diff = Standard - pseudopvalue)%>%
  filter(prior_sigma == "Standard") %>%
  mutate(prior_range = factor(prior_range, levels = c("Short","Long","Weak"))) %>%
  ggplot(aes(x = avg_cluster, y = diff)) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex = 0.5) +
  facet_grid(Distribution ~ prior_range ) +
  theme(legend.position = "none") +
  labs(x = "Average expression",
       y= "Differences in pseudopvalues\n(Standard - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/Final/SuppSPDE_PriorRange_pvalueavgexpr_densemesh.jpeg"), height = 120, width = 250, units = "mm")


results_inla_priors%>%
  select(Gene,Distribution, meshscenario, prior_range, prior_sigma, Estimate,avg_cluster) %>%
  pivot_wider(names_from =prior_range, values_from = Estimate  ) %>%
  pivot_longer(c(Long,Short,Weak), names_to = "prior_range", values_to = "Estimate") %>%
  mutate(diff = Standard - Estimate)%>%
  filter(prior_sigma == "Standard") %>%
  mutate(prior_range = factor(prior_range, levels = c("Short","Long","Weak"))) %>%
  ggplot(aes(x = avg_cluster, y = diff)) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex=0.5) +
  facet_grid(Distribution ~ prior_range, scales = "free" ) +
  theme(legend.position = "none")+
  labs(x = "Average expression",
       y= "Differences in estimates\n(Standard - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/Final/SuppSPDE_PriorRange_estimateavgexpr_densemesh.jpeg"), height = 120, width = 250, units = "mm")

## Prior sigma -------------------------------------------------------------

(results_inla_priors %>%
   select(Gene,Distribution, meshscenario, prior_range, prior_sigma, pseudopvalue) %>%
   pivot_wider(names_from =prior_sigma, values_from = pseudopvalue  ) %>%
   pivot_longer(c(Conservative,Moderate,Permissive,Flat),
                names_to = "prior_sigma", values_to = "pseudopvalue") %>%
   filter(prior_range == "Standard") %>%
   mutate(prior_sigma = factor(prior_sigma, levels = c("Conservative","Moderate","Standard","Permissive","Flat"))) %>%
   ggplot(aes(x = -log10(Standard), y = -log10(pseudopvalue),label = Gene)) +
   geom_abline(slope=1,intercept=0, lty = 2, col = "darkgrey")+
   geom_point(cex=0.5) +
   facet_grid(Distribution ~ prior_sigma)) +
  theme_bw() +
  labs(x = "-log10(pseudopvalue Standard)",y = "-log10(pseudopvalue other)" )+
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/Final/SuppSPDE_PriorSigma_densemesh.jpeg"), height = 120, width = 250, units = "mm")


# varying prior_range
results_inla_priors %>%
  select(Gene,Distribution, meshscenario, prior_range, prior_sigma, pseudopvalue,avg_cluster) %>%
  pivot_wider(names_from =prior_sigma, values_from = pseudopvalue  ) %>%
  pivot_longer(c(Conservative,Moderate,Permissive,Flat),
               names_to = "prior_sigma", values_to = "pseudopvalue") %>%
  mutate(diff = Standard - pseudopvalue)%>%
  filter(prior_range == "Standard") %>%
  mutate(prior_sigma = factor(prior_sigma, levels = c("Conservative","Moderate","Standard","Permissive","Flat"))) %>%
  ggplot(aes(x = avg_cluster, y = diff)) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex = 0.5) +
  facet_grid(Distribution ~ prior_sigma,scales ="free" ) +
  theme(legend.position = "none") +
  labs(x = "Average expression",
       y= "Differences in pseudopvalues\n(Standard - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/Final/SuppSPDE_PriorSigma_pvalueavgexpr_densemesh.jpeg"), height = 120, width = 250, units = "mm")


results_inla_priors%>%
  select(Gene,Distribution, meshscenario, prior_range, prior_sigma, Estimate,avg_cluster) %>%
  pivot_wider(names_from =prior_sigma, values_from = Estimate  ) %>%
  pivot_longer(c(Conservative,Moderate,Permissive,Flat),
               names_to = "prior_sigma", values_to = "Estimate") %>%
  mutate(diff = Standard - Estimate)%>%
  filter(prior_range == "Standard") %>%
  mutate(prior_sigma = factor(prior_sigma, levels = c("Conservative","Moderate","Standard","Permissive","Flat"))) %>%
  ggplot(aes(x = avg_cluster, y = diff)) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex=0.5) +
  facet_grid(Distribution ~ prior_sigma,scales = "free" ) +
  theme(legend.position = "none")+
  labs(x = "Average expression",
       y= "Differences in estimates\n(Standard - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/Final/SuppSPDE_PriorSigma_estimateavgexpr_densemesh.jpeg"), height = 120, width = 250, units = "mm")

# BYM2 --------------------------------------------------------------------
avg_neurons <- overlap_metrics_allen$result %>%
  filter(class == "01 IT-ET Glut") %>%
  select(target, avg_cluster)

results_bym2_wb %>% filter(failed == T)
results_bym2_wb %>% filter(Distribution %in% c("negbin","gaussian"),
                           !is.na(padj)) %>%
  dplyr::group_by(Distribution, prior_phi, prior_prec, kprop ) %>%
  summarise(nDE = sum(padj == 0, na.rm = T),
            n = n()) %>%
  mutate(label = paste(nDE,n,sep = "/")) %>%
  select(!c(nDE,n)) %>%
  pivot_wider(names_from = prior_prec, values_from = label)


results_bym2_priors <- results_bym2_wb %>% left_join(avg_neurons, by = c("Gene"="target")) %>%
  filter(Distribution %in% c("negbin","gaussian"),
         kprop %in% c(0.250)) %>%
  mutate(prior_phi = factor(prior_phi, levels = c("iid","balanced","spatial"),
                            labels = c("IID","Balanced","Spatial")),
         prior_prec = factor(prior_prec, levels = c("conservative","moderate","standard","permissive","flat"),
                             labels =c("Conservative","Moderate","Standard","Permissive","Flat")),
         Distribution = factor(Distribution, levels = c("negbin","gaussian"),
                               labels = c("Negative binomial","Gaussian")))


## Prior phi -------------------------------------------------------------

(results_bym2_priors %>%
   select(Gene,Distribution, prior_phi, prior_prec, pseudopvalue) %>%
   pivot_wider(names_from =prior_phi, values_from = pseudopvalue  ) %>%
   pivot_longer(c(IID, Spatial), names_to = "prior_phi", values_to = "pseudopvalue") %>%
   filter(prior_prec == "Standard") %>%
   mutate(prior_phi = factor(prior_phi, levels = c("IID", "Spatial"))) %>%
   ggplot(aes(x = -log10(Balanced), y = -log10(pseudopvalue),label = Gene)) +
   geom_abline(slope = 1, intercept = 0, lty = 2, col = "darkgrey")+
   geom_point(cex=0.5) +
   facet_grid(Distribution ~ prior_phi)) +
  theme_bw() +
  labs(x = "-log10(pseudopvalue Balanced)",y = "-log10(pseudopvalue other)" )+
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/SuppBYM2_PriorPhi_kprop25.jpeg"),height = 120, width = 250,units = "mm")



# varying prior_phi
results_bym2_priors %>%
  select(Gene,Distribution, prior_phi, prior_prec, pseudopvalue,avg_cluster) %>%
  pivot_wider(names_from =prior_phi, values_from = pseudopvalue  ) %>%
  pivot_longer(c(IID, Spatial), names_to = "prior_phi", values_to = "pseudopvalue") %>%
  mutate(diff = Balanced - pseudopvalue)%>%
  filter(prior_prec == "Standard") %>%
  mutate(prior_phi = factor(prior_phi, levels = c("IID", "Spatial"))) %>%
  ggplot(aes(x = avg_cluster, y = diff
  )) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex = 0.5) +
  facet_grid(Distribution ~ prior_phi ,scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Average expression",
       y= "Differences in pseudopvalues\n(Balanced - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/SuppBYM2_PriorPhi_pvalueavgexpr_kprop25.jpeg"), height = 120, width = 250, units = "mm")


results_bym2_priors%>%
  select(Gene,Distribution, prior_phi, prior_prec, Estimate,avg_cluster) %>%
  pivot_wider(names_from =prior_phi, values_from = Estimate  ) %>%
  pivot_longer(c(IID, Spatial), names_to = "prior_phi", values_to = "Estimate") %>%
  mutate(diff = Balanced - Estimate)%>%
  filter(prior_prec == "Standard") %>%
  mutate(prior_phi = factor(prior_phi, levels = c("IID", "Spatial"))) %>%
  ggplot(aes(x = avg_cluster, y = diff)) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex=0.5) +
  facet_grid(Distribution ~ prior_phi, scales = "free" ) +
  theme(legend.position = "none")+
  labs(x = "Average expression",
       y= "Differences in estimates\n(Balanced - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/SuppBYM2_PriorPhi_estimateavgexpr_kprop25.jpeg"), height = 120, width = 250, units = "mm")

## Prior sigma -------------------------------------------------------------

(results_bym2_priors %>%
   select(Gene,Distribution, prior_phi, prior_prec, pseudopvalue) %>%
   pivot_wider(names_from =prior_prec, values_from = pseudopvalue  ) %>%
   pivot_longer(c(Conservative,Moderate,Permissive,Flat),
                names_to = "prior_prec", values_to = "pseudopvalue") %>%
   filter(prior_phi == "Balanced") %>%
   mutate(prior_prec = factor(prior_prec, levels = c("Conservative","Moderate","Standard","Permissive","Flat"))) %>%
   ggplot(aes(x = -log10(Standard), y = -log10(pseudopvalue),label = Gene)) +
   geom_abline(slope=1,intercept=0, lty = 2, col = "darkgrey")+
   geom_point(cex=0.5) +
   facet_grid(Distribution ~ prior_prec)) +
  theme_bw() +
  labs(x = "-log10(pseudopvalue Standard)",y = "-log10(pseudopvalue other)" )+
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/SuppBYM2_PriorPrec_kprop25.jpeg"), height = 120, width = 250, units = "mm")


# varying prior_phi
results_bym2_priors %>%
  select(Gene,Distribution, prior_phi, prior_prec, pseudopvalue,avg_cluster) %>%
  pivot_wider(names_from =prior_prec, values_from = pseudopvalue  ) %>%
  pivot_longer(c(Conservative,Moderate,Permissive,Flat),
               names_to = "prior_prec", values_to = "pseudopvalue") %>%
  mutate(diff = Standard - pseudopvalue)%>%
  filter(prior_phi == "Balanced") %>%
  mutate(prior_prec = factor(prior_prec, levels = c("Conservative","Moderate","Standard","Permissive","Flat"))) %>%
  ggplot(aes(x = avg_cluster, y = diff)) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex = 0.5) +
  facet_grid(Distribution ~ prior_prec,scales ="free" ) +
  theme(legend.position = "none") +
  labs(x = "Average expression",
       y= "Differences in pseudopvalues\n(Standard - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/SuppBYM2_PriorPrec_pvalueavgexpr_kprop25.jpeg"), height = 120, width = 250, units = "mm")


results_bym2_priors%>%
  select(Gene,Distribution, prior_phi, prior_prec, Estimate,avg_cluster) %>%
  pivot_wider(names_from =prior_prec, values_from = Estimate  ) %>%
  pivot_longer(c(Conservative,Moderate,Permissive,Flat),
               names_to = "prior_prec", values_to = "Estimate") %>%
  mutate(diff = Standard - Estimate)%>%
  filter(prior_phi == "Balanced") %>%
  mutate(prior_phi = factor(prior_phi, levels = c("Conservative","Moderate","Standard","Permissive","Flat"))) %>%
  ggplot(aes(x = avg_cluster, y = diff)) +
  geom_hline(yintercept = 0, lty = 2 , col = "darkgrey") +
  geom_point(cex=0.5) +
  facet_grid(Distribution ~ prior_prec,scales = "free" ) +
  theme(legend.position = "none")+
  labs(x = "Average expression",
       y= "Differences in estimates\n(Standard - other)") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave(paste0(dir, "Images/SuppBYM2_PriorPrec_estimateavgexpr_kprop25.jpeg"), height = 120, width = 250, units = "mm")




