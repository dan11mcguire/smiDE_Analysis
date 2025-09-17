#------------------------#
#-- Results Simulation --#
#------------------------#

dir.out = "~/Simulation/"
library(tidyverse)
library(cowplot)

dataplot <- read.csv(paste0(dir.out, "data/ResultsSimulation.csv"))
fulldata <- read.csv(paste0(dir.out,"data/Lung5-3data.csv"))

# Figures -----------------------------------------------------------------
fulldata <- fulldata %>% 
  filter(cell_type == "macrophage",
         niche %in% c("myeloid-enriched stroma","stroma"))

my_colors <- c( "#5F4690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05",
                "#CC503E", "#94346E", "#6F4070", "#994E95","#666666")

# Different clusters ------------------------------------------------------
fulldata$clusterexample <- kmeans(fulldata[,c('sdimx','sdimy')], 59, iter.max = 50)$cluster
plotclusters1<- fulldata %>% mutate(niche = factor(niche,levels = c("myeloid-enriched stroma","stroma"),
                                                   labels = c("Myeloid-enriched stroma","Stroma"))) %>%
  ggplot(aes(x = sdimx, y = sdimy, col = as.factor(clusterexample))) +
  geom_point( cex = 0.4) +
  scale_color_manual(values = colorRampPalette(my_colors)(59))  +
  coord_fixed() + theme_bw() +
  labs(col = "") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  ggtitle("58 clusters (K/n=1%)")+
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 6),
        legend.box.margin=margin(-10,-10,-10,-10),plot.title = element_text(hjust = 0.5, size =7))

fulldata$clusterexample2 <- kmeans(fulldata[,c('sdimx','sdimy')], 294 , iter.max = 50)$cluster
plotclusters2<- fulldata %>% mutate(niche = factor(niche,levels = c("myeloid-enriched stroma","stroma"),
                                                   labels = c("Myeloid-enriched stroma","Stroma"))) %>%
  ggplot(aes(x = sdimx, y = sdimy, col = as.factor(clusterexample2))) +
  geom_point( cex = 0.4) +
  scale_color_manual(values = colorRampPalette(my_colors)(294))  +
  coord_fixed() + theme_bw() +
  labs(col = "") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  ggtitle("294 clusters (K/n=5%)")+
  theme(legend.position = "none",
        text = element_text(size = 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 6),
        legend.box.margin=margin(-10,-10,-10,-10),plot.title = element_text(hjust = 0.5, size = 7))

# Different mesh ------------------------------------------------------
library(fmesher)
# save mesh and add with photoshop
plotmeshdense<- fulldata %>%
  ggplot(aes(x = sdimx, y = sdimy, col = as.factor(clusterexample))) +
  theme_void() +
  ggtitle("Dense mesh (m=0.5)")+
  theme(plot.title = element_text(hjust = 0.5, size = 7))

plotmeshcoarse<- fulldata %>%
  ggplot(aes(x = sdimx, y = sdimy, col = as.factor(clusterexample))) +
  theme_void() +
  ggtitle("Coarse mesh (m=0.5)")+
  theme(plot.title = element_text(hjust = 0.5, size = 7))

range <- max(max(fulldata$sdimx) - min(fulldata$sdimx), max(fulldata$sdimy) - min(fulldata$sdimy))
cutoff <- range/sqrt(nrow(fulldata))
mesh0.5 <- fm_mesh_2d(loc = fulldata[, c("sdimx","sdimy")],
                      offset = c(0.01 * range,  0.05 * range),
                      cutoff = cutoff * 0.5
)
mesh2 <- fm_mesh_2d(loc = fulldata[, c("sdimx","sdimy")],
                    offset = c(0.01 * range,  0.05 * range),
                    cutoff = cutoff * 3
)
png(paste0(dir.out,"Images/Densemesh.jpeg"), width = 800, height = 600)   # open graphics device
plot(mesh0.5)
dev.off()

png(paste0(dir.out,"Images/Coarsemesh.jpeg"), width = 800, height = 600)   # open graphics device
plot(mesh2)
dev.off()


# Macrophages -------------------------------------------------------------
macrophages <- fulldata %>% mutate(niche = factor(niche,levels = c("myeloid-enriched stroma","stroma"),
                                                  labels = c("Myeloid-enriched\nstroma","Stroma"))) %>%
  ggplot(aes(x = sdimx, y = sdimy, col = niche)) +
  geom_point( cex = 0.4) +
  scale_color_manual(values = c("dark grey", my_colors[1]))  +
  coord_fixed() + theme_bw() +
  labs(col = "") +
  guides(colour = guide_legend(override.aes = list(size=2),
                               keywidth     = unit(2, "mm"),
                               keyheight    = unit(2, "mm"))) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 6),
        plot.margin       = margin(3, 20, 3, 3),  
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin     = margin(0, 0, 0, 0),
        legend.key.width     = unit(2, "mm"),
        legend.key.height    = unit(2, "mm"),
        legend.key           = element_blank(),
        legend.spacing.x     = unit(10, "mm"),
        legend.spacing.y     = unit(1, "mm"))


# Examples ----------------------------------------------------------------
plotexample <- function(sim,cexpoints=0.5, labels){
  coefs <- read.csv(paste0(dir.out, "Examples/coefs",sim,".csv"))
  dataranef <- read.csv(paste0(dir.out, "Examples/dataranef",sim,".csv"))
  
  p1 <- dataranef %>% arrange(Ycount) %>%
    ggplot(aes(x = sdimx, y = sdimy, col = Ycount)) +
    geom_point(alpha = 0.8,cex=cexpoints) +
    theme_bw()+
    labs(col = "")+
    scale_color_viridis_c()+
    theme_light() +
    coord_fixed()+
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    theme(legend.position = "top",
          legend.text = element_text(size = 6),
          legend.box.margin=margin(-10,-10,-10,-10)) +
    guides(col = guide_colorbar(title.vjust = 0.7,
                                barheight = 0.5,
                                barwidth = 3))
  
  
  p2 <- coefs %>%
    filter(Model %in% c("SPDEdense","Naive","IndepClust5","GP25","BYM2_25")) %>%
    mutate(Model = factor(Model, levels = c("Naive","IndepClust5","GP25","BYM2_25","SPDEdense"),
                          labels = c("Naive","Independent\nClusters","GP","BYM2","SPDE-INLA"))) %>%
    ggplot(aes(x = fct_rev(Model), y = exp(Estimate), col = Model)) +
    geom_hline(yintercept = 1) +
    geom_point(position = position_dodge(width=0.9)) +
    geom_errorbar(aes(ymin = exp(LB), ymax = exp(UB)), position = position_dodge2(width=0.9)) +
    scale_color_manual(values = my_colors[c(12,3,5,7,8)]) +
    theme_bw()+
    labs(x = "", y = "Fold change")+
    coord_flip()+
    theme(legend.position = "none",
          text = element_text(size = 8)) +
    theme(axis.text.x=element_text(angle = -90, hjust = 0));p2
  
  
  #-2.42  3.04
  ranef <- dataranef %>% 
    pivot_longer(c(glmres:bym250),
                 names_to = "Model",
                 values_to = "ranef") %>%
    filter(Model %in% c("glmres","indep5","gp25","bym225","spdedense")) %>%
    mutate(Model = factor(Model, levels = c("glmres","indep5","gp25","bym225","spdedense"),
                          labels = c("Naive","Independent\nClusters","GP","BYM2","SPDE-INLA"))) %>%
    filter(!Model == "Naive") %>%
    arrange(abs(ranef)) %>%
    ggplot(aes(x = sdimx, y = sdimy, col = ranef)) +
    geom_point(alpha = 0.8,cex=0.03) +
    facet_grid(.~Model)+
    scale_colour_gradient2(low = my_colors[3], high = my_colors[7],
                           mid = "light grey",
    ) +
    theme_bw()+
    coord_fixed()+
    labs(col = "Random effects")+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 6, angle = -90),
          legend.text = element_text(size = 5),
          legend.box.margin=margin(-10,-10,-10,-10))+
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) +
    theme(strip.background =element_blank(),
          strip.text = element_text(colour = "black", size = 8))+
    guides(col = guide_colorbar(title.vjust = 0,
                                barwidth = 0.5,
                                barheight = 2))
  
  cowplot::plot_grid(p1,p2,ranef,nrow=1,rel_widths = c(1/7,1.2/7,3/7), labels = labels)
}
examples <- cowplot::plot_grid(plotexample(1, labels = c("d","f","h")),
                               plotexample(99, labels = c("e","g","i")),
                               ncol =1)


# Results -----------------------------------------------------------------
library(scales)
library(ggrepel)
expm1_trans <-  function() trans_new("expm1", "expm1", "log1p")
plotnew_t1err <- dataplot %>%
  filter((prior_range == "standard" | is.na(prior_range)),
         (prior_sigma == "standard" | is.na(prior_sigma)),
         (prior_phi == "balanced" | is.na(prior_phi)),
         (prior_prec == "standard" | is.na(prior_prec))) %>%
  filter(prop %in% c(0.01,0.05,0.10,0.10,0.30,0.5,1) | is.na(prop)) %>%
  ggplot(aes(x = (Time), y = Type1Error, col = Model,
             group = paste(Model,Distribution),
             label = prop)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.05) +
  geom_label_repel(data = dataplot %>%
                     filter((prop %in% c(0.05,0.3) &Model %in% c("Independent\nclusters")) |
                              (prop %in% c(1) & Model %in% c("BYM2")) |
                              (prop %in% c(0.01,0.3,1)& Model %in% c("GP")) |
                              (Model %in% c("SPDE-INLA"))
                     ),
                   aes(label = ifelse(Model == "SPDE-INLA",ifelse(mesh == 0.5, "Dense","Coarse"),paste0(prop*100,"%"))
                   ),show.legend = FALSE, size = 2.5, max.overlaps = Inf,
                   box.padding = 0.5, label.padding = 0.1, min.segment.length = 0.1)+
  scale_x_continuous(trans=log_trans(), breaks = c(0,1, 10, 60,60*5,60*20, 60*60,7*60*60),
                     labels = c("0s","1s","10s","1m","5m","20m","1h","7h")) +
  scale_y_continuous(breaks = c(0,0.05,0.25,0.5,0.75,1))+
  labs(x = "Time", y = "Type 1 Error", col ="", lty = "", pch = "")+
  scale_color_manual(limits = c("Wilcox","MAST","NB GLM","DESeq2","C-SIDE",
                                "Independent\nclusters","BYM2","GP","SPDE-INLA"),
                     values = c("black","#666666","#264653","#1D6996","#5F4690",
                                "#38A6A5","#73AF48","#E17C05","#CC503E"))+
  scale_shape_manual(values = c( 3,19)) +
  theme_bw() +
  facet_wrap(Distribution~.)+
  theme(legend.position = "top",
        strip.background = element_blank(),
        text = element_text(size = 9),
        strip.text = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.margin = margin(-0.1,0,0,0, unit="cm"),
        legend.key.width = unit(6,"cm"))+
  guides(colour = guide_legend(nrow = 1,
                               keywidth=0.15,
                               keyheight=0.1,
                               default.unit="inch",
                               byrow = T),
         pch = guide_legend(nrow = 1,
                            keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch",
                            byrow = T),
         cex = guide_legend(nrow = 1,
                            keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch",
                            byrow = T),
         lty = guide_legend(nrow = 1,
                            keywidth=0.5,
                            keyheight=0.1,
                            default.unit="inch",
                            byrow = T));plotnew_t1err
plotnew_rank <- dataplot %>% 
  filter((prior_range == "standard" | is.na(prior_range)),
         (prior_sigma == "standard" | is.na(prior_sigma)),
         (prior_phi == "balanced" | is.na(prior_phi)),
         (prior_prec == "standard" | is.na(prior_prec))) %>%
  filter(prop %in% c(0.01,0.05,0.10,0.10,0.30,0.5,1) | is.na(prop)) %>%
  ggplot(aes(x = (Time), y = Rank, col = Model,
             group = paste(Model,Distribution),
             label = prop)) +
  geom_point() + geom_line() +
  geom_label_repel(data = dataplot %>%
                     filter((prop %in% c(0.05,0.3) &Model %in% c("Independent\nclusters")) |
                              (prop %in% c(1) & Model %in% c("BYM2")) |
                              (prop %in% c(0.01,0.3,1)& Model %in% c("GP")) |
                              (Model %in% c("SPDE-INLA"))
                     ),
                   aes(label = ifelse(Model == "SPDE-INLA",ifelse(mesh == 0.5, "Dense","Coarse"),paste0(prop*100,"%"))
                   ),show.legend = FALSE, size = 2.5, max.overlaps = Inf,
                   box.padding = 0.5, label.padding = 0.1, min.segment.length = 0.1)+
  scale_x_continuous(trans=log_trans(), breaks = c(0,1, 10, 60,60*5,60*20, 60*60,7*60*60),
                     labels = c("0s","1s","10s","1m","5m","20m","1h","7h")) +
  scale_y_continuous(limits= c(0.75, 1))+
  facet_wrap(Distribution~.)+
  labs(x = "Time", y = "Rank", col ="", lty = "", pch = "")+
  scale_color_manual(limits = c("Wilcox","MAST","NB GLM","DESeq2","C-SIDE",
                                "Independent\nclusters","BYM2","GP","SPDE-INLA"),
                     values = c("black","#666666","#264653","#1D6996","#5F4690",
                                "#38A6A5","#73AF48","#E17C05","#CC503E"))+
  scale_shape_manual(values = c( 3,19)) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 9),
        strip.text = element_text(size = 9));plotnew_rank
plot_results<-cowplot::plot_grid(legend = cowplot::get_plot_component(plotnew_t1err, 'guide-box-top', return_all = TRUE),
                                 cowplot::plot_grid(plotnew_t1err + theme(legend.position = "none"),
                                                    plotnew_rank + theme(legend.position = "none"),
                                                    nrow = 1 ,labels = c("j","k")
                                 ),
                                 nrow = 2,rel_heights = c(0.15,0.85))
# Plot --------------------------------------------------------------------

plot_grid(plot_grid(plotclusters1,plotclusters2, plotmeshdense,plotmeshcoarse,macrophages,nrow=1,
                    labels = c("a","","b","", "c"),
                    align = "v",axis = "l",
                    rel_widths = c(1,1,1,1,1.5)),
          examples,
          plot_results,
          ncol = 1, rel_heights = c(1.2,2.5,3))
ggsave(paste0(dir.out,"Images/Figuresimulations_revised.jpeg"), width = 250, height = 200, units = "mm")

# Supplement --------------------------------------------------------------
## Power -------
datapower <- read.csv(paste0(dir.out, "data/ResultsSimulation_power.csv"))
datapower %>% filter(beta <= 0.5,beta > 0,
  prop %in% c(0.05,0.3) | is.na(prop)) %>%
  filter(Model %in% c("BYM2","SPDE-INLA","GP","Independent\nclusters")) %>%
  mutate(label = ifelse(Model == "SPDE-INLA", paste0("m: ",mesh), paste0("K/n:",prop))) %>%
  ggplot(aes(x = beta, y = Power, col = as.factor(label),label = prop, lty =Distribution,
             group=paste(Model,mesh,prop,Distribution))) +
  geom_smooth(se=F,span = 0.5, cex=.8 )+
  geom_hline(yintercept = 0.05)+
  scale_linetype_manual(values = c(2,1))+
  facet_grid(.~Model) +
  labs(x = "Beta", y = "Power", col = "") +
  theme_bw() +
  guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm"))
ggsave(paste0(dir.out,"Images/Figure_power.jpeg"), width = 250, height = 120, units = "mm")

# Priors ------------------------------------------------------------------
# SPDE
dataplot %>% filter(mesh == 0.5,
                    Model == "SPDE-INLA") %>%
  select(Distribution, prior_range, prior_sigma, Type1Error) %>%
  pivot_wider(names_from = prior_sigma, values_from = Type1Error)

# BYM2
dataplot %>% filter(prop == 1,
                    Model == "BYM2") %>%
  select(Distribution, prior_phi, prior_prec, Type1Error) %>%
  pivot_wider(names_from = prior_prec, values_from = Type1Error)
