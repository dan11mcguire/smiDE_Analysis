# devtools::install_github("Winnie09/Palo")
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
setwd("/home/rstudio/data/dan/smiDE_package/smiDE_Analysis/data/")
orm_metrics <- fread("colon_de_results/orm_metrics_colon.csv")
load("coloncancer.RData")

RhpcBLASctl::blas_set_num_threads(1)
metainfo <- data.table(cbind(annot, clust))

#targset <- as.numeric(args[1])

### pre de
pre_de_bcell <- 
  smiDE::pre_de(counts = Matrix::t(raw)
                ,normalized_data = Matrix::t(raw)
                ,weight_colname = NULL
                ,metadata = metainfo
                ,ref_celltype = "B-cell"
                ,cell_type_metadata_colname = "clust"
                ,sdimx_colname = "x_slide_mm"
                ,sdimy_colname = "y_slide_mm"
                ,mm_radius = 0.015 ### 15 microns
                ,aggregation = "sum"
                ,adjacencies_only = FALSE)

#pre_de_bcell$nblist$adjacency_counts_by_ct[1:10] ### count # of cells of each 
#                                                 ### neighbor cell type

bcell_cellid <- metainfo[clust=="B-cell", cell_ID]
pre_de_bcell$nblist$adjacency_counts_by_ct[match(bcell_cellid
                                                 ,cell_ID)
                                           ,summary(rowSums(.SD))
                                           ,.SDcols=setdiff(colnames(pre_de_bcell$nblist$adjacency_counts_by_ct), "cell_ID")]

immune_types <- c("B-cell", "macrophage", "mDC", "neutrophil"
                  ,"pDC", "mast", "monocyte", "NK", "plasmablast"
                  ,"T CD8 memory", "plasma.cell", "plasmablast"
                  ,"T CD4 naive", "Treg")

setdiff(metainfo[,unique(clust)], immune_types)

celltype_nbrs <- 
  copy(pre_de_bcell$nblist$adjacency_counts_by_ct)[,immune_ct_neighbors:=rowSums(.SD)
                                                   ,.SDcols=c(immune_types)
  ]

metainfo_bcell <- merge(metainfo[clust=="B-cell"]
                        ,celltype_nbrs[,.(cell_ID, immune_ct_neighbors)]
                        , by="cell_ID")

ggplot(metainfo, aes(x_slide_mm, y_slide_mm, color=factor(fov))) + 
  geom_point(size=0.1) + 
  geom_point(data=metainfo[clust=="B-cell"], size=0.1) + 
  coord_fixed() + 
  theme_bw() + 
  scale_color_manual(values=rep(unname(pals::alphabet()),3)
#  scale_color_manual(values=ctpal[levs]
 #                    ,breaks=levs
                     ,guide=guide_legend(override.aes=list(size=4))) + 
  geom_text(data=metainfo[,lapply(.SD, mean),.SDcols=c("x_slide_mm", "y_slide_mm"), by=fov]
            ,aes(x=x_slide_mm, y=y_slide_mm, label=fov), inherit.aes = FALSE, color='black') 


metainfo[,tissgrp:="ll"]
metainfo[fov >=60, tissgrp:="ur"]
metainfo[,bcell_nbr:=0]
metainfo[match(pre_de_bcell$cell_adjacency_dt[from %in% metainfo_bcell[["cell_ID"]] | 
                                                to %in% metainfo_bcell[["cell_ID"]]
                                              ,unique(c(from, to))]
               ,cell_ID), bcell_nbr:=1]
metainfo[,.N,by=.(bcell_nbr)]

metainfo[,ct:="other"]
metainfo[clust %in% immune_types,ct:="other immune cell type"]
metainfo[clust == "B-cell",ct:="B cell"]
metainfo[,ct:=factor(ct, levels=c("other", "other immune cell type", "B cell"))]
ctpal2 <- c("grey50", "#33A02C", "black")
names(ctpal2) <- c("other", "other immune cell type", "B cell")


plist <- 
lapply(split(metainfo, by="tissgrp"), function(xx){
  ggplot(xx[clust!="B-cell"][order(ct)], aes(x_slide_mm, y_slide_mm, color=ct)) + 
    geom_point(size=0.05, alpha=0.5) + 
    geom_point(data=xx[clust=="B-cell"], size=0.3) + 
    coord_fixed() + 
    theme_bw() + 
    scale_color_manual(values=rev(ctpal2)
                       ,breaks=rev(names(ctpal2))
                       ,name="cell type"
                       ,guide=guide_legend(override.aes=list(size=4, alpha=1))) + 
    theme(legend.position="bottom")
    
})
cowplot::plot_grid(plist[[1]] + theme(plot.margin = unit(c(0.1,-2,0,0), "cm")) + theme(legend.position="none")
                   ,plist[[2]] + theme(plot.margin = unit(c(0,0.1,0,-2), "cm")) + 
                     theme(legend.position = "none")
                   ,rel_widths = c(8,4))


pll <- plist[[1]] + theme(plot.margin = unit(c(0.1,-2,0,0), "cm")) + theme(legend.position="none")
pll2 <- plist[[2]] + theme(plot.margin = unit(c(0,0.1,0,-2), "cm")) + theme(legend.position = "none")
pcelltype <- 
cowplot::plot_grid(pll + theme(plot.margin = unit(c(0,-2,0,0), "cm"))
                   ,cowplot::plot_grid(pll2,NULL, rel_heights = c(3, 4.2), nrow=2) + 
                      theme(plot.margin = unit(c(0,0.1,0,-2), "cm"))
                   ,rel_widths = c(3, 2.5))

pcelltype

metainfo_bcell[,.N,by=.(immune_ct_neighbors)]
metainfo_bcell[,imm_cut:=cut(immune_ct_neighbors
                       ,quantile(immune_ct_neighbors, seq(0,1,1/7))
                       ,include.lowest=TRUE)]

metainfo_bcell[,.N,by=.(imm_cut)]
metainfo_bcell[,tissgrp:="ll"]
metainfo_bcell[fov >=60, tissgrp:="ur"]
cutpal <- colorRampPalette(brewer.pal(11, "RdYlBu"))(7)
names(cutpal) <- metainfo_bcell[,levels(imm_cut)]
plist_immcut <- 
lapply(split(metainfo_bcell, by="tissgrp"), function(xx){
  ggplot(xx[order(imm_cut)]
         ,aes(x_slide_mm, y_slide_mm, fill=imm_cut)) + 
    geom_point(size=2
               ,pch=21
               ,alpha=1
               ) + 
    coord_fixed() + 
    theme_bw() + 
    scale_fill_manual(name = "# of immune cell neighbors"
                      ,values=cutpal, guide=guide_legend(
                                                         override.aes=list(size=4,alpha=1))) + 
    theme(legend.position="right")
})

rightside3 <- cowplot::plot_grid(plist_immcut$ll + theme(legend.position="none") +
                                   theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title = element_blank()) + 
                                         theme(plot.margin = unit(c(0.1,0,0,0), "cm"))
                                       ,cowplot::plot_grid(plist_immcut$ur + theme(legend.position="none") +
                                         theme(axis.text = element_blank(), axis.ticks=element_blank() ,axis.title = element_blank()) + 
                                         theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, rel_widths = c(2.5,0.1))
                                       , nrow=2
                                       ,rel_heights = c(6, 3.3)
                                ,align = 'v'
          )
rightside3


leftside <- 
cowplot::plot_grid(cowplot::plot_grid(plist$ll + theme(legend.position="none") 
                                      ,plist$ur + theme(legend.position="none")
                                      ,nrow=2, rel_heights = c(6,3))
                   )

cowplot::plot_grid(leftside, rightside3, nrow=1)

leg1 <- cowplot::get_plot_component(plist[[1]] + theme(legend.position = "right"), 'guide-box-right', return_all=TRUE)
leg2 <- cowplot::get_plot_component(plist_immcut[[1]] , 'guide-box-right', return_all=TRUE)


combleg <- cowplot::plot_grid(cowplot::ggdraw(leg1) + theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              ,cowplot::ggdraw(leg2)+ theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              , nrow=2)

plot_ab <- 
cowplot::plot_grid(leftside
                   ,combleg
                   ,rightside3
                   ,nrow=1
                   ,labels = c("a", "", "b")
                   ,rel_widths = c(1, 0.4, 1)) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
print(plot_ab)
#ggsave("colon_data/deplot1.png"
#       ,width=938*sqrt(25)
#       ,height=904*sqrt(25)
#       ,dpi=450
#       ,units = "px")
#
dev1f <- list.files("colon_de_results", pattern="naive_de",full.names=TRUE)
dev2f <- list.files("colon_de_results", pattern="de_seminaive",full.names=TRUE)
dev3f <- list.files("colon_de_results", pattern="de_spatial",full.names=TRUE)

deobj <- readRDS(dev2f[1])
for(ff in 2:length(dev2f)){
  message(ff / length(ff))
  deobjff <- readRDS(dev2f[ff])
  deobj$results <- c(deobj$results, deobjff$results)
  deobj$targets <- c(deobj$target, deobjff$targets)
}
saveRDS(deobj, file="colon_data/seminaive_de_results.rds")

deobj <- readRDS(dev1f[1])
for(ff in 2:length(dev1f)){
  message(ff / length(ff))
  deobjff <- readRDS(dev1f[ff])
  deobj$results <- c(deobj$results, deobjff$results)
  deobj$targets <- c(deobj$target, deobjff$targets)
}
saveRDS(deobj, file="colon_data/naive_de_results.rds")

deobj <- readRDS(dev3f[1])
for(ff in 2:length(dev3f)){
  message(ff / length(ff))
  deobjff <- readRDS(dev3f[ff])
  deobj$results <- c(deobj$results, deobjff$results)
  deobj$targets <- c(deobj$target, deobjff$targets)
}
saveRDS(deobj, file="colon_data/inlaspatialaware_de_results.rds")



deobj1 <- readRDS("/home/rstudio/data/dan/nsclc_de_analysis/colon_data/naive_de_results.rds")
deobj2 <- readRDS("/home/rstudio/data/dan/nsclc_de_analysis/colon_data/seminaive_de_results.rds")
deobj3 <- readRDS("/home/rstudio/data/dan/nsclc_de_analysis/colon_data/inlaspatialaware_de_results.rds")

ms1 <- smiDE::results(deobj1, "model_summary")[[1]]
pw1 <- smiDE::results(deobj1, "pairwise")[[1]]
emm1 <-smiDE::results(deobj1, "emmeans")[[1]]

ms2 <- smiDE::results(deobj2, "model_summary")[[1]]
pw2 <- smiDE::results(deobj2, "pairwise")[[1]]
emm2 <-smiDE::results(deobj2, "emmeans")[[1]]

ms3 <- smiDE::results(deobj3, "model_summary")[[1]]
pw3 <- smiDE::results(deobj3, "pairwise")[[1]]
emm3 <-smiDE::results(deobj3, "emmeans")[[1]]
sre <- smiDE::results(deobj3, "spatial_random_effect")[[1]]

pw1p <- 
merge(pw1[term=="immune_ct_neighbors"][!grepl("Custom", target)]
      ,orm_metrics[clust=="B-cell",.(target, contam_ratio=ratio)]
      ,by="target"
      )

pw2p <- 
merge(pw1p, pw2[term=="immune_ct_neighbors"], by=c("target"), suffixes=c("_v1", "_v2"))

pw3p <- merge(pw2p, pw3[term=="immune_ct_neighbors"], by="target")

pw3p <- 
merge(pw3p, ms3[term=="Stdev for s"][,.(spatial_re_sd_mean = mean,target)]
      ,by="target")

pw3p[,lcl_bonf:=`9.255831e-06quant`]
pw3p[,ucl_bonf:=`0.9999907quant`]
pw3p[,lcl:=`0.025quant`]
pw3p[,ucl:=`0.975quant`]

pw3p[,sig_bonf:=as.numeric(sign(log(lcl_bonf))==sign(log(ucl_bonf)))]
pw3p[,sig:=as.numeric(sign(log(lcl))==sign(log(ucl)))]
pw3p[,showbound:=`0.5quant`]
pw3p[sig==1 & `0.975quant` > 1,showbound:=`0.025quant`]
pw3p[sig==1 & `0.975quant` < 1,showbound:=`0.975quant`]
pw3p[`0.975quant` < 1, showbound:=`0.975quant`]
pw3p[,sigcat:="nonsig"]
pw3p[sig==1,sigcat:="significant at 95% CI"]
pw3p[sig_bonf==1,sigcat:="significant at 95% Bonferroni adjusted CI"]

pw3p[,.(sig*abs(`0.5quant`-1), `0.5quant`, log2(`0.5quant`), sig*abs(log2(`0.5quant`)-0))]
pw3p[,yval1:=sig*abs(log2(`0.5quant`)-0)]
pw3p[,yval2:=sig*abs(`0.5quant`-1)]
pw3p[,yval3:=sig*abs(`0.5quant`*(14.54-9.458)-1)]

sigpal <- c("grey", brewer.pal(9,"Pastel1")[2],scales::muted("red"))
names(sigpal) <- c("nonsig", "significant at 95% CI", "significant at 95% Bonferroni adjusted CI")
pw3p[,version:="spatial-aware (INLA)"]
vinla <- 
ggplot(pw3p, aes(x=log2(showbound)
                 ,y = yval2
                 ,color=sigcat)) + 
  geom_point() + 
  theme_bw(base_size=20) + 
  scale_color_manual(values=sigpal
                       ,guide=guide_legend(override.aes=list(size=4),nrow=2)
                     ,breaks=names(sigpal)[2:3]
                     ,name="Significance") + 
  geom_label_repel(data=pw3p[sigcat=="significant at 95% Bonferroni adjusted CI"]
                   ,aes(label=target)
                   ,show.legend=FALSE
                   ) + 
  labs(x="log2 (Fold Change Credible Interval Bound closest to 1)"
       ,y="Fold Change (+ Positive, - Negative)") + 
  scale_y_continuous(breaks=seq(0,0.8, 0.1)
                     ,labels=paste0("+ ", c(1 + seq(0,0.8, 0.1)), "\n-", c(1 - seq(0,0.8, 0.1)))) + 
  facet_wrap(~version, labeller=label_both)

vinla

pw12 <- 
rbindlist(list(pw1[,version:="naive"]
               ,pw2[,version:="semi-naive"])
          ,use.names=TRUE,fill=TRUE)

v12 <- 
smiDE::volcano(pw12[term=="immune_ct_neighbors"][!grepl("Custom", target)]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]] + 
  facet_wrap(~version,labeller=label_both)

v12 + labs()
g1 = ggplotGrob(v12)
gtable_show_layout(g1)


xaxis <- g1[c(12),]
ggplotify::as.ggplot(xaxis)

g1_am0 <- g1[,c(5:9)]
g1p <- ggplotify::as.ggplot(g1_am0)
g2_am0 <- g1[,c(9:12)]

g2p <- ggplotify::as.ggplot(g2_am0)

gab <- 
cowplot::plot_grid(g1p  #+ labs(title="Naive model with no control variable adjustment or filtering") + theme(plot.title=element_text(hjust=0.5, vjust=-8, size=20))
                   ,g2p #+ labs(title="Semi-Naive")
                   ,rel_widths = c(1,0.8)
                   ,labels=c("a", "b")
                   ,label_size = 20)
gab
gabp <- 
cowplot::plot_grid(gab + theme(plot.margin = unit(c(0,0,-1,0), "cm"))
                   ,xaxis #+ theme(plot.margin = unit(c(0,0,1,0), "cm"))
                   , rel_heights = c(1, 0.05), nrow=2)



cowplot::plot_grid(gabp
                   ,cowplot::plot_grid(NULL
                                       ,vinla + theme(legend.position = "bottom")
                                       ,NULL
                                       ,rel_heights = c(0.2, 1, 0.1)
                                       ,nrow=3
                                       )
                   ,nrow=1, rel_widths = c(2, 0.8))

gabcp <- 
cowplot::plot_grid(gabp + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
                   ,vinla + theme(legend.position = "top"
                                  ,axis.title.x = element_text(size=16)
                                  )
                   ,labels=c("", "c"), label_size = 20
                   ,nrow=1, rel_widths = c(2, 1))

gabcp + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)

#ggsave("colon_data/deplot2.png"
#       ,dpi=450
#       ,width=7.54 * sqrt(80)#sqrt(100)
#       ,height=700 * sqrt(80)#sqrt(100)
#       ,units = "px"
#       )
#ggsave("colon_data/deplot22.png"
#       ,dpi=450
#       ,width=19*sqrt(15)
#       ,height=8*2.54/3 *sqrt(15)
#       ,units = "cm"
#       )



pmiddle <- 
cowplot::plot_grid(g1p  #+ labs(title="Naive")
                   ,g2p #+ labs(title="Semi-Naive")
                   ,rel_widths = c(1,0.8)
                   ,labels=c("b", "c")
                   ,label_size = 20)

pmiddle <- 
cowplot::plot_grid(pmiddle + theme(plot.margin = unit(c(0,0,-1,0), "cm"))
                   ,xaxis #+ theme(plot.margin = unit(c(0,0,1,0), "cm"))
                   , rel_heights = c(1, 0.05), nrow=2)


pmiddle <- 
cowplot::plot_grid(pmiddle + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
                   ,vinla + theme(legend.position = "top"
                                  ,axis.title.x = element_text(size=14))
                   ,labels=c("", "d"), label_size = 20
                   ,nrow=1, rel_widths = c(2, 1))

pmiddle <- pmiddle + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)




lbld <- rbindlist(list(pw3p[sig==1][order(`0.5quant`)][,head(.SD,10)]
                       ,pw3p[sig==1][order(`0.5quant`)][,tail(.SD,10)]
                       ,pw3p[sig_bonf==1]
                       ,pw3p[p.value_v2 < 0.05/.N][order(-abs(fold_change_v2))]
                       )
                  )[,unique(.SD)]

lbld <- rbindlist(list(pw3p[sig==1][order(`0.5quant`)][,head(.SD,10)]
                       ,pw3p[sig==1][order(`0.5quant`)][,tail(.SD,10)]
                       ,pw3p[sig_bonf==1]
                       ,pw3p[p.value_v2 < 0.05/.N][order(-abs(fold_change_v2))][1:10]
                       )
                  )[,unique(.SD)]
lbld[fold_change_v2 < `0.5quant`,ndg:=0.1 - fold_change_v2]
lbld[fold_change_v2 > `0.5quant`,ndg:=1.8 - fold_change_v2]


spatialp <- 
ggplot(pw3p, aes(x=fold_change_v2
                 ,y=`0.5quant`,
                 fill=spatial_re_sd_mean, shape=)) + 
  geom_point(data=pw3p[sig_bonf==1]
             ,aes(shape=0),size=6,stroke=1.2) + 
  geom_point(data=pw3p[p.value_v2 < 0.05/.N]
             ,aes(shape=8),size=3,stroke=1.2) + 
  geom_point(aes(shape=21),size=3) + 
  scale_shape_identity(guide=guide_legend(override.aes=list(size=6))
                       ,breaks=c(0,8)
                       ,name=""
                       ,labels=c("9 Significant genes\n(INLA spatial RE model)"
                                 ,"114 Significant genes\n(non-spatial semi-naive model)")) + 
  theme_bw(base_size = 20) + 
  labs(x="Fold change (without spatial RE)"
       ,y ="Fold change (INLA, with spatial RE)") + 
#  scale_color_gradient2() + 
  scale_fill_gradient2(name = "Std Dev of\nSpatial\nRandom Effect"
                       ,high=scales::muted("red"), mid="white", low=scales::muted("blue")
                       #,midpoint = median(pw3p[["spatial_re_sd_mean"]])
                       ,midpoint = mean(pw3p[["spatial_re_sd_mean"]])
                       ) + 
  geom_abline(slope=1,intercept=0,lty=2) + 
  scale_x_continuous( n.breaks=20) + 
  scale_y_continuous( n.breaks=20) + 
  geom_label_repel(data=lbld
                   ,aes(label=target)
                   ,max.overlaps=20
                   ,nudge_x = lbld[,ndg]
                   ,direction="y"
                   ) + 
  theme(legend.position="bottom") + 
  coord_fixed() 

pw3p[,islbl:=as.numeric(target %in% lbld[["target"]])]
spatialp <- 
ggplot(pw3p[order(islbl)], aes(x=fold_change_v2
                 ,y=`0.5quant`,
                 fill=spatial_re_sd_mean, shape=)) + 
  geom_point(data=pw3p[p.value_v2 < 0.05/.N]
             ,aes(shape=8),size=3,stroke=1.2) + 
  geom_point(aes(shape=21),size=3) + 
  geom_point(data=pw3p[sig_bonf==1]
             ,aes(shape=0),size=6,stroke=1.2) + 
  scale_shape_identity(guide=guide_legend(override.aes=list(size=6))
                       ,breaks=c(0,8)
                       ,name=""
                       ,labels=c("9 Significant genes\n(INLA spatial RE model)"
                                 ,"114 Significant genes\n(non-spatial semi-naive model)")) + 
  theme_bw(base_size = 20) + 
  labs(x="Fold change (without spatial RE)"
       ,y ="Fold change (INLA, with spatial RE)") + 
#  scale_color_gradient2() + 
  scale_fill_gradient2(name = "Std Dev of\nSpatial\nRandom Effect"
                       ,high=scales::muted("red"), mid="white", low=scales::muted("blue")
                       #,midpoint = median(pw3p[["spatial_re_sd_mean"]])
                       ,midpoint = mean(pw3p[["spatial_re_sd_mean"]])
                       ) + 
  geom_abline(slope=1,intercept=0,lty=2) + 
  scale_x_continuous( n.breaks=20) + 
  scale_y_continuous( n.breaks=20) + 
  geom_label_repel(data=lbld
                   ,aes(label=target)
                   ,max.overlaps=20
                   ,nudge_x = lbld[,ndg]
                   ,direction="y"
                   ) + 
  theme(legend.position="bottom") + 
  theme(axis.text.x=element_text(size=18)) + 
  coord_fixed() 

spatialp <- spatialp + 
  geom_hline(yintercept = 1,lty=2) + 
  geom_vline(xintercept = 1,lty=2)

pw1p[p.value < 0.05/.N]
pw1p[p.value < 0.05/.N][,.N,by=.(fold_change < 1)]
pw1p[p.value < 0.05/.N][,.N,by=.(fold_change < 1)][,prop:=N/sum(N)][]

pw2p[p.value_v2 < 0.05/.N]
pw2p[p.value_v2 < 0.05/.N][,.N,by=.(fold_change_v2 < 1)]
pw2p[p.value_v2 < 0.05/.N][,.N,by=.(fold_change_v2 < 1)][,prop:=N/sum(N)][]
pw2p[p.value_v2 < 0.05/.N][p.value_v1 > 0.05/.N]

pd <- merge(metainfo[,.(x_slide_mm, y_slide_mm,tissgrp)]
            ,sre[target == "WNT7A"]
            ,by.x=c("x_slide_mm", "y_slide_mm")
            ,by.y=c("x", "y"))#[order(mean)]

wnt7ap <- 
ggplot(pd[order(mean)], aes(x_slide_mm, y_slide_mm, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw(base_size=20) + 
  scale_color_gradient2(name = "Predicted\nSpatial\nRandom Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = mean(pd[["mean"]]))  + 
  coord_fixed()  + 
  theme(legend.position="bottom") + 
  labs(title="WNT7A Predicted\nSpatial Random Effects")

wnt7ap
pd <- merge(metainfo[,.(x_slide_mm, y_slide_mm,tissgrp)]
            ,sre[target == "CD79A"]
            ,by.x=c("x_slide_mm", "y_slide_mm")
            ,by.y=c("x", "y"))#[order(mean)]

cd79p <- 
ggplot(pd[order(mean)], aes(x_slide_mm, y_slide_mm, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw(base_size=20) + 
  scale_color_gradient2(name = "Predicted\nSpatial\nRandom Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = median(pw2p[["contam_ratio"]]))  + 
  coord_fixed()  + 
  theme(legend.position="bottom", legend.text=element_text(size=8)) + 
  labs(title="CD79A Predicted\nSpatial Random Effects")


cowplot::plot_grid(spatialp + theme(plot.margin = unit(c(0.1,0,0,0.1), "cm")
                                    ,legend.box = "vertical"
                                    ) + 
                     labs(title="Impact of Spatial Random Effect Variation Fold Change Estimation and Inference ")
                   ,cd79p + theme(plot.margin = unit(c(0,0,0,0), "cm"))
                   ,wnt7ap + theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
                   ,rel_widths = c(1, 0.5, 0.5)
                   ,labels = c("a", "b", "c")
                   ,nrow = 1) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)


ggsave("colon_data/deplot3.png"
       ,width=1858*sqrt(25)
       ,height=919*sqrt(25)
       ,dpi=450
       ,units = "px")

plot3efg <- 
cowplot::plot_grid(spatialp + theme(plot.margin = unit(c(0.1,0,0,0.1), "cm")
                                    ,legend.box = "vertical"
                                    ,axis.text.x = element_text(size=12)
                                    ) + 
                     labs(title="Impact of Spatial Random Effect Variation\non Fold Change Estimation and Inference ")
                   ,cd79p + theme(plot.margin = unit(c(0.1,0,0,0), "cm")) 
                   ,wnt7ap + theme(plot.margin = unit(c(0.1,0.1,0,0), "cm"))
                   ,rel_widths = c(1, 0.5, 0.5)
                   ,labels = c("e", "f", "g")
                   ,label_size=20
                   ,nrow = 1) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)

plot3efg






plist2 <- 
lapply(split(metainfo, by="tissgrp"), function(xx){
  ggplot(xx[clust!="B-cell"][order(ct)], aes(y_slide_mm, x_slide_mm, color=ct)) + 
    geom_point(size=0.05, alpha=0.5) + 
    geom_point(data=xx[clust=="B-cell"], size=0.3) + 
    coord_fixed() + 
    theme_bw() + 
    scale_color_manual(values=rev(ctpal2)
                       ,breaks=rev(names(ctpal2))
                       ,name="cell type"
                       ,guide=guide_legend(override.aes=list(size=4, alpha=1))) + 
    theme(legend.position="bottom")
    
})

leg1 <- cowplot::get_plot_component(plist[[1]] + 
                                      #theme(legend.position = "right")
                                      theme(legend.position = "right", legend.text = element_text(size=20), legend.title = element_text(size=20))
                                    , 'guide-box-right', return_all=TRUE)
upside <- 
cowplot::plot_grid(cowplot::plot_grid(plist2$ll + theme(legend.position="none") 
                                      ,cowplot::plot_grid(leg1
                                                          ,plist$ur + theme(legend.position="none")
                                                          ,nrow=2
                                                          ,rel_heights = c(1,3))
                                      ,nrow=1, rel_widths = c(6, 2.8))
                   ,labels=c("a", "")
                   ,label_size = 20) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
upside
ggsave("colon_data/colon_deplot_pt1.png"
       ,height =576*sqrt(30)
       ,width=1163*sqrt(30)
       ,units="px"
       ,dpi=450
       )
cowplot::plot_grid(upside, plot3cde, nrow=2, rel_heights = c(1, 1))




png("colon_data/deplot_3row.png"
    ,res=450
    ,width=19*sqrt(15)
    ,height=8*2.54*sqrt(15)
    ,units = "cm"
    )
cowplot::plot_grid(upside, pmiddle, plot3efg,nrow=3)
dev.off()
cowplot::plot_grid(upside, pmiddle,nrow=2)

png("colon_data/deplot_3row_v2.png"
    ,res=450
    ,width=19*sqrt(9)
    ,height=8*2.54*sqrt(9)
    ,units = "cm"
    )
cowplot::plot_grid(upside, pmiddle, plot3efg,nrow=3)
dev.off()

png("colon_data/deplot_3row_v3.png"
    ,res=450
    ,width=19*sqrt(4)
    ,height=8*2.54*sqrt(4)
    ,units = "cm"
    )
cowplot::plot_grid(upside, pmiddle, plot3efg,nrow=3)
dev.off()













cowplot::plot_grid(pmiddle,plot3efg,nrow=2)
ggsave("colon_data/deplot_2row.png"
       ,dpi=450
       ,width=19*sqrt(15)
       ,height=8*2.54/3*2 *sqrt(15)
       ,units = "cm"
       )




cowplot::plot_grid(upside, plot3cde, nrow=2, rel_heights = c(1, 1))
#ggsave("colon_data/deplottest.png"
#       ,height = 1482*20
#       ,width=1482*40
#       ,units="px"
#       ,dpi=450)

cowplot::plot_grid(leftside, rightside3, nrow=1)

leg2 <- cowplot::get_plot_component(plist_immcut[[1]] , 'guide-box-right', return_all=TRUE)


combleg <- cowplot::plot_grid(cowplot::ggdraw(leg1) + theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              ,cowplot::ggdraw(leg2)+ theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              , nrow=2)



ggsave("colon_data/deplot3cde.png"
       ,width=1858*sqrt(25)
       ,height=919*sqrt(25)
       ,dpi=450
       ,units = "px")


ggsave("colon_data/deplot3.png"
       ,width=1858*sqrt(25)
       ,height=919*sqrt(25)
       ,dpi=450
       ,units = "px")

hbbp <- 
ggplot(sre[target == "HBB"][order(mean)], aes(x, y, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw() + 
  scale_color_gradient2(name = "Predicted Spatial Random Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = median(pw2p[["contam_ratio"]]))  + 
  coord_fixed()  + 
  theme(legend.position="bottom") + 
  labs(title="WNT7A Predicted Spatial Random Effects")

sre_sub <- merge(sre[target=="CCKBR"], metainfo_bcell, by.x=c("x", "y"), by.y=c("x_slide_mm", "y_slide_mm"))
ggplot(sre_sub[target == "CCKBR"][order(mean)], aes(x, y, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw() + 
  scale_color_gradient2(name = "Predicted Spatial Random Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = median(pw2p[["contam_ratio"]]))  + 
#  coord_fixed()  + 
  theme(legend.position="bottom") + 
  labs(title="WNT7A Predicted Spatial Random Effects") + 
  facet_wrap(~tissgrp,scales="free")




ggplot(pw2p, aes(x=-log10(p.value_v2)
                 ,y=-log10(p.value_v1),
                 fill=contam_ratio)) + 
  geom_point(pch=21,size=3) + 
  theme_bw() + 
#  scale_color_gradient2() + 
  scale_fill_gradient2(high=muted("red"), mid="white", low=muted("blue"), midpoint = median(pw2p[["contam_ratio"]])
                       ,guide=guide_legend(override.aes=list(size=4))) + 
  geom_abline(slope=1,intercept=0,lty=2) + 
  scale_x_continuous(trans='log10', n.breaks=10) + 
  scale_y_continuous(trans='log10', n.breaks=10)  


smiDE::volcano(rbindlist(list(pw1[term=="immune_ct_neighbors"][,typ:="naive"]
                              ,pw2[term=="immune_ct_neighbors"][,typ:="v2"]
                              ))
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]] + 
  facet_wrap(~typ,ncol=2)

pdcomb <- 
rbindlist(list(pw1[term=="immune_ct_neighbors"][,typ:="naive"]
                              ,pw2[term=="immune_ct_neighbors"][,typ:="v2"]
                              ))
smiDE::volcano(pdcomb
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]] + 
  geom_line(data=pdcomb, aes(x=log2(fold_change), y= -log10(p.value), group=target)
            ,inherit.aes=FALSE)

ggplot(pw1, aes())
pw1[term=="immune_ct_neighbors"]


ggplot(pw1[term=="immune_ct_neighbors"], aes(x=log2(fold_change),y=-log10(p.value))) + 
  geom_point() + 
  scale_y_continuous(n.breaks=)


vnaive$`immune_ct_neighbors 14.54 / immune_ct_neighbors 9.458`

vp2 <- 
smiDE::volcano(pw2[term=="immune_ct_neighbors"]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)

vp2$`immune_ct_neighbors 14.54 / immune_ct_neighbors 9.458`

inlapw <- pw3[term=="immune_ct_neighbors"]
inlapw[,lcl_bonf:=`9.255831e-06quant`]
inlapw[,ucl_bonf:=`0.9999907quant`]
inlapw[,lcl:=`0.025quant`]
inlapw[,ucl:=`0.975quant`]

inlapw[,sig_bonf:=as.numeric(sign(log(lcl_bonf))==sign(log(ucl_bonf)))]
inlapw[,sig:=as.numeric(sign(log(lcl))==sign(log(ucl)))]
inlapw[,sig:=as.numeric(sign(log(`9.255831e-06quant`))==sign(log(`0.9999907quant`)))]


inlapw

inlapw[,.(sign(log(`9.255831e-06quant`)),sign(log(`0.9999907quant`)))]


vp3 <- 
smiDE::volcano(pw3[term=="immune_ct_neighbors"]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)
vp3$`immune_ct_neighbors 14.54 / immune_ct_neighbors 9.458`




pw
smiDE::volcano(pw[term=="immune_ct_neighbors"]
               ,fold_change_cutoff = 2^0.5
               ,interactive=FALSE
               ,limit_labels = TRUE
               )



dev3 <- readRDS("colon_data/de_results/de_v3.targset__0.rds")
smiDE::results(dev3, "pairwise")



##
ggplot(metainfo[clust %in% c("B-cell", other_immune_types)], aes(x_slide_mm, y_slide_mm, color=clust)) + 
  geom_point(size=0.1) + 
  theme_bw() + 
  coord_fixed() + 
  scale_color_manual(values=rep(unname(pals::alphabet2()),2) #
                     ,guide=guide_legend(override.aes=list(size=4))
  )

ggplot(metainfo_bcell[clust %in% c("B-cell", other_immune_types)], aes(x_slide_mm, y_slide_mm, color=clust)) + 
  geom_point(size=0.1) + 
  theme_bw() + 
  coord_fixed() + 
  scale_color_manual(values=rep(unname(pals::alphabet2()),2) #
                     ,guide=guide_legend(override.aes=list(size=4))
  )



celltype_nbrs <- 
  copy(pre_de_bcell$nblist$adjacency_counts_by_ct)[,other_immune_ct_neighbors:=rowSums(.SD)
                                                   ,.SDcols=c("B-cell", other_immune_types)
  ]



ggplot(celltype_nbrs[match(bcell_cellid, cell_ID)]
       ,aes(x=other_immune_ct_neighbors)) + 
  geom_histogram(bins = 21)

pd <- celltype_nbrs[match(bcell_cellid, cell_ID),hist(other_immune_ct_neighbors)]
pd <- celltype_nbrs[match(bcell_cellid, cell_ID),hist(other_immune_ct_neighbors)]

ggplot(pd, aes(x_slide_mm, y_slide_mm, color=other_immune_ct_neighbors)) + 
  geom_point(size=0.1) + 
  theme_bw() + 
  coord_fixed() 
#  scale_color_manual(values=rep(unname(pals::alphabet2()),2) #
#                     ,guide=guide_legend(override.aes=list(size=4))
#  )
