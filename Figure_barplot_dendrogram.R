setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/LongGene_code/src/')
suppressPackageStartupMessages({
  source('libraries.R')
  library(Seurat)
  library(scrattch.hicat)
  library(dendextend)
  library(dplyr)
  library(matrixStats)
  library(Matrix)
  library(scrattch.io)
  library(clustree)
  library(stringr)
  library(here)
  library(ggpubr)
  library(tidyverse)
  library(feather)
})

setwd('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/')
data_path <- ('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/')
load("processed_data/human_10x_10XALS_integration/integrated_seurat/Cv3_tenxals.final.rda")
#load("processed_data/human_10x_10XALS_integration/integrated_seurat/all.cell.integrated.rda")
#load("processed_data/human_10x_10XALS_integration/integrated_seurat/Cv3_dendrogram.rda")
load("processed_data/human_10x_10XALS_integration/Figure_1/aucell_metadata.3.rda")


#### 
aucell = aucell[aucell$Barcode_new %in% merged_combined@meta.data$sample_id,]
aucell = aucell[,c("Barcode_new","mouse_facs")]

aucell <- aucell[ order(match(aucell$Barcode_new, merged_combined@meta.data$sample_id)), ]
merged_combined_added <- AddMetaData(merged_combined, metadata = aucell)

merged_combined_added@meta.data$mouse_facs.colors = ifelse(merged_combined_added@meta.data$mouse_facs == "0", "#D9D9D9",ifelse(merged_combined_added@meta.data$mouse_facs == "1","#6A51A3", "#41AB5D" ))

########
colors.cluster = unique(merged_combined_added@meta.data$cluster_label.dendro.color)
names(colors.cluster) = unique(merged_combined_added@meta.data$cluster_label.dendro)

colors.facs = unique(merged_combined_added@meta.data$mouse_facs.colors)
names(colors.facs) = unique(merged_combined_added@meta.data$mouse_facs)

p1 <- DimPlot(merged_combined_added, reduction = "umap", group.by = "cluster_label.dendro", cols = colors.cluster , raster = F)  + NoLegend()
p2 <- DimPlot(merged_combined_added, reduction = "umap", group.by = "mouse_facs", cols = colors.facs , raster = F)  + NoLegend()

library(ggpubr)
pdf(file = "processed_data/figures_raw/umap_cluster.pdf", width = 11, height = 11)
ggarrange(p1, p2, ncol = 2, nrow = 1)
dev.off()

# integrated@meta.data$cluster_label.dendro = ifelse(integrated@meta.data$cluster_label.dendro == "Sncg-2", "Lamp5-13", ifelse(integrated@meta.data$cluster_label.dendro == "Sncg-4", "Lamp5-14", ifelse(integrated@meta.data$cluster_label.dendro == "Sncg-5", "Lamp5-15", ifelse(integrated@meta.data$cluster_label.dendro == "Sncg-6", "Lamp5-16", integrated@meta.data$cluster_label.dendro))))
# integrated@meta.data$subclass_label = ifelse(integrated@meta.data$cluster_label.dendro == "Lamp5-13", "Lamp5", ifelse(integrated@meta.data$cluster_label.dendro ==  "Lamp5-14", "Lamp5" , ifelse(integrated@meta.data$cluster_label.dendro ==  "Lamp5-15","Lamp5", ifelse(integrated@meta.data$cluster_label.dendro == "Lamp5-16","Lamp5",  integrated@meta.data$subclass_label))))
# integrated@meta.data$cluster_label.dendro = ifelse(integrated@meta.data$cluster_label.dendro == "Sst-21", "Sst-1", integrated@meta.data$cluster_label.dendro)
# integrated@meta.data$cluster_label.dendro = ifelse(integrated@meta.data$cluster_label.dendro == "OPC-2", "Micro-PVM-3", integrated@meta.data$cluster_label.dendro)
# integrated@meta.data$subclass_label = ifelse(integrated@meta.data$cluster_label.dendro == "Micro-PVM-3", "Micro-PVM", integrated@meta.data$subclass_label)
# 
# 
# cluster = integrated@meta.data[!duplicated(integrated@meta.data$cluster_label.dendro),]
# cluster = cluster[,c(13,14)]
# cluster <- dplyr::arrange(cluster, cluster_label.dendro)
# cluster$cluster_label.dendro.id = c(1:nrow(cluster))
# cluster = cluster[,c(1,3)]
# integrated@meta.data = merge(integrated@meta.data,cluster,by=c("cluster_label.dendro")) 
# 
# color = merged_combined@meta.data[!duplicated(merged_combined@meta.data$cluster_label.dendro),]
# color  = color[,c(2,23,25,26)]
# integrated@meta.data = merge(integrated@meta.data,color,by=c("cluster_label.dendro")) 
# row.names(integrated@meta.data) = integrated@meta.data$sample_id
# 
# 
# pdf(file = "processed_data/figures/umap_cluster.pdf", width = 11, height = 11)
# DimPlot(integrated, reduction = "umap", group.by = "cluster_label.dendro", raster = F, cols = ucols) + NoLegend()
# dev.off()

###### Split object ######
# integrated.tenx <- subset(integrated, subset = orig.ident == "tenxals")
# integrated.sm <- subset(integrated, subset = orig.ident == "sm")
# 
# pdf(file = "processed_data/figures/umap_clusterSS4.pdf", width = 11, height = 11)
# DimPlot(integrated.sm, reduction = "umap", group.by = "cluster_label.dendro", raster = F, cols = ucols) + NoLegend()
# dev.off()
# 
# 
# colors <- RColorBrewer::brewer.pal(11, "RdBu")[c(2,3,8,9)]
# nb.cols <- 4
# mycolors <- colorRampPalette(colors)(nb.cols)
# ucols.layer = mycolors
# pdf(file = "processed_data/figures/umap_layer.pdf", width = 11, height = 11)
# DimPlot(integrated, reduction = "umap", group.by = "predicted.id", raster = F,cols = ucols.layer) 
# dev.off()
# 
# 
# ucols.subclass = unique(integrated@meta.data$subclass_color)
# pdf(file = "processed_data/figures/umap_subclass.pdf", width = 11, height = 11)
# DimPlot(integrated, reduction = "umap", group.by = "subclass_label", raster = F, cols = ucols.subclass) 
# dev.off()

####### MAke dendrograme and save to file ############
load(paste0(data_path,"processed_data/human_10x_10XALS_integration/Figure_1/human_meta.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/Figure_1/human_als.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/Figure_1/aucell_metadata.3.rda"))



als = als[,-c(28,29)]
aucell = aucell[,c(16,23)]
als = merge(als,aucell, by=c("row.names"))

pdf(file = "processed_data/figures/dendrogram.pdf", width = 16, height = 7 )
dend %>%
  rotate(as.character(dend_order)) %>% #rotate to match labels new order
  plot(main = "Motor cortex taxonomy human")
dev.off()


#Plot 1: layer heatmap
ssv4 = meta %>% drop_na(layer_label)
#check all SS clusters are represented
ss_cl <- ssv4 %>%
  filter(orig.ident  %>% str_detect("sm")) %>%
  distinct(cluster_label.dendro) %>% pull() 
dend_order[which(dend_order %in% ss_cl == FALSE)]

not_ss_cl = setdiff(dend_order,ss_cl)

#check layer info
ssv4 %>% distinct(layer_label) 

#make plotting object
to_plot_1 <- ssv4 %>%
  filter(orig.ident %>% str_detect("sm")) %>%
  
  group_by(cluster_label.dendro) %>%
  mutate(cluster_n = n()) %>%
  ungroup() %>%
  
  group_by(cluster_label.dendro, layer_label) %>%
  mutate(cluster_layer_n = n()) %>%
  ungroup() %>%
  
  mutate(layer_prop = cluster_layer_n / cluster_n) %>%
  
  distinct(cluster_label.dendro, layer_label, layer_prop) %>%
  
  #factorize
  mutate(layer = layer_label %>% as_factor() %>% fct_relevel(c("1", "2", "3", "5", "6")) %>% fct_rev(),
         cluster = cluster_label.dendro %>% as_factor() %>% fct_relevel(ss_cl))


#### Add layer proportions for missing cluster by taking the avergae over each subclass ####
to_plot_tmp = to_plot_1
to_plot_tmp$cluster_label.dendro = sub("-[^-]+$", "", to_plot_tmp$cluster_label.dendro)

to_plot_2 <- to_plot_tmp %>%
  group_by(cluster_label.dendro,layer_label) %>%
  mutate(layer_prop = mean(layer_prop))  %>% 
  distinct(layer_prop) %>% 
  mutate(cluster_label.dendro  = case_when(
  cluster_label.dendro == "L6 CT" ~ "L6 CT-5",
  cluster_label.dendro == "L2/3 IT" ~ "L2/3 IT-8",
  cluster_label.dendro == "Astro" ~ "Astro-2",
  cluster_label.dendro == "Micro-PVM" ~ "Micro-PVM-2"
)) %>% drop_na(cluster_label.dendro) %>%
  #factorize
  mutate(layer = layer_label %>% as_factor() %>% fct_relevel(c("1", "2", "3", "5", "6")) %>% fct_rev(),
         cluster = cluster_label.dendro %>% as_factor() %>% fct_relevel(not_ss_cl))
  
to_plot = rbind(to_plot_1,to_plot_2)


order = as.data.frame(dend_order)
colnames(order)[1] = "cluster"
order$position = row.names(order)

to_plot = merge(to_plot,order,by=c("cluster"))
to_plot= to_plot[order(as.numeric(as.character(to_plot$position))),]


p1 <- to_plot %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = layer, fill = layer_prop)) +
  
  scale_fill_gradient(low = "white", high = "black") +
  
  labs(fill = "Proportion", 
       y = "Layer",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

p1

#Plot 1b: CSN and CPN heatmap
#check all Cv3 clusters are represented
cv3_cl <- als %>%
  filter(orig.ident  %>% str_detect("tenxals")) %>%
  distinct(cluster_label.dendro) %>% pull() 
dend_order[which(dend_order %in% cv3_cl == FALSE)]


#make plotting object
to_plot_3 <- als %>%
  filter(orig.ident %>% str_detect("tenxals")) %>%
  
  group_by(cluster_label.dendro) %>%
  mutate(cluster_n = n()) %>%
  ungroup() %>%
  
  group_by(cluster_label.dendro, mouse_facs) %>%
  mutate(cluster_facs_n = n()) %>%
  ungroup() %>%
  
  mutate(facs_prop = cluster_facs_n / cluster_n) %>%
  
  distinct(cluster_label.dendro, mouse_facs, facs_prop) %>%
  
  #factorize
  mutate(facs = mouse_facs %>% as_factor() %>% fct_relevel(c("0", "1", "2")) %>% fct_rev(),
         cluster = cluster_label.dendro %>% as_factor() %>% fct_relevel(cv3_cl))


order = as.data.frame(dend_order)
colnames(order)[1] = "cluster"
order$position = row.names(order)

to_plot_3 = merge(to_plot_3,order,by=c("cluster"))
to_plot_3= to_plot_3[order(as.numeric(as.character(to_plot_3$position))),]
to_plot_3 = to_plot_3[to_plot_3$mouse_facs != "0",]

p1b <- to_plot_3 %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = mouse_facs, fill = facs_prop)) +
  
  scale_fill_gradient(low = "white", high = "black") +
  
  labs(fill = "Proportion", 
       y = "Layer",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

pdf(p1b,file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/Heatmap_csn.pdf")
p1b
dev.off()

#Plot 2: proportions by class
to_plot_4 <- als %>%
  filter(orig.ident %>% str_detect("tenxals")) %>%
  
  group_by(Batch, class_label) %>%
  mutate(donor_n = n()) %>%
  ungroup() %>%
  
  group_by(Sample_ID, cluster_label.dendro) %>%
  mutate(donor_cluster_n = n()) %>%
  ungroup() %>%
  
  mutate(cluster_prop = donor_cluster_n / donor_n * 100) %>%
  
  group_by(cluster_label.dendro) %>%
  summarise(mean_prop = mean(cluster_prop),
            sd_prop = sd(cluster_prop)) %>%
  ungroup() %>%
  
  mutate(cluster = cluster_label.dendro %>% as_factor() %>% fct_relevel(dend_order)) 

order = as.data.frame(dend_order)
colnames(order)[1] = "cluster"
order$position = row.names(order)

to_plot_4 = merge(to_plot_4,order,by=c("cluster"))

tmp <- als %>% arrange(cluster_label.dendro) %>% distinct(cluster_label.dendro, cluster_label.dendro.color)


to_plot_4 = merge(to_plot_4,tmp,by=c("cluster_label.dendro"))
to_plot_4= to_plot_4[order(as.numeric(as.character(to_plot_4$position))),]

colors_use <- to_plot_4$cluster_label.dendro.color
names(colors_use) <- to_plot_4$cluster_label.dendro.color

p2 <- to_plot_4 %>%
  ggplot() +
  geom_errorbar(aes(x = cluster, ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop)) +
  geom_point(aes(x = cluster, y = mean_prop, color = colors_use), size = 2, alpha = 0.6) +
  
  scale_color_manual(values = colors_use) +
  scale_y_log10() +
  
  labs(y = "Percent of class",
       x = "") +  
  
  theme(aspect.ratio = 12/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 1, linetype = 2),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

p2

#Plot 3: Stacked barplot dataset
to_plot_6<- meta %>% 
  
  mutate(
    dataset_label = case_when(
      orig.ident == "tenx" ~ "Cv3 all layers",
      orig.ident %>% str_detect("sm") ~ "SSv4 layer dissected",
      TRUE ~ "Cv3 all layers ALS/FTD"),
    
    dataset_label  = dataset_label %>% as_factor() %>% fct_relevel(c("Cv3 all layers", "Cv3 all layers ALS/FTD", "SSv4 layer dissected")),
    
    cluster = cluster_label.dendro %>% as_factor() %>% fct_relevel(dend_order))


p4 <- to_plot_6 %>%
  ggplot() +
  geom_bar(aes(x = cluster, fill = dataset_label), position = "fill") +
  
  scale_fill_manual(values = c("#FF6600", "#66CCFF", "#CC0099")) +
  
  labs(x = "",
       y = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
  )

p4



#Plot 5: Violin plot genes detected
to_plot_7 <- meta %>% 
  filter(orig.ident %>% str_detect("tenx")) %>%
  mutate(cluster = cluster_label.dendro %>% as_factor() %>% fct_relevel(dend_order))


order = as.data.frame(dend_order)
colnames(order)[1] = "cluster"
order$position = row.names(order)

to_plot_7 = merge(to_plot_7,order,by=c("cluster"))

tmp <- als %>% arrange(cluster_label.dendro) %>% distinct(cluster_label.dendro, cluster_label.dendro.color)
colnames(tmp)[2] <- "new_color"

to_plot_7 = merge(to_plot_7,tmp,by=c("cluster_label.dendro"))
to_plot_7= to_plot_7[order(as.numeric(as.character(to_plot_7$position))),]

colors_use <- to_plot_7$cluster_label.dendro.color
names(colors_use) <- to_plot_7$cluster_label.dendro.color

p5 <- to_plot_7 %>%
  ggplot() +
  geom_violin(aes(x = cluster, y = nFeature_RNA, fill = colors_use), scale = "width") +
  geom_abline(slope = 0, intercept = 3000, linetype = 2, color = "lightgrey") +
  geom_abline(slope = 0, intercept = 6000, linetype = 2, color = "lightgrey") +
  geom_abline(slope = 0, intercept = 9000, linetype = 2, color = "lightgrey") +
  
  
  scale_fill_manual(values = colors_use) +
  
  labs(x = "",
       y = "Genes detected (10x)") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
  )

p5

library(data.table)
load("processed_data/human_10x_10XALS_integration/Figure_3/seurat_combat.rda")
median = sce.to.seurat@meta.data[,c("nFeature_RNA","subclass_label")]
stats = setDT(median)[,list(Mean=mean(nFeature_RNA), Max=max(nFeature_RNA), Min=min(nFeature_RNA), Median=as.numeric(median(nFeature_RNA)), Std=sd(nFeature_RNA)), by=subclass_label]

#load("processed_data/human_10x_10XALS_integration/integration_exc.final.rda")
load("processed_data/human_10x_10XALS_integration/integrated_seurat/all.cell.integrated.rda")
load("processed_data/human_10x_10XALS_integration/integrated_seurat/Cv3_tenxals.final.rda")
meta_data_cl <- merged_combined@meta.data[,c("cluster_label.dendro", "cluster_label.dendro.deg", "subclass_label", "cluster_label.dendro.color", "subclass_color", "cluster_label.dendro.deg.colors")]
meta_data_cl <- meta_data_cl[!duplicated(meta_data_cl$cluster_label.dendro),]


metadata_new <- integrated@meta.data
metadata_new <- merge(metadata_new,meta_data_cl,by=c("cluster_label.dendro"))
row.names(metadata_new) <- metadata_new$sample_id

integrated <- subset(integrated, sample_id %in% metadata_new$sample_id)
integrated <- AddMetaData(object = integrated,metadata = metadata_new)

Idents(integrated) <- "cluster_label.dendro.deg"

rownames(integrated) %>% duplicated() %>% table()
Idents(integrated) %>% anyNA()

p<- DotPlot(integrated, features = c("MALAT1"),dot.scale = 4) + RotatedAxis() +scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")+
  theme(text = element_text(face = "bold"),
        axis.text.x=element_text(angle=45, hjust=1, size=5),
        legend.text=element_text(size=5),
        legend.title=element_text(size=5))

colors.cluster = unique(integrated@meta.data$cluster_label.dendro.deg.colors)
names(colors.cluster) = unique(integrated@meta.data$cluster_label.dendro.deg)


p1 <- DimPlot(integrated, reduction = "umap", group.by = "cluster_label.dendro.deg", cols = colors.cluster , raster = F)  + NoLegend()
p2 <- FeaturePlot(object = integrated,features = "MALAT1", raster = F)
#p3 <- VlnPlot(object = integrated,features = "MALAT1")

ggsave(filename = 'processed_data/figures_raw/bubble_plot_MALAT1.pdf', p)
ggsave(filename = 'processed_data/figures_raw/dimplot_plot_MALAT1.pdf', p1)
ggsave(filename = 'processed_data/figures_raw/feature_plot_MALAT1.pdf', p2)


tmp <- SplitObject(integrated, split.by = "orig.ident")
p1<-DimPlot(tmp$sm, reduction = "tsne", group.by = "cluster_label.dendro.deg", cols = colors.cluster , raster = F)  + NoLegend()
p2<-DimPlot(tmp$tenx, reduction = "tsne", group.by = "cluster_label.dendro.deg", cols = colors.cluster , raster = F)  + NoLegend()
p3<-DimPlot(tmp$tenxals, reduction = "tsne", group.by = "cluster_label.dendro.deg", cols = colors.cluster , raster = F)  + NoLegend()

p1<-FeaturePlot(object = tmp$sm,features = "MALAT1")
p2<-FeaturePlot(object = tmp$tenx,features = "MALAT1")
p3<-FeaturePlot(object = tmp$tenxals,features = "MALAT1")



colors.cluster = unique(integrated@meta.data$subclass_color)
names(colors.cluster) = unique(integrated@meta.data$subclass_label)


p1 <- DimPlot(integrated, reduction = "umap", group.by = "subclass_label", cols = colors.cluster , raster = F)  + NoLegend()
p2 <- FeaturePlot(integrated, features = "nFeature_RNA", raster = F)  
ggsave(filename = 'processed_data/figures_raw/dimplot_subclass_nGenes.pdf', p1)
ggsave(filename = 'processed_data/figures_raw/feature_plot_nGenes.pdf', p2)



