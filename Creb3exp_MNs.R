setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/LongGene_code/src/')
suppressPackageStartupMessages({
  source('libraries.R')
  library(Seurat)
  library(SeuratDisk)
  library(dendextend)
  library(dplyr)
  library(matrixStats)
  library(stringr)
  library(hdf5r)
  library(rhdf5)
  library(DESeq2)
  library(here)
  library(DEFormats)
  library(WGCNA)
  library(plotrix)
  library(tidyverse)
  library(Matrix)
  library(ggrepel)
  library(reshape2)
})
setwd('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/')
source('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/scripts_github/helper_sc_functions.R')
library(fastMatMR)
# Intermediate save directory for .RDS files:
load.dir <- '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/mouse_human_spinal/ms_spinal_cord/'

# Save directory for plot output:
save.dir <- '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/mouse_human_spinal/ms_spinal_cord/'
dir.create(file.path(save.dir))

graph.dir <- paste(save.dir, '/graphs/', sep ='')
dir.create(file.path(graph.dir))

# enable high maximum RAM
options(future.globals.maxSize= 4194304000)
### Commented out, but this is how to read in saved RDS file (RAM intensive process, best to only do once then save though saving is time-intensive)
mouse_spinal_cord <- readRDS(paste0(load.dir,'lepichon_gitler_integrated.RDS'))
mouse_spinal_cord <- subset(mouse_spinal_cord, subset = group == 'LePichon')
mouse_spinal_cord_UMAP <- DimPlot(mouse_spinal_cord, reduction = "umap", group.by = "cholinergictypes", label = T) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
mouse_spinal_cord_UMAP <- rasterize(mouse_spinal_cord_UMAP, layers='Point', dpi=1000)
mouse_spinal_cord_UMAP
#ggsave(mouse_spinal_cord_UMAP, filename = paste(graph.dir, 'UMAP_spinal_cord.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")

Idents(mouse_spinal_cord) <- mouse_spinal_cord$cholinergictypes
markers_cholinergictypes <- FindAllMarkers(mouse_spinal_cord, only.pos = TRUE)
p_dot<- DotPlot(mouse_spinal_cord,features = c("Pitx2", "C1qtnf7", "Pax8", "Gpc4", "Pax2","Ebf2", "Nxph2", "Bcl11a", "Pou6f2", "Aox1", "Grin3b","Esrrb", "Rreb1", "Chst9","Pdgfd", "Ahnak2", "Pde8a", "Vwa5b1", "Mme", "Fbn2", "Gnb4", "Cpa6", "Crebrf" ,"Creb3l2", 
                                              "Creb3", "Specc1"))+NoLegend()
#ggsave(p_dot, filename = paste(graph.dir, 'DotPlot_mouse.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")

mn.markers <- FindMarkers(mouse_spinal_cord, ident.1 = "Skeletal Motor Neurons", ident.2 = "Visceral Motor Neurons", max.cells.per.ident = 500, logfc.threshold = 0, min.pct = 0)
mn.markers <- FindMarkers(mouse_spinal_cord, ident.1 = "Cholinergic Interneurons", ident.2 = "Visceral Motor Neurons", max.cells.per.ident = 500, logfc.threshold = 0, min.pct = 0)



alphamns = sceasy::convertFormat("mouse_human_spinal/ms_spinal_cord/h5_files/alphamns.h5ad", from = "anndata", to = "seurat")
alphamns@meta.data$skeletal_subtype <- "Alpha Motor Neurons"
gammamns = sceasy::convertFormat("mouse_human_spinal/ms_spinal_cord/h5_files/gammamns.h5ad", from = "anndata", to = "seurat")


mns.combined <- merge(alphamns, y = gammamns, add.cell.ids = c("alphamns", "gammamns"), project = "skeletal_subtype")
mns.combined
#mns.combined <- subset(mns.combined, subset = skeletal_subtype == "Beta(?) Motor Neurons", invert = TRUE)
mns.combined <- FindVariableFeatures(mns.combined, selection.method = "vst")
all.genes <- rownames(mns.combined)
mns.combined <- ScaleData(mns.combined, features = all.genes)
# PCA
mns.combined <- RunPCA(mns.combined, features = VariableFeatures(object = mns.combined))
# Cluster
mns.combined <- FindNeighbors(mns.combined, dims = 1:10)
mns.combined <- FindClusters(mns.combined, resolution = 2)
# UMAP
mns.combined <- RunUMAP(mns.combined, dims = 1:10)
DimPlot(mns.combined, reduction = "umap", group.by = "skeletal_subtype", label = TRUE)+NoLegend()

Idents(mns.combined) <- mns.combined@meta.data$skeletal_subtype
DotPlot(mns.combined,features="Creb3l2")
VlnPlot(object = mns.combined,features = "Creb3l2")

Idents(mns.combined) <- mns.combined$skeletal_subtype
mn.markers <- FindMarkers(mns.combined, ident.1 = "Alpha Motor Neurons", ident.2 = "Gamma Motor Neurons", max.cells.per.ident = 500, logfc.threshold = 0, min.pct = 0)
mn.markers <- FindMarkers(mns.combined, ident.1 = "Gamma Motor Neurons", ident.2 = "Beta(?) Motor Neurons", max.cells.per.ident = 1000, logfc.threshold = 0, min.pct = 0)

write.csv(mn.markers, file=paste(graph.dir, 'Motor_neurons_mouse_markers.csv', sep='/'))

mn.markers$gene <- rownames(mn.markers)
rownames(mn.markers) <- NULL
mouse_mn_volcano_plot <- make_volcano_plot(mn.markers, 3, labeled_genes = c('Tpd52l1', 'Creb3', "Creb3l2"), p_val_cutoff = 0.05)
mouse_mn_volcano_plot <- rasterize(mouse_mn_volcano_plot, layers="Point", dpi=500) + NoLegend()
mouse_mn_volcano_plot



#########################################
#########################################
###### Analysis of human data ###########
#########################################
#########################################
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(gplots)
library(svglite)
library(ggrastr)
library(matrixStats)

library(future)
plan("multisession", workers = 12)
options(future.globals.maxSize= +Inf)
options(future.rng.onMisuse="ignore")


fig_dir <- '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/mouse_human_spinal/hs_spinal_cord/figures/'
dir.create(file.path(fig_dir))

csv_dir <- '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/mouse_human_spinal/hs_spinal_cord/csv_files/'
dir.create(file.path(csv_dir))


gautier_mns <- readRDS('mouse_human_spinal/hs_spinal_cord/GSE228778_gautier_mns.rds')
gautier_prop_intronic <- readRDS('mouse_human_spinal/hs_spinal_cord/GSE228778_gautier_prop_intronic.rds')

Idents(gautier_prop_intronic) <- gautier_prop_intronic@meta.data$top_level_annotation
p <- DimPlot(object = gautier_prop_intronic, reduction = 'umap', group.by = "top_level_annotation", label = TRUE) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
p <- rasterize(p, layers='Point', dpi=500)
p

FeaturePlot(gautier_prop_intronic, features = c("CSF1R"))
DotPlot(gautier_prop_intronic,features="ST3GAL6")

###### Subcluster neurons #####
gautier_neurons <- subset(gautier_prop_intronic, idents = "Neurons")
gautier_neurons <- FindVariableFeatures(gautier_neurons, selection.method = "vst", nfeatures = 33000)
all.genes <- rownames(gautier_neurons)
gautier_neurons <- ScaleData(gautier_neurons, features = all.genes)
# PCA
gautier_neurons <- RunPCA(gautier_neurons, features = VariableFeatures(object = gautier_neurons))
# Cluster
gautier_neurons <- FindNeighbors(gautier_neurons, dims = 1:50)
gautier_neurons <- FindClusters(gautier_neurons, resolution = 0.5)
# UMAP
gautier_neurons <- RunUMAP(gautier_neurons, dims = 1:50)
DimPlot(gautier_neurons, reduction = "umap", label = TRUE)
# Glial marker genes
FeaturePlot(gautier_neurons, features = c("MBP"))
FeaturePlot(gautier_neurons, features = c("GFAP"))
FeaturePlot(gautier_neurons, features = c("CSF1R"))

# Remove cluster 9 (glial/neuronal doublets) and re-analyze
gautier_neurons <- subset(gautier_neurons, idents = c("10"), invert = TRUE)

gautier_neurons <- FindVariableFeatures(gautier_neurons, selection.method = "vst", nfeatures = 33000)
all.genes <- rownames(gautier_neurons)
gautier_neurons <- ScaleData(gautier_neurons, features = all.genes)
# PCA
gautier_neurons <- RunPCA(gautier_neurons, features = VariableFeatures(object = gautier_neurons))
# Cluster
gautier_neurons <- FindNeighbors(gautier_neurons, dims = 1:50)
gautier_neurons <- FindClusters(gautier_neurons, resolution = 2.5)
# UMAP
gautier_neurons <- RunUMAP(gautier_neurons, dims = 1:50)
gautier_neurons_UMAP <- DimPlot(gautier_neurons, reduction = "umap", group.by = "seurat_clusters", label = T) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
gautier_neurons_UMAP <- rasterize(gautier_neurons_UMAP, layers='Point', dpi=1000)
gautier_neurons_UMAP
ggsave(gautier_neurons_UMAP, filename = paste(fig_dir, 'UMAP_spinal_cord.svg', sep='/'), device='svg', width = 3, height = 3, units = "in")

# Marker gene
Idents(gautier_neurons) <- gautier_neurons$seurat_clusters
cluster50.markers <- FindMarkers(gautier_neurons, ident.1 = 52, only.pos = TRUE)
head(cluster50.markers, n = 52)

gautier_neurons
gautier_neurons@meta.data$motor_neurons <- ifelse(gautier_neurons$seurat_clusters  == "52", yes = "Motor Neurons", no = "Other Neurons")


p2 <- DimPlot(gautier_neurons, reduction = "umap", group.by = "motor_neurons", label = T) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
p2
#ggsave(p2, filename = paste(fig_dir, 'UMAP_MNs_vs_Other_CREB3_MNs.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")

p3 <- DotPlot(object = gautier_neurons,group.by = "motor_neurons", features = c("SLC5A7", "CREB3", "CREB3L2","CREBRF"))
p3

FeaturePlot(object = gautier_neurons,features = "CREB5")
#ggsave(p3, filename = paste(fig_dir, 'DotPlot_CREB3_MNs.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")

Idents(gautier_neurons) <- gautier_neurons$motor_neurons
mn.markers <- FindMarkers(gautier_neurons, ident.1 = "Motor Neurons", ident.2 = "Other Neurons", max.cells.per.ident = 500, logfc.threshold = 0, min.pct = 0)
write.csv(mn.markers, file=paste(csv_dir, 'Motor_neurons_markers.csv', sep='/'))

mn.markers$gene <- rownames(mn.markers)
rownames(mn.markers) <- NULL

gautier_mn_volcano_plot <- make_volcano_plot(mn.markers, 3, labeled_genes = c('CREB3L2', 'CHAT', "CREBRF"), p_val_cutoff = 0.05)
gautier_mn_volcano_plot <- rasterize(gautier_mn_volcano_plot, layers="Point", dpi=500) + NoLegend()
gautier_mn_volcano_plot

Idents(gautier_neurons) <- gautier_neurons$motor_neurons
gautier_mns <- subset(gautier_neurons, idents = "Motor Neurons")

gautier_mns <- FindVariableFeatures(gautier_mns, selection.method = "vst", nfeatures = 40)
all.genes <- rownames(gautier_mns)
gautier_mns <- ScaleData(gautier_mns, features = all.genes)
# PCA
gautier_mns <- RunPCA(gautier_mns, features = VariableFeatures(object = gautier_mns), approx = FALSE)
# Cluster
gautier_mns <- FindNeighbors(gautier_mns, dims = 1:2)
gautier_mns <- FindClusters(gautier_mns, resolution = 0.4, graph.name = "RNA_snn")
DimPlot(gautier_mns, reduction = "pca") + ggtitle("Subtype") + xlim(-6,12) + ylim(-3,7.5) + theme(plot.title = element_text(hjust = 0.5))

gautier_mns@meta.data$motor_neuron_subtype <- ifelse(gautier_mns$seurat_clusters  == "0", yes = "Gamma", no = "Alpha")
Idents(gautier_mns) <- gautier_mns$motor_neuron_subtype

subtype_pca <- DimPlot(gautier_mns, reduction = "pca") + ggtitle("Subtype") + xlim(-6,12) + ylim(-3,7.5) + theme(plot.title = element_text(hjust = 0.5))
subtype_pca
ggsave(subtype_pca, filename = paste(fig_dir, 'subptype_pca.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")

creb5_pca <- FeaturePlot(gautier_mns, reduction = "pca", features = c("CREB5")) + xlim(-6,12) + ylim(-3,7.5)
creb5_pca
ggsave(creb5_pca, filename = paste(fig_dir, 'creb5_pca.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")

tpd52l1_pca <- FeaturePlot(gautier_mns, reduction = "pca", features = c("TPD52L1")) + xlim(-6,12) + ylim(-3,7.5)
tpd52l1_pca
ggsave(tpd52l1_pca, filename = paste(fig_dir, 'stpd52l1_pca.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")

creb3l2_pca <- FeaturePlot(gautier_mns, reduction = "pca", features = c("CREB3L2")) + xlim(-6,12) + ylim(-3,7.5)
creb3l2_pca
ggsave(creb3l2_pca, filename = paste(fig_dir, 'creb3l2_pca.pdf', sep='/'), device='pdf', width = 3, height = 3, units = "in")


#########################################################
#########################################################
###### Analysis of human developemental  data ###########
#########################################################
#########################################################









