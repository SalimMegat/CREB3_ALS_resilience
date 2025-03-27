
setwd('~/Desktop/Analysis_fus/LongGene_code/src/')
suppressPackageStartupMessages({
  source('libraries.R')
  library(Seurat)
  library(dendextend)
  library(dplyr)
  library(matrixStats)
  library(stringr)
  library(here)
  library(DESeq2)
  library(DEFormats)
  library(multiGSEA)
  library("scProportionTest")
  library(tidyverse)
  library(feather)
})
setwd('~/Desktop/single_cell/sc_als/')
data_path <- here('~/Desktop/single_cell/sc_als/')
here(data_path)

load(here(data_path,"processed_data/human_10x_10XALS_integration/Figure_2/Cv3_tenxals.final.noSST.rda"))
load(here(data_path, "processed_data/human_10x_10XALS_integration/integrated_seurat/markers.cells.rda"))
load(here(data_path, "processed_data/human_10x_10XALS_integration/Figure_2/AggregatedCountsCombat.rda"))
load(here(data_path, "processed_data/human_10x_10XALS_integration/Figure_1/aucell_metadata.3.rda"))
markers_all = unique(c(markers_exc,markers_glia,markers_inh))


## transpose ###
cts = rna$RNA
cts.t = t(cts)

### convert to dataframe ###
cts.t = as.data.frame(cts.t)
cts.t[1:10,1:10]

### get value to split !!!!!
splitRows <- gsub('_.*', '',row.names(cts.t))

cts.split <- split.data.frame(cts.t, f = factor(splitRows))

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- sub("^[^_]*_", "", rownames(x))
  t(x)
})

dge_groups = meta.sub[!duplicated(meta.sub$cluster_label.dendro.deg),]
#whole$DGE_group = gsub("_", "-", whole$DGE_group)
dge_groups = dge_groups$cluster_label.dendro.deg


creb3.bound = read.table('processed_data/human_10x_10XALS_integration/Figure_3/CREB3_bound.genes', header = T)

###########3
dge_array <- lapply(dge_groups, function(cell) {
  message('*******************************************************************************')
  message(paste(Sys.time(), 'processing', cell))
  
  dir = paste0('processed_data/celltype_counts/edgeR_results_Group/')
  
  file1 = read.table(file = paste0(dir,"edgeR_dge_results_sALSvsControl_",cell , '.txt'), header = T)
  file1$cell = cell 
  file1$group = "sALS"
  file1$fcsign = sign(file1$logFC)
  file1$metric= file1$mlog10PValue/file1$fcsign
  file2 = read.table(file = paste0(dir,"edgeR_dge_results_C9ALSvsControl_",cell , '.txt'), header = T)
  file2$cell = cell 
  file2$group = "c9ALS"
  file2$fcsign = sign(file2$logFC)
  file2$metric= file2$mlog10PValue/file2$fcsign
  file3 = read.table(file = paste0(dir,"edgeR_dge_results_sFTLDvsControl_",cell , '.txt'), header = T)
  file3$cell = cell 
  file3$group = "sFTLD"
  file3$fcsign = sign(file3$logFC)
  file3$metric= file3$mlog10PValue/file3$fcsign
  file4 = read.table(file = paste0(dir,"edgeR_dge_results_C9FTLDvsControl_",cell , '.txt'), header = T)
  file4$cell = cell 
  file4$group = "c9FTLD"
  file4$fcsign = sign(file4$logFC)
  file4$metric= file4$mlog10PValue/file4$fcsign
  
  tmp <- rbind(file1,file2,file3,file4)
  
})

final_deg <- do.call(rbind,dge_array)
final_deg$group <- factor(final_deg$group, levels = c( "sALS", "c9ALS", "sFTLD", "c9FTLD"))
final_deg_creb3 = final_deg[final_deg$gene_name %in% creb3.bound$Gene,]


final_deg_creb3 = aggregate(final_deg_creb3[, 11], list(final_deg_creb3$cell,final_deg_creb3$group), mean)
colnames(final_deg_creb3) = c("cell","Group","Zscore")

to_plot_1 = final_deg_creb3


dend_order <- labels(dend.labeled)
order = as.data.frame(dend_order)
colnames(order)[1] = "cell"
order$position = row.names(order)

to_plot_1 = merge(to_plot_1,order,by=c("cell"))
to_plot_1= to_plot_1[order(as.numeric(as.character(to_plot_1$position))),]

to_plot_1$Group <- factor(to_plot_1$Group,levels = c("c9ALS","sALS", "c9FTLD", "sFTLD"))

p1 <- to_plot_1 %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = Group, fill = Zscore)) +
  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  
  scale_fill_gradientn(
    colors=c("blue","white","red"),
    values=rescale(c(-0.3,0.13,0.5)),
    limits=c(-0.3,0.5)
  ) +
  
  labs(fill = "Creb3_Zscore", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

pdf(p1,file = 'processed_data/figures_raw/creb3_zscore_hetmap.pdf')
p1
dev.off()

######### Creb3 targets genes in CSN mouse RNAseq #############
csn_edgeR_deg.post = read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/edgeR_dge_results.txtCSN_post.txt', header = T)
cpn_edgeR_deg.post = read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/edgeR_dge_results.txtCPN_post.txt', header = T)
csn_edgeR_deg.pre = read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/edgeR_dge_results.txtCSN_pre.txt', header = T)
cpn_edgeR_deg.pre = read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/edgeR_dge_results.txtCPN_pre.txt', header = T)

csn_edgeR_deg.post$cell = "csn"
cpn_edgeR_deg.post$cell = "cpn"
cpn_edgeR_deg.post$time = "post"
csn_edgeR_deg.post$time = "post"

csn_edgeR_deg.pre$cell = "csn"
cpn_edgeR_deg.pre$cell = "cpn"
csn_edgeR_deg.pre$time = "pre"
cpn_edgeR_deg.pre$time = "pre"

deg_file.mouse = rbind(csn_edgeR_deg.post,cpn_edgeR_deg.post,csn_edgeR_deg.pre,cpn_edgeR_deg.pre)


human_homologous =  read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/orthologous.cleaned', header = T)

deg_file.mouse = merge(deg_file.mouse,human_homologous,by=c("gene_name"))
deg_file.mouse = deg_file.mouse[deg_file.mouse$human_genes %in% creb3.bound$Gene,]
deg_file.mouse$fcsign = sign(deg_file.mouse$logFC)
deg_file.mouse$metric= deg_file.mouse$mlog10PValue/deg_file.mouse$fcsign

deg_file.mouse.csn = deg_file.mouse[deg_file.mouse$cell == "csn",]

pdf('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/cpn_aging.pdf', height = 4,width = 7 )


dodge <- position_dodge(width = 0.4)
deg_file.mouse.csn$time <- factor(deg_file.mouse.csn$time, levels = c( "pre", "post"))

p<-ggplot(deg_file.mouse.csn, aes(x=time, y=metric, fill=cell))
p + geom_violin(
  aes(color = cell), trim = FALSE,
  position = dodge
) +
  geom_boxplot(width=.1, outlier.colour=NA, position = dodge)
dev.off()

##Compute the analysis of variance
res.aov3 <- aov(metric ~ time, data = deg_file.mouse.csn)
summary(res.aov3)
TukeyHSD(res.aov3, which = "age:genotype")
# Summary of the analysis
plot(res.aov3, 2)
library(car)
leveneTest(value ~ age*genotype, data = cpn_count)

#######################################################
#######################################################
######### Creb3 mRNA in human single cell #############
#######################################################
#######################################################

dge_array <- lapply(dge_groups, function(cell) {
  message('*******************************************************************************')
  message(paste(Sys.time(), 'processing', cell))
  
  dir = paste0('processed_data/celltype_counts/edgeR_results_Condition_unfiltered/')
  
  file1 = read.table(file = paste0(dir,"edgeR_dge_results_ALSvsControl_",cell , '.txt'), header = T)
  file1$cell = cell 
  file1$group = "ALS"
  file1$fcsign = sign(file1$logFC)
  file1$metric= file1$mlog10PValue/file1$fcsign
  file3 = read.table(file = paste0(dir,"edgeR_dge_results_FTLDvsControl_",cell , '.txt'), header = T)
  file3$cell = cell 
  file3$group = "FTLD"
  file3$fcsign = sign(file3$logFC)
  file3$metric= file3$mlog10PValue/file3$fcsign
  
  tmp <- rbind(file1,file3)
  
})

final_deg <- do.call(rbind,dge_array)
final_deg$group <- factor(final_deg$group, levels = c( "ALS", "FTLD"))
final_deg_creb3 = final_deg[final_deg$gene_name == "CREB3" ,]


final_deg_creb3 = aggregate(final_deg_creb3[, 11], list(final_deg_creb3$cell,final_deg_creb3$group), mean)
colnames(final_deg_creb3) = c("cell","Group","Zscore")

to_plot_1 = final_deg_creb3


dend_order <- labels(dend.labeled)
order = as.data.frame(dend_order)
colnames(order)[1] = "cell"
order$position = row.names(order)

to_plot_1 = merge(to_plot_1,order,by=c("cell"))
to_plot_1= to_plot_1[order(as.numeric(as.character(to_plot_1$position))),]

to_plot_1$Group <- factor(to_plot_1$Group,levels = c("c9ALS","sALS", "c9FTLD", "sFTLD"))

p1 <- to_plot_1 %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = Group, fill = Zscore)) +
  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  
  scale_fill_gradientn(
    colors=c("blue","white","red"),
    values=rescale(c(-0.3,0.13,0.5)),
    limits=c(-0.3,0.5)
  ) +
  
  labs(fill = "Creb3_Zscore", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

pdf(p1,file = 'processed_data/figures_raw/creb3_zscore_hetmap.pdf')
p1
dev.off()
