
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
  library(plotrix)
  library(tidyverse)
  library(reshape2)
})
setwd('~/Desktop/single_cell/sc_als/')
data_path <- here('~/Desktop/single_cell/sc_als/')
here(data_path)

load(here(data_path,"processed_data/human_10x_10XALS_integration/Figure_5/merged_combined_rep.rda"))
load(here(data_path,"processed_data/human_10x_10XALS_integration/Figure_2/Cv3_tenxals.final.noSST.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/integrated_seurat/markers.cells.rda"))
markers_all = unique(c(markers_exc,markers_glia,markers_inh))


p1 <- DimPlot(merged_combined, reduction = "umap", group.by = "subclass_label", label = TRUE, label.size = 3,repel = TRUE,raster = F) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(merged_combined_rep, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE, raster = F) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2


merged_combined_rep@meta.data <- merged_combined_rep@meta.data %>% mutate(predicted.id_colors = case_when(
  predicted.id == "Astro" ~ "#8C6BB1",
  predicted.id == "Endo" ~ "#FDE0EF" ,
  predicted.id == "L2/3 IT" ~ "#C7E9C0",
  predicted.id == "L5 IT" ~ "#238B45",
  predicted.id == "L6 IT" ~ "#00441B",
  predicted.id == "L5 ET" ~ "#C6DBEF",
  predicted.id == "L5/6 NP" ~  "#9EBCDA" ,
  predicted.id == "L6 CT" ~ "#08306B",
  predicted.id == "L6b" ~  "#08519C",
  predicted.id == "L6 IT Car3" ~ "#9ECAE1" ,
  predicted.id == "Lamp5" ~ "#E5C494",
  predicted.id == "Vip" ~ "#A6761D",
  predicted.id == "Micro-PVM" ~ "#41B6C4",
  predicted.id == "VLMC" ~ "#FFFFD9",
  predicted.id == "OPC" ~ "#BDBDBD",
  predicted.id == "Oligo" ~ "#F781BF",
  predicted.id == "Pvalb" ~ "#EF3B2C",
  predicted.id == "Sst" ~ "#FB6A4A",
  predicted.id == "Chandelier" ~ "#FCBBA1"
))


colors = unique(merged_combined_rep@meta.data$predicted.id_colors)
names(colors) = unique(merged_combined_rep@meta.data$predicted.id)
DimPlot(merged_combined_rep, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3,repel = TRUE,raster = F, cols = colors) + NoLegend() + ggtitle("Reference annotations")


####### Build dendrogram for subclass label ############
anno = as.data.frame(merged_combined_rep@meta.data)

select.cells <- row.names(anno)
length(select.cells)

ref_anno <- anno %>%
  dplyr::filter(row.names(anno) %in% select.cells)

# Make a data.frame of unique cluster id, type, color, and broad type
ref.cl.df <- ref_anno %>%
  dplyr::select(predicted.id, 
                predicted.id_colors) %>%
  unique()

#standardize cluster annoation with cluster_id, cluster_label and cluster_color. These are the required fields to visualize clusters properly.
colnames(ref.cl.df)[1:2] <- c("cluster_label", "cluster_color")
row.names(ref.cl.df) <- ref.cl.df$cluster_label

# Set up the ref.cl factor object
ref.cl <- setNames(factor(ref_anno$predicted.id),row.names(ref_anno))


new.marker.genes = markers_all
sample.3.data_t.df = as.data.frame(merged_combined_rep@assays$RNA@data)

library(Matrix)
num.breaks <- (round(ncol(sample.3.data_t.df) / 10000))
print(num.breaks)
sample.3.data_t.df_tmp <- list()
d <- 1
while(d < num.breaks + 1){
  last.break <- d * 10000
  first.break <- last.break - 9999
  if(last.break < ncol(sample.3.data_t.df)){
    matrix.to.add <- as.matrix(sample.3.data_t.df[ , first.break:last.break])
    matrix.to.add <- Matrix(matrix.to.add, sparse = TRUE)
    matrix.to.add@x <- log2(matrix.to.add@x + 1)
    sample.3.data_t.df_tmp[[d]] <- matrix.to.add
  }
  if(last.break > ncol(sample.3.data_t.df)){
    matrix.to.add <- as.matrix(sample.3.data_t.df[ , first.break:ncol(sample.3.data_t.df)])
    matrix.to.add <- Matrix(matrix.to.add, sparse = TRUE)
    matrix.to.add@x <- log2(matrix.to.add@x + 1)
    sample.3.data_t.df_tmp[[d]] <- matrix.to.add
  }
  d <- d + 1
}
gc()

#Stitches subdivided log2 normalized chunks back together into single matrix
norm.dat <- sample.3.data_t.df_tmp[[1]]
d <- 2
while(d < num.breaks + 1){   
  norm.dat <- cbind(norm.dat, sample.3.data_t.df_tmp[[d]])
  d <- d + 1
}
gc()


library(scrattch.hicat)
library(scrattch.io)
library(clustree)
source('build_dend.R')

idx = markers_all %in% row.names(norm.dat)
markers_all = markers_all[idx]

cl.med <- get_cl_medians(norm.dat[markers_all,], 
                         ref.cl)
##The prefered order for the leaf nodes.
l.rank <- setNames(1:nrow(ref.cl.df), 
                   row.names(ref.cl.df))

##Color of the leaf nodes.
l.color <- setNames(as.character(ref.cl.df$cluster_color), row.names(ref.cl.df))

# Build the dendrogram

cl.dat = cl.med[,levels(ref.cl)]
dend.result <- build_dend(cl.med[,levels(ref.cl)],
                          l.rank, 
                          l.color,
                          nboot = 100)

dend <- dend.result$dend


###attach cluster labels to the leaves of the tree 
dend.labeled <- dend
labels(dend.labeled) <- ref.cl.df[labels(dend), "cluster_label"]
dend.labeled = dend.labeled %>% set("labels_col", l.color[labels(dend.labeled)])
dend.labeled = dend.labeled %>% set("leaves_col", l.color[labels(dend.labeled)])


pdf(file = "processed_data/figures_raw/dendrogram_sublcass_replication.pdf", width = 15, height = 7 )
dend.labeled %>%
  plot(main = "Sublcass taxonomy in motor and frontal cortex")
dev.off()


meta.sub_rep = merged_combined_rep@meta.data
cts <- AggregateExpression(merged_combined_rep,group.by=c('predicted.id', 'region', 'subject'), assays = 'RNA', slot = 'counts', return.seurat = FALSE)


save(cts,dend.labeled,meta.sub_rep, file = 'processed_data/human_10x_10XALS_integration/Figure_5/AggregatedCountsSublclass.rda')

load(here(data_path, "processed_data/human_10x_10XALS_integration/Figure_5/AggregatedCountsSublclass.rda"))


## transpose ###
cts = cts$RNA
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

cts.split.modified_final = lapply(cts.split.modified, as.data.frame)

vector_names = names(cts.split.modified_final)
vector_names = gsub('/', '-', vector_names)


cts.split.modified_final <- lapply(cts.split.modified_final, function(x){
  colnames(x) <- sub('medial frontal cortex', "medial_frontal_cortex", colnames(x))
  colnames(x) <- sub('motor cortex', "motor_cortex", colnames(x))
  x
})

meta.sub_rep$region = gsub(' ', '_', meta.sub_rep$region)

meta.sub_rep$Sample_ID = paste0(meta.sub_rep$region, '_', meta.sub_rep$subject)
meta.sub_rep$disease = gsub('-', '', meta.sub_rep$disease)


dge_groups = meta.sub_rep[!duplicated(meta.sub_rep$predicted.id),]

#whole$DGE_group = gsub("_", "-", whole$DGE_group)
dge_groups = dge_groups$predicted.id


#######################################################################
#######################################################################
#######################################################################
################ Write table DEG for the paper ########################
setwd('~/Desktop/Analysis_fus/LongGene_code/src/')
source('libraries.R')
source('ComBat-seq/ComBat_seq.R')
source('ComBat-seq/helper_seq.R')
setwd('~/Desktop/single_cell/sc_als/')

dge_array <- lapply(dge_groups, function(cell) {
  message('*******************************************************************************')
  message(paste(Sys.time(), 'processing', cell))
  
  dir = paste0('processed_data/human_10x_10XALS_integration/Figure_5/DGE_analysis/', cell)
  dir.create(file.path(dir), showWarnings = TRUE)
  
  counts = cts.split.modified_final[[cell]]
  edgeR_metadata = meta.sub_rep[!duplicated(meta.sub_rep$Sample_ID),]
  idx = colnames(counts) %in% edgeR_metadata$Sample_ID
  counts = counts[,idx]
  idx2 = edgeR_metadata$Sample_ID %in% colnames(counts)
  edgeR_metadata = edgeR_metadata[idx2,]
  edgeR_metadata$Group = paste0(edgeR_metadata$disease, '_', edgeR_metadata$region)
  edgeR_metadata$seq_batch = gsub(".*_","",edgeR_metadata$seq_batch)
  counts <- ComBat_seq(counts, batch=edgeR_metadata$seq_batch, group = edgeR_metadata$Group)
  cols <- c("Sample_ID","Group","sex")
  edgeR_metadata[cols] <- lapply(edgeR_metadata[cols], factor)  
  edgeR_metadata = edgeR_metadata[!duplicated(edgeR_metadata$Sample_ID),]
  rownames(edgeR_metadata) = edgeR_metadata$Sample_ID
  
  
  reorder_idx <- match(colnames(counts),row.names(edgeR_metadata))
  edgeR_metadata = edgeR_metadata[reorder_idx,]
  
  dds <- DESeqDataSetFromMatrix(round(counts), 
                                colData = edgeR_metadata, 
                                design = ~Group)
  relevel(dds$Group, ref = "Control_motor_cortex")
  vsd <- DESeq2::vst(dds, blind=FALSE)
  z <- DESeq2::plotPCA(vsd, intgroup=c("disease"))
  z + geom_label(aes(label = name))
  plot(z)

  
  
  ### Create DGE objects ####
  
  de.exon = as.DGEList(dds)
  de.exon$samples$Group = de.exon$samples$group
  
  keep <- filterByExpr(de.exon, group = de.exon$samples$Group)
  de.exon <- de.exon[keep,,keep.lib.sizes=FALSE]
  de.exon <- calcNormFactors(de.exon)
  design <- model.matrix(~0+Group+sex, de.exon$samples)
  de.exon <- estimateDisp(de.exon,design)
  log.exon <- edgeR::cpm(de.exon, log = T)
  
  norm.exon <- edgeR::cpm(de.exon, log = F)
  norm.exon =  as.data.frame(norm.exon)
  dim(norm.exon)
  
  log.exon = as.data.frame(log.exon)
  log.exon$gene_name = row.names(log.exon)
  dim(de.exon)
  
  #```{r edgeR-estimate-disp}
  ## Estimate dispersion and fit model
  qlfit <- glmQLFit(de.exon, design = design)
  
  # Define contrasts 
  
  #Before testing for differences in gene expression, we define the contrasts
  #we wish to test for. Here we represent the constrasts as a numeric matrix:
  
  
  #contrasts <- as.data.frame(makeContrasts(C9ALSvsControl=Groupc9ALS-GroupPN,C9FTLDvsControl=Groupc9FTLD-GroupPN,sFTLDvsControl=GroupsFTLD-GroupPN,sALSvsControl=GroupsALS-GroupPN, levels = design))
  #contrasts <- as.data.frame(makeContrasts(C9ALSvsControl=Groupc9ALS-GroupPN,C9FTLDvsControl=Groupc9FTLD-GroupPN,sFTLDvsControl=GroupsFTLD-GroupPN,sALSvsControl=GroupsALS-GroupPN, levels = design))
  contrasts <- as.data.frame(makeContrasts(ALSvsControl_motor=GroupC9ALS_motor_cortex-GroupControl_motor_cortex,FTDvsControl_motor=GroupC9FTD_motor_cortex-GroupControl_motor_cortex,ALSvsControl_frontal=GroupC9ALS_medial_frontal_cortex-GroupControl_medial_frontal_cortex,FTDvsControl_frontal=GroupC9FTD_medial_frontal_cortex-GroupControl_medial_frontal_cortex, levels = design))
  
  signif3 <- function(x) signif(x, digits = 3)
  edgeR_res <- lapply(contrasts, function(cm) {
    qlf <- glmQLFTest(qlfit, contrast = cm)
    tt <- topTags(qlf, n = Inf, sort.by = "none")$table
    tt %>%
      dplyr::mutate(mlog10PValue = -log10(PValue)) %>% 
      dplyr::mutate_at(vars(one_of(c("logFC", "logCPM", "F", 
                                     "PValue", "FDR", "mlog10PValue"))), 
                       list(signif3))
  })
  
  edgeR_res[[1]]$gene_name <- row.names(edgeR_res[[1]])
  edgeR_res[[2]]$gene_name <- row.names(edgeR_res[[2]])
  edgeR_res[[3]]$gene_name <- row.names(edgeR_res[[3]])
  edgeR_res[[4]]$gene_name <- row.names(edgeR_res[[4]])
  
  
  edgeR_res[[1]]$cell <- gsub(' ', '_', cell)
  edgeR_res[[2]]$cell <- gsub(' ', '_', cell)
  edgeR_res[[3]]$cell <- gsub(' ', '_', cell)
  edgeR_res[[4]]$cell <- gsub(' ', '_', cell)
  
  edgeR_res[[1]]$tissue <- "motor"
  edgeR_res[[2]]$tissue <- "motor"
  edgeR_res[[3]]$tissue <- "frontal"
  edgeR_res[[4]]$tissue <- "frontal"
  
  edgeR_res[[1]]$test <- "als_control"
  edgeR_res[[2]]$test <- "ftd_control"
  edgeR_res[[3]]$test <- "als_control"
  edgeR_res[[4]]$test <- "ftd_control"
  
  
  
  for (nm in names(edgeR_res)) {
    write.table(edgeR_res[[nm]] %>% dplyr::arrange(PValue), 
                file = paste0(dir, '/', "edgeR_dge_results_", nm, '_', cell, ".txt"), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  # tagwise = de.exon$tagwise.dispersion
  # counts = de.exon$AveLogCPM
  # 
  # df = data_frame(tagwise,counts)
  # df$cell = cell
  # write.table(df, 
  #             file = paste0(dir, '/', "tagwise_dispersion", '_', cell, ".txt"), 
  #             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
})


creb3.bound = read.table('processed_data/human_10x_10XALS_integration/Figure_3/CREB3_bound.genes', header = T)

###########
dge_array <- lapply(dge_groups, function(cell) {
  message('*******************************************************************************')
  message(paste(Sys.time(), 'processing', cell))
  
  dir = paste0('processed_data/human_10x_10XALS_integration/Figure_5/DGE_analysis/files_degs/')
  
  file1 = read.table(file = paste0(dir,"edgeR_dge_results_ALSvsControl_frontal_",cell , '.txt'), header = T)
  file1$group = "ALS"
  file1$fcsign = sign(file1$logFC)
  file1$metric= file1$mlog10PValue/file1$fcsign
  file2 = read.table(file = paste0(dir,"edgeR_dge_results_ALSvsControl_motor_",cell , '.txt'), header = T)
  file2$group = "ALS"
  file2$fcsign = sign(file2$logFC)
  file2$metric= file2$mlog10PValue/file2$fcsign
  file3 = read.table(file = paste0(dir,"edgeR_dge_results_FTDvsControl_frontal_",cell , '.txt'), header = T)
  file3$group = "FTLD"
  file3$fcsign = sign(file3$logFC)
  file3$metric= file3$mlog10PValue/file3$fcsign
  file4 = read.table(file = paste0(dir,"edgeR_dge_results_FTDvsControl_motor_",cell , '.txt'), header = T)
  file4$group = "FTLD"
  file4$fcsign = sign(file4$logFC)
  file4$metric= file4$mlog10PValue/file4$fcsign
  
  tmp <- rbind(file1,file2,file3,file4)
  
})

final_deg <- do.call(rbind,dge_array)
final_deg$group <- factor(final_deg$group, levels = c( "ALS", "FTLD"))
final_deg_creb3 = final_deg[final_deg$gene_name %in% creb3.bound$Gene,]


final_deg_creb3 = aggregate(final_deg_creb3[, 13], list(final_deg_creb3$cell,final_deg_creb3$tissue,final_deg_creb3$group), mean)
colnames(final_deg_creb3) = c("cell","tissue","group","Zscore")

to_plot_1 = final_deg_creb3


dend_order <- labels(dend.labeled)
order = as.data.frame(dend_order)
colnames(order)[1] = "cell"
order$position = row.names(order)
order$cell = gsub(' ', '_', order$cell)


to_plot_1 = merge(to_plot_1,order,by=c("cell"))
to_plot_1= to_plot_1[order(as.numeric(as.character(to_plot_1$position))),]
to_plot_1$group <- factor(to_plot_1$group,levels = c("ALS", "FTLD"))
to_plot_1$Zscore = to_plot_1$Zscore * 10

to_plot_1_motor = to_plot_1[to_plot_1$tissue == "motor",]
to_plot_1_frontal = to_plot_1[to_plot_1$tissue != "motor",]


p1 <- to_plot_1_motor %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = group, fill = Zscore)) +
  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  
#  scale_fill_gradientn(
 #   colors=c("blue","white","red"),
 #   values=rescale(c(-0.3,0.13,0.5)),
 #   limits=c(-0.3,0.5)
 # ) +
  
  labs(fill = "Creb3_Zscore", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

p2 <- to_plot_1_frontal %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = group, fill = Zscore)) +
  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  
  #  scale_fill_gradientn(
  #   colors=c("blue","white","red"),
  #   values=rescale(c(-0.3,0.13,0.5)),
  #   limits=c(-0.3,0.5)
  # ) +
  
  labs(fill = "Creb3_Zscore", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())


pdf(p1,file = 'processed_data/figures_raw/creb3_zscore_heatmap_replication.pdf')
p1+p2
dev.off()




####################################################### 
#######################################################
######### Creb3 mRNA in human single cell #############
#######################################################
#######################################################
final_deg_creb3 = final_deg[final_deg$gene_name %in% creb3.bound$Gene,]

motor_creb3 = final_deg_creb3[final_deg_creb3$tissue == "motor",]
frontal_creb3 = final_deg_creb3[final_deg_creb3$tissue != "motor",]

##Compute the analysis of variance
res.aov_motor <- aov(metric ~ cell*group, data = motor_creb3)
summary(res.aov_motor)
pwc_motor = TukeyHSD(res.aov_motor, which = "cell:group")
pwc_motor = as.data.frame(pwc_motor[[1]])
pwc_motor$comp = row.names(pwc_motor)
pwc_motor = pwc_motor[order(pwc_motor$comp),]

comp <- data.frame(do.call('rbind', strsplit(as.character(pwc_motor$comp),':',fixed=TRUE)))
comp$X2 = sub("-", ":", comp$X2)
comp2 <- data.frame(do.call('rbind', strsplit(as.character(comp$X2),':',fixed=TRUE)))
comp3 = cbind(comp,comp2)
comp3 = comp3[,c(1,4,5,3)]
colnames(comp3) <- c("cell1", "group1", "cell2", "group2")
pwc_motor = cbind(pwc_motor,comp3)

idx = (pwc_motor$cell1 == pwc_motor$cell2 & pwc_motor$group1 != pwc_motor$group2)
pwc_motor = pwc_motor[idx,]
pwc_motor$tissue = "motor_cortex"

### Anova cortex frontal 
res.aov_frontal <- aov(metric ~ cell*group, data =frontal_creb3)
summary(res.aov_frontal)
pwc_frontal = TukeyHSD(res.aov_frontal, which = "cell:group")
pwc_frontal = pwc_frontal[[1]]
pwc_frontal = as.data.frame(pwc_frontal)
pwc_frontal$comp = row.names(pwc_frontal)
pwc_frontal = pwc_frontal[order(pwc_frontal$comp),]

comp <- data.frame(do.call('rbind', strsplit(as.character(pwc_frontal$comp),':',fixed=TRUE)))
comp$X2 = sub("-", ":", comp$X2)
comp2 <- data.frame(do.call('rbind', strsplit(as.character(comp$X2),':',fixed=TRUE)))
comp3 = cbind(comp,comp2)
comp3 = comp3[,c(1,4,5,3)]
colnames(comp3) <- c("cell1", "group1", "cell2", "group2")
pwc_frontal = cbind(pwc_frontal,comp3)

idx = (pwc_frontal$cell1 == pwc_frontal$cell2 & pwc_frontal$group1 != pwc_frontal$group2)
pwc_frontal = pwc_frontal[idx,]
pwc_frontal$tissue = "frontal_cortex"

all_stats_motor_frontal = rbind(pwc_motor,pwc_frontal)

write.table(all_stats_motor_frontal,file = "processed_data/article_marques/stats_heatmap_crosstissue.txt", quote = F, row.names = F, col.names = T, sep = "\t")


dodge <- position_dodge(width = 0.4)
deg_file.mouse.csn$time <- factor(deg_file.mouse.csn$time, levels = c( "pre", "post"))

p<-ggplot(deg_file.mouse.csn, aes(x=time, y=metric, fill=cell))
p + geom_violin(
  aes(color = cell), trim = FALSE,
  position = dodge
) +
  geom_boxplot(width=.1, outlier.colour=NA, position = dodge)
dev.off()



#######################################################################
#######################################################################
#######################################################################
############## CREB3 target genes exp and survival ####################
#######################################################################
#######################################################################
#######################################################################

# load(here(data_path, "processed_data/human_10x_10XALS_integration/Figure_5/AggregatedCountsSublclass.rda"))
# ## transpose ###
# cts = cts$RNA
# cts.t = t(cts)
# 
# ### convert to dataframe ###
# cts.t = as.data.frame(cts.t)
# cts.t[1:10,1:10]
# 
# ### get value to split !!!!!
# splitRows <- gsub('_.*', '',row.names(cts.t))
# 
# cts.split <- split.data.frame(cts.t, f = factor(splitRows))
# 
# cts.split.modified <- lapply(cts.split, function(x){
#   rownames(x) <- sub("^[^_]*_", "", rownames(x))
#   t(x)
# })
# 
# cts.split.modified_final = lapply(cts.split.modified, as.data.frame)
# 
# names(cts.split.modified_final) = gsub(' ', '-', names(cts.split.modified_final))
# 
# 
# cts.split.modified_final <- lapply(cts.split.modified_final, function(x){
#   colnames(x) <- sub('medial frontal cortex', "medial_frontal_cortex", colnames(x))
#   colnames(x) <- sub('motor cortex', "motor_cortex", colnames(x))
#   x
# })
# 
# meta.sub_rep$region = gsub(' ', '_', meta.sub_rep$region)
# 
# meta.sub_rep$Sample_ID = paste0(meta.sub_rep$region, '_', meta.sub_rep$subject)
# meta.sub_rep$disease = gsub('-', '', meta.sub_rep$disease)
# 
# 
# dge_groups = meta.sub_rep[!duplicated(meta.sub_rep$predicted.id),]
# 
# #whole$DGE_group = gsub("_", "-", whole$DGE_group)
# dge_groups = dge_groups$predicted.id
# dge_groups = gsub(' ', '-', dge_groups)
# 
# turquoise_genes = read.table('processed_data/article_marques/gene_list_modules.txt', header = T)
# turquoise_genes = turquoise_genes[turquoise_genes$module == "turquoise",]
# 
# 
# 
# dge_array <- lapply(dge_groups, function(cell) {
#   message('*******************************************************************************')
# message(paste(Sys.time(), 'processing', cell))
# 
# dir = paste0('processed_data/human_10x_10XALS_integration/Figure_5/CREB3_target_progression/', cell)
# dir.create(file.path(dir), showWarnings = TRUE)
# 
# counts = cts.split.modified_final[[cell]]
# edgeR_metadata = meta.sub_rep[!duplicated(meta.sub_rep$Sample_ID),]
# idx = colnames(counts) %in% edgeR_metadata$Sample_ID
# counts = counts[,idx]
# idx2 = edgeR_metadata$Sample_ID %in% colnames(counts)
# edgeR_metadata = edgeR_metadata[idx2,]
# edgeR_metadata$Group = paste0(edgeR_metadata$disease, '_', edgeR_metadata$region)
# edgeR_metadata$seq_batch = gsub(".*_","",edgeR_metadata$seq_batch)
# counts <- ComBat_seq(counts, batch=edgeR_metadata$seq_batch, group = edgeR_metadata$sex)
# cols <- c("Sample_ID","Group","sex")
# edgeR_metadata[cols] <- lapply(edgeR_metadata[cols], factor)  
# edgeR_metadata = edgeR_metadata[!duplicated(edgeR_metadata$Sample_ID),]
# rownames(edgeR_metadata) = edgeR_metadata$Sample_ID
# 
# survival = read.table('processed_data/human_10x_10XALS_integration/Figure_5/clinical_data_singlecell_final.txt', header = T)
# edgeR_metadata = merge(edgeR_metadata,survival,by=c("subject"))
# row.names(edgeR_metadata) = edgeR_metadata$Sample_ID
# 
# counts = counts[,colnames(counts) %in% row.names(edgeR_metadata)]
# 
# 
# reorder_idx <- match(colnames(counts),row.names(edgeR_metadata))
# edgeR_metadata = edgeR_metadata[reorder_idx,]
# edgeR_metadata$Disease_Duration = edgeR_metadata$Disease_Duration * 12
# 
# 
# dds <- DESeqDataSetFromMatrix(round(counts), 
#                               colData = edgeR_metadata, 
#                               design = ~Group)
# 
# vsd <- vst(dds, blind=FALSE)
# z <- plotPCA(vsd, intgroup=c("sex"))
# z + geom_label(aes(label = name))
# 
# 
# de.exon = as.DGEList(dds)
# de.exon$samples$Group = de.exon$samples$group
# 
# keep <- filterByExpr(de.exon, group = de.exon$samples$Group)
# de.exon <- de.exon[keep,,keep.lib.sizes=FALSE]
# de.exon <- calcNormFactors(de.exon)
# design <- model.matrix(~0+Group+sex, de.exon$samples)
# de.exon <- estimateDisp(de.exon,design)
# log.exon <- edgeR::cpm(de.exon, log = T)
# 
# log.exon = log.exon[row.names(log.exon) %in% creb3.bound$Gene,]
# #log.exon = log.exon[row.names(log.exon) %in% turquoise_genes$gene,]
# 
# 
# matrix = t(scale(t(log.exon)))
# 
# 
# log.exon_melt = reshape2::melt(as.matrix(matrix))
# colnames(log.exon_melt) = c("gene", "Sample_ID", "counts")
# log.exon_melt = merge(log.exon_melt,edgeR_metadata,by=c("Sample_ID"))
# log.exon_melt = log.exon_melt[,c(1,2,3,4,8,9,10,46,47,48,49)]
# 
# log.exon_melt = log.exon_melt[log.exon_melt$Sample_ID != "motor_cortex_F2" & log.exon_melt$subject != "motor_cortex_F1" & log.exon_melt$subject != "motor_cortex_F5" & log.exon_melt$subject != "medial_frontal_cortex_F3" ,]
# 
# x = do.call(data.frame,aggregate(log.exon_melt[, 3], list(log.exon_melt$subject,log.exon_melt$region), function(x) c(mean = median(x), sd = sd(x))))
# colnames(x) <- c("subject", "tissue", "mean_Z", "se_Z")
# x = x[x$subject != "F2",]
# x$cell = cell
# 
# y = log.exon_melt[,c(4,9,10,11)]
# y = y[!duplicated(y$subject),]
# z = merge(x,y,by=c('subject'))
# 
# write.table(z, 
#             file = paste0(dir, "/", "CREB3_targetGenes", '_', cell, ".txt"), 
#             sep = "\t", row.names = T, col.names = TRUE, quote = FALSE)
# 
# })
# 
# 
# merged_bcv <- 
#   do.call(rbind,
#           lapply(list.files(path = "processed_data/human_10x_10XALS_integration/Figure_5/CREB3_target_progression/creb3_files/", recursive = T, pattern = "\\.txt$", full.names = T), read.table, header=T, row.names=NULL))
# 
# 
# colors = meta.sub_rep[!duplicated(meta.sub_rep$predicted.id),]
# colors = colors[,c("predicted.id", "predicted.id_colors")]
# colnames(colors) = c("cell", "colors")
# colors$cell = gsub(' ', '-', colors$cell)
# 
# merged_bcv = merge(merged_bcv,colors,by=c("cell"))
# 
# 
# merged_bcv <- merged_bcv %>% mutate(class = case_when(
#   cell == "Astro" ~ "Glia",
#   cell == "Endo" ~ "Vascular",
#   cell == "Micro-PVM" ~ "Glia",
#   cell == "Oligo" ~ "Glia",
#   cell == "OPC" ~ "Glia",
#   cell == "VLMC" ~ "Vascular",
#   cell == "Chandelier" ~ "Inhibitory",
#   cell == "Sst" ~ "Inhibitory",
#   cell == "Vip" ~ "Inhibitory",
#   cell == "Pvalb" ~ "Inhibitory",
#   cell == "Lamp5" ~ "Inhibitory",
#   cell == "L2-3-IT" ~ "Excitatory",
#   cell == "L5-6-NP" ~ "Excitatory",
#   cell == "L5-ET" ~ "Excitatory",
#   cell == "L5-IT" ~ "Excitatory",
#   cell == "L6-IT" ~ "Excitatory",
#   cell == "L6-CT" ~ "Excitatory",
#   cell == "L6b" ~ "Excitatory",
#   cell == "L6-IT-Car3" ~ "Excitatory",
# 
# ))
# 
# 
# library(plyr)
# corfun<-function(x, y) {
#   corr=(cor.test(x, y,
#                  alternative="two.sided", method="pearson"))
# }
# cor_table = ddply(merged_bcv, .(class), summarise,z=corfun(mean_Z,Disease_Duration)$statistic,
#                   pval=corfun(mean_Z,Disease_Duration)$p.value,
#                   tau.est=corfun(mean_Z,Disease_Duration)$estimate,
#                   alt=corfun(mean_Z,Disease_Duration)$alternative
# ) 
# 
# 
# 
# ggplot(merged_bcv, aes(x=Disease_Duration, y=mean_Z)) + geom_point()
# 
# ggplot(mtcars, aes(x=wt, y=mpg, color=cyl, shape=cyl)) +
#   geom_point() + 
#   geom_smooth(method=lm, se=T, fullrange=TRUE)+
#   scale_shape_manual(values=c(3, 16, 17))+ 
#   scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
#   theme(legend.position="top")
# 
# merged_bcv = aggregate(merged_bcv$tagwise, list(merged_bcv$cell,merged_bcv$new_cpm), FUN=mean) 
# colnames(avg_bcv) = c("cluster", "bin", "dispersion")
# 
# 
# 
# 
# 
