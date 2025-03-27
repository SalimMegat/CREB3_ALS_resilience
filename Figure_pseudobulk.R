
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
  library("scProportionTest")
  library(tidyverse)
  library(feather)
})
setwd('~/Desktop/single_cell/sc_als/')
data_path <- ('~/Desktop/single_cell/sc_als/')

load(paste0(data_path,"processed_data/human_10x_10XALS_integration/Figure_2/Cv3_tenxals.final.noSST.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/integrated_seurat/markers.cells.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/Figure_2/AggregatedCountsCombat.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/Figure_1/aucell_metadata.3.rda"))
markers_all = unique(c(markers_exc,markers_glia,markers_inh))

# meta = na.omit(merged_combined@meta.data)
# meta = meta[,c(1,17,18,19,20,21,22,16)]
# meta$Barcode = gsub("\\-.*","",meta$Barcode)
# meta$Barcode = paste0(meta$Barcode,"-",1)
# write.table(meta, file = "~/Desktop/whole_metadata_new",quote = F,sep = "\t",row.names = F,col.names = T)


#Plot 1b: CSN and CPN heatmap
#check all Cv3 clusters are represented
dend_order <- labels(dend.labeled)
order.dendo = as.data.frame(dend_order)
colnames(order.dendo)[1] = "cluster"
order.dendo$position = row.names(order.dendo)

pdf(file = "processed_data/figures_raw/dendrogram_degGroup.pdf", width = 16, height = 7 )
dend.labeled %>%
  rotate(as.character(dend_order)) %>% #rotate to match labels new order
  plot(main = "Motor cortex taxonomy human")
dev.off()



als <- merged_combined@meta.data
aucell_tmp = aucell[,c(16,23)]
als = merge(als,aucell_tmp, by=c("row.names"))
row.names(als) = als$Row.names

cv3_cl <- als %>%
  filter(orig.ident  %>% str_detect("tenxals")) %>%
  distinct(cluster_label.dendro.deg) %>% pull() 
dend_order[which(dend_order %in% cv3_cl == FALSE)]

#make plotting object
to_plot_3 <- als %>%
  filter(orig.ident %>% str_detect("tenxals")) %>%
  
  group_by(cluster_label.dendro.deg) %>%
  mutate(cluster_n = n()) %>%
  ungroup() %>%
  
  group_by(cluster_label.dendro.deg, mouse_facs) %>%
  mutate(cluster_facs_n = n()) %>%
  ungroup() %>%
  
  mutate(facs_prop = cluster_facs_n / cluster_n) %>%
  
  distinct(cluster_label.dendro.deg, mouse_facs, facs_prop) %>%
  
  #factorize
  mutate(facs = mouse_facs %>% as_factor() %>% fct_relevel(c("0", "1", "2")) %>% fct_rev(),
         cluster = cluster_label.dendro.deg %>% as_factor() %>% fct_relevel(cv3_cl))




to_plot_3 = merge(to_plot_3,order.dendo,by=c("cluster"))
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

pdf(p1b,file = "processed_data/figures_raw/Heatmap_csn_deg.pdf")
p1b
dev.off()

#Plot 2: proportions of DEG cell group per disease group
prop_test <- sc_utils(merged_combined)

prop_test.subclass.sALS <- permutation_test(
  prop_test, cluster_identity = "cluster_label.dendro.deg",
  sample_1 = "PN", sample_2 = "sALS",
  sample_identity = "Group"
)

prop_test.subclass.sFTLD <- permutation_test(
  prop_test, cluster_identity = "cluster_label.dendro.deg",
  sample_1 = "PN", sample_2 = "sFTLD",
  sample_identity = "Group"
)

prop_test.subclass.c9ALS <- permutation_test(
  prop_test, cluster_identity = "cluster_label.dendro.deg",
  sample_1 = "PN", sample_2 = "c9ALS",
  sample_identity = "Group"
)

prop_test.subclass.c9FTLD <- permutation_test(
  prop_test, cluster_identity = "cluster_label.dendro.deg",
  sample_1 = "PN", sample_2 = "c9FTLD",
  sample_identity = "Group"
)


permutation_plot(prop_test.subclass.sALS, FDR_threshold = 0.01, log2FD_threshold = 0.58)
permutation_plot(prop_test.subclass.sFTLD, FDR_threshold = 0.01, log2FD_threshold = 0.58)
permutation_plot(prop_test.subclass.c9ALS, FDR_threshold = 0.01, log2FD_threshold = 0.58)
permutation_plot(prop_test.subclass.c9FTLD, FDR_threshold = 0.01, log2FD_threshold = 0.58)


tmp1 = prop_test.subclass.c9ALS@results$permutation
tmp1$Group = "c9ALS"
tmp1 = tmp1[,-c(2,3)]
tmp2 = prop_test.subclass.c9FTLD@results$permutation
tmp2$Group = "c9FTLD"
tmp2 = tmp2[,-c(2,3)]
tmp3 = prop_test.subclass.sFTLD@results$permutation
tmp3$Group = "sFTLD"
tmp3 = tmp3[,-c(2,3)]
tmp4 = prop_test.subclass.sALS@results$permutation
tmp4$Group = "sALS"
tmp4 = tmp4[,-c(2,3)]


to_plot_5 <- rbind(tmp1,tmp2,tmp3,tmp4)
colnames(to_plot_5)[1] <- "cluster"

order = as.data.frame(dend_order)
colnames(order)[1] = "cluster"
order$position = row.names(order)

to_plot_5 = merge(to_plot_5,order,by=c("cluster"))


tmp <- als %>% arrange(cluster_label.dendro.deg) %>% distinct(cluster_label.dendro.deg, cluster_label.dendro.deg.colors)
colnames(tmp)[1] <- "cluster"

to_plot_5 = merge(to_plot_5,tmp,by=c("cluster"))
to_plot_5= to_plot_5[order(as.numeric(as.character(to_plot_5$position))),]
colors_use <- to_plot_5$cluster_label.dendro.deg.colors
names(colors_use) <- to_plot_5$cluster_label.dendro.deg.colors


p3 <- to_plot_5 %>%
  ggplot() +
  geom_point(aes(x = reorder(position, sort(as.numeric(position))), y = obs_log2FD, color = colors_use), size = 2, alpha = 0.6) +
  scale_color_manual(values = colors_use) +
  
  labs(y = "Cell Proportion (Control/ALS_FTD)",
       x = "") +  
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        legend.position = "none",
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 1, linetype = 2),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

p3
pdf(p3,file = "processed_data/figures_raw/cell_proportion.pdf")
p3
dev.off()

#Plot 3: transcriptome wide shift heatmap ######
dge_groups = meta.sub[!duplicated(meta.sub$cluster_label.dendro.deg),]
dge_groups = dge_groups$cluster_label.dendro.deg

dir = 'processed_data/celltype_counts/pca_distance/'
files = list.files('processed_data/celltype_counts/pca_distance/')
big.list.of.data.frames <- lapply(paste0(dir,'/' , files), read.table, 
                                  header = TRUE,
                                  stringsAsFactors = FALSE)
big.data.frame <- do.call(rbind,big.list.of.data.frames)

title = as.data.frame(files)
title$files = gsub("transcriptome_wide_shift_", "",title$files)
title$files = gsub(".txt", "",title$files)
title = as.data.frame(title[rep(seq_len(nrow(title)), each = 4), ])
colnames(title) <- 'cluster'
to_plot_4 = cbind(title,big.data.frame)

order = as.data.frame(dend_order)
colnames(order)[1] = "cluster"
order$position = row.names(order)

to_plot_4 = merge(to_plot_4,order,by=c("cluster"))
to_plot_4= to_plot_4[order(as.numeric(as.character(to_plot_4$position))),]


p4 <- to_plot_4 %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = Group, fill = NormalizedDistance)) +
  
  #scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  
  scale_fill_gradientn(
    colors=c("blue","white","red"),
    values=rescale(c(-1,0,1)),
    limits=c(-2,2)
  ) +
  
  labs(fill = "Normalized Genome-Wide Distance", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

p4

pdf(p4,file = "processed_data/figures_raw/transcriptome_wide_distance.pdf")
p4
dev.off()

#Plot 4: CSN DE gene overlap with human DE genes ######
csn = read.table('aucell/csn_gene_human', header = F, stringsAsFactors = F)
cpn = read.table('aucell/cpn_gene_human', header = F, stringsAsFactors = F)
csn = as.character(as.factor(csn$V1))
cpn = as.character(as.factor(cpn$V1))

files <- list.files(path = "processed_data/celltype_counts/edgeR_results_Condition/", full.names = T, recursive = F,pattern = "*.txt")
names(files) <- dir('processed_data/celltype_counts/edgeR_results_Condition/', '*.txt$')

pvalues_csn <- vector()
enrichment_csn <- vector()
pvalues_csn.up <- vector()
pvalues_csn.dn <- vector()

pvalues_cpn <- vector()
enrichment_cpn <- vector()
enrichment_csn.up <- vector()
enrichment_csn.dn <- vector()


for (i in 1:60) {
  message('*******************************************************************************')
  message(paste(Sys.time(), 'processing', i))
  
  single_cell = read.table(file = files[i], header = T)
  single_cell = single_cell[single_cell$FDR< 0.05,]

  single_cell.up = single_cell[single_cell$FDR< 0.05 & single_cell$logFC > 0,]
  single_cell.dn = single_cell[single_cell$FDR< 0.05 & single_cell$logFC < 0,]
  
  
    
  if (nrow(single_cell)>0){
    
    csn_genes = sum(single_cell$gene_name %in% csn == "TRUE")
    cpn_genes = sum(single_cell$gene_name %in% cpn == "TRUE")

    csn_genes.up = sum(single_cell.up$gene_name %in% csn == "TRUE")
    csn_genes.dn = sum(single_cell.dn$gene_name %in% csn == "TRUE")
    
    
    fold.enrichment.up <-  (csn_genes.up / nrow(single_cell.up) ) / (324 / 13325)
    fold.enrichment.dn <-  (csn_genes.dn / nrow(single_cell.dn) ) / (324 / 13325)
    
    
    p.values1.up <- phyper(q=csn_genes.up-1, m=324, n=13325, k=nrow(single_cell.up), lower.tail=F)
    p.values2.dn <- phyper(q=csn_genes.dn-1, m=324, n=13325, k=nrow(single_cell.dn), lower.tail=F)
    
    
    pvalues_csn.up <- c(pvalues_csn.up, p.values1.up)
    pvalues_csn.dn <- c(pvalues_csn.dn, p.values2.dn)
  
    
    enrichment_csn.up <- c(enrichment_csn.up, fold.enrichment.up)
    enrichment_csn.dn <- c(enrichment_csn.dn, fold.enrichment.dn)
    
  }
  
}

csn_matrix.up = as.data.frame(cbind(enrichment_csn.up, pvalues_csn.up))
csn_matrix.dn = as.data.frame(cbind(enrichment_csn.dn, pvalues_csn.dn))

names = gsub("edgeR_dge_results_","", names(files))
names = gsub("*.txt","",names)
names = gsub("*vsControl","",names)

csn_matrix.up$cell_type_condition <- names
csn_matrix.dn$cell_type_condition <- names

library(tidyverse)
library(ggpubr)

csn_matrix.up = csn_matrix.up %>% 
  separate(cell_type_condition, into = c("Group", "CellType"), sep="_(?=[^_]+$)")
csn_matrix.dn = csn_matrix.dn %>% 
  separate(cell_type_condition, into = c("Group", "CellType"), sep="_(?=[^_]+$)")

#csn_matrix$pvalues_csn = -log(csn_matrix$pvalues_csn)
#cpn_matrix$pvalues_cpn = -log(cpn_matrix$pvalues_cpn)
colnames(csn_matrix.up)[4] <- "cluster"
colnames(csn_matrix.dn)[4] <- "cluster"

colnames(csn_matrix.dn) <- c("enrichment_csn","pvalues_csn","Group","cluster","direction")
colnames(csn_matrix.up) <- c("enrichment_csn","pvalues_csn","Group","cluster","direction")

csn_matrix.dn$direction = 'down-regulated_in_mouseCSN'
csn_matrix.up$direction = 'up-regulated_in_mouseCSN'

combined = rbind(csn_matrix.up,csn_matrix.dn)

write.table(combined,file = 'processed_data/article_marques/DEG_enrichment_crossspecies.txt', quote = F, row.names = F, col.names = T)

csn_matrix.up$z_score = qnorm(csn_matrix.up$pvalues_csn,lower.tail = F)
is.na(csn_matrix.up)<-sapply(csn_matrix.up, is.infinite)
csn_matrix.up[is.na(csn_matrix.up)]<-0

csn_matrix.dn$z_score = qnorm(csn_matrix.dn$pvalues_csn,lower.tail = F)
is.na(csn_matrix.dn)<-sapply(csn_matrix.dn, is.infinite)
csn_matrix.dn[is.na(csn_matrix.dn)]<-0

deg_cell = meta.sub[!duplicated(meta.sub$cluster_label.dendro.deg),]
deg_cell = deg_cell[,c("cluster_label.dendro.deg" ,"class_label")]
colnames(deg_cell)[1]  = "cluster"

csn_matrix.up = merge(csn_matrix.up,deg_cell,by=c('cluster'))

res.aov_csn_matrix.up <- aov(z_score ~ class_label, data = csn_matrix.up)
summary(res.aov_csn_matrix.up)
pwc <- TukeyHSD(res.aov_csn_matrix.up, which = "class_label")
pwc <- as.data.frame(pwc[[1]])


write.table(pwc,file = 'processed_data/article_marques/post-hoc_DEG_enrichment_crossspecies.txt', quote = F, row.names = F, col.names = T)


min = min(csn_matrix.up$z_score) - 0.1
max = max(csn_matrix.up$z_score) + 0.1

pdf("~/Desktop/single_cell/sc_als/processed_data/figures_raw/Creb3_luciferase.pdf")
csn_matrix.up %>%
  mutate(class_label = factor(class_label, levels = c("Glutamatergic","GABAergic","Non-Neuronal"))) %>%
  ggplot(aes(fill=class_label, y=z_score, x=class_label)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Z-score CSN-ALS signature enrichement in human neurons") +
  ylim(min,max)
dev.off()



to_plot_2 = merge(csn_matrix.up,order.dendo,by=c("cluster"))
to_plot_2= to_plot_2[order(as.numeric(as.character(to_plot_2$position))),]


p2 <- to_plot_2 %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = Group, fill = z_score)) +
  
  #scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
   scale_fill_gradientn(
    colors=c("blue","white","red"),
    values=rescale(c(-3,3,10)),
    limits=c(-3,10)
  ) +
  
  
  labs(fill = "CSN DEG enrichement in human ALS", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

p2
pdf(p2,file = "processed_data/figures_raw/csn_deg_enrichment.pdf")
p2
dev.off()


##### Calculate cumulative rank based on proportion and overlap DEG #########
overlap = to_plot_2
overlap = overlap[overlap$Group != "FTLD",]
overlap = overlap[,c(1,5)]

proportion = to_plot_3
proportion = proportion[!duplicated(proportion$cluster_label.dendro.deg),]
proportion = proportion[,c(1,4)]

final_tab = merge(overlap,proportion)
final_tab$facs_prop = final_tab$facs_prop * 10
final_tab$cumsum = final_tab$z_score + final_tab$facs_prop

to_plot_6 = merge(final_tab,order.dendo,by=c("cluster"))
to_plot_6= to_plot_6[order(as.numeric(as.character(to_plot_6$position))),]
to_plot_6$Rank = c("Cumulative Rank")

p6 <- to_plot_6 %>%
  ggplot() +
  geom_tile(aes(x = reorder(position, sort(as.numeric(position))), y = Rank, fill = cumsum)) +
  
  #scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  scale_fill_gradientn(
    colors=c("blue","white","red"),
    values=rescale(c(-3,3,10)),
    limits=c(-3,11)
  ) +
  
  
  labs(fill = "CSN DEG enrichement in human ALS", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 6/length(dend_order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())

p6
pdf(p6,file = "processed_data/figures_raw/cumsum.pdf")
p6
dev.off()




####### Violin plot of NEFH expression across neurons #######
aucell_tmp2 = aucell[,c(5,23)]
merged_combined_added <- AddMetaData(merged_combined, metadata = aucell_tmp2)
merged_combined_added@meta.data$mouse_facs.colors = ifelse(merged_combined_added@meta.data$mouse_facs == "0", "#D9D9D9",ifelse(merged_combined_added@meta.data$mouse_facs == "1","#6A51A3", "#41AB5D" ))


# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(merged_combined_added, features = "NEFH", ncol = 2)

########
colors.cluster = unique(merged_combined_added@meta.data$cluster_label.dendro.deg.colors)
names(colors.cluster) = unique(merged_combined_added@meta.data$cluster_label.dendro.deg)

colors.facs = unique(merged_combined_added@meta.data$mouse_facs.colors)
names(colors.facs) = unique(merged_combined_added@meta.data$mouse_facs)

p3 <- FeaturePlot(merged_combined, features = "NEFH",raster = F)
p1 <- DimPlot(merged_combined_added, reduction = "umap", group.by = "cluster_label.dendro.deg", cols = colors.cluster , raster = F)  + NoLegend()
p2 <- DimPlot(merged_combined_added, reduction = "umap", group.by = "mouse_facs", cols = colors.facs , raster = F)  + NoLegend()

library(ggpubr)
pdf(file = "processed_data/figures_raw/umap_clusterNEFH.pdf", width = 11, height = 11)
ggarrange(p1, p2, p3, ncol = 2, nrow = 2)
dev.off()

####### 
gene_exp = read.table('processed_data/celltype_counts/allgenes_countsincell.txt', header = T)
gene_exp = aggregate(gene_exp$logCPM, list(gene_exp$cluster,gene_exp$gene), FUN=mean) 
colnames(gene_exp) = c("cluster", 'gene','logCPM')

gene_exp = merge(gene_exp,to_plot_3,by=c("cluster"))
gene_exp = gene_exp[gene_exp$facs == "2",]


tmp <- als %>% arrange(cluster_label.dendro.deg) %>% distinct(cluster_label.dendro.deg, cluster_label.dendro.deg.colors)
colnames(tmp)[1] <- "cluster"

gene_exp = merge(gene_exp,tmp,by=c("cluster"))
gene_exp= gene_exp[order(as.numeric(as.character(gene_exp$position))),]
colors_use <- gene_exp$cluster_label.dendro.deg.colors
names(colors_use) <- gene_exp$cluster_label.dendro.deg.colors

require(plyr)
func <- function(nefh_exp)
{
  return(data.frame(COR = cor(nefh_exp$logCPM, nefh_exp$facs_prop)))
}

cor_matrix = ddply(gene_exp, .(gene), func)
order.scores<-order(-cor_matrix$COR)
cor_matrix$rank <- NA
cor_matrix$rank[order.scores] <- 1:nrow(cor_matrix)

# filter dataframe to get data to be highligheted
highlight_df <- cor_matrix %>% 
  filter(gene=="NEFH")
pdf("processed_data/figures_raw/gene_correlation_NEFH_mCSN.pdf")
cor_matrix %>% 
  ggplot(aes(x=rank,y=COR)) + 
  geom_point(alpha=0.5,size=1) +
  geom_point(data=highlight_df, 
             aes(x=rank,y=COR), 
             color='red',
             size=3) + 
  ylab("Gene correlation to mCSN/human enrichment") +
  xlab("Gene rank") 
dev.off()




nefh_exp = gene_exp[gene_exp$gene == "NEFH",]
colors_use <- nefh_exp$cluster_label.dendro.deg.colors
names(colors_use) <- nefh_exp$cluster_label.dendro.deg.colors
pdf("processed_data/figures_raw/correlation_NEFH_mCSN.pdf")
nefh_exp %>%
  ggplot(aes(y=logCPM, x=facs_prop)) + 
  geom_smooth(method='lm',) +
  geom_point(color=colors_use, alpha = 0.5,size=2.5) +
  ylab("logCPM NEFH") +
  xlab("mCSN geneset expression in human cells (%)") 

dev.off()

cor.test(nefh_exp$logCPM,nefh_exp$facs_prop, method = "pearson")

