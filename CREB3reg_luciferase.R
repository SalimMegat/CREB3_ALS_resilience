library(WGCNA)
library(edgeR)
library(tidyverse)
library(SummarizedExperiment)
library(limma)
library(devtools)
library(tidyr)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(EnhancedVolcano)
library(ggrepel)
library(pheatmap)
library(ggplot2)
library(sva)
library(ggpubr)
setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/LongGene_code/src/')
source('libraries.R')
source('ComBat-seq/ComBat_seq.R')
source('ComBat-seq/helper_seq.R')
library(enrichR)
library(gplots)
setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/Exon_Intron/sod_christine/')
##### Load Laura's data
load("~/Desktop/iMac_U1118/Desktop/Analysis_fus/Exon_Intron/sod_christine/Counts/Intron.RData")
magda_intron_count <- as.data.frame(counts$counts[,c('CNRX57', 'CNRX77','CNRX83', 'CNRX103', 'CNRX79', 'CNRX81', 'CNRX105', 'CNRX67', 'CNRX69', 'CNRX87', 'CNRX133', 'CNRX71', 'CNRX91', 'CNRX95', 'CNRX131', 'CNRX63', 'CNRX117', 'CNRX137', 'CNRX143', 'CNRX65', 'CNRX119', 'CNRX139','CNRX153', 'CNRX97', 'CNRX101', 'CNRX107', 'CNRX125', 'CNRX109', 'CNRX157', 'CNRX58', 'CNRX78', 'CNRX84','CNRX104','CNRX62', 'CNRX80', 'CNRX82', 'CNRX68', 'CNRX70', 'CNRX92', 'CNRX96', 'CNRX132', 'CNRX64', 'CNRX118', 'CNRX138', 'CNRX144', 'CNRX66', 'CNRX120', 'CNRX140','CNRX154', 'CNRX98', 'CNRX102', 'CNRX108', 'CNRX126','CNRX110', 'CNRX128', 'CNRX152', 'CNRX158')])
colnames(magda_intron_count) = c("CSN_WT1_30d","CSN_WT2_30d","CSN_WT3_30d","CSN_WT4_30d","CSN_SOD1_30d", "CSN_SOD2_30d", "CSN_SOD3_30d","CSN_WT1_60d","CSN_WT2_60d","CSN_WT3_60d","CSN_WT4_60d","CSN_SOD1_60d", "CSN_SOD2_60d", "CSN_SOD3_60d", "CSN_SOD4_60d","CSN_WT1_90d","CSN_WT2_90d","CSN_WT3_90d","CSN_WT4_90d","CSN_SOD1_90d", "CSN_SOD2_90d", "CSN_SOD3_90d", "CSN_SOD4_90d", "CSN_WT1_105d","CSN_WT2_105d","CSN_WT3_105d","CSN_WT4_105d","CSN_SOD1_105d", "CSN_SOD2_105d","CPN_WT1_30d","CPN_WT2_30d","CPN_WT3_30d","CPN_WT4_30d","CPN_SOD1_30d", "CPN_SOD2_30d", "CPN_SOD3_30d","CPN_WT1_60d","CPN_WT2_60d","CPN_SOD1_60d", "CPN_SOD2_60d", "CPN_SOD3_60d","CPN_WT1_90d","CPN_WT2_90d","CPN_WT3_90d","CPN_WT4_90d","CPN_SOD1_90d", "CPN_SOD2_90d", "CPN_SOD3_90d", "CPN_SOD4_90d","CPN_WT1_105d","CPN_WT2_105d","CPN_WT3_105d","CPN_WT4_105d", "CPN_SOD1_105d", "CPN_SOD2_105d", "CPN_SOD3_105d", "CPN_SOD4_105d")
row.names(magda_intron_count) = gsub('\\..*', '', row.names(magda_intron_count))
##### Load Laura's data
load("~/Desktop/iMac_U1118/Desktop/Analysis_fus/Exon_Intron/sod_christine/Counts/Exon.RData")
magda_exon_count <- as.data.frame(counts$counts[,c('CNRX57', 'CNRX77','CNRX83', 'CNRX103', 'CNRX79', 'CNRX81', 'CNRX105', 'CNRX67', 'CNRX69', 'CNRX87', 'CNRX133', 'CNRX71', 'CNRX91', 'CNRX95', 'CNRX131', 'CNRX63', 'CNRX117', 'CNRX137', 'CNRX143', 'CNRX65', 'CNRX119', 'CNRX139','CNRX153', 'CNRX97', 'CNRX101', 'CNRX107', 'CNRX125', 'CNRX109', 'CNRX157', 'CNRX58', 'CNRX78', 'CNRX84','CNRX104','CNRX62', 'CNRX80', 'CNRX82', 'CNRX68', 'CNRX70', 'CNRX92', 'CNRX96', 'CNRX132', 'CNRX64', 'CNRX118', 'CNRX138', 'CNRX144', 'CNRX66', 'CNRX120', 'CNRX140','CNRX154', 'CNRX98', 'CNRX102', 'CNRX108', 'CNRX126','CNRX110', 'CNRX128', 'CNRX152', 'CNRX158')])
colnames(magda_exon_count) = c("CSN_WT1_30d","CSN_WT2_30d","CSN_WT3_30d","CSN_WT4_30d","CSN_SOD1_30d", "CSN_SOD2_30d", "CSN_SOD3_30d","CSN_WT1_60d","CSN_WT2_60d","CSN_WT3_60d","CSN_WT4_60d","CSN_SOD1_60d", "CSN_SOD2_60d", "CSN_SOD3_60d", "CSN_SOD4_60d","CSN_WT1_90d","CSN_WT2_90d","CSN_WT3_90d","CSN_WT4_90d","CSN_SOD1_90d", "CSN_SOD2_90d", "CSN_SOD3_90d", "CSN_SOD4_90d", "CSN_WT1_105d","CSN_WT2_105d","CSN_WT3_105d","CSN_WT4_105d","CSN_SOD1_105d", "CSN_SOD2_105d","CPN_WT1_30d","CPN_WT2_30d","CPN_WT3_30d","CPN_WT4_30d","CPN_SOD1_30d", "CPN_SOD2_30d", "CPN_SOD3_30d","CPN_WT1_60d","CPN_WT2_60d","CPN_SOD1_60d", "CPN_SOD2_60d", "CPN_SOD3_60d","CPN_WT1_90d","CPN_WT2_90d","CPN_WT3_90d","CPN_WT4_90d","CPN_SOD1_90d", "CPN_SOD2_90d", "CPN_SOD3_90d", "CPN_SOD4_90d","CPN_WT1_105d","CPN_WT2_105d","CPN_WT3_105d","CPN_WT4_105d", "CPN_SOD1_105d", "CPN_SOD2_105d", "CPN_SOD3_105d", "CPN_SOD4_105d")
row.names(magda_exon_count) = gsub('\\..*', '', row.names(magda_exon_count))
##### Load Laura's data
load("~/Desktop/iMac_U1118/Desktop/Analysis_fus/Exon_Intron/sod_christine/Counts/Genebody.RData")
magda_gene_count <- as.data.frame(counts$counts[,c('CNRX57', 'CNRX77','CNRX83', 'CNRX103', 'CNRX79', 'CNRX81', 'CNRX105', 'CNRX67', 'CNRX69', 'CNRX87', 'CNRX133', 'CNRX71', 'CNRX91', 'CNRX95', 'CNRX131', 'CNRX63', 'CNRX117', 'CNRX137', 'CNRX143', 'CNRX65', 'CNRX119', 'CNRX139','CNRX153', 'CNRX97', 'CNRX101', 'CNRX107', 'CNRX125', 'CNRX109', 'CNRX157', 'CNRX58', 'CNRX78', 'CNRX84','CNRX104','CNRX62', 'CNRX80', 'CNRX82', 'CNRX68', 'CNRX70', 'CNRX92', 'CNRX96', 'CNRX132', 'CNRX64', 'CNRX118', 'CNRX138', 'CNRX144', 'CNRX66', 'CNRX120', 'CNRX140','CNRX154', 'CNRX98', 'CNRX102', 'CNRX108', 'CNRX126','CNRX110', 'CNRX128', 'CNRX152', 'CNRX158')])
colnames(magda_gene_count) = c("CSN_WT1_30d","CSN_WT2_30d","CSN_WT3_30d","CSN_WT4_30d","CSN_SOD1_30d", "CSN_SOD2_30d", "CSN_SOD3_30d","CSN_WT1_60d","CSN_WT2_60d","CSN_WT3_60d","CSN_WT4_60d","CSN_SOD1_60d", "CSN_SOD2_60d", "CSN_SOD3_60d", "CSN_SOD4_60d","CSN_WT1_90d","CSN_WT2_90d","CSN_WT3_90d","CSN_WT4_90d","CSN_SOD1_90d", "CSN_SOD2_90d", "CSN_SOD3_90d", "CSN_SOD4_90d", "CSN_WT1_105d","CSN_WT2_105d","CSN_WT3_105d","CSN_WT4_105d","CSN_SOD1_105d", "CSN_SOD2_105d","CPN_WT1_30d","CPN_WT2_30d","CPN_WT3_30d","CPN_WT4_30d","CPN_SOD1_30d", "CPN_SOD2_30d", "CPN_SOD3_30d","CPN_WT1_60d","CPN_WT2_60d","CPN_SOD1_60d", "CPN_SOD2_60d", "CPN_SOD3_60d","CPN_WT1_90d","CPN_WT2_90d","CPN_WT3_90d","CPN_WT4_90d","CPN_SOD1_90d", "CPN_SOD2_90d", "CPN_SOD3_90d", "CPN_SOD4_90d","CPN_WT1_105d","CPN_WT2_105d","CPN_WT3_105d","CPN_WT4_105d", "CPN_SOD1_105d", "CPN_SOD2_105d", "CPN_SOD3_105d", "CPN_SOD4_105d")
row.names(magda_gene_count) = gsub('\\..*', '', row.names(magda_gene_count))


genotype <- data.frame(c(rep("WT", 4), rep("SOD", 3),rep("WT", 4),rep("SOD", 4),rep("WT", 4),rep("SOD", 4),rep("WT", 4),rep("SOD", 2),rep("WT", 4),rep("SOD", 3),rep("WT", 2),rep("SOD", 3),rep("WT", 4),rep("SOD", 4),rep("WT", 4),rep("SOD", 4)))
genotype$tissue <- c(rep("CSN", 29), rep("CPN", 28))
genotype$status <- c(rep("pre", 7), rep("pre", 8),rep("post", 8),rep("post", 6),rep("pre", 7),rep("pre", 5),rep("post", 8),rep("post", 8))
genotype$age <- c(rep("30d", 7), rep("60d", 8),rep("90d", 8),rep("105d", 6),rep("30d", 7),rep("60d", 5),rep("90d", 8),rep("105d", 8))
genotype$genotype_tissue_geno = paste0(genotype$tissue, '_', genotype$c.rep..WT...4...rep..SOD...3...rep..WT...4...rep..SOD...4...rep..WT...)
genotype$group = paste0(genotype$genotype_tissue_geno,'_',genotype$age)

row.names(genotype) <- colnames(magda_gene_count)
colnames(genotype)[1]= c("genotypes")

#genotype$batch <- c(1,1,2,2,2,1,3,1,2,3,2,2,2,1,1,3,2,3,2,3,3,3,3,3,2,1,1,1,1,1,1,2,3,3,3,3,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,1,1,1,1,1)
genotype$age_class = c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4)

genotype$batch <- runif(nrow(genotype), min = 0, max = 4)
genotype$batch <- round(genotype$batch)
genotype$batch <- paste0("batch", genotype$batch)

##### Load Salim's data
## get the gene ID info 
genebody_info <- read.table("~/Desktop/iMac_U1118/Desktop/Analysis_fus/Exon_Intron/tdp43KI_white/Gene_information.txt", header = TRUE, row.names = 1)
row.names(genebody_info) = gsub('\\..*', '', row.names(genebody_info))
genebody_info$id = row.names(genebody_info)
genebody_info = genebody_info[,-c(1,3)]
#write.table(genebody_info, '~/Desktop/Analysis_fus/Exon_Intron/wgcna_intron/annotation_mouse.txt', row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
## merge intron with genebody info to get the gene name
intron_count_batch<- merge(genebody_info,magda_gene_count, by=c("row.names"))
intron_count_batch = intron_count_batch[!duplicated(intron_count_batch$Symbol),]
row.names(intron_count_batch) <- intron_count_batch[,2]
intron_count_batch <- intron_count_batch[,-c(1:3)]


gene = genebody_info[!duplicated(genebody_info$Symbol),]



dds <- DESeqDataSetFromMatrix(round(intron_count_batch), 
                              colData = genotype, 
                              design = ~group)
relevel(dds$group, ref = "CSN_WT_30d")
# Transform counts for data visualization
vsd <- DESeq2::vst(dds, blind=FALSE,nsub = )
z <- DESeq2::plotPCA(vsd, intgroup=c("tissue"))
z + geom_label(aes(label = name))
plot(z)




vsd.csn <- vsd[ , vsd$tissue %in% c("CSN") ]
vsd.cpn <- vsd[ , vsd$tissue %in% c("CPN") ]


a<-plotPCA(vsd.csn, "age")
b<-plotPCA(vsd.csn, "genotypes")
e<-plotPCA(vsd.csn, "batch")

c<-plotPCA(vsd.cpn, "age")
d<-plotPCA(vsd.cpn, "genotypes")
f<- plotPCA(vsd.cpn, "batch")

p<- ggarrange(a, b, e, c,d,f, 
              ncol = 3, nrow = 2)



ggsave(filename = '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/figures_raw/PCA_CSN_allgroup.pdf', plot = p)




##### Pipeline edgeR differential expression ######
## merge intron with genebody info to get the gene name
de.exon <- DGEList(counts  = adjusted,samples  = genotype )
de.exon <- calcNormFactors(de.exon)
keep <- filterByExpr(de.exon, group = de.exon$samples$genotypes )
de.exon <- de.exon[keep,,keep.lib.sizes=FALSE]
design <- model.matrix(~0+group,de.exon$samples)
de.exon <- estimateDisp(de.exon,design)
norm.exon <- edgeR::cpm(de.exon, normalized.lib.sizes = TRUE)
log.exon <- edgeR::cpm(de.exon, log = TRUE, prior.count = 1)
log.exon = as.data.frame(log.exon)
#log.exon$gene_name = row.names(log.exon)
dim(de.exon)
## plot PCA plot 
### Raw PCA with ggplot 
#We model the count data using a quasi-likelihood (QL) negative binomial (NB)
#generalized log-linear model, which accounts for gene-specific variability from
#both biological and technical sources. Before fitting the model, we estimate
#the NB dispersion (overall biological variability across all genes), and the QL
#dispersion (gene-specific) using the `estimateDisp()` function.
#```{r edgeR-estimate-disp}
## Estimate dispersion and fit model
qlfit <- glmQLFit(de.exon, design = design)
## Plot dispersions
plotBCV(de.exon)


# Define contrasts 

#Before testing for differences in gene expression, we define the contrasts
#we wish to test for. Here we represent the constrasts as a numeric matrix:

contrasts <- as.data.frame(makeContrasts(CSN_pre=groupSOD.pre.CSN-groupWT.pre.CSN, CSN_post=groupSOD.post.CSN-groupWT.post.CSN,CPN_pre=groupSOD.pre.CPN-groupWT.pre.CPN,CPN_post=groupSOD.post.CPN-groupWT.post.CPN, CSN_enriched=groupWT.pre.CSN-groupWT.pre.CPN,CPN_enriched=groupWT.pre.CPN-groupWT.pre.CSN, levels = design))
#contrasts <- as.data.frame(makeContrasts(CSN_enriched=groupCSN_WT-groupCPN_WT,CPN_enriched=groupCPN_WT-groupCSN_WT, levels = mm))

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


#### Split data and look grouped of genes variantion across time 
mat_melted = reshape2::melt(as.matrix(log.exon))
colnames(mat_melted) = c("gene", "sampleID", "counts")

mat_melted = merge(mat_melted,genotype,by=c("sampleID"))
mat_melted$genotypes <- factor(mat_melted$genotype,levels = c("WT","SOD"))

gene_interest = "Creb3"

creb3_mat = mat_melted[mat_melted$gene == gene_interest,]

res.aov3 <- aov(counts ~ status*genotypes, data = creb3_mat)
summary(res.aov3)
TukeyHSD(res.aov3, which = "status:genotypes")
# Summary of the analysis
plot(res.aov3, 2)

library(viridis)
library(hrbrthemes)

min = min(creb3_mat$counts) - 0.4
max = max(creb3_mat$counts) + 0.4

pdf("~/Desktop/single_cell/sc_als/processed_data/figures_raw/Creb3_mouse_CSN.pdf")
creb3_mat %>%
  mutate(status = factor(status, levels = c("pre","post"))) %>%
  ggplot(aes(fill=genotypes, y=counts, x=status)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()


###############################################################################
###############################################################################
############ Load mouse CSN dataset and target TF gene expression ##############
###############################################################################
###############################################################################
tf.targets = read.table('~/Desktop/single_cell/sc_als/processed_data/WGCNA/TFs_turquoise_target', header = T)
tf.list = unique(tf.targets$TF)

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

deg_file.mouse = rbind(csn_edgeR_deg.post,csn_edgeR_deg.pre)
deg_file.mouse = merge(deg_file.mouse,human_homologous,by=c("gene_name"))

tf_mouse_array <- lapply(tf.list, function(tf) {
  target = tf.targets[tf.targets$TF == tf,]
  deg_file.mouse = deg_file.mouse[deg_file.mouse$human_genes %in% target$Target,]
  deg_file.mouse$fcsign = sign(deg_file.mouse$logFC)
  deg_file.mouse$metric= deg_file.mouse$mlog10PValue/deg_file.mouse$fcsign
  deg_file.mouse$TF  = tf
  tmp <- rbind(deg_file.mouse)
})

final_mouse_tf <- do.call(rbind,tf_mouse_array)

final_mouse_tf$time <- factor(final_mouse_tf$time,levels = c("pre", "post"))

### Select CREB3 TF #####
creb3_tf = final_mouse_tf[final_mouse_tf$TF == "CREB3",]
creb3_tf = creb3_tf[creb3_tf$gene_name != "Sod1",]
write.table(creb3_tf,'~/Desktop/single_cell/sc_als/processed_data/article_marques/creb3_tfMouse.txt', sep = "\t", quote = F, col.names = T, row.names = F)


# Build the linear model
library(tidyverse)
library(ggpubr)
library(rstatix)
model  <- lm(metric ~ time, data = creb3_tf)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
### If data not normal #####
min = min(creb3_tf$metric) - 0.1
max = max(creb3_tf$metric) + 0.1

pdf("~/Desktop/single_cell/sc_als/processed_data/figures_raw/Creb3_mouse_CSN_targets.pdf")
creb3_tf %>%
  mutate(time = factor(time, levels = c("pre","post"))) %>%
  ggplot(aes(fill=time, y=metric, x=time)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()

res <- wilcox.test(metric ~ time, data = creb3_tf,
                   exact = FALSE)
res


#######################################################################################################################
#######################################################################################################################
##### Check expression of the CREB3 regulons in mouse motorneuros SOD1 model single cell (Commuinication Biology) #####
#######################################################################################################################
#######################################################################################################################
#creb3_bound_gene = tf.targets[tf.targets$TF == "CREB3",]

library(viridis)
library(hrbrthemes)



creb3.bound = read.table('~/Desktop/single_cell/sc_als/processed_data/human_10x_10XALS_integration/Figure_3/CREB3_bound.genes', header = T)

tmp1 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/motor_neurons_ALS.csv', header = T)
tmp1$cell = "motorneurons"
tmp2 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/neurons_ALS.csv', header = T)
tmp2$cell = "other_neurons"
tmp3 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/astrocytes_ALS.csv', header = T)
tmp3$cell = "astrocytes"
tmp4 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/mircroglia_ALS.csv', header = T)
tmp4$cell = "microglia"
tmp5 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/endothelial_ALS.csv', header = T)
tmp5$cell = "endothelial"


tmp= rbind(tmp1,tmp2,tmp3,tmp4,tmp5)
tmp = merge(tmp,human_homologous,by=c("gene_name"))
tmp = tmp[tmp$human_genes %in% creb3.bound$Gene,]
tmp$fcsign = sign(tmp$logFC)
tmp$metric= -log10(tmp$P.Value)/tmp$fcsign

mean = aggregate(tmp[, 13], list(tmp$cell), mean)


min = min(tmp$metric) - 0.1
max = max(tmp$metric) + 0.1

pdf("~/Desktop/single_cell/sc_als/processed_data/figures_raw/Creb3_mouse_MN_targets.pdf")
tmp %>%
  mutate(cell = factor(cell, levels = c("astrocytes","microglia","endothelial", "motorneurons", "other_neurons"))) %>%
  ggplot(aes(fill=cell, y=metric, x=cell)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()

##### Non parametric to test average z score across cell in the SOD1 spinal cord ####
# Pairwise comparisons
pwc2 <- tmp %>% 
  wilcox_test(metric ~ cell, p.adjust.method = "fdr")
pwc2


#######################################################################################################################
#######################################################################################################################
##### Check expression of the turquoise module  mouse motorneuros SOD1 model single cell (Commuinication Biology) #####
#######################################################################################################################
#######################################################################################################################

library(viridis)
library(hrbrthemes)

module = read.table('~/Desktop/single_cell/sc_als/processed_data/article_marques/gene_list_modules.txt', header = T)
module = module[module$module == "turquoise",]

tmp1 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/motor_neurons_ALS.csv', header = T)
tmp1$cell = "motorneurons"
tmp2 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/neurons_ALS.csv', header = T)
tmp2$cell = "other_neurons"
tmp3 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/astrocytes_ALS.csv', header = T)
tmp3$cell = "astrocytes"
tmp4 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/mircroglia_ALS.csv', header = T)
tmp4$cell = "microglia"
tmp5 = read.csv('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/endothelial_ALS.csv', header = T)
tmp5$cell = "endothelial"


tmp= rbind(tmp1,tmp2,tmp3,tmp4,tmp5)
tmp = merge(tmp,human_homologous,by=c("gene_name"))
tmp = tmp[tmp$human_genes %in% module$gene,]
tmp$fcsign = sign(tmp$logFC)
tmp$metric= -log10(tmp$P.Value)/tmp$fcsign

mean = aggregate(tmp[, 13], list(tmp$cell), mean)


min = min(tmp$metric) - 0.1
max = max(tmp$metric) + 0.1

pdf("~/Desktop/single_cell/sc_als/processed_data/figures_raw/Turquoise_mouse_MN_targets.pdf")
tmp %>%
  mutate(cell = factor(cell, levels = c("astrocytes","microglia","endothelial", "other_neurons", "motorneurons"))) %>%
  ggplot(aes(fill=cell, y=metric, x=cell)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()

# Pairwise comparisons
pwc2 <- tmp %>% 
  wilcox_test(metric ~ cell, p.adjust.method = "fdr")
pwc2

#######################################################################################################################
#######################################################################################################################
##### Check expression of the turquoise module in mouse corticospinal neuron (Moya et al 2022)#########################
#######################################################################################################################
#######################################################################################################################

library(viridis)
library(hrbrthemes)


module = read.table('~/Desktop/single_cell/sc_als/processed_data/article_marques/gene_list_modules.txt', header = T)
module = module[module$module == "turquoise",]

moya = read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Moya/moya_all_cell.txt', header = T,sep = ' ')

human_homologous =  read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/orthologous.cleaned', header = T)
moya = merge(moya,human_homologous,by=c("gene_name"))

moya = moya[moya$human_genes %in% module$gene,]
moya$fcsign = sign(moya$log2fc)
moya$metric= -log10(moya$pvalue)/moya$fcsign
moya = moya[complete.cases(moya$metric),]


min = min(moya$metric) - 0.1
max = max(moya$metric) + 0.1

pdf("~/Desktop/single_cell/sc_als/processed_data/figures_raw/Turquoise_mouse_Moya_targets.pdf")
moya %>%
  mutate(cell = factor(cell, levels = c("M1","Cogalt2","Gprin3"))) %>%
  ggplot(aes(fill=cell, y=metric, x=cell)) + 
  geom_violin(alpha=0.5) +
  geom_boxplot(width = 0.15) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()

aggregate(moya$metric, list(moya$cell), FUN=mean) 



#######################################################################################################################
#######################################################################################################################
##### Analysis of Luciferase Assay with CREB3 R119G mutation  ########################################################
#######################################################################################################################
#######################################################################################################################
library(tidyverse)
library(data.table)
library(ggpubr)
library(rstatix)
library(viridis)
library(multcomp)
luciferase = read.csv('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/article_marques/rapporteur_luciferase/table_luciferase.csv')
luciferase$manip = paste0("m", luciferase$manip)
luciferase = luciferase[order(luciferase$manip),]
model  <- lm(luciferase ~ plasmide, data = luciferase)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
luciferase %>%
  group_by(plasmide) %>%
  shapiro_test(luciferase)
library(car)


# Levene's test with one independent variable
res.aov_creb3 <- aov(luciferase ~ plasmide, data = luciferase)
summary(res.aov_creb3)
TukeyHSD(res.aov_creb3, which = "plasmide")


library(lme4)
library(lmerTest)
summary(aov(luciferase~plasmide + Error(interaction(plasmide,manip)),data=luciferase))
model <- lmer(luciferase~plasmide + (1|manip), data = luciferase)
summary(model)
summary(glht(model, linfct = mcp(plasmide = "Tukey")), test = adjusted("bonferroni"))


ggplot(luciferase, aes(x=factor(manip), y=luciferase, fill=plasmide)) +
  geom_boxplot()


# Summary of the analysis
aggregate(luciferase$luciferase, list(luciferase$plasmide), FUN=mean) 

min = min(luciferase$luciferase) - 0.1
max = max(luciferase$luciferase) + 0.1

pdf("~/Desktop/single_cell/sc_als/processed_data/figures_raw/Creb3_luciferase.pdf")
luciferase %>%
  mutate(plasmide = factor(plasmide, levels = c("pCMV6","pCMV6_CREB3","pCMV6_CREB3_R119G"))) %>%
  ggplot(aes(fill=plasmide, y=luciferase, x=plasmide)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Luciferase activity") +
  ylim(min,max)
dev.off()
