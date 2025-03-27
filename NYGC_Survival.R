####

### libraries
url <- "https://cran.r-project.org/src/contrib/Archive/NanoStringNorm/NanoStringNorm_1.2.1.1.tar.gz"
pkgFile <- "NanoStringNorm_1.2.1.1.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)


# creb3.bound = read.table('~/Desktop/single_cell/sc_als/processed_data/human_10x_10XALS_integration/Figure_3/CREB3_bound.genes', header = T)
# gene1 = creb3.bound$Gene
# 
# turquoise_genes = read.table('~/Desktop/single_cell/sc_als/processed_data/article_marques/gene_list_modules.txt', header = T)
# turquoise_genes = turquoise_genes[turquoise_genes$module == "turquoise",]
# gene2 = turquoise_genes$gene
# 
# #creb3.reg = read.table('~/Desktop/single_cell/sc_als/processed_data/human_10x_10XALS_integration/Figure_3/CREB3_regulon', header = T)
# #gene3 = creb3.reg$target
# 
# 
# genes = c(gene1,gene2)
# genes = genes[!duplicated(genes)]
#write.table(genes,file = "~/Desktop/single_cell/sc_als/processed_data/human_10x_10XALS_integration/Figure_3/CREB3_regulon", quote = F, row.names = F, col.names = F)
creb3.reg = read.table('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/human_10x_10XALS_integration/Figure_5/CREB3_regulons.txt', header = T)

fus.reg = read.table('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_Raph/wgcna_intron/genesTurquoise_sorted.txt', header = F)
fus.reg = fus.reg$V1

setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/LongGene_code/src/')
source('libraries.R')
source('ComBat-seq/ComBat_seq.R')
source('ComBat-seq/helper_seq.R')
library(jaffelab)
# Load `SummarizedExperiment` object
setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/Human_NYGC/')
nygc <- read.table("ALS.NYGC.all.counts.txt",sep = "\t", 
                          stringsAsFactors=F, header = TRUE, row.names = 1)

cov = read.table('covariates.ALS_All', header = T, stringsAsFactors = FALSE)
clinical = read.table('covariates.ALS_Clinical', header = T, stringsAsFactors = FALSE)
geno = read.table('covariates.ALS_Genotype', header = T, stringsAsFactors = FALSE)
cov = merge(cov,clinical, by=c('Subject_ID'))
cov = cov[!duplicated(cov$Sample_ID),]
cov = merge(cov,geno, by=c("Subject_ID"))
cov = cov[!duplicated(cov$Sample_ID.x),]
cov = subset(cov, Sample_ID.x %in% colnames(nygc))
cortex_idx = colnames(nygc) %in% cov$Sample_ID.x
nygc = nygc[,cortex_idx]
order = match(cov$Sample_ID.x,colnames(nygc),nomatch = 0)
nygc = nygc[,order]

row.names(cov) = cov$Sample_ID.x

### perform DESeq normalization on count and PCA
dds <- DESeqDataSetFromMatrix(nygc, cov, design = ~Sex)
dds$Sex <- relevel(factor(dds$Sex), "F")
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
z <- plotPCA(vsd, intgroup=c("Sex"))
z + geom_label(aes(label = group))


#####
cov$pca = z$data$PC2
cov$batch <- ifelse(cov$pca>=0, 1, 2)

# genotype_sns <- data.frame(c(rep("WT", 12), rep("KI", 12)))
# rownames(genotype_sns) <- colnames(allsns)
# colnames(genotype_sns) = "genotypes"
# genotype_sns$genotypes <- relevel(genotype_sns$genotypes, "WT")

batch = cov$batch
adjusted <- ComBat_seq(nygc, batch=batch, group=cov$AAO)
adjusted = as.data.frame(adjusted)

# Create DGEList and include average transcript length offsets
de.exon <- DGEList(counts  = adjusted,samples  = cov)
keep <- filterByExpr(de.exon, group = de.exon$samples$Group)
de.exon <- de.exon[keep,,keep.lib.sizes=FALSE]
de.exon <- calcNormFactors(de.exon)
design <- model.matrix(~0+AAO, de.exon$samples)
de.exon <- estimateDisp(de.exon,design)

log.exon <- edgeR::cpm(de.exon, log = TRUE, prior.count = 1)
norm.exon <- edgeR::cpm(de.exon, normalized.lib.sizes = TRUE)

################
## load expression

#de.exon$samples$Group = factor(de.exon$samples$Group,
#              levels = c("Control", "ALS"))

mod = model.matrix(~0+AAO + Sex + RIN + PMI + Tissue + DD,
                   data = de.exon$samples)

v_exon <- voom(de.exon, design, plot = TRUE)
fit_exon <- lmFit(v_exon, design)
exprsNorm_exon <- v_exon$E
fit_exon <- eBayes(fit_exon)


## Extract top results
# ortholog = read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/csn_gene_orthologs', header = F)
# top_exon <- topTable(fit_exon, coef = 1, n = nrow(de.exon), sort.by = 'none',adjust = "fdr")

## residualize expression  
gExprs = log2(norm.exon+1)        
gExprs = jaffelab::cleaningY(gExprs, mod, P=1)


gExprs = gExprs[row.names(gExprs) %in% fus.reg,]
mat = apply(gExprs,2, function(x) c( "Stand dev" = sd(x), 
                         "Mean"= mean(x,na.rm=TRUE),
                         "n" = length(x),
                         "Median" = median(x),
                         "CoeffofVariation" = sd(x)/mean(x,na.rm=TRUE),
                         "Minimum" = min(x),
                         "Maximun" = max(x),
                         "Upper Quantile" = quantile(x,1),
                         "LowerQuartile" = quantile(x,0)
)
)


mat = t(mat)


final_mat = merge(cov,mat,by=c("row.names"))
final_mat$Tissue = gsub("Cortex_Frontal", "Cortex_Motor", final_mat$Tissue)
final_mat$Tissue = gsub("Cortex_Motor_Lateral", "Cortex_Motor", final_mat$Tissue)
final_mat = final_mat[final_mat$Sample_ID.x != "CGND_HRA_00245",]
final_mat = final_mat[final_mat$Tissue != "Cortex_Motor_Medial",]

median = median(final_mat$Median)
final_mat$group = ifelse(final_mat$Median >=0.013, "high_CREB3","low_CREB3")



final_mat = final_mat[final_mat$DD < 120,]
colors <- RColorBrewer::brewer.pal(9, "PuRd")[5]

pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_survival_transcriptomic.pdf")
final_mat %>%
  ggplot(aes(x=AAO, y=Mean, color=Tissue)) +
  geom_point(alpha=0.5, size=3) + 
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  scale_color_manual(values=colors)+
  theme(legend.position="top")
dev.off()
aggregate(final_mat[, c(16)], list(final_mat$group), mean)

corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="pearson"))
}
cor_table = ddply(final_mat, .(Tissue), summarise,z=corfun(Median,AAO)$statistic,
                  pval=corfun(Median,AAO)$p.value,
                  tau.est=corfun(Median,AAO)$estimate,
                  alt=corfun(Median,AAO)$alternative
) 


library(outliers)
library(ggsurvfit)
library(survival)

final_mat$time = final_mat$AAO
final_mat$status = 1
final_mat = final_mat[!duplicated(final_mat$Sample_ID.y),]

km_trt_fit <- survfit(Surv(time, status) ~ group, data=final_mat)

summary(survfit(Surv(time, status) ~ group, data = final_mat), times = 25)

surv_diff <- survdiff(Surv(time, status) ~ group, data = final_mat)
surv_diff

pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_survival_transcriptomic_kaplan.pdf")
survfit2(Surv(time, status) ~ group, data = final_mat) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) 
dev.off()
aggregate(final_mat[, c(16)], list(final_mat$group), mean)




dodge <- position_dodge(width = 0.4)

pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_aao_transcriptomic.pdf", width = 8, height = 12)
final_mat %>%
  ggplot(aes(x=group, y=AAO, color=group)) +
  geom_violin(aes(color = group),position = dodge) + 
  geom_boxplot(width=.1, outlier.colour=NA, position = dodge) +
  theme(legend.position="top") + ylim(0,max(final_mat$AAO))
dev.off()
aggregate(final_mat[, c(15)], list(final_mat$group), mean)







