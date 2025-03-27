####

### libraries
url <- "https://cran.r-project.org/src/contrib/Archive/NanoStringNorm/NanoStringNorm_1.2.1.1.tar.gz"
pkgFile <- "NanoStringNorm_1.2.1.1.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)


creb3.reg = read.table('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/human_10x_10XALS_integration/Figure_3/CREB3_bound.genes', header = T)
#write.table(genes,file = "~/Desktop/single_cell/sc_als/processed_data/human_10x_10XALS_integration/Figure_3/CREB3_regulon", quote = F, row.names = F, col.names = F)

setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/LongGene_code/src/')
source('libraries.R')
source('ComBat-seq/ComBat_seq.R')
source('ComBat-seq/helper_seq.R')
library(jaffelab)
library(org.Hs.eg.db)
library(annotate)
# Load `SummarizedExperiment` object
setwd('~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/')
nygc <- read.table("human_10x_10XALS_integration/Figure_6/sALS_RNAseq.blood.txt",sep = "\t", 
                   stringsAsFactors=F, header = TRUE)

map = as.data.frame(getSYMBOL(as.character(nygc$EntrezGeneID), data='org.Hs.eg'))
nygc = cbind(map,nygc)
colnames(nygc)[1] <- 'genes'
nygc = nygc[!duplicated(nygc$genes),]
nygc = na.omit(nygc)
row.names(nygc) = nygc[,1]
nygc = nygc[,-c(1,2)]

cov = read.csv('human_10x_10XALS_integration/Figure_6/covariates.csv')
cov = cov[cov$Disease_status == "sALS",]

cortex_idx = colnames(nygc) %in% cov$Collection_ID
nygc = nygc[,cortex_idx]
order = match(cov$Collection_ID,colnames(nygc),nomatch = 0)
nygc = nygc[,order]
row.names(cov) = cov$Collection_ID

### perform DESeq normalization on count and PCA

dds <- DESeqDataSetFromMatrix(round(nygc), cov, design = ~Sex)
dds$Sex <- relevel(factor(dds$Sex), "F")
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
z <- plotPCA(vsd, intgroup=c("Site_of_onset"))
z + geom_label(aes(label = group))


# #####
# cov$pca = z$data$PC2
# cov$batch <- ifelse(cov$pca>=0, 1, 2)
# 
# # genotype_sns <- data.frame(c(rep("WT", 12), rep("KI", 12)))
# # rownames(genotype_sns) <- colnames(allsns)
# # colnames(genotype_sns) = "genotypes"
# # genotype_sns$genotypes <- relevel(genotype_sns$genotypes, "WT")
# 
# batch = cov$batch
# adjusted <- ComBat_seq(nygc, batch=batch, group=cov$AAO)
# adjusted = as.data.frame(adjusted)
cov$Age_of_onset_code = ifelse(cov$Age_of_onset >= 61, "old", "young")

# Create DGEList and include average transcript length offsets
de.exon <- DGEList(counts  = nygc,samples  = cov)
keep <- filterByExpr(de.exon, group = de.exon$samples$group)
de.exon <- de.exon[keep,,keep.lib.sizes=FALSE]
de.exon <- calcNormFactors(de.exon)
design <- model.matrix(~0+Disease_duration_months, de.exon$samples)
de.exon <- estimateDisp(de.exon,design)

log.exon <- edgeR::cpm(de.exon, log = TRUE, prior.count = 1)
norm.exon <- edgeR::cpm(de.exon, normalized.lib.sizes = TRUE)

################
## load expression

#de.exon$samples$Group = factor(de.exon$samples$Group,
#              levels = c("Control", "ALS"))


mod = model.matrix(~0+Disease_duration_months + Sex,
                   data = de.exon$samples)

v_exon <- voom(de.exon, design, plot = TRUE)
fit_exon <- lmFit(v_exon, design)
exprsNorm_exon <- v_exon$E
fit_exon <- eBayes(fit_exon)


## Extract top results
# ortholog = read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/csn_gene_orthologs', header = F)
# top_exon <- topTable(fit_exon, coef = 1, n = nrow(de.exon), sort.by = 'none',adjust = "fdr")

## residualize expression  
library(jaffelab)
gExprs = log2(norm.exon+1)        
gExprs = cleaningY(gExprs, mod, P=1)


gExprs = gExprs[row.names(gExprs) %in% creb3.reg$Gene,]
mat = apply(gExprs,2, function(x) c( "Stand dev" = sd(x), 
                                     "Mean"= mean(x,na.rm=TRUE),
                                     "n" = length(x),
                                     "Median" = median(x),
                                     "CoeffofVariation" = sd(x)/mean(x,na.rm=TRUE),
                                     "Minimum" = min(x),
                                     "Maximun" = max(x),
                                     "Upper Quantile" = quantile(x,.10),
                                     "LowerQuartile" = quantile(x,.90)
)
)


mat = t(mat)

colors <- RColorBrewer::brewer.pal(9, "PuRd")[5]
final_mat = merge(cov,mat,by=c("row.names"))
median = median(final_mat$Median)
final_mat$group = ifelse(final_mat$Median >=median, "high_CREB3","low_CREB3")


pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_survival_transcriptomic_blood.pdf")
final_mat %>%
  ggplot(aes(x=Disease_duration_months, y=Median, color=Disease_status)) +
  geom_point(alpha=0.5, size=3) + 
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  scale_color_manual(values=colors)+
  theme(legend.position="top")
dev.off()

dodge <- position_dodge(width = 0.4)

pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_aao_transcriptomic_blood.pdf", width = 8, height = 12)
final_mat %>%
  ggplot(aes(x=group, y=Age_of_onset, color=group)) +
  geom_violin(aes(color = group),position = dodge) + 
  geom_boxplot(width=.1, outlier.colour=NA, position = dodge) +
  theme(legend.position="top") + ylim(0,max(final_mat$Age_of_onset))
dev.off()
aggregate(final_mat[, c(6)], list(final_mat$group), mean)

library(tidyverse)
library(ggpubr)
library(rstatix)
final_mat %>%
  group_by(group) %>%
  get_summary_stats(Age_of_onset, type = "mean_sd")

final_mat %>%
  group_by(group) %>%
  shapiro_test(Age_of_onset)
final_mat %>% levene_test(Age_of_onset ~ group)

stat.test2 <- final_mat %>%
  t_test(Age_of_onset ~ group, var.equal = F) %>%
  add_significance()
stat.test2



cor(final_mat$Disease_duration_months,final_mat$Median)
corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="pearson"))
}

test = corfun(x = final_mat$Disease_duration_months,final_mat$Median)

library(outliers)
library(ggsurvfit)

final_mat$time = final_mat$Disease_duration_months
final_mat$status = final_mat$Deceased

km_trt_fit <- survfit(Surv(time, status) ~ group, data=final_mat)
summary(survfit(Surv(time, status) ~ group, data = final_mat), times = 25)

surv_diff <- survdiff(Surv(time, status) ~ group, data = final_mat)
surv_diff

pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_survival_transcriptomic_blood_kaplan.pdf")
survfit2(Surv(time, status) ~ group, data = final_mat) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) 
dev.off()
aggregate(final_mat[, c(9)], list(final_mat$group), mean)




final_mat_blood = final_mat






