library(WGCNA)
library(edgeR)
library(tidyverse)
library(SummarizedExperiment)
library(limma)
library(GenomicRanges)
library(devtools)
library(tidyr)
library(DESeq2)
library(sva)
library(RColorBrewer)
library(scales)
library(data.table)
library(magrittr)
library(dplyr)
library(EnhancedVolcano)
library(ggrepel)
library(pheatmap)
library(plyr)
library(ggplot2)
setwd('~/Desktop/iMac_U1118/Desktop/Analysis_fus/LongGene_code/src/')
source('libraries.R')
source('ComBat-seq/ComBat_seq.R')
source('ComBat-seq/helper_seq.R')
library(enrichR)
library(gplots)
library(igraph)
##### Load Laura's data

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
#  lnames = load(file = "~/Desktop/Analysis_fus/Exon_Intron/sod_christine/L5ET_mouseInputWGCNA.RData");
#  human_homologous =  read.table('~/Desktop/Analysis_fus/Exon_Intron/sod_christine/orthologous.cleaned', header = T)
#  human_homologous = human_homologous[!duplicated(human_homologous$human_genes),]
#  row.names(human_homologous) = human_homologous$gene_name
# 
#  datExpr = t(datExpr)
#  datExpr = merge(datExpr,human_homologous,by=c("row.names"))
#  row.names(datExpr) = datExpr$human_genes
#  datExpr = datExpr[-c(1,31,32)]
# 
#  fus_nls = t(datExpr)
#  fus_nls_traits = datTraits
# 
# # #The variable lnames contains the names of loaded variables.
# lnames = load(file = "~/Desktop/Analysis_fus/Exon_Intron/sod_christine/L2IT_humanInputWGCNA.RData");
# human_l5et = datExpr
# human_l5et_traits = datTraits
# 
# lnames = load(file = "~/Desktop/Analysis_fus/Exon_Intron/sod_christine/L5ET_humanInputWGCNA.RData");
# human_l2it = datExpr
# human_l2it_traits = datTraits
# 
# 
# 
# 
# 
# # Display the current working directory
# # If necessary, change the path below to the directory where the data files are stored.
# # "." means current directory. On Windows use a forward slash / instead of the usual \.
# workingDir = "~/Desktop/Analysis_fus/Exon_Intron/sod_christine/";
# setwd(workingDir);
# # Load the package
# library(WGCNA);
# # The following setting is important, do not omit.
# options(stringsAsFactors = FALSE);
# #Read in the female liver data set
# idx = colnames(fus_nls)
# idx2 = colnames(human_l5et)
# idx3 = colnames(human_l2it)
# 
# final = Reduce(intersect, list(idx,idx2,idx3))
# 
# fus_nls = fus_nls[,colnames(fus_nls) %in% final]
# human_l5et = human_l5et[,colnames(human_l5et) %in% final]
# human_l2it = human_l2it[,colnames(human_l2it) %in% final]
# 
# 
# human_l2it = as.data.frame(human_l2it)
# human_l5et = as.data.frame(human_l5et)
# fus_nls = as.data.frame(fus_nls)
# 
# dim(fus_nls)
# dim(human_l5et)
# dim(human_l2it)
# 
# genes=colnames(fus_nls)
# # #=====================================================================================
# # #
# # #  Code chunk 2
# # #
# # #=====================================================================================
# #
# #
# # We work with two sets:
# nSets = 3;
# # For easier labeling of plots, create a vector holding descriptive names of the two sets.
# setLabels = c("Human_L5ET", "Human_L2IT", "Mouse_L5ET")
# shortLabels = c("Human_L5ET", "Human_L2IT", "Mouse_L5ET")
# # Form multi-set expression data: columns starting from 9 contain actual expression data.
# multiExpr = vector(mode = "list", length = nSets)
# 
# multiExpr[[1]] = list(data = as.data.frame(human_l5et))
# names(multiExpr[[1]]$data) = names(human_l5et)
# rownames(multiExpr[[1]]$data) = row.names(human_l5et)
# 
# multiExpr[[2]] = list(data = as.data.frame(human_l2it))
# names(multiExpr[[2]]$data) = names(human_l2it)
# rownames(multiExpr[[2]]$data) = row.names(human_l2it)
# 
# 
# multiExpr[[3]] = list(data = as.data.frame(fus_nls))
# names(multiExpr[[3]]$data) = names(fus_nls)
# rownames(multiExpr[[3]]$data) = row.names(fus_nls)
# 
# # Check that the data has the correct format for many functions operating on multiple sets:
# exprSize = checkSets(multiExpr)
# 
# # #=====================================================================================
# # #
# # #  Code chunk 3
# # #
# # #=====================================================================================
# #
# #
# # Check that all genes and samples have sufficiently low numbers of missing values.
# gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
# gsg$allOK
# 
# #
# # #=====================================================================================
# # #
# # #  Code chunk 4
# # #
# # #=====================================================================================
# #
# 
# if (!gsg$allOK)
# {
#   # Print information about the removed genes:
#   if (sum(!gsg$goodGenes) > 0)
#     printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
#                                               collapse = ", ")))
#   for (set in 1:exprSize$nSets)
#   {
#     if (sum(!gsg$goodSamples[[set]]))
#       printFlush(paste("In set", setLabels[set], "removing samples",
#                        paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
#     # Remove the offending genes and samples
#     multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
#   }
#   # Update exprSize
#   exprSize = checkSets(multiExpr)
# }
# 
# 
# #=====================================================================================
# # #
# # #  Code chunk 5
# # #
# # #=====================================================================================
# #
# #
# sampleTrees = list()
# for (set in 1:nSets)
# {
#   sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
# }
# 
# #
# # #=====================================================================================
# # #
# # #  Code chunk 6
# # #
# # #=====================================================================================
# #
# #
# #pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
# par(mfrow=c(2,1))
# par(mar = c(0, 4, 2, 0))
# for (set in 1:nSets)
#   plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
#        xlab="", sub="", cex = 0.7);
# 
# 
# # #=====================================================================================
# # #
# # #  Code chunk 9
# # #
# # #=====================================================================================
# #
# human_l5et_traits = human_l5et_traits[,c(2,5)]
# human_l5et_traits$Sex = "1"
# colnames(human_l5et_traits) = c("genotypes", "batch")
# 
# human_l2it_traits = human_l2it_traits[,c(2,5)]
# human_l2it_traits$Sex = "1"
# colnames(human_l2it_traits) = c("genotypes", "batch")
# 
# 
# traitData = rbind(human_l5et_traits,human_l2it_traits,fus_nls_traits)
# traitData$Samples = row.names(traitData)
# traitData$Species = c(rep("Human",68), rep("Mouse", 29))
# 
# # remove columns that hold information we do not need.
# allTraits = traitData
# # See how big the traits are and what are the trait and sample names
# dim(allTraits)
# names(allTraits)
# allTraits$Samples
# # Form a multi-set structure that will hold the clinical traits.
# Traits = vector(mode="list", length = nSets);
# for (set in 1:nSets)
# {
#   setSamples = rownames(multiExpr[[set]]$data);
#   traitRows = match(setSamples, allTraits$Samples);
#   Traits[[set]] = list(data = allTraits[traitRows, -2]);
#   rownames(Traits[[set]]$data) = allTraits[traitRows, 3];
# }
# collectGarbage();
# # Define data set dimensions
# nGenes = exprSize$nGenes;
# nSamples = exprSize$nSamples;
# #
# #
# # #=====================================================================================
# # #
# # #  Code chunk 10
# # #
# # #=====================================================================================
# #
# #
# save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
#       file = "~/Desktop/Analysis_fus/Exon_Intron/sod_christine/Consensus_matrix_L2IT_L5ET_crossSpecies.rda");
# 

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

setwd("~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/")
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "WGCNA/Consensus_matrix_L2IT_crossSpecies.rda");
#The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets
genes = colnames(multiExpr[[1]]$data)

load("WGCNA/TOMmatrix_L2IT.rda")
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
allowWGCNAThreads()

# powers = c(seq(1,10,by=1), seq(12,20, by=2))
# # Initialize a list to hold the results of scale-free analysis
# powerTables = vector(mode = "list", length = nSets)
# # Call the network topology analysis function for each set in turn
# for (set in 1:nSets)
#   powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, networkType = "signed",
#                                                      verbose = 2)[[2]])
# 
# collectGarbage();
# # Plot the results:
# colors = c("black", "red")
# # Will plot these columns of the returned scale free analysis tables
# plotCols = c(2,5,6,7)
# colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
#              "Max connectivity");
# # Get the minima and maxima of the plotted points
# ylim = matrix(NA, nrow = 2, ncol = 4);
# for (set in 1:nSets)
# {
#   for (col in 1:length(plotCols))
#   {
#     ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
#     ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
#   }
# }
# # Plot the quantities in the chosen columns vs. the soft thresholding power
# sizeGrWindow(8, 6)
# par(mfcol = c(2,2));
# par(mar = c(4.2, 4.2 , 2.2, 0.5))
# cex1 = 0.7;
# for (col in 1:length(plotCols)) for (set in 1:nSets)
# {
#   if (set==1)
#   {
#     plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
#          xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
#          main = colNames[col]);
#     addGrid();
#   }
#   if (col==1)
#   {
#     text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
#          labels=powers,cex=cex1,col=colors[set]);
#   } else
#     text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
#          labels=powers,cex=cex1,col=colors[set]);
#   if (col==1)
#   {
#     legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
#   } else
#     legend("topright", legend = setLabels, col = colors, pch = 20) ;
# }
# 
# consTOM <- consensusTOM(multiExpr,
#                         checkMissingData = TRUE,
#                         maxBlockSize = Inf,
#                         nPreclusteringCenters = NULL,
#                         randomSeed = 12345,
#                         corType = "pearson",
#                         maxPOutliers = 0.05,
#                         quickCor = 0,
#                         pearsonFallback = "individual",
#                         power = 14,
#                         networkType = "signed",
#                         checkPower = TRUE,
#                         TOMType = "signed",
#                         TOMDenom = "min",
#                         saveIndividualTOMs = TRUE,
#                         individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
#                         networkCalibration = "full quantile",
#                         saveCalibratedIndividualTOMs = TRUE,
#                         calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set%s-Block%b.RData",
#                         calibrationQuantile = 0.95,
#                         sampleForCalibration = TRUE, sampleForCalibrationFactor = 5000,
#                         getNetworkCalibrationSamples = FALSE,
#                         consensusQuantile = 0,
#                         useMean = FALSE,
#                         setWeights = NULL,
#                         saveConsensusTOMs = TRUE,
#                         consensusTOMFilePattern = "consensusTOM-Block%b.RData",
#                         returnTOMs = TRUE,
#                         nThreads = 10,
#                         verbose = 5)
# 
# 
# 
# TOM<-as.matrix(consTOM[["consensusTOM"]][[1]]) # get the consensus TOM
# rownames(TOM)<-colnames(multiExpr[[1]]$data)
# colnames(TOM)<-colnames(multiExpr[[2]]$data)

consTree = hclust(as.dist(1-TOM), method = "average");

minModuleSize = 50;

dynamicMods = cutreeDynamic(dendro = consTree, distM = 1-TOM,
                       deepSplit = 2, cutHeight = 0.99,
                       minClusterSize = minModuleSize,
                       pamRespectsDendro = FALSE );

table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

pdf('~/Desktop/single_cell/sc_als/processed_data/figures_raw/dendroclors.pdf')
plotDendroAndColors(consTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs0 = moduleEigengenes(multiExpr[[1]]$data, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(Traits[[1]]$data$genotypes);
names(weight) = "weight"
weight = as.numeric(weight$weight)
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(multiExpr[[1]]$data, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(multiExpr[[1]]$data, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# Create the starting data frame
geneInfo0 = data.frame(geneSymbol = genes,
                       moduleColor = dynamicColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.));
geneInfo = geneInfo0[geneOrder, ]



##### Calculate module membeship for each of the transcription factors #########
source('WGCNA/topHubs.R')
tf <- c("YBX1", "SP2", "CREB1", "ZHX2", "DR1", "BRF2", "CREB3", "MAFG", "JUN", "ZBTB24", "HMG20B", "KLF6", "E2F5", "THAP11", "MLX", "SNAPC5", "CEBPG", "DBP")
#tf_turquoise = geneInfo0[geneInfo0$geneSymbol %in% tf,]
#tf_turquoise = geneInfo0[geneInfo0$moduleColor == "turquoise",]
#tf_turquoise = tf_turquoise[,c("MM.turquoise", "p.MM.turquoise")]
modulekME = signedKME(multiExpr[[1]]$data,MEs)
modulekME$gene = row.names(modulekME)


tf_turquoise = modulekME[modulekME$gene %in% tf,]
gene_sig_turquoise = colnames(multiExpr[[1]]$data[,dynamicColors=="turquoise"])

tf_turquoise = modulekME[modulekME$gene %in% gene_sig_turquoise,]
tf_turquoise = tf_turquoise[,c("kMEturquoise", "gene")]
tf_turquoise <- tf_turquoise[order(tf_turquoise$kMEturquoise,decreasing=TRUE),]
tf_turquoise$rank <- 1:nrow(tf_turquoise)
tf_turquoise$colors <- ifelse(tf_turquoise$gene == "CREB3", "#8DD3C7","#969696")

library(ggrepel)

tf_turquoise$gene <- factor(tf_turquoise$gene,levels = unique(tf_turquoise$gene))

p0 <- ggplot(tf_turquoise, aes(x=rank, y=kMEturquoise)) + geom_point()

ggsave(filename = 'figures_raw/module_membership_overmodule.pdf',p0)



p<- ggplot(tf_turquoise, aes(x=gene, y=kMEturquoise, fill=colors)) +
  geom_col() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

write.csv(tf_turquoise,file = "article_marques/mme_data.csv", quote = F)
ggsave(filename = 'figures_raw/module_membership.pdf',p)



#new_modules = mergeCloseModules(exprData =multiExpr[[1]]$data, MEs = MEs,colors = dynamicColors,corFnc = "cor", corOptions = list(use = 'p'),cutHeight = 0.05 )
#new_colors = new_modules$colors 

#####
#MEs_l5et = moduleEigengenes(multiExpr[[1]]$data, dynamicColors)$eigengenes
MEs_l2it = moduleEigengenes(multiExpr[[1]]$data, dynamicColors)$eigengenes
MEs_mouse = moduleEigengenes(multiExpr[[2]]$data, dynamicColors)$eigengenes
#moduleTraitCor_l5et = cor(MEs_l5et, as.numeric(Traits[[1]]$data$genotypes), use = "p")
moduleTraitCor_l2it = cor(MEs_l2it,  as.numeric(Traits[[1]]$data$genotypes), use = "p")
moduleTraitCor_mouse = cor(MEs_mouse,  as.numeric(Traits[[2]]$data$genotypes), use = "p")
#l5et_nSamples = nrow(multiExpr[[1]]$data)
l2it_nSamples = nrow(multiExpr[[1]]$data)
mouse_nSamples = nrow(multiExpr[[2]]$data)
#moduleTraitPvalue_l5et = corPvalueStudent(moduleTraitCor_l5et, l5et_nSamples)
moduleTraitPvalue_l2it = corPvalueStudent(moduleTraitCor_l2it, l2it_nSamples)
moduleTraitPvalue_mouse = corPvalueStudent(moduleTraitCor_mouse, mouse_nSamples)
#moduleTraitPadj_l5et = as.data.frame(p.adjust(moduleTraitPvalue_l5et, method="fdr"))
moduleTraitPadj_l2it = as.data.frame(p.adjust(moduleTraitPvalue_l2it, method="fdr"))
moduleTraitPadj_mouse = as.data.frame(p.adjust(moduleTraitPvalue_mouse, method="fdr"))
#row.names(moduleTraitPadj_l5et) = row.names(moduleTraitPvalue_l5et)
row.names(moduleTraitPadj_l2it) = row.names(moduleTraitPvalue_l2it)
row.names(moduleTraitPadj_mouse) = row.names(moduleTraitPvalue_mouse)
#colnames(moduleTraitPadj_l5et) = "V1"
colnames(moduleTraitPadj_l2it) = "V1"
colnames(moduleTraitPadj_mouse) = "V1"

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor_l2it), ncol(moduleTraitCor_l2it));
consensusPvalue = matrix(NA, nrow(moduleTraitCor_mouse), ncol(moduleTraitCor_mouse));
# Find consensus negative correlations
negative = moduleTraitCor_l2it < 0 & moduleTraitCor_mouse < 0;
consensusCor[negative] = pmax(moduleTraitCor_l2it[negative], moduleTraitCor_mouse[negative]);
consensusPvalue[negative] = pmax( moduleTraitPadj_l2it[negative],moduleTraitPadj_mouse[negative]);
# Find consensus positive correlations
positive = moduleTraitCor_l2it > 0 & moduleTraitCor_mouse > 0;
consensusCor[positive] = pmin(moduleTraitCor_l2it[positive],moduleTraitCor_mouse[positive]);
consensusPvalue[positive] = pmax( moduleTraitPadj_l2it[positive],moduleTraitPadj_mouse[positive]);

textMatrix = paste(signif(consensusCor, 2), "\n(",signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor_l2it)

sizeGrWindow(10,7)
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = c("Consensus"),
               yLabels = names(MEs_l2it),
               ySymbols = names(MEs_l2it),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))


df_cor = data.frame(row.names = row.names(moduleTraitCor_l2it),consensusCor,consensusPvalue)
df_cor = na.omit(df_cor)
df_cor$Module = row.names(df_cor)
df_cor$Module <- factor(df_cor$Module, levels = df_cor$Module[order(df_cor$consensusCor,decreasing = TRUE)])
write.table(df_cor,file = "processed_data/article_marques/WGCNA_correlationToPheno.txt", sep = "\t", quote = F, col.names = T, row.names = F)

pdf('~/Desktop/single_cell/sc_als/processed_data/figures_raw/histogram_module.pdf')
p<- ggplot(df_cor, aes(x=reorder(Module,-consensusCor), y=consensusCor, fill="Neuonal Module")) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +  ylim(-0.8, 0.8)
p + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


dev.off()

gene_names = colnames(multiExpr[[1]]$data)

gene_sig_turquoise = colnames(multiExpr[[1]]$data[,dynamicColors=="turquoise"])
gene_sig_yellow = colnames(multiExpr[[1]]$data[,dynamicColors=="lightyellow"])


turquoise= as.data.frame(colnames(multiExpr[[1]]$data[,dynamicColors=="turquoise"]))
turquoise$module = "turquoise"
colnames(turquoise)[1] <- "gene"
lightyellow= as.data.frame(colnames(multiExpr[[1]]$data[,dynamicColors=="lightyellow"]))
lightyellow$module = "lightyellow"
colnames(lightyellow)[1] <- "gene"

combined_module = rbind(turquoise,lightyellow)
write.table(combined_module,file = '~/Desktop/single_cell/sc_als/processed_data/article_marques/gene_list_modules.txt', quote = F,col.names = T,row.names = F)

#write.table(gene_sig_blue,file = '~/Desktop/BLUE', quote = F, row.names = F, col.names = F)
############################################################                    
library(gplots)
col_scale = colorpanel(100, 'purple','black','yellow')
heatmap.2(t(multiExpr[[1]]$data[,dynamicColors=="turquoise"]), scale = "row", 
          col=col_scale, density.info ="none", trace="none", 
          cexCol=0.7, cexRow=0.8, margin=c(19,11), main = "LoF Turquoise Module", Colv = T)
heatmap.2(t(multiExpr[[1]]$data[,dynamicColors=="lightyellow"]), scale = "row", 
          col=col_scale, density.info ="none", trace="none", 
          cexCol=0.7, cexRow=0.8, margin=c(19,11), main = "LoF Turquoise Module", Colv = T)


#######################
nSelect = 1000
set.seed(10);
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM1 = 1-TOM  #Calculate topological overlap matrix for datExpr1
select = sample(nGenes, size = nSelect);
selectTOM1 = dissTOM1[select, select];
diag(selectTOM1)=NA
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM1), method = "average")
selectColors = dynamicColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
myheatcol = colorpanel(250,'red',"yellow",'black')
plotDiss = selectTOM1^7;

pdf('~/Desktop/single_cell/sc_als/processed_data/figures_raw/tomplot.pdf')
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=myheatcol)
dev.off()
#########################################################################
#########################################################################
################## Perform enrichment analysis for modules ##############
#########################################################################
#########################################################################

library(enrichR)
dbs <- c("KEGG_2016", "Reactome_2016","GO_Biological_Process_2019")
dbs2 <- c("GO_Biological_Process_2023", "Reactome_2016","GO_Molecular_Function_2023","KEGG_2016")
enrichrT <- enrichr(gene_sig_turquoise, dbs2)
enrichrY <- enrichr(gene_sig_yellow, dbs2)
#printEnrich(enrichr5, "~/Desktop/Analysis_fus/LongGene_code/Imma_paper/GO_analysis/enrichr5.txt", sep = "\t", columns = c(1,2,4,9))
#printEnrich(enrichr22, "~/Desktop/Analysis_fus/LongGene_code/Imma_paper/GO_analysis/enrichr22.txt", sep = "\t", columns = c(1,2,4,9))

## convert list to dataframe for enrichR sets
enrichrT <- rbindlist(enrichrT, fill=TRUE)
enrichrY <- rbindlist(enrichrY, fill=TRUE)

enrichrT = enrichrT[order(enrichrT$Adjusted.P.value, decreasing = FALSE),]
enrichrY = enrichrY[order(enrichrY$Adjusted.P.value, decreasing = FALSE),]

enrichrT = enrichrT[1:30,]
enrichrY = enrichrY[1:30,]

enrichrT$Adjusted.P.value = -log(enrichrT$Adjusted.P.value)
enrichrY$Adjusted.P.value = -log(enrichrY$Adjusted.P.value)


enrichrT$Term = gsub("\\s*\\([^\\)]+\\)","",enrichrT$Term)
enrichrY$Term = gsub("\\s*\\([^\\)]+\\)","",enrichrY$Term)

# Basic histogram for turquoise module
enrichrT$Term <- factor(enrichrT$Term, levels = enrichrT$Term[order(enrichrT$Adjusted.P.value)])
p<- ggplot(enrichrT, aes(x=Term, y=Adjusted.P.value, fill="Turquoise Module")) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +  ylim(0, 50)
p + theme_bw()
p + coord_flip()
# Basic histogram for yellow module
enrichrY$Term <- factor(enrichrY$Term, levels = enrichrY$Term[order(enrichrY$Adjusted.P.value)])

p<- ggplot(enrichrY, aes(x=Term, y=Adjusted.P.value, fill="Yellow Module")) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +  ylim(0, 5)
p + theme_bw()
p + coord_flip()


##################################################################
##################################################################
##################### Network ####################################
##################################################################
##################################################################
library(hdWGCNA)
library(Seurat)
library(Matrix)
data_path <- '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/'
load(paste0(data_path,"processed_data/human_10x_10XALS_integration/Figure_3/seurat_combat.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/integrated_seurat/markers.cells.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/Figure_1/aucell_metadata.3.rda"))
load(paste0(data_path, "processed_data/human_10x_10XALS_integration/Figure_2/AggregatedCountsCombat.rda"))

sce.to.seurat
#sce.to.seurat.als = subset(sce.to.seurat, subset = Condition != 'FTLD')

seurat_wgcna_ALS <- SetupForWGCNA(
  sce.to.seurat,
  gene_select = "custom", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "ALS_Control", features = geneInfo0$geneSymbol  # the name of the hdWGCNA experiment
)

seurat_wgcna_ALS <- MetacellsByGroups(
  seurat_obj = seurat_wgcna_ALS,
  group.by = c("class_label"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 50, # nearest-neighbors parameter
  ident.group = 'class_label', # set the Idents of the metacell seurat object
  target_metacells = 1000,
  min_cells = 100, verbose = TRUE
)

### normalize metacell expression matrix:
#seurat_wgcna_ALS <- NormalizeMetacells(seurat_wgcna_ALS)
#seurat_wgcna_ALS <- ScaleMetacells(seurat_wgcna_ALS, features=VariableFeatures(seurat_wgcna_ALS))
seurat_wgcna_ALS_sub <- SetDatExpr(
  seurat_wgcna_ALS,
  group_name = c("Glutamatergic"), # the name of the group of interest in the group.by column
  group.by='class_label', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)


####
datExpr = seurat_wgcna_ALS_sub@misc$ALS_Control$datExpr
TOM = TOMsimilarityFromExpr(datExpr = datExpr,power = 12)
colnames(TOM) = colnames(datExpr)
row.names(TOM) = colnames(datExpr)

#### Eigengenes + Connectivity ######
colors_filtered = geneInfo0$geneSymbol %in% colnames(datExpr)
color_final = geneInfo0[colors_filtered,2]
eigengenes = moduleEigengenes(datExpr, 
                              color_final, 
                              impute = TRUE, 
                              nPC = 1, 
                              align = "along average", 
                              excludeGrey = FALSE, 
                              grey = if (is.numeric(colors)) 0 else "grey",
                              subHubs = TRUE,
                              trapErrors = FALSE, 
                              softPower = 12,
                              scale = TRUE,
                              verbose = 0, indent = 0)
eigengenes_exp =eigengenes$eigengenes

target = geneInfo0[,c(1,2,2)]
colnames(target) = c("gene_name", "module", "color")
modules <- target[match(colnames(seurat_wgcna_ALS_sub@misc$ALS_Control$datExpr), target$gene_name),]
module_keep = df_cor$Module
module_keep  = gsub("ME", "", module_keep)
module_keep = c("blue", "turquoise", "lightyellow", "darkgreen", "pink", "darkgrey")

kMEs <- WGCNA::signedKME(
  datExpr,
  eigengenes_exp,
  outputColumnName = "kME",
  corFnc = 'cor', corOptions = "use='p'"
)


modules <- modules[,1:3]
colnames(kMEs) <- colnames(eigengenes_exp)
colnames(kMEs) <- paste0("kME_", colnames(kMEs))
kMEs <- cbind(modules, kMEs)
kMEs <- subset(kMEs, module %in% module_keep)


modules$module = as.factor(modules$module)
mods <- levels(modules$module)
mods <- mods[mods %in% module_keep]
modules <- subset(modules, module %in% module_keep)



# get hub genes:
connectivity_table= topHubs(multiExpr[[1]]$data, colorh=  dynamicColors, power=2, type="signed")
connectivity_table = connectivity_table[,c(1,2)]
colnames(connectivity_table) = c("gene_name", "connectivity")
modules = merge(modules, connectivity_table,by=c("gene_name"))
hub_list <- lapply(mods, function(cur_mod) {
  cur <- subset(modules, module == cur_mod)
  cur = cur %>% group_by(module) %>% top_n( 2, wt = connectivity)
  cur = cur$gene_name
  
})
names(hub_list) <- mods

# get all genes that aren't in gray mod
selected_genes <- modules[modules$module %in% mods,'gene_name']

# subset the TOM for umap
# keep all genes as rows, and keep only hubs as cols (features)
feature_mat <- TOM[selected_genes,unlist(hub_list)]

# print('running supervised UMAP:')
# hub_umap <-  uwot::umap(
#   X = feature_mat,
#   min_dist = 0.4,
#   n_neighbors= 25,
#   metric = "cosine",
#   spread=1)

print('running supervised UMAP:')
hub_umap <-  uwot::umap(
    X = feature_mat,
    min_dist = 0.4,
    n_neighbors= 50,
    metric = "cosine",
    spread=0.8,
    y = modules$module, target_weight = 0.60,n_epochs = 200# for supervised UMAP
  )


# set up plotting df
plot_df <- as.data.frame(hub_umap)
colnames(plot_df) <- c("UMAP1", "UMAP2")
plot_df$gene <- rownames(feature_mat)

# add module color, and hub gene status to the plotting df:
ix <- match(plot_df$gene, modules$gene_name)
plot_df$module <- modules$module[ix]
plot_df$color <- modules$color[ix]
plot_df$hub <- ifelse(
  plot_df$gene %in% as.character(unlist(hub_list)), 'hub', 'other'
)

# get kME values for each gene
kMEs <- do.call(rbind, lapply(mods, function(cur_mod){
  cur <- subset(kMEs, module == cur_mod)
  cur <- cur[,c('gene_name', paste0('kME_',"ME", cur_mod))]
  colnames(cur) <- c('gene_name', 'kME')
  
  # scale kMEs between 0 & 1:
  cur$kME <- scale01(cur$kME)
  cur
}))

ix <- kMEs$gene_name[match(plot_df$gene, kMEs$gene_name)]
plot_df$kME <- kMEs[ix, 'kME']

umap_df = plot_df

#### Change light yellow to darkyellow for visibility #####
umap_df$module = ifelse(umap_df$module == "lightyellow", "yellow",umap_df$module )
umap_df$color = ifelse(umap_df$color == "lightyellow", "yellow",umap_df$color )


# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()


save(TOM,modules,umap_df,hub_list, file = "WGCNA/network_files.rda")
pdf(file = "figures_raw/network2.pdf")
umap_plot(currentTOM = TOM,currentModules = modules,currentUMAP = umap_df,currentHub = hub_list,edge_prop = 0.2,edge.alpha = 3,vertex.label.cex = 0.01,sample_edges = TRUE)
dev.off()


load("WGCNA/network_files.rda")
##################################################################
##################################################################
##################### Module preservation #########################
##################################################################
##################################################################
names(multiExpr) <- c("Human_L2IT", "Mouse_GFP")
multiColor = list(Human_L2IT = dynamicColors);
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
# Save the results
save(mp, file = "modulePreservation_mouse_human.RData");

load('WGCNA/modulePreservation_mouse_human.RData')


ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);


# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )


# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];



plotMods = !(modColors %in% c("grey", "gold"));


keep = !(colnames(MEs) %in% row.names(df_cor));
MEs_to_plot = MEs[,keep]
MEs_to_plot = colnames(MEs_to_plot)
MEs_to_plot = gsub("^ME", "", MEs_to_plot)
MEs_to_plot = c(MEs_to_plot,"gold")

plotMods = !(modColors %in% MEs_to_plot);


# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(5, 10);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))

pdf("figures_raw/Module_conservation.pdf", width=6, height=7)
for (p in 2:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off()


# If plotting into a file, close it
################################################################
######### dataset ALS-NYGC eigengenes in mouse and human  ######
################################################################
human_l2it_eigengenes = moduleEigengenes(multiExpr[[1]]$data, 
                                         dynamicColors, 
                            impute = TRUE, 
                            nPC = 1, 
                            align = "along average", 
                            excludeGrey = FALSE, 
                            grey = if (is.numeric(colors)) 0 else "grey",
                            subHubs = TRUE,
                            trapErrors = FALSE, 
                            softPower = 5,
                            scale = TRUE,
                            verbose = 0, indent = 0)


mouse_gfp_eigengenes = moduleEigengenes(multiExpr[[2]]$data, 
                                        dynamicColors, 
                                         impute = TRUE, 
                                         nPC = 1, 
                                         align = "along average", 
                                         excludeGrey = FALSE, 
                                         grey = if (is.numeric(colors)) 0 else "grey",
                                         subHubs = TRUE,
                                         trapErrors = FALSE, 
                                         softPower = 5,
                                         scale = TRUE,
                                         verbose = 0, indent = 0)


cell = "L5-ET"
load(paste0( "~/Desktop/iMac_U1118/Desktop/Analysis_fus/Exon_Intron/sod_christine/Consensus_WGNA/input_wgcna/", cell,"/","human_", cell,".WGCNA.rda"))
datExpr = datExpr[,colnames(datExpr) %in% genes]


colors_filtered = geneInfo0$geneSymbol %in% colnames(datExpr)
color_final = geneInfo0[colors_filtered,2]


eigengenes = moduleEigengenes(datExpr, 
                              color_final, 
                              impute = TRUE, 
                              nPC = 1, 
                              align = "along average", 
                              excludeGrey = FALSE, 
                              grey = if (is.numeric(colors)) 0 else "grey",
                              subHubs = TRUE,
                              trapErrors = FALSE, 
                              softPower = 5,
                              scale = TRUE,
                              verbose = 0, indent = 0)

eigengenes =eigengenes$averageExpr[,c("AEturquoise","AElightyellow", "AEdarkgrey", "AEdarkgreen")]
human_l2it_eigengenes =human_l2it_eigengenes$averageExpr[,c("AEturquoise","AElightyellow", "AEdarkgrey", "AEdarkgreen")]
mouse_gfp_eigengenes =mouse_gfp_eigengenes$averageExpr[,c("AEturquoise","AElightyellow", "AEdarkgrey", "AEdarkgreen")]




# Change box plot colors by groups
# Turquoise module
human_l2it_eigengenes = melt(as.matrix(human_l2it_eigengenes))
colnames(human_l2it_eigengenes) = c("group", "module", "eigengenes")
human_l2it_eigengenes$group = str_sub(human_l2it_eigengenes$group, 8, 9)
human_l2it_eigengenes$group = gsub("AL", "ALS", human_l2it_eigengenes$group)
human_l2it_eigengenes$cell = "L2/3-IT"
human_l2it_eigengenes$species = "human"


mouse_gfp_eigengenes = melt(as.matrix(mouse_gfp_eigengenes))
colnames(mouse_gfp_eigengenes) = c("group", "module", "eigengenes")
mouse_gfp_eigengenes$group  = gsub("\\_.*","",mouse_gfp_eigengenes$group)
mouse_gfp_eigengenes$group = gsub("^[0-9]|[0-9]$", "", mouse_gfp_eigengenes$group)
mouse_gfp_eigengenes$group = ifelse(mouse_gfp_eigengenes$group == "WT", "PN", ifelse(mouse_gfp_eigengenes$group == "SOD", "ALS", mouse_gfp_eigengenes$group))
mouse_gfp_eigengenes$cell = "CSN"
mouse_gfp_eigengenes$species = "mouse"


eigengenes = melt(as.matrix(eigengenes))
colnames(eigengenes) = c("group", "module", "eigengenes")
eigengenes$group = str_sub(eigengenes$group, 8, 9)
eigengenes$group = gsub("AL", "ALS", eigengenes$group)
colnames(eigengenes) = c("group", "module", "eigengenes")
eigengenes$cell = "L5-ET"
eigengenes$species = "human"


all_eigengenes = rbind(human_l2it_eigengenes,eigengenes,mouse_gfp_eigengenes)
all_eigengenes$group <- factor(all_eigengenes$group,levels = c("PN","ALS"))
all_eigengenes$cell <- factor(all_eigengenes$cell,levels = c("L2/3-IT","L5-ET", "CSN"))
write.table(all_eigengenes,file = "processed_data/article_marques/all.eigengenes_human", sep = "\t", quote = F, col.names = T, row.names = F)

turquoise = all_eigengenes[all_eigengenes$module == "AEturquoise",]
yellow = all_eigengenes[all_eigengenes$module == "AElightyellow",]

# Grouped
extrafont::loadfonts()

library(viridis)
library(hrbrthemes)

min = min(turquoise$eigengenes) - 0.1
max = max(turquoise$eigengenes) + 0.1

pdf("figures_raw/CellEigengenes_Turquoise.pdf")
turquoise %>%
  mutate(cell = factor(cell, levels = c("L2/3-IT","L5-ET","CSN"))) %>%
  ggplot(aes(fill=group, y=eigengenes, x=cell)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()
### perform 2 way ANOVA
library(multcomp)
library(rstatix)
library(ggpubr)
## get summary statistic
turquoise %>%
  group_by(group, cell) %>%
  get_summary_stats(eigengenes, type = "mean_sd")

res.aov2 <- aov(eigengenes ~ group*cell, data = turquoise)
summary(res.aov2)
TukeyHSD(res.aov2)


min = min(turquoise$eigengenes) - 0.1
max = max(turquoise$eigengenes) + 0.1

pdf("figures_raw/CellEigengenes_Yellow.pdf")
yellow %>%
  mutate(cell = factor(cell, levels = c("L2/3-IT","L5-ET","CSN"))) %>%
  ggplot(aes(fill=group, y=eigengenes, x=cell)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()
### perform 2 way ANOVA
library(multcomp)
library(rstatix)
library(ggpubr)
## get summary statistic
yellow %>%
  group_by(group, cell) %>%
  get_summary_stats(eigengenes, type = "mean_sd")

res.aov3 <- aov(eigengenes ~ group*cell, data = yellow)
summary(res.aov3)
TukeyHSD(res.aov3)



######### Mouse GFP eigengenes time course ########
mouse_gfp_eigengene_course = moduleEigengenes(multiExpr[[2]]$data, 
                                        dynamicColors, 
                                        impute = TRUE, 
                                        nPC = 1, 
                                        align = "along average", 
                                        excludeGrey = FALSE, 
                                        grey = if (is.numeric(colors)) 0 else "grey",
                                        subHubs = TRUE,
                                        trapErrors = FALSE, 
                                        softPower = 5,
                                        scale = TRUE,
                                        verbose = 0, indent = 0)

mouse_gfp_eigengene_course =mouse_gfp_eigengene_course$averageExpr[,c("AEturquoise","AElightyellow")]
mouse_gfp_eigengene_course = melt(as.matrix(mouse_gfp_eigengene_course))
colnames(mouse_gfp_eigengene_course) = c("group", "module", "eigengenes")
mouse_gfp_eigengene_course$time = ifelse(grepl("batch1", mouse_gfp_eigengene_course$group), "30", ifelse(grepl("batch2", mouse_gfp_eigengene_course$group),"60", ifelse(grepl("batch3", mouse_gfp_eigengene_course$group),"90", "105")))
mouse_gfp_eigengene_course$group  = gsub("\\_.*","",mouse_gfp_eigengene_course$group)
mouse_gfp_eigengene_course$group = gsub("^[0-9]|[0-9]$", "", mouse_gfp_eigengene_course$group)
mouse_gfp_eigengene_course$cell = "CSN"
mouse_gfp_eigengene_course$species = "mouse"

mouse_gfp_eigengene_course$time <- factor(mouse_gfp_eigengene_course$time,levels = c("30","60", "90", "105"))
#write.table(mouse_gfp_eigengene_course,file = "processed_data/article_marques/mouse.eigengenes_turquoise", sep = "\t", quote = F, col.names = T, row.names = F)


min = min(mouse_gfp_eigengene_course$eigengenes) - 0.1
max = max(mouse_gfp_eigengene_course$eigengenes) + 0.1

mouse_gfp_eigengene_course$module = gsub("AElightyellow", "AEturquoise", mouse_gfp_eigengene_course$module)


mouse_gfp_eigengene_course$color = ifelse(mouse_gfp_eigengene_course$group == "WT", "#6A51A3", "#8DD3C7")

pdf("figures_raw/CellEigengenes_mouseCourse.pdf", width = 16, height = 8)
mouse_gfp_eigengene_course$time = as.numeric(mouse_gfp_eigengene_course$time)
ggplot(mouse_gfp_eigengene_course, aes(x=time, y=eigengenes, color=group)) +
  geom_point(alpha = 0.5, color = mouse_gfp_eigengene_course$color) + 
  geom_smooth(method=lm, level=0.5,aes(color=group,fill=group)) + scale_color_manual(values=c("#8DD3C7","#6A51A3")) +
  scale_fill_manual(values=c("#8DD3C7","#6A51A3"))
dev.off()

### perform 2 way ANOVA
library(multcomp)
library(rstatix)
library(ggpubr)
## get summary statistic
mouse_gfp_eigengene_course %>%
  group_by(group, time) %>%
  get_summary_stats(eigengenes, type = "mean_sd")

mouse_gfp_eigengene_course$time = as.factor(mouse_gfp_eigengene_course$time)
res.aov2 <- aov(eigengenes ~ group*time, data = mouse_gfp_eigengene_course)
summary(res.aov2)
res <- TukeyHSD(res.aov2, "group:time", ordered = T)
res <- as.data.frame(res$`group:time`)
res$comp = row.names(res)
write.table(res,file = "article_marques/manuscrit_v9/mouse_stats_timecourse.txt", quote = F, sep = "\t", row.names = F)


###############################################
###############################################
###### Merge eigengenes in cell class #########
###############################################
###############################################
merged_eigengenes <- 
  do.call(rbind,
          lapply(list.files(path = "WGCNA/output_wgcna_class/", recursive = T, pattern = "\\.txt$", full.names = T), read.table))
colnames(merged_eigengenes) = c("group", "module", "eigengenes", "cell")
merged_eigengenes$group <- factor(merged_eigengenes$group,levels = c("PN","ALS"))


####### Grouped for turquoise module ###########
library(viridis)
library(hrbrthemes)
merged_eigengenes_turquoise = merged_eigengenes[merged_eigengenes$module == "AEturquoise",]
min = min(merged_eigengenes_turquoise$eigengenes) - 0.1
max = max(merged_eigengenes_turquoise$eigengenes) + 0.1
p1 <- merged_eigengenes_turquoise %>%
  mutate(module = factor(cell, levels = c("Astrocyte","Endothelial","Excitatory", "Inhibitory", "Microglia", "Oligodendrocyte", "OPC"))) %>%
  ggplot(aes(fill=group, y=eigengenes, x=cell)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylim(-1.5,1.5) +
  ylab("Eigengenes") 

### perform 2 way ANOVA
library(tidyverse)
library(rstatix)
library(ggpubr)
## get summary statistic
merged_eigengenes_turquoise %>%
  group_by(group, cell) %>%
  get_summary_stats(eigengenes, type = "mean_sd")

res.aov_turquoise <- aov(eigengenes ~ group*cell, data = merged_eigengenes_turquoise)
summary(res.aov_turquoise)

pwc <- merged_eigengenes_turquoise %>% tukey_hsd(eigengenes ~ group * cell)


pwc <- merged_eigengenes %>%
  group_by(cell) %>%
  wilcox_test(eigengenes ~ group, p.adjust.method = "Bonferroni") 


write.table(pwc,file = "~/Desktop/single_cell/sc_als/processed_data/article_marques/eigengenes_stats_CellClass_turquoise.txt", quote = F, row.names = F, col.names = T)


####### Grouped for lightyellow module ###########
merged_eigengenes_yellow = merged_eigengenes[merged_eigengenes$module == "AElightyellow",]
min = min(merged_eigengenes_yellow$eigengenes) - 0.1
max = max(merged_eigengenes_yellow$eigengenes) + 0.1
pdf("figures_raw/ClassEigengenes_Turquoise.pdf")
p2 <- merged_eigengenes_yellow %>%
  mutate(module = factor(cell, levels = c("Astrocyte","Endothelial","Excitatory", "Inhibitory", "Microglia", "Oligodendrocyte", "OPC"))) %>%
  ggplot(aes(fill=group, y=eigengenes, x=cell)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylim(-1.5,1.5) +
  ylab("Eigengenes") 

### perform 2 way ANOVA
library(tidyverse)
library(rstatix)
library(ggpubr)
## get summary statistic
merged_eigengenes_yellow %>%
  group_by(group, cell) %>%
  get_summary_stats(eigengenes, type = "mean_sd")

res.aov_turquoise <- aov(eigengenes ~ group*cell, data = merged_eigengenes_yellow)
summary(res.aov_turquoise)

pwc <- merged_eigengenes_yellow %>% tukey_hsd(eigengenes ~ group * cell)
write.table(pwc,file = "~/Desktop/single_cell/sc_als/processed_data/article_marques/eigengenes_stats_CellClass_yellow.txt", quote = F, row.names = F, col.names = T)

library(ggpubr)
pdf(file = "figures_raw/ClassEigengenes.pdf", width = 11, height = 11)
ggarrange(p1, p2, ncol = 1, nrow = 2)
dev.off()


###############################################
###############################################
###### Merge eigengenes in cell population ####
###############################################
###############################################
merged_eigengenes <- 
  do.call(rbind,
          lapply(list.files(path = "WGCNA/output_wgcna/", recursive = T, pattern = "\\.txt$", full.names = T), read.table))
colnames(merged_eigengenes) = c("group", "module", "eigengenes", "cell")
merged_eigengenes$group <- factor(merged_eigengenes$group,levels = c("PN","ALS"))

library(tidyverse)
library(rstatix)
library(ggpubr)
library(agricolae)

# q3 <- quantile(merged_eigengenes$eigengenes, 0.75)
# iqr <- IQR(merged_eigengenes$eigengenes)
# upper_bound <- q3 + 1.5*iqr
# q1 <- quantile(merged_eigengenes$eigengenes, 0.25)
# lower_bound <- q1 - 1.5*iqr
# 
# merged_eigengenes <- merged_eigengenes[merged_eigengenes$eigengenes < upper_bound & merged_eigengenes$eigengenes > lower_bound,]



summary_all = merged_eigengenes %>%
  group_by(group, cell) %>%
  get_summary_stats(eigengenes, type = "mean_sd")


res.aov <- aov(eigengenes ~ group*cell, data = merged_eigengenes)
summary(res.aov)

pwc <- merged_eigengenes %>%
  group_by(cell) %>%
  wilcox_test(eigengenes ~ group, p.adjust.method = "Bonferroni") 


#res <- TukeyHSD(res.aov, "group:cell", ordered = T)
#res <- as.data.frame(res$`group:cell`)
#res$comp = row.names(res)
write.table(pwc,file = 'article_marques/eigengenes_cell_group_tukey.txt', quote = F, sep = "\t", row.names = F)

dend_order <- labels(dend.labeled)
order.dendo = as.data.frame(dend_order)
colnames(order.dendo)[1] = "cluster"
order.dendo$position = row.names(order.dendo)
order = order.dendo
colnames(order) <- c("cell","position")

to_plot_1 = merged_eigengenes
to_plot_1 = merge(to_plot_1,order,by=c("cell"))
to_plot_1= to_plot_1[order(as.numeric(as.character(to_plot_1$position))),]
#to_plot_1 = to_plot_1[to_plot_1$Group != "FTLD",]
to_plot_1$group <- factor(to_plot_1$group,levels = c("PN","ALS"))
to_plot_1$cell <- factor(to_plot_1$cell,levels = unique(to_plot_1$cell))

# pwc <- merged_eigengenes %>%
#    group_by(cell) %>%
#   wilcox_test(eigengenes ~ group, p.adjust.method = "Bonferroni") 

final_eigen_class = aggregate(merged_eigengenes[, 3], list(merged_eigengenes$cell,merged_eigengenes$group), mean)
colnames(final_eigen_class) = c("cell","Group","Zscore")

# clust = final_eigen_class
# clust = as.data.frame.matrix(xtabs(Zscore ~ cell + Group, clust))
# row.order <- hclust(dist(clust))$order
# clust = clust[row.order, ]
# clust$cell = row.names(clust)
# rownames(clust) <- NULL
# clust$position = row.names(clust)
# order = clust[,c(3,4)]

dend_order <- labels(dend.labeled)
order.dendo = as.data.frame(dend_order)
colnames(order.dendo)[1] = "cluster"
order.dendo$position = row.names(order.dendo)
order = order.dendo
colnames(order) <- c("cell","position")

to_plot_1 = final_eigen_class
to_plot_1 = merge(to_plot_1,order,by=c("cell"))
to_plot_1= to_plot_1[order(as.numeric(as.character(to_plot_1$position))),]
#to_plot_1 = to_plot_1[to_plot_1$Group != "FTLD",]
to_plot_1$Group <- factor(to_plot_1$Group,levels = c("PN","ALS"))
to_plot_1$cell <- factor(to_plot_1$cell,levels = unique(to_plot_1$cell))

p1b <- to_plot_1 %>%
  ggplot(aes(x = cell, y = Group, fill = Zscore)) +
  geom_tile(color = "black", lwd=0.5,linetype=1) +
  #scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0,limits=range(to_plot_1$Zscore)) +
  
  # scale_fill_gradientn(
  #   colors=c("blue","white","red"),
  #   values=scales::rescale(c(-0.2,0,0.2)),
  #   limits=c(-0.5,0.5)
  # ) +
  
  labs(fill = "Turquoise Z-score", 
       y = "Group",
       x = "") +
  
  theme(aspect.ratio = 2/6*length(order),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank())


p1b+ scale_fill_gradientn(colours = rev(brewer.pal(n=5,name ="RdBu")),values=scales::rescale(c(-0.3,-0.1,0,0.1,0.2)),limits=c(-0.5,0.3))
ggsave(filename = 'figures_raw/eigengenes_cell_groups.pdf',p1b+ scale_fill_gradientn(colours = rev(brewer.pal(n=5,name ="RdBu")),values=scales::rescale(c(-0.5,-0.1,0,0.1,0.3)),limits=c(-0.5,0.3)))




#############################################################################
#############################################################################
########  Variance in Aggregating Seurat all DEG-cell group #################
#############################################################################
#############################################################################
sce.to.seurat

object = subset(sce.to.seurat,subset = Condition ==  "ALS")
pb.method = "average"
assays = NULL
features = NULL
return.seurat = FALSE
group.by = c('cluster_label.dendro.deg')
add.ident = NULL
slot = "counts"
verbose = TRUE

var = data.return$RNA
var = reshape2::melt(var)


library(viridis)
library(hrbrthemes)

min = min(var$value)
max = max(var$value) + 10
pdf("processed_data/figures_raw/ClassEigengenes_Turquoise.pdf")
var %>%
  ggplot(aes(fill=Var2, y=value, x=Var2)) + 
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(min,max)
dev.off()
#group.by=c('cluster_label.dendro.deg', 'Sample_ID'), assays = 'RNA', slot = 'counts', return.seurat = FALSE


#############################################################################
#############################################################################
############## Get value for boxplot (median,quaritles....) #################
#############################################################################
#############################################################################
ggplot2_boxplot <- function(x){
  
  quartiles <- as.numeric(quantile(x,
                                   probs = c(0.25, 0.5, 0.75)))
  
  names(quartiles) <- c("25th percentile",
                        "50th percentile\n(median)",
                        "75th percentile")
  
  IQR <- diff(quartiles[c(1,3)])
  
  upper_whisker <- max(x[x < (quartiles[3] + 1.5 * IQR)])
  lower_whisker <- min(x[x > (quartiles[1] - 1.5 * IQR)])
  
  upper_dots <- x[x > (quartiles[3] + 1.5*IQR)]
  lower_dots <- x[x < (quartiles[1] - 1.5*IQR)]
  
  return(list("quartiles" = quartiles,
              "25th percentile" = as.numeric(quartiles[1]),
              "50th percentile\n(median)" = as.numeric(quartiles[2]),
              "75th percentile" = as.numeric(quartiles[3]),
              "IQR" = IQR,
              "upper_whisker" = upper_whisker,
              "lower_whisker" = lower_whisker,
              "upper_dots" = upper_dots,
              "lower_dots" = lower_dots))
}



data_1mWT <- ggplot2_boxplot(neuronMod$AEturquoise[1:6])
data_1mKI <- ggplot2_boxplot(neuronMod$AEturquoise[7:12])
data_6mWT <- ggplot2_boxplot(neuronMod$AEturquoise[13:18])
data_6mKI <- ggplot2_boxplot(neuronMod$AEturquoise[19:24])
data_22mWT <- ggplot2_boxplot(neuronMod$AEturquoise[25:27])
data_22mKI <- ggplot2_boxplot(neuronMod$AEturquoise[28:30])

all_turquoise = list(data_1mWT,data_1mKI,data_6mWT,data_6mKI,data_22mWT,data_22mKI)
all_turquoise = data.frame(matrix(unlist(all_turquoise), nrow=length(all_turquoise), byrow=T),stringsAsFactors=FALSE)
colnames(all_turquoise) = c('quartiles','25th_percentile','50th_percentile',"75th_percentile", 'IQR','upper_whisker','lower_whisker','upper_dots','lower_dots')
row.names(all_turquoise) = c('data_1mWT','data_1mKI','data_6mWT','data_6mKI','data_22mWT','data_22mKI')
write.table(all_turquoise, '~/Desktop/Analysis_fus/Exon_Intron/wgcna_intron/turquoiseMod.boxplot.txt', sep = '\t', quote = FALSE)

data_1mWT <- ggplot2_boxplot(YellowMod$AEyellow[1:6])
data_1mKI <- ggplot2_boxplot(YellowMod$AEyellow[7:12])
data_6mWT <- ggplot2_boxplot(YellowMod$AEyellow[13:18])
data_6mKI <- ggplot2_boxplot(YellowMod$AEyellow[19:24])
data_22mWT <- ggplot2_boxplot(YellowMod$AEyellow[25:27])
data_22mKI <- ggplot2_boxplot(YellowMod$AEyellow[28:30])

all_yellow = list(data_1mWT,data_1mKI,data_6mWT,data_6mKI,data_22mWT,data_22mKI)
all_yellow = data.frame(matrix(unlist(all_yellow), nrow=length(all_yellow), byrow=T),stringsAsFactors=FALSE)
all_yellow = all_yellow[,c(1:9)]
colnames(all_yellow)[1:9] = c('quartiles','25th_percentile','50th_percentile',"75th_percentile", 'IQR','upper_whisker','lower_whisker','upper_dots','lower_dots')
row.names(all_yellow) = c('data_1mWT','data_1mKI','data_6mWT','data_6mKI','data_22mWT','data_22mKI')
write.table(all_yellow, '~/Desktop/Analysis_fus/Exon_Intron/wgcna_intron/yellowMod.boxplot.txt', sep = '\t', quote = FALSE)






