## Packages needed to run code
## ----------------------------
library("dplyr")     ## data manipulation
library("data.table") ## read data, fast!
library("forcats")   ## dealing with factors
library("ggplot2")   ## dataviz
library("magrittr")  ## piping
library("metafolio")  ## colorpalette
library("skimr")     ## summarising data
library("qqman")     ## Manhattan plot
library("patchwork") ## assembling plots
library("CMplot")
library('plyr')
library(rstatix)
library(Vennerable)
library("ggpubr")
library(survival)
library(ggplot2)
library(reshape)
library(data.table)
library(dplyr)
library("ggpubr")
library(survival)
library(tidyverse)
library(viridis)
library(survivalAnalysis)
setwd('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_progression/clinic_timeseries/')

save.dir <- '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/figures_raw'

#### Load Df2 Mine data ######
surv_final <- read.table('survival_data_CREB3_pooled_cohort',header = T)
#write.table(surv_final,file = 'survival_data_CREB3_pooled_cohort')

lower_bound <- quantile(surv_final$MOD, 0.005)
upper_bound <- quantile(surv_final$MOD, 0.995)

outlier_ind <- which(surv_final$MOD < lower_bound & surv_final$MOD > upper_bound)
keep_ind <- which(surv_final$MOD > lower_bound & surv_final$MOD< upper_bound)

surv_final_nooutlier <- surv_final[keep_ind,]

library(ggsurvfit)
library(ggpubr)
m <- coxph(Surv(MOD, Status) ~ Genotype+AAO+Sex+Origin, data = surv_final_nooutlier) 
df_coxph <- cox_as_data_frame(
  m,
  unmangle_dict = NULL,
  factor_id_sep = ":",
  sort_by = NULL
)
write.table(df_coxph,file = '~/Desktop/iMac_U1118/Desktop/single_cell/sc_als/processed_data/article_marques/coxph_suppfigure11.txt', sep = '\t', row.names = F, quote = F)



p_geno <- survfit2(Surv(MOD, Status) ~ Genotype, data = surv_final_nooutlier) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) 

p_origin <- survfit2(Surv(MOD, Status) ~ Origin, data = surv_final_nooutlier) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) 


p_sex <- survfit2(Surv(MOD, Status) ~ Sex, data = surv_final_nooutlier) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) 

ggsave(p_geno, filename = paste(save.dir, 'survival_plot_geno.pdf', sep='/'), width = 5, height = 3.5)
ggsave(p_sex, filename = paste(save.dir, 'survival_plot_sex.pdf', sep='/'), width = 5, height = 3.5)
ggsave(p_origin, filename = paste(save.dir, 'survival_plot_origin.pdf', sep='/'), width = 5, height = 3.5)


aggregate(surv_final_nooutlier[, c("MOD")], list(surv_final_nooutlier$Genotype), mean)


surv_perm <- surv_final_nooutlier
surv_perm$Genotype = as.factor(surv_perm$Genotype)
print(ggplot(surv_perm,aes(Genotype,MOD))
      + geom_boxplot(fill="lightgray")
      + stat_sum(alpha=0.7)
      
      + scale_size(breaks=1:2, range=c(3,6))
)

surv_final_nooutlier %>% wilcox_test(MOD ~ Genotype) 



set.seed(101) ## for reproducibility
nsim <- 100000
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(surv_perm))
  bdat <- transform(surv_perm,MOD=MOD[perm])
  ## compute & store difference in means; store the value
  res[i] <- mean(bdat$MOD[bdat$Genotype=="1"])-
    mean(bdat$MOD[bdat$Genotype=="0"])
}
obs <- mean(surv_perm$MOD[surv_perm$Genotype=="1"])-
  mean(surv_perm$MOD[surv_perm$Genotype=="0"])
## append the observed value to the list of results
res <- c(res,obs)

hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")

permuted_p <- 2*mean(res>=obs)    

mean<- aggregate(surv_final_nooutlier[, c("MOD")], list(surv_final_nooutlier$Genotype), mean)
sd <- sd(surv_final_nooutlier$MOD)
d_cohen <- (mean$x[2] - mean$x[1]) / sd

p3<-surv_final_nooutlier %>%
  ggplot(aes(fill=as.factor(Genotype), y=MOD, x=as.factor(Genotype))) +
  geom_violin(position="dodge", alpha=0.5,) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(0,250)
p3

ggsave(p3, filename = paste(save.dir, 'violin_plot_survival_CREB3R119G.pdf', sep='/'), width = 3, height = 3)

table(surv_final_nooutlier$Genotype)






surv_control = surv[surv$Genotype == "0",]
surv_als = surv[surv$Genotype == "1",]

surv_control_list <- lapply(1:100, function(x) surv_control[sample(nrow(surv_control), size = 200, replace = TRUE),])
files <- paste("list", seq(1:100), sep = "_")
names(surv_control_list) <- files


library(broom)
dge_array <- lapply(files, function(list) {
  message('*******************************************************************************')
  message(paste(Sys.time(), 'processing', list))
  
  surv = rbind(surv_control_list[[list]],surv_als)
  surv$time = surv$MOD
  km_trt_fit <- survfit(Surv(time, status) ~ Genotype, data=surv) 
  summary(survfit(Surv(time, status) ~ Genotype, data = surv), times = 40)
  surv_diff <- survdiff(Surv(time, status) ~ Genotype, data = surv)
  s <- glance(surv_diff)
})

perm <- as.data.frame(do.call(rbind, dge_array))
perm <- -log10(perm$p.value)
sum(perm<1.3) / 100
hist(perm,col="gray",las=1,main="",breaks = 20)
abline(v=1.3,col="red")


