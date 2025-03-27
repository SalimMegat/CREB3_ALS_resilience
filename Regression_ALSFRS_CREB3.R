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
setwd('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_progression/clinic_timeseries/')

############################################################################
############################################################################
######## Plot progression for rare variant and C9 and CREB3 carriers #######
############################################################################
############################################################################

alsfrs_prog = read.csv('~/Desktop/iMac_U1118/Desktop/MGoP/full_data_anchorALSFRS.csv', header = T)
alsfrs_prog = setDT(alsfrs_prog)
colnames(alsfrs_prog) <- c("row", "ID", "time", "alsfrst")
alsfrs_prog = alsfrs_prog[alsfrs_prog$time >=0,]
alsfrs_prog$time = alsfrs_prog$time*12
alsfrs_prog  <- separate(alsfrs_prog,ID, into = c("IID", "dataset"), sep="_(?=[^_]+$)")
alsfrs_prog$IID = gsub("CASE-", "CASE_", alsfrs_prog$IID)
#alsfrs_prog = alsfrs_prog[alsfrs_prog$time>0,]


# ### Load genotype  C9orf72 carriers #####
# genotype = fread('../geno_rare/expansion_C9orf72_all')
# c9_sample <- genotype[genotype$Genotype == "1",1]
# 
# c9orf72_alsfrsVisit = merge(alsfrs_prog, genotype, by=c("IID"))
# c9orf72_alsfrsVisit = c9orf72_alsfrsVisit[c9orf72_alsfrsVisit$time < 50,]
# tab=table(c9orf72_alsfrsVisit$IID)
# c9orf72_alsfrsVisit = c9orf72_alsfrsVisit[ifelse(tab[c9orf72_alsfrsVisit$IID]==1, FALSE, TRUE),] 
# 
# #pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/C9expansion_progression.pdf")
# ggplot(c9orf72_alsfrsVisit, aes(x = time, y = alsfrst, group = Genotype,color = Genotype )) +
# geom_smooth(method = "lm", level = 0.58, aes(fill=Genotype))
# #dev.off()
# 
# model <- lm(alsfrst ~ Genotype*time, data = c9orf72_alsfrsVisit)
# summary(model)
# 
# res_c9<- c9orf72_alsfrsVisit %>% 
#   group_by(IID) %>%
#   dplyr::summarise(across(starts_with("alsfrst"),
#                    list(slope = ~lm(. ~ time)$coef[2])))
# 
# # 
# # res_c9 <- c9orf72_alsfrsVisit[, 
# #                             {
# #                               ux <- mean(time)
# #                               uy <- mean(alsfrst)
# #                               slope <- sum((time - ux) * (alsfrst - uy)) / sum((time - ux) ^ 2)
# #                               list(slope=slope, intercept=uy - slope * ux)
# #                             }, by=IID]
# 
# 
# res_c9 = merge(c9orf72_alsfrsVisit,res_c9,by=c("IID"))
# res_c9 = res_c9[!duplicated(res_c9$IID),]
# res_c9$genotype = "C9"
### Load genotype  CREB3 carriers #####

genotype = fread('../geno_common/all_cohort_creb3_var.raw')

creb3_alsfrsVisit = merge(alsfrs_prog, genotype, by=c("IID"))
#creb3_alsfrsVisit = creb3_alsfrsVisit[creb3_alsfrsVisit$time < 30,]
creb3_alsfrsVisit = creb3_alsfrsVisit[creb3_alsfrsVisit$IID != "c36-2729-2729",]
tab=table(creb3_alsfrsVisit$IID)
creb3_alsfrsVisit = creb3_alsfrsVisit[ifelse(tab[creb3_alsfrsVisit$IID]==1, FALSE, TRUE),] 
#dev.off()
model <- lm(alsfrst ~ Genotype*time, data = creb3_alsfrsVisit)
summary(model)

res_creb3<- creb3_alsfrsVisit %>% 
  group_by(IID) %>%
  dplyr::summarise(across(starts_with("alsfrst"),
                          list(slope = ~lm(. ~ time)$coef[2])))

# res_creb3 <- creb3_alsfrsVisit[, 
#                               {
#                                 ux <- mean(time)
#                                 uy <- mean(alsfrst)
#                                 slope <- sum((time - ux) * (alsfrst - uy)) / sum((time - ux) ^ 2)
#                                 list(slope=slope, intercept=uy - slope * ux)
#                               }, by=IID
# ]

res_creb3 = merge(creb3_alsfrsVisit,res_creb3,by=c("IID"))
res_creb3 = res_creb3[!duplicated(res_creb3$IID),]
res_creb3$genotype = "CREB3"
res_creb3$exp = exp(res_creb3$alsfrst_slope)
#all_slope <- all_slope[all_slope$IID != "RES03916",]

#itals <- d50[d50$dataset == "itals",]
q3 <- quantile(res_creb3$exp, 0.75)
iqr <- IQR(res_creb3$exp)
upper_bound <- q3 + 1.5*iqr
q1 <- quantile(res_creb3$exp, 0.25)
lower_bound <- q1 - 1.5*iqr

res_creb3 <- res_creb3[res_creb3$exp<1.3,]

c9_sample = fread('../geno_rare/expansion_C9orf72_all')
c9_sample <- c9_sample[c9_sample$Genotype == "1",1]
id_to_removed <- read.csv('ID_itals_to_removed')
id_to_removed = id_to_removed$x

res_creb3 = res_creb3[!(res_creb3$IID %in% id_to_removed),]
res_creb3 = res_creb3[!(res_creb3$IID %in% c9_sample$IID),]


pwc <- res_creb3 %>%
  wilcox_test(exp ~ Genotype) 
pwc
aggregate(res_creb3$alsfrst_slope,list(res_creb3$Genotype), FUN=mean)
sd <- sd(res_creb3$alsfrst_slope)
d_cohen <- (-0.59 - (-0.80)) / sd


#pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_carriers.pdf", width = 10, height = 10)
time_course = creb3_alsfrsVisit[creb3_alsfrsVisit$time < 50,]
ggplot(time_course, aes(x = time, y = as.numeric(alsfrst), group = Genotype,color = Genotype )) +
  geom_smooth(method = "lm", level = 0.78, aes(fill=Genotype)) + ylim(20,50)

#pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_carriers.pdf")
# library(viridis)
# library(hrbrthemes)
# all_slope %>%
#   mutate(Genotype = factor(Genotype, levels = c("0", "1"))) %>%
#   mutate(genotype = factor(genotype, levels = c("C9", "CREB3"))) %>%
#   ggplot(aes(fill=Genotype, y=exp, x=genotype)) + 
#   geom_violin(position="dodge", alpha=0.5) +
#   geom_boxplot(width = 0.15,alpha=0.5,position = position_dodge(0.9)) +
#   scale_fill_viridis(discrete=T, name="") +
#   xlab("") +
#   ylab("ALSFRS Slope")
#dev.off()

set.seed(101) ## for reproducibility
nsim <- 9999
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(res_creb3))
  bdat <- transform(res_creb3,exp=exp[perm])
  ## compute & store difference in means; store the value
  res[i] <- mean(bdat$exp[bdat$Genotype=="1"])-
    mean(bdat$exp[bdat$Genotype=="0"])
}
obs <- mean(res_creb3$exp[res_creb3$Genotype=="1"])-
  mean(res_creb3$exp[res_creb3$Genotype=="0"])
## append the observed value to the list of results
res <- c(res,obs)

hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")

permuted_p <- 2*mean(res>=obs)    


library(viridis)
library(hrbrthemes)
min = min(res_creb3$exp) - 0.01
max = max(res_creb3$exp) + 0.01
pdf(file = "~/Desktop/single_cell/sc_als/processed_data/figures_raw/CREB3_carriers_eALSFRS.pdf")
res_creb3 %>%
  mutate(Genotype = factor(Genotype, levels = c("0", "1"))) %>%
  ggplot(aes(fill=Genotype, y=exp, x=Genotype)) + 
  geom_violin(position="dodge", alpha=0.5) +
  geom_boxplot(width = 0.15,alpha=0.5,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("ALSFRS Slope")
dev.off()


#### Load survival data and CREB3 genotype #####
library(tidyr)
library(multcomp)
library(rstatix)
library(ggpubr)
#### Process FTD GRN  ######
library(viridis)
library(hrbrthemes)
d50 = read.csv('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_progression/clinic_timeseries/all_d50_sept2024_final.csv', header = F)
colnames(d50) <- c("SubjectUID", "d50")
d50  <- d50 %>% 
  separate(SubjectUID, into = c("SubjectUID", "dataset"), sep="_(?=[^_]+$)")
d50$SubjectUID = gsub("CASE-", "CASE_", d50$SubjectUID)

surv_orig = read.csv('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_progression/clinic_timeseries/all_survival_sept2024_final.csv', header = T)
colnames(surv_orig)[1] <- "SubjectUID"
surv <- merge(surv_orig,d50,by=c("SubjectUID"))
surv = na.omit(surv)

surv %>%
  group_by(dataset) %>%
  summarize(cor=cor(d50, MOD))

d50 %>%
  group_by(dataset) %>%
  get_summary_stats(d50, type = "mean_sd")

surv %>%
  group_by(dataset) %>%
  get_summary_stats(MOD, type = "mean_sd")

res.aov2 <- aov(d50 ~ dataset, data = d50)
summary(res.aov2)
res <- TukeyHSD(res.aov2)
res <- res$dataset

surv %>%
  group_by(dataset) %>%
  get_summary_stats(MOD, type = "mean_sd")
df <- d50 %>% group_by(dataset) %>% dplyr::mutate(N=n())
p2<-d50 %>%
  ggplot(aes(fill=dataset, y=d50, x=dataset)) +
  geom_violin(position="dodge", alpha=0.5,) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(0,30)


q3 <- quantile(d50$d50, 0.75)
iqr <- IQR(d50$d50)
upper_bound <- q3 + 1.5*iqr
q1 <- quantile(d50$d50, 0.25)
lower_bound <- q1 - 1.5*iqr

tmp <- d50[d50$d50<8,]
itals_to_removed <- d50[d50$d50 > 8,1]


#d50 <- d50[d50$d50<10.8,]

p2<-tmp %>%
  ggplot(aes(fill=dataset, y=d50, x=dataset)) +
  geom_violin(position="dodge", alpha=0.5,) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(0,15)

tmp %>%
  group_by(dataset) %>%
  get_summary_stats(d50, type = "mean_sd")


#write.csv(itals$SubjectUID,"../clinic_timeseries/ID_itals.csv",quote = F, row.names = F)
#write.csv(itals_to_removed,"../clinic_timeseries/ID_itals_to_removed",quote = F, row.names = F)
#### Load survival data and CREB3 genotype #####
library(data.table)
id_to_removed <- read.csv('ID_itals_to_removed')
id_to_removed = id_to_removed$x
genotype = fread('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_progression/geno_common/all_cohort_creb3_r119g.raw')
#surv = read.csv('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_progression/als_italian_clinic/')
surv = fread('~/Desktop/iMac_U1118/Desktop/U1118_experiment/Article_progression/clinic_timeseries/all_survival_sept2024_final.csv')
surv = merge(genotype,surv, by=c("IID"))
#surv = surv[!(surv$IID %in% id_to_removed),]
surv = surv[surv$IID != "RES03091",]

sd <- 3*sd(surv$MOD)
surv = surv[surv$MOD<120,]

#aad <- read.table('covariates_ALS.PROG', header = T)

library(ggsurvfit)
library(broom)
surv$time = surv$MOD
surv$status = surv$status
km_trt_fit <- survfit(Surv(time, status) ~ Genotype, data=surv) 
summary(survfit(Surv(time, status) ~ Genotype, data = surv), times = 40)
surv_diff <- survdiff(Surv(time, status) ~ Genotype, data = surv)
s <- glance(surv_diff)

m <- coxph(Surv(time, status) ~ Genotype, data = surv) 


survfit2(Surv(time, status) ~ Genotype, data = surv) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) 
add_confidence_interval()
mean <- aggregate(surv[, c("MOD")], list(surv$Genotype), mean)
sd <- sd(surv$MOD)
d_cohen <- (55.4 - 37.7) / sd

obs <- -log10(s$p.value)
#surv = surv %>% mutate(points_bin = cut(MOD, breaks=c(0, 100, 200, 300)))
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



##### Compare mean survival with permutation test #######
pwc <- surv %>%
wilcox_test(MOD ~ Genotype, p.adjust.method = "BH") 


p3<-surv %>%
  ggplot(aes(fill=as.factor(Genotype), y=MOD, x=as.factor(Genotype))) +
  geom_violin(position="dodge", alpha=0.5,) +
  geom_boxplot(width = 0.15,position = position_dodge(0.9)) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("Eigengenes") +
  ylim(0,130)
p3

library("ggplot2"); theme_set(theme_bw())
library("lmPerm")
library("gtools")

surv_perm <- surv
surv_perm$Genotype = as.factor(surv_perm$Genotype)
print(ggplot(surv_perm,aes(Genotype,MOD))
      + geom_boxplot(fill="lightgray")
      + stat_sum(alpha=0.7)
      
      + scale_size(breaks=1:2, range=c(3,6))
)

set.seed(101) ## for reproducibility
nsim <- 9999
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


##### Cox proportional hazard and survival ######



