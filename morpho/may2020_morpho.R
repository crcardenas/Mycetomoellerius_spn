##################################
# NMDS of fovouros sp n vs zeteki#
#                                #                            
# Cody Raul Cardenas May 2020    #
##################################
library(tidyverse)
library(readr)
library(vegan)
library(ggplot2)
library(ggord)
library(rcompanion)
library(ggpubr)
library(broom)
library(purrr)
library(psych)
#library(MASS) #overwrites dplyr select option, use "dplyr::" before functions; ex dplyr::select(df, -c(column_name))
#import data
d <- read_csv("all_measurments.csv") %>% drop_na() %>% filter(!duplicated(bc_num))
dw <- d %>% filter(caste != "gyne", caste != "male") #177 workers
dw_BR <- dw %>% filter(species !="fovouros") 
dw_LB <- dw %>% filter(species !="zeteki") 
dg <- d %>% filter(caste != "worker", caste != "male") #53 gynes
dg_BR <- dg %>% filter(species !="fovouros")
dg_LB <- dg %>% filter(species !="zeteki")
dm <- d %>% filter(caste != "gyne", caste != "worker") #43 males
dm_BR <- dm %>% filter(species !="fovouros")
dm_LB <- dm %>% filter(species !="zeteki")
############
#! Summary #
############
# worker
# BR
wBRmorph_tab <- describe(dw_BR, na.rm=T,omit=T,quant=c(.05,.95)) %>% rownames_to_column() %>% 
  mutate(range=paste(Q0.25,Q0.75,sep="-")) %>% column_to_rownames() %>% 
  select(mean,se,range)
# LB
wLBmorph_tab <- describe(dw_LB, na.rm=T,omit=T,quant=c(.05,.95)) %>% rownames_to_column() %>% 
  mutate(range=paste(Q0.25,Q0.75,sep="-")) %>% column_to_rownames() %>% 
  select(mean,se,range)

# gyne
# BR
gBRmorph_tab <- describe(dg_BR, na.rm=T,omit=T,quant=c(.25,.75)) %>% rownames_to_column() %>% 
  mutate(range=paste(Q0.25,Q0.75,sep="-")) %>% column_to_rownames() %>% 
  select(mean,se,range)
# LB
gLBmorph_tab <- describe(dg_LB, na.rm=T,omit=T,quant=c(.25,.75)) %>% rownames_to_column() %>% 
  mutate(range=paste(Q0.25,Q0.75,sep="-")) %>% column_to_rownames() %>% 
  select(mean,se,range)

# male
# BR
mBRmorph_tab <- describe(dm_BR, na.rm=T,omit=T,quant=c(.25,.75)) %>% rownames_to_column() %>% 
  mutate(range=paste(Q0.25,Q0.75,sep="-")) %>% column_to_rownames() %>% 
  select(mean,se,range)
# LB
mLBmorph_tab <- describe(dm_LB, na.rm=T,omit=T,quant=c(.25,.75)) %>% rownames_to_column() %>% 
  mutate(range=paste(Q0.25,Q0.75,sep="-")) %>% column_to_rownames() %>% 
  select(mean,se,range)
#write files
write.table(wBRmorph_tab, file="morph_Mzet_workers.csv", sep=",", row.names=T, col.names = T)
write.table(wLBmorph_tab, file="morph_Mfov_workers.csv", sep=",", row.names=T, col.names = T)
write.table(gBRmorph_tab, file="morph_Mzet_gynes.csv", sep=",", row.names=T, col.names = T)
write.table(gLBmorph_tab, file="morph_Mfov_gynes.csv", sep=",", row.names=T, col.names = T)
write.table(mBRmorph_tab, file="morph_Mzet_males.csv", sep=",", row.names=T, col.names = T)
write.table(mLBmorph_tab, file="morph_Mfov_males.csv", sep=",", row.names=T, col.names = T)

##################################
#! normality test (shapiro wilks)#
##################################
# workers
shapw_w <- sapply(dw[,3:17],shapiro.test,simplify=T)
shapw_wBR <- sapply(dw_BR[,3:17],shapiro.test,simplify=T)
shapw_wLB <- sapply(dw_LB[,3:17],shapiro.test,simplify=T)
print(shapw_w)
# All samples; Not normal (bimodal, generate plots)
print(shapw_wBR)
print(shapw_wLB)
# normal BR: HW, HL, SL, FL, ML, PL, PPL, GL, CI, FLI, WL,TL 
# normal LB: HW, ML, GL, TL
# worker normal camparisons: HW, ML, GL, TL
# workers nonnormal compariosons: HL, SL, EL, FL, PL, PPL, CI, EI, SI, FLI, WL

#gynes
shapw_g <- sapply(dg[,3:17],shapiro.test,simplify=T)
shapw_gBR <- sapply(dg_BR[,3:17],shapiro.test,simplify=T)
shapw_gLB <- sapply(dg_LB[,3:17],shapiro.test,simplify=T)
print(shapw_g)
# CI, EI, FLI, are normal when comp all samples
print(shapw_gBR)
print(shapw_gLB)
# normal BR: PL, CI, SI, FLI, WL
# normal LB: PL, PPL, EI, WL, TL
# gyne normal camparisons: PL, WL
# gyne nonnormal compariosons: HW, HL, SL, EL, FL, ML, PPL, GL, CI, EI, SI, FLI, TL

#males
shapw_m <- sapply(dm[,3:17],shapiro.test,simplify=T)
shapw_mBR <- sapply(dm_BR[,3:17],shapiro.test,simplify=T)
shapw_mLB <- sapply(dm_LB[,3:17],shapiro.test,simplify=T)
print(shapw_m)
# GL, CI, EI, SI, FLI, TL, are all normally distributed for both
print(shapw_mBR)
print(shapw_mLB)
# normal BR: HW, HL, SL, ML, PL, SI, FLI, WL, TL
# normal LB: SL, ML, GL, CI, EI, SI, WL, TL
# male normal camparisons:  SL, ML, SI, WL, TL
# male nonnormal compariosons: HW, HL, EL, FL, PL, PPL, GL, CI, EI, FLI, 


#####################
#!test of difference#
#####################
# for normal data use Welch T's
# worker normal camparisons: ML, GL
# gyne normal camparisons: PL, FLI, WL
# male normal camparisons:  SL, ML, SI, WL, TL

welch_wHW <- t.test(HW~species,dw)                        
welch_wML <- t.test(dw$ML~dw$species) 
welch_wGL <- t.test(dw$GL~dw$species)
welch_wTL <- t.test(TL~species,dw)
welch_gPL <- t.test(PL~species,dg)
welch_gWL <- t.test(WL~species,dg)
welch_mSL <- t.test(dm$SL~dm$species)
welch_mML<- t.test(dm$ML~dm$species)
welch_mSI <- t.test(SI~species,dw)
welch_mWL <- t.test(dm$WL~dm$species)
welch_mTL <- t.test(dm$TL~dm$species)

welch_tab <- purrr::map_df(list(welch_wHW, welch_wML, welch_wGL, welch_wTL,
                                welch_gPL, welch_gWL,
                                welch_mSL, welch_mML, welch_mSI, welch_mWL, welch_mTL), broom::tidy) %>% 
  mutate("fovouros mean"=estimate1,"zeteki mean"=estimate2,"t-statistic"=statistic,DF=parameter,
         "confidence intverval low" = conf.low, "confidence interval high" = conf.high) %>%
  select("fovouros mean","zeteki mean", "t-statistic",DF,p.value) %>% 
  as.matrix() 
rownames(welch_tab) <- c("worker HW", "worker ML", "worker GL", "worker TL", 
                         "gyne PL", "gyne WL",
                         "male SL", "male ML", "male SI", "male WL", "male TL")

# For non-normal use Wilcoxon Rank-Sum test (requires that the population distribution of 
# differneces be symmetric about the unknown median)
# workers

# workers nonnormal compariosons: HL, SL, EL, FL, PL, PPL, CI, EI, SI, FLI, WL, TL
cox_wHL <- wilcox.test(dw_BR$HL, dw_LB$HL, conf.int = T, exact=F)
cox_wSL <- wilcox.test(dw_BR$SL, dw_LB$SL, conf.int = T, exact=F)
cox_wEL <- wilcox.test(dw_BR$EL, dw_LB$EL, conf.int = T, exact=F)
cox_wFL <- wilcox.test(dw_BR$FL, dw_LB$FL, conf.int = T, exact=F)
cox_wPL <- wilcox.test(dw_BR$PL, dw_LB$PL, conf.int = T, exact=F)
cox_wPPL <- wilcox.test(dw_BR$PPL, dw_LB$PPL, conf.int = T, exact=F)
cox_wCI <- wilcox.test(dw_BR$CI, dw_LB$CI, conf.int = T, exact=F)
cox_wEI <- wilcox.test(dw_BR$EI, dw_LB$EI, conf.int = T, exact=F)
cox_wSI <- wilcox.test(dw_BR$SI, dw_LB$SI, conf.int = T, exact=F)
cox_wFLI <- wilcox.test(dw_BR$FLI, dw_LB$FLI, conf.int = T, exact=F)
cox_wWL <- wilcox.test(dw_BR$WL, dw_LB$WL, conf.int = T, exact=F)

# gynes
# gyne nonnormal compariosons: HW, HL, SL, EL, FL, ML, PPL, GL, CI, EI, SI, FLI, TL
cox_gHW <- wilcox.test(dg_BR$HW, dg_LB$HW, conf.int = T, exact=F) 
cox_gHL <- wilcox.test(dg_BR$HL, dg_LB$HL, conf.int = T, exact=F) 
cox_gSL <- wilcox.test(dg_BR$SL, dg_LB$SL, conf.int = T, exact=F)  
cox_gEL <- wilcox.test(dg_BR$EL, dg_LB$EL, conf.int = T, exact=F)   
cox_gFL <- wilcox.test(dg_BR$FL, dg_LB$FL, conf.int = T, exact=F) 
cox_gML <- wilcox.test(dg_BR$ML, dg_LB$ML, conf.int = T, exact=F) 
cox_gPPL <- wilcox.test(dg_BR$PPL, dg_LB$PPL, conf.int = T, exact=F) 
cox_gGL <- wilcox.test(dg_BR$GL, dg_LB$GL, conf.int = T, exact=F) 
cox_gCI <- wilcox.test(dg_BR$CI, dg_LB$CI, conf.int = T, exact=F) 
cox_gEI <- wilcox.test(dg_BR$EI, dg_LB$EI, conf.int = T, exact=F) 
cox_gSI <- wilcox.test(dg_BR$SI, dg_LB$SI, conf.int = T, exact=F) 
cox_gFLI <- wilcox.test(dg_BR$FLI, dg_LB$FLI, conf.int = T, exact = F) 
cox_gTL <- wilcox.test(dg_BR$TL, dg_LB$TL, conf.int = T, exact=F) 

# males
# male nonnormal compariosons: HW, HL, EL, FL, PL, PPL, GL, CI, EI, FLI, 
cox_mHW <- wilcox.test(dm_BR$HW, dm_LB$HW, conf.int = T, exact=F) 
cox_mHL <- wilcox.test(dm_BR$HL, dm_LB$HL, conf.int = T, exact=F) 
cox_mEL <- wilcox.test(dm_BR$EL, dm_LB$EL, conf.int = T, exact=F) 
cox_mFL <- wilcox.test(dm_BR$FL, dm_LB$FL, conf.int = T, exact=F) 
cox_mPL <- wilcox.test(dm_BR$PL, dm_LB$PL, conf.int = T, exact=F) 
cox_mPPL <- wilcox.test(dm_BR$PPL, dm_LB$PPL, conf.int = T, exact=F) 
cox_mGL <- wilcox.test(dm_BR$GL, dm_LB$GL, conf.int = T, exact=F) 
cox_mCI <- wilcox.test(dm_BR$CI, dm_LB$CI, conf.int = T, exact=F) 
cox_mEI <- wilcox.test(dm_BR$EI, dm_LB$EI, conf.int = T, exact=F) 
cox_mFLI <- wilcox.test(dm_BR$FLI, dm_LB$FLI, conf.int = T, exact=F) 

wilcox_tab <- purrr::map_df(list(cox_wHL,cox_wSL,cox_wEL,cox_wFL,cox_wPL,cox_wPPL,cox_wCI,
                                 cox_wEI,cox_wSI,cox_wFLI,cox_wWL,
                                 cox_gHW,cox_gHL,cox_gSL,cox_gEL,cox_gFL,cox_gML,cox_gPPL,cox_gGL,
                                 cox_gCI,cox_gEI,cox_gSI,cox_gFLI,cox_gTL,
                                 cox_mHW,cox_mHL,cox_mEL,cox_mFL,cox_mPL,cox_mPPL,cox_mGL,
                                 cox_mCI,cox_mEI,cox_mFLI), 
                            broom::tidy) %>% 
  mutate("W-Stat"=statistic, "confidence low" = conf.low, "confidence high" = conf.high) %>% 
  select("W-Stat","confidence low","confidence high",p.value) %>% 
  as.matrix()

rownames(wilcox_tab) <- c("worker HL","worker SL", "worker EL", "worker FL", "worker PL", 
                          "worker PPL","worker CI", "worker EI", "worker SI", "worker FLI", "worker WL",
                          "gyne HW", "gyne HL", "gyne SL", "gyne EL", "gyne FL", "gyne ML", "gyne PPL",
                          "gyne GL", "gyne CI", "gyne EI", "gyne SI", "gyne FLI", "gyne TL",
                          "male HW", "male HL", "male EL", "male FL", "male PL", "male PPL",
                          "male GL", "male CI", "male EI", "male FLI")

###############
#!test results#
###############
#welch t-test 
print(welch_tab)

# Wilcoxon Rank Sum Test
print(wilcox_tab)
# cannot use exact = True because there are ties in the data (IE some data has the same value)
# thus the P-value is not an exact computation

# W is the number of ranks that can be compared.
# account for bonferroni correction to adjust P-value according to
# the number of tests performed for each partition of our data
# P value < 0.003? 
# All but Gyne FLI, & Male FLI are significant. This makes sense
# these data were  normially distributed when looking at the 
# shapiro-wilks test

write.table(welch_tab, file="welch_ttest.csv", sep=",", row.names=T, col.names = T)
write.table(wilcox_tab, file="wilcox_ranktest.csv", sep=",", row.names=T, col.names = T)

#remove samples
dg2 <- dg %>% select(!FLI) #rm FLI
dm2 <- dm %>% select(!FLI) #rm FLI

############
#!  NMDS   #
############
set.seed(666)
#calculate disimilarity
#workers
dw.dis <- vegan::vegdist(dw[,3:17])
dw.mds1 <- vegan::metaMDS(dw[,3:17],try=1000,k=2,distance = "bray")
dw.mds1 # Stress: 0.1291963  
#gynes
dg2.dis <- vegdist(dg2[,3:16])
dg2.mds1 <- metaMDS(dg2[,3:16],try=1000,k=2,distance = "bray")
dg2.mds1 # Stress: 0.1190165 
#males
dm2.dis <- vegdist(dm2[,3:16])
dm2.mds1 <- metaMDS(dm2[,3:16],try=1000,k=2,distance = "bray")
dm2.mds1 # Stress: 0.1554117  

# shepards plots
stressplot(dw.mds1,dw.dis, main="Worker NMDS stress plot")
stressplot(dg2.mds1,dg2.dis, main="Gyne NMDS stress plot")
stressplot(dm2.mds1,dm2.dis, main="Male NMDS stress plot")

#plots ### CHECK YOUR CHICKEN PLOTS YA DONGK
#get types & paratype names and points
wtype_para <- head(as.data.frame(scores(dw.mds1)),10)
row.names(wtype_para) <- head(dw$bc_num,10)
wtype_para$spp<-head(dw$species,10)
#get all other points and names
wpoints <- as.data.frame(scores(dw.mds1))
row.names(wpoints) <- dw$bc_num
wpoints$spp <-dw$species
#plot using ggord & ggplot functions
worker_NMDS <- ggord(dw.mds1, axes=c("1","2"), grp_in=dw$species, 
                     ellipse_pro = 0.95, cols=c("#E69F00","#CC79A7"), 
                     size=2) +
  annotate("text", x=0.01, y=-0.03, label="n=171 stress=0.1291963") +
  scale_shape_manual('Groups', values = c(16, 17)) +
  ggtitle("Worker NMDS")+
  theme(plot.title = element_text(hjust = 0.5))
worker_NMDS = worker_NMDS + geom_text(data=wtype_para, #labels for paratypes & types
                     aes(x=NMDS1, y=NMDS2,
                         label=row.names(wtype_para)), 
                     alpha=0.60)
worker_NMDS
  


gtype_para <- head(as.data.frame(scores(dg2.mds1)),8)
row.names(gtype_para) <- head(dg2$bc_num,8)
gtype_para$spp<-head(dg2$species,8)
#get all other points and names
gpoints <- as.data.frame(scores(dg2.mds1))
row.names(gpoints) <- dg2$bc_num
gpoints$spp <-dg2$species
#... uhhh

gyne_NMDS <- ggord(dg2.mds1, axes=c("1","2"), grp_in=dg2$species, 
                     ellipse_pro = 0.95, cols=c("#E69F00","#CC79A7"), 
                     size=3) +
  annotate("text", x=0.03, y=-0.0275, label="n=53 stress=0.1190165") +
  scale_shape_manual('Groups', values = c(16, 17)) +
  ggtitle("gyne NMDS")+
  theme(plot.title = element_text(hjust = 0.5))
gyne_NMDS = gyne_NMDS + geom_text(data=gtype_para, #labels for paratypes & types
                                      aes(x=NMDS1, y=NMDS2,
                                          label=row.names(gtype_para)), 
                                      alpha=0.60)
gyne_NMDS



mtype_para <- head(as.data.frame(scores(dm2.mds1)),7)
row.names(mtype_para) <- head(dm2$bc_num,7)
mtype_para$spp<-head(dm2$species,7)
#get all other points and names
mpoints <- as.data.frame(scores(dm2.mds1))
row.names(mpoints) <- dm2$bc_num
mpoints$spp <-dm2$species

male_NMDS <- ggord(dm2.mds1, axes=c("1","2"), grp_in=dm2$species, 
                   ellipse_pro = 0.95, cols=c("#E69F00","#CC79A7"), 
                   size=3) +
  annotate("text", x=0.01, y=-0.02, label="n=43 stress=0.1554117") +
  scale_shape_manual('Groups', values = c(16, 17)) +
  ggtitle("male NMDS")+
  theme(plot.title = element_text(hjust = 0.5))
male_NMDS = male_NMDS + geom_text(data=mtype_para, #labels for paratypes & types
                                  aes(x=NMDS1, y=NMDS2,
                                      label=row.names(mtype_para)), 
                                  alpha=0.60)
male_NMDS