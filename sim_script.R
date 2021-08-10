#left entorhinal ST24CV
#left inferior temporal ST32CV
#left middle  temporal ST40CV
#left fusiform         ST26CV

#right entorhinal ST83CV
#right inferior temporal ST91CV
#right middle  temporal ST99CV
#right fusiform         ST85CV



#LeftLateralVentricle (ST37); 
#RightLateralVentricle (ST96); 
#LeftInferiorLateralVentricle (ST30); 
#RightInferiorLateralVentricle (ST89); 
#ThirdVentricle (ST127); 
#FourthVentricle (ST9); 
#FifthVentrical (ST8); 
#LeftChoroidPlexus (ST21);
#RightChoroidPlexus (ST80)

#CSFABETA cutpoint <= 254
#CSFPTAU cutpoint >= 24
#PETAMY cutpooint >= .78
install.packages("simr")

vol.ims <- c("ST103CV", "ST44CV", "ST29SV", "ST88SV", "ST10CV", 
             "ST24CV", "ST32CV", "ST40CV", 
             "ST26CV", "ST83CV", "ST91CV", 
             "ST99CV", "ST85CV", "ST37SV", "ST96SV",
             "ST30SV", "ST89SV", "ST127SV", "ST9SV", "ST21SV", "ST80SV")

source.script <- FALSE

#library(simstudy)
library(ADNIMERGE)
library(plyr)
library(dplyr)
library(furniture)
library(lme4)
library(nlme)
library(simr)
#library(MargCond)
library(stringr)

test.adas       <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADASSCORES.csv")
adni_imaging    <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/unharmonized_freesurfer_imaging.csv")
adas_scores_1   <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADASSCORES.csv")
adas_scores_23  <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADAS_ADNIGO23.csv")
adni_neuropsych <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Neuropsychological (1)/UWNPSYCHSUM_03_09_21.csv")
csf.upenn9      <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK9_04_19_17.csv")
csf.upenn10     <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK10_07_29_19.csv")
csf.upenn12     <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK12_01_04_21.csv")
neuropath.data  <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/NEUROPATH_05_17_21.csv")

adas_scores_1   <- adas_scores_1[,c("RID", "VISCODE",  "TOTAL11")]
adas_scores_23  <- adas_scores_23[,c("RID", "VISCODE2", "TOTSCORE")]
adni_neuropsych <- adni_neuropsych[,c("RID", "VISCODE2", "ADNI_MEM", "ADNI_EF")]
adni_imaging    <- adni_imaging[,c("RID", "VISCODE2", vol.ims)]
csf.upenn9      <- csf.upenn9[,c("RID", "VISCODE2", "EXAMDATE",  "ABETA", "TAU", "PTAU", "COMMENT")]
csf.upenn10      <- csf.upenn10[,c("RID", "VISCODE2", "DRAWDATE",  "ABETA42", "TAU", "PTAU", "NOTE")]
csf.upenn12      <- csf.upenn12[,c("RID", "VISCODE2", "EXAMDATE",  "ABETA", "TAU", "PTAU", "NOTE")]

colnames(adas_scores_23) <- c("RID", "VISCODE", "TOTAL11")
colnames(adni_neuropsych)<- c("RID", "VISCODE", "ADNI_MEM", "ADNI_EF")
colnames(adni_imaging)   <- c("RID", "VISCODE", vol.ims)
colnames(csf.upenn9) <- colnames(csf.upenn10) <- colnames(csf.upenn12) <- c("RID", "VISCODE", "EXAMDATE", "ABETA", "TAU", "PTAU", "COMMENT")
adas_scores              <- rbind(adas_scores_1, adas_scores_23)
csf.data                 <- rbind(csf.upenn9, csf.upenn10, csf.upenn12)
csf.data$EXAMDATE <- as.POSIXct(csf.data$EXAMDATE)
csf.data                 <- csf.data[order(csf.data$RID, csf.data$EXAMDATE, decreasing = FALSE), ]
rownames(csf.data)       <- 1:nrow(csf.data)
csf.drop <- which(duplicated(csf.data[,c("RID", "VISCODE")]))
image.drop <- which(duplicated(adni_imaging[,c("RID", "VISCODE")]))
adni_imaging <- adni_imaging[-image.drop,]
csf.data <- csf.data[-csf.drop,]
adas_scores              <- adas_scores[order(adas_scores$RID, 
                                              adas_scores$VISCODE, 
                                              decreasing = FALSE), ]
column.keeps             <- c("RID", "DX", "COLPROT", "ORIGPROT", 
                               "VISCODE", "EXAMDATE", "AGE", 
                               "PTGENDER", "PTEDUCAT", 
                               "APOE4", "CDRSB", "MMSE", 
                               "M", "mPACCtrailsB")
adas_demog       <- adnimerge[,column.keeps]
drops1           <- which(adas_scores$VISCODE == "")
drops2           <- which(adas_scores$VISCODE == "uns1")
adas_scores      <- adas_scores[-unique(union(drops1, drops2)), ]
adas_merge_demog <- merge(adas_demog, adas_scores, by=c("RID", "VISCODE"))
adas_merge_demog <- na.omit(adas_merge_demog)
adas_merge_demog <- merge(adas_merge_demog, adni_neuropsych, by=c("RID", "VISCODE"), all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, adni_imaging, by=c("RID", "VISCODE"), all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, csf.data, by=c("RID", "VISCODE"), all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, neuropath.data, by="RID", all = TRUE)
adas_merge_demog$EXAMDATE.x <- as.POSIXct(adas_merge_demog$EXAMDATE.x)
adas_merge_demog <- adas_merge_demog[order(adas_merge_demog$RID, adas_merge_demog$EXAMDATE.x, decreasing = FALSE), ]


av45                 <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/UCBERKELEYAV45_01_14_21.csv")
av45.keeps           <- c("RID", "EXAMDATE", "VISCODE2", "SUMMARYSUVR_COMPOSITE_REFNORM")
av45                 <- av45[,av45.keeps]
colnames(av45)       <- c("RID", "EXAMDATE", "VISCODE", "SUMMARYSUVR_COMPOSITE_REFNORM")
av45$AmyPos          <- rep(NA, nrow(av45))
adas_merge_demog$RID <- factor(adas_merge_demog$RID)
adas_merge_demog$EXAMDATE.x <- as.POSIXct(adas_merge_demog$EXAMDATE.x)
av45$EXAMDATE <- as.POSIXct(av45$EXAMDATE)
av45$RID             <- factor(av45$RID)
adas_merge_demog     <- merge(adas_merge_demog, av45, by=c("RID", "VISCODE"), all = TRUE)
names(adas_merge_demog)[names(adas_merge_demog) == 'EXAMDATE.x'] <- 'EXAMDATE_adnimerge'
names(adas_merge_demog)[names(adas_merge_demog) == 'EXAMDATE.y'] <- 'EXAMDATE_csf'
names(adas_merge_demog)[names(adas_merge_demog) == 'EXAMDATE'] <- 'EXAMDATE_pet'
adas_merge_demog$AGE <- round(adas_merge_demog$AGE + (adas_merge_demog$M / 12), 1)
adas_merge_demog$diff_adas_csf <- abs(difftime(adas_merge_demog$EXAMDATE_adnimerge, adas_merge_demog$EXAMDATE_csf, units = "weeks"))
adas_merge_demog$diff_adas_pet <- abs(difftime(adas_merge_demog$EXAMDATE_adnimerge, adas_merge_demog$EXAMDATE_pet, units = "weeks"))
adas_merge_demog$adas_csf_valid <- adas_merge_demog$adas_pet_valid <- rep(NA, nrow(adas_merge_demog))
adas_merge_demog["adas_csf_valid"][which(adas_merge_demog$diff_adas_csf/52 <= .5), ] <- 1
adas_merge_demog["adas_csf_valid"][which(adas_merge_demog$diff_adas_csf/52 > .5), ] <- 0
adas_merge_demog["adas_pet_valid"][which(adas_merge_demog$diff_adas_pet/52 <= .5), ] <- 1
adas_merge_demog["adas_pet_valid"][which(adas_merge_demog$diff_adas_pet/52 > .5), ] <- 0

adas.outcome.data <- adas_merge_demog[,c("RID", "VISCODE", "DX", "COLPROT", "ORIGPROT", 
                                         "VISCODE", "EXAMDATE_adnimerge", "AGE", 
                                         "PTGENDER", "PTEDUCAT", 
                                         "APOE4", "CDRSB", "MMSE", "TOTAL11", 
                                         "adas_pet_valid", "adas_csf_valid", "EXAMDATE_pet", "M",
                                         "ABETA", "TAU", "PTAU", "SUMMARYSUVR_COMPOSITE_REFNORM")]
rids.adas.outcome <- which(!is.na(adas.outcome.data$adas_pet_valid) | !is.na(adas.outcome.data$adas_csf_valid))
rids.adas.outcome <- unique(adas.outcome.data["RID"][rids.adas.outcome,])
adas.outcome.data <- subset(adas.outcome.data, RID %in% rids.adas.outcome)
adas.outcome.data <- adas.outcome.data[order(adas.outcome.data$RID, adas.outcome.data$EXAMDATE_adnimerge, decreasing = FALSE), ]
adas.drops <- which(is.na(adas.outcome.data$TOTAL11))
adas.outcome.data <- adas.outcome.data[-adas.drops, ]
adas.outcome.data$ABETA <- as.numeric(adas.outcome.data$ABETA)
adas.outcome.data$AmyPos <- rep(NA, nrow(adas.outcome.data))
rownames(adas.outcome.data) <- 1:nrow(adas.outcome.data)
adas.outcome.amypos <- which(adas.outcome.data$SUMMARYSUVR_COMPOSITE_REFNORM >= 0.78 | adas.outcome.data$ABETA <= 254)
adas.outcome.data["AmyPos"][adas.outcome.amypos, ] <- 1
adas.outcome.amyneg <- which(adas.outcome.data$SUMMARYSUVR_COMPOSITE_REFNORM < 0.78 & adas.outcome.data$ABETA > 254)
adas.outcome.data["AmyPos"][adas.outcome.amyneg, ] <- 0

adas.outcome.data$RID <- factor(adas.outcome.data$RID)
adas.outcome.data <- TimeSinceBaselineValidAmy(adas.outcome.data, "M")
adas.outcome.data <- subset(adas.outcome.data, new_time <= 24 & new_time != 18)
adas.outcome.data$RID <- factor(adas.outcome.data$RID)
adas.outcome.data <- CreateBaselineVar(adas.outcome.data, "AmyPos")
#Baseline cohort

adas_baseline           <- subset(adas.outcome.data, new_time == 0)
adas_baseline$DX        <- factor(adas_baseline$DX)
adas_baseline$APOE4     <- factor(adas_baseline$APOE4)
adas_baseline$PTGENDER  <- factor(adas_baseline$PTGENDER)
adas_baseline$COLPROT   <- factor(adas_baseline$COLPROT)
adas_baseline$AmyPos    <- factor(adas_baseline$AmyPos)
adas_baseline$new_time  <- factor(adas_baseline$new_time)

adas_desc <- adas_baseline[,c("DX", "APOE4", "PTGENDER", "AGE", "TOTAL11", 
                              "MMSE", "CDRSB", "PTEDUCAT", "SUMMARYSUVR_COMPOSITE_REFNORM",
                              "AmyPos")]
desc.table <- table1(adas_desc, splitby = ~ DX)$Table1
desc.table <- desc.table[c(1, 6:nrow(desc.table)), ]

adas_desc_time <- adas.outcome.data[,c("TOTAL11", "new_time", "DX")]
adas_desc_time$new_time <- factor(adas_desc_time$new_time)
desc.table.time <- table1(adas_desc_time, splitby = ~ new_time)$Table1
desc.table.time <- desc.table.time[c(1:3, 9:12), ]

#### model fitting


base.mci1     <- subset(adas_baseline, DX=="MCI" & CDRSB>=.5 & MMSE >= 20 & MMSE <=26 & AGE>=50 & AGE <=85 & AmyPos==1)
desc.table.base.mci <- base.mci1[,c("DX", "APOE4", "PTGENDER", "AGE", "TOTAL11", 
                                      "MMSE", "CDRSB", "PTEDUCAT",
                                      "AmyPos")]
desc.table.base.mci <- table1(desc.table.base.mci)$Table1

data.mci1     <- PullLongData(base.mci1, adas.outcome.data)



data.mci1$new_time <- as.numeric(data.mci1$new_time)
data.mci1$APOE4 <- factor(data.mci1$APOE4)
all.rids <- levels(data.mci1[["RID"]])
treat.rids <- sample(all.rids, 90, prob = rep(.5, length(all.rids)), replace = FALSE)
base.rids <- subset(all.rids, all.rids %notin% treat.rids)
data.mci1$treat <- rep(NA, nrow(data.mci1))
treat.rows <- which(data.mci1$RID %in% treat.rids)
base.rows <- which(data.mci1$RID %in% base.rids)
data.mci1["treat"][treat.rows,] <- 1
data.mci1["treat"][base.rows,] <- 0
data.mci1$treat <- factor(data.mci1$treat)
mod2 <- lmer(TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID), data = data.mci1)
mod2_large <- mod2
summary(mod2_large)
fixef(mod2_large)["treat1"] <- 0
fixef(mod2_large)["treat1:new_time"] <- (0.136017*.5)
summary(mod2_large)
mod2_ext_class <- extend(mod2_large, along="RID", n=650)
p_curve_treat_mod2 <- powerCurve(mod2_ext_class, test=fcompare(TOTAL11~new_time), along="RID", breaks = c(100, 200, 300, 400, 500, 600, 650))
sim_time           <-  powerSim(mod2_large, fcompare(TOTAL11 ~ new_time), nsim=100)

summary(sim_time)

plot(p_curve_treat_mod2)
this<- summary(p_curve_treat_mod2)
this$mean <- this$mean*100
this$lower <- this$lower*100
this$upper <- this$upper*100
mci.enrich.plot.1 <- ggplot(data = this, aes(x=nlevels, y=mean)) +geom_errorbar(aes(ymin=lower, ymax=upper)) + geom_line() + geom_point() + scale_x_discrete(limits=this$nlevels)
mci.enrich.plot.1 <- mci.enrich.plot.1 + xlab("Sample Size per Arm") + ylab("Statistical Power (%)")
mci.enrich.plot.1
all.rids <- levels(adas.outcome.data$RID)
nrids    <- nlevels(adas.outcome.data$RID)
keeps.rids  <- sample(all.rids, 300, prob = rep(.5, length(all.rids)))
unenriched.data<- subset(adas.outcome.data, RID %in% keeps.rids)
unenriched.data$RID <- factor(unenriched.data$RID)

unenriched.totall1<- SampleSizeSimulation(sim.data = unenriched.data, formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                          fcompare_str = "TOTAL11~new_time", breaks = seq(100, 1000, by=100))
data.mci1.totall11 <- SampleSizeSimulation(sim.data = data.mci1, formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                          fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100))
mci1.plot.total11 <- CombineSimPlots(nonenrich = unenriched.totall1, enrich = data.mci1.totall11, limits = seq(100, 1000, by=100))

unenriched.cdr <- SampleSizeSimulation(sim.data = unenriched.data, formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                          fcompare_str = "CDRSB~new_time", breaks = seq(100, 800, by=100))
data.mci1.cdr <- SampleSizeSimulation(sim.data = data.mci1, formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                          fcompare_str = "CDRSB~new_time", breaks = seq(100, 800, by=100))
mci1.plot.cdr <- CombineSimPlots(nonenrich = unenriched.cdr, enrich = data.mci1.cdr, limits = seq(100, 800, by=100))

firstplot <- mci1.plot.total11$plot
firstplot <- firstplot + labs(title = "ADAS-COG11") +geom_hline(yintercept = 80, linetype="dashed")
firstplot

secondplot <- mci1.plot.cdr$plot
secondplot <- secondplot + labs(title = "CDRSB") +geom_hline(yintercept = 80, linetype="dashed")
secondplot





amci1.analyzed <- PlotObsData(data.mci1, formula.fixed = "TOTAL11 ~ new_time", ylab= "ADAS-COG11")
mci1.analyzed
mci1.analyzed.cdr<- PlotObsData(data.mci1, formula.fixed = "CDRSB ~ new_time", ylab= "CDRSB")
mci1.analyzed.cdr


base.mci2     <- subset(adas_baseline, DX=="MCI" & CDRSB>=.5 & MMSE >= 20 & MMSE <=26 & AGE>=50 & AGE <=65 & AmyPos==1)
data.mci2     <- PullLongData(base.mci2, adas.outcome.data)
mci2.analyzed <- PlotObsData(data.mci2, formula.fixed = "TOTAL11 ~ new_time", ylab= "ADAS-COG11")
mci2.analyzed
mci2.analyzed.cdr <- PlotObsData(data.mci2, formula.fixed = "CDRSB ~ new_time", ylab= "CDRSB")
mci2.analyzed.cdr

base.ad1     <- subset(adas_baseline, DX == "Dementia" & CDRSB >= 1 & MMSE >= 20 & MMSE <= 26 & AGE >=50 & AGE <= 85) 
data.ad1     <- PullLongData(base.ad1, adas.outcome.data)
ad1.analyzed <- PlotObsData(data.ad1,formula.fixed = "TOTAL11 ~ new_time", ylab= "ADAS-COG11")
ad1.analyzed
ad1.analyzed.cdr <- PlotObsData(data.ad1,formula.fixed = "CDRSB ~ new_time", ylab= "CDRSB")
ad1.analyzed.cdr

base.ad2     <- subset(adas_baseline, DX == "Dementia" & CDRSB >= 1 & MMSE >= 20 & MMSE <= 26 & AGE >=50 & AGE <= 65) 
data.ad2     <- PullLongData(base.ad2, adas.outcome.data)
ad2.analyzed <- PlotObsData(data.ad2,formula.fixed = "TOTAL11 ~ new_time", ylab= "ADAS-COG11")
ad2.analyzed
ad2.analyzed.cdr <- PlotObsData(data.ad2,formula.fixed = "CDRSB ~ new_time", ylab= "CDRSB")
ad2.analyzed.cdr

intervals(mod1)

