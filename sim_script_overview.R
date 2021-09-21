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

#CSFABETA cutpoint <= 980
#CSFPTAU cutpoint  >= 24
#PETAMY cutpooint  >= .78



feature.correction <- function(training.data,  data, formula, cr.feat1, feat) {
  model <- lm(formula = as.formula(formula), data = training.data)
  mean.val1 <- mean(training.data[[cr.feat1]])
  mean.val1 <- rep(mean.val1, nrow(data))
  coef.correction1 <- model$coefficients[[cr.feat1]]
  new.feat <- data[[feat]]
  new.feat <- (new.feat - (coef.correction1*data[[cr.feat1]]))
  new.feat <- new.feat  + (coef.correction1*mean.val1)
  return(new.feat)
}



source.script <- FALSE
library(zoo)
#library(simstudy)
library(ADNIMERGE)
library(plyr)
library(dplyr)
library(furniture)
library(lme4)
library(nlme)
library(simr)
library(stringr)
library(matrixStats)
library(lmerTest)
library(splitstackshape)
library(purrr)
library(ggplot2)
#CDR
cdr_global           <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/CDR (1).csv")
cdr_global$EXAMDATE  <- as.POSIXct(cdr_global$EXAMDATE, format="%Y-%M-%D")
cdr_global           <- cdr_global[order(cdr_global$RID, cdr_global$EXAMDATE, decreasing = FALSE), ]
cdr_global$VISCODE2  <- as.character(cdr_global$VISCODE2)
cdr_global           <- cdr_global[,c("RID", "VISCODE2", "CDGLOBAL")]
colnames(cdr_global) <- c("RID", "VISCODE", "CDGLOBAL")
cdr_global["VISCODE"][which(cdr_global$VISCODE=="sc"),] <- "bl"

#imaging
# harmonized in harmonize_imaging_for_simulation

imaging_simulation_data <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/harmed_freesurfer_imaging.csv")



#ADAS
adas_scores_1   <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADASSCORES.csv")
adas_scores_23  <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADAS_ADNIGO23.csv")
adas_scores_1   <- adas_scores_1[,c("RID", "VISCODE",  "TOTAL11", "TOTALMOD")]
adas_scores_23  <- adas_scores_23[,c("RID", "VISCODE2", "TOTSCORE", "TOTAL13")]
colnames(adas_scores_1) <- colnames(adas_scores_23) <- c("RID", "VISCODE", "TOTAL11", "TOTAL13")
adas_scores              <- rbind(adas_scores_1, adas_scores_23)
adas_scores              <- adas_scores[order(adas_scores$RID, 
                                              adas_scores$VISCODE, 
                                              decreasing = FALSE), ]
drops1           <- which(adas_scores$VISCODE == "")
adas_scores      <- adas_scores[-drops1, ]


#Neuropsych
adni_neuropsych <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Neuropsychological (1)/UWNPSYCHSUM_03_09_21.csv")
adni_neuropsych <- adni_neuropsych[,c("RID", "VISCODE2", "ADNI_MEM", "ADNI_EF")]
colnames(adni_neuropsych)<- c("RID", "VISCODE", "ADNI_MEM", "ADNI_EF")


#CSF
csf.upenn9      <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK9_04_19_17.csv")
csf.upenn10     <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK10_07_29_19.csv")
csf.upenn12     <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK12_01_04_21.csv")
csf.upenn9      <- csf.upenn9[,c("RID", "VISCODE2", "EXAMDATE",  "ABETA", "TAU", "PTAU", "COMMENT")]
csf.upenn10      <- csf.upenn10[,c("RID", "VISCODE2", "DRAWDATE",  "ABETA42", "TAU", "PTAU", "NOTE")]
csf.upenn12      <- csf.upenn12[,c("RID", "VISCODE2", "EXAMDATE",  "ABETA", "TAU", "PTAU", "NOTE")]
csf.upenn9$ABETA <-           as.numeric(as.character(csf.upenn9$ABETA))
csf.upenn9$PTAU <-           as.numeric(as.character(csf.upenn9$PTAU))
csf.upenn10$ABETA42 <-           as.numeric(as.character(csf.upenn10$ABETA42))
csf.upenn10$PTAU <-           as.numeric(as.character(csf.upenn10$PTAU))
csf.upenn12$ABETA <-           as.numeric(as.character(csf.upenn12$ABETA))
csf.upenn12$PTAU <-           as.numeric(as.character(csf.upenn12$PTAU))
csf.upenn9$transform <- rep(TRUE, nrow(csf.upenn9))
csf.upenn10$transform <- rep(TRUE, nrow(csf.upenn10))
csf.upenn12$transform <- rep(FALSE, nrow(csf.upenn12))
colnames(adni_imaging)   <- c("RID", "VISCODE", vol.ims)
colnames(csf.upenn9) <- colnames(csf.upenn10) <- colnames(csf.upenn12) <- c("RID", "VISCODE", "EXAMDATE", "ABETA", "TAU", "PTAU", "COMMENT", "TRANSFORM")
csf.data                 <- rbind(csf.upenn9, csf.upenn10, csf.upenn12)
csf.data$EXAMDATE <- as.POSIXct(csf.data$EXAMDATE)
csf.data                 <- csf.data[order(csf.data$RID, csf.data$EXAMDATE, decreasing = FALSE), ]
csf.data$ABETA <-           as.numeric(as.character(csf.data$ABETA))
csf.data$PTAU <-           as.numeric(as.character(csf.data$PTAU))
csf.data["ABETA"][which(csf.data$TRANSFORM == TRUE), ] <- (csf.data["ABETA"][which(csf.data$TRANSFORM == TRUE), ]*1.014) + 29.25
csf.data["PTAU"][which(csf.data$TRANSFORM == TRUE), ] <- (csf.data["PTAU"][which(csf.data$TRANSFORM == TRUE), ]*.961) - .694
amyloid.cutpoint <- (980 *1.014) + 29.25
ptau.cutpoint  <-  (24*.961) - .694
rownames(csf.data)       <- 1:nrow(csf.data)
csf.drop <- which(duplicated(csf.data[,c("RID", "VISCODE")]))
csf.data <- csf.data[-csf.drop,]


#Neuropathology Data
neuropath.data  <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/NEUROPATH_05_17_21.csv")
non.ad.imputation <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADNI-nonADimputations-class.csv")
non.ad.imputation["TDP43"][which(non.ad.imputation$TDP43 == TRUE), ] <- 1
non.ad.imputation["TDP43"][which(non.ad.imputation$TDP43 == FALSE), ] <- 0
non.ad.imputation["LEWY"][which(non.ad.imputation$LEWY == TRUE), ] <- 1
non.ad.imputation["LEWY"][which(non.ad.imputation$LEWY == FALSE), ] <- 0
non.ad.imputation["CAA"][which(non.ad.imputation$CAA == TRUE), ] <- 1
non.ad.imputation["CAA"][which(non.ad.imputation$CAA == FALSE), ] <- 0



#PET Amyloid
av45                 <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/UCBERKELEYAV45_01_14_21.csv")
av45.keeps           <- c("RID", "EXAMDATE", "VISCODE2", "SUMMARYSUVR_COMPOSITE_REFNORM")
av45                 <- av45[,av45.keeps]
colnames(av45)       <- c("RID", "EXAMDATE", "VISCODE", "SUMMARYSUVR_COMPOSITE_REFNORM")


#ADNI Merge
column.keeps             <- c("RID", "DX", "COLPROT", "ORIGPROT", 
                               "VISCODE", "EXAMDATE", "AGE", 
                               "PTGENDER", "PTEDUCAT", 
                               "APOE4", "CDRSB", "MMSE", 
                               "M", "mPACCtrailsB")
adas_demog       <- adnimerge[,column.keeps]



#Merging Data
adas_merge_demog <- merge(adas_demog, adas_scores, by=c("RID", "VISCODE"), all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, adni_neuropsych, by=c("RID", "VISCODE"), all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, csf.data, by=c("RID", "VISCODE"), all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, imaging_simulation_data, by=c("RID", "VISCODE"), all = TRUE)


names(adas_merge_demog)[names(adas_merge_demog) == 'EXAMDATE.x'] <- 'EXAMDATE_adnimerge'
names(adas_merge_demog)[names(adas_merge_demog) == "EXAMDATE.y"] <- "EXAMDATE_csf"
adas_merge_demog <- merge(adas_merge_demog, neuropath.data, by="RID", all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, non.ad.imputation, by="RID", all = TRUE)
adas_merge_demog <- merge(adas_merge_demog, cdr_global, by=c("RID", "VISCODE"), all = TRUE)
adas_merge_demog            <- merge(adas_merge_demog, av45, by=c("RID", "VISCODE"), all = TRUE)
names(adas_merge_demog)[names(adas_merge_demog) == 'EXAMDATE'] <- 'EXAMDATE_pet'
adas_merge_demog$AGE <- round(adas_merge_demog$AGE + (adas_merge_demog$M / 12), 1)
adas_merge_demog$EXAMDATE_adnimerge <- as.POSIXct(adas_merge_demog$EXAMDATE_adnimerge)
adas_merge_demog$EXAMDATE_csf<- as.POSIXct(adas_merge_demog$EXAMDATE_csf)
adas_merge_demog$EXAMDATE_pet<- as.POSIXct(adas_merge_demog$EXAMDATE_pet)

adas_merge_demog$diff_csf <- abs(difftime(adas_merge_demog$EXAMDATE_adnimerge, adas_merge_demog$EXAMDATE_csf, units = "weeks"))
adas_merge_demog$diff_pet <- abs(difftime(adas_merge_demog$EXAMDATE_adnimerge, adas_merge_demog$EXAMDATE_pet, units = "weeks"))
adas_merge_demog$csf_valid <- adas_merge_demog$pet_valid <- rep(NA, nrow(adas_merge_demog))
adas_merge_demog["csf_valid"][which(adas_merge_demog$diff_csf/52 <= .5), ] <- 1
adas_merge_demog["csf_valid"][which(adas_merge_demog$diff_csf/52 > .5), ] <- 0
adas_merge_demog["pet_valid"][which(adas_merge_demog$diff_pet/52 <= .5), ] <- 1
adas_merge_demog["pet_valid"][which(adas_merge_demog$diff_pet/52 > .5), ] <- 0
adas_merge_demog$csf_pos <- adas_merge_demog$suvr_pos<- adas_merge_demog$ptau_pos<- rep(NA, nrow(adas_merge_demog))
adas_merge_demog$ABETA   <- as.numeric(as.character(adas_merge_demog$ABETA))
adas_merge_demog["suvr_pos"][which(adas_merge_demog$pet_valid == 1 & adas_merge_demog$SUMMARYSUVR_COMPOSITE_REFNORM >= .78), ] <- 1
adas_merge_demog["suvr_pos"][which(adas_merge_demog$pet_valid == 1 & adas_merge_demog$SUMMARYSUVR_COMPOSITE_REFNORM < .78), ] <- 0
adas_merge_demog["csf_pos"][which(adas_merge_demog$csf_valid == 1 & adas_merge_demog$ABETA < amyloid.cutpoint), ] <- 1
adas_merge_demog["csf_pos"][which(adas_merge_demog$csf_valid == 1 & adas_merge_demog$ABETA >= amyloid.cutpoint), ] <- 0
adas_merge_demog["ptau_pos"][which(adas_merge_demog$csf_valid == 1 & adas_merge_demog$PTAU >= ptau.cutpoint), ] <- 1
adas_merge_demog["ptau_pos"][which(adas_merge_demog$csf_valid == 1 & adas_merge_demog$PTAU < ptau.cutpoint), ] <- 0
adas_merge_demog <- adas_merge_demog[order(adas_merge_demog$RID, adas_merge_demog$VISCODE, decreasing = FALSE), ]
adas_merge_demog$AmyPos <- rep(NA, nrow(adas_merge_demog))
adas_merge_demog$AmyPos_full <- rep(NA, nrow(adas_merge_demog))
adas_merge_demog$TauPos_full <- rep(NA, nrow(adas_merge_demog))

for(i in 1:nrow(adas_merge_demog)) {
  if(!is.na(adas_merge_demog["suvr_pos"][i,]) & adas_merge_demog["suvr_pos"][i,]==1) {
    adas_merge_demog["AmyPos"][i,] <- 1
  } else if(!is.na(adas_merge_demog["csf_pos"][i,]) & adas_merge_demog["csf_pos"][i,]==1) {
    adas_merge_demog["AmyPos"][i,] <- 1
  } else if(!is.na(adas_merge_demog["csf_pos"][i,]) | !is.na(adas_merge_demog["suvr_pos"][i,])) {
    adas_merge_demog["AmyPos"][i,] <- 0
  }
}



adas_merge_demog$M_vis <- substr(adas_merge_demog$VISCODE, 2, 10)
adas_merge_demog$M_vis <- as.numeric(adas_merge_demog$M_vis)
adas_merge_demog["M_vis"][which(adas_merge_demog$VISCODE=="bl"), ] <- 0
adas_merge_demog <- adas_merge_demog[order(adas_merge_demog$RID, adas_merge_demog$M_vis, decreasing = FALSE), ]
adas_merge_demog$fulllewy <- adas_merge_demog$fulltdp43 <- adas_merge_demog$fullcaa <- rep(NA, nrow(adas_merge_demog)) 
baseline.var.list <- c("DX","AGE", "PTGENDER", "PTEDUCAT", 
                       "APOE4", "CDRSB", "MMSE", "mPACCtrailsB", "TOTAL11", "TOTAL13", 
                       "ABETA", "TAU", "PTAU", "AmyPos","ptau_pos", "CDGLOBAL", "ADNI_MEM", "ADNI_EF")
adas_merge_demog$DX <- as.character(adas_merge_demog$DX)
adas_merge_demog <- adas_merge_demog[order(adas_merge_demog$RID, adas_merge_demog$M_vis, decreasing = FALSE), ]
adas_merge_demog$RID <- factor(as.character(adas_merge_demog$RID))

for(i in baseline.var.list) {
  adas_merge_demog <- CreateBaselineVar(adas_merge_demog, "M_vis", i)
}
baseline.var.list <- paste(baseline.var.list, "_bl", sep="")



adas_merge_demog$image_remerging <- 1:nrow(adas_merge_demog)
image_columns <- grep("_harmonized", colnames(adas_merge_demog), value = TRUE)
image_data_for_adj <- adas_merge_demog[,c("RID", "DX_bl", "AGE_bl", "PTGENDER", "AmyPos_bl", image_columns, "image_remerging")]
image_data_for_adj <- na.omit(image_data_for_adj)
image_data_for_adj$Baseline <- !duplicated(image_data_for_adj$RID)
image_data_for_adj["x1"] <- 5
icv.training <- subset(image_data_for_adj, Baseline==TRUE & DX_bl=="CN" & AmyPos_bl==0)
image_columns <- image_columns[!image_columns=="ST10CV_harmonized"]


for(i in image_columns) {
  correction_formula <- paste(i, "~ AGE_bl + PTGENDER + ST10CV_harmonized", sep=" ")
  corrected_name <- paste(i, "icv_adj", sep="_")
  corrected_feature  <- feature.correction(training.data = icv.training,
                                           data = image_data_for_adj,
                                           formula = correction_formula,
                                           cr.feat1 = "ST10CV_harmonized",
                                           feat = i)
  image_data_for_adj[corrected_name] <- corrected_feature
}

image_columns.adj <- grep("_harmonized_icv_adj", colnames(image_data_for_adj), value = TRUE)
imaging_to_remerge <- image_data_for_adj[,c(image_columns.adj, "image_remerging")]
adas_merge_demog <- merge(adas_merge_demog, imaging_to_remerge, by="image_remerging", all.x = TRUE)















# All neuropathology data for checking effect size of neurpats on outcome 

write.csv(adas_merge_demog, "/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/final_merged_data_4sim.csv")


aut.rows <- which(complete.cases(adas_merge_demog[,c("NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                                    "NPLBOD", "NPAMY")]))



imputed.rows <- which(complete.cases(adas_merge_demog[,c("TDP43", "LEWY", "CAA")]))
keeps <- unique(union(aut.rows, imputed.rows))
all.neuropat <- adas_merge_demog[keeps,]
all.neuropat <- SetNeuroData(all.neuropat)






for(i in 1:nrow(all.neuropat)) {
  if(!is.na(all.neuropat["Amy_pos_path"][i,]) & all.neuropat["Amy_pos_path"][i,]==1) {
    all.neuropat["AmyPos_full"][i,] <- 1
  } else if(!is.na(all.neuropat["AmyPos"][i,]) & all.neuropat["AmyPos"][i,]==1) {
    all.neuropat["AmyPos_full"][i,] <- 1
  } else if(!is.na(all.neuropat["Amy_pos_path"][i,]) | !is.na(all.neuropat["AmyPos"][i,])) {
    all.neuropat["AmyPos_full"][i,] <- 0
  }
}

for(i in 1:nrow(all.neuropat)) {
  if(!is.na(all.neuropat["TAU_pos_path"][i,]) & all.neuropat["TAU_pos_path"][i,]==1) {
    all.neuropat["TauPos_full"][i,] <- 1
  } else if(!is.na(all.neuropat["ptau_pos"][i,]) & all.neuropat["ptau_pos"][i,]==1) {
    all.neuropat["TauPos_full"][i,] <- 1
  } else if(!is.na(all.neuropat["TAU_pos_path"][i,]) | !is.na(all.neuropat["ptau_pos"][i,])) {
    all.neuropat["TauPos_full"][i,] <- 0
  }
}

all.neuropat$RID <- factor(all.neuropat$RID)
all.neuropat <- CreateBaselineVar(all.neuropat, "M_vis", "AmyPos_full")
all.neuropat <- CreateBaselineVar(all.neuropat, "M_vis", "TauPos_full")
all.neuropat$imputed <- rep(NA, nrow(all.neuropat))





all.neuropat <- all.neuropat[complete.cases(all.neuropat[,c("DX_bl", "PTEDUCAT_bl", "PTGENDER_bl", "fulllewy", 
                                                            "fulltdp43", "fullcaa", "CDGLOBAL_bl", 
                                                            "AGE_bl", "MMSE_bl", "AmyPos_full_bl", "TauPos_full_bl")]), ]

all.neuropat.bl <- subset(all.neuropat, M_vis == 0)
all.neuropat.bl$RID <- factor(all.neuropat.bl$RID)

testneurodata <- all.neuropat.bl[,c("RID", "PTEDUCAT_bl", "PTGENDER_bl", "fulllewy", 
                                    "fulltdp43", "fullcaa", 
                                    "AGE_bl")]
testneurodata <- StratifyContVar(testneurodata, stratcols =  c("PTEDUCAT_bl",  "AGE_bl"))
testneurodata.rids <- RandomizeTreatment2(testneurodata, stratcolumns = c("PTEDUCAT_bl_strat", "PTGENDER_bl", "fulllewy", 
                                                                          "fulltdp43", "fullcaa", 
                                                                          "AGE_bl_strat"), longdata = all.neuropat)







autopsy.data <- all.neuropat[complete.cases(all.neuropat[,c("NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                            "NPLBOD", "NPAMY")]),]
imputed.data <- all.neuropat[complete.cases(all.neuropat[,c("TDP43", "LEWY", "CAA")]),]


autopsy.data.bl <- subset(autopsy.data, M_vis==0)
imputed.data.bl <- subset(imputed.data, M_vis==0)
autopsy.data.bl$RID         <- factor(autopsy.data.bl$RID)
autopsy.data.bl$fulllewy    <- factor(autopsy.data.bl$fulllewy)
autopsy.data.bl$fulltdp43   <- factor(autopsy.data.bl$fulltdp43)
autopsy.data.bl$fullcaa     <- factor(autopsy.data.bl$fullcaa)
autopsy.data.bl$DX_bl       <- factor(autopsy.data.bl$DX_bl)
autopsy.data.bl$PTGENDER    <- factor(autopsy.data.bl$PTGENDER)
autopsy.data.bl$AmyPos_full    <- factor(autopsy.data.bl$AmyPos_full)
autopsy.data.bl$TauPos_full    <- factor(autopsy.data.bl$TauPos_full)
autopsy.data$L
autopsy.data.bl$PTGENDER    <- factor(autopsy.data.bl$PTGENDER)
autopsy.desc        <- table1(autopsy.data.bl[,c("DX_bl", "AGE_bl", "PTGENDER", "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl", 
                                           "mPACCtrailsB_bl", "TauPos_full", "AmyPos_full","fullcaa", "fulltdp43", "fulllewy", "VentricalSum")], splitby = ~DX_bl)$Table1

View(autopsy.desc)

cns  <- which(autopsy.data.bl$DX_bl=="CN")
mcis <- which(autopsy.data.bl$DX_bl=="MCI")
ads  <- which(autopsy.data.bl$DX_bl=="Dementia")

cns <- cns[1:10]
length(cns)
mcis <- mcis[1:8]
ads  <- ads[1:11]

balanced.autopsy <- autopsy.data.bl[c(cns, mcis, ads),]
balanced.autopsy$RID <- factor(balanced.autopsy$RID)
balanced.autopsy$Amy_pos_path <- factor(balanced.autopsy$Amy_pos_path)
balanced.autopsy$TAU_pos_path <- factor(balanced.autopsy$TAU_pos_path)
unname(quantile(balanced.autopsy$AGE_bl)[2:5])




balanced.autopsy.desc <- table1(balanced.autopsy[,c("DX_bl", "AGE_bl", "PTGENDER", "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl", 
                                                  "mPACCtrailsB_bl", "TAU_pos_path", "Amy_pos_path","fullcaa", "fulltdp43", "fulllewy")], 
                                                                                                               splitby = ~DX_bl)$Table1
balanced.autopsy.long <- subset(autopsy.data, RID %in% balanced.autopsy$RID)
balanced.autopsy.long <- balanced.autopsy.long[complete.cases(balanced.autopsy.long[,c("RID", "M_vis", "MMSE", "TOTAL11", "CDRSB", "mPACCtrailsB", "ADNI_MEM", "ADNI_EF", "VentricalSum")]), ]
balanced.autopsy.long <- TimeSinceBaseline(balanced.autopsy.long, timecol = "M_vis")
balanced.autopsy.long <- ZscoreAdj(balanced.autopsy.long, c("MMSE", "TOTAL11", "CDRSB", "mPACCtrailsB", "ADNI_MEM", "ADNI_EF", "VentricalSum"))
balanced.autopsy.long$COMPOSITEZSCORE <- (balanced.autopsy.long$MMSE_zscore + balanced.autopsy.long$mPACCtrailsB_zscore+ balanced.autopsy.long$ADNI_MEM+ balanced.autopsy.long$ADNI_EF) / 4
balanced.autopsy.long$new_time  <- (balanced.autopsy.long$new_time / 12)
balanced.autopsy.long$fulllewy  <- factor(balanced.autopsy.long$fulllewy)
balanced.autopsy.long$fulltdp43 <- factor(balanced.autopsy.long$fulltdp43)
balanced.autopsy.long$fullcaa   <- factor(balanced.autopsy.long$fullcaa)
balanced.autopsy.long$Amy_pos_path   <- factor(balanced.autopsy.long$Amy_pos_path)
balanced.autopsy.long$TAU_pos_path   <- factor(balanced.autopsy.long$TAU_pos_path)
balanced.autopsy.long$PTGENDER_bl    <- factor(balanced.autopsy.long$PTGENDER_bl, levels = c("Female", "Male"))
contrasts(balanced.autopsy.long$PTGENDER_bl) <- contr.treatment(levels(balanced.autopsy.long$PTGENDER_bl),
                                           base=which(levels(balanced.autopsy.long$PTGENDER_bl) == 'Female'))



autopsy.lme    <- lmerTest::lmer(LeftMeta ~ new_time + new_time*AGE_bl  + 
                                   new_time*PTGENDER_bl + new_time*PTEDUCAT_bl  + fulllewy*new_time  + 
                                   fullcaa*new_time     + fulltdp43*new_time    + fulltdp43*new_time + 
                                   AmyPos_full*new_time + TauPos_full*new_time  + (1 | RID),  data = balanced.autopsy.long)

summary(autopsy.lme)







autopsy.list <- list("unbalanced.autopsy.table" = autopsy.desc,
                     "balanced.autopsy.table" = balanced.autopsy.desc,
                     "balanced.autopsy.lme"   = autopsy.lme)




summary(autopsy.list$balanced.autopsy.lme)
imputed.data.bl$RID            <- factor(imputed.data.bl$RID)
imputed.data.bl$fulllewy       <- factor(imputed.data.bl$fulllewy)
imputed.data.bl$fulltdp43      <- factor(imputed.data.bl$fulltdp43)
imputed.data.bl$fullcaa        <- factor(imputed.data.bl$fullcaa)
imputed.data.bl$DX_bl          <- factor(imputed.data.bl$DX_bl)
imputed.data.bl$PTGENDER       <- factor(imputed.data.bl$PTGENDER)
imputed.data.bl$AmyPos_full    <- factor(imputed.data.bl$AmyPos_full)
imputed.data.bl$TauPos_full    <- factor(imputed.data.bl$TauPos_full)

imputed.data.bl$PTGENDER    <- factor(imputed.data.bl$PTGENDER)
imputed.desc        <- table1(imputed.data.bl[,c("DX_bl", "AGE_bl", "PTGENDER", "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl", 
                                                 "mPACCtrailsB_bl", "TauPos_full", "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")], splitby = ~DX_bl)$Table1


cns.imputed  <- which(imputed.data.bl$DX_bl=="CN")
mcis.imputed <- which(imputed.data.bl$DX_bl=="MCI")
ads.imputed  <- which(imputed.data.bl$DX_bl=="Dementia")

cns.imputed <- cns.imputed[1:190]
mcis.imputed <- mcis.imputed[1:175]
ads.imputed <- ads.imputed[1:200]

balanced.imputed <- imputed.data.bl[c(cns.imputed, mcis.imputed, ads.imputed),]
balanced.imputed$RID <- factor(balanced.imputed$RID)
balanced.imputed$AmyPos_full <- factor(balanced.imputed$AmyPos_full)
balanced.imputed$TauPos_full <- factor(balanced.imputed$TauPos_full)

balanced.imputed.desc <- table1(balanced.imputed[,c("DX_bl", "AGE_bl", "PTGENDER", "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl", 
                                                    "mPACCtrailsB_bl", "AmyPos_full", "TauPos_full","fullcaa", "fulltdp43", "fulllewy")], 
                                splitby = ~DX_bl)$Table1



balanced.imputed.long <- subset(imputed.data, RID %in% balanced.imputed$RID)
balanced.imputed.long <- balanced.imputed.long[complete.cases(balanced.imputed.long[,c("RID", "M_vis", "VentricalSum")]), ]
balanced.imputed.long <- TimeSinceBaseline(balanced.imputed.long, timecol = "M_vis")
balanced.imputed.long <- ZscoreAdj(balanced.imputed.long, c( "VentricalSum"))
#balanced.imputed.long$COMPOSITEZSCORE <- (balanced.imputed.long$MMSE_zscore + balanced.imputed.long$mPACCtrailsB_zscore+ balanced.imputed.long$ADNI_MEM+ balanced.imputed.long$ADNI_EF) / 4
balanced.imputed.long$new_time   <- (balanced.imputed.long$new_time / 12)
balanced.imputed.long$fulllewy  <- factor(balanced.imputed.long$fulllewy)
balanced.imputed.long$fulltdp43 <- factor(balanced.imputed.long$fulltdp43)
balanced.imputed.long$fullcaa   <- factor(balanced.imputed.long$fullcaa)
balanced.imputed.long$AmyPos_full <- factor(balanced.imputed.long$AmyPos_full)
balanced.imputed.long$TauPos_full <- factor(balanced.imputed.long$TauPos_full)
balanced.imputed.long$RID <- factor(balanced.imputed.long$RID)
which(is.na(balanced.imputed.long[,c("DX_bl", "AGE_bl", "PTGENDER", "PTEDUCAT_bl",
                                "VentricalSum", "AmyPos_full", "TauPos_full","fullcaa", "fulltdp43", "fulllewy")]), arr.ind = TRUE)

imputed.lme    <- lmerTest::lmer(VentricalSum ~    new_time + new_time*AGE_bl   + new_time*PTGENDER_bl +  new_time*PTEDUCAT_bl  + fulllewy*new_time  + 
                                   fullcaa*new_time     + fulltdp43*new_time    + 
                                   AmyPos_full*new_time + TauPos_full*new_time  +
                                     (1 + new_time|RID)  ,  data = balanced.imputed.long, control = lmerControl(optimizer = "Nelder_Mead"))


summary(imputed.lme)

rem.names <- names(which(unlist(Map(nrow, split(balanced.imputed.long, balanced.imputed.long$RID))) <= 2))
test.remove <- subset(balanced.imputed.long, RID %notin% rem.names)
test.remove$RID <- factor(test.remove$RID)
unlist(Map(nrow, split(test.remove, test.remove$RID)))
nrow(test.remove)
View(test.remove[1:50,])
imputed.lme2    <-  lmerTest::lmer(VentricalSum    ~  new_time + new_time*AGE_bl  + new_time*PTGENDER_bl +  new_time*PTEDUCAT_bl  + fulllewy*new_time  + 
                                     fullcaa*new_time     + fulltdp43*new_time    + 
                                     AmyPos_full*new_time + TauPos_full*new_time  +
                                     (1 + new_time|RID)  ,  data = test.remove)

summary(imputed.lme2)


imputed.list <- list("unbalanced.imputed.table" = imputed.desc,
                     "balanced.imputed.table"   = balanced.imputed.desc,
                     "balanced.imputed.lme"     = imputed.lme)



saveRDS(list("autopsy" = autopsy.list,
             "imputed" = imputed.list), "/Users/adamgabriellang/Desktop/clinical_trial_sim/save_ptau_and_total/aut_imp.rds")



necc.cols.t11   <- c("RID", "M_vis", "TOTAL11")
necc.cols.cdr   <- c("RID", "M_vis", "CDRSB")
necc.cols.mmse  <- c("RID", "M_vis", "MMSE")
necc.cols.mpacc <- c("RID", "M_vis", "mPACCtrailsB")
necc.cols.mem   <- c("RID", "M_vis", "ADNI_MEM", "ADNI_EF")

all.neuropat$RID               <- factor(all.neuropat$RID)
all.neuropat$fulllewy          <- factor(all.neuropat$fulllewy)
all.neuropat$fulltdp43         <- factor(all.neuropat$fulltdp43)
all.neuropat$fullcaa           <- factor(all.neuropat$fullcaa)
all.neuropat$AmyPos_full_bl    <- factor(all.neuropat$AmyPos_full_bl)
all.neuropat$TauPos_full_bl    <- factor(all.neuropat$TauPos_full_bl)
all.neuropat$Amy_pos_path      <- factor(all.neuropat$Amy_pos_path)
all.neuropat$TAU_pos_path      <- factor(all.neuropat$TAU_pos_path)
all.neuropat$PTGENDER          <- factor(all.neuropat$PTGENDER)


combined.11    <- all.neuropat[complete.cases(all.neuropat[,necc.cols.t11]),]
combined.cdr   <- all.neuropat[complete.cases(all.neuropat[,necc.cols.cdr]),]
combined.mmse  <- all.neuropat[complete.cases(all.neuropat[,necc.cols.mmse]),]
combined.mpacc <- all.neuropat[complete.cases(all.neuropat[,necc.cols.mpacc]),]
combined.mem   <- all.neuropat[complete.cases(all.neuropat[,necc.cols.mem]),]

full.data.list <- list("combined.11"      = combined.11,
                       "combined.cdr"     = combined.cdr,
                       "combined.mmse"    = combined.mmse,
                       "combined.mpacc"   = combined.mpacc,
                       "combined.mem"     = combined.mem)

full.data.list <- Map(QuickAdjust, full.data.list)



control.cohort                     <- lapply(full.data.list, function(x)  subset(x, new_time==0 &  DX_bl=="CN" & AmyPos_full_bl==0 & TauPos_full_bl == 0 & CDGLOBAL_bl == 0 & fulllewy == 0 & fulltdp43 == 0 & fullcaa == 0))
control.cohort.long                <- purrr::map2(control.cohort, full.data.list, PullLongData)

mci.scen1.generic                  <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30))
mci.scen1.generic.tplus            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & TauPos_full_bl==1 &  AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30))
mci.scen1.earlyage                 <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & CDGLOBAL_bl == 0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & AGE >=50 & AGE <= 65))

mci.scen1.generic.long             <- purrr::map2(mci.scen1.generic, full.data.list, PullLongData)
mci.scen1.generic.long.tplus       <- purrr::map2(mci.scen1.generic.tplus, full.data.list, PullLongData)
mci.scen1.earlyage.long            <- purrr::map2(mci.scen1.earlyage, full.data.list, PullLongData)

early.ad.scen1.generic            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.generic.tplus      <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & TauPos_full_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.earlyage           <- lapply(full.data.list, function(x)  subset(x,  new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & AGE_bl >=55 & AGE_bl <= 65 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))

early.ad.scen1.generic.long        <- purrr::map2(early.ad.scen1.generic, full.data.list, PullLongData)
early.ad.scen1.generic.long.tplus  <- purrr::map2(early.ad.scen1.generic.tplus, full.data.list, PullLongData)
early.ad.scen1.earlyage.long       <- purrr::map2(early.ad.scen1.earlyage, full.data.list, PullLongData)


ad.scen1.generic                   <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 28 & AGE_bl >=65 & AGE_bl <= 85))
ad.scen1.generic.tplus             <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & TauPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 28 & AGE_bl >=65 & AGE_bl <= 85))
ad.scen1.earlyage                  <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 28 & AGE_bl >=55 & AGE_bl <= 65))

ad.scen1.generic.long              <- purrr::map2(ad.scen1.generic, full.data.list, PullLongData)
ad.scen1.generic.long.tplus        <- purrr::map2(ad.scen1.generic.tplus, full.data.list, PullLongData)
ad.scen1.earlyage.long             <- purrr::map2(ad.scen1.earlyage, full.data.list, PullLongData)


#Controls
lme.control.combined.total11 <- lmerTest::lmer(TOTAL11      ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl+  (1|RID), data = control.cohort.long[["combined.11"]])
lme.control.combined.cdrsb   <- lmerTest::lmer(CDRSB        ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl+  (1|RID), data = control.cohort.long[["combined.cdr"]])
lme.control.combined.mmse    <- lmerTest::lmer(MMSE         ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl+  (1|RID), data = control.cohort.long[["combined.mmse"]])
lme.control.combined.mpacc   <- lmerTest::lmer(mPACCtrailsB ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl+  (1|RID), data = control.cohort.long[["combined.mpacc"]])
lme.control.combined.mem     <- lmerTest::lmer(ADNI_MEM     ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl+  (1|RID), data = control.cohort.long[["combined.mem"]])
lme.control.combined.ef      <- lmerTest::lmer(ADNI_MEM     ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl+  (1|RID), data = control.cohort.long[["combined.mem"]])
lme.control.combined.total13 <- lmerTest::lmer(TOTAL13      ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl+  (1|RID), data = control.cohort.long[["combined.13"]])


#MCI
lme.mci.total11  <- lmerTest::lmer(TOTAL11      ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.generic.long[["combined.11"]])
lme.mci.cdrsb    <- lmerTest::lmer(CDRSB        ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.generic.long[["combined.cdr"]])
lme.mci.mmse     <- lmerTest::lmer(MMSE         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.generic.long[["combined.mmse"]])
lme.mci.mpacc    <- lmerTest::lmer(mPACCtrailsB ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.generic.long[["combined.mpacc"]])


lme.mci.tplus.total11  <- lmerTest::lmer(TOTAL11      ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID),  data = mci.scen1.generic.long.tplus[["combined.11"]])
lme.mci.tplus.cdrsb    <- lmerTest::lmer(CDRSB        ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID),  data = mci.scen1.generic.long.tplus[["combined.cdr"]])
lme.mci.tplus.mmse     <- lmerTest::lmer(MMSE         ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID),  data = mci.scen1.generic.long.tplus[["combined.mmse"]])
lme.mci.tplus.mpacc    <- lmerTest::lmer(mPACCtrailsB ~ new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID),  data = mci.scen1.generic.long.tplus[["combined.mpacc"]])

lme.mci.earlyage.total11  <- lmerTest::lmer(TOTAL11      ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.earlyage.long[["combined.11"]])
lme.mci.earlyage.cdrsb    <- lmerTest::lmer(CDRSB        ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.earlyage.long[["combined.cdr"]])
lme.mci.earlyage.mmse     <- lmerTest::lmer(MMSE         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.earlyage.long[["combined.mmse"]])
lme.mci.earlyage.mpacc    <- lmerTest::lmer(mPACCtrailsB ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = mci.scen1.earlyage.long[["combined.mpacc"]])

#AD

lme.ad.total11  <- lmerTest::lmer(TOTAL11      ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long[["combined.11"]])
lme.ad.cdrsb    <- lmerTest::lmer(CDRSB        ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long[["combined.cdr"]])
lme.ad.mmse     <- lmerTest::lmer(MMSE         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long[["combined.mmse"]])
lme.ad.mpacc    <- lmerTest::lmer(mPACCtrailsB ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long[["combined.mpacc"]])

lme.ad.tplus.total11  <- lmerTest::lmer(TOTAL11         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long.tplus[["combined.11"]])
lme.ad.tplus.cdrsb    <- lmerTest::lmer(CDRSB           ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long.tplus[["combined.cdr"]])
lme.ad.tplus.mmse        <- lmerTest::lmer(MMSE         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long.tplus[["combined.mmse"]])
lme.ad.tplus.mpacc       <- lmerTest::lmer(mPACCtrailsB ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.generic.long.tplus[["combined.mpacc"]])

lme.ad.earlyage.total11  <- lmerTest::lmer(TOTAL11      ~   new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.earlyage.long[["combined.11"]])
lme.ad.earlyage.cdrsb    <- lmerTest::lmer(CDRSB        ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.earlyage.long[["combined.cdr"]])
lme.ad.earlyage.mmse     <- lmerTest::lmer(MMSE         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.earlyage.long[["combined.mmse"]])
lme.ad.earlyage.mpacc    <- lmerTest::lmer(mPACCtrailsB ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = ad.scen1.earlyage.long[["combined.mpacc"]])




#Early AD

lme.early.ad.total11  <- lmerTest::lmer(TOTAL11      ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long[["combined.11"]])
lme.early.ad.cdrsb    <- lmerTest::lmer(CDRSB        ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long[["combined.cdr"]])
lme.early.ad.mmse     <- lmerTest::lmer(MMSE         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long[["combined.mmse"]])
lme.early.ad.mpacc    <- lmerTest::lmer(mPACCtrailsB ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long[["combined.mpacc"]])

lme.early.ad.tplus.total11     <- lmerTest::lmer(TOTAL11       ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long.tplus[["combined.11"]])
lme.early.ad.tplus.cdrsb       <- lmerTest::lmer(CDRSB         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long.tplus[["combined.cdr"]])
lme.early.ad.tplus.mmse        <- lmerTest::lmer(MMSE          ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long.tplus[["combined.mmse"]])
lme.early.ad.tplus.mpacc       <- lmerTest::lmer(mPACCtrailsB  ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.generic.long.tplus[["combined.mpacc"]])

lme.early.ad.earlyage.total11  <- lmerTest::lmer(TOTAL11       ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.earlyage.long[["combined.11"]])
lme.early.ad.earlyage.cdrsb    <- lmerTest::lmer(CDRSB         ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.earlyage.long[["combined.cdr"]])
lme.early.ad.earlyage.mmse     <- lmerTest::lmer(MMSE          ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.earlyage.long[["combined.mmse"]])
lme.early.ad.earlyage.mpacc    <- lmerTest::lmer(mPACCtrailsB  ~  new_time + PTEDUCAT_bl + AGE_bl +PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)  ,  data = early.ad.scen1.earlyage.long[["combined.mpacc"]])


g1.mci <- "MCI A+ CDR = .5 \nAge 65-85  MMMSE 24-30 \n"
g2.mci <- "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n"
g3.mci <- "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n"

g1.ad <- "AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n"
g2.ad <- "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n"
g3.ad <- "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n"


g1.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n"
g2.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"
g3.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n"



total11.lmes <- list( "MCI A+ CDR = .5 \nAge 65-85  MMMSE 24-30 \n" = lme.mci.total11, 
                      "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n" = lme.mci.tplus.total11,
                      "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n" = lme.mci.earlyage.total11,
                      "AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n"= lme.ad.total11, 
                      "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n" = lme.ad.tplus.total11,
                     "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n" = lme.ad.earlyage.total11,
                     "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n" = lme.early.ad.total11, 
                     "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n" = lme.early.ad.tplus.total11,
                     "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"  = lme.early.ad.earlyage.total11)


cdrsb.lmes <- list( "MCI A+ CDR = .5 \nAge 65-85  MMMSE 24-30 \n" = lme.mci.cdrsb, 
                      "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n" = lme.mci.tplus.cdrsb,
                      "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n" = lme.mci.earlyage.cdrsb,
                      "AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n"= lme.ad.cdrsb, 
                      "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n" = lme.ad.tplus.cdrsb,
                      "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n" = lme.ad.earlyage.cdrsb,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n" = lme.early.ad.cdrsb, 
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n" = lme.early.ad.tplus.cdrsb,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"  = lme.early.ad.earlyage.cdrsb)

mmse.lmes <- list( "MCI A+ CDR = .5 \nAge 65-85  MMMSE 24-30 \n" = lme.mci.mmse, 
                      "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n" = lme.mci.tplus.mmse,
                      "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n" = lme.mci.earlyage.mmse,
                      "AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n"= lme.ad.mmse, 
                      "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n" = lme.ad.tplus.mmse,
                      "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n" = lme.ad.earlyage.mmse,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n" = lme.early.ad.mmse, 
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n" = lme.early.ad.tplus.mmse,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"  = lme.early.ad.earlyage.mmse)

mpacc.lmes <- list( "MCI A+ CDR = .5 \nAge 65-85  MMMSE 24-30 \n" = lme.mci.mpacc, 
                      "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n" = lme.mci.tplus.mpacc,
                      "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n" = lme.mci.earlyage.mpacc,
                      "AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n"= lme.ad.mpacc, 
                      "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n" = lme.ad.tplus.mpacc,
                      "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n" = lme.ad.earlyage.mpacc,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n" = lme.early.ad.mpacc, 
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n" = lme.early.ad.tplus.mpacc,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"  = lme.early.ad.earlyage.mpacc)




totall11.data <- list(g1.mci = mci.scen1.generic.long[["combined.11"]], 
                      "mci.tplus" = mci.scen1.generic.long.tplus[["combined.11"]],
                      "mci.earlyage" = mci.scen1.earlyage.long[["combined.11"]],
                      "ad" = ad.scen1.generic.long[["combined.11"]], 
                      "ad.tplus" =ad.scen1.generic.long.tplus[["combined.11"]],
                      "ad.earlyage" = ad.scen1.earlyage.long[["combined.11"]],
                      "early.ad" = early.ad.scen1.generic.long[["combined.11"]], 
                      "early.ad.tplus" = early.ad.scen1.generic.long.tplus[["combined.11"]],
                      "early.ad.earlyage" = early.ad.scen1.earlyage.long[["combined.11"]])


cdrsb.data <- list("mci" = mci.scen1.generic.long[["combined.cdr"]], 
                      "mci.tplus" = mci.scen1.generic.long.tplus[["combined.cdr"]],
                      "mci.earlyage" = mci.scen1.earlyage.long[["combined.cdr"]],
                      "ad" = ad.scen1.generic.long[["combined.cdr"]], 
                      "ad.tplus" =ad.scen1.generic.long.tplus[["combined.cdr"]],
                      "ad.earlyage" = ad.scen1.earlyage.long[["combined.cdr"]],
                      "early.ad" = early.ad.scen1.generic.long[["combined.cdr"]], 
                      "early.ad.tplus" = early.ad.scen1.generic.long.tplus[["combined.cdr"]],
                      "early.ad.earlyage" = early.ad.scen1.earlyage.long[["combined.cdr"]])

mmse.data <- list("mci" = mci.scen1.generic.long[["combined.mmse"]], 
                      "mci.tplus" = mci.scen1.generic.long.tplus[["combined.mmse"]],
                      "mci.earlyage" = mci.scen1.earlyage.long[["combined.mmse"]],
                      "ad" = ad.scen1.generic.long[["combined.mmse"]], 
                      "ad.tplus" =ad.scen1.generic.long.tplus[["combined.mmse"]],
                      "ad.earlyage" = ad.scen1.earlyage.long[["combined.mmse"]],
                      "early.ad" = early.ad.scen1.generic.long[["combined.mmse"]], 
                      "early.ad.tplus" = early.ad.scen1.generic.long.tplus[["combined.mmse"]],
                      "early.ad.earlyage" = early.ad.scen1.earlyage.long[["combined.mmse"]])

mpacc.data <- list("mci" = mci.scen1.generic.long[["combined.mpacc"]], 
                      "mci.tplus" = mci.scen1.generic.long.tplus[["combined.mpacc"]],
                      "mci.earlyage" = mci.scen1.earlyage.long[["combined.mpacc"]],
                      "ad" = ad.scen1.generic.long[["combined.mpacc"]], 
                      "ad.tplus" =ad.scen1.generic.long.tplus[["combined.mpacc"]],
                      "ad.earlyage" = ad.scen1.earlyage.long[["combined.mpacc"]],
                      "early.ad" = early.ad.scen1.generic.long[["combined.mpacc"]], 
                      "early.ad.tplus" = early.ad.scen1.generic.long.tplus[["combined.mpacc"]],
                      "early.ad.earlyage" = early.ad.scen1.earlyage.long[["combined.mpacc"]])

total11sig <- Map(BuildSignificanceTable, total11.lmes)
cdrsbsig   <- Map(BuildSignificanceTable, cdrsb.lmes)
mmsesig    <- Map(BuildSignificanceTable, mmse.lmes)
mpaccsig   <- Map(BuildSignificanceTable, mpacc.lmes)


#ADAS 11

Randomize11 <- Map(RandomizeTreatment, totall11.data)
BuildModels11 <- purrr::map2(total11.lmes, Randomize11, BuildSimulationModel, formula="TOTAL11 ~ new_time*treat + (1|RID)" )
DiseasePlots11 <- purrr::map2(Randomize11, BuildModels11, SampleSizeSimulation2, formula = "TOTAL11 ~ new_time*treat + (1|RID)", fcompare = "TOTAL11 ~ new_time + (1|RID)",
                                                          efficacy=.5, breaks=c(100,200,300), yaxislab_dpm="ADAS11", return_dpm=TRUE)
names(DiseasePlots11) <- names(total11.lmes)
mcis.11 <- GroupDiseaseTraj(DiseasePlots11[1:3],  yaxislab_dpm="ADAS11")
ads.11 <- GroupDiseaseTraj(DiseasePlots11[4:6],  yaxislab_dpm="ADAS11")
early.ads.11 <- GroupDiseaseTraj(DiseasePlots11[7:9],  yaxislab_dpm="ADAS11")


#CDRSB

RandomizeCDR <- Map(RandomizeTreatment, cdrsb.data)
BuildModelsCDR <- purrr::map2(cdrsb.lmes, RandomizeCDR, BuildSimulationModel, formula="CDRSB ~ new_time*treat + (1|RID)" )
DiseasePlotsCDR <- purrr::map2(RandomizeCDR, BuildModelsCDR, SampleSizeSimulation2, formula = "CDRSB ~ new_time*treat + (1|RID)", fcompare = "CDRSB ~ new_time + (1|RID)",
                              efficacy=.5, breaks=c(100,200,300), yaxislab_dpm="CDR Sum of Boxes", return_dpm=TRUE)
names(DiseasePlotsCDR) <- names(cdrsb.lmes)
mcis.cdr <- GroupDiseaseTraj(DiseasePlotsCDR[1:3],  yaxislab_dpm="CDR Sum of Boxes")
ads.cdr <- GroupDiseaseTraj(DiseasePlotsCDR[4:6],  yaxislab_dpm="CDR Sum of Boxes")
early.ads.cdr <- GroupDiseaseTraj(DiseasePlotsCDR[7:9],  yaxislab_dpm="CDR Sum of Boxes")



#MMSE

Randomizemmse <- Map(RandomizeTreatment, mmse.data)
BuildModelsmmse <- purrr::map2(mmse.lmes, Randomizemmse, BuildSimulationModel, formula="MMSE ~ new_time*treat + (1|RID)" )
DiseasePlotsmmse <- purrr::map2(Randomizemmse, BuildModelsmmse, SampleSizeSimulation2, formula = "MMSE ~ new_time*treat + (1|RID)", fcompare = "MMSE ~ new_time + (1|RID)",
                               efficacy=.5, breaks=c(100,200,300), yaxislab_dpm="MMSE", return_dpm=TRUE)
names(DiseasePlotsmmse) <- names(mmse.lmes)
mcis.mmse <- GroupDiseaseTraj(DiseasePlotsmmse[1:3],  yaxislab_dpm="MMSE")
ads.mmse <- GroupDiseaseTraj(DiseasePlotsmmse[4:6],  yaxislab_dpm="MMSE")
early.ads.mmse <- GroupDiseaseTraj(DiseasePlotsmmse[7:9],  yaxislab_dpm="MMSE")


#mPACC

Randomizempacc <- Map(RandomizeTreatment, mpacc.data)
BuildModelsmpacc <- purrr::map2(mpacc.lmes, Randomizempacc, BuildSimulationModel, formula="mPACCtrailsB ~ new_time*treat + (1|RID)" )
DiseasePlotsmpacc <- purrr::map2(Randomizempacc, BuildModelsmpacc, SampleSizeSimulation2, formula = "mPACCtrailsB ~ new_time*treat + (1|RID)", fcompare = "mPACCtrailsB ~ new_time + (1|RID)",
                                efficacy=.5, breaks=c(100,200,300), yaxislab_dpm="mPACCtrailsB", return_dpm=TRUE)
names(DiseasePlotsmpacc) <- names(mpacc.lmes)
mcis.mpacc <- GroupDiseaseTraj(DiseasePlotsmpacc[1:3],  yaxislab_dpm="mPACCtrailsB")
ads.mpacc <- GroupDiseaseTraj(DiseasePlotsmpacc[4:6],  yaxislab_dpm="mPACCtrailsB")
early.ads.mpacc <- GroupDiseaseTraj(DiseasePlotsmpacc[7:9],  yaxislab_dpm="mPACCtrailsB")



savedpmlist <- list("ADAS" = list("MCI" = mcis.11,
                                  "AD" = ads.11,
                                  "Early AD" = early.ads.11),
                    "CDRSB"= list("MCI" = mcis.cdr,
                                  "AD" = ads.cdr,
                                  "Early AD" = early.ads.cdr),
                    "MMSE"= list("MCI" = mcis.mmse,
                                  "AD" = ads.mmse,
                                  "Early AD" = early.ads.mmse),
                    "mPACCtrailsB"= list("MCI" = mcis.mpacc,
                                 "AD" = ads.mpacc,
                                 "Early AD" = early.ads.mpacc))


saveRDS(savedpmlist, "/Users/adamgabriellang/Desktop/clinical_trial_sim/save_ptau_and_total/save_dpms.rds")









#ADAS 11
adas.data.complete <- adas_merge_demog[,c("RID","M_vis", "VISCODE", "DX", "COLPROT", "ORIGPROT", 
                                         "VISCODE", "EXAMDATE_adnimerge", "AGE", 
                                         "PTGENDER", "PTEDUCAT", 
                                         "APOE4", "CDRSB", "MMSE", "TOTAL11", "mPACCtrailsB",baseline.var.list,
                                         "adas_pet_valid", "adas_csf_valid", "EXAMDATE_pet", "M_vis",
                                         "ABETA", "TAU", "PTAU", "SUMMARYSUVR_COMPOSITE_REFNORM", "AmyPos", "ptau_pos","TDP43", "LEWY", "CAA", "NPBRAAK", "NPNEUR","NPTDPA", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                         "NPLBOD", "NPAMY", "fulllewy", "fulltdp43", "fullcaa")]
necc.cols.t11 <- c("RID", "VISCODE", "TOTAL11", "AmyPos_bl", "ptau_pos_bl", "CDGLOBAL_bl", "AGE_bl", "MMSE_bl")
adas.neuropath.outcome.rows <- which(complete.cases(adas.data.complete[,c(necc.cols.t11,"NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                         "NPLBOD", "NPAMY")]) | complete.cases(adas.data.complete[,c(necc.cols.t11,"TDP43", "LEWY", "CAA")]))
adas.neuropath.outcome <- adas.data.complete[adas.neuropath.outcome.rows,]



#CDRSB
cdr.data.complete <- adas_merge_demog[,c("RID","M_vis", "VISCODE", "DX", "COLPROT", "ORIGPROT", 
                                          "VISCODE", "EXAMDATE_adnimerge", "AGE", 
                                          "PTGENDER", "PTEDUCAT", 
                                          "APOE4", "CDRSB", "MMSE", "TOTAL11", "mPACCtrailsB", baseline.var.list,
                                          "adas_pet_valid", "adas_csf_valid", "EXAMDATE_pet", "M_vis",
                                          "ABETA", "TAU", "PTAU", "SUMMARYSUVR_COMPOSITE_REFNORM", "AmyPos", "ptau_pos","TDP43", "LEWY", "CAA", "NPBRAAK", "NPNEUR","NPTDPA", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                          "NPLBOD", "NPAMY", "fulllewy", "fulltdp43", "fullcaa")]
necc.cols.cdr <- c("RID", "VISCODE", "CDRSB", "AmyPos_bl", "ptau_pos_bl", "CDGLOBAL_bl", "AGE_bl",  "MMSE_bl")
cdr.neuropath.outcome.rows <- which(complete.cases(cdr.data.complete[,c(necc.cols.cdr,"NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                                          "NPLBOD", "NPAMY")]) | complete.cases(cdr.data.complete[,c(necc.cols.cdr,"TDP43", "LEWY", "CAA")]))
cdr.neuropath.outcome <- cdr.data.complete[cdr.neuropath.outcome.rows,]


#MMSE
mmse.data.complete <- adas_merge_demog[,c("RID","M_vis", "VISCODE", "DX", "COLPROT", "ORIGPROT", 
                                         "VISCODE", "EXAMDATE_adnimerge", "AGE", 
                                         "PTGENDER", "PTEDUCAT", 
                                         "APOE4", "CDRSB", "MMSE", "TOTAL11", "mPACCtrailsB",baseline.var.list,
                                         "adas_pet_valid", "adas_csf_valid", "EXAMDATE_pet", "M_vis",
                                         "ABETA", "TAU", "PTAU", "SUMMARYSUVR_COMPOSITE_REFNORM", "AmyPos", "ptau_pos","TDP43", "LEWY", "CAA", "NPBRAAK", "NPNEUR","NPTDPA", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                         "NPLBOD", "NPAMY", "fulllewy", "fulltdp43", "fullcaa")]
necc.cols.mmse <- c("RID", "VISCODE", "MMSE", "AmyPos_bl", "ptau_pos_bl", "CDGLOBAL_bl", "AGE_bl",  "MMSE_bl")
mmse.neuropath.outcome.rows <- which(complete.cases(mmse.data.complete[,c(necc.cols.t11,"NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                                        "NPLBOD", "NPAMY")]) | complete.cases(mmse.data.complete[,c(necc.cols.t11,"TDP43", "LEWY", "CAA")]))
mmse.neuropath.outcome <- mmse.data.complete[mmse.neuropath.outcome.rows,]


mpacc.data.complete <- adas_merge_demog[,c("RID","M_vis", "VISCODE", "DX", "COLPROT", "ORIGPROT", 
                                          "VISCODE", "EXAMDATE_adnimerge", "AGE", 
                                          "PTGENDER", "PTEDUCAT", 
                                          "APOE4", "CDRSB", "MMSE", "TOTAL11", "mPACCtrailsB",baseline.var.list,
                                          "adas_pet_valid", "adas_csf_valid", "EXAMDATE_pet", "M_vis",
                                          "ABETA", "TAU", "PTAU", "SUMMARYSUVR_COMPOSITE_REFNORM", "AmyPos", "ptau_pos","TDP43", "LEWY", "CAA", "NPBRAAK", "NPNEUR","NPTDPA", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                          "NPLBOD", "NPAMY", "fulllewy", "fulltdp43", "fullcaa")]
necc.cols.mpacc <- c("RID", "VISCODE", "mPACCtrailsB", "AmyPos_bl", "ptau_pos_bl", "CDGLOBAL_bl", "AGE_bl",  "MMSE_bl")
mpacc.neuropath.outcome.rows <- which(complete.cases(mpacc.data.complete[,c(necc.cols.t11,"NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                                          "NPLBOD", "NPAMY")]) | complete.cases(mpacc.data.complete[,c(necc.cols.t11,"TDP43", "LEWY", "CAA")]))
mpacc.neuropath.outcome <- mpacc.data.complete[mpacc.neuropath.outcome.rows,]






#ADNIMEM
mem.data.complete <- adas_merge_demog[,c("RID", "VISCODE", "DX", "COLPROT", "ORIGPROT", 
                                             "VISCODE", "EXAMDATE_adnimerge", "AGE", 
                                             "PTGENDER", "PTEDUCAT", 
                                             "APOE4", "ADNI_MEM", "ADNI_EF",
                                             "adas_pet_valid", "adas_csf_valid", "CDRSB", "MMSE","EXAMDATE_pet", "M_vis","mPACCtrailsB",baseline.var.list,
                                             "ABETA", "TAU", "PTAU", "SUMMARYSUVR_COMPOSITE_REFNORM", "AmyPos", "ptau_pos","TDP43", "LEWY", "CAA",
                                             "NPBRAAK", "NPNEUR","NPTDPA", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                             "NPLBOD", "NPAMY", "fulllewy", "fulltdp43", "fullcaa")]

necc.cols.mem <- c("RID", "VISCODE", "ADNI_MEM", "ADNI_EF", "AmyPos_bl", "ptau_pos_bl", "CDGLOBAL_bl", "AGE_bl", "mPACCtrailsB")
mem.outcome.rows <- which(complete.cases(mem.data.complete[,c(necc.cols.mem,"NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                                         "NPLBOD", "NPAMY")]) | complete.cases(mem.data.complete[,c(necc.cols.mem,"TDP43", "LEWY", "CAA")]))
mem.neuropath.outcome <- mem.data.complete[mem.outcome.rows,]
adas.neuropath.outcome  <- SetNeuroData(adas.neuropath.outcome)
cdr.neuropath.outcome   <- SetNeuroData(cdr.neuropath.outcome)
mmse.neuropath.outcome  <- SetNeuroData(mmse.neuropath.outcome)
mpacc.neuropath.outcome <- SetNeuroData(mpacc.neuropath.outcome)
mem.neuropath.outcome   <- SetNeuroData(mem.neuropath.outcome)




adas.neuropath.outcome <- adas.neuropath.outcome[order(adas.neuropath.outcome$RID, adas.neuropath.outcome$M_vis, decreasing = FALSE),]
adas.neuropath.outcome <- QuickAdjust(adas.neuropath.outcome)


cdr.neuropath.outcome <- cdr.neuropath.outcome[order(cdr.neuropath.outcome$RID, cdr.neuropath.outcome$M_vis, decreasing = FALSE),]
cdr.neuropath.outcome <- QuickAdjust(cdr.neuropath.outcome)


mmse.neuropath.outcome <- mmse.neuropath.outcome[order(mmse.neuropath.outcome$RID, mmse.neuropath.outcome$M_vis, decreasing = FALSE),]
mmse.neuropath.outcome <- QuickAdjust(mmse.neuropath.outcome)


mpacc.neuropath.outcome <- mpacc.neuropath.outcome[order(mpacc.neuropath.outcome$RID, mpacc.neuropath.outcome$M_vis, decreasing = FALSE),]
mpacc.neuropath.outcome <- QuickAdjust(mpacc.neuropath.outcome)


mem.neuropath.outcome <- mem.neuropath.outcome[order(mem.neuropath.outcome$RID, mem.neuropath.outcome$M_vis, decreasing = FALSE),]
mem.neuropath.outcome <- QuickAdjust(mem.neuropath.outcome)



all.data.list <- list(adas.neuropath.outcome,cdr.neuropath.outcome,
                      mmse.neuropath.outcome,mpacc.neuropath.outcome,
                      mem.neuropath.outcome)


mci.scen1.generic          <- lapply(all.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_bl==1 & AGE >= 65 & AGE <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30))
mci.scen1.earlyage         <- lapply(all.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_bl==1 & CDGLOBAL_bl == 0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & AGE >=50 & AGE <= 65))
mci.scen1.tplus            <- lapply(all.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_bl==1 & CDGLOBAL_bl == 0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & AGE >=65 & AGE <= 85 & ptau_pos_bl==1))
mci.scen1.neur             <- lapply(all.data.list[c(2,4)], function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_bl==1 & CDGLOBAL_bl == 0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & AGE >=65 & AGE <= 85 & fullcaa==0 & fulltdp43==0 & fulllewy==0))
mci.scen1.neur.tplus       <- lapply(all.data.list[c(2,4)], function(x)subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_bl==1 & CDGLOBAL_bl == 0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & AGE >=65 & AGE <= 85 & fullcaa==0 & fulltdp43==0 & fulllewy==0 & ptau_pos_bl==1))###############)

ad.scen1.generic           <- lapply(all.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 28 & AGE >=65 & AGE <= 85))
ad.scen1.earlyage          <- lapply(all.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 28 & AGE >=55 & AGE <= 65))
ad.scen1.tplus             <- lapply(all.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 28 & AGE >=65 & AGE <= 85  & ptau_pos_bl==1))
ad.scen1.neur              <- lapply(all.data.list[c(2,4)], function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE >=65 & AGE <= 85 & fullcaa==0 & fulltdp43==0 & fulllewy==0))
ad.scen1.neur.tplus        <- lapply(all.data.list[c(2,4)], function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE >=65 & AGE <= 85 & fullcaa==0 & fulltdp43==0 & fulllewy==0 & ptau_pos_bl==1))

early.ad.scen1.generic      <- lapply(all.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_bl==1 & AGE>=65 & AGE <= 85 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.earlyage     <- lapply(all.data.list, function(x)  subset(x,  new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_bl==1 & AGE>=55 & AGE <= 65 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.tplus        <- lapply(all.data.list, function(x)  subset(x,  new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_bl==1 & AGE>=65 & AGE <= 85 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1) & ptau_pos_bl==1))
early.ad.scen1.neur         <- lapply(all.data.list[c(2,4)], function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_bl==1 & AGE>=65 & AGE <= 85 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1) & fullcaa==0 & fulltdp43==0 & fulllewy==0))
early.ad.scen1.neur.tplus   <- lapply(all.data.list[c(2,4)], function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_bl==1 & AGE>=65 & AGE <= 85 & MMSE_bl >=20 & MMSE_bl <= 30 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1) & fullcaa==0 & fulltdp43==0 & fulllewy==0 & ptau_pos_bl==1))##############


g1.mci <- "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 \n"
g2.mci <- "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n"
g3.mci <- "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n"
g4.mci <- "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 \nNo Copathologies (Lewy- TDP43- CAA-)\n"
g5.mci <- "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 \nNo Copathologies (Lewy- TDP43- CAA-) Tau+ (PTAU)\n"

g1.ad <- "AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n"
g2.ad <- "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n"
g3.ad <- "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n"
g4.ad <- "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 \nNo Copathologies (Lewy- TDP43- CAA-)\n"
g5.ad <- "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 \nNo Copathologies (Lewy- TDP43- CAA-) Tau+ (PTAU)\n"

g1.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n"
g2.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"
g3.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n"
g4.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 \nNo Copathologies (Lewy- TDP43- CAA-)\n"
g5.early.ad <- "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 \nNo Copathologies (Lewy- TDP43- CAA-) Tau+ (PTAU)\n"

mci.scen1.generic.long     <- purrr::map2(mci.scen1.generic, all.data.list, PullLongData)
mci.scen1.earlyage.long  <- purrr::map2(mci.scen1.earlyage, all.data.list, PullLongData)
mci.scen1.tplus.long <- purrr::map2(mci.scen1.tplus, all.data.list, PullLongData)
mci.scen1.neur.long <- purrr::map2(mci.scen1.neur, all.data.list[c(2,4)], PullLongData)
mci.scen1.neur.tplus.long <- purrr::map2(mci.scen1.neur.tplus, all.data.list[c(2,4)], PullLongData)

ad.scen1.generic.long    <- purrr::map2(ad.scen1.generic, all.data.list, PullLongData)
ad.scen1.earlyage.long   <- purrr::map2(ad.scen1.earlyage, all.data.list, PullLongData)
ad.scen1.tplus.long      <- purrr::map2(ad.scen1.tplus, all.data.list, PullLongData)
ad.scen1.neur.long       <- purrr::map2(ad.scen1.neur, all.data.list[c(2,4)], PullLongData)
ad.scen1.neur.tplus.long <- purrr::map2(ad.scen1.neur.tplus, all.data.list[c(2,4)], PullLongData)

early.ad.scen1.generic.long<- purrr::map2(early.ad.scen1.generic, all.data.list, PullLongData)
early.ad.scen1.earlyage.long<- purrr::map2(early.ad.scen1.earlyage, all.data.list, PullLongData)
early.ad.scen1.tplus.long<- purrr::map2(early.ad.scen1.tplus, all.data.list, PullLongData)
early.ad.scen1.neur.long<- purrr::map2(early.ad.scen1.neur, all.data.list[c(2,4)], PullLongData)
early.ad.scen1.neur.tplus.long<- purrr::map2(early.ad.scen1.neur.tplus, all.data.list[c(2,4)], PullLongData)

testdata1     <- mci.scen1.generic.long[[1]]
testdata1$new_PTAU_bl <- testdata1$PTAU_bl
testdata1.bl <- subset(testdata1, new_time==0)
testdata1$new_PTAU_bl <- (testdata1$new_PTAU_bl - mean(testdata1.bl$PTAU_bl))
testdata1$RID <- factor(testdata1$RID)
all.rids      <- levels(testdata1[["RID"]])
nrids <- nlevels(testdata1[["RID"]])
n.sample <- round(nrids/2)
treat.rids <- sample(all.rids, n.sample, prob = rep(.5, length(all.rids)), replace = FALSE)
base.rids <- subset(all.rids, all.rids %notin% treat.rids)
testdata1$treat <- rep(NA, nrow(testdata1))
treat.rows <- which(testdata1$RID %in% treat.rids)
base.rows <- which(testdata1$RID %in% base.rids)
testdata1["treat"][treat.rows,] <- 1
testdata1["treat"][base.rows,] <- 0
testdata1$treat <- factor(testdata1$treat)
testdata1$fullcaa <- factor(testdata1$fullcaa)
testdata1$fulllewy <- factor(testdata1$fulllewy)
testdata1$fulltdp43 <- factor(testdata1$fulltdp43)
testdata1$treat <- factor(testdata1$treat)
testmodel1 <- lmer(TOTAL11 ~new_time*PTAU_bl+new_time*fulllewy+new_time*fulltdp43+new_time*fulltdp43+(1|RID), data = testdata1)
fixd       <- fixef(testmodel1)
fixd["treat1"] <- 0
fixd["new_time:treat1"] <- 0
fixd["treat1:PTAU_bl"] <- 0
fixd["new_time:treat1:PTAU_bl"] <- (-1*(fixd["new_time:PTAU_bl"] /2))
check<- summary(testmodel1)$sigma
varcor<- as.numeric(summary(testmodel1)$varcor[[1]])
constr.lme                                          <- makeLmer(TOTAL11 ~ new_time*treat*PTAU_bl + new_time*fulllewy + new_time*fulltdp43 + new_time*fulltdp43 + (1|RID),fixef = fixd, VarCorr=list(varcor), sigma = check, data = testdata1)
sim_ext_class1     <- extend(constr.lme, along="RID", n=max(breaks))
testpower          <- powerCurve(sim_ext_class1, test = compare(TOTAL11 ~ new_time*PTAU_bl + new_time*fulllewy + new_time*fulltdp43 + new_time*fulltdp43 + (1|RID)), along="RID", breaks = breaks) 
testpower



adas.neuropath.outcome$CDGLOBAL_bl <- factor(adas.neuropath.outcome$CDGLOBAL_bl)
adas.neuropath.outcome$PTGENDER_bl <- factor(adas.neuropath.outcome$PTGENDER_bl)
adas.neuropath.outcome$DX_bl <- factor(adas.neuropath.outcome$DX_bl)
adas.neuropath.outcome$AmyPos_bl <- factor(adas.neuropath.outcome$AmyPos_bl)
adas.neuropath.outcome$fullcaa <- factor(adas.neuropath.outcome$fullcaa)
adas.neuropath.outcome$fulltdp43 <- factor(adas.neuropath.outcome$fulltdp43)
adas.neuropath.outcome$fulllewy <- factor(adas.neuropath.outcome$fulllewy)
adas.neuropath.outcome.bl <- subset(adas.neuropath.outcome, new_time==0)
adas.desc.table <- table1(adas.neuropath.outcome.bl[c("DX_bl","AmyPos_bl","PTGENDER_bl","ptau_pos_bl","MMSE_bl","CDGLOBAL_bl","TOTAL11", "fullcaa", "fulltdp43", "fulllewy")], splitby = ~DX_bl)$Table1


View(cdr.neuropath.outcome)
cdr.neuropath.outcome$CDGLOBAL_bl <- factor(cdr.neuropath.outcome$CDGLOBAL_bl)
cdr.neuropath.outcome$PTGENDER_bl <- factor(cdr.neuropath.outcome$PTGENDER_bl)
cdr.neuropath.outcome$DX_bl <- factor(cdr.neuropath.outcome$DX_bl)
cdr.neuropath.outcome$AmyPos_bl <- factor(cdr.neuropath.outcome$AmyPos_bl)
cdr.neuropath.outcome$ptau_pos_bl <- factor(cdr.neuropath.outcome$ptau_pos_bl)
cdr.neuropath.outcome$fullcaa <- factor(cdr.neuropath.outcome$fullcaa)
cdr.neuropath.outcome$fulltdp43 <- factor(cdr.neuropath.outcome$fulltdp43)
cdr.neuropath.outcome$fulllewy <- factor(cdr.neuropath.outcome$fulllewy)
cdr.neuropath.outcome.bl <- subset(cdr.neuropath.outcome, new_time==0)
cdr.desc.table <- table1(cdr.neuropath.outcome.bl[c("DX_bl","AmyPos_bl","PTGENDER_bl","ptau_pos_bl","MMSE_bl","CDGLOBAL_bl","CDRSB_bl", "fullcaa", "fulltdp43", "fulllewy")], splitby = ~DX_bl)$Table1

mmse.neuropath.outcome$CDGLOBAL_bl <- factor(mmse.neuropath.outcome$CDGLOBAL_bl)
mmse.neuropath.outcome$PTGENDER_bl <- factor(mmse.neuropath.outcome$PTGENDER_bl)
mmse.neuropath.outcome$DX_bl <- factor(mmse.neuropath.outcome$DX_bl)
mmse.neuropath.outcome$AmyPos_bl <- factor(mmse.neuropath.outcome$AmyPos_bl)
mmse.neuropath.outcome$ptau_pos_bl <- factor(mmse.neuropath.outcome$ptau_pos_bl)
mmse.neuropath.outcome$fullcaa <- factor(mmse.neuropath.outcome$fullcaa)
mmse.neuropath.outcome$fulltdp43 <- factor(mmse.neuropath.outcome$fulltdp43)
mmse.neuropath.outcome$fulllewy <- factor(mmse.neuropath.outcome$fulllewy)
mmse.neuropath.outcome.bl <- subset(mmse.neuropath.outcome, new_time==0)
mmse.desc.table <- table1(mmse.neuropath.outcome.bl[c("DX_bl","AmyPos_bl","PTGENDER_bl","ptau_pos_bl","MMSE_bl","CDGLOBAL_bl", "fullcaa", "fulltdp43", "fulllewy")], splitby = ~DX_bl)$Table1

mpacc.neuropath.outcome$CDGLOBAL_bl <- factor(mpacc.neuropath.outcome$CDGLOBAL_bl)
mpacc.neuropath.outcome$PTGENDER_bl <- factor(mpacc.neuropath.outcome$PTGENDER_bl)
mpacc.neuropath.outcome$DX_bl <- factor(mpacc.neuropath.outcome$DX_bl)
mpacc.neuropath.outcome$AmyPos_bl <- factor(mpacc.neuropath.outcome$AmyPos_bl)
mpacc.neuropath.outcome$ptau_pos_bl <- factor(mpacc.neuropath.outcome$ptau_pos_bl)
mpacc.neuropath.outcome$fullcaa <- factor(mpacc.neuropath.outcome$fullcaa)
mpacc.neuropath.outcome$fulltdp43 <- factor(mpacc.neuropath.outcome$fulltdp43)
mpacc.neuropath.outcome$fulllewy <- factor(mpacc.neuropath.outcome$fulllewy)
mpacc.neuropath.outcome.bl <- subset(mpacc.neuropath.outcome, new_time==0)
mpacc.desc.table <- table1(mpacc.neuropath.outcome.bl[c("DX_bl","AmyPos_bl","PTGENDER_bl","ptau_pos_bl","MMSE_bl","CDGLOBAL_bl","mPACCtrailsB_bl", "fullcaa", "fulltdp43", "fulllewy")], splitby = ~DX_bl)$Table1


mem.neuropath.outcome$CDGLOBAL_bl <- factor(mem.neuropath.outcome$CDGLOBAL_bl)
mem.neuropath.outcome$PTGENDER_bl <- factor(mem.neuropath.outcome$PTGENDER_bl)
mem.neuropath.outcome$DX_bl <- factor(mem.neuropath.outcome$DX_bl)
mem.neuropath.outcome$AmyPos_bl <- factor(mem.neuropath.outcome$AmyPos_bl)
mem.neuropath.outcome$ptau_pos_bl <- factor(mem.neuropath.outcome$ptau_pos_bl)
mem.neuropath.outcome$fullcaa <- factor(mem.neuropath.outcome$fullcaa)
mem.neuropath.outcome$fulltdp43 <- factor(mem.neuropath.outcome$fulltdp43)
mem.neuropath.outcome$fulllewy <- factor(mem.neuropath.outcome$fulllewy)
mem.neuropath.outcome.bl <- subset(mem.neuropath.outcome, new_time==0)
mem.desc.table <- table1(mem.neuropath.outcome.bl[c("DX_bl","AmyPos_bl","PTGENDER_bl","ptau_pos_bl","MMSE_bl","CDGLOBAL_bl", "ADNI_MEM_bl", "ADNI_EF_bl", "fullcaa", "fulltdp43", "fulllewy")], splitby = ~DX_bl)$Table1






if(source.script) {



saveRDS(adas.desc.table, "/Users/adamgabriellang/Desktop/clinical_trial_sim/save_ptau_and_total/adas_descr.rds")
saveRDS(adas.desc.table.neuro, "/Users/adamgabriellang/Desktop/clinical_trial_sim/save_ptau_and_total/adas_descr_nero.rds")


if(source.script) { ########################


#MCI fitting
unenriched.mci.totall11.model<- SampleSizeSimulation(sim.data = mci.gen.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                          fcompare_str = "TOTAL11~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.total11.model    <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                         fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.earlyage.total11.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.tplus.total11.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.neuro.total11.model  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")



unenriched.mci.cdrsb.model<- SampleSizeSimulation(sim.data = mci.gen.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.cdrsb.model    <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.earlyage.cdrsb.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

mci.scen1.tplus.cdrsb.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.neuro.cdrsb.model  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

mci.scen1.neuro.cdrsb.model.1000  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")




unenriched.mci.cdrsb.model<- SampleSizeSimulation(sim.data = mci.gen.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.cdrsb.model    <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.earlyage.cdrsb.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

mci.scen1.tplus.cdrsb.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.neuro.cdrsb.model  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")


unenriched.mci.mem.model<- SampleSizeSimulation(sim.data = mci.gen.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "ADNI_MEM~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.mem.model    <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.earlyage.mem.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.tplus.mem.model  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                    fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.neuro.mem.model  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[2]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")


unenriched.mci.ef.model<- SampleSizeSimulation(sim.data = mci.gen.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                fcompare_str = "ADNI_EF~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
mci.scen1.ef.model              <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                               fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
mci.scen1.earlyage.ef.model     <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
mci.scen1.tplus.ef.model        <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
mci.scen1.neuro.ef.model        <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[2]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")




#AD fitting
unenriched.ad.totall11.model<- SampleSizeSimulation(sim.data = ad.gen.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "TOTAL11~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
ad.scen1.total11.model    <- SampleSizeSimulation(sim.data = ad.scen1.i.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
ad.scen1.earlyage.total11.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.earlyage.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                         fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")

ad.scen1.tplus.total11.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
ad.scen1.neuro.total11.model  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")



unenriched.ad.cdrsb.model<- SampleSizeSimulation(sim.data = ad.gen.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
ad.scen1.cdrsb.model    <- SampleSizeSimulation(sim.data = ad.scen1.i.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
ad.scen1.earlyage.cdrsb.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

ad.scen1.tplus.cdrsb.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
ad.scen1.neuro.cdrsb.model  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

ad.scen1.neuro.cdrsb.model.1000  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")




unenriched.ad.cdrsb.model<- SampleSizeSimulation(sim.data = ad.gen.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
ad.scen1.cdrsb.model    <- SampleSizeSimulation(sim.data = ad.scen1.i.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
ad.scen1.earlyage.cdrsb.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

ad.scen1.tplus.cdrsb.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
ad.scen1.neuro.cdrsb.model  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")


unenriched.ad.mem.model<- SampleSizeSimulation(sim.data = ad.gen.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                fcompare_str = "ADNI_MEM~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
ad.scen1.mem.model    <- SampleSizeSimulation(sim.data = ad.scen1.i.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                               fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
ad.scen1.earlyage.mem.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.earlyage.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
ad.scen1.tplus.mem.model  <- SampleSizeSimulation(sim.data = ad.scen1.i.tplus.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
ad.scen1.neuro.mem.model  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[2]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")


unenriched.ad.ef.model<- SampleSizeSimulation(sim.data = ad.gen.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                               fcompare_str = "ADNI_EF~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
ad.scen1.ef.model              <- SampleSizeSimulation(sim.data = ad.scen1.i.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
ad.scen1.earlyage.ef.model     <- SampleSizeSimulation(sim.data = ad.scen1.i.earlyage.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
ad.scen1.tplus.ef.model        <- SampleSizeSimulation(sim.data = ad.scen1.i.tplus.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
ad.scen1.neuro.ef.model        <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[2]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")


# Early AD fitting
unenriched.early.early.ad.totall11.model<- SampleSizeSimulation(sim.data = early.ad.gen.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "TOTAL11~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.total11.model    <- SampleSizeSimulation(sim.data = early.ad.scen1.i.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.earlyage.total11.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.earlyage.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                         fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.tplus.total11.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.neuro.total11.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")



unenriched.early.early.ad.cdrsb.model<- SampleSizeSimulation(sim.data = early.ad.gen.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
early.ad.scen1.cdrsb.model    <- SampleSizeSimulation(sim.data = early.ad.scen1.i.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
early.ad.scen1.earlyage.cdrsb.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

early.ad.scen1.tplus.cdrsb.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
early.ad.scen1.neuro.cdrsb.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

early.ad.scen1.neuro.cdrsb.model.1000  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                         fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")




unenriched.early.early.ad.cdrsb.model<- SampleSizeSimulation(sim.data = early.ad.gen.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
early.ad.scen1.cdrsb.model    <- SampleSizeSimulation(sim.data = early.ad.scen1.i.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
early.ad.scen1.earlyage.cdrsb.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

early.ad.scen1.tplus.cdrsb.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
early.ad.scen1.neuro.cdrsb.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")


unenriched.early.early.ad.mem.model<- SampleSizeSimulation(sim.data = early.ad.gen.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                               fcompare_str = "ADNI_MEM~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
early.ad.scen1.mem.model    <- SampleSizeSimulation(sim.data = early.ad.scen1.i.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                              fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
early.ad.scen1.earlyage.mem.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.earlyage.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
early.ad.scen1.tplus.mem.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.i.tplus.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
early.ad.scen1.neuro.mem.model  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[2]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")


unenriched.early.early.ad.ef.model   <- SampleSizeSimulation(sim.data = early.ad.gen.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                              fcompare_str = "ADNI_EF~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
early.ad.scen1.ef.model              <- SampleSizeSimulation(sim.data = early.ad.scen1.i.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
early.ad.scen1.earlyage.ef.model     <- SampleSizeSimulation(sim.data = early.ad.scen1.i.earlyage.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
early.ad.scen1.tplus.ef.model        <- SampleSizeSimulation(sim.data = early.ad.scen1.i.tplus.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
early.ad.scen1.neuro.ef.model        <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[2]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")

































} ####################################################
#if(source.script){
### rerun everything with 100 simulations

simrOptions(nsim=100)
unenriched.mci.total11.model100<- SampleSizeSimulation(sim.data = mci.scen1.generic.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "TOTAL11~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.total11.model100    <- SampleSizeSimulation(sim.data = mci.scen1.earlyage.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.earlyage.total11.model100  <- SampleSizeSimulation(sim.data = mci.scen1.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.tplus.total11.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
mci.scen1.neuro.total11.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")



mci.total11.list <- list(unenriched.mci.total11.model100,
                        mci.scen1.total11.model100,
                        mci.scen1.earlyage.total11.model100,
                        mci.scen1.tplus.total11.model100,
                        mci.scen1.neuro.total11.model100)


names(mci.total11.list) <- c(g1.mci, g2.mci, g3.mci, g4.mci, g5.mci)
total11.powerlines.mci<-CombineSimPlots(mci.total11.list, limits = seq(100, 1000, by=100))
total11.mci.dtm <- GroupDiseaseTraj(mci.total11.list, yaxislab_dpm = "ADAS-11")

save.all.list[["mci.total11"]] <- list("powerlines" = total11.powerlines.mci,
                                     "dpm"        = total11.mci.dtm)




unenriched.mci.cdr.model100<- SampleSizeSimulation(sim.data = mci.scen1.generic.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
mci.scen1.cdr.model100    <- SampleSizeSimulation(sim.data = mci.scen1.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
mci.scen1.earlyage.cdr.model100  <- SampleSizeSimulation(sim.data = mci.scen1.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
mci.scen1.tplus.cdr.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
mci.scen1.neuro.cdr.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")



mci.cdr.list <- list(unenriched.mci.cdr.model100,
                         mci.scen1.cdr.model100,
                         mci.scen1.earlyage.cdr.model100,
                         mci.scen1.tplus.cdr.model100,
                         mci.scen1.neuro.cdr.model100)


names(mci.cdr.list) <- c(g1.mci, g2.mci, g3.mci, g4.mci, g5.mci)
cdr.powerlines.mci<-CombineSimPlots(mci.cdr.list, limits = seq(100, 1000, by=100))
cdr.mci.dtm <- GroupDiseaseTraj(mci.cdr.list, yaxislab_dpm = "CDR-SOB")

save.all.list[["mci.cdr"]] <- list("powerlines" = cdr.powerlines.mci,
                                    "dpm"        = cdr.mci.dtm)





unenriched.mci.mmse.model100<- SampleSizeSimulation(sim.data = mci.scen1.generic.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "MMSE~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
mci.scen1.mmse.model100    <- SampleSizeSimulation(sim.data = mci.scen1.earlyage.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
mci.scen1.earlyage.mmse.model100  <- SampleSizeSimulation(sim.data = mci.scen1.tplus.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                             fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
mci.scen1.tplus.mmse.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
mci.scen1.neuro.mmse.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.tplus.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")

mci.mmse.list <- list(unenriched.mci.mmse.model100,
                     mci.scen1.mmse.model100,
                     mci.scen1.earlyage.mmse.model100,
                     mci.scen1.tplus.mmse.model100,
                     mci.scen1.neuro.mmse.model100)


names(mci.mmse.list) <- c(g1.mci, g2.mci, g3.mci, g4.mci, g5.mci)
mmse.powerlines.mci<-CombineSimPlots(mci.mmse.list, limits = seq(100, 1000, by=100))
mmse.mci.dtm <- GroupDiseaseTraj(mci.mmse.list, yaxislab_dpm = "MMSE")

save.all.list[["mci.mmse"]] <- list("powerlines" = mmse.powerlines.mci,
                                   "dpm"        = mmse.mci.dtm)






unenriched.mci.mpacc.model100<- SampleSizeSimulation(sim.data = mci.scen1.generic.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "mPACCtrailsB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.mpacc.model100    <- SampleSizeSimulation(sim.data = mci.scen1.earlyage.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.earlyage.mpacc.model100  <- SampleSizeSimulation(sim.data = mci.scen1.tplus.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.tplus.mpacc.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.neuro.mpacc.model100  <- SampleSizeSimulation(sim.data = mci.scen1.neur.tplus.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")

mci.mpacc.list <- list(unenriched.mci.mpacc.model100,
                      mci.scen1.mpacc.model100,
                      mci.scen1.earlyage.mpacc.model100,
                      mci.scen1.tplus.mpacc.model100,
                      mci.scen1.neuro.mpacc.model100)

names(mci.mpacc.list) <- c(g1.mci, g2.mci, g3.mci, g4.mci, g5.mci)
mpacc.powerlines.mci      <-CombineSimPlots(mci.mpacc.list, limits = seq(100, 1000, by=100))
mpacc.mci.dtm         <- GroupDiseaseTraj(mci.mpacc.list, yaxislab_dpm = "mPACCtrailsB")
save.all.list[["mci.mpacc"]] <- list("powerlines" = mpacc.powerlines.mci,
                                   "dpm"        = mpacc.mci.dtm)









unenriched.ad.total11.model100<- SampleSizeSimulation(sim.data = ad.scen1.generic.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "TOTAL11~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
ad.scen1.total11.model100    <- SampleSizeSimulation(sim.data = ad.scen1.earlyage.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
ad.scen1.earlyage.total11.model100  <- SampleSizeSimulation(sim.data = ad.scen1.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                             fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
ad.scen1.tplus.total11.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
ad.scen1.neuro.total11.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")



ad.total11.list <- list(unenriched.ad.total11.model100,
                         ad.scen1.total11.model100,
                         ad.scen1.earlyage.total11.model100,
                         ad.scen1.tplus.total11.model100,
                         ad.scen1.neuro.total11.model100)


names(ad.total11.list) <- c(g1.ad, g2.ad, g3.ad, g4.ad, g5.ad)
total11.powerlines.ad<-CombineSimPlots(ad.total11.list, limits = seq(100, 1000, by=100))
total11.ad.dtm <- GroupDiseaseTraj(ad.total11.list, yaxislab_dpm = "ADAS-11")

save.all.list[["ad.total11"]] <- list("powerlines" = total11.powerlines.ad,
                                     "dpm"        = total11.ad.dtm)




unenriched.ad.cdr.model100<- SampleSizeSimulation(sim.data = ad.scen1.generic.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
ad.scen1.cdr.model100    <- SampleSizeSimulation(sim.data = ad.scen1.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
ad.scen1.earlyage.cdr.model100  <- SampleSizeSimulation(sim.data = ad.scen1.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                         fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
ad.scen1.tplus.cdr.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
ad.scen1.neuro.cdr.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")





ad.cdr.list <- list(unenriched.ad.cdr.model100,
                        ad.scen1.cdr.model100,
                        ad.scen1.earlyage.cdr.model100,
                        ad.scen1.tplus.cdr.model100,
                        ad.scen1.neuro.cdr.model100)


names(ad.cdr.list) <- c(g1.ad, g2.ad, g3.ad, g4.ad, g5.ad)
cdr.powerlines.ad<-CombineSimPlots(ad.cdr.list, limits = seq(100, 1000, by=100))
cdr.ad.dtm <- GroupDiseaseTraj(ad.cdr.list, yaxislab_dpm = "CDR-SOB")
save.all.list[["ad.cdr"]] <- list("powerlines" = cdr.powerlines.ad,
                                      "dpm"        = cdr.ad.dtm)







unenriched.ad.mmse.model100<- SampleSizeSimulation(sim.data = ad.scen1.generic.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "MMSE~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
ad.scen1.mmse.model100    <- SampleSizeSimulation(sim.data = ad.scen1.earlyage.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
ad.scen1.earlyage.mmse.model100  <- SampleSizeSimulation(sim.data = ad.scen1.tplus.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
ad.scen1.tplus.mmse.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
ad.scen1.neuro.mmse.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.tplus.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")




ad.mmse.list <- list(unenriched.ad.mmse.model100,
                    ad.scen1.mmse.model100,
                    ad.scen1.earlyage.mmse.model100,
                    ad.scen1.tplus.mmse.model100,
                    ad.scen1.neuro.mmse.model100)


names(ad.mmse.list) <- c(g1.ad, g2.ad, g3.ad, g4.ad, g5.ad)
mmse.powerlines.ad<-CombineSimPlots(ad.mmse.list, limits = seq(100, 1000, by=100))
mmse.ad.dtm <- GroupDiseaseTraj(ad.mmse.list, yaxislab_dpm = "MMSE")
save.all.list[["ad.mmse"]] <- list("powerlines" = mmse.powerlines.ad,
                                  "dpm"        = mmse.ad.dtm)








unenriched.ad.mpacc.model100<- SampleSizeSimulation(sim.data = ad.scen1.generic.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "mPACCtrailsB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
ad.scen1.mpacc.model100    <- SampleSizeSimulation(sim.data = ad.scen1.earlyage.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
ad.scen1.earlyage.mpacc.model100  <- SampleSizeSimulation(sim.data = ad.scen1.tplus.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                           fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
ad.scen1.tplus.mpacc.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
ad.scen1.neuro.mpacc.model100  <- SampleSizeSimulation(sim.data = ad.scen1.neur.tplus.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")


ad.mpacc.list <- list(unenriched.ad.mpacc.model100,
                     ad.scen1.mpacc.model100,
                     ad.scen1.earlyage.mpacc.model100,
                     ad.scen1.tplus.mpacc.model100,
                     ad.scen1.neuro.mpacc.model100)
names(ad.mpacc.list) <- c(g1.ad, g2.ad, g3.ad, g4.ad, g5.ad)
mpacc.powerlines.ad<-CombineSimPlots(ad.mpacc.list, limits = seq(100, 1000, by=100))
mpacc.ad.dtm <- GroupDiseaseTraj(ad.mpacc.list, yaxislab_dpm = "mPACCtrailsB")
save.all.list[["ad.mpacc"]] <- list("powerlines" = mpacc.powerlines.ad,
                                   "dpm"        = mpacc.ad.dtm)



unenriched.early.ad.total11.model100<- SampleSizeSimulation(sim.data = early.ad.scen1.generic.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "TOTAL11~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.total11.model100    <- SampleSizeSimulation(sim.data = early.ad.scen1.earlyage.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.earlyage.total11.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                            fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.tplus.total11.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                         fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")
early.ad.scen1.neuro.total11.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.tplus.long[[1]], formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                                         fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-11")


early.ad.total11.list <- list(unenriched.early.ad.total11.model100,
                              early.ad.scen1.total11.model100,
                              early.ad.scen1.earlyage.total11.model100,
                              early.ad.scen1.tplus.total11.model100,
                              early.ad.scen1.neuro.total11.model100)


names(early.ad.total11.list) <- c(g1.early.ad, g2.early.ad, g3.early.ad, g4.early.ad, g5.early.ad)
total11.powerlines.early.ad  <- CombineSimPlots(early.ad.total11.list, limits = seq(100, 1000, by=100))
early.ad.total11.dtm <- GroupDiseaseTraj(early.ad.total11.list, yaxislab_dpm = "ADAS-11")
save.all.list[["early.ad.total11"]] <- list("powerlines" = total11.powerlines.early.ad,
                                         "dpm"        = early.ad.total11.dtm)









unenriched.early.ad.cdr.model100<- SampleSizeSimulation(sim.data = early.ad.scen1.generic.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
early.ad.scen1.cdr.model100    <- SampleSizeSimulation(sim.data = early.ad.scen1.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
early.ad.scen1.earlyage.cdr.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
early.ad.scen1.tplus.cdr.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")
early.ad.scen1.neuro.cdr.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                                                                           fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDR-SOB")


early.ad.cdr.list <- list(unenriched.early.ad.cdr.model100,
                              early.ad.scen1.cdr.model100,
                              early.ad.scen1.earlyage.cdr.model100,
                              early.ad.scen1.tplus.cdr.model100,
                              early.ad.scen1.neuro.cdr.model100)


names(early.ad.cdr.list) <- c(g1.early.ad, g2.early.ad, g3.early.ad, g4.early.ad, g5.early.ad)
cdr.powerlines.early.ad  <- CombineSimPlots(early.ad.cdr.list, limits = seq(100, 1000, by=100))
early.ad.cdr.dtm <- GroupDiseaseTraj(early.ad.cdr.list, yaxislab_dpm = "CDR-SOB")
save.all.list[["early.ad.cdr"]] <- list("powerlines" = cdr.powerlines.early.ad,
                                            "dpm"        = early.ad.cdr.dtm)



unenriched.early.ad.mmse.model100<- SampleSizeSimulation(sim.data = early.ad.scen1.generic.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "MMSE~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
early.ad.scen1.mmse.model100    <- SampleSizeSimulation(sim.data = early.ad.scen1.earlyage.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
early.ad.scen1.earlyage.mmse.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.tplus.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                           fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
early.ad.scen1.tplus.mmse.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")
early.ad.scen1.neuro.mmse.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.tplus.long[[1]], formula = "MMSE ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "MMSE~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "MMSE")



early.ad.mmse.list <- list(unenriched.early.ad.mmse.model100,
                          early.ad.scen1.mmse.model100,
                          early.ad.scen1.earlyage.mmse.model100,
                          early.ad.scen1.tplus.mmse.model100,
                          early.ad.scen1.neuro.mmse.model100)


names(early.ad.mmse.list) <- c(g1.early.ad, g2.early.ad, g3.early.ad, g4.early.ad, g5.early.ad)
mmse.powerlines.early.ad<-CombineSimPlots(early.ad.mmse.list, limits = seq(100, 1000, by=100))
early.ad.mmse.dtm <- GroupDiseaseTraj(early.ad.mmse.list, yaxislab_dpm = "MMSE")
save.all.list[["early.ad.mmse"]] <- list("powerlines" = mmse.powerlines.early.ad,
                                        "dpm"        = early.ad.mmse.dtm)





unenriched.early.ad.mpacc.model100<- SampleSizeSimulation(sim.data = early.ad.scen1.generic.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "mPACCtrailsB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
early.ad.scen1.mpacc.model100    <- SampleSizeSimulation(sim.data = early.ad.scen1.earlyage.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
early.ad.scen1.earlyage.mpacc.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.tplus.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
early.ad.scen1.tplus.mpacc.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
early.ad.scen1.neuro.mpacc.model100  <- SampleSizeSimulation(sim.data = early.ad.scen1.neur.tplus.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                       fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")


early.ad.mpacc.list <- list(unenriched.early.ad.mpacc.model100,
                           early.ad.scen1.mpacc.model100,
                           early.ad.scen1.earlyage.mpacc.model100,
                           early.ad.scen1.tplus.mpacc.model100,
                           early.ad.scen1.neuro.mpacc.model100)


names(early.ad.mpacc.list)  <- c(g1.early.ad, g2.early.ad, g3.early.ad, g4.early.ad, g5.early.ad)
mpacc.powerlines.early.ad   <- CombineSimPlots(early.ad.mpacc.list, limits = seq(100, 1000, by=100))
early.ad.mpacc.dtm <- GroupDiseaseTraj(early.ad.mpacc.list, yaxislab_dpm = "mPACCtrailsB")
save.all.list[["early.ad.mpacc"]] <- list("powerlines" = mpacc.powerlines.early.ad,
                                        "dpm"        = early.ad.mpacc.dtm)



saveRDS(save.all.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/power_curves_dpm.rds")




if(source.script) {
  mmse.powerlines.mci <- CombineSimPlots(list( "All MCI A+" = unenriched.mci.mmse.model100,
                                               "MCI A+ CDR >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.mmse.model100,
                                               "Early-Age" = mci.scen1.earlyage.mmse.model100,
                                               "Tau-Pos" = mci.scen1.tplus.mmse.model100,
                                               "No Copathologies"= mci.scen1.neuro.mmse.model100), limits = seq(100, 1000, by=100))
  
  
  mmse.powerlines.ad <- CombineSimPlots(list( "All MCI A+" = unenriched.ad.mmse.model100,
                                              "MCI A+ CDR >= 0.5 \n MMSE 24-30 Age 50-85" = ad.scen1.mmse.model100,
                                              "Early-Age" = ad.scen1.earlyage.mmse.model100,
                                              "Tau-Pos" = ad.scen1.tplus.mmse.model100,
                                              "No Copathologies"= ad.scen1.neuro.mmse.model100), limits = seq(100, 1000, by=100))
  
  
  
  
  unenriched.mci.mmse.model100$disease_progression_plot
  
  
  
  
  
  
  mmse.powerlines.mci$plot
  mmse.powerlines.ad$plot

unenriched.mci.cdrsb.model500<- SampleSizeSimulation(sim.data = mci.gen.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "CDRSB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.cdrsb.model500    <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                 fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.earlyage.cdrsb.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.earlyage500.cdrsb.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                           fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

mci.scen1.tplus.cdrsb.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")
mci.scen1.neuro.cdrsb.model500  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")

mci.scen1.neuro.cdrsb.model500.1000  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                          fcompare_str = "CDRSB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "CDRSB")



unenriched.mci.mpacc.model500<- SampleSizeSimulation(sim.data = mci.gen.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "mPACCtrailsB~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.mpacc.model500    <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                    fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.earlyage.mpacc.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                           fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.earlyage500.mpacc.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                              fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")

mci.scen1.tplus.mpacc.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")
mci.scen1.neuro.mpacc.model500  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[1]], formula = "mPACCtrailsB ~ treat + new_time+ new_time*treat + (1|RID)",
                                                        fcompare_str = "mPACCtrailsB~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "mPACCtrailsB")





unenriched.mci.mem.model500      <- SampleSizeSimulation(sim.data = mci.gen.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                fcompare_str = "ADNI_MEM~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.mem.model500           <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                               fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.earlyage.mem.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                      fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.tplus.mem.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[3]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")
mci.scen1.neuro.mem.model500  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[2]], formula = "ADNI_MEM ~ treat + new_time+ new_time*treat + (1|RID)",
                                                   fcompare_str = "ADNI_MEM~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_MEM")



unenriched.mci.ef.model500<- SampleSizeSimulation(sim.data = mci.gen.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                               fcompare_str = "ADNI_EF~new_time", breaks = seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")

mci.scen1.ef.model500    <- SampleSizeSimulation(sim.data = mci.scen1.i.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                              fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")
mci.scen1.earlyage.ef.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.earlyage.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                     fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")

mci.scen1.tplus.ef.model500  <- SampleSizeSimulation(sim.data = mci.scen1.i.tplus.long[[3]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")

mci.scen1.neuro.ef.model500  <- SampleSizeSimulation(sim.data = mci.scen1.neur.long[[2]], formula = "ADNI_EF ~ treat + new_time+ new_time*treat + (1|RID)",
                                                  fcompare_str = "ADNI_EF~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADNI_EF")




total11.powerlines<-CombineSimPlots(list("All MCI A+" = unenriched.mci.totall11.model500,
                                     "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.total11.model500,
                                     "Early-Age" = mci.scen1.earlyage.total11.model500,
                                     "Tau-Pos" = mci.scen1.tplus.total11.model500,
                                     "No Copathologies"= mci.scen1.neuro.total11.model500), limits = seq(100, 1000, by=100))

cdr.powerlines <-CombineSimPlots(list("All MCI A+" = unenriched.mci.cdrsb.model500,
                                "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.cdrsb.model500,
                                "Early-Age" = mci.scen1.earlyage.cdrsb.model500,
                                "Tau-Pos" = mci.scen1.tplus.cdrsb.model500,
                                "No Copathologies"= mci.scen1.neuro.cdrsb.model500), limits = seq(100, 1000, by=100))

mpacc.powerlines <- CombineSimPlots(list("All MCI A+" = unenriched.mci.mpacc.model500,
                                         "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.mpacc.model500,
                                         "Early-Age" = mci.scen1.earlyage.mpacc.model500,
                                         "Tau-Pos" = mci.scen1.tplus.mpacc.model500,
                                         "No Copathologies"= mci.scen1.neuro.mpacc.model500), limits = seq(100, 1000, by=100))

mem.powerlines  <- CombineSimPlots(list("All MCI A+" = unenriched.mci.mem.model500,
                                     "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.mem.model500,
                                     "Early-Age" = mci.scen1.earlyage.mpacc.model500,
                                     "Tau-Pos" = mci.scen1.tplus.mem.model500,
                                     "No Copathologies"= mci.scen1.neuro.mem.model500), limits = seq(100, 1000, by=100))

ef.powerlines  <- CombineSimPlots(list( "All MCI A+" = unenriched.mci.ef.model500,
                                        "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.ef.model500,
                                        "Early-Age" = mci.scen1.earlyage.ef.model500,
                                        "Tau-Pos" = mci.scen1.tplus.ef.model500,
                                        "No Copathologies"= mci.scen1.neuro.ef.model500), limits = seq(100, 1000, by=100))




save_list <- list("total11" = list(total11.powerlines,
                                   list("All MCI A+" = unenriched.mci.totall11.model500,
                                        "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.total11.model500,
                                        "Early-Age" = mci.scen1.earlyage.total11.model500,
                                        "Tau-Pos" = mci.scen1.tplus.total11.model500,
                                        "No Copathologies"= mci.scen1.neuro.total11.model500)),
                 "cdr"=list(cdr.powerlines,
                            list("All MCI A+" = unenriched.mci.cdrsb.model500,
                                 "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.cdrsb.model500,
                                 "Early-Age" = mci.scen1.earlyage.cdrsb.model500,
                                 "Tau-Pos" = mci.scen1.tplus.cdrsb.model500,
                                 "No Copathologies"= mci.scen1.neuro.cdrsb.model500)),
                 "mpacc" = list(mpacc.powerlines, list("All MCI A+" = unenriched.mci.mpacc.model500,
                                                       "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.mpacc.model500,
                                                       "Early-Age" = mci.scen1.earlyage.mpacc.model500,
                                                       "Tau-Pos" = mci.scen1.tplus.mpacc.model500,
                                                       "No Copathologies"= mci.scen1.neuro.mpacc.model500)),
                 "mem" = list(mem.powerlines, list("All MCI A+" = unenriched.mci.mem.model500,
                                                   "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.mem.model500,
                                                   "Early-Age" = mci.scen1.earlyage.mpacc.model500,
                                                   "Tau-Pos" = mci.scen1.tplus.mem.model500,
                                                   "No Copathologies"= mci.scen1.neuro.mem.model500)),
                 "ef" = list(ef.powerlines, list( "All MCI A+" = unenriched.mci.ef.model500,
                                                  "MCI A+ CDRSB >= 0.5 \n MMSE 24-30 Age 50-85" = mci.scen1.ef.model500,
                                                  "Early-Age" = mci.scen1.earlyage.ef.model500,
                                                  "Tau-Pos" = mci.scen1.tplus.ef.model500,
                                                  "No Copathologies"= mci.scen1.neuro.ef.model500)))



saveRDS(save_list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/powerlines.rds")



rm(check)
























test.out <-  SampleSizeSimulation(sim.data = data.mci1, formula = "TOTAL11 ~ treat + new_time+ new_time*treat + (1|RID)",
                                    fcompare_str = "TOTAL11~new_time", breaks =  seq(100, 1000, by=100), yaxislab_dpm = "ADAS-COG11")






mci1.plot.total11 <- CombineSimPlots(nonenrich = unenriched.totall1, enrich = data.mci1.totall11, limits = seq(100, 1000, by=100))

unenriched.cdr <- SampleSizeSimulation(sim.data = unenriched.data, formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                          fcompare_str = "CDRSB~new_time", breaks = seq(100, 800, by=100))
data.mci1.cdr <- SampleSizeSimulation(sim.data = data.mci1, formula = "CDRSB ~ treat + new_time+ new_time*treat + (1|RID)",
                                          fcompare_str = "CDRSB~new_time", breaks = seq(100, 800, by=100))



mci1.plot.cdr <- CombineSimPlots(nonenrich = unenriched.cdr, enrich = data.mci1.cdr, limits = seq(100, 800, by=100))

firstplot <- mci1.plot.total11$plot
firstplot <- firstplot + labs(title = "ADAS-COG11") + geom_hline(yintercept = 80, linetype="dashed")
firstplot

secondplot <- mci1.plot.cdr$plot
secondplot <- secondplot + labs(title = "CDRSB") +geom_hline(yintercept = 80, linetype="dashed")
secondplot





amci1.analyzed    <- PlotObsData(data.mci1, formula.fixed = "TOTAL11 ~ new_time", ylab= "ADAS-COG11")
mci1.analyzed
mci1.analyzed.cdr <- PlotObsData(data.mci1, formula.fixed = "CDRSB ~ new_time", ylab= "CDRSB")
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
}
}