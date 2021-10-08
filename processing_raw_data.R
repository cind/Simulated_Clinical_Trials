########## Left Meta ROI #########

#Left entorhinal ST24CV
#Left inferior temporal ST32CV
#Left middle  temporal ST40CV
#Left fusiform         ST26CV

##################################

######### Right Meta ROI #########

#Right entorhinal ST83CV
#Right inferior temporal ST91CV
#Right middle  temporal ST99CV
#Right fusiform         ST85CV

##################################

######### Ventricular region #########

#LeftLateralVentricle (ST37); 
#RightLateralVentricle (ST96); 
#LeftInferiorLateralVentricle (ST30); 
#RightInferiorLateralVentricle (ST89); 
#ThirdVentricle (ST127); 
#FourthVentricle (ST9); 
#FifthVentrical (ST8); 
#LeftChoroidPlexus (ST21);
#RightChoroidPlexus (ST80)

##################################

######### Cutpoints #########

#CSFABETA cutpoint <= 980
#CSFPTAU cutpoint  >= 24
#PETAMY cutpooint  >= .78

##################################

vol.ims <- c("ST103CV", "ST44CV",  "ST29SV",  "ST88SV",  
             "ST10CV",  "ST24CV",  "ST32CV",  "ST40CV",  
             "ST26CV",  "ST83CV",  "ST91CV",  "ST99CV",  
             "ST85CV",  "ST37SV",  "ST96SV",  "ST30SV",  
             "ST89SV",  "ST127SV", "ST9SV",  "ST21SV",  "ST80SV")

#defaults to not writing data
#check if it looks okay first

writecsv <- FALSE

# reading all data
cdr_global        <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/CDR (1).csv")
adas_scores_1     <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADASSCORES.csv")
adas_scores_23    <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADAS_ADNIGO23.csv")
adni_neuropsych   <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Neuropsychological (1)/UWNPSYCHSUM_03_09_21.csv")
csf.upenn9        <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK9_04_19_17.csv")
csf.upenn10       <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK10_07_29_19.csv")
csf.upenn12       <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/Biospecimen_Results/UPENNBIOMK12_01_04_21.csv")
neuropath.data    <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/NEUROPATH_05_17_21.csv")
non.ad.imputation       <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADNI-nonADimputations-class.csv")
imaging_simulation_data <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/harmed_unharmed_freesurfer_imaging.csv")
av45                    <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/UCBERKELEYAV45_01_14_21.csv")



#CDR
cdr_global           <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/CDR (1).csv")
cdr_global$EXAMDATE  <- as.POSIXct(cdr_global$EXAMDATE, format="%Y-%M-%D")
cdr_global           <- cdr_global[order(cdr_global$RID, cdr_global$EXAMDATE, decreasing = FALSE), ]
cdr_global$VISCODE2  <- as.character(cdr_global$VISCODE2)
cdr_global           <- cdr_global[,c("RID", "VISCODE2", "CDGLOBAL")]
colnames(cdr_global) <- c("RID", "VISCODE", "CDGLOBAL")
cdr_global["VISCODE"][which(cdr_global$VISCODE=="sc"),] <- "bl"

#imaging

imaging_simulation_data$VISCODE.1 <- NULL
imaging_simulation_data$NA..1 <-imaging_simulation_data$NA..2 <- imaging_simulation_data$NA. <- NULL

#ADAS
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
adni_neuropsych <- adni_neuropsych[,c("RID", "VISCODE2", "ADNI_MEM", "ADNI_EF")]
colnames(adni_neuropsych)<- c("RID", "VISCODE", "ADNI_MEM", "ADNI_EF")


#CSF
csf.upenn9       <- csf.upenn9[ ,c("RID", "VISCODE2", "EXAMDATE",  "ABETA",   "TAU", "PTAU", "COMMENT")]
csf.upenn10      <- csf.upenn10[,c("RID", "VISCODE2", "DRAWDATE",  "ABETA42", "TAU", "PTAU", "NOTE")]
csf.upenn12      <- csf.upenn12[,c("RID", "VISCODE2", "EXAMDATE",  "ABETA",   "TAU", "PTAU", "NOTE")]

csf.upenn9$ABETA      <- as.numeric(as.character(csf.upenn9$ABETA))
csf.upenn9$PTAU       <- as.numeric(as.character(csf.upenn9$PTAU))
csf.upenn10$ABETA42   <- as.numeric(as.character(csf.upenn10$ABETA42))
csf.upenn10$PTAU      <- as.numeric(as.character(csf.upenn10$PTAU))
csf.upenn12$ABETA     <- as.numeric(as.character(csf.upenn12$ABETA))
csf.upenn12$PTAU      <- as.numeric(as.character(csf.upenn12$PTAU))
csf.upenn9$transform  <- rep(TRUE, nrow(csf.upenn9))
csf.upenn10$transform <- rep(TRUE, nrow(csf.upenn10))
csf.upenn12$transform <- rep(FALSE, nrow(csf.upenn12))
colnames(csf.upenn9)     <- colnames(csf.upenn10) <- colnames(csf.upenn12) <- c("RID", "VISCODE", "EXAMDATE", "ABETA", "TAU", "PTAU", "COMMENT", "TRANSFORM")
csf.data                 <- rbind(csf.upenn9, csf.upenn10, csf.upenn12)
csf.data$EXAMDATE        <- as.POSIXct(csf.data$EXAMDATE)
csf.data                 <- csf.data[order(csf.data$RID, csf.data$EXAMDATE, decreasing = FALSE), ]
csf.data$ABETA           <- as.numeric(as.character(csf.data$ABETA))
csf.data$PTAU            <- as.numeric(as.character(csf.data$PTAU))
csf.data$ABETA_untransformed <- csf.data$ABETA
csf.data$PTAU_untransformed  <- csf.data$PTAU
csf.data["ABETA"][which(csf.data$TRANSFORM == TRUE), ] <- (csf.data["ABETA"][which(csf.data$TRANSFORM == TRUE), ] * 1.014) + 29.25
csf.data["PTAU"][which(csf.data$TRANSFORM == TRUE), ] <- (csf.data["PTAU"][which(csf.data$TRANSFORM == TRUE), ] * .961) - .694
amyloid.cutpoint   <- (980 *1.014) + 29.25
ptau.cutpoint      <-  (24*.961) - .694
rownames(csf.data) <- 1:nrow(csf.data)
csf.drop <- which(duplicated(csf.data[,c("RID", "VISCODE")]))
csf.data <- csf.data[-csf.drop,]


#Neuropathology Data
non.ad.imputation["TDP43"][which(non.ad.imputation$TDP43 == TRUE),  ]  <- 1
non.ad.imputation["TDP43"][which(non.ad.imputation$TDP43 == FALSE), ]  <- 0
non.ad.imputation["LEWY"][which(non.ad.imputation$LEWY   == TRUE),  ]  <- 1
non.ad.imputation["LEWY"][which(non.ad.imputation$LEWY   == FALSE), ]  <- 0
non.ad.imputation["CAA"][which(non.ad.imputation$CAA     == TRUE),  ]  <- 1
non.ad.imputation["CAA"][which(non.ad.imputation$CAA     == FALSE), ]  <- 0



#PET Amyloid
av45.keeps           <- c("RID", "EXAMDATE", "VISCODE2", "SUMMARYSUVR_COMPOSITE_REFNORM")
av45                 <- av45[,av45.keeps]
colnames(av45)       <- c("RID", "EXAMDATE", "VISCODE", "SUMMARYSUVR_COMPOSITE_REFNORM")


#ADNI Merge
column.keeps             <- c("RID", "DX", "COLPROT", "ORIGPROT", 
                              "VISCODE", "EXAMDATE", "AGE", 
                              "PTGENDER", "PTEDUCAT", 
                              "APOE4", "CDRSB", "MMSE", 
                              "M", "mPACCtrailsB")

fulldata_demog           <- adnimerge[,column.keeps]



#Merging Data
all_outcomes_demog_merged <- merge(fulldata_demog, adas_scores, by=c("RID", "VISCODE"), all = TRUE)
all_outcomes_demog_merged <- merge(all_outcomes_demog_merged, adni_neuropsych, by=c("RID", "VISCODE"), all = TRUE)
all_outcomes_demog_merged <- merge(all_outcomes_demog_merged, csf.data, by=c("RID", "VISCODE"), all = TRUE)
all_outcomes_demog_merged <- merge(all_outcomes_demog_merged, imaging_simulation_data, by=c("RID", "VISCODE"), all = TRUE)

names(all_outcomes_demog_merged)[names(all_outcomes_demog_merged) == 'EXAMDATE.x'] <- 'EXAMDATE_adnimerge'
names(all_outcomes_demog_merged)[names(all_outcomes_demog_merged) == "EXAMDATE.y"] <- "EXAMDATE_csf"
all_outcomes_demog_merged <- merge(all_outcomes_demog_merged, neuropath.data, by="RID", all = TRUE)
all_outcomes_demog_merged <- merge(all_outcomes_demog_merged, non.ad.imputation, by="RID", all = TRUE)
all_outcomes_demog_merged <- merge(all_outcomes_demog_merged, cdr_global, by=c("RID", "VISCODE"), all = TRUE)
all_outcomes_demog_merged <- merge(all_outcomes_demog_merged, av45, by=c("RID", "VISCODE"), all = TRUE)
names(all_outcomes_demog_merged)[names(all_outcomes_demog_merged) == 'EXAMDATE'] <- 'EXAMDATE_pet'
all_outcomes_demog_merged$AGE <- round(all_outcomes_demog_merged$AGE.x + (all_outcomes_demog_merged$M.x / 12), 1)
all_outcomes_demog_merged$EXAMDATE_adnimerge <- as.POSIXct(all_outcomes_demog_merged$EXAMDATE_adnimerge)
all_outcomes_demog_merged$EXAMDATE_csf<- as.POSIXct(all_outcomes_demog_merged$EXAMDATE_csf)
all_outcomes_demog_merged$EXAMDATE_pet<- as.POSIXct(all_outcomes_demog_merged$EXAMDATE_pet)

all_outcomes_demog_merged$diff_csf <- abs(difftime(all_outcomes_demog_merged$EXAMDATE_adnimerge, all_outcomes_demog_merged$EXAMDATE_csf, units = "weeks"))
all_outcomes_demog_merged$diff_pet <- abs(difftime(all_outcomes_demog_merged$EXAMDATE_adnimerge, all_outcomes_demog_merged$EXAMDATE_pet, units = "weeks"))
all_outcomes_demog_merged$csf_valid <- all_outcomes_demog_merged$pet_valid <- rep(NA, nrow(all_outcomes_demog_merged))
all_outcomes_demog_merged["csf_valid"][which(all_outcomes_demog_merged$diff_csf / 52 <= .5), ] <- 1
all_outcomes_demog_merged["csf_valid"][which(all_outcomes_demog_merged$diff_csf / 52 >  .5), ] <- 0

all_outcomes_demog_merged["pet_valid"][which(all_outcomes_demog_merged$diff_pet / 52 <= .5), ] <- 1
all_outcomes_demog_merged["pet_valid"][which(all_outcomes_demog_merged$diff_pet / 52 >  .5), ] <- 0

all_outcomes_demog_merged$csf_pos <- all_outcomes_demog_merged$suvr_pos<- all_outcomes_demog_merged$ptau_pos<- rep(NA, nrow(all_outcomes_demog_merged))
all_outcomes_demog_merged$ABETA   <- as.numeric(as.character(all_outcomes_demog_merged$ABETA))

all_outcomes_demog_merged["suvr_pos"][which(all_outcomes_demog_merged$pet_valid == 1 & all_outcomes_demog_merged$SUMMARYSUVR_COMPOSITE_REFNORM >= .78), ] <- 1
all_outcomes_demog_merged["suvr_pos"][which(all_outcomes_demog_merged$pet_valid == 1 & all_outcomes_demog_merged$SUMMARYSUVR_COMPOSITE_REFNORM < .78), ] <- 0

all_outcomes_demog_merged["csf_pos"][which(all_outcomes_demog_merged$csf_valid  == 1 & all_outcomes_demog_merged$ABETA < amyloid.cutpoint), ] <- 1
all_outcomes_demog_merged["csf_pos"][which(all_outcomes_demog_merged$csf_valid  == 1 & all_outcomes_demog_merged$ABETA >= amyloid.cutpoint), ] <- 0

all_outcomes_demog_merged["ptau_pos"][which(all_outcomes_demog_merged$csf_valid == 1 & all_outcomes_demog_merged$PTAU >= ptau.cutpoint), ] <- 1
all_outcomes_demog_merged["ptau_pos"][which(all_outcomes_demog_merged$csf_valid == 1 & all_outcomes_demog_merged$PTAU < ptau.cutpoint), ] <- 0

all_outcomes_demog_merged <- all_outcomes_demog_merged[order(all_outcomes_demog_merged$RID, all_outcomes_demog_merged$VISCODE, decreasing = FALSE), ]
all_outcomes_demog_merged$AmyPos      <- rep(NA, nrow(all_outcomes_demog_merged))
all_outcomes_demog_merged$AmyPos_full <- rep(NA, nrow(all_outcomes_demog_merged))
all_outcomes_demog_merged$TauPos_full <- rep(NA, nrow(all_outcomes_demog_merged))

for(i in 1:nrow(all_outcomes_demog_merged)) {
  if(!is.na(all_outcomes_demog_merged["suvr_pos"][i,]) & all_outcomes_demog_merged["suvr_pos"][i,]==1) {
    all_outcomes_demog_merged["AmyPos"][i,] <- 1
  } else if(!is.na(all_outcomes_demog_merged["csf_pos"][i,]) & all_outcomes_demog_merged["csf_pos"][i,]==1) {
    all_outcomes_demog_merged["AmyPos"][i,] <- 1
  } else if(!is.na(all_outcomes_demog_merged["csf_pos"][i,]) | !is.na(all_outcomes_demog_merged["suvr_pos"][i,])) {
    all_outcomes_demog_merged["AmyPos"][i,] <- 0
  }
}



all_outcomes_demog_merged$M_vis <- substr(all_outcomes_demog_merged$VISCODE, 2, 10)
all_outcomes_demog_merged$M_vis <- as.numeric(all_outcomes_demog_merged$M_vis)
all_outcomes_demog_merged["M_vis"][which(all_outcomes_demog_merged$VISCODE=="bl"), ] <- 0
all_outcomes_demog_merged          <- all_outcomes_demog_merged[order(all_outcomes_demog_merged$RID, all_outcomes_demog_merged$M_vis, decreasing = FALSE), ]
all_outcomes_demog_merged$fulllewy <- all_outcomes_demog_merged$fulltdp43 <- all_outcomes_demog_merged$fullcaa <- rep(NA, nrow(all_outcomes_demog_merged)) 
baseline.var.list <- c("DX","AGE", "PTGENDER", "PTEDUCAT", 
                       "APOE4", "CDRSB", "MMSE", "mPACCtrailsB", "TOTAL11", "TOTAL13", 
                       "ABETA", "PTAU", "TAU","ABETA_untransformed", "PTAU_untransformed", "AmyPos","ptau_pos", "CDGLOBAL", "ADNI_MEM", "ADNI_EF")

all_outcomes_demog_merged$DX <- as.character(all_outcomes_demog_merged$DX.x)
all_outcomes_demog_merged <- all_outcomes_demog_merged[order(all_outcomes_demog_merged$RID, all_outcomes_demog_merged$M_vis, decreasing = FALSE), ]
all_outcomes_demog_merged$RID <- factor(as.character(all_outcomes_demog_merged$RID))

########################## RUNS SLOW ##########################

for(i in baseline.var.list) {
  all_outcomes_demog_merged <- CreateBaselineVar(all_outcomes_demog_merged, "M_vis", i)
}

##############################################################


baseline.var.list <- paste(baseline.var.list, "_bl", sep="")


all_outcomes_demog_merged$image_remerging <- 1:nrow(all_outcomes_demog_merged)
image_columns                             <- grep("_harmonized", 
                                                  colnames(all_outcomes_demog_merged), 
                                                  value = TRUE)
image_data_for_adj                        <- all_outcomes_demog_merged[,c("RID", "DX_bl", "AGE_bl", 
                                                        "PTGENDER.x", "AmyPos_bl", 
                                                        image_columns, "image_remerging")]
image_data_for_adj          <- na.omit(image_data_for_adj)

image_data_for_adj$Baseline <- !duplicated(image_data_for_adj$RID)

icv.training                <- subset(image_data_for_adj, 
                                      Baseline==TRUE & DX_bl=="CN" & AmyPos_bl==0)

image_columns               <- image_columns[!image_columns=="ST10CV_harmonized"]


for(i in image_columns) {
  correction_formula <- paste(i, "~ AGE_bl + PTGENDER.x + ST10CV_harmonized", sep=" ")
  corrected_name     <- paste(i, "icv_adj", sep="_")
  corrected_feature  <- feature.correction(training.data = icv.training,
                                           data          = image_data_for_adj,
                                           formula       = correction_formula,
                                           cr.feat1      = "ST10CV_harmonized",
                                           feat          = i)
  image_data_for_adj[corrected_name] <- corrected_feature
}

image_columns.adj           <- grep("_harmonized_icv_adj", colnames(image_data_for_adj), value = TRUE)
imaging_to_remerge          <- image_data_for_adj[,c(image_columns.adj, "image_remerging")]
all_outcomes_demog_merged   <- merge(all_outcomes_demog_merged, imaging_to_remerge, by="image_remerging", all.x = TRUE)


# WRITE OUTPUT OF MERGED FILE
if(writecsv) {
write.csv(all_outcomes_demog_merged, "/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/final_merged_data_4sim_remerged.csv")
}



