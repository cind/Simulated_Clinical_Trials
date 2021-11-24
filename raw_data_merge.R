#define amyloid and tau positive thresholds
#980 defined in upennbmk9 for amyloid
#24 defined in upennbmk9 for ptau

abeta.cutpoint     <- (980 *1.014) + 29.25
ptau.cutpoint      <-  (24*.961) - .694


#get all columns available from adnimerge dataframe
adni.mergecols <- adnimerge[,c("RID","COLPROT", "ORIGPROT","VISCODE", "EXAMDATE","DX",  
                               "AGE", "M", "PTGENDER", "PTEDUCAT", "APOE4", 
                               "ADAS11", "ADAS13","CDRSB", "MMSE", "mPACCtrailsB")]
#SUVR data
pet_amy.adni   <- ADNIMERGE::ucberkeleyav45
pet_amy.adni   <- pet_amy.adni[,c("RID", "VISCODE", "SUMMARYSUVR_COMPOSITE_REFNORM", "SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF")]

#CDR data
cdr.adni <- ADNIMERGE::cdr
cdr.adni <- cdr.adni[,c("RID", "VISCODE", "CDGLOBAL")]
cdr.adni["VISCODE"][which(cdr.adni$VISCODE == "sc"), ] <- "bl" #no bl level in cdr. bl and sc always a few days apart


#Merge
adni.fulldata <- merge(adni.mergecols,  pet_amy.adni, by = c("RID", "VISCODE"), all=TRUE)
adni.fulldata <- merge(adni.fulldata,   cdr.adni, by=c("RID", "VISCODE"), all=TRUE)

#CSF data
penn.9  <- ADNIMERGE::upennbiomk9
penn.10 <- ADNIMERGE::upennbiomk10
penn.12 <- ADNIMERGE::upennbiomk12_2020


#Scale upenn9 and upenn10 to match upenn12
penn.9       <- penn.9[, c("RID", "VISCODE", "ABETA", "PTAU")]
penn.9$TABLE <- "b9"
penn.9$ABETA <- as.numeric(penn.9$ABETA) #generates NA for >1700 pg/ml
penn.9$PTAU  <- as.numeric(penn.9$PTAU) #generates NA for < 3 pg/ml
penn.9$ABETA <- (penn.9$ABETA * 1.014) + 29.25
penn.9$PTAU  <- (penn.9$PTAU * .961)   - .694


penn.10           <- penn.10[, c("RID", "VISCODE", "ABETA42", "PTAU")]
penn.10$TABLE     <- "b10"
penn.10$ABETA42   <- as.numeric(penn.10$ABETA42)
penn.10$PTAU      <- as.numeric(penn.10$PTAU)
penn.10$ABETA42   <- (penn.10$ABETA42 * 1.014) + 29.25
penn.10$PTAU      <- (penn.10$PTAU * .961)   - .694
colnames(penn.10) <- c("RID", "VISCODE", "ABETA", "PTAU", "TABLE")

penn.12       <- penn.12[, c("RID", "VISCODE", "ABETA", "PTAU")]
penn.12$TABLE <- "b12"
penn.12$ABETA <- as.numeric(penn.12$ABETA)
penn.12$PTAU  <- as.numeric(penn.12$PTAU)

csf.adni           <- rbind(penn.9, penn.10, penn.12)

#initialize CSF Amyloid and Tau Positive columns
csf.adni$PTAU_pos  <- NA
csf.adni$ABETA_pos <- NA

#Create Amy Tau pos col based on CSF
csf.adni["PTAU_pos"][which(csf.adni$PTAU >= ptau.cutpoint), ]    <- 1
csf.adni["PTAU_pos"][which(csf.adni$PTAU < ptau.cutpoint), ]     <- 0
csf.adni["ABETA_pos"][which(csf.adni$ABETA <= abeta.cutpoint), ] <- 1
csf.adni["ABETA_pos"][which(csf.adni$ABETA > abeta.cutpoint), ]  <- 0

#merge
adni.fulldata <- merge(adni.fulldata, csf.adni, by=c("RID", "VISCODE"), all = TRUE)

## neuropath 
neuropath.data     <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/NEUROPATH_05_17_21.csv")
neuropath.data     <- SetNeuroData(neuropath.data)
neuropath.data     <- neuropath.data[,c("RID", "CAA_path", "Lewy_pos_path", "TDP_pos_path",  "Amy_pos_path",  "TAU_pos_path")]
non.ad.imputation  <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADNI-nonADimputations-class.csv")

#merge autopsy and imputed data
adni.fulldata <- merge(adni.fulldata, neuropath.data,    by = "RID", all= TRUE)
adni.fulldata <- merge(adni.fulldata, non.ad.imputation, by = "RID", all= TRUE)

#imaging
imaging_simulation_data <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/data_processed/harmed_and_unharmed_freesurfer_imaging.csv")
imaging_simulation_data <- imaging_simulation_data[,c(2,3, 9:ncol(imaging_simulation_data)-1)]
adni.fulldata           <- merge(adni.fulldata, imaging_simulation_data, by=c("RID", "VISCODE"), all = TRUE)
adni.fulldata$AGE       <- adni.fulldata$AGE + (adni.fulldata$M / 12)

#combining pos columns
adni.fulldata$AmyloidPos <- adni.fulldata$TauPos <-  adni.fulldata$LewyPos <-  adni.fulldata$TDP43Pos <- adni.fulldata$CAAPos <- NA

#Amyloid Order if subject has multiple: 1-Autopsy, 2-CSF, 3-SUVR
for(i in 1:nrow(adni.fulldata)) {
  if(!is.na(adni.fulldata["Amy_pos_path"][i,])) {
    adni.fulldata["AmyloidPos"][i,] <-  adni.fulldata["Amy_pos_path"][i,]
  } else if(!is.na(adni.fulldata["ABETA_pos"][i,])) {
    adni.fulldata["AmyloidPos"][i,] <- adni.fulldata["ABETA_pos"][i,]
  } else if(!is.na(adni.fulldata["SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF"][i,])) {
    adni.fulldata["AmyloidPos"][i,] <- adni.fulldata["SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF"][i,]
  } 
}

#TAU Order if subject has multiple: 1-Autopsy, 2-CSF
for(i in 1:nrow(adni.fulldata)) {
  if(!is.na(adni.fulldata["TAU_pos_path"][i,])) {
    adni.fulldata["TauPos"][i,] <- adni.fulldata["TAU_pos_path"][i,]
  } else if(!is.na(adni.fulldata["PTAU_pos"][i,])) {
    adni.fulldata["TauPos"][i,] <- adni.fulldata["PTAU_pos"][i,]
  } 
}


#Neurpath order 1-Autopsy, 2-Imputed
for(i in 1:nrow(adni.fulldata)) {
  if(!is.na(adni.fulldata["Lewy_pos_path"][i,])) {
    adni.fulldata["LewyPos"][i,] <- adni.fulldata["Lewy_pos_path"][i,]
  } else if(!is.na(adni.fulldata["LEWY"][i,])) {
    adni.fulldata["LewyPos"][i,] <- adni.fulldata["LEWY"][i,]
  } 
}

for(i in 1:nrow(adni.fulldata)) {
  if(!is.na(adni.fulldata["TDP_pos_path"][i,])) {
    adni.fulldata["TDP43Pos"][i,] <- adni.fulldata["TDP_pos_path"][i,]
  } else if(!is.na(adni.fulldata["TDP43"][i,])) {
    adni.fulldata["TDP43Pos"][i,] <- adni.fulldata["TDP43"][i,]
  } 
}

for(i in 1:nrow(adni.fulldata)) {
  if(!is.na(adni.fulldata["CAA_path"][i,])) {
    adni.fulldata["CAAPos"][i,] <- adni.fulldata["CAA_path"][i,]
  } else if(!is.na(adni.fulldata["CAA"][i,])) {
    adni.fulldata["CAAPos"][i,] <- adni.fulldata["CAA"][i,]
  } 
}



#Creating baseline variables from longitudinal variables
baseline.var.list <- c("DX", "AGE", "PTEDUCAT", "APOE4", 
                       "ADAS11", "ADAS13","CDRSB", "MMSE", "mPACCtrailsB",
                       "CDGLOBAL", "AmyloidPos", "TauPos")

for(i in baseline.var.list) {
  adni.fulldata <- CreateBaselineVar(adni.fulldata, "M", i)
}



#Adjust harmonized imaging columns for ICV
vol.ims <- c("ST103CV", "ST44CV",  "ST29SV",  "ST88SV",  
             "ST24CV",  "ST32CV",  "ST40CV",  "ST80SV",
             "ST26CV",  "ST83CV",  "ST91CV",  "ST99CV",  
             "ST85CV",  "ST37SV",  "ST96SV",  "ST30SV",  
             "ST89SV",  "ST127SV", "ST9SV",  "ST21SV")

vol.ims <- paste(vol.ims, "_harmonized", sep="")


image_data_for_adj <- adni.fulldata[,c("RID", "VISCODE", "DX_bl", "PTGENDER", "AmyloidPos_bl", "AGE_bl", vol.ims, "ST10CV_harmonized")]
image_data_for_adj <- image_data_for_adj[complete.cases(image_data_for_adj[,c(vol.ims)]),]
image_data_for_adj$Baseline <- !duplicated(image_data_for_adj$RID)

icv.training                <- subset(image_data_for_adj, 
                                      Baseline==TRUE & DX_bl==1 & AmyloidPos_bl==0)


for(i in vol.ims) {
  correction_formula <- paste(i, "~ AGE_bl + PTGENDER + ST10CV_harmonized", sep=" ")
  corrected_name     <- paste(i, "icv_adj", sep="_")
  corrected_feature  <- feature.correction(training.data = icv.training,
                                           data          = image_data_for_adj,
                                           formula       = correction_formula,
                                           cr.feat1      = "ST10CV_harmonized",
                                           feat          = i)
  image_data_for_adj[corrected_name] <- corrected_feature
}


vol.ims.icv.adj <- paste(vol.ims, "icv_adj", sep="_")
image_data_for_adj <- image_data_for_adj[,c("RID", "VISCODE", vol.ims.icv.adj)]

#merge imaging columns
adni.fulldata <- merge(adni.fulldata, image_data_for_adj, by = c("RID", "VISCODE"), all.x=TRUE)
adni.fulldata["DX_bl"][which(adni.fulldata$DX_bl == 1), ] <- "CN"
adni.fulldata["DX_bl"][which(adni.fulldata$DX_bl == 2), ] <- "MCI"
adni.fulldata["DX_bl"][which(adni.fulldata$DX_bl == 3), ] <- "Dementia"
write.csv(adni.fulldata, "/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/adni_fulldata.csv")





