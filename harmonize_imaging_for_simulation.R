source("/Users/adamgabriellang/Desktop/ComGamPackage/ComBatPackage/helper_functions.R")
source("/Users/adamgabriellang/Desktop/ComGamPackage/ComBatPackage/main.R")

adni_imaging         <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/unharmonized_freesurfer_imaging.csv")
adni_imaging["DXCURREN"][which(is.na(adni_imaging$DXCURREN)), ]   <- 0
adni_imaging["DXCHANGE"][which(is.na(adni_imaging$DXCHANGE)), ]   <- 0
adni_imaging["DIAGNOSIS"][which(is.na(adni_imaging$DIAGNOSIS)), ] <- 0
adni_imaging$diagnosis_image <- rowSums2(as.matrix(adni_imaging[,c("DXCURREN", "DXCHANGE", "DIAGNOSIS")]))
adni_imaging$age_image       <- adni_imaging$AGE
adni_imaging$ptgender_image  <- adni_imaging$PTGENDER
adni_imaging$RID <- substr(adni_imaging$PTID, 7, 15)
adni_imaging$RID <- as.numeric(adni_imaging$RID)
adni_imaging$RID <- as.character(adni_imaging$RID)
adni_imaging$RID <- factor(adni_imaging$RID)
adni_imaging$VISCODE <- adni_imaging$VISCODE2
adni_imaging <- subset(adni_imaging, Phase != "ADNI1 3T")
adni_imaging$RID     <- factor(adni_imaging$RID)

harm.demogs <- adnimerge[,c("RID", "VISCODE", "PTGENDER", "COLPROT", "APOE4", "DX", "PTEDUCAT", "AGE", "M")]
full.data.set        <- merge(harm.demogs, adni_imaging, by=c("RID", "VISCODE"))

vol.ims <- c("ST103CV", "ST44CV",  "ST29SV",  "ST88SV",  
             "ST10CV",  "ST24CV",  "ST32CV",  "ST40CV",  
             "ST26CV",  "ST83CV",  "ST91CV",  "ST99CV",  
             "ST85CV",  "ST37SV",  "ST96SV",  "ST30SV",  
             "ST89SV",  "ST127SV", "ST9SV",  "ST21SV",  "ST80SV")

imaging_data_set     <- full.data.set[,c("RID", "VISCODE", vol.ims, "PTGENDER.x", "COLPROT", "APOE4", "DX", "PTEDUCAT", "AGE.x", "M")]
imaging_data_set     <- na.omit(imaging_data_set)
imaging_data_set$AGE <- round(imaging_data_set$AGE + (imaging_data_set$M / 12), 1)
imaging_data_set$STUDY <- imaging_data_set$COLPROT
imaging_data_set <- imaging_data_set[order(imaging_data_set$RID, imaging_data_set$M, decreasing = FALSE), ]
imaging_data_set <- imaging_data_set[-which(duplicated(imaging_data_set[,c("RID", "VISCODE")])),]
imaging_data_set$STUDY <- factor(imaging_data_set$STUDY)
imaging_data_set$Baseline <- !duplicated(imaging_data_set[,c("RID", "STUDY")])
imaging_data_set_baseline <- subset(imaging_data_set, Baseline==TRUE)
colnames(imaging_data_set) <- c("RID", "VISCODE", vol.ims, "PTGENDER", "STUDY", "APOE4", "DX", "PTEDUCAT", "AGE", "M")
colnames(imaging_data_set_baseline) <- c("RID", "VISCODE", vol.ims, "PTGENDER", "STUDY", "APOE4", "DX", "PTEDUCAT", "AGE", "M")


im.cols      <- vol.ims
drop         <- which(im.cols == "ST10CV")
im.cols      <- im.cols[-drop]
icv.col      <- as.data.frame(imaging_data_set_baseline$ST10CV)
colnames(icv.col) <- c("ST10CV")
cov.cols     <- c("APOE4", "PTEDUCAT", "AGE", "PTGENDER", "STUDY", "DX")
cov.cols.icv <- c("PTGENDER", "STUDY")



ims_preharm  <- imaging_data_set_baseline[,c(im.cols)]
covs_preharm <- imaging_data_set_baseline[,c(cov.cols)]
covs_icv_preharm <- imaging_data_set_baseline[,c(cov.cols.icv)] 

covs_preharm$APOE4 <- factor(covs_preharm$APOE4)
covs_preharm$PTGENDER <- factor(covs_preharm$PTGENDER)
covs_preharm$STUDY <- factor(covs_preharm$STUDY)

covs_icv_preharm$PTGENDER <- factor(covs_icv_preharm$PTGENDER)
covs_icv_preharm$STUDY <- factor(covs_icv_preharm$STUDY)
training <- FindIndices(covs_preharm, "DX", "CN")
covs_preharm$DX <- NULL



#### harmonize ICV
harmedICV <- ComGamHarm(feature.data = as.data.frame(icv.col),
                        covar.data = covs_icv_preharm,
                        training.indicies = as.numeric(training))
harm.icv.feat         <- as.data.frame(t(harmedICV$harm.results))
harm.icv.shift.params <- harmedICV$shift.scale.params
harm.icv.model <- harmedICV$models.list
covs_preharm <- cbind(covs_preharm, harm.icv.feat)
colnames(covs_preharm) <- c("APOE4", "PTEDUCAT", "AGE", "PTGENDER", "STUDY", "ICV_harmonized")

#### harmonize Feats
harmedFeatures <- ComGamHarm(feature.data = ims_preharm,
                             covar.data   = covs_preharm,
                             training.indicies = training,
                             smooth.terms      = c("AGE"),
                             k.val             = 5)
harmed_feats <- as.data.frame(t(harmedFeatures$harm.results))
harm.shift   <- harmedFeatures$shift.scale.params
harm.models  <- harmedFeatures$models.list


#harmonize ICV longitudinally
icv.long <- imaging_data_set$ST10CV
icv.long <- as.data.frame(icv.long)
cov.cols.icv.long <- imaging_data_set[ ,cov.cols.icv]
colnames(icv.long) <- "ST10CV"
long.icv.harm <- ApplyHarm(long.im     = icv.long,
                           long.cov    = cov.cols.icv.long,
                           site.params = harm.icv.shift.params,
                           model.list  = harm.icv.model)
long.icv.harm           <- as.data.frame(t(long.icv.harm))
rownames(long.icv.harm) <- 1:nrow(long.icv.harm)

#harmonize imaging features longitudinally 
features.long <- imaging_data_set[,im.cols]
cov.cols.long <- imaging_data_set[,cov.cols]
cov.cols.long <- cbind(cov.cols.long, long.icv.harm)
colnames(cov.cols.long) <- c(cov.cols, "ICV_harmonized")
cov.cols.long$DX <- NULL
long.feats.harm  <- ApplyHarm(long.im   = features.long,
                           long.cov     = cov.cols.long,
                           site.params  = harm.shift,
                           model.list   = harm.models)
long.feats.harm          <- as.data.frame(t(long.feats.harm))
harmed.imaging           <- cbind(long.feats.harm, long.icv.harm)
colnames(harmed.imaging) <- paste(colnames(harmed.imaging), "_harmonized", sep="")
harmed.imaging <- cbind(imaging_data_set[,c("RID", "VISCODE")], harmed.imaging)
ventric.region <- c("ST37SV_harmonized", "ST96SV_harmonized", "ST30SV_harmonized", "ST89SV_harmonized", 
                    "ST127SV_harmonized", "ST9SV_harmonized", "ST21SV_harmonized",
                    "ST80SV_harmonized")
left.meta  <- c("ST24CV_harmonized", "ST32CV_harmonized", "ST40CV_harmonized", "ST26CV_harmonized")
right.meta <- c("ST83CV_harmonized", "ST91CV_harmonized", "ST99CV_harmonized", "ST85CV_harmonized")




harmed.imaging$VentricalSum <- matrixStats::rowSums2(as.matrix(harmed.imaging[,ventric.region]))
harmed.imaging$LeftMeta     <- matrixStats::rowSums2(as.matrix(harmed.imaging[,left.meta]))
harmed.imaging$RightMeta    <- matrixStats::rowSums2(as.matrix(harmed.imaging[,right.meta]))
total.imaging <- cbind(imaging_data_set, harmed.imaging[,3:ncol(harmed.imaging)])


if(FALSE) {
write.csv(total.imaging,"/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/harmed_unharmed_freesurfer_imaging.csv")
}
  

