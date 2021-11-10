adni.fulldata     <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/Data/adni_fulldata.csv")

#seperate dataframes by outcome
necc.cols.t11     <- c("RID", "M", "ADAS11","CAAPos", "LewyPos", "TDP43Pos")
necc.cols.t13     <- c("RID", "M", "ADAS13","CAAPos", "LewyPos", "TDP43Pos")
necc.cols.cdr     <- c("RID", "M", "CDRSB","CAAPos", "LewyPos", "TDP43Pos")
necc.cols.mmse    <- c("RID", "M", "MMSE","CAAPos", "LewyPos", "TDP43Pos")
necc.cols.mpacc   <- c("RID", "M", "mPACCtrailsB","CAAPos", "LewyPos", "TDP43Pos")
necc.cols.image   <- c("RID", "M", "ST29SV_harmonized_icv_adj", "ST88SV_harmonized_icv_adj","CAAPos", "LewyPos", "TDP43Pos")

adni.fulldata$RID            <- factor(adni.fulldata$RID)
adni.fulldata$LewyPos        <- factor(adni.fulldata$LewyPos)
adni.fulldata$TDP43Pos       <- factor(adni.fulldata$TDP43Pos)
adni.fulldata$CAAPos         <- factor(adni.fulldata$CAAPos)
adni.fulldata$AmyloidPos     <- factor(adni.fulldata$AmyloidPos)
adni.fulldata$TauPos         <- factor(adni.fulldata$TauPos)
adni.fulldata$PTGENDER       <- factor(adni.fulldata$PTGENDER)
adni.fulldata$APOE4          <- factor(adni.fulldata$APOE4)

combined.11    <- adni.fulldata[complete.cases(adni.fulldata[,necc.cols.t11]),]
combined.13    <- adni.fulldata[complete.cases(adni.fulldata[,necc.cols.t13]),]
combined.cdr   <- adni.fulldata[complete.cases(adni.fulldata[,necc.cols.cdr]),]
combined.mmse  <- adni.fulldata[complete.cases(adni.fulldata[,necc.cols.mmse]),]
combined.mpacc <- adni.fulldata[complete.cases(adni.fulldata[,necc.cols.mpacc]),]
combined.image <- adni.fulldata[complete.cases(adni.fulldata[,necc.cols.image]),]

full.data.list <- list("ADAS11"      = combined.11,
                       "ADAS13"      = combined.13,
                       "CDRSB"       = combined.cdr,
                       "MMSE"        = combined.mmse,
                       "MPACC"       = combined.mpacc,
                       "Imaging"     = combined.image) 


#subsetting to first two years per subject
full.data.list <- Map(QuickAdjust, full.data.list)

#Enriching into groups

#Within MCI
mci.scen1.generic                  <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyloidPos_bl==1 & AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30))
mci.scen1.generic.tplus            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyloidPos_bl==1 & TauPos_bl==1 &  AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30))
mci.scen1.earlyage                 <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyloidPos_bl==1 & CDGLOBAL_bl == 0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & AGE_bl >=50 & AGE_bl <= 65))
mci.scen1.no.path.tplus            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyloidPos_bl==1 & TauPos_bl==1 & AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & LewyPos==0 & CAAPos==0 & TDP43Pos==0))

mci.scen1.generic.long             <- purrr::map2(mci.scen1.generic, full.data.list, PullLongData)
mci.scen1.generic.long.tplus       <- purrr::map2(mci.scen1.generic.tplus, full.data.list, PullLongData)
mci.scen1.earlyage.long            <- purrr::map2(mci.scen1.earlyage, full.data.list, PullLongData)
mci.scen1.no.path.tplus.long       <- purrr::map2(mci.scen1.no.path.tplus, full.data.list, PullLongData)

#Within Early AD
early.ad.scen1.generic             <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyloidPos_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.generic.tplus       <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyloidPos_bl==1 & TauPos_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.earlyage            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyloidPos_bl==1 & AGE_bl >=55 & AGE_bl <= 65 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.nopath.tplus        <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyloidPos_bl==1 & TauPos_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1) & LewyPos==0 & CAAPos==0 & TDP43Pos==0))

early.ad.scen1.generic.long        <- purrr::map2(early.ad.scen1.generic, full.data.list, PullLongData)
early.ad.scen1.generic.long.tplus  <- purrr::map2(early.ad.scen1.generic.tplus, full.data.list, PullLongData)
early.ad.scen1.earlyage.long       <- purrr::map2(early.ad.scen1.earlyage, full.data.list, PullLongData)
early.ad.scen1.nopath.tplus.long   <- purrr::map2(early.ad.scen1.nopath.tplus, full.data.list, PullLongData)

#Within AD
ad.scen1.generic                   <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyloidPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=65 & AGE_bl <= 85))
ad.scen1.generic.tplus             <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyloidPos_bl==1 & TauPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=65 & AGE_bl <= 85))
ad.scen1.earlyage                  <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyloidPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=55 & AGE_bl <= 65))
ad.scen1.nopath.tplus              <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyloidPos_bl==1 & TauPos_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=65 & AGE_bl <= 85 & LewyPos==0 & CAAPos==0 & TDP43Pos==0))

ad.scen1.generic.long              <- purrr::map2(ad.scen1.generic, full.data.list, PullLongData)
ad.scen1.generic.long.tplus        <- purrr::map2(ad.scen1.generic.tplus, full.data.list, PullLongData)
ad.scen1.earlyage.long             <- purrr::map2(ad.scen1.earlyage, full.data.list, PullLongData)
ad.scen1.nopath.tplus.long         <- purrr::map2(ad.scen1.nopath.tplus, full.data.list, PullLongData)


saveRDS(list("cs" = early.ad.scen1.generic,       "long" = early.ad.scen1.generic.long),       "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/earlyadcohorts.rds")
saveRDS(list("cs" = early.ad.scen1.generic.tplus, "long" = early.ad.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/earlyadcohorts_tplus.rds")
saveRDS(list("cs" = early.ad.scen1.nopath.tplus,  "long" = early.ad.scen1.nopath.tplus.long),  "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/earlyadcohorts_neuroenriched_tplus.rds")


saveRDS(list("cs" = ad.scen1.generic,       "long" = ad.scen1.generic.long),       "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/adcohorts.rds")
saveRDS(list("cs" = ad.scen1.generic.tplus, "long" = ad.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/adcohorts_tplus.rds")
saveRDS(list("cs" = ad.scen1.nopath.tplus,  "long" = ad.scen1.nopath.tplus.long),  "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/adcohorts_neuroenriched_tplus.rds")


saveRDS(list("cs" = mci.scen1.generic,       "long" = mci.scen1.generic.long),       "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/mcicohorts.rds")
saveRDS(list("cs" = mci.scen1.generic.tplus, "long" = mci.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/mcicohorts_tplus.rds")
saveRDS(list("cs" = mci.scen1.no.path.tplus, "long" = mci.scen1.no.path.tplus.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/mcicohorts_neuroenriched_tplus.rds")


