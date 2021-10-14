full.merged.data <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/final_merged_data_4sim.csv")
aut.rows <- which(complete.cases(full.merged.data[,c("NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                     "NPLBOD", "NPAMY")]))
imputed.rows <- which(complete.cases(full.merged.data[,c("TDP43", "LEWY", "CAA")]))
keeps <- unique(union(aut.rows, imputed.rows))
all.neuropat <- full.merged.data[keeps,]
all.neuropat <- SetNeuroData(all.neuropat)
#
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
all.neuropat$autopsy <- rep(NA, nrow(all.neuropat))
all.neuropat <- all.neuropat[complete.cases(all.neuropat[,c("DX_bl", "PTEDUCAT_bl", "PTGENDER_bl", "fulllewy", 
                                                            "fulltdp43", "fullcaa", "CDGLOBAL_bl", 
                                                            "AGE_bl", "MMSE_bl", "AmyPos_full_bl", "TauPos_full_bl")]), ]
all.neuropat.bl <- subset(all.neuropat, M_vis == 0)
all.neuropat.bl$RID <- factor(all.neuropat.bl$RID)

all.neuropat$autopsy <- complete.cases(all.neuropat[,c("NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                       "NPLBOD", "NPAMY")])
necc.cols.t11     <- c("RID", "M_vis", "TOTAL11")
necc.cols.t13     <- c("RID", "M_vis", "TOTAL13")
necc.cols.cdr     <- c("RID", "M_vis", "CDRSB")
necc.cols.mmse    <- c("RID", "M_vis", "MMSE")
necc.cols.mpacc   <- c("RID", "M_vis", "mPACCtrailsB")
necc.cols.mem     <- c("RID", "M_vis", "ADNI_MEM", "ADNI_EF")
necc.cols.image   <- c("RID", "M_vis", "ST29SV_harmonized_icv_adj", "ST88SV_harmonized_icv_adj")

all.neuropat$RID            <- factor(all.neuropat$RID)
all.neuropat$fulllewy       <- factor(all.neuropat$fulllewy)
all.neuropat$fulltdp43      <- factor(all.neuropat$fulltdp43)
all.neuropat$fullcaa        <- factor(all.neuropat$fullcaa)
all.neuropat$AmyPos_full_bl <- factor(all.neuropat$AmyPos_full_bl)
all.neuropat$TauPos_full_bl <- factor(all.neuropat$TauPos_full_bl)
all.neuropat$Amy_pos_path   <- factor(all.neuropat$Amy_pos_path)
all.neuropat$TAU_pos_path   <- factor(all.neuropat$TAU_pos_path)
all.neuropat$PTGENDER       <- factor(all.neuropat$PTGENDER)


combined.11    <- all.neuropat[complete.cases(all.neuropat[,necc.cols.t11]),]
combined.13    <- all.neuropat[complete.cases(all.neuropat[,necc.cols.t13]),]
combined.cdr   <- all.neuropat[complete.cases(all.neuropat[,necc.cols.cdr]),]
combined.mmse  <- all.neuropat[complete.cases(all.neuropat[,necc.cols.mmse]),]
combined.mpacc <- all.neuropat[complete.cases(all.neuropat[,necc.cols.mpacc]),]
combined.mem   <- all.neuropat[complete.cases(all.neuropat[,necc.cols.mem]),]
combined.image <- all.neuropat[complete.cases(all.neuropat[,necc.cols.image]),]

full.data.list <- list("ADAS11"      = combined.11,
                       "ADAS13"      = combined.13,
                       "CDRSB"       = combined.cdr,
                       "MMSE"        = combined.mmse,
                       "MPACC"       = combined.mpacc,
                       "MEM"         = combined.mem,
                       "EF"          = combined.mem,
                       "LH"          = combined.image,
                       "RH"          = combined.image) #use same dataset for ADNI EF and ADNI MEM

full.data.list <- Map(QuickAdjust, full.data.list)
control.cohort                     <- lapply(full.data.list, function(x)  subset(x, new_time==0 &  AGE_bl >= 65 & AGE_bl <= 85 & DX_bl=="CN" & AmyPos_full_bl==0 & TauPos_full_bl == 0 & CDGLOBAL_bl == 0 & fulllewy == 0 & fulltdp43 == 0 & fullcaa == 0 & MMSE_bl >=25))
control.cohort.long                <- purrr::map2(control.cohort, full.data.list, PullLongData)


full.data.list.zscored <- pmap(list(full.data.list, c("TOTAL11", "TOTAL13", "CDRSB","MMSE", "mPACCtrailsB", "ADNI_MEM", "ADNI_EF", "ST29SV_harmonized_icv_adj", "ST88SV_harmonized_icv_adj"), control.cohort.long), ZscoreAdj)
zscore.no.controls     <- pmap(list(full.data.list, c("TOTAL11", "TOTAL13", "CDRSB","MMSE", "mPACCtrailsB", "ADNI_MEM", "ADNI_EF", "ST29SV_harmonized_icv_adj", "ST88SV_harmonized_icv_adj"), full.data.list), ZscoreAdj)
control.cohort.long    <- pmap(list(control.cohort.long, c("TOTAL11", "TOTAL13", "CDRSB","MMSE", "mPACCtrailsB", "ADNI_MEM", "ADNI_EF", "ST29SV_harmonized_icv_adj", "ST88SV_harmonized_icv_adj"), control.cohort.long), ZscoreAdj)
full.data.list         <- full.data.list.zscored


mci.scen1.generic                  <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30))
mci.scen1.generic.tplus            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & TauPos_full_bl==1 &  AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30))
mci.scen1.earlyage                 <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & CDGLOBAL_bl == 0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & AGE >=50 & AGE <= 65))
mci.scen1.no.path                  <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & fulllewy==0 & fullcaa==0 & fulltdp43==0))
mci.scen1.no.path.tplus            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="MCI" & AmyPos_full_bl==1 & TauPos_full_bl==1 & AGE_bl >= 65 & AGE_bl <= 85 & CDGLOBAL_bl==0.5 & MMSE_bl >=24 & MMSE_bl <= 30 & fulllewy==0 & fullcaa==0 & fulltdp43==0))

mci.scen1.generic.long             <- purrr::map2(mci.scen1.generic, full.data.list, PullLongData)
mci.scen1.generic.long.tplus       <- purrr::map2(mci.scen1.generic.tplus, full.data.list, PullLongData)
mci.scen1.earlyage.long            <- purrr::map2(mci.scen1.earlyage, full.data.list, PullLongData)
mci.scen1.no.path.long             <- purrr::map2(mci.scen1.no.path, full.data.list, PullLongData)
mci.scen1.no.path.tplus.long       <- purrr::map2(mci.scen1.no.path.tplus, full.data.list, PullLongData)


early.ad.scen1.generic            <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.generic.tplus      <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & TauPos_full_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.earlyage           <- lapply(full.data.list, function(x)  subset(x,  new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & AGE_bl >=55 & AGE_bl <= 65 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1)))
early.ad.scen1.nopath             <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1) & fulllewy==0 & fullcaa==0 & fulltdp43==0))
early.ad.scen1.nopath.tplus       <- lapply(full.data.list, function(x)  subset(x, new_time==0 & (DX_bl=="Dementia" | DX_bl=="MCI") & AmyPos_full_bl==1 & TauPos_full_bl==1 & AGE_bl>=65 & AGE_bl <= 85 & MMSE_bl >=20 & MMSE_bl <= 28 & (CDGLOBAL_bl==.5 | CDGLOBAL_bl==1) & fulllewy==0 & fullcaa==0 & fulltdp43==0))

early.ad.scen1.generic.long        <- purrr::map2(early.ad.scen1.generic, full.data.list, PullLongData)
early.ad.scen1.generic.long.tplus  <- purrr::map2(early.ad.scen1.generic.tplus, full.data.list, PullLongData)
early.ad.scen1.earlyage.long       <- purrr::map2(early.ad.scen1.earlyage, full.data.list, PullLongData)
early.ad.scen1.nopath.long         <- purrr::map2(early.ad.scen1.nopath, full.data.list, PullLongData)
early.ad.scen1.nopath.tplus.long   <- purrr::map2(early.ad.scen1.nopath.tplus, full.data.list, PullLongData)

ad.scen1.generic                   <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=65 & AGE_bl <= 85))
ad.scen1.generic.tplus             <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & TauPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=65 & AGE_bl <= 85))
ad.scen1.earlyage                  <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=55 & AGE_bl <= 65))
ad.scen1.nopath                    <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=65 & AGE_bl <= 85 & fulllewy==0 & fullcaa==0 & fulltdp43==0))
ad.scen1.nopath.tplus              <- lapply(full.data.list, function(x)  subset(x, new_time==0 & DX_bl=="Dementia" & AmyPos_full_bl==1 & TauPos_full_bl==1 & CDGLOBAL_bl >= 1 & MMSE_bl >=20 & MMSE_bl <= 26 & AGE_bl >=65 & AGE_bl <= 85 & fulllewy==0 & fullcaa==0 & fulltdp43==0))

ad.scen1.generic.long              <- purrr::map2(ad.scen1.generic, full.data.list, PullLongData)
ad.scen1.generic.long.tplus        <- purrr::map2(ad.scen1.generic.tplus, full.data.list, PullLongData)
ad.scen1.earlyage.long             <- purrr::map2(ad.scen1.earlyage, full.data.list, PullLongData)
ad.scen1.nopath.long               <- purrr::map2(ad.scen1.nopath, full.data.list, PullLongData)
ad.scen1.nopath.tplus.long         <- purrr::map2(ad.scen1.nopath.tplus, full.data.list, PullLongData)

baseline.control.data <- list("control" = control.cohort)
full.control.data     <- list("control" = control.cohort.long)

baseline.enriched.data  <- list("mci"              = mci.scen1.generic,
                                "mci.tplus"        = mci.scen1.generic.tplus,
                                "mci.earlyage"     = mci.scen1.earlyage, 
                                "earlyad"          = early.ad.scen1.generic,
                                "earlyad.tplus"    = early.ad.scen1.generic.tplus,
                                "earlyad.earlyage" = early.ad.scen1.earlyage,
                                "ad"               = ad.scen1.generic,
                                "ad.tplus"         = ad.scen1.generic.tplus,
                                "ad.earlyage"      = ad.scen1.earlyage)

full.enriched.data <- list("mci"              = mci.scen1.generic.long,
                           "mci.tplus"        = mci.scen1.generic.long.tplus,
                           "mci.earlyage"     = mci.scen1.earlyage.long, 
                           "earlyad"          = early.ad.scen1.generic.long,
                           "earlyad.tplus"    = early.ad.scen1.generic.long.tplus,
                           "earlyad.earlyage" = early.ad.scen1.earlyage.long,
                           "ad"               = ad.scen1.generic.long,
                           "ad.tplus"         = ad.scen1.generic.long.tplus,
                           "ad.earlyage"      = ad.scen1.earlyage.long)


full.nopath.data.baseline   <- list("mci"     = mci.scen1.no.path,
                                    "earlyad" = early.ad.scen1.nopath,
                                    "ad"      = ad.scen1.nopath,
                                    "earlyad.tplus" = early.ad.scen1.nopath.tplus,
                                    "ad.tplus" = ad.scen1.nopath.tplus)

full.nopath.data.long   <- list("mci"     = mci.scen1.no.path.long,
                                "earlyad" = early.ad.scen1.nopath.long,
                                "ad"      = ad.scen1.nopath.long,
                                "earlyad.tplus" = early.ad.scen1.nopath.tplus.long,
                                "ad.tplus" = ad.scen1.nopath.tplus.long)


saveRDS(list("cs"=early.ad.scen1.generic,       "long" = early.ad.scen1.generic.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts.rds")
saveRDS(list("cs"=early.ad.scen1.generic.tplus, "long" = early.ad.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_tplus.rds")
saveRDS(list("cs"=early.ad.scen1.nopath,        "long" = early.ad.scen1.nopath.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_neuroenriched.rds")
saveRDS(list("cs"=early.ad.scen1.nopath.tplus,  "long" = early.ad.scen1.nopath.tplus.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_neuroenriched_tplus.rds")


saveRDS(list("cs"=ad.scen1.generic, "long"=ad.scen1.generic.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts.rds")
saveRDS(list("cs"=ad.scen1.generic.tplus, "long"=ad.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts_tplus.rds")
saveRDS(list("cs"=ad.scen1.nopath, "long"=ad.scen1.nopath.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts_neuroenriched.rds")
saveRDS(list("cs"=ad.scen1.nopath.tplus, "long"=ad.scen1.nopath.tplus.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts_neuroenriched_tplus.rds")


saveRDS(list("cs"=mci.scen1.generic, "long"=mci.scen1.generic.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts.rds")
saveRDS(list("cs"=mci.scen1.generic.tplus, "long"=mci.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts_tplus.rds")
saveRDS(list("cs"=mci.scen1.no.path, "long"=mci.scen1.no.path.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts_neuroenriched.rds")
saveRDS(list("cs"=mci.scen1.no.path.tplus, "long"=mci.scen1.no.path.tplus.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts_neuroenriched_tplus.rds")

