full.merged.data <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/final_merged_data_4sim.csv")
aut.rows <- which(complete.cases(full.merged.data[,c("NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                     "NPLBOD", "NPAMY")]))

imputed.rows <- which(complete.cases(full.merged.data[,c("TDP43", "LEWY", "CAA")]))
keeps <- unique(union(aut.rows, imputed.rows))
all.neuropat <- full.merged.data[keeps,]
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


strat.continuous <- replicate(7, c("PTEDUCAT_bl", "AGE_bl"), simplify = FALSE)

control.formulas <- list("ADAS11" = "TOTAL11_zscore      ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  (1|RID)",
                         "ADAS13" = "TOTAL13_zscore      ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  (1|RID)",
                         "CDRSB"  = "CDRSB_zscore        ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  (1|RID)",
                         "MMSE"   = "MMSE_zscore         ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  (1|RID)",
                         "MPACC"  = "mPACCtrailsB_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  (1|RID)",
                         "MEM"    = "ADNI_MEM_zscore     ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  (1|RID)",
                         "EF"     = "ADNI_EF_zscore      ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  (1|RID)")

enrichment.formulas <-  list(  "ADAS11"  = "TOTAL11_zscore_neuradj      ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                               "ADAS13"  = "TOTAL13_zscore_neuradj      ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                               "CDRSB"   = "CDRSB_zscore_neuradj        ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                               "MMSE"    = "MMSE_zscore_neuradj         ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                               "MPACC"   = "mPACCtrailsB_zscore_neuradj ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                               "MEM"     = "ADNI_MEM_zscore_neuradj     ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                               "EF"      = "ADNI_EF_zscore_neuradj      ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)")

simulation.formulas <-  list(  "ADAS11"  = "TOTAL11_zscore_neuradj      ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "ADAS13"  = "TOTAL13_zscore_neuradj      ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "CDRSB"   = "CDRSB_zscore_neuradj        ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "MMSE"    = "MMSE_zscore_neuradj         ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "MPACC"   = "mPACCtrailsB_zscore_neuradj ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "MEM"     = "ADNI_MEM_zscore_neuradj     ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "EF"      = "ADNI_EF_zscore_neuradj      ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time +  (1|RID)")

comparison.formulas <-  list(  "ADAS11"  = "TOTAL11_zscore_neuradj      ~ new_time +  PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "ADAS13"  = "TOTAL13_zscore_neuradj      ~ new_time +  PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "CDRSB"   = "CDRSB_zscore_neuradj        ~ new_time +  PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "MMSE"    = "MMSE_zscore_neuradj         ~ new_time +  PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "MPACC"   = "mPACCtrailsB_zscore_neuradj ~ new_time +  PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "MEM"     = "ADNI_MEM_zscore_neuradj     ~ new_time +  PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)",
                               "EF"      = "ADNI_EF_zscore_neuradj      ~ new_time +  PTEDUCAT_bl + AGE_bl + PTGENDER_bl +  fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time +  (1|RID)")

simulation.formulas.nopath <-  list("ADAS11"  = "TOTAL11_zscore      ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                                    "ADAS13"  = "TOTAL13_zscore      ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                                    "CDRSB"   = "CDRSB_zscore        ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                                    "MMSE"    = "MMSE_zscore         ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                                    "MPACC"   = "mPACCtrailsB        ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                                    "MEM"     = "ADNI_MEM_zscore     ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)",
                                    "EF"      = "ADNI_EF_zscore      ~ new_time + new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + (1|RID)")




boyle.es.incr <- replicate(3, c(0.592, 0.265, 0.347, 0.062, 0.018, 0.032), simplify = FALSE)
boyle.es.decr <- replicate(4, c(-0.592, -0.265, -0.347, -0.062, -0.018, -0.032), simplify = FALSE)
boyle.es <- do.call(c, list(boyle.es.incr, boyle.es.decr))
outcomes.neuro.regress <- list("TOTAL11_zscore", "TOTAL13_zscore",   "CDRSB_zscore", 
                               "MMSE_zscore", "mPACCtrailsB_zscore", "ADNI_MEM_zscore", "ADNI_EF_zscore")

subtract.control      <- replicate(7, "controlled", simplify = FALSE)
not.subtract.control  <- replicate(7, "not-subtracting-controls", simplify = FALSE)

# stratify continuous variables and assign treatment based on block randomization
treatment.enrichment.list <- list()
for(i in 1:length(baseline.enriched.data)) {
  sublist <- baseline.enriched.data[[i]]
  sublist.long <- full.enriched.data[[i]]
  enrich.data.stratcont <- map2(sublist, strat.continuous, StratifyContVar)
  full.treatment        <- map2(enrich.data.stratcont, sublist.long, RandomizeTreatment2, stratcolumns = c("PTGENDER", "fulllewy", "fullcaa", "fulltdp43", "PTEDUCAT_bl_strat", "AGE_bl_strat"))
  treatment.enrichment.list[[i]]        <- full.treatment
  names(treatment.enrichment.list[[i]]) <- names(full.data.list)
}
names(treatment.enrichment.list) <- names(baseline.enriched.data)

saveRDS(treatment.enrichment.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/stratification_example_again.rds")



# fit all models for control data
control.models.list <- list()
for(i in 1:length(full.control.data)) {
  sublist <- full.control.data[[i]]
  control.models.list[[i]]         <- map2(sublist, control.formulas, MapLmer)
  names(control.models.list[[i]])  <- names(control.formulas)
}
names(control.models.list) <- names(full.control.data)

saveRDS(control.models.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/control_models_example.rds")


# regress out neuropathology values to calculate mean rate accurately
neur.regr <- list()
for(i in 1:length(treatment.enrichment.list)) {
  sublist <- treatment.enrichment.list[[i]]
  neur.regr[[i]] <- pmap(list(sublist, boyle.es, outcomes.neuro.regress), RegressNeuro)
  names(neur.regr[[i]]) <- names(sublist)
}
names(neur.regr) <- names(treatment.enrichment.list)
treatment.enrichment.list <- neur.regr

saveRDS(neur.regr, "/Users/adamgabriellang/Desktop/clinical_trial_sim/neuro_regression_example.rds")

#fit all models for enriched data
enriched.models.list <- list()
for(i in 1:length(treatment.enrichment.list)) {
  sublist <- treatment.enrichment.list[[i]]
  enriched.models.list[[i]] <- map2(sublist, enrichment.formulas, MapLmer)
  names(enriched.models.list[[i]]) <- names(enrichment.formulas)
}
names(enriched.models.list) <- names(treatment.enrichment.list)
saveRDS(enriched.models.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/enriched_models_example.rds")


#calculate effect of mean decline -  healthy aging
adjusted.lmes <- list()
for(i in 1:length(treatment.enrichment.list)) {
  sublist <- treatment.enrichment.list[[i]]
  outcome.list <- list()
  for(j in 1:length(sublist)) {
    df <- sublist[[j]]
    newlmes <- RemoveNormalAging(control.model = control.models.list[[1]][[j]],
                      enriched.model = enriched.models.list[[i]][[j]],
                      formula.model = enrichment.formulas[[j]],
                      data = df)
    outcome.list[[j]] <- newlmes
     
  }
adjusted.lmes[[i]] <- outcome.list
names(adjusted.lmes[[i]]) <- names(enrichment.formulas)
}
names(adjusted.lmes) <- names(treatment.enrichment.list)
saveRDS(adjusted.lmes, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adjusted_lmes_example.rds")


#change effect sizes for neuropathologies
#boyle.adjusted.list <- list()
#for(i in 1:length(adjusted.lmes)) {
#  sublist <- adjusted.lmes[[i]]
#  boyle.adjusted.lmes <- map2(sublist, boyle.es, ChangeNeuroFixEf)
#  boyle.adjusted.list[[i]] <- boyle.adjusted.lmes
#  names(boyle.adjusted.list[[i]]) <- names(enrichment.formulas)
#}
#names(boyle.adjusted.list) <- names(full.enriched.data)

#create models for simulation
simulation.model.list <- list()
for(i in 1:length(adjusted.lmes)) {
  sublist <- adjusted.lmes[[i]]
  subdata <- treatment.enrichment.list[[i]] 
  simlist <- pmap(list(sublist, simulation.formulas, subdata, subtract.control, boyle.es), BuildSimulationModel)
  simulation.model.list[[i]] <- simlist
  names(simulation.model.list[[i]]) <- names(enrichment.formulas)
  }
names(simulation.model.list) <- names(full.enriched.data)
saveRDS(simulation.model.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/simulation_models_example.rds")


list_for_modeling <- list("treatment.enrichment.list"  = treatment.enrichment.list,
                          "simulation.formulas"        = simulation.formulas,
                          "enrichment.formulas"        = comparison.formulas,
                          "simulation.model.list"      = simulation.model.list)


saveRDS(list_for_modeling, "/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling_minus_control_2.rds")


#####################################################################  NO NEUROPATHOLOGIES!!!!!!!!!!!!! ################################################################# 
# stratify continuous variables and assign treatment based on block randomization
nopath.treatment.enrichment.list <- list()
for(i in 1:length(full.nopath.data.baseline)) {
  sublist <- full.nopath.data.baseline[[i]]
  sublist.long <- full.nopath.data.long[[i]]
  enrich.data.stratcont <- map2(sublist, strat.continuous, StratifyContVar)
  full.treatment        <- map2(enrich.data.stratcont, sublist.long, RandomizeTreatment2NoPath, stratcolumns = c("PTGENDER", "PTEDUCAT_bl_strat", "AGE_bl_strat"))
  nopath.treatment.enrichment.list[[i]]        <- full.treatment
  names(nopath.treatment.enrichment.list[[i]]) <- names(full.data.list)
}
names(nopath.treatment.enrichment.list) <- names(full.nopath.data.baseline)


#fit all models for enriched data
nopath.enriched.models.list <- list()
for(i in 1:length(nopath.treatment.enrichment.list)) {
  sublist <- nopath.treatment.enrichment.list[[i]]
  nopath.enriched.models.list[[i]] <- map2(sublist, control.formulas, MapLmer)
  names(nopath.enriched.models.list[[i]]) <- names(control.formulas)
}
names(nopath.enriched.models.list) <- names(nopath.treatment.enrichment.list)

#calculate effect of mean decline -  healthy aging
nopath.adjusted.lmes <- list()
for(i in 1:length(nopath.treatment.enrichment.list)) {
  sublist <- nopath.treatment.enrichment.list[[i]]
  outcome.list <- list()
  for(j in 1:length(sublist)) {
    df <- sublist[[j]]
    newlmes <- RemoveNormalAging(control.model  = control.models.list[[1]][[j]],
                                 enriched.model = nopath.enriched.models.list[[i]][[j]],
                                 formula.model  = control.formulas[[j]],
                                 data = df)
    outcome.list[[j]] <- newlmes
    
  }
  nopath.adjusted.lmes[[i]] <- outcome.list
  names(nopath.adjusted.lmes[[i]]) <- names(control.formulas)
}
names(nopath.adjusted.lmes) <- names(nopath.treatment.enrichment.list)


#create models for simulation
nopath.simulation.model.list <- list()
for(i in 1:length(nopath.adjusted.lmes)) {
  sublist <- nopath.adjusted.lmes[[i]]
  subdata <- nopath.treatment.enrichment.list[[i]] 
  simlist <- pmap(list(sublist, simulation.formulas.nopath, subdata, subtract.control), BuildSimulationModelNoPath)
  nopath.simulation.model.list[[i]] <- simlist
  names(nopath.simulation.model.list[[i]]) <- names(control.formulas)
}
names(nopath.simulation.model.list) <- names(nopath.treatment.enrichment.list)



list_for_modeling.nopath <- list("treatment.enrichment.list"         = nopath.treatment.enrichment.list,
                                 "simulation.formulas"               = simulation.formulas.nopath,
                                 "enrichment.formulas"               = control.formulas,
                                 "simulation.model.list"             = nopath.simulation.model.list)


saveRDS(list_for_modeling.nopath, "/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling_minus_control_nopath_2.rds")





#saveRDS(list_for_modeling, "/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling.rds")
#saveRDS(list_for_modeling, "/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling_minus_control.rds")
#saveRDS(control.models.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/control_models.rds")

objects.for.sim           <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling_minus_control.rds")
treatment.enrichment.list <- objects.for.sim$treatment.enrichment.list
simulation.formulas       <- objects.for.sim$simulation.formulas
simulation.model.list     <- objects.for.sim$simulation.model.list
enrichment.formulas       <- objects.for.sim$enrichment.formulas

# to simulate all models
simrOptions(nsim=10)
all.sims <- list()
for(i in 1:length(treatment.enrichment.list)) {
  sublist <- treatment.enrichment.list[[i]]
  sub.sim <- list()
  for(j in 1:length(sublist)) {
    df  <- sublist[[j]]
    print(paste(names(treatment.enrichment.list[i]), names(enrichment.formulas)[j], sep="_"))
    sim <- SampleSizeSimulation2(sim.data = df, 
                          formula      = simulation.formulas[[j]],
                          compare_str  = enrichment.formulas[[j]],
                          breaks       = seq(100, 1000, by=100), 
                          yaxislab_dpm = names(enrichment.formulas)[j],
                          model = simulation.model.list[[i]][[j]],
                          return_dpm = TRUE)
    sub.sim[[j]] <- sim
    
  }
  all.sims[[i]] <- sub.sim
  names(all.sims[[i]]) <- names(enrichment.formulas)
}
names(all.sims) <- names(treatment.enrichment.list)

allsimulationspt1 <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/allsimspt1.rds")
allsimulationspt2 <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/allsimspt2.rds")
allsimulationspt3 <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/allsimspt3.rds")
allsimulationspt4 <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/allsimspt4.rds")

allsimulationspt1 <- allsimulationspt1[1:2]
allsimulationspt2 <- allsimulationspt2[c(3,4)]
allsimulationspt3 <- allsimulationspt3[c(5,6)]
allsimulationspt4 <- allsimulationspt4[c(7,8,9)]

#allsimulations <- list()
allsimulations[[1]] <- allsimulationspt1[[1]]
allsimulations[[2]] <- allsimulationspt1[[2]]

allsimulations[[3]] <- allsimulationspt2[[1]]
allsimulations[[4]] <- allsimulationspt2[[2]]

allsimulations[[5]] <- allsimulationspt3[[1]]
allsimulations[[6]] <- allsimulationspt3[[2]]

allsimulations[[7]] <- allsimulationspt4[[1]]
allsimulations[[8]] <- allsimulationspt4[[2]]
allsimulations[[9]] <- allsimulationspt4[[3]]

names(allsimulations) <- names(full.enriched.data)
saveRDS(allsimulations, "/Users/adamgabriellang/Desktop/clinical_trial_sim/allsimulations.rds")
saveRDS(all.sims, "/Users/adamgabriellang/Desktop/clinical_trial_sim/allsimulations_dpm.rds")
testallsim2 <- readRDS("~/Desktop/clinical_trial_sim/allsimspt2_2.rds")
testallsim2<-testallsim2[4:6]
testallsim2[[1]][[7]]


saveRDS(list("cs"=early.ad.scen1.generic, "long"=early.ad.scen1.generic.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts.rds")
saveRDS(list("cs"=early.ad.scen1.generic.tplus, "long"=early.ad.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_tplus.rds")
saveRDS(list("cs"=early.ad.scen1.nopath, "long"=early.ad.scen1.nopath.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_neuroenriched.rds")
saveRDS(list("cs"=early.ad.scen1.nopath.tplus, "long"=early.ad.scen1.nopath.tplus.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_neuroenriched_tplus.rds")


saveRDS(list("cs"=ad.scen1.generic, "long"=ad.scen1.generic.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts.rds")
saveRDS(list("cs"=ad.scen1.generic.tplus, "long"=ad.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts_tplus.rds")
saveRDS(list("cs"=ad.scen1.nopath, "long"=ad.scen1.nopath.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts_neuroenriched.rds")
saveRDS(list("cs"=ad.scen1.nopath.tplus, "long"=ad.scen1.nopath.tplus.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/adcohorts_neuroenriched_tplus.rds")


saveRDS(list("cs"=mci.scen1.generic, "long"=mci.scen1.generic.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts.rds")
saveRDS(list("cs"=mci.scen1.generic.tplus, "long"=mci.scen1.generic.long.tplus), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts_tplus.rds")
saveRDS(list("cs"=mci.scen1.no.path, "long"=mci.scen1.no.path.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts_neuroenriched.rds")
saveRDS(list("cs"=mci.scen1.no.path.tplus, "long"=mci.scen1.no.path.tplus.long), "/Users/adamgabriellang/Desktop/clinical_trial_sim/mcicohorts_neuroenriched_tplus.rds")

