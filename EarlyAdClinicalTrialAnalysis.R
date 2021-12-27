save.rds <- FALSE

earlyadcohorts                     <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/earlyadcohorts.rds")
earlyadcohorts.tplus               <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/earlyadcohorts_tplus.rds")
earlyadcohorts.neuroenriched.tplus <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/earlyadcohorts_neuroenriched_tplus.rds")

#keeping subjects with 3 or more time points
#update descriptive statistics and a few plots in powerpoint as this will remove a few subjects
earlyadcohorts                     <- Keep1YearorMore(earlyadcohorts)
earlyadcohorts.tplus               <- Keep1YearorMore(earlyadcohorts.tplus)
earlyadcohorts.neuroenriched.tplus <- Keep1YearorMore(earlyadcohorts.neuroenriched.tplus)

#Reading data

######################## ADAS13 ######################## 
cs.earlyad.adas13     <- earlyadcohorts$cs$ADAS13
long.earlyad.adas13   <-  earlyadcohorts$long$ADAS13

#Tau+
cs.earlyad.tplus.adas13    <- earlyadcohorts.tplus$cs$ADAS13
long.earlyad.tplus.adas13  <-  earlyadcohorts.tplus$long$ADAS13

# Neuroenriched Tau+
cs.earlyad.neuroenriched.tplus.adas13    <- earlyadcohorts.neuroenriched.tplus$cs$ADAS13
long.earlyad.neuroenriched.tplus.adas13  <-  earlyadcohorts.neuroenriched.tplus$long$ADAS13



######################### HIPPOCAMPUS ######################## 
cs.earlyad.hipp    <- earlyadcohorts$cs$Imaging
long.earlyad.hipp  <-  earlyadcohorts$long$Imaging

cs.earlyad.hipp$hipp_average    <- (cs.earlyad.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.hipp$hipp_average  <- (long.earlyad.hipp$ST29SV_harmonized_icv_adj + long.earlyad.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right

#Tau+
cs.earlyad.tplus.hipp           <- earlyadcohorts.tplus$cs$Imaging
long.earlyad.tplus.hipp         <-  earlyadcohorts.tplus$long$Imaging


cs.earlyad.tplus.hipp$hipp_average   <- (cs.earlyad.tplus.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.tplus.hipp$hipp_average <- (long.earlyad.tplus.hipp$ST29SV_harmonized_icv_adj + long.earlyad.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right

# Neuroenriched Tau+
cs.earlyad.neuroenriched.tplus.hipp    <- earlyadcohorts.neuroenriched.tplus$cs$Imaging
long.earlyad.neuroenriched.tplus.hipp <-  earlyadcohorts.neuroenriched.tplus$long$Imaging

cs.earlyad.neuroenriched.tplus.hipp$hipp_average   <- (cs.earlyad.neuroenriched.tplus.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.neuroenriched.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.neuroenriched.tplus.hipp$hipp_average <- (long.earlyad.neuroenriched.tplus.hipp$ST29SV_harmonized_icv_adj + long.earlyad.neuroenriched.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right



######################## MPACC######################## 
cs.earlyad.mpacc     <- earlyadcohorts$cs$MPACC
long.earlyad.mpacc   <-  earlyadcohorts$long$MPACC


#Tau+
cs.earlyad.tplus.mpacc    <- earlyadcohorts.tplus$cs$MPACC
long.earlyad.tplus.mpacc  <-  earlyadcohorts.tplus$long$MPACC


# Neuroenriched Tau+
cs.earlyad.neuroenriched.tplus.mpacc    <- earlyadcohorts.neuroenriched.tplus$cs$MPACC
long.earlyad.neuroenriched.tplus.mpacc  <-  earlyadcohorts.neuroenriched.tplus$long$MPACC




#Descriptive statistics tables
######################## ADAS13 ########################
earlyad.adas13.desc <- cs.earlyad.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                            "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                            "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                            "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]


earlyad.adas13.desc$TauPos_bl     <- factor(earlyad.adas13.desc$TauPos_bl)
earlyad.adas13.desc$AmyloidPos_bl <- factor(earlyad.adas13.desc$AmyloidPos_bl)
colnames(earlyad.adas13.desc)     <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                     "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                     "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                     "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.adas13.desc <- table1(earlyad.adas13.desc)
earlyad.adas13.desc <- earlyad.adas13.desc$Table1
View(earlyad.adas13.desc)
earlyad.adas13.desc <- earlyad.adas13.desc[c(1, 6:9, 11:nrow(earlyad.adas13.desc)),]


#Enriched for TAU
earlyad.tplus.adas13.desc <- cs.earlyad.tplus.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                        "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                        "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

earlyad.tplus.adas13.desc$TauPos_bl <- factor(earlyad.tplus.adas13.desc$TauPos_bl)
earlyad.tplus.adas13.desc$AmyloidPos_bl <- factor(earlyad.tplus.adas13.desc$AmyloidPos_bl)

colnames(earlyad.tplus.adas13.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                         "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                         "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                         "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.tplus.adas13.desc <- table1(earlyad.tplus.adas13.desc, splitby = ~ Diagnosis)
earlyad.tplus.adas13.desc <- earlyad.tplus.adas13.desc$Table1
earlyad.tplus.adas13.desc <- earlyad.tplus.adas13.desc[c(1, 6:9, 11:nrow(earlyad.tplus.adas13.desc)),]


#Neuroenriched Tau+
earlyad.neuroenriched.tplus.adas13.desc <- cs.earlyad.neuroenriched.tplus.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                    "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                    "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                                                    "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

earlyad.neuroenriched.tplus.adas13.desc$TauPos_bl <- factor(earlyad.neuroenriched.tplus.adas13.desc$TauPos_bl)
earlyad.neuroenriched.tplus.adas13.desc$AmyloidPos_bl <- factor(earlyad.neuroenriched.tplus.adas13.desc$AmyloidPos_bl)

colnames(earlyad.neuroenriched.tplus.adas13.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                       "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.tplus.adas13.desc <- table1(earlyad.neuroenriched.tplus.adas13.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.tplus.adas13.desc <- earlyad.neuroenriched.tplus.adas13.desc$Table1
earlyad.neuroenriched.tplus.adas13.desc <- earlyad.neuroenriched.tplus.adas13.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.tplus.adas13.desc)),]



#Descriptive statistics tables
######################## Hippocampus ########################
earlyad.hipp.desc <- cs.earlyad.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                            "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                            "hipp_average" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                            "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]


earlyad.hipp.desc$TauPos_bl <- factor(earlyad.hipp.desc$TauPos_bl)
earlyad.hipp.desc$AmyloidPos_bl <- factor(earlyad.hipp.desc$AmyloidPos_bl)
colnames(earlyad.hipp.desc)   <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                     "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                     "Hippocampus (L/R Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                     "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.hipp.desc <- table1(earlyad.hipp.desc, splitby = ~ Diagnosis)
earlyad.hipp.desc <- earlyad.hipp.desc$Table1
earlyad.hipp.desc <- earlyad.hipp.desc[c(1, 6:9, 11:nrow(earlyad.hipp.desc)),]


#Enriched for TAU
earlyad.tplus.hipp.desc <- cs.earlyad.tplus.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                        "hipp_average" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                        "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

earlyad.tplus.hipp.desc$TauPos_bl <- factor(earlyad.tplus.hipp.desc$TauPos_bl)
earlyad.tplus.hipp.desc$AmyloidPos_bl <- factor(earlyad.tplus.hipp.desc$AmyloidPos_bl)

colnames(earlyad.tplus.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                         "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                         "Hippocampus (L/R Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                         "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.tplus.hipp.desc <- table1(earlyad.tplus.hipp.desc, splitby = ~ Diagnosis)
earlyad.tplus.hipp.desc <- earlyad.tplus.hipp.desc$Table1
earlyad.tplus.hipp.desc <- earlyad.tplus.hipp.desc[c(1, 6:9, 11:nrow(earlyad.tplus.hipp.desc)),]


#Neuroenriched Tau+
earlyad.neuroenriched.tplus.hipp.desc <- cs.earlyad.neuroenriched.tplus.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                    "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                    "hipp_average" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                                                    "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

earlyad.neuroenriched.tplus.hipp.desc$TauPos_bl <- factor(earlyad.neuroenriched.tplus.hipp.desc$TauPos_bl)
earlyad.neuroenriched.tplus.hipp.desc$AmyloidPos_bl <- factor(earlyad.neuroenriched.tplus.hipp.desc$AmyloidPos_bl)

colnames(earlyad.neuroenriched.tplus.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                       "Hippocampus (L/R Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.tplus.hipp.desc <- table1(earlyad.neuroenriched.tplus.hipp.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.tplus.hipp.desc <- earlyad.neuroenriched.tplus.hipp.desc$Table1
earlyad.neuroenriched.tplus.hipp.desc <- earlyad.neuroenriched.tplus.hipp.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.tplus.hipp.desc)),]



#Descriptive statistics tables
######################## mPACCtrailsB ########################
earlyad.mpacc.desc <- cs.earlyad.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                            "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                            "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                            "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

earlyad.mpacc.desc$TauPos_bl <- factor(earlyad.mpacc.desc$TauPos_bl)
earlyad.mpacc.desc$AmyloidPos_bl <- factor(earlyad.mpacc.desc$AmyloidPos_bl)
colnames(earlyad.mpacc.desc)   <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                     "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                     "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                     "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.mpacc.desc <- table1(earlyad.mpacc.desc, splitby = ~ Diagnosis)
earlyad.mpacc.desc <- earlyad.mpacc.desc$Table1
earlyad.mpacc.desc <- earlyad.mpacc.desc[c(1, 6:9, 11:nrow(earlyad.mpacc.desc)),]


#Enriched for TAU
earlyad.tplus.mpacc.desc <- cs.earlyad.tplus.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                        "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                        "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

earlyad.tplus.mpacc.desc$TauPos_bl <- factor(earlyad.tplus.mpacc.desc$TauPos_bl)
earlyad.tplus.mpacc.desc$AmyloidPos_bl <- factor(earlyad.tplus.mpacc.desc$AmyloidPos_bl)

colnames(earlyad.tplus.mpacc.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                         "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                         "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                         "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.tplus.mpacc.desc <- table1(earlyad.tplus.mpacc.desc, splitby = ~ Diagnosis)
earlyad.tplus.mpacc.desc <- earlyad.tplus.mpacc.desc$Table1
earlyad.tplus.mpacc.desc <- earlyad.tplus.mpacc.desc[c(1, 6:9, 11:nrow(earlyad.tplus.mpacc.desc)),]


#Neuroenriched Tau+
earlyad.neuroenriched.tplus.mpacc.desc <- cs.earlyad.neuroenriched.tplus.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                    "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                    "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                                                    "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

earlyad.neuroenriched.tplus.mpacc.desc$TauPos_bl <- factor(earlyad.neuroenriched.tplus.mpacc.desc$TauPos_bl)
earlyad.neuroenriched.tplus.mpacc.desc$AmyloidPos_bl <- factor(earlyad.neuroenriched.tplus.mpacc.desc$AmyloidPos_bl)

colnames(earlyad.neuroenriched.tplus.mpacc.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                       "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.tplus.mpacc.desc <- table1(earlyad.neuroenriched.tplus.mpacc.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.tplus.mpacc.desc <- earlyad.neuroenriched.tplus.mpacc.desc$Table1
earlyad.neuroenriched.tplus.mpacc.desc <- earlyad.neuroenriched.tplus.mpacc.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.tplus.mpacc.desc)),]


#fit repeated measures model on longitudinal data
# want observed scores, just use time and dont control for other variables


long.earlyad.adas13 <- MMRMTime(long.earlyad.adas13)
long.earlyad.tplus.adas13 <- MMRMTime(long.earlyad.tplus.adas13)
long.earlyad.neuroenriched.tplus.adas13 <- MMRMTime(long.earlyad.neuroenriched.tplus.adas13)


long.earlyad.hipp <-  MMRMTime(long.earlyad.hipp)
long.earlyad.tplus.hipp <- MMRMTime(long.earlyad.tplus.hipp)
long.earlyad.neuroenriched.tplus.hipp <- MMRMTime(long.earlyad.neuroenriched.tplus.hipp)


long.earlyad.mpacc <- MMRMTime(long.earlyad.mpacc)
long.earlyad.tplus.mpacc <- MMRMTime(long.earlyad.tplus.mpacc)
long.earlyad.neuroenriched.tplus.mpacc <- MMRMTime(long.earlyad.neuroenriched.tplus.mpacc)


######################## ADAS13 ######################## 
mmrm.earlyad.adas13        <- gls(ADAS13~factor(new_time_mmrm),
                                  na.action=na.omit, data=long.earlyad.adas13,
                                  correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                  weights=nlme::varIdent(form=~1|new_time_mmrm))

#Enriched for TAU
mmrm.earlyad.tplus.adas13<- gls(ADAS13~factor(new_time_mmrm),
                                na.action=na.omit, data=long.earlyad.tplus.adas13,
                                correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                weights=nlme::varIdent(form=~1|new_time_mmrm))


#Neuroenriched Tau+
mmrm.earlyad.neuroenriched.tplus.adas13<- gls(ADAS13~factor(new_time_mmrm),
                                              na.action=na.omit, data=long.earlyad.neuroenriched.tplus.adas13,
                                              correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                              weights=nlme::varIdent(form=~1|new_time_mmrm))



######################## HIPPOCAMPUS ######################## 
mmrm.earlyad.hipp   <- gls(hipp_average~factor(new_time_mmrm),
                           na.action=na.omit, data=long.earlyad.hipp,
                           correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                           weights=nlme::varIdent(form=~1|new_time_mmrm))
mmrm.earlyad.hipp.test   <- gls(hipp_average~factor(new_time_mmrm),
                           na.action=na.omit, data=long.earlyad.hipp,
                           correlation=nlme::corSymm(form=~1|RID),
                           weights=nlme::varIdent(form=~1|new_time_mmrm))


#Enriched for TAU
mmrm.earlyad.tplus.hipp <- gls(hipp_average~factor(new_time_mmrm),
                               na.action=na.omit, data=long.earlyad.tplus.hipp,
                               correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                               weights=nlme::varIdent(form=~1|new_time_mmrm))
#Neuroenriched Tau+
mmrm.earlyad.neuroenriched.tplus.hipp <- gls(hipp_average~factor(new_time_mmrm),
                                             na.action=na.omit, data=long.earlyad.neuroenriched.tplus.hipp,
                                             correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                             weights=nlme::varIdent(form=~1|new_time_mmrm))





######################## MPACC ######################## 
mmrm.earlyad.mpacc        <- gls(mPACCtrailsB~factor(new_time_mmrm),
                                 na.action=na.omit, data=long.earlyad.mpacc,
                                 correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                 weights=nlme::varIdent(form=~1|new_time_mmrm))


#Enriched for TAU
mmrm.earlyad.tplus.mpacc<- gls(mPACCtrailsB~factor(new_time_mmrm),
                               na.action=na.omit, data=long.earlyad.tplus.mpacc,
                               correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                               weights=nlme::varIdent(form=~1|new_time_mmrm))

#Neuroenriched Tau+
mmrm.earlyad.neuroenriched.tplus.mpacc<- gls(mPACCtrailsB~factor(new_time_mmrm),
                                             na.action=na.omit, data=long.earlyad.neuroenriched.tplus.mpacc,
                                             correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                             weights=nlme::varIdent(form=~1|new_time_mmrm))



#get coefficients and standard errors

######################## ADAS13 ######################## 
coefs.earlyad.adas13     <- unname(coef(summary(mmrm.earlyad.adas13))[, "Value"])
stderror.earlyad.adas13  <- unname(coef(summary(mmrm.earlyad.adas13))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.adas13     <- unname(coef(summary(mmrm.earlyad.tplus.adas13))[, "Value"])
stderror.earlyad.tplus.adas13  <- unname(coef(summary(mmrm.earlyad.tplus.adas13))[, "Std.Error"])


#Neuroenriched Tau+
coefs.earlyad.neuroenriched.tplus.adas13     <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.adas13))[, "Value"])
stderror.earlyad.neuroenriched.tplus.adas13  <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.adas13))[, "Std.Error"])


######################## HIPPOCAMPUS ######################## 
coefs.earlyad.hipp     <- unname(coef(summary(mmrm.earlyad.hipp))[, "Value"])
stderror.earlyad.hipp  <- unname(coef(summary(mmrm.earlyad.hipp))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.hipp     <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Value"])
stderror.earlyad.tplus.hipp  <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Std.Error"])



#Neuroenriched Tau+
coefs.earlyad.neuroenriched.tplus.hipp     <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.hipp))[, "Value"])
stderror.earlyad.neuroenriched.tplus.hipp  <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.hipp))[, "Std.Error"])



######################## MPACC ######################## 
coefs.earlyad.mpacc     <- unname(coef(summary(mmrm.earlyad.mpacc))[, "Value"])
stderror.earlyad.mpacc  <- unname(coef(summary(mmrm.earlyad.mpacc))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.mpacc     <- unname(coef(summary(mmrm.earlyad.tplus.mpacc))[, "Value"])
stderror.earlyad.tplus.mpacc  <- unname(coef(summary(mmrm.earlyad.tplus.mpacc))[, "Std.Error"])


#Neuroenriched Tau+
coefs.earlyad.neuroenriched.tplus.mpacc     <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.mpacc))[, "Value"])
stderror.earlyad.neuroenriched.tplus.mpacc  <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.mpacc))[, "Std.Error"])



#create dataset for plotting
######################## ADAS13 ######################## 
earlyad.adas13.plotting.data <- data.frame("coefs" = coefs.earlyad.adas13,
                                           "stderr" = stderror.earlyad.adas13)
earlyad.adas13.plotting.data$ci_low     <-  -1.96* earlyad.adas13.plotting.data$stderr
earlyad.adas13.plotting.data$ci_hi      <-   1.96* earlyad.adas13.plotting.data$stderr
earlyad.adas13.plotting.data$means      <-  earlyad.adas13.plotting.data$coefs
earlyad.adas13.plotting.data$means[2:5] <-  earlyad.adas13.plotting.data$means[2:5] + earlyad.adas13.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.adas13.plotting.data$ci_low     <-  earlyad.adas13.plotting.data$ci_low+ earlyad.adas13.plotting.data$means
earlyad.adas13.plotting.data$ci_hi      <-  earlyad.adas13.plotting.data$ci_hi+ earlyad.adas13.plotting.data$means
earlyad.adas13.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.adas13.plotting.data$Enrichment <- "No Enrichment"

#Enriched for TAU
earlyad.tplus.adas13.plotting.data <- data.frame("coefs" = coefs.earlyad.tplus.adas13,
                                                 "stderr" = stderror.earlyad.tplus.adas13)
earlyad.tplus.adas13.plotting.data$ci_low     <-  -1.96* earlyad.tplus.adas13.plotting.data$stderr
earlyad.tplus.adas13.plotting.data$ci_hi      <-   1.96* earlyad.tplus.adas13.plotting.data$stderr
earlyad.tplus.adas13.plotting.data$means      <-  earlyad.tplus.adas13.plotting.data$coefs
earlyad.tplus.adas13.plotting.data$means[2:5] <-  earlyad.tplus.adas13.plotting.data$means[2:5] + earlyad.tplus.adas13.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.tplus.adas13.plotting.data$ci_low     <-  earlyad.tplus.adas13.plotting.data$ci_low+ earlyad.tplus.adas13.plotting.data$means
earlyad.tplus.adas13.plotting.data$ci_hi      <-  earlyad.tplus.adas13.plotting.data$ci_hi+ earlyad.tplus.adas13.plotting.data$means
earlyad.tplus.adas13.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.tplus.adas13.plotting.data$Enrichment <- "Tau+"

#Neuroenriched Tau+
earlyad.neuroenriched.tplus.adas13.plotting.data <- data.frame("coefs" = coefs.earlyad.neuroenriched.tplus.adas13,
                                                               "stderr" = stderror.earlyad.neuroenriched.tplus.adas13)
earlyad.neuroenriched.tplus.adas13.plotting.data$ci_low     <-  -1.96* earlyad.neuroenriched.tplus.adas13.plotting.data$stderr
earlyad.neuroenriched.tplus.adas13.plotting.data$ci_hi      <-   1.96* earlyad.neuroenriched.tplus.adas13.plotting.data$stderr
earlyad.neuroenriched.tplus.adas13.plotting.data$means      <-  earlyad.neuroenriched.tplus.adas13.plotting.data$coefs
earlyad.neuroenriched.tplus.adas13.plotting.data$means[2:5] <-  earlyad.neuroenriched.tplus.adas13.plotting.data$means[2:5] + earlyad.neuroenriched.tplus.adas13.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.neuroenriched.tplus.adas13.plotting.data$ci_low     <-  earlyad.neuroenriched.tplus.adas13.plotting.data$ci_low+ earlyad.neuroenriched.tplus.adas13.plotting.data$means
earlyad.neuroenriched.tplus.adas13.plotting.data$ci_hi      <-  earlyad.neuroenriched.tplus.adas13.plotting.data$ci_hi+ earlyad.neuroenriched.tplus.adas13.plotting.data$means
earlyad.neuroenriched.tplus.adas13.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.neuroenriched.tplus.adas13.plotting.data$Enrichment <- "No Copathologies Tau+"




######################## HIPPOCAMPUS ######################## 
earlyad.hipp.plotting.data <- data.frame("coefs" = coefs.earlyad.hipp,
                                         "stderr" = stderror.earlyad.hipp)
earlyad.hipp.plotting.data$ci_low     <-  -1.96* earlyad.hipp.plotting.data$stderr
earlyad.hipp.plotting.data$ci_hi      <-   1.96* earlyad.hipp.plotting.data$stderr
earlyad.hipp.plotting.data$means      <-  earlyad.hipp.plotting.data$coefs
earlyad.hipp.plotting.data$means[2:5] <-  earlyad.hipp.plotting.data$means[2:5] + earlyad.hipp.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.hipp.plotting.data$ci_low     <-  earlyad.hipp.plotting.data$ci_low+ earlyad.hipp.plotting.data$means
earlyad.hipp.plotting.data$ci_hi      <-  earlyad.hipp.plotting.data$ci_hi+ earlyad.hipp.plotting.data$means

earlyad.hipp.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.hipp.plotting.data$Enrichment <- "No Enrichment"

#Enriched for TAU
earlyad.tplus.hipp.plotting.data <- data.frame("coefs" = coefs.earlyad.tplus.hipp,
                                               "stderr" = stderror.earlyad.tplus.hipp)
earlyad.tplus.hipp.plotting.data$ci_low     <-  -1.96* earlyad.tplus.hipp.plotting.data$stderr
earlyad.tplus.hipp.plotting.data$ci_hi      <-   1.96* earlyad.tplus.hipp.plotting.data$stderr
earlyad.tplus.hipp.plotting.data$means      <-  earlyad.tplus.hipp.plotting.data$coefs
earlyad.tplus.hipp.plotting.data$means[2:5] <-  earlyad.tplus.hipp.plotting.data$means[2:5] + earlyad.tplus.hipp.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.tplus.hipp.plotting.data$ci_low     <-  earlyad.tplus.hipp.plotting.data$ci_low+ earlyad.tplus.hipp.plotting.data$means
earlyad.tplus.hipp.plotting.data$ci_hi      <-  earlyad.tplus.hipp.plotting.data$ci_hi+ earlyad.tplus.hipp.plotting.data$means
earlyad.tplus.hipp.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.tplus.hipp.plotting.data$Enrichment <- "Tau+"



#Neuroenriched Tau+
earlyad.neuroenriched.tplus.hipp.plotting.data <- data.frame("coefs" = coefs.earlyad.neuroenriched.tplus.hipp,
                                                             "stderr" = stderror.earlyad.neuroenriched.tplus.hipp)
earlyad.neuroenriched.tplus.hipp.plotting.data$ci_low     <-  -1.96* earlyad.neuroenriched.tplus.hipp.plotting.data$stderr
earlyad.neuroenriched.tplus.hipp.plotting.data$ci_hi      <-   1.96* earlyad.neuroenriched.tplus.hipp.plotting.data$stderr
earlyad.neuroenriched.tplus.hipp.plotting.data$means      <-  earlyad.neuroenriched.tplus.hipp.plotting.data$coefs
earlyad.neuroenriched.tplus.hipp.plotting.data$means[2:5] <-  earlyad.neuroenriched.tplus.hipp.plotting.data$means[2:5] + earlyad.neuroenriched.tplus.hipp.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.neuroenriched.tplus.hipp.plotting.data$ci_low     <-  earlyad.neuroenriched.tplus.hipp.plotting.data$ci_low+ earlyad.neuroenriched.tplus.hipp.plotting.data$means
earlyad.neuroenriched.tplus.hipp.plotting.data$ci_hi      <-  earlyad.neuroenriched.tplus.hipp.plotting.data$ci_hi+ earlyad.neuroenriched.tplus.hipp.plotting.data$means
earlyad.neuroenriched.tplus.hipp.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.neuroenriched.tplus.hipp.plotting.data$Enrichment <- "No Copathologies Tau+"


######################## MPACC ######################## 
earlyad.mpacc.plotting.data <- data.frame("coefs" = coefs.earlyad.mpacc,
                                          "stderr" = stderror.earlyad.mpacc)
earlyad.mpacc.plotting.data$ci_low     <-  -1.96* earlyad.mpacc.plotting.data$stderr
earlyad.mpacc.plotting.data$ci_hi      <-   1.96* earlyad.mpacc.plotting.data$stderr
earlyad.mpacc.plotting.data$means      <-  earlyad.mpacc.plotting.data$coefs
earlyad.mpacc.plotting.data$means[2:5] <-  earlyad.mpacc.plotting.data$means[2:5] + earlyad.mpacc.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.mpacc.plotting.data$ci_low     <-  earlyad.mpacc.plotting.data$ci_low+ earlyad.mpacc.plotting.data$means
earlyad.mpacc.plotting.data$ci_hi      <-  earlyad.mpacc.plotting.data$ci_hi+ earlyad.mpacc.plotting.data$means
earlyad.mpacc.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.mpacc.plotting.data$Enrichment <- "No Enrichment"

#Enriched for TAU
earlyad.tplus.mpacc.plotting.data <- data.frame("coefs" = coefs.earlyad.tplus.mpacc,
                                                "stderr" = stderror.earlyad.tplus.mpacc)
earlyad.tplus.mpacc.plotting.data$ci_low     <-  -1.96* earlyad.tplus.mpacc.plotting.data$stderr
earlyad.tplus.mpacc.plotting.data$ci_hi      <-   1.96* earlyad.tplus.mpacc.plotting.data$stderr
earlyad.tplus.mpacc.plotting.data$means      <-  earlyad.tplus.mpacc.plotting.data$coefs
earlyad.tplus.mpacc.plotting.data$means[2:5] <-  earlyad.tplus.mpacc.plotting.data$means[2:5] + earlyad.tplus.mpacc.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.tplus.mpacc.plotting.data$ci_low     <-  earlyad.tplus.mpacc.plotting.data$ci_low+ earlyad.tplus.mpacc.plotting.data$means
earlyad.tplus.mpacc.plotting.data$ci_hi      <-  earlyad.tplus.mpacc.plotting.data$ci_hi+ earlyad.tplus.mpacc.plotting.data$means
earlyad.tplus.mpacc.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.tplus.mpacc.plotting.data$Enrichment <- "Tau+"


#Neuroenriched Tau+
earlyad.neuroenriched.tplus.mpacc.plotting.data <- data.frame("coefs" = coefs.earlyad.neuroenriched.tplus.mpacc,
                                                              "stderr" = stderror.earlyad.neuroenriched.tplus.mpacc)
earlyad.neuroenriched.tplus.mpacc.plotting.data$ci_low     <-  -1.96* earlyad.neuroenriched.tplus.mpacc.plotting.data$stderr
earlyad.neuroenriched.tplus.mpacc.plotting.data$ci_hi      <-   1.96* earlyad.neuroenriched.tplus.mpacc.plotting.data$stderr
earlyad.neuroenriched.tplus.mpacc.plotting.data$means      <-  earlyad.neuroenriched.tplus.mpacc.plotting.data$coefs
earlyad.neuroenriched.tplus.mpacc.plotting.data$means[2:5] <-  earlyad.neuroenriched.tplus.mpacc.plotting.data$means[2:5] + earlyad.neuroenriched.tplus.mpacc.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.neuroenriched.tplus.mpacc.plotting.data$ci_low     <-  earlyad.neuroenriched.tplus.mpacc.plotting.data$ci_low+ earlyad.neuroenriched.tplus.mpacc.plotting.data$means
earlyad.neuroenriched.tplus.mpacc.plotting.data$ci_hi      <-  earlyad.neuroenriched.tplus.mpacc.plotting.data$ci_hi+ earlyad.neuroenriched.tplus.mpacc.plotting.data$means
earlyad.neuroenriched.tplus.mpacc.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.neuroenriched.tplus.mpacc.plotting.data$Enrichment <- "No Copathologies Tau+"



## Plots
# Combine non-enriched and tau enriched data for plotting

######################## ADAS13 ######################## 
earlyad.adas13.plotting.data <- rbind(earlyad.adas13.plotting.data,
                                      earlyad.tplus.adas13.plotting.data,
                                      earlyad.neuroenriched.tplus.adas13.plotting.data)
earlyad.adas13.mmrm.plot <- ggplot(earlyad.adas13.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.adas13.mmrm.plot <- earlyad.adas13.mmrm.plot + xlab("Time (Years)") + ylab("Mean ADAS-13 (95% CI)") + ylim(20, 40)
tiff('mpaccmrmm.tiff', units="in", width=6, height=3, res=500)
earlyad.mpacc.mmrm.plot
dev.off()
######################## HIPPOCAMPUS ######################## 
earlyad.hipp.plotting.data <- rbind(earlyad.hipp.plotting.data,
                                    earlyad.tplus.hipp.plotting.data,
                                    earlyad.neuroenriched.tplus.hipp.plotting.data)
earlyad.hipp.mmrm.plot <- ggplot(earlyad.hipp.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.hipp.mmrm.plot <- earlyad.hipp.mmrm.plot + xlab("Time (Years)") + ylab("Mean Hippocampus (95% CI)") + ylim(2700, 3800)

######################## MPACC ######################## 
earlyad.mpacc.plotting.data <- rbind(earlyad.mpacc.plotting.data,
                                     earlyad.tplus.mpacc.plotting.data,
                                     earlyad.neuroenriched.tplus.mpacc.plotting.data)
earlyad.mpacc.mmrm.plot <- ggplot(earlyad.mpacc.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.mpacc.mmrm.plot <- earlyad.mpacc.mmrm.plot + xlab("Time (Years)") + ylab("Mean mPACCtrailsB (95% CI)") 


######################## ADAS13 ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.earlyad.adas13 <- StratifyContVar(cs.earlyad.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.adas13 <- StratifyContVar(cs.earlyad.tplus.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.tplus.adas13 <- StratifyContVar(cs.earlyad.neuroenriched.tplus.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))




#Assign treatment and placebo by blocks
long.earlyad.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.adas13, long.earlyad.adas13, no.prop = TRUE)
long.earlyad.tplus.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.adas13, long.earlyad.tplus.adas13, no.prop = TRUE)
long.earlyad.neuroenriched.tplus.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.tplus.adas13, long.earlyad.neuroenriched.tplus.adas13, no.prop = TRUE)

######################## HIPPOCAMPUS ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.earlyad.hipp <- StratifyContVar(cs.earlyad.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.hipp <- StratifyContVar(cs.earlyad.tplus.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.tplus.hipp <- StratifyContVar(cs.earlyad.neuroenriched.tplus.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))


#Assign treatment and placebo by blocks
long.earlyad.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.hipp, long.earlyad.hipp, no.prop = TRUE)
long.earlyad.tplus.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.hipp, long.earlyad.tplus.hipp, no.prop = TRUE)
long.earlyad.neuroenriched.tplus.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.tplus.hipp, long.earlyad.neuroenriched.tplus.hipp, no.prop = TRUE)


######################## MPACC ######################## 

# Categorize continuous variables for block randomization
cs.earlyad.mpacc <- StratifyContVar(cs.earlyad.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.mpacc <- StratifyContVar(cs.earlyad.tplus.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.tplus.mpacc <- StratifyContVar(cs.earlyad.neuroenriched.tplus.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))


#Assign treatment and placebo by blocks
long.earlyad.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.mpacc, long.earlyad.mpacc, no.prop = TRUE)
long.earlyad.tplus.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.mpacc, long.earlyad.tplus.mpacc, no.prop = TRUE)
long.earlyad.neuroenriched.tplus.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.tplus.mpacc, long.earlyad.neuroenriched.tplus.mpacc, no.prop = TRUE)



######################## ADAS13 ######################## 

formula.earlyad.adas13.rs                              <- "ADAS13~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
long.earlyad.adas13_with_treatment                     <- StratifyContinuous(long.earlyad.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.tplus.adas13_with_treatment               <- StratifyContinuous(long.earlyad.tplus.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.neuroenriched.tplus.adas13_with_treatment <- StratifyContinuous(long.earlyad.neuroenriched.tplus.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))

model.earlyad.adas13.rs                     <- MapLmer(newdata = long.earlyad.adas13_with_treatment,
                                                       formula.model = formula.earlyad.adas13.rs)

model.earlyad.tplus.adas13.rs               <- MapLmer(newdata = long.earlyad.tplus.adas13_with_treatment,
                                                       formula.model = formula.earlyad.adas13.rs)
model.earlyad.neuroenriched.tplus.adas13.rs <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.adas13_with_treatment,
                                                       formula.model = formula.earlyad.adas13.rs)

formula.earlyad.adas13.simulation.rs                   <- "ADAS13 ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
simulation.model.earlyad.adas13.rs                     <- BuildSimulationModelNoPath(model.earlyad.adas13.rs, formula.earlyad.adas13.simulation.rs, long.earlyad.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.adas13.rs               <- BuildSimulationModelNoPath(model.earlyad.tplus.adas13.rs, formula.earlyad.adas13.simulation.rs, long.earlyad.tplus.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.adas13.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.adas13.rs, formula.earlyad.adas13.simulation.rs, long.earlyad.neuroenriched.tplus.adas13_with_treatment, "not-controlled")

relcontr.earlyad.adas13                            <- GetRelContributions(model.earlyad.adas13.rs, long.earlyad.adas13_with_treatment)
relcontr.earlyad.tplus.adas13                      <- GetRelContributions(model.earlyad.tplus.adas13.rs, long.earlyad.tplus.adas13_with_treatment)
relcontr.earlyad.neuroenriched.tplus.adas13        <- GetRelContributions(model.earlyad.neuroenriched.tplus.adas13.rs, long.earlyad.neuroenriched.tplus.adas13_with_treatment)


relcontrcount.earlyad.adas13        <- BuildNeuroCountPlot(relcontr.earlyad.adas13)
relcontrcount.earlyad.tplus.adas13  <- BuildNeuroCountPlot(relcontr.earlyad.tplus.adas13)





simulation.model.earlyad.adas13.tau.rs                     <- BuildSimulationModelNoPath(model.earlyad.adas13.rs, 
                                                                                         formula.earlyad.adas13.simulation.rs, 
                                                                                         long.earlyad.adas13_with_treatment,  
                                                                                         relcontr.earlyad.adas13$Rates_Data$Mean_Rate[[4]])

simulation.model.earlyad.tplus.adas13.tau.rs               <- BuildSimulationModelNoPath(model.earlyad.tplus.adas13.rs, 
                                                                                         formula.earlyad.adas13.simulation.rs, 
                                                                                         long.earlyad.tplus.adas13_with_treatment,  
                                                                                         relcontr.earlyad.tplus.adas13$Rates_Data$Mean_Rate[[4]])

simulation.model.earlyad.neuroenriched.tplus.adas13.tau.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.adas13.rs, 
                                                                                         formula.earlyad.adas13.simulation.rs, 
                                                                                         long.earlyad.neuroenriched.tplus.adas13_with_treatment, 
                                                                                         relcontr.earlyad.neuroenriched.tplus.adas13$Rates_Data$Mean_Rate[[4]])



plotting.simulation.model.earlyad.adas13.tau.rs <- simulation.model.earlyad.adas13.tau.rs
fixef(plotting.simulation.model.earlyad.adas13.tau.rs)["new_time"] <- relcontr.earlyad.adas13$Rates_Data$Mean_Rate[[4]]
plotting.simulation.model.earlyad.tplus.adas13.tau.rs <- simulation.model.earlyad.tplus.adas13.tau.rs
fixef(plotting.simulation.model.earlyad.tplus.adas13.tau.rs)["new_time"] <- relcontr.earlyad.tplus.adas13$Rates_Data$Mean_Rate[[4]]



dpms.earlyad.adas13 <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.adas13.rs,
                                     "Tau+" = simulation.model.earlyad.tplus.adas13.rs,
                                     "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.adas13.rs), ylab="ADAS13", ylim.low = 19, ylim.high = 35.5)

dpms.earlyad.adas13.tau <- DPMPlots(list("No Enrichment" = plotting.simulation.model.earlyad.adas13.tau.rs,
                                     "Tau+" = plotting.simulation.model.earlyad.tplus.adas13.tau.rs,
                                     "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.adas13.tau.rs), ylab="ADAS13", ylim.low = 19, ylim.high = 35.5)


adas.overall.vs.tau <- OverallvsTau(dpms.earlyad.adas13, dpms.earlyad.adas13.tau, ylab = "ADAS13", ymin = 19, ymax =33.5)
tiff('adastplussdpm.tiff', units="in", width=7, height=3, res=400)
adas.overall.vs.tau[[2]]
dev.off()
######################## Hippcampus ######################## 



formula.earlyad.hipp.rs                              <- "hipp_average~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
long.earlyad.hipp_with_treatment                     <- StratifyContinuous(long.earlyad.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.tplus.hipp_with_treatment               <- StratifyContinuous(long.earlyad.tplus.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.neuroenriched.tplus.hipp_with_treatment <- StratifyContinuous(long.earlyad.neuroenriched.tplus.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
model.earlyad.hipp.rs                     <- MapLmer(newdata = long.earlyad.hipp_with_treatment,
                                                       formula.model = formula.earlyad.hipp.rs)


model.earlyad.tplus.hipp.rs               <- MapLmer(newdata = long.earlyad.tplus.hipp_with_treatment,
                                                       formula.model = formula.earlyad.hipp.rs)

model.earlyad.neuroenriched.tplus.hipp.rs <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.hipp_with_treatment,
                                                       formula.model = formula.earlyad.hipp.rs)

formula.earlyad.hipp.simulation.rs                   <- "hipp_average ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
simulation.model.earlyad.hipp.rs                     <- BuildSimulationModelNoPath(model.earlyad.hipp.rs, formula.earlyad.hipp.simulation.rs, long.earlyad.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.hipp.rs               <- BuildSimulationModelNoPath(model.earlyad.tplus.hipp.rs, formula.earlyad.hipp.simulation.rs, long.earlyad.tplus.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.hipp.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.hipp.rs, formula.earlyad.hipp.simulation.rs, long.earlyad.neuroenriched.tplus.hipp_with_treatment, "not-controlled")


relcontr.earlyad.hipp                            <- GetRelContributions(model.earlyad.hipp.rs, long.earlyad.hipp_with_treatment)
relcontr.earlyad.tplus.hipp                      <- GetRelContributions(model.earlyad.tplus.hipp.rs, long.earlyad.tplus.hipp_with_treatment)
relcontr.earlyad.neuroenriched.tplus.hipp        <- GetRelContributions(model.earlyad.neuroenriched.tplus.hipp.rs, long.earlyad.neuroenriched.tplus.hipp_with_treatment)


relcontrcount.earlyad.hipp        <- BuildNeuroCountPlot(relcontr.earlyad.hipp)
relcontrcount.earlyad.tplus.hipp  <- BuildNeuroCountPlot(relcontr.earlyad.tplus.hipp)



simulation.model.earlyad.hipp.tau.rs                     <- BuildSimulationModelNoPath(model.earlyad.hipp.rs, 
                                                                                         formula.earlyad.hipp.simulation.rs, 
                                                                                         long.earlyad.hipp_with_treatment,  
                                                                                         relcontr.earlyad.hipp$Rates_Data$Mean_Rate[[4]])

simulation.model.earlyad.tplus.hipp.tau.rs               <- BuildSimulationModelNoPath(model.earlyad.tplus.hipp.rs, 
                                                                                         formula.earlyad.hipp.simulation.rs, 
                                                                                         long.earlyad.tplus.hipp_with_treatment,  
                                                                                         relcontr.earlyad.tplus.hipp$Rates_Data$Mean_Rate[[4]])

simulation.model.earlyad.neuroenriched.tplus.hipp.tau.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.hipp.rs, 
                                                                                         formula.earlyad.hipp.simulation.rs, 
                                                                                         long.earlyad.neuroenriched.tplus.hipp_with_treatment, 
                                                                                         relcontr.earlyad.neuroenriched.tplus.hipp$Rates_Data$Mean_Rate[[4]])
plotting.simulation.model.earlyad.hipp.tau.rs <- simulation.model.earlyad.hipp.tau.rs
fixef(plotting.simulation.model.earlyad.hipp.tau.rs)["new_time"] <- relcontr.earlyad.hipp$Rates_Data$Mean_Rate[[4]]
plotting.simulation.model.earlyad.tplus.hipp.tau.rs <- simulation.model.earlyad.tplus.hipp.tau.rs
fixef(plotting.simulation.model.earlyad.tplus.hipp.tau.rs)["new_time"] <- relcontr.earlyad.tplus.hipp$Rates_Data$Mean_Rate[[4]]


dpms.earlyad.hipp <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.hipp.rs,
                                     "Tau+" = simulation.model.earlyad.tplus.hipp.rs,
                                     "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.hipp.rs), ylab="Hippocampus", ylim.low = 2800, ylim.high = 4000)

dpms.earlyad.hipp.tau <- DPMPlots(list("No Enrichment" = plotting.simulation.model.earlyad.hipp.tau.rs,
                                         "Tau+" = plotting.simulation.model.earlyad.tplus.hipp.tau.rs,
                                         "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.hipp.tau.rs), ylab="Hippocampus", ylim.low = 2800, ylim.high = 4000)


hipp.overall.vs.tau <- OverallvsTau(dpms.earlyad.hipp, dpms.earlyad.hipp.tau, ylab = "Hippocampus", ymin = 2800, ymax =4100)
tiff('hipptplussdpm.tiff', units="in", width=7, height=3, res=400)
hipp.overall.vs.tau[[2]]
dev.off()
######################## mPACCtrailsB ######################## 
formula.earlyad.mpacc.rs                              <- "mPACCtrailsB~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
long.earlyad.mpacc_with_treatment                     <- StratifyContinuous(long.earlyad.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.tplus.mpacc_with_treatment               <- StratifyContinuous(long.earlyad.tplus.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.neuroenriched.tplus.mpacc_with_treatment <- StratifyContinuous(long.earlyad.neuroenriched.tplus.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))



model.earlyad.mpacc.rs                     <- MapLmer(newdata = long.earlyad.mpacc_with_treatment,
                                                     formula.model = formula.earlyad.mpacc.rs)
model.earlyad.tplus.mpacc.rs               <- MapLmer(newdata = long.earlyad.tplus.mpacc_with_treatment,
                                                     formula.model = formula.earlyad.mpacc.rs)
model.earlyad.neuroenriched.tplus.mpacc.rs <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.mpacc_with_treatment,
                                                     formula.model = formula.earlyad.mpacc.rs)

formula.earlyad.mpacc.simulation.rs                   <- "mPACCtrailsB ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
simulation.model.earlyad.mpacc.rs                     <- BuildSimulationModelNoPath(model.earlyad.mpacc.rs, formula.earlyad.mpacc.simulation.rs, long.earlyad.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.mpacc.rs               <- BuildSimulationModelNoPath(model.earlyad.tplus.mpacc.rs, formula.earlyad.mpacc.simulation.rs, long.earlyad.tplus.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.mpacc.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.mpacc.rs, formula.earlyad.mpacc.simulation.rs, long.earlyad.neuroenriched.tplus.mpacc_with_treatment, "not-controlled")


relcontr.earlyad.mpacc                            <- GetRelContributions(model.earlyad.mpacc.rs, long.earlyad.mpacc_with_treatment)
relcontr.earlyad.tplus.mpacc                      <- GetRelContributions(model.earlyad.tplus.mpacc.rs, long.earlyad.tplus.mpacc_with_treatment)
relcontr.earlyad.neuroenriched.tplus.mpacc        <- GetRelContributions(model.earlyad.neuroenriched.tplus.mpacc.rs, long.earlyad.neuroenriched.tplus.mpacc_with_treatment)


relcontrcount.earlyad.mpacc        <- BuildNeuroCountPlot(relcontr.earlyad.mpacc)
relcontrcount.earlyad.tplus.mpacc  <- BuildNeuroCountPlot(relcontr.earlyad.tplus.mpacc)



simulation.model.earlyad.mpacc.tau.rs                     <- BuildSimulationModelNoPath(model.earlyad.mpacc.rs, 
                                                                                       formula.earlyad.mpacc.simulation.rs, 
                                                                                       long.earlyad.mpacc_with_treatment,  
                                                                                       relcontr.earlyad.mpacc$Rates_Data$Mean_Rate[[4]])

simulation.model.earlyad.tplus.mpacc.tau.rs               <- BuildSimulationModelNoPath(model.earlyad.tplus.mpacc.rs, 
                                                                                       formula.earlyad.mpacc.simulation.rs, 
                                                                                       long.earlyad.tplus.mpacc_with_treatment,  
                                                                                       relcontr.earlyad.tplus.mpacc$Rates_Data$Mean_Rate[[4]])

simulation.model.earlyad.neuroenriched.tplus.mpacc.tau.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.mpacc.rs, 
                                                                                       formula.earlyad.mpacc.simulation.rs, 
                                                                                       long.earlyad.neuroenriched.tplus.mpacc_with_treatment, 
                                                                                       relcontr.earlyad.neuroenriched.tplus.mpacc$Rates_Data$Mean_Rate[[4]])



plotting.simulation.model.earlyad.mpacc.tau.rs <- simulation.model.earlyad.mpacc.tau.rs
fixef(plotting.simulation.model.earlyad.mpacc.tau.rs)["new_time"] <- relcontr.earlyad.mpacc$Rates_Data$Mean_Rate[[4]]
plotting.simulation.model.earlyad.tplus.mpacc.tau.rs <- simulation.model.earlyad.tplus.mpacc.tau.rs
fixef(plotting.simulation.model.earlyad.tplus.mpacc.tau.rs)["new_time"] <- relcontr.earlyad.tplus.mpacc$Rates_Data$Mean_Rate[[4]]


dpms.earlyad.mpacc <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.mpacc.rs,
                                   "Tau+" = simulation.model.earlyad.tplus.mpacc.rs,
                                   "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.mpacc.rs), ylab="mpaccocampus", ylim.low = 2800, ylim.high = 4000)

dpms.earlyad.mpacc.tau <- DPMPlots(list("No Enrichment" = plotting.simulation.model.earlyad.mpacc.tau.rs,
                                       "Tau+" = plotting.simulation.model.earlyad.tplus.mpacc.tau.rs,
                                       "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.mpacc.tau.rs), ylab="mpaccocampus", ylim.low = 2800, ylim.high = 4000)


mpacc.overall.vs.tau <- OverallvsTau(dpms.earlyad.mpacc, dpms.earlyad.mpacc.tau, ylab = "mPACCtrailsB", ymin = -18, ymax =-8)
tiff('mpacctplussdpm.tiff', units="in", width=7, height=3, res=400)
mpacc.overall.vs.tau[[2]]
dev.off()

if(save.rds) {



#ADAS#################################### #################################### ####################################

earlyadadas13.sim.list <- list("formula_largemodel" = formula.earlyad.adas13.simulation.rs,
                              "largemodel" = simulation.model.earlyad.adas13.rs,
                              "formula_smallmodel" = formula.earlyad.adas13.rs,
                              "smallmodel" = model.earlyad.adas13.rs,
                              "sample_sizes" = seq(100, 800, by=100),
                              "nsim"=500,
                              "data" = long.earlyad.adas13_with_treatment)



earlyadadas13.tplus.sim.list <- list("formula_largemodel" = formula.earlyad.adas13.simulation.rs,
                                    "largemodel" = simulation.model.earlyad.tplus.adas13.rs,
                                    "formula_smallmodel" = formula.earlyad.adas13.rs,
                                    "smallmodel" = model.earlyad.tplus.adas13.rs,
                                    "sample_sizes" = seq(100, 800, by=100),
                                    "nsim"=500,
                                    "data" = long.earlyad.tplus.adas13_with_treatment)

earlyadadas13.neuro.tplus.sim.list <- list("formula_largemodel" = formula.earlyad.adas13.simulation.rs,
                                          "largemodel" = simulation.model.earlyad.neuroenriched.tplus.adas13.rs,
                                          "formula_smallmodel" = formula.earlyad.adas13.rs,
                                          "smallmodel" = model.earlyad.neuroenriched.tplus.adas13.rs,
                                          "sample_sizes" = seq(100, 800, by=100),
                                          "nsim"=500,
                                          "data" = long.earlyad.neuroenriched.tplus.adas13_with_treatment)



earlyadadas13.sim.tau.list <- list("formula_largemodel" = formula.earlyad.adas13.simulation.rs,
                                  "largemodel" = simulation.model.earlyad.adas13.tau.rs,
                                  "formula_smallmodel" = formula.earlyad.adas13.rs,
                                  "smallmodel" = model.earlyad.adas13.rs,
                                  "sample_sizes" = seq(100, 800, by=100),
                                  "nsim"=500,
                                  "data" = long.earlyad.adas13_with_treatment)


earlyadadas13.tplus.sim.tau.list <- list("formula_largemodel" = formula.earlyad.adas13.simulation.rs,
                                        "largemodel" = simulation.model.earlyad.tplus.adas13.tau.rs,
                                        "formula_smallmodel" = formula.earlyad.adas13.rs,
                                        "smallmodel" = model.earlyad.tplus.adas13.rs,
                                        "sample_sizes" = seq(100, 800, by=100),
                                        "nsim"=500,
                                        "data" = long.earlyad.tplus.adas13_with_treatment)

earlyadadas13.neuro.tplus.sim.tau.list <- list("formula_largemodel" = formula.earlyad.adas13.simulation.rs,
                                              "largemodel" = simulation.model.earlyad.neuroenriched.tplus.adas13.tau.rs,
                                              "formula_smallmodel" = formula.earlyad.adas13.rs,
                                              "smallmodel" = model.earlyad.neuroenriched.tplus.adas13.rs,
                                              "sample_sizes" = seq(100, 800, by=100),
                                              "nsim"=500,
                                              "data" = long.earlyad.neuroenriched.tplus.adas13_with_treatment)


saveRDS(earlyadadas13.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/earlyadadas13.rds")
saveRDS(earlyadadas13.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/earlyadadas13_tplus.rds")
saveRDS(earlyadadas13.neuro.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/earlyadadas13_neuro_tplus.rds")


saveRDS(earlyadadas13.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/earlyadadas13_tau.rds")
saveRDS(earlyadadas13.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/earlyadadas13_tplus_tau.rds")
saveRDS(earlyadadas13.neuro.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/earlyadadas13_neuro_tplus_tau.rds")




#Hipp#################################### #################################### ####################################


earlyadhipp.sim.list <- list("formula_largemodel" = formula.earlyad.hipp.simulation.rs,
                               "largemodel" = simulation.model.earlyad.hipp.rs,
                               "formula_smallmodel" = formula.earlyad.hipp.rs,
                               "smallmodel" = model.earlyad.hipp.rs,
                               "sample_sizes" = seq(100, 800, by=100),
                               "nsim"=500,
                               "data" = long.earlyad.hipp_with_treatment)


earlyadhipp.tplus.sim.list <- list("formula_largemodel" = formula.earlyad.hipp.simulation.rs,
                                     "largemodel" = simulation.model.earlyad.tplus.hipp.rs,
                                     "formula_smallmodel" = formula.earlyad.hipp.rs,
                                     "smallmodel" = model.earlyad.tplus.hipp.rs,
                                     "sample_sizes" = seq(100, 800, by=100),
                                     "nsim"=500,
                                     "data" = long.earlyad.tplus.hipp_with_treatment)

earlyadhipp.neuro.tplus.sim.list <- list("formula_largemodel" = formula.earlyad.hipp.simulation.rs,
                                           "largemodel" = simulation.model.earlyad.neuroenriched.tplus.hipp.rs,
                                           "formula_smallmodel" = formula.earlyad.hipp.rs,
                                           "smallmodel" = model.earlyad.neuroenriched.tplus.hipp.rs,
                                           "sample_sizes" = seq(100, 800, by=100),
                                           "nsim"=500,
                                           "data" = long.earlyad.neuroenriched.tplus.hipp_with_treatment)



earlyadhipp.sim.tau.list <- list("formula_largemodel" = formula.earlyad.hipp.simulation.rs,
                                   "largemodel" = simulation.model.earlyad.hipp.tau.rs,
                                   "formula_smallmodel" = formula.earlyad.hipp.rs,
                                   "smallmodel" = model.earlyad.hipp.rs,
                                   "sample_sizes" = seq(100, 800, by=100),
                                   "nsim"=500,
                                   "data" = long.earlyad.hipp_with_treatment)


earlyadhipp.tplus.sim.tau.list <- list("formula_largemodel" = formula.earlyad.hipp.simulation.rs,
                                         "largemodel" = simulation.model.earlyad.tplus.hipp.tau.rs,
                                         "formula_smallmodel" = formula.earlyad.hipp.rs,
                                         "smallmodel" = model.earlyad.tplus.hipp.rs,
                                         "sample_sizes" = seq(100, 800, by=100),
                                         "nsim"=500,
                                         "data" = long.earlyad.tplus.hipp_with_treatment)

earlyadhipp.neuro.tplus.sim.tau.list <- list("formula_largemodel" = formula.earlyad.hipp.simulation.rs,
                                               "largemodel" = simulation.model.earlyad.neuroenriched.tplus.hipp.tau.rs,
                                               "formula_smallmodel" = formula.earlyad.hipp.rs,
                                               "smallmodel" = model.earlyad.neuroenriched.tplus.hipp.rs,
                                               "sample_sizes" = seq(100, 800, by=100),
                                               "nsim"=500,
                                               "data" = long.earlyad.neuroenriched.tplus.hipp_with_treatment)

saveRDS(earlyadhipp.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/earlyadhipp.rds")
saveRDS(earlyadhipp.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/earlyadhipp_tplus.rds")
saveRDS(earlyadhipp.neuro.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/earlyadhipp_neuro_tplus.rds")


saveRDS(earlyadhipp.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/earlyadhipp_tau.rds")
saveRDS(earlyadhipp.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/earlyadhipp_tplus_tau.rds")
saveRDS(earlyadhipp.neuro.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/earlyadhipp_neuro_tplus_tau.rds")




#MPACC #################################### #################################### ####################################

earlyadmpacc.sim.list <- list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs,
                               "largemodel" = simulation.model.earlyad.mpacc.rs,
                               "formula_smallmodel" = formula.earlyad.mpacc.rs,
                               "smallmodel" = model.earlyad.mpacc.rs,
                               "sample_sizes" = seq(100, 800, by=100),
                               "nsim"=500,
                               "data" = long.earlyad.mpacc_with_treatment)


earlyadmpacc.tplus.sim.list <- list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs,
                                     "largemodel" = simulation.model.earlyad.tplus.mpacc.rs,
                                     "formula_smallmodel" = formula.earlyad.mpacc.rs,
                                     "smallmodel" = model.earlyad.tplus.mpacc.rs,
                                     "sample_sizes" = seq(100, 800, by=100),
                                     "nsim"=500,
                                     "data" = long.earlyad.tplus.mpacc_with_treatment)

earlyadmpacc.neuro.tplus.sim.list <- list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs,
                                           "largemodel" = simulation.model.earlyad.neuroenriched.tplus.mpacc.rs,
                                           "formula_smallmodel" = formula.earlyad.mpacc.rs,
                                           "smallmodel" = model.earlyad.neuroenriched.tplus.mpacc.rs,
                                           "sample_sizes" = seq(100, 800, by=100),
                                           "nsim"=500,
                                           "data" = long.earlyad.neuroenriched.tplus.mpacc_with_treatment)



earlyadmpacc.sim.tau.list <- list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs,
                                   "largemodel" = simulation.model.earlyad.mpacc.tau.rs,
                                   "formula_smallmodel" = formula.earlyad.mpacc.rs,
                                   "smallmodel" = model.earlyad.mpacc.rs,
                                   "sample_sizes" = seq(100, 800, by=100),
                                   "nsim"=500,
                                   "data" = long.earlyad.mpacc_with_treatment)


earlyadmpacc.tplus.sim.tau.list <- list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs,
                                         "largemodel" = simulation.model.earlyad.tplus.mpacc.tau.rs,
                                         "formula_smallmodel" = formula.earlyad.mpacc.rs,
                                         "smallmodel" = model.earlyad.tplus.mpacc.rs,
                                         "sample_sizes" = seq(100, 800, by=100),
                                         "nsim"=500,
                                         "data" = long.earlyad.tplus.mpacc_with_treatment)

earlyadmpacc.neuro.tplus.sim.tau.list <- list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs,
                                               "largemodel" = simulation.model.earlyad.neuroenriched.tplus.mpacc.tau.rs,
                                               "formula_smallmodel" = formula.earlyad.mpacc.rs,
                                               "smallmodel" = model.earlyad.neuroenriched.tplus.mpacc.rs,
                                               "sample_sizes" = seq(100, 800, by=100),
                                               "nsim"=500,
                                               "data" = long.earlyad.neuroenriched.tplus.mpacc_with_treatment)


saveRDS(earlyadmpacc.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/earlyadmpacc.rds")
saveRDS(earlyadmpacc.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/earlyadmpacc_tplus.rds")
saveRDS(earlyadmpacc.neuro.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/earlyadmpacc_neuro_tplus.rds")


saveRDS(earlyadmpacc.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/earlyadmpacc_tau.rds")
saveRDS(earlyadmpacc.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/earlyadmpacc_tplus_tau.rds")
saveRDS(earlyadmpacc.neuro.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/earlyadmpacc_neuro_tplus_tau.rds")

#################################### #################################### ####################################

}


checkadas    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/fitted_adas13/earlyad_adas13_unenriched_overall_20.rds")
checkadastau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/fitted_adas13/earlyad_adas13_tplus_overall_20.rds")
checkadasnocopathtau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/fitted_adas13/earlyadadas13_neuro_tplus_overall_20.rds")
checkhipp    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/fitted_hipp/earlyad_hipp_unenriched_overall_20.rds")
checkhipptau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/fitted_hipp/earlyad_hipp_tplus_overall_20.rds")
checkhippnocopathtau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/fitted_hipp/earlyadhipp_neuro_tplus_overall_20.rds")
checkmpacc    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/fitted_mpacc/earlyad_mpacc_unenriched_overall_20.rds")
checkmpacctau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/fitted_mpacc/earlyad_mpacc_tplus_overall_20.rds")
checkmpaccnocopathtau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/fitted_mpacc/earlyadmpacc_neuro_tplus_overall_20.rds")


checkadas.t    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/fitted_adas13/earlyad_adas13_unenriched_tau_20.rds")
checkadastau.t <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/fitted_adas13/earlyad_adas13_tplus_tau_20.rds")
checkhipp.t    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/fitted_hipp/earlyad_hipp_unenriched_tau_20.rds")
checkhipptau.t <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/fitted_hipp/earlyad_hipp_tplus_tau_20.rds")
checkmpacc.t    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/fitted_mpacc/earlyad_mpacc_unenriched_tau_20.rds")
checkmpacctau.t <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/fitted_mpacc/earlyad_mpacc_tplus_tau_20.rds")


all.new.list <- list("checkadas" = checkadas,
                     "checkadastau" = checkadastau,
                     "checkadasnocopathtau" = checkadasnocopathtau,
                     "checkhipp" = checkhipp,
                     "checkhipptau" = checkhipptau,
                     "checkhippnocopathtau" = checkhippnocopathtau,
                     "checkmpacc" = checkmpacc,
                     "checkmpacctau" = checkmpacctau,
                     "checkmpaccnocopathtau" = checkmpaccnocopathtau,
                     "checkadas.t" = checkadas.t,
                     "checkadastau.t" = checkadastau.t,
                     "checkhipp.t" = checkhipp.t,
                     "checkhipptau.t" = checkhipptau.t,
                     "checkmpacc.t" = checkmpacc.t,
                     "checkmpacctau.t" = checkmpacctau.t)


checkall <- Map(function(x){CalculateSampleAtPower_markdown(x$power.data)}, all.new.list)

checkadas_new    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/fitted_adas13/earlyad_adas13_unenriched_tau_20.rds")

