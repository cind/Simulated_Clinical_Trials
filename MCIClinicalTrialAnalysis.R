save.rds <- FALSE

MCIcohorts                     <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/mcicohorts.rds")
MCIcohorts.tplus               <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/mcicohorts_tplus.rds")
MCIcohorts.neuroenriched.tplus <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adni_full_enriched/mcicohorts_neuroenriched_tplus.rds")

#keeping subjects with 3 or more time points
#update descriptive statistics and a few plots in powerpoint as this will remove a few subjects

MCIcohorts                     <- Keep1YearorMore(MCIcohorts)
MCIcohorts.tplus               <- Keep1YearorMore(MCIcohorts.tplus)
MCIcohorts.neuroenriched.tplus <- Keep1YearorMore(MCIcohorts.neuroenriched.tplus)

#Reading data

######################## ADAS13 ######################## 
cs.MCI.adas13     <- MCIcohorts$cs$ADAS13
long.MCI.adas13   <-  MCIcohorts$long$ADAS13

#Tau+
cs.MCI.tplus.adas13    <- MCIcohorts.tplus$cs$ADAS13
long.MCI.tplus.adas13  <-  MCIcohorts.tplus$long$ADAS13

# Neuroenriched Tau+
cs.MCI.neuroenriched.tplus.adas13    <- MCIcohorts.neuroenriched.tplus$cs$ADAS13
long.MCI.neuroenriched.tplus.adas13  <-  MCIcohorts.neuroenriched.tplus$long$ADAS13



######################### HIPPOCAMPUS ######################## 
cs.MCI.hipp    <- MCIcohorts$cs$Imaging
long.MCI.hipp  <-  MCIcohorts$long$Imaging

cs.MCI.hipp$hipp_average    <- (cs.MCI.hipp$ST29SV_harmonized_icv_adj + cs.MCI.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.MCI.hipp$hipp_average  <- (long.MCI.hipp$ST29SV_harmonized_icv_adj + long.MCI.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right

#Tau+
cs.MCI.tplus.hipp           <- MCIcohorts.tplus$cs$Imaging
long.MCI.tplus.hipp         <-  MCIcohorts.tplus$long$Imaging


cs.MCI.tplus.hipp$hipp_average   <- (cs.MCI.tplus.hipp$ST29SV_harmonized_icv_adj + cs.MCI.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.MCI.tplus.hipp$hipp_average <- (long.MCI.tplus.hipp$ST29SV_harmonized_icv_adj + long.MCI.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right

# Neuroenriched Tau+
cs.MCI.neuroenriched.tplus.hipp    <- MCIcohorts.neuroenriched.tplus$cs$Imaging
long.MCI.neuroenriched.tplus.hipp <-  MCIcohorts.neuroenriched.tplus$long$Imaging

cs.MCI.neuroenriched.tplus.hipp$hipp_average   <- (cs.MCI.neuroenriched.tplus.hipp$ST29SV_harmonized_icv_adj + cs.MCI.neuroenriched.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.MCI.neuroenriched.tplus.hipp$hipp_average <- (long.MCI.neuroenriched.tplus.hipp$ST29SV_harmonized_icv_adj + long.MCI.neuroenriched.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right



######################## MPACC######################## 
cs.MCI.mpacc     <- MCIcohorts$cs$MPACC
long.MCI.mpacc   <-  MCIcohorts$long$MPACC


#Tau+
cs.MCI.tplus.mpacc    <- MCIcohorts.tplus$cs$MPACC
long.MCI.tplus.mpacc  <-  MCIcohorts.tplus$long$MPACC


# Neuroenriched Tau+
cs.MCI.neuroenriched.tplus.mpacc    <- MCIcohorts.neuroenriched.tplus$cs$MPACC
long.MCI.neuroenriched.tplus.mpacc  <-  MCIcohorts.neuroenriched.tplus$long$MPACC




#Descriptive statistics tables
######################## ADAS13 ########################
MCI.adas13.desc <- cs.MCI.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                            "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                            "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                            "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.adas13.desc$TauPos_bl     <- factor(MCI.adas13.desc$TauPos_bl)
MCI.adas13.desc$AmyloidPos_bl <- factor(MCI.adas13.desc$AmyloidPos_bl)
colnames(MCI.adas13.desc)     <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                       "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.adas13.desc <- table1(MCI.adas13.desc, splitby = ~ Diagnosis)
MCI.adas13.desc <- MCI.adas13.desc$Table1
MCI.adas13.desc <- MCI.adas13.desc[c(1, 6:9, 11:nrow(MCI.adas13.desc)),]

#Enriched for TAU
MCI.tplus.adas13.desc <- cs.MCI.tplus.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                        "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                        "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.tplus.adas13.desc$TauPos_bl <- factor(MCI.tplus.adas13.desc$TauPos_bl)
MCI.tplus.adas13.desc$AmyloidPos_bl <- factor(MCI.tplus.adas13.desc$AmyloidPos_bl)

colnames(MCI.tplus.adas13.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                         "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                         "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                         "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.tplus.adas13.desc <- table1(MCI.tplus.adas13.desc, splitby = ~ Diagnosis)
MCI.tplus.adas13.desc <- MCI.tplus.adas13.desc$Table1
MCI.tplus.adas13.desc <- MCI.tplus.adas13.desc[c(1, 6:9, 11:nrow(MCI.tplus.adas13.desc)),]


#Neuroenriched Tau+
MCI.neuroenriched.tplus.adas13.desc <- cs.MCI.neuroenriched.tplus.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                    "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                    "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                                                    "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.neuroenriched.tplus.adas13.desc$TauPos_bl <- factor(MCI.neuroenriched.tplus.adas13.desc$TauPos_bl)
MCI.neuroenriched.tplus.adas13.desc$AmyloidPos_bl <- factor(MCI.neuroenriched.tplus.adas13.desc$AmyloidPos_bl)

colnames(MCI.neuroenriched.tplus.adas13.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                       "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.neuroenriched.tplus.adas13.desc <- table1(MCI.neuroenriched.tplus.adas13.desc, splitby = ~ Diagnosis)
MCI.neuroenriched.tplus.adas13.desc <- MCI.neuroenriched.tplus.adas13.desc$Table1
MCI.neuroenriched.tplus.adas13.desc <- MCI.neuroenriched.tplus.adas13.desc[c(1, 6:9, 11:nrow(MCI.neuroenriched.tplus.adas13.desc)),]



#Descriptive statistics tables
######################## Hippocampus ########################
MCI.hipp.desc <- cs.MCI.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                        "hipp_average" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                        "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]


MCI.hipp.desc$TauPos_bl <- factor(MCI.hipp.desc$TauPos_bl)
MCI.hipp.desc$AmyloidPos_bl <- factor(MCI.hipp.desc$AmyloidPos_bl)
colnames(MCI.hipp.desc)   <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                   "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                   "Hippocampus (L/R Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                   "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.hipp.desc <- table1(MCI.hipp.desc, splitby = ~ Diagnosis)
MCI.hipp.desc <- MCI.hipp.desc$Table1
MCI.hipp.desc <- MCI.hipp.desc[c(1, 6:9, 11:nrow(MCI.hipp.desc)),]


#Enriched for TAU
MCI.tplus.hipp.desc <- cs.MCI.tplus.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                    "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                    "hipp_average" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                    "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.tplus.hipp.desc$TauPos_bl <- factor(MCI.tplus.hipp.desc$TauPos_bl)
MCI.tplus.hipp.desc$AmyloidPos_bl <- factor(MCI.tplus.hipp.desc$AmyloidPos_bl)

colnames(MCI.tplus.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                       "Hippocampus (L/R Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.tplus.hipp.desc <- table1(MCI.tplus.hipp.desc, splitby = ~ Diagnosis)
MCI.tplus.hipp.desc <- MCI.tplus.hipp.desc$Table1
MCI.tplus.hipp.desc <- MCI.tplus.hipp.desc[c(1, 6:9, 11:nrow(MCI.tplus.hipp.desc)),]


#Neuroenriched Tau+
MCI.neuroenriched.tplus.hipp.desc <- cs.MCI.neuroenriched.tplus.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                "hipp_average" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                                                "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.neuroenriched.tplus.hipp.desc$TauPos_bl <- factor(MCI.neuroenriched.tplus.hipp.desc$TauPos_bl)
MCI.neuroenriched.tplus.hipp.desc$AmyloidPos_bl <- factor(MCI.neuroenriched.tplus.hipp.desc$AmyloidPos_bl)

colnames(MCI.neuroenriched.tplus.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                     "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                     "Hippocampus (L/R Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                     "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.neuroenriched.tplus.hipp.desc <- table1(MCI.neuroenriched.tplus.hipp.desc, splitby = ~ Diagnosis)
MCI.neuroenriched.tplus.hipp.desc <- MCI.neuroenriched.tplus.hipp.desc$Table1
MCI.neuroenriched.tplus.hipp.desc <- MCI.neuroenriched.tplus.hipp.desc[c(1, 6:9, 11:nrow(MCI.neuroenriched.tplus.hipp.desc)),]



#Descriptive statistics tables
######################## mPACCtrailsB ########################
MCI.mpacc.desc <- cs.MCI.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                          "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                          "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                          "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.mpacc.desc$TauPos_bl <- factor(MCI.mpacc.desc$TauPos_bl)
MCI.mpacc.desc$AmyloidPos_bl <- factor(MCI.mpacc.desc$AmyloidPos_bl)
colnames(MCI.mpacc.desc)   <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                    "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                    "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                    "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.mpacc.desc <- table1(MCI.mpacc.desc, splitby = ~ Diagnosis)
MCI.mpacc.desc <- MCI.mpacc.desc$Table1
MCI.mpacc.desc <- MCI.mpacc.desc[c(1, 6:9, 11:nrow(MCI.mpacc.desc)),]


#Enriched for TAU
MCI.tplus.mpacc.desc <- cs.MCI.tplus.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                      "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                      "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                      "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.tplus.mpacc.desc$TauPos_bl <- factor(MCI.tplus.mpacc.desc$TauPos_bl)
MCI.tplus.mpacc.desc$AmyloidPos_bl <- factor(MCI.tplus.mpacc.desc$AmyloidPos_bl)

colnames(MCI.tplus.mpacc.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                        "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                        "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                        "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.tplus.mpacc.desc <- table1(MCI.tplus.mpacc.desc, splitby = ~ Diagnosis)
MCI.tplus.mpacc.desc <- MCI.tplus.mpacc.desc$Table1
MCI.tplus.mpacc.desc <- MCI.tplus.mpacc.desc[c(1, 6:9, 11:nrow(MCI.tplus.mpacc.desc)),]


#Neuroenriched Tau+
MCI.neuroenriched.tplus.mpacc.desc <- cs.MCI.neuroenriched.tplus.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                  "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                  "ADAS13" ,"mPACCtrailsB_bl", "TauPos_bl", 
                                                                                  "AmyloidPos_bl","CAAPos", "TDP43Pos", "LewyPos")]

MCI.neuroenriched.tplus.mpacc.desc$TauPos_bl <- factor(MCI.neuroenriched.tplus.mpacc.desc$TauPos_bl)
MCI.neuroenriched.tplus.mpacc.desc$AmyloidPos_bl <- factor(MCI.neuroenriched.tplus.mpacc.desc$AmyloidPos_bl)

colnames(MCI.neuroenriched.tplus.mpacc.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                      "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                      "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                      "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

MCI.neuroenriched.tplus.mpacc.desc <- table1(MCI.neuroenriched.tplus.mpacc.desc, splitby = ~ Diagnosis)
MCI.neuroenriched.tplus.mpacc.desc <- MCI.neuroenriched.tplus.mpacc.desc$Table1
MCI.neuroenriched.tplus.mpacc.desc <- MCI.neuroenriched.tplus.mpacc.desc[c(1, 6:9, 11:nrow(MCI.neuroenriched.tplus.mpacc.desc)),]


#fit repeated measures model on longitudinal data
# want observed scores, just use time and dont control for other variables


long.MCI.adas13 <- MMRMTime(long.MCI.adas13)
long.MCI.tplus.adas13 <- MMRMTime(long.MCI.tplus.adas13)
long.MCI.neuroenriched.tplus.adas13 <- MMRMTime(long.MCI.neuroenriched.tplus.adas13)


long.MCI.hipp <-  MMRMTime(long.MCI.hipp)
long.MCI.tplus.hipp <- MMRMTime(long.MCI.tplus.hipp)
long.MCI.neuroenriched.tplus.hipp <- MMRMTime(long.MCI.neuroenriched.tplus.hipp)


long.MCI.mpacc <- MMRMTime(long.MCI.mpacc)
long.MCI.tplus.mpacc <- MMRMTime(long.MCI.tplus.mpacc)
long.MCI.neuroenriched.tplus.mpacc <- MMRMTime(long.MCI.neuroenriched.tplus.mpacc)


######################## ADAS13 ######################## 
mmrm.MCI.adas13        <- gls(ADAS13~factor(new_time_mmrm),
                                  na.action=na.omit, data=long.MCI.adas13,
                                  correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                  weights=nlme::varIdent(form=~1|new_time_mmrm))

#Enriched for TAU
mmrm.MCI.tplus.adas13<- gls(ADAS13~factor(new_time_mmrm),
                                na.action=na.omit, data=long.MCI.tplus.adas13,
                                correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                weights=nlme::varIdent(form=~1|new_time_mmrm))


#Neuroenriched Tau+
mmrm.MCI.neuroenriched.tplus.adas13<- gls(ADAS13~factor(new_time_mmrm),
                                              na.action=na.omit, data=long.MCI.neuroenriched.tplus.adas13,
                                              correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                              weights=nlme::varIdent(form=~1|new_time_mmrm))



######################## HIPPOCAMPUS ######################## 
mmrm.MCI.hipp   <- gls(hipp_average~factor(new_time_mmrm),
                           na.action=na.omit, data=long.MCI.hipp,
                           correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                           weights=nlme::varIdent(form=~1|new_time_mmrm))
mmrm.MCI.hipp.test   <- gls(hipp_average~factor(new_time_mmrm),
                                na.action=na.omit, data=long.MCI.hipp,
                                correlation=nlme::corSymm(form=~1|RID),
                                weights=nlme::varIdent(form=~1|new_time_mmrm))


#Enriched for TAU
mmrm.MCI.tplus.hipp <- gls(hipp_average~factor(new_time_mmrm),
                               na.action=na.omit, data=long.MCI.tplus.hipp,
                               correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                               weights=nlme::varIdent(form=~1|new_time_mmrm))
#Neuroenriched Tau+
mmrm.MCI.neuroenriched.tplus.hipp <- gls(hipp_average~factor(new_time_mmrm),
                                             na.action=na.omit, data=long.MCI.neuroenriched.tplus.hipp,
                                             correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                             weights=nlme::varIdent(form=~1|new_time_mmrm))





######################## MPACC ######################## 
mmrm.MCI.mpacc        <- gls(mPACCtrailsB~factor(new_time_mmrm),
                                 na.action=na.omit, data=long.MCI.mpacc,
                                 correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                 weights=nlme::varIdent(form=~1|new_time_mmrm))


#Enriched for TAU
mmrm.MCI.tplus.mpacc<- gls(mPACCtrailsB~factor(new_time_mmrm),
                               na.action=na.omit, data=long.MCI.tplus.mpacc,
                               correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                               weights=nlme::varIdent(form=~1|new_time_mmrm))

#Neuroenriched Tau+
mmrm.MCI.neuroenriched.tplus.mpacc<- gls(mPACCtrailsB~factor(new_time_mmrm),
                                             na.action=na.omit, data=long.MCI.neuroenriched.tplus.mpacc,
                                             correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                             weights=nlme::varIdent(form=~1|new_time_mmrm))



#get coefficients and standard errors

######################## ADAS13 ######################## 
coefs.MCI.adas13     <- unname(coef(summary(mmrm.MCI.adas13))[, "Value"])
stderror.MCI.adas13  <- unname(coef(summary(mmrm.MCI.adas13))[, "Std.Error"])

#Enriched for TAU
coefs.MCI.tplus.adas13     <- unname(coef(summary(mmrm.MCI.tplus.adas13))[, "Value"])
stderror.MCI.tplus.adas13  <- unname(coef(summary(mmrm.MCI.tplus.adas13))[, "Std.Error"])


#Neuroenriched Tau+
coefs.MCI.neuroenriched.tplus.adas13     <- unname(coef(summary(mmrm.MCI.neuroenriched.tplus.adas13))[, "Value"])
stderror.MCI.neuroenriched.tplus.adas13  <- unname(coef(summary(mmrm.MCI.neuroenriched.tplus.adas13))[, "Std.Error"])


######################## HIPPOCAMPUS ######################## 
coefs.MCI.hipp     <- unname(coef(summary(mmrm.MCI.hipp))[, "Value"])
stderror.MCI.hipp  <- unname(coef(summary(mmrm.MCI.hipp))[, "Std.Error"])

#Enriched for TAU
coefs.MCI.tplus.hipp     <- unname(coef(summary(mmrm.MCI.tplus.hipp))[, "Value"])
stderror.MCI.tplus.hipp  <- unname(coef(summary(mmrm.MCI.tplus.hipp))[, "Std.Error"])



#Neuroenriched Tau+
coefs.MCI.neuroenriched.tplus.hipp     <- unname(coef(summary(mmrm.MCI.neuroenriched.tplus.hipp))[, "Value"])
stderror.MCI.neuroenriched.tplus.hipp  <- unname(coef(summary(mmrm.MCI.neuroenriched.tplus.hipp))[, "Std.Error"])



######################## MPACC ######################## 
coefs.MCI.mpacc     <- unname(coef(summary(mmrm.MCI.mpacc))[, "Value"])
stderror.MCI.mpacc  <- unname(coef(summary(mmrm.MCI.mpacc))[, "Std.Error"])

#Enriched for TAU
coefs.MCI.tplus.mpacc     <- unname(coef(summary(mmrm.MCI.tplus.mpacc))[, "Value"])
stderror.MCI.tplus.mpacc  <- unname(coef(summary(mmrm.MCI.tplus.mpacc))[, "Std.Error"])


#Neuroenriched Tau+
coefs.MCI.neuroenriched.tplus.mpacc     <- unname(coef(summary(mmrm.MCI.neuroenriched.tplus.mpacc))[, "Value"])
stderror.MCI.neuroenriched.tplus.mpacc  <- unname(coef(summary(mmrm.MCI.neuroenriched.tplus.mpacc))[, "Std.Error"])



#create dataset for plotting
######################## ADAS13 ######################## 
MCI.adas13.plotting.data <- data.frame("coefs" = coefs.MCI.adas13,
                                           "stderr" = stderror.MCI.adas13)
MCI.adas13.plotting.data$ci_low     <-  -1.96* MCI.adas13.plotting.data$stderr
MCI.adas13.plotting.data$ci_hi      <-   1.96* MCI.adas13.plotting.data$stderr
MCI.adas13.plotting.data$means      <-  MCI.adas13.plotting.data$coefs
MCI.adas13.plotting.data$means[2:5] <-  MCI.adas13.plotting.data$means[2:5] + MCI.adas13.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.adas13.plotting.data$ci_low     <-  MCI.adas13.plotting.data$ci_low+ MCI.adas13.plotting.data$means
MCI.adas13.plotting.data$ci_hi      <-  MCI.adas13.plotting.data$ci_hi+ MCI.adas13.plotting.data$means
MCI.adas13.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.adas13.plotting.data$Enrichment <- "No Enrichment"

#Enriched for TAU
MCI.tplus.adas13.plotting.data <- data.frame("coefs" = coefs.MCI.tplus.adas13,
                                                 "stderr" = stderror.MCI.tplus.adas13)
MCI.tplus.adas13.plotting.data$ci_low     <-  -1.96* MCI.tplus.adas13.plotting.data$stderr
MCI.tplus.adas13.plotting.data$ci_hi      <-   1.96* MCI.tplus.adas13.plotting.data$stderr
MCI.tplus.adas13.plotting.data$means      <-  MCI.tplus.adas13.plotting.data$coefs
MCI.tplus.adas13.plotting.data$means[2:5] <-  MCI.tplus.adas13.plotting.data$means[2:5] + MCI.tplus.adas13.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.tplus.adas13.plotting.data$ci_low     <-  MCI.tplus.adas13.plotting.data$ci_low+ MCI.tplus.adas13.plotting.data$means
MCI.tplus.adas13.plotting.data$ci_hi      <-  MCI.tplus.adas13.plotting.data$ci_hi+ MCI.tplus.adas13.plotting.data$means
MCI.tplus.adas13.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.tplus.adas13.plotting.data$Enrichment <- "Tau+"

#Neuroenriched Tau+
MCI.neuroenriched.tplus.adas13.plotting.data <- data.frame("coefs" = coefs.MCI.neuroenriched.tplus.adas13,
                                                               "stderr" = stderror.MCI.neuroenriched.tplus.adas13)
MCI.neuroenriched.tplus.adas13.plotting.data$ci_low     <-  -1.96* MCI.neuroenriched.tplus.adas13.plotting.data$stderr
MCI.neuroenriched.tplus.adas13.plotting.data$ci_hi      <-   1.96* MCI.neuroenriched.tplus.adas13.plotting.data$stderr
MCI.neuroenriched.tplus.adas13.plotting.data$means      <-  MCI.neuroenriched.tplus.adas13.plotting.data$coefs
MCI.neuroenriched.tplus.adas13.plotting.data$means[2:5] <-  MCI.neuroenriched.tplus.adas13.plotting.data$means[2:5] + MCI.neuroenriched.tplus.adas13.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.neuroenriched.tplus.adas13.plotting.data$ci_low     <-  MCI.neuroenriched.tplus.adas13.plotting.data$ci_low+ MCI.neuroenriched.tplus.adas13.plotting.data$means
MCI.neuroenriched.tplus.adas13.plotting.data$ci_hi      <-  MCI.neuroenriched.tplus.adas13.plotting.data$ci_hi+ MCI.neuroenriched.tplus.adas13.plotting.data$means
MCI.neuroenriched.tplus.adas13.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.neuroenriched.tplus.adas13.plotting.data$Enrichment <- "No Copathologies Tau+"




######################## HIPPOCAMPUS ######################## 
MCI.hipp.plotting.data <- data.frame("coefs" = coefs.MCI.hipp,
                                         "stderr" = stderror.MCI.hipp)
MCI.hipp.plotting.data$ci_low     <-  -1.96* MCI.hipp.plotting.data$stderr
MCI.hipp.plotting.data$ci_hi      <-   1.96* MCI.hipp.plotting.data$stderr
MCI.hipp.plotting.data$means      <-  MCI.hipp.plotting.data$coefs
MCI.hipp.plotting.data$means[2:5] <-  MCI.hipp.plotting.data$means[2:5] + MCI.hipp.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.hipp.plotting.data$ci_low     <-  MCI.hipp.plotting.data$ci_low+ MCI.hipp.plotting.data$means
MCI.hipp.plotting.data$ci_hi      <-  MCI.hipp.plotting.data$ci_hi+ MCI.hipp.plotting.data$means

MCI.hipp.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.hipp.plotting.data$Enrichment <- "No Enrichment"

#Enriched for TAU
MCI.tplus.hipp.plotting.data <- data.frame("coefs" = coefs.MCI.tplus.hipp,
                                               "stderr" = stderror.MCI.tplus.hipp)
MCI.tplus.hipp.plotting.data$ci_low     <-  -1.96* MCI.tplus.hipp.plotting.data$stderr
MCI.tplus.hipp.plotting.data$ci_hi      <-   1.96* MCI.tplus.hipp.plotting.data$stderr
MCI.tplus.hipp.plotting.data$means      <-  MCI.tplus.hipp.plotting.data$coefs
MCI.tplus.hipp.plotting.data$means[2:5] <-  MCI.tplus.hipp.plotting.data$means[2:5] + MCI.tplus.hipp.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.tplus.hipp.plotting.data$ci_low     <-  MCI.tplus.hipp.plotting.data$ci_low+ MCI.tplus.hipp.plotting.data$means
MCI.tplus.hipp.plotting.data$ci_hi      <-  MCI.tplus.hipp.plotting.data$ci_hi+ MCI.tplus.hipp.plotting.data$means
MCI.tplus.hipp.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.tplus.hipp.plotting.data$Enrichment <- "Tau+"



#Neuroenriched Tau+
MCI.neuroenriched.tplus.hipp.plotting.data <- data.frame("coefs" = coefs.MCI.neuroenriched.tplus.hipp,
                                                             "stderr" = stderror.MCI.neuroenriched.tplus.hipp)
MCI.neuroenriched.tplus.hipp.plotting.data$ci_low     <-  -1.96* MCI.neuroenriched.tplus.hipp.plotting.data$stderr
MCI.neuroenriched.tplus.hipp.plotting.data$ci_hi      <-   1.96* MCI.neuroenriched.tplus.hipp.plotting.data$stderr
MCI.neuroenriched.tplus.hipp.plotting.data$means      <-  MCI.neuroenriched.tplus.hipp.plotting.data$coefs
MCI.neuroenriched.tplus.hipp.plotting.data$means[2:5] <-  MCI.neuroenriched.tplus.hipp.plotting.data$means[2:5] + MCI.neuroenriched.tplus.hipp.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.neuroenriched.tplus.hipp.plotting.data$ci_low     <-  MCI.neuroenriched.tplus.hipp.plotting.data$ci_low+ MCI.neuroenriched.tplus.hipp.plotting.data$means
MCI.neuroenriched.tplus.hipp.plotting.data$ci_hi      <-  MCI.neuroenriched.tplus.hipp.plotting.data$ci_hi+ MCI.neuroenriched.tplus.hipp.plotting.data$means
MCI.neuroenriched.tplus.hipp.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.neuroenriched.tplus.hipp.plotting.data$Enrichment <- "No Copathologies Tau+"


######################## MPACC ######################## 
MCI.mpacc.plotting.data <- data.frame("coefs" = coefs.MCI.mpacc,
                                          "stderr" = stderror.MCI.mpacc)
MCI.mpacc.plotting.data$ci_low     <-  -1.96* MCI.mpacc.plotting.data$stderr
MCI.mpacc.plotting.data$ci_hi      <-   1.96* MCI.mpacc.plotting.data$stderr
MCI.mpacc.plotting.data$means      <-  MCI.mpacc.plotting.data$coefs
MCI.mpacc.plotting.data$means[2:5] <-  MCI.mpacc.plotting.data$means[2:5] + MCI.mpacc.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.mpacc.plotting.data$ci_low     <-  MCI.mpacc.plotting.data$ci_low+ MCI.mpacc.plotting.data$means
MCI.mpacc.plotting.data$ci_hi      <-  MCI.mpacc.plotting.data$ci_hi+ MCI.mpacc.plotting.data$means
MCI.mpacc.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.mpacc.plotting.data$Enrichment <- "No Enrichment"

#Enriched for TAU
MCI.tplus.mpacc.plotting.data <- data.frame("coefs" = coefs.MCI.tplus.mpacc,
                                                "stderr" = stderror.MCI.tplus.mpacc)
MCI.tplus.mpacc.plotting.data$ci_low     <-  -1.96* MCI.tplus.mpacc.plotting.data$stderr
MCI.tplus.mpacc.plotting.data$ci_hi      <-   1.96* MCI.tplus.mpacc.plotting.data$stderr
MCI.tplus.mpacc.plotting.data$means      <-  MCI.tplus.mpacc.plotting.data$coefs
MCI.tplus.mpacc.plotting.data$means[2:5] <-  MCI.tplus.mpacc.plotting.data$means[2:5] + MCI.tplus.mpacc.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.tplus.mpacc.plotting.data$ci_low     <-  MCI.tplus.mpacc.plotting.data$ci_low+ MCI.tplus.mpacc.plotting.data$means
MCI.tplus.mpacc.plotting.data$ci_hi      <-  MCI.tplus.mpacc.plotting.data$ci_hi+ MCI.tplus.mpacc.plotting.data$means
MCI.tplus.mpacc.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.tplus.mpacc.plotting.data$Enrichment <- "Tau+"


#Neuroenriched Tau+
MCI.neuroenriched.tplus.mpacc.plotting.data <- data.frame("coefs" = coefs.MCI.neuroenriched.tplus.mpacc,
                                                              "stderr" = stderror.MCI.neuroenriched.tplus.mpacc)
MCI.neuroenriched.tplus.mpacc.plotting.data$ci_low     <-  -1.96* MCI.neuroenriched.tplus.mpacc.plotting.data$stderr
MCI.neuroenriched.tplus.mpacc.plotting.data$ci_hi      <-   1.96* MCI.neuroenriched.tplus.mpacc.plotting.data$stderr
MCI.neuroenriched.tplus.mpacc.plotting.data$means      <-  MCI.neuroenriched.tplus.mpacc.plotting.data$coefs
MCI.neuroenriched.tplus.mpacc.plotting.data$means[2:5] <-  MCI.neuroenriched.tplus.mpacc.plotting.data$means[2:5] + MCI.neuroenriched.tplus.mpacc.plotting.data$coefs[1] #add intercept to each estimate for plotting
MCI.neuroenriched.tplus.mpacc.plotting.data$ci_low     <-  MCI.neuroenriched.tplus.mpacc.plotting.data$ci_low+ MCI.neuroenriched.tplus.mpacc.plotting.data$means
MCI.neuroenriched.tplus.mpacc.plotting.data$ci_hi      <-  MCI.neuroenriched.tplus.mpacc.plotting.data$ci_hi+ MCI.neuroenriched.tplus.mpacc.plotting.data$means
MCI.neuroenriched.tplus.mpacc.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
MCI.neuroenriched.tplus.mpacc.plotting.data$Enrichment <- "No Copathologies Tau+"



## Plots
# Combine non-enriched and tau enriched data for plotting

######################## ADAS13 ######################## 
MCI.adas13.plotting.data <- rbind(MCI.adas13.plotting.data,
                                      MCI.tplus.adas13.plotting.data,
                                      MCI.neuroenriched.tplus.adas13.plotting.data)
MCI.adas13.mmrm.plot <- ggplot(MCI.adas13.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
MCI.adas13.mmrm.plot <- MCI.adas13.mmrm.plot + xlab("Time (Years)") + ylab("Mean ADAS-13 (95% CI)") + ylim(10, 30)

######################## HIPPOCAMPUS ######################## 
MCI.hipp.plotting.data <- rbind(MCI.hipp.plotting.data,
                                    MCI.tplus.hipp.plotting.data,
                                    MCI.neuroenriched.tplus.hipp.plotting.data)
MCI.hipp.mmrm.plot <- ggplot(MCI.hipp.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
MCI.hipp.mmrm.plot <- MCI.hipp.mmrm.plot + xlab("Time (Years)") + ylab("Mean Hippocampus (95% CI)") + ylim(2700, 3800)

######################## MPACC ######################## 
MCI.mpacc.plotting.data <- rbind(MCI.mpacc.plotting.data,
                                     MCI.tplus.mpacc.plotting.data,
                                     MCI.neuroenriched.tplus.mpacc.plotting.data)
MCI.mpacc.mmrm.plot <- ggplot(MCI.mpacc.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
MCI.mpacc.mmrm.plot <- MCI.mpacc.mmrm.plot + xlab("Time (Years)") + ylab("Mean mPACCtrailsB (95% CI)") 


######################## ADAS13 ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.MCI.adas13 <- StratifyContVar(cs.MCI.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.MCI.tplus.adas13 <- StratifyContVar(cs.MCI.tplus.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.MCI.neuroenriched.tplus.adas13 <- StratifyContVar(cs.MCI.neuroenriched.tplus.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))




#Assign treatment and placebo by blocks
long.MCI.adas13_with_treatment <- RandomizeTreatment2(cs.MCI.adas13, long.MCI.adas13, no.prop = TRUE)
long.MCI.tplus.adas13_with_treatment <- RandomizeTreatment2(cs.MCI.tplus.adas13, long.MCI.tplus.adas13, no.prop = TRUE)
long.MCI.neuroenriched.tplus.adas13_with_treatment <- RandomizeTreatment2(cs.MCI.neuroenriched.tplus.adas13, long.MCI.neuroenriched.tplus.adas13, no.prop = TRUE)

######################## HIPPOCAMPUS ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.MCI.hipp <- StratifyContVar(cs.MCI.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.MCI.tplus.hipp <- StratifyContVar(cs.MCI.tplus.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.MCI.neuroenriched.tplus.hipp <- StratifyContVar(cs.MCI.neuroenriched.tplus.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))


#Assign treatment and placebo by blocks
long.MCI.hipp_with_treatment <- RandomizeTreatment2(cs.MCI.hipp, long.MCI.hipp, no.prop = TRUE)
long.MCI.tplus.hipp_with_treatment <- RandomizeTreatment2(cs.MCI.tplus.hipp, long.MCI.tplus.hipp, no.prop = TRUE)
long.MCI.neuroenriched.tplus.hipp_with_treatment <- RandomizeTreatment2(cs.MCI.neuroenriched.tplus.hipp, long.MCI.neuroenriched.tplus.hipp, no.prop = TRUE)


######################## MPACC ######################## 

# Categorize continuous variables for block randomization
cs.MCI.mpacc <- StratifyContVar(cs.MCI.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.MCI.tplus.mpacc <- StratifyContVar(cs.MCI.tplus.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.MCI.neuroenriched.tplus.mpacc <- StratifyContVar(cs.MCI.neuroenriched.tplus.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))


#Assign treatment and placebo by blocks
long.MCI.mpacc_with_treatment <- RandomizeTreatment2(cs.MCI.mpacc, long.MCI.mpacc, no.prop = TRUE)
long.MCI.tplus.mpacc_with_treatment <- RandomizeTreatment2(cs.MCI.tplus.mpacc, long.MCI.tplus.mpacc, no.prop = TRUE)
long.MCI.neuroenriched.tplus.mpacc_with_treatment <- RandomizeTreatment2(cs.MCI.neuroenriched.tplus.mpacc, long.MCI.neuroenriched.tplus.mpacc, no.prop = TRUE)



######################## ADAS13 ######################## 

formula.MCI.adas13.rs                              <- "ADAS13~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + (1 + new_time|RID)"
long.MCI.adas13_with_treatment                     <- StratifyContinuous(long.MCI.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.MCI.tplus.adas13_with_treatment               <- StratifyContinuous(long.MCI.tplus.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.MCI.neuroenriched.tplus.adas13_with_treatment <- StratifyContinuous(long.MCI.neuroenriched.tplus.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))



model.MCI.adas13.rs                     <- MapLmer(newdata = long.MCI.adas13_with_treatment,
                                                       formula.model = formula.MCI.adas13.rs)

model.MCI.tplus.adas13.rs               <- MapLmer(newdata = long.MCI.tplus.adas13_with_treatment,
                                                       formula.model = formula.MCI.adas13.rs)
model.MCI.neuroenriched.tplus.adas13.rs <- MapLmer(newdata = long.MCI.neuroenriched.tplus.adas13_with_treatment,
                                                       formula.model = formula.MCI.adas13.rs)

formula.MCI.adas13.simulation.rs                   <- "ADAS13 ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + (1 + new_time|RID)"
simulation.model.MCI.adas13.rs                     <- BuildSimulationModelNoPath(model.MCI.adas13.rs, formula.MCI.adas13.simulation.rs, long.MCI.adas13_with_treatment, "not-controlled")
simulation.model.MCI.tplus.adas13.rs               <- BuildSimulationModelNoPath(model.MCI.tplus.adas13.rs, formula.MCI.adas13.simulation.rs, long.MCI.tplus.adas13_with_treatment, "not-controlled")
simulation.model.MCI.neuroenriched.tplus.adas13.rs <- BuildSimulationModelNoPath(model.MCI.neuroenriched.tplus.adas13.rs, formula.MCI.adas13.simulation.rs, long.MCI.neuroenriched.tplus.adas13_with_treatment, "not-controlled")

relcontr.MCI.adas13                            <- GetRelContributions(model.MCI.adas13.rs, long.MCI.adas13_with_treatment)
relcontr.MCI.tplus.adas13                      <- GetRelContributions(model.MCI.tplus.adas13.rs, long.MCI.tplus.adas13_with_treatment)
relcontr.MCI.neuroenriched.tplus.adas13        <- GetRelContributions(model.MCI.neuroenriched.tplus.adas13.rs, long.MCI.neuroenriched.tplus.adas13_with_treatment)


relcontrcount.MCI.adas13        <- BuildNeuroCountPlot(relcontr.MCI.adas13)
relcontrcount.MCI.tplus.adas13  <- BuildNeuroCountPlot(relcontr.MCI.tplus.adas13)





simulation.model.MCI.adas13.tau.rs                     <- BuildSimulationModelNoPath(model.MCI.adas13.rs, 
                                                                                         formula.MCI.adas13.simulation.rs, 
                                                                                         long.MCI.adas13_with_treatment,  
                                                                                         relcontr.MCI.adas13$Rates_Data$Mean_Rate[[4]])

simulation.model.MCI.tplus.adas13.tau.rs               <- BuildSimulationModelNoPath(model.MCI.tplus.adas13.rs, 
                                                                                         formula.MCI.adas13.simulation.rs, 
                                                                                         long.MCI.tplus.adas13_with_treatment,  
                                                                                         relcontr.MCI.tplus.adas13$Rates_Data$Mean_Rate[[4]])

simulation.model.MCI.neuroenriched.tplus.adas13.tau.rs <- BuildSimulationModelNoPath(model.MCI.neuroenriched.tplus.adas13.rs, 
                                                                                         formula.MCI.adas13.simulation.rs, 
                                                                                         long.MCI.neuroenriched.tplus.adas13_with_treatment, 
                                                                                         relcontr.MCI.neuroenriched.tplus.adas13$Rates_Data$Mean_Rate[[4]])


dpms.MCI.adas13 <- DPMPlots(list("No Enrichment" = simulation.model.MCI.adas13.rs,
                                     "Tau+" = simulation.model.MCI.tplus.adas13.rs,
                                     "No Copathologies Tau+" = simulation.model.MCI.neuroenriched.tplus.adas13.rs), ylab="ADAS13", ylim.low = 19, ylim.high = 35.5)

dpms.MCI.adas13.tau <- DPMPlots(list("No Enrichment" = simulation.model.MCI.adas13.tau.rs,
                                         "Tau+" = simulation.model.MCI.tplus.adas13.tau.rs,
                                         "No Copathologies Tau+" = simulation.model.MCI.neuroenriched.tplus.adas13.tau.rs), ylab="ADAS13", ylim.low = 19, ylim.high = 35.5)

######################## Hippcampus ######################## 
formula.MCI.hipp.rs                              <- "hipp_average~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl + (1 + new_time|RID)"
long.MCI.hipp_with_treatment                     <- StratifyContinuous(long.MCI.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.MCI.tplus.hipp_with_treatment               <- StratifyContinuous(long.MCI.tplus.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.MCI.neuroenriched.tplus.hipp_with_treatment <- StratifyContinuous(long.MCI.neuroenriched.tplus.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))

model.MCI.hipp.rs                     <- MapLmer(newdata = long.MCI.hipp_with_treatment,
                                                     formula.model = formula.MCI.hipp.rs)


model.MCI.tplus.hipp.rs               <- MapLmer(newdata = long.MCI.tplus.hipp_with_treatment,
                                                     formula.model = formula.MCI.hipp.rs)

model.MCI.neuroenriched.tplus.hipp.rs <- MapLmer(newdata = long.MCI.neuroenriched.tplus.hipp_with_treatment,
                                                     formula.model = formula.MCI.hipp.rs)

formula.MCI.hipp.simulation.rs                   <- "hipp_average ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl +  (1 + new_time|RID)"
simulation.model.MCI.hipp.rs                     <- BuildSimulationModelNoPath(model.MCI.hipp.rs, formula.MCI.hipp.simulation.rs, long.MCI.hipp_with_treatment, "not-controlled")
simulation.model.MCI.tplus.hipp.rs               <- BuildSimulationModelNoPath(model.MCI.tplus.hipp.rs, formula.MCI.hipp.simulation.rs, long.MCI.tplus.hipp_with_treatment, "not-controlled")
simulation.model.MCI.neuroenriched.tplus.hipp.rs <- BuildSimulationModelNoPath(model.MCI.neuroenriched.tplus.hipp.rs, formula.MCI.hipp.simulation.rs, long.MCI.neuroenriched.tplus.hipp_with_treatment, "not-controlled")


relcontr.MCI.hipp                            <- GetRelContributions(model.MCI.hipp.rs, long.MCI.hipp_with_treatment)
relcontr.MCI.tplus.hipp                      <- GetRelContributions(model.MCI.tplus.hipp.rs, long.MCI.tplus.hipp_with_treatment)
relcontr.MCI.neuroenriched.tplus.hipp        <- GetRelContributions(model.MCI.neuroenriched.tplus.hipp.rs, long.MCI.neuroenriched.tplus.hipp_with_treatment)


relcontrcount.MCI.hipp        <- BuildNeuroCountPlot(relcontr.MCI.hipp)
relcontrcount.MCI.tplus.hipp  <- BuildNeuroCountPlot(relcontr.MCI.tplus.hipp)



simulation.model.MCI.hipp.tau.rs                     <- BuildSimulationModelNoPath(model.MCI.hipp.rs, 
                                                                                       formula.MCI.hipp.simulation.rs, 
                                                                                       long.MCI.hipp_with_treatment,  
                                                                                       relcontr.MCI.hipp$Rates_Data$Mean_Rate[[4]])

simulation.model.MCI.tplus.hipp.tau.rs               <- BuildSimulationModelNoPath(model.MCI.tplus.hipp.rs, 
                                                                                       formula.MCI.hipp.simulation.rs, 
                                                                                       long.MCI.tplus.hipp_with_treatment,  
                                                                                       relcontr.MCI.tplus.hipp$Rates_Data$Mean_Rate[[4]])

simulation.model.MCI.neuroenriched.tplus.hipp.tau.rs <- BuildSimulationModelNoPath(model.MCI.neuroenriched.tplus.hipp.rs, 
                                                                                       formula.MCI.hipp.simulation.rs, 
                                                                                       long.MCI.neuroenriched.tplus.hipp_with_treatment, 
                                                                                       relcontr.MCI.neuroenriched.tplus.hipp$Rates_Data$Mean_Rate[[4]])


dpms.MCI.hipp <- DPMPlots(list("No Enrichment" = simulation.model.MCI.hipp.rs,
                                   "Tau+" = simulation.model.MCI.tplus.hipp.rs,
                                   "No Copathologies Tau+" = simulation.model.MCI.neuroenriched.tplus.hipp.rs), ylab="Hippocampus", ylim.low = 2800, ylim.high = 4000)

dpms.MCI.hipp.tau <- DPMPlots(list("No Enrichment" = simulation.model.MCI.hipp.tau.rs,
                                       "Tau+" = simulation.model.MCI.tplus.hipp.tau.rs,
                                       "No Copathologies Tau+" = simulation.model.MCI.neuroenriched.tplus.hipp.tau.rs), ylab="Hippocampus", ylim.low = 2800, ylim.high = 4000)



######################## mPACCtrailsB ######################## 
formula.MCI.mpacc.rs                              <- "mPACCtrailsB~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl +  (1 + new_time|RID)"
long.MCI.mpacc_with_treatment                     <- StratifyContinuous(long.MCI.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.MCI.tplus.mpacc_with_treatment               <- StratifyContinuous(long.MCI.tplus.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.MCI.neuroenriched.tplus.mpacc_with_treatment <- StratifyContinuous(long.MCI.neuroenriched.tplus.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))


model.MCI.mpacc.rs                     <- MapLmer(newdata = long.MCI.mpacc_with_treatment,
                                                      formula.model = formula.MCI.mpacc.rs)

model.MCI.tplus.mpacc.rs               <- MapLmer(newdata = long.MCI.tplus.mpacc_with_treatment,
                                                      formula.model = formula.MCI.mpacc.rs)
model.MCI.neuroenriched.tplus.mpacc.rs <- MapLmer(newdata = long.MCI.neuroenriched.tplus.mpacc_with_treatment,
                                                      formula.model = formula.MCI.mpacc.rs)

formula.MCI.mpacc.simulation.rs                   <- "mPACCtrailsB ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER + MMSE_bl +  (1 + new_time|RID)"
simulation.model.MCI.mpacc.rs                     <- BuildSimulationModelNoPath(model.MCI.mpacc.rs, formula.MCI.mpacc.simulation.rs, long.MCI.mpacc_with_treatment, "not-controlled")
simulation.model.MCI.tplus.mpacc.rs               <- BuildSimulationModelNoPath(model.MCI.tplus.mpacc.rs, formula.MCI.mpacc.simulation.rs, long.MCI.tplus.mpacc_with_treatment, "not-controlled")
simulation.model.MCI.neuroenriched.tplus.mpacc.rs <- BuildSimulationModelNoPath(model.MCI.neuroenriched.tplus.mpacc.rs, formula.MCI.mpacc.simulation.rs, long.MCI.neuroenriched.tplus.mpacc_with_treatment, "not-controlled")


relcontr.MCI.mpacc                            <- GetRelContributions(model.MCI.mpacc.rs, long.MCI.mpacc_with_treatment)
relcontr.MCI.tplus.mpacc                      <- GetRelContributions(model.MCI.tplus.mpacc.rs, long.MCI.tplus.mpacc_with_treatment)
relcontr.MCI.neuroenriched.tplus.mpacc        <- GetRelContributions(model.MCI.neuroenriched.tplus.mpacc.rs, long.MCI.neuroenriched.tplus.mpacc_with_treatment)


relcontrcount.MCI.mpacc        <- BuildNeuroCountPlot(relcontr.MCI.mpacc)
relcontrcount.MCI.tplus.mpacc  <- BuildNeuroCountPlot(relcontr.MCI.tplus.mpacc)



simulation.model.MCI.mpacc.tau.rs                     <- BuildSimulationModelNoPath(model.MCI.mpacc.rs, 
                                                                                        formula.MCI.mpacc.simulation.rs, 
                                                                                        long.MCI.mpacc_with_treatment,  
                                                                                        relcontr.MCI.mpacc$Rates_Data$Mean_Rate[[4]])

simulation.model.MCI.tplus.mpacc.tau.rs               <- BuildSimulationModelNoPath(model.MCI.tplus.mpacc.rs, 
                                                                                        formula.MCI.mpacc.simulation.rs, 
                                                                                        long.MCI.tplus.mpacc_with_treatment,  
                                                                                        relcontr.MCI.tplus.mpacc$Rates_Data$Mean_Rate[[4]])

simulation.model.MCI.neuroenriched.tplus.mpacc.tau.rs <- BuildSimulationModelNoPath(model.MCI.neuroenriched.tplus.mpacc.rs, 
                                                                                        formula.MCI.mpacc.simulation.rs, 
                                                                                        long.MCI.neuroenriched.tplus.mpacc_with_treatment, 
                                                                                        relcontr.MCI.neuroenriched.tplus.mpacc$Rates_Data$Mean_Rate[[4]])



dpms.MCI.mpacc <- DPMPlots(list("No Enrichment" = simulation.model.MCI.mpacc.rs,
                                    "Tau+" = simulation.model.MCI.tplus.mpacc.rs,
                                    "No Copathologies Tau+" = simulation.model.MCI.neuroenriched.tplus.mpacc.rs), ylab="mPACCtrailsB", ylim.low = -19, ylim.high = -7)

dpms.MCI.mpacc.tau <- DPMPlots(list("No Enrichment" = simulation.model.MCI.mpacc.tau.rs,
                                        "Tau+" = simulation.model.MCI.tplus.mpacc.tau.rs,
                                        "No Copathologies Tau+" = simulation.model.MCI.neuroenriched.tplus.mpacc.tau.rs), ylab="mPACCtrailsB", ylim.low = -19, ylim.high = -7)


if(save.rds) {
  
  
  
  #ADAS#################################### #################################### ####################################
  
  MCIadas13.sim.list <- list("formula_largemodel" = formula.MCI.adas13.simulation.rs,
                                 "largemodel" = simulation.model.MCI.adas13.rs,
                                 "formula_smallmodel" = formula.MCI.adas13.rs,
                                 "smallmodel" = model.MCI.adas13.rs,
                                 "sample_sizes" = seq(100, 800, by=100),
                                 "nsim"=500,
                                 "data" = long.MCI.adas13_with_treatment)
  
  
  MCIadas13.tplus.sim.list <- list("formula_largemodel" = formula.MCI.adas13.simulation.rs,
                                       "largemodel" = simulation.model.MCI.tplus.adas13.rs,
                                       "formula_smallmodel" = formula.MCI.adas13.rs,
                                       "smallmodel" = model.MCI.tplus.adas13.rs,
                                       "sample_sizes" = seq(100, 800, by=100),
                                       "nsim"=500,
                                       "data" = long.MCI.tplus.adas13_with_treatment)
  
  MCIadas13.neuro.tplus.sim.list <- list("formula_largemodel" = formula.MCI.adas13.simulation.rs,
                                             "largemodel" = simulation.model.MCI.neuroenriched.tplus.adas13.rs,
                                             "formula_smallmodel" = formula.MCI.adas13.rs,
                                             "smallmodel" = model.MCI.neuroenriched.tplus.adas13.rs,
                                             "sample_sizes" = seq(100, 800, by=100),
                                             "nsim"=500,
                                             "data" = long.MCI.neuroenriched.tplus.adas13_with_treatment)
  
  
  
  MCIadas13.sim.tau.list <- list("formula_largemodel" = formula.MCI.adas13.simulation.rs,
                                     "largemodel" = simulation.model.MCI.adas13.tau.rs,
                                     "formula_smallmodel" = formula.MCI.adas13.rs,
                                     "smallmodel" = model.MCI.adas13.rs,
                                     "sample_sizes" = seq(100, 800, by=100),
                                     "nsim"=500,
                                     "data" = long.MCI.adas13_with_treatment)
  
  
  MCIadas13.tplus.sim.tau.list <- list("formula_largemodel" = formula.MCI.adas13.simulation.rs,
                                           "largemodel" = simulation.model.MCI.tplus.adas13.tau.rs,
                                           "formula_smallmodel" = formula.MCI.adas13.rs,
                                           "smallmodel" = model.MCI.tplus.adas13.rs,
                                           "sample_sizes" = seq(100, 800, by=100),
                                           "nsim"=500,
                                           "data" = long.MCI.tplus.adas13_with_treatment)
  
  MCIadas13.neuro.tplus.sim.tau.list <- list("formula_largemodel" = formula.MCI.adas13.simulation.rs,
                                                 "largemodel" = simulation.model.MCI.neuroenriched.tplus.adas13.tau.rs,
                                                 "formula_smallmodel" = formula.MCI.adas13.rs,
                                                 "smallmodel" = model.MCI.neuroenriched.tplus.adas13.rs,
                                                 "sample_sizes" = seq(100, 800, by=100),
                                                 "nsim"=500,
                                                 "data" = long.MCI.neuroenriched.tplus.adas13_with_treatment)
  
  
  saveRDS(MCIadas13.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/MCIadas13.rds")
  saveRDS(MCIadas13.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/MCIadas13_tplus.rds")
  saveRDS(MCIadas13.neuro.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/MCIadas13_neuro_tplus.rds")
  
  
  saveRDS(MCIadas13.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/MCIadas13_tau.rds")
  saveRDS(MCIadas13.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/MCIadas13_tplus_tau.rds")
  saveRDS(MCIadas13.neuro.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/MCIadas13_neuro_tplus_tau.rds")
  
  
  
  
  #Hipp#################################### #################################### ####################################
  
  
  MCIhipp.sim.list <- list("formula_largemodel" = formula.MCI.hipp.simulation.rs,
                               "largemodel" = simulation.model.MCI.hipp.rs,
                               "formula_smallmodel" = formula.MCI.hipp.rs,
                               "smallmodel" = model.MCI.hipp.rs,
                               "sample_sizes" = seq(100, 800, by=50),
                               "nsim"=500,
                           "data" = long.MCI.hipp_with_treatment
                           )
  
  
  MCIhipp.tplus.sim.list <- list("formula_largemodel" = formula.MCI.hipp.simulation.rs,
                                     "largemodel" = simulation.model.MCI.tplus.hipp.rs,
                                     "formula_smallmodel" = formula.MCI.hipp.rs,
                                     "smallmodel" = model.MCI.tplus.hipp.rs,
                                     "sample_sizes" = seq(100, 800, by=50),
                                     "nsim"=500,
                                
                                 "data" = long.MCI.tplus.hipp_with_treatment
                                 )
  
  MCIhipp.neuro.tplus.sim.list <- list("formula_largemodel" = formula.MCI.hipp.simulation.rs,
                                           "largemodel" = simulation.model.MCI.neuroenriched.tplus.hipp.rs,
                                           "formula_smallmodel" = formula.MCI.hipp.rs,
                                           "smallmodel" = model.MCI.neuroenriched.tplus.hipp.rs,
                                           "sample_sizes" = seq(100, 800, by=50),
                                           "nsim"=500,
                                       
                                       "data" = long.MCI.neuroenriched.tplus.hipp_with_treatment)
  
  
  
  MCIhipp.sim.tau.list <- list("formula_largemodel" = formula.MCI.hipp.simulation.rs,
                                   "largemodel" = simulation.model.MCI.hipp.tau.rs,
                                   "formula_smallmodel" = formula.MCI.hipp.rs,
                                   "smallmodel" = model.MCI.hipp.rs,
                                   "sample_sizes" = seq(100, 800, by=50),
                                   "nsim"=500,
                               "data" = long.MCI.hipp_with_treatment
                               )
  
  
  MCIhipp.tplus.sim.tau.list <- list("formula_largemodel" = formula.MCI.hipp.simulation.rs,
                                         "largemodel" = simulation.model.MCI.tplus.hipp.tau.rs,
                                         "formula_smallmodel" = formula.MCI.hipp.rs,
                                         "smallmodel" = model.MCI.tplus.hipp.rs,
                                         "sample_sizes" = seq(100, 800, by=50),
                                         "nsim"=500,
                                    
                                     "data" = long.MCI.tplus.hipp_with_treatment
                                     )
  
  MCIhipp.neuro.tplus.sim.tau.list <- list("formula_largemodel" = formula.MCI.hipp.simulation.rs,
                                               "largemodel" = simulation.model.MCI.neuroenriched.tplus.hipp.tau.rs,
                                               "formula_smallmodel" = formula.MCI.hipp.rs,
                                               "smallmodel" = model.MCI.neuroenriched.tplus.hipp.rs,
                                               "sample_sizes" = seq(100, 800, by=50),
                                               "nsim"=500,
                                           
                                           "data" = long.MCI.neuroenriched.tplus.hipp_with_treatment)
  
  
  saveRDS(MCIhipp.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/MCIhipp.rds")
  saveRDS(MCIhipp.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/MCIhipp_tplus.rds")
  saveRDS(MCIhipp.neuro.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/MCIhipp_neuro_tplus.rds")
  
  
  saveRDS(MCIhipp.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/MCIhipp_tau.rds")
  saveRDS(MCIhipp.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/MCIhipp_tplus_tau.rds")
  saveRDS(MCIhipp.neuro.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/hipp_power_data/MCIhipp_neuro_tplus_tau.rds")
  
  
  

  
  
  #MPACC #################################### #################################### ####################################
  
  MCImpacc.sim.list <- list("formula_largemodel" = formula.MCI.mpacc.simulation.rs,
                                "largemodel" = simulation.model.MCI.mpacc.rs,
                                "formula_smallmodel" = formula.MCI.mpacc.rs,
                                "smallmodel" = model.MCI.mpacc.rs,
                                "sample_sizes" = seq(100, 500, by=20),
                                "nsim"=500,
                            
                            "data" = long.MCI.mpacc_with_treatment
                           
  )
  
  
  MCImpacc.tplus.sim.list <- list("formula_largemodel" = formula.MCI.mpacc.simulation.rs,
                                      "largemodel" = simulation.model.MCI.tplus.mpacc.rs,
                                      "formula_smallmodel" = formula.MCI.mpacc.rs,
                                      "smallmodel" = model.MCI.tplus.mpacc.rs,
                                      "sample_sizes" = seq(100, 500, by=20),
                                      "nsim"=500,
                                  
                                 
                                  "data" = long.MCI.tplus.mpacc_with_treatment
                                 
  )
  
  MCImpacc.neuro.tplus.sim.list <- list("formula_largemodel" = formula.MCI.mpacc.simulation.rs,
                                            "largemodel" = simulation.model.MCI.neuroenriched.tplus.mpacc.rs,
                                            "formula_smallmodel" = formula.MCI.mpacc.rs,
                                            "smallmodel" = model.MCI.neuroenriched.tplus.mpacc.rs,
                                            "sample_sizes" = seq(100, 500, by=20),
                                            "nsim"=500, 
                                        
                                       
                                        "data" = long.MCI.neuroenriched.tplus.mpacc_with_treatment
  )
  
  
  
  MCImpacc.sim.tau.list <- list("formula_largemodel" = formula.MCI.mpacc.simulation.rs,
                                    "largemodel" = simulation.model.MCI.mpacc.tau.rs,
                                    "formula_smallmodel" = formula.MCI.mpacc.rs,
                                    "smallmodel" = model.MCI.mpacc.rs,
                                    "sample_sizes" = seq(100, 500, by=20),
                                    "nsim"=500,
                                
                                "data" = long.MCI.mpacc_with_treatment)
                          
  
  
  MCImpacc.tplus.sim.tau.list <- list("formula_largemodel" = formula.MCI.mpacc.simulation.rs,
                                          "largemodel" = simulation.model.MCI.tplus.mpacc.tau.rs,
                                          "formula_smallmodel" = formula.MCI.mpacc.rs,
                                          "smallmodel" = model.MCI.tplus.mpacc.rs,
                                          "sample_sizes" = seq(100, 500, by=20),
                                          "nsim"=500,
                                      
                                    
                                      "data" = long.MCI.tplus.mpacc_with_treatment
                                      
  )
  
  MCImpacc.neuro.tplus.sim.tau.list <- list("formula_largemodel" = formula.MCI.mpacc.simulation.rs,
                                                "largemodel" = simulation.model.MCI.neuroenriched.tplus.mpacc.tau.rs,
                                                "formula_smallmodel" = formula.MCI.mpacc.rs,
                                                "smallmodel" = model.MCI.neuroenriched.tplus.mpacc.rs,
                                                "sample_sizes" = seq(100, 500, by=20),
                                                "nsim"= 500,
                                            
                                           
                                            "data" = long.MCI.neuroenriched.tplus.mpacc_with_treatment
  )
  
  
  saveRDS(MCImpacc.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/MCImpacc.rds")
  saveRDS(MCImpacc.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/MCImpacc_tplus.rds")
  saveRDS(MCImpacc.neuro.tplus.sim.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/MCImpacc_neuro_tplus.rds")
  
  
  saveRDS(MCImpacc.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/MCImpacc_tau.rds")
  saveRDS(MCImpacc.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/MCImpacc_tplus_tau.rds")
  saveRDS(MCImpacc.neuro.tplus.sim.tau.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mpacc_power_data/MCImpacc_neuro_tplus_tau.rds")
  
  #################################### #################################### ####################################
  
}

mci.adas13.longpower.est <- longpower::lmmpower(model.MCI.adas13.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.hipp.longpower.est   <- longpower::lmmpower(model.MCI.hipp.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.mpacc.longpower.est  <- longpower::lmmpower(model.MCI.mpacc.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]

mci.tplus.adas13.longpower.est <- longpower::lmmpower(model.MCI.tplus.adas13.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.tplus.hipp.longpower.est <- longpower::lmmpower(model.MCI.tplus.hipp.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.tplus.mpacc.longpower.est <- longpower::lmmpower(model.MCI.tplus.mpacc.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]

mci.neuroenriched.tplus.adas13.longpower.est <-longpower::lmmpower(model.MCI.neuroenriched.tplus.adas13.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.neuroenriched.tplus.hipp.longpower.est <-longpower::lmmpower(model.MCI.neuroenriched.tplus.hipp.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.neuroenriched.tplus.mpacc.longpower.est <-longpower::lmmpower(model.MCI.neuroenriched.tplus.mpacc.rs, pct.change=.5, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]


# Tau Related Decline


mci.adas13.longpower.tau.est <- longpower::lmmpower(model.MCI.adas13.rs, delta= -0.770000, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.hipp.longpower.tau.est   <- longpower::lmmpower(model.MCI.hipp.rs, delta=43.66, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.mpacc.longpower.tau.est  <- longpower::lmmpower(model.MCI.mpacc.rs, delta=.5555, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]

mci.tplus.adas13.longpower.tau.est <- longpower::lmmpower(model.MCI.tplus.adas13.rs, delta = -1, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.tplus.hipp.longpower.tau.est   <- longpower::lmmpower(model.MCI.tplus.hipp.rs, delta   = 46, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]
mci.tplus.mpacc.longpower.tau.est  <- longpower::lmmpower(model.MCI.tplus.mpacc.rs, delta=0.71000, sig.level=0.05, power=.8, t=c(0, .5, 1, 1.5, 2))["n"]


longpower.estimates <- list( "Overall" = list("ADAS" = mci.adas13.longpower.est,
                                              "HIPP" = mci.hipp.longpower.est,
                                              "MPACC" = mci.mpacc.longpower.est,
                                              "ADAS_tplus" = mci.tplus.adas13.longpower.est,
                                              "HIPP_tplus" = mci.tplus.hipp.longpower.est,
                                              "MPACC_tplus" = mci.tplus.mpacc.longpower.est,
                                              "ADAS_neuro_tplus" = mci.neuroenriched.tplus.adas13.longpower.est,
                                              "HIPP_neuro_tplus" = mci.neuroenriched.tplus.hipp.longpower.est,
                                              "MPACC_neuro_tplus" = mci.neuroenriched.tplus.mpacc.longpower.est),
                             "Tau" = list("ADAS" = mci.adas13.longpower.tau.est,
                                          "HIPP" = mci.hipp.longpower.tau.est,
                                          "MPACC" = mci.mpacc.longpower.tau.est,
                                          "ADAS_tplus" = mci.tplus.adas13.longpower.tau.est,
                                          "HIPP_tplus" = mci.tplus.hipp.longpower.tau.est,
                                          "MPACC_tplus" = mci.tplus.mpacc.longpower.tau.est))


saveRDS(longpower.estimates, "/Users/adamgabriellang/Desktop/clinical_trial_sim/longpower_estimates.rds")
