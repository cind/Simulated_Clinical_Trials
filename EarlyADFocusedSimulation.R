library(nlme)
library(sjPlot)
#read in early ad data
earlyadcohorts       <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts.rds")
earlyadcohorts.tplus <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_tplus.rds")
earlyadcohorts.neuroenriched <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_neuroenriched.rds")

View(earlyadcohorts)
######################## ADAS13 ######################## 
cs.earlyad.adas13     <- earlyadcohorts$cs$ADAS13
long.earlyad.adas13   <-  earlyadcohorts$long$ADAS13

#Tau+
cs.earlyad.tplus.adas13    <- earlyadcohorts.tplus$cs$ADAS13
long.earlyad.tplus.adas13  <-  earlyadcohorts.tplus$long$ADAS13


# Neuroenriched
cs.earlyad.neuroenriched.adas13    <- earlyadcohorts.neuroenriched$cs$ADAS13
long.earlyad.neuroenriched.adas13  <-  earlyadcohorts.neuroenriched$long$ADAS13

######################### HIPPOCAMPUS ######################## 
#LH and RH dataframes are identical and have both left/right hippocampus regions (just use LH)
cs.earlyad.hipp    <- earlyadcohorts$cs$LH
long.earlyad.hipp  <-  earlyadcohorts$long$LH

cs.earlyad.hipp$hipp_average    <- (cs.earlyad.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.hipp$hipp_average  <- (long.earlyad.hipp$ST29SV_harmonized_icv_adj + long.earlyad.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right

#Tau+
cs.earlyad.tplus.hipp    <- earlyadcohorts.tplus$cs$LH
long.earlyad.tplus.hipp  <-  earlyadcohorts.tplus$long$LH

cs.earlyad.tplus.hipp$hipp_average   <- (cs.earlyad.tplus.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.tplus.hipp$hipp_average <- (long.earlyad.tplus.hipp$ST29SV_harmonized_icv_adj + long.earlyad.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right


# Neuroenriched
cs.earlyad.neuroenriched.hipp    <- earlyadcohorts.neuroenriched$cs$LH
long.earlyad.neuroenriched.hipp  <-  earlyadcohorts.neuroenriched$long$LH

cs.earlyad.neuroenriched.hipp$hipp_average   <- (cs.earlyad.neuroenriched.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.neuroenriched.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.neuroenriched.hipp$hipp_average <- (long.earlyad.neuroenriched.hipp$ST29SV_harmonized_icv_adj + long.earlyad.neuroenriched.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right



######################## MPACC######################## 
cs.earlyad.mpacc     <- earlyadcohorts$cs$MPACC
long.earlyad.mpacc   <-  earlyadcohorts$long$MPACC

#Tau+
cs.earlyad.tplus.mpacc    <- earlyadcohorts.tplus$cs$MPACC
long.earlyad.tplus.mpacc  <-  earlyadcohorts.tplus$long$MPACC


# Neuroenriched
cs.earlyad.neuroenriched.mpacc    <- earlyadcohorts.neuroenriched$cs$MPACC
long.earlyad.neuroenriched.mpacc  <-  earlyadcohorts.neuroenriched$long$MPACC





#Descriptive table
######################## ADAS13 ########################
earlyad.adas13.desc <- cs.earlyad.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                            "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                            "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                            "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.adas13.desc$TauPos_full <- factor(earlyad.adas13.desc$TauPos_full)
earlyad.adas13.desc$AmyPos_full <- factor(earlyad.adas13.desc$AmyPos_full)
colnames(earlyad.adas13.desc)   <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                   "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                   "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                   "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.adas13.desc <- table1(earlyad.adas13.desc, splitby = ~ Diagnosis)
earlyad.adas13.desc <- earlyad.adas13.desc$Table1
earlyad.adas13.desc <- earlyad.adas13.desc[c(1, 6:9, 11:nrow(earlyad.adas13.desc)),]


#Enriched for TAU
earlyad.tplus.adas13.desc <- cs.earlyad.tplus.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                        "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                        "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.tplus.adas13.desc$TauPos_full <- factor(earlyad.tplus.adas13.desc$TauPos_full)
earlyad.tplus.adas13.desc$AmyPos_full <- factor(earlyad.tplus.adas13.desc$AmyPos_full)

colnames(earlyad.tplus.adas13.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                         "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                         "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                         "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.tplus.adas13.desc <- table1(earlyad.tplus.adas13.desc, splitby = ~ Diagnosis)
earlyad.tplus.adas13.desc <- earlyad.tplus.adas13.desc$Table1
earlyad.tplus.adas13.desc <- earlyad.tplus.adas13.desc[c(1, 6:9, 11:nrow(earlyad.tplus.adas13.desc)),]


#Neuroenriched
earlyad.neuroenriched.adas13.desc <- cs.earlyad.neuroenriched.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                        "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                        "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.neuroenriched.adas13.desc$TauPos_full <- factor(earlyad.neuroenriched.adas13.desc$TauPos_full)
earlyad.neuroenriched.adas13.desc$AmyPos_full <- factor(earlyad.neuroenriched.adas13.desc$AmyPos_full)

colnames(earlyad.neuroenriched.adas13.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                         "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                         "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                         "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.adas13.desc <- table1(earlyad.neuroenriched.adas13.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.adas13.desc <- earlyad.neuroenriched.adas13.desc$Table1
earlyad.neuroenriched.adas13.desc <- earlyad.neuroenriched.adas13.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.adas13.desc)),]






######################## HIPPOCAMPUS ######################## 
earlyad.hipp.desc <- cs.earlyad.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                        "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl",
                                        "hipp_average" ,"mPACCtrailsB_bl", "TauPos_full", 
                                        "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.hipp.desc$TauPos_full <- factor(earlyad.hipp.desc$TauPos_full)
earlyad.hipp.desc$AmyPos_full <- factor(earlyad.hipp.desc$AmyPos_full)

colnames(earlyad.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                 "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                 "Hippocampus (L.R average_ICV adjusted_Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                 "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")
earlyad.hipp.desc <- table1(earlyad.hipp.desc, splitby = ~ Diagnosis)
earlyad.hipp.desc <- earlyad.hipp.desc$Table1
earlyad.hipp.desc <- earlyad.hipp.desc[c(1, 6:9, 11:nrow(earlyad.hipp.desc)),]

#Enriched for TAU
earlyad.tplus.hipp.desc <- cs.earlyad.tplus.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                    "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl",
                                                    "hipp_average" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                    "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.tplus.hipp.desc$TauPos_full <- factor(earlyad.tplus.hipp.desc$TauPos_full)
earlyad.tplus.hipp.desc$AmyPos_full <- factor(earlyad.tplus.hipp.desc$AmyPos_full)

colnames(earlyad.tplus.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                       "Hippocampus (L.R average_ICV adjusted_Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.tplus.hipp.desc <- table1(earlyad.tplus.hipp.desc, splitby = ~ Diagnosis)
earlyad.tplus.hipp.desc <- earlyad.tplus.hipp.desc$Table1
earlyad.tplus.hipp.desc <- earlyad.tplus.hipp.desc[c(1, 6:9, 11:nrow(earlyad.tplus.hipp.desc)),]


# Neuroenriched
earlyad.neuroenriched.hipp.desc <- cs.earlyad.neuroenriched.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                    "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl",
                                                    "hipp_average" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                    "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.neuroenriched.hipp.desc$TauPos_full <- factor(earlyad.neuroenriched.hipp.desc$TauPos_full)
earlyad.neuroenriched.hipp.desc$AmyPos_full <- factor(earlyad.neuroenriched.hipp.desc$AmyPos_full)

colnames(earlyad.neuroenriched.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                       "Hippocampus (L.R average_ICV adjusted_Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.hipp.desc <- table1(earlyad.neuroenriched.hipp.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.hipp.desc <- earlyad.neuroenriched.hipp.desc$Table1
earlyad.neuroenriched.hipp.desc <- earlyad.neuroenriched.hipp.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.hipp.desc)),]






######################## MPACC ######################## 


earlyad.mpacc.desc <- cs.earlyad.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                            "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                            "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                            "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.mpacc.desc$TauPos_full <- factor(earlyad.mpacc.desc$TauPos_full)
earlyad.mpacc.desc$AmyPos_full <- factor(earlyad.mpacc.desc$AmyPos_full)
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
                                                        "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                        "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.tplus.mpacc.desc$TauPos_full <- factor(earlyad.tplus.mpacc.desc$TauPos_full)
earlyad.tplus.mpacc.desc$AmyPos_full <- factor(earlyad.tplus.mpacc.desc$AmyPos_full)

colnames(earlyad.tplus.mpacc.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                         "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                         "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                         "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.tplus.mpacc.desc <- table1(earlyad.tplus.mpacc.desc, splitby = ~ Diagnosis)
earlyad.tplus.mpacc.desc <- earlyad.tplus.mpacc.desc$Table1
earlyad.tplus.mpacc.desc <- earlyad.tplus.mpacc.desc[c(1, 6:9, 11:nrow(earlyad.tplus.mpacc.desc)),]


#Neuroenriched
earlyad.neuroenriched.mpacc.desc <- cs.earlyad.neuroenriched.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                        "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                        "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                                        "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.neuroenriched.mpacc.desc$TauPos_full <- factor(earlyad.neuroenriched.mpacc.desc$TauPos_full)
earlyad.neuroenriched.mpacc.desc$AmyPos_full <- factor(earlyad.neuroenriched.mpacc.desc$AmyPos_full)

colnames(earlyad.neuroenriched.mpacc.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                 "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                 "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                 "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.mpacc.desc <- table1(earlyad.neuroenriched.mpacc.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.mpacc.desc <- earlyad.neuroenriched.mpacc.desc$Table1
earlyad.neuroenriched.mpacc.desc <- earlyad.neuroenriched.mpacc.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.mpacc.desc)),]


View(earlyad.neuroenriched.mpacc.desc)


#fit repeated measures model on longitudinal data
# want observed scores, just use time and dont control for other variables


long.earlyad.adas13 <- MMRMTime(long.earlyad.adas13)
long.earlyad.tplus.adas13 <- MMRMTime(long.earlyad.tplus.adas13)
long.earlyad.neuroenriched.adas13 <- MMRMTime(long.earlyad.neuroenriched.adas13)


long.earlyad.hipp <-  MMRMTime(long.earlyad.hipp)
long.earlyad.tplus.hipp <- MMRMTime(long.earlyad.tplus.hipp)
long.earlyad.neuroenriched.hipp <- MMRMTime(long.earlyad.neuroenriched.hipp)


long.earlyad.mpacc <- MMRMTime(long.earlyad.mpacc)
long.earlyad.tplus.mpacc <- MMRMTime(long.earlyad.tplus.mpacc)
long.earlyad.neuroenriched.mpacc <- MMRMTime(long.earlyad.neuroenriched.mpacc)



######################## ADAS13 ######################## 
mmrm.earlyad.adas13        <- gls(TOTAL13~factor(new_time_mmrm),
                                    na.action=na.omit, data=long.earlyad.adas13,
                                    correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                    weights=nlme::varIdent(form=~1|new_time_mmrm))


#Enriched for TAU
mmrm.earlyad.tplus.adas13<- gls(TOTAL13~factor(new_time_mmrm),
                                na.action=na.omit, data=long.earlyad.tplus.adas13,
                                correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                weights=nlme::varIdent(form=~1|new_time_mmrm))
#Neuroenriched
mmrm.earlyad.neuroenriched.adas13<- gls(TOTAL13~factor(new_time_mmrm),
                                na.action=na.omit, data=long.earlyad.neuroenriched.adas13,
                                correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                weights=nlme::varIdent(form=~1|new_time_mmrm))


######################## HIPPOCAMPUS ######################## 
mmrm.earlyad.hipp<- gls(hipp_average~factor(new_time_mmrm),
                                na.action=na.omit, data=long.earlyad.hipp,
                                correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                weights=nlme::varIdent(form=~1|new_time_mmrm))
#Enriched for TAU
mmrm.earlyad.tplus.hipp <- gls(hipp_average~factor(new_time_mmrm),
                        na.action=na.omit, data=long.earlyad.tplus.hipp,
                        correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                        weights=nlme::varIdent(form=~1|new_time_mmrm))


#Neuroenriched
mmrm.earlyad.neuroenriched.hipp <- gls(hipp_average~factor(new_time_mmrm),
                               na.action=na.omit, data=long.earlyad.neuroenriched.hipp,
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
#Neuroenriched
mmrm.earlyad.neuroenriched.mpacc<- gls(mPACCtrailsB~factor(new_time_mmrm),
                                        na.action=na.omit, data=long.earlyad.neuroenriched.mpacc,
                                        correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                        weights=nlme::varIdent(form=~1|new_time_mmrm))








#get coefficients and standard errors

######################## ADAS13 ######################## 
coefs.earlyad.adas13     <- unname(coef(summary(mmrm.earlyad.adas13))[, "Value"])
stderror.earlyad.adas13  <- unname(coef(summary(mmrm.earlyad.adas13))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.adas13     <- unname(coef(summary(mmrm.earlyad.tplus.adas13))[, "Value"])
stderror.earlyad.tplus.adas13  <- unname(coef(summary(mmrm.earlyad.tplus.adas13))[, "Std.Error"])

#Neuroenriched
coefs.earlyad.neuroenriched.adas13     <- unname(coef(summary(mmrm.earlyad.neuroenriched.adas13))[, "Value"])
stderror.earlyad.neuroenriched.adas13  <- unname(coef(summary(mmrm.earlyad.neuroenriched.adas13))[, "Std.Error"])

######################## HIPPOCAMPUS ######################## 
coefs.earlyad.hipp     <- unname(coef(summary(mmrm.earlyad.hipp))[, "Value"])
stderror.earlyad.hipp  <- unname(coef(summary(mmrm.earlyad.hipp))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.hipp     <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Value"])
stderror.earlyad.tplus.hipp  <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Std.Error"])

coefs.earlyad.neuroenriched.hipp     <- unname(coef(summary(mmrm.earlyad.neuroenriched.hipp))[, "Value"])
stderror.earlyad.neuroenriched.hipp  <- unname(coef(summary(mmrm.earlyad.neuroenriched.hipp))[, "Std.Error"])


######################## MPACC ######################## 
coefs.earlyad.mpacc     <- unname(coef(summary(mmrm.earlyad.mpacc))[, "Value"])
stderror.earlyad.mpacc  <- unname(coef(summary(mmrm.earlyad.mpacc))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.mpacc     <- unname(coef(summary(mmrm.earlyad.tplus.mpacc))[, "Value"])
stderror.earlyad.tplus.mpacc  <- unname(coef(summary(mmrm.earlyad.tplus.mpacc))[, "Std.Error"])

#Neuroenriched
coefs.earlyad.neuroenriched.mpacc     <- unname(coef(summary(mmrm.earlyad.neuroenriched.mpacc))[, "Value"])
stderror.earlyad.neuroenriched.mpacc  <- unname(coef(summary(mmrm.earlyad.neuroenriched.mpacc))[, "Std.Error"])








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


earlyad.neuroenriched.adas13.plotting.data <- data.frame("coefs" = coefs.earlyad.neuroenriched.adas13,
                                                 "stderr" = stderror.earlyad.neuroenriched.adas13)
earlyad.neuroenriched.adas13.plotting.data$ci_low     <-  -1.96* earlyad.neuroenriched.adas13.plotting.data$stderr
earlyad.neuroenriched.adas13.plotting.data$ci_hi      <-   1.96* earlyad.neuroenriched.adas13.plotting.data$stderr
earlyad.neuroenriched.adas13.plotting.data$means      <-  earlyad.neuroenriched.adas13.plotting.data$coefs
earlyad.neuroenriched.adas13.plotting.data$means[2:5] <-  earlyad.neuroenriched.adas13.plotting.data$means[2:5] + earlyad.neuroenriched.adas13.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.neuroenriched.adas13.plotting.data$ci_low     <-  earlyad.neuroenriched.adas13.plotting.data$ci_low+ earlyad.neuroenriched.adas13.plotting.data$means
earlyad.neuroenriched.adas13.plotting.data$ci_hi      <-  earlyad.neuroenriched.adas13.plotting.data$ci_hi+ earlyad.neuroenriched.adas13.plotting.data$means
earlyad.neuroenriched.adas13.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.neuroenriched.adas13.plotting.data$Enrichment <- "No Copathologies"





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

earlyad.neuroenriched.hipp.plotting.data <- data.frame("coefs" = coefs.earlyad.neuroenriched.hipp,
                                                         "stderr" = stderror.earlyad.neuroenriched.hipp)
earlyad.neuroenriched.hipp.plotting.data$ci_low     <-  -1.96* earlyad.neuroenriched.hipp.plotting.data$stderr
earlyad.neuroenriched.hipp.plotting.data$ci_hi      <-   1.96* earlyad.neuroenriched.hipp.plotting.data$stderr
earlyad.neuroenriched.hipp.plotting.data$means      <-  earlyad.neuroenriched.hipp.plotting.data$coefs
earlyad.neuroenriched.hipp.plotting.data$means[2:5] <-  earlyad.neuroenriched.hipp.plotting.data$means[2:5] + earlyad.neuroenriched.hipp.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.neuroenriched.hipp.plotting.data$ci_low     <-  earlyad.neuroenriched.hipp.plotting.data$ci_low+ earlyad.neuroenriched.hipp.plotting.data$means
earlyad.neuroenriched.hipp.plotting.data$ci_hi      <-  earlyad.neuroenriched.hipp.plotting.data$ci_hi+ earlyad.neuroenriched.hipp.plotting.data$means
earlyad.neuroenriched.hipp.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.neuroenriched.hipp.plotting.data$Enrichment <- "No Copathologies"


#create dataset for plotting
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


earlyad.neuroenriched.mpacc.plotting.data <- data.frame("coefs" = coefs.earlyad.neuroenriched.mpacc,
                                                         "stderr" = stderror.earlyad.neuroenriched.mpacc)
earlyad.neuroenriched.mpacc.plotting.data$ci_low     <-  -1.96* earlyad.neuroenriched.mpacc.plotting.data$stderr
earlyad.neuroenriched.mpacc.plotting.data$ci_hi      <-   1.96* earlyad.neuroenriched.mpacc.plotting.data$stderr
earlyad.neuroenriched.mpacc.plotting.data$means      <-  earlyad.neuroenriched.mpacc.plotting.data$coefs
earlyad.neuroenriched.mpacc.plotting.data$means[2:5] <-  earlyad.neuroenriched.mpacc.plotting.data$means[2:5] + earlyad.neuroenriched.mpacc.plotting.data$coefs[1] #add intercept to each estimate for plotting
earlyad.neuroenriched.mpacc.plotting.data$ci_low     <-  earlyad.neuroenriched.mpacc.plotting.data$ci_low+ earlyad.neuroenriched.mpacc.plotting.data$means
earlyad.neuroenriched.mpacc.plotting.data$ci_hi      <-  earlyad.neuroenriched.mpacc.plotting.data$ci_hi+ earlyad.neuroenriched.mpacc.plotting.data$means
earlyad.neuroenriched.mpacc.plotting.data$time       <-  c(0, .5, 1, 1.5, 2)
earlyad.neuroenriched.mpacc.plotting.data$Enrichment <- "No Copathologies"





## Plots
# Combine non-enriched and tau enriched data for plotting

######################## ADAS13 ######################## 
earlyad.adas13.plotting.data <- rbind(earlyad.adas13.plotting.data,
                                      earlyad.tplus.adas13.plotting.data,
                                      earlyad.neuroenriched.adas13.plotting.data)
earlyad.adas13.mmrm.plot <- ggplot(earlyad.adas13.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.adas13.mmrm.plot <- earlyad.adas13.mmrm.plot + xlab("Time (Years)") + ylab("Mean ADAS-13 (95% CI)") + ylim(20, 35)
earlyad.adas13.mmrm.plot


######################## HIPPOCAMPUS ######################## 
earlyad.hipp.plotting.data <- rbind(earlyad.hipp.plotting.data,
                                    earlyad.tplus.hipp.plotting.data,
                                    earlyad.neuroenriched.hipp.plotting.data)
earlyad.hipp.mmrm.plot <- ggplot(earlyad.hipp.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.hipp.mmrm.plot <- earlyad.hipp.mmrm.plot + xlab("Time (Years)") + ylab("Mean Hippocampus (95% CI)") + ylim(2700, 3800)
earlyad.hipp.mmrm.plot



######################## MPACC ######################## 
earlyad.mpacc.plotting.data <- rbind(earlyad.mpacc.plotting.data,
                                      earlyad.tplus.mpacc.plotting.data,
                                      earlyad.neuroenriched.mpacc.plotting.data)
earlyad.mpacc.mmrm.plot <- ggplot(earlyad.mpacc.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.mpacc.mmrm.plot <- earlyad.mpacc.mmrm.plot + xlab("Time (Years)") + ylab("Mean mPACCtrailsB (95% CI)") 
earlyad.mpacc.mmrm.plot


######################################## Simulations #################################################


######################## ADAS13 ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.earlyad.adas13 <- StratifyContVar(cs.earlyad.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.adas13 <- StratifyContVar(cs.earlyad.tplus.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.adas13 <- StratifyContVar(cs.earlyad.neuroenriched.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))





#Assign treatment and placebo by blocks
long.earlyad.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.adas13, long.earlyad.adas13)
long.earlyad.tplus.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.adas13, long.earlyad.tplus.adas13)
long.earlyad.neuroenriched.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.adas13, long.earlyad.neuroenriched.adas13)

######################## HIPPOCAMPUS ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.earlyad.hipp <- StratifyContVar(cs.earlyad.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.hipp <- StratifyContVar(cs.earlyad.tplus.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.hipp <- StratifyContVar(cs.earlyad.neuroenriched.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))





#Assign treatment and placebo by blocks
long.earlyad.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.hipp, long.earlyad.hipp)
long.earlyad.tplus.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.hipp, long.earlyad.tplus.hipp)
long.earlyad.neuroenriched.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.hipp, long.earlyad.neuroenriched.hipp)


######################## MPACC ######################## 

# Categorize continuous variables for block randomization
cs.earlyad.mpacc <- StratifyContVar(cs.earlyad.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.mpacc <- StratifyContVar(cs.earlyad.tplus.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.mpacc <- StratifyContVar(cs.earlyad.neuroenriched.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))





#Assign treatment and placebo by blocks
long.earlyad.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.mpacc, long.earlyad.mpacc)
long.earlyad.tplus.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.mpacc, long.earlyad.tplus.mpacc)
long.earlyad.neuroenriched.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.mpacc, long.earlyad.neuroenriched.mpacc)




#long.earlyad.adas13_with_treatment_example <-RandomizeTreatment2(cs.earlyad.adas13, long.earlyad.adas13)
#treatment.ex <- long.earlyad.adas13_with_treatment_example[[1]]
#placebo.ex <- long.earlyad.adas13_with_treatment_example[[2]]
#treatment.ex <- treatment.ex[,c("CDGLOBAL_bl", "PTEDUCAT_bl_strat", "AGE_bl_strat", "MMSE_bl_strat", "fulllewy", "fulltdp43", "fullcaa", "TauPos_full_bl")]
#placebo.ex <- placebo.ex[,c("CDGLOBAL_bl", "PTEDUCAT_bl_strat", "AGE_bl_strat", "MMSE_bl_strat", "fulllewy", "fulltdp43", "fullcaa", "TauPos_full_bl")]
#colnames(treatment.ex) <- colnames(placebo.ex) <- c("CDR", "Education (Cat.)", "Age (Cat.)", "MMSE (Cat.)", "Lewy", "TDP43", "CAA", "Tau")
#treatment.ex$CDR <- factor(treatment.ex$CDR)
#placebo.ex$CDR <- factor(placebo.ex$CDR)
#treat.desc <- table1(treatment.ex)$Table1
#placebo.desc <- table1(placebo.ex)$Table1


# Calculate outcome for overall decline
# Not controlling for neuropathologies or AD to get overall decline, controlling for demographic variables


######################## ADAS13 ######################## 
formula.earlyad.adas13 <- "TOTAL13 ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
model.earlyad.adas13 <- MapLmer(newdata = long.earlyad.adas13_with_treatment,
                                formula.model = formula.earlyad.adas13)
model.earlyad.tplus.adas13 <- MapLmer(newdata = long.earlyad.tplus.adas13_with_treatment,
                                formula.model = formula.earlyad.adas13)
model.earlyad.neuroenriched.adas13 <- MapLmer(newdata = long.earlyad.neuroenriched.adas13_with_treatment,
                                      formula.model = formula.earlyad.adas13)




formula.earlyad.adas13.simulation             <- "TOTAL13 ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
simulation.model.earlyad.adas13               <- BuildSimulationModelNoPath(model.earlyad.adas13, formula.earlyad.adas13.simulation, long.earlyad.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.adas13         <- BuildSimulationModelNoPath(model.earlyad.tplus.adas13, formula.earlyad.adas13.simulation, long.earlyad.tplus.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.adas13 <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.adas13, formula.earlyad.adas13.simulation, long.earlyad.neuroenriched.adas13_with_treatment, "not-controlled")


dpm.data.earlyad.adas13 <- get_model_data(simulation.model.earlyad.adas13, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.adas13 <- get_model_data(simulation.model.earlyad.tplus.adas13, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.adas13 <- get_model_data(simulation.model.earlyad.neuroenriched.adas13, terms = "new_time", type = "int") 

dpm.data.earlyad.adas13$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.adas13$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.adas13$Enrichment <- "No Copathologies"

fulldpm.earlyad.adas13 <- rbind(dpm.data.earlyad.adas13,dpm.data.earlyad.tplus.adas13,dpm.data.earlyad.neuroenriched.adas13)
fulldpm.earlyad.adas13$Treatment <- fulldpm.earlyad.adas13$group



######################## HIPPOCAMPUS ######################## 

formula.earlyad.hipp <- "hipp_average ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
model.earlyad.hipp <- MapLmer(newdata = long.earlyad.hipp_with_treatment,
                                formula.model = formula.earlyad.hipp)
model.earlyad.tplus.hipp <- MapLmer(newdata = long.earlyad.tplus.hipp_with_treatment,
                                      formula.model = formula.earlyad.hipp)
model.earlyad.neuroenriched.hipp <- MapLmer(newdata = long.earlyad.neuroenriched.hipp_with_treatment,
                                              formula.model = formula.earlyad.hipp)




formula.earlyad.hipp.simulation             <- "hipp_average ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
simulation.model.earlyad.hipp               <- BuildSimulationModelNoPath(model.earlyad.hipp, formula.earlyad.hipp.simulation, long.earlyad.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.hipp         <- BuildSimulationModelNoPath(model.earlyad.tplus.hipp, formula.earlyad.hipp.simulation, long.earlyad.tplus.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.hipp <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.hipp, formula.earlyad.hipp.simulation, long.earlyad.neuroenriched.hipp_with_treatment, "not-controlled")


dpm.data.earlyad.hipp <- get_model_data(simulation.model.earlyad.hipp, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.hipp <- get_model_data(simulation.model.earlyad.tplus.hipp, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.hipp <- get_model_data(simulation.model.earlyad.neuroenriched.hipp, terms = "new_time", type = "int") 

dpm.data.earlyad.hipp$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.hipp$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.hipp$Enrichment <- "No Copathologies"

fulldpm.earlyad.hipp <- rbind(dpm.data.earlyad.hipp,dpm.data.earlyad.tplus.hipp,dpm.data.earlyad.neuroenriched.hipp)
fulldpm.earlyad.hipp$Treatment <- fulldpm.earlyad.hipp$group


######################## MPACC ######################## 
formula.earlyad.mpacc <- "mPACCtrailsB ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
model.earlyad.mpacc <- MapLmer(newdata = long.earlyad.mpacc_with_treatment,
                                formula.model = formula.earlyad.mpacc)
model.earlyad.tplus.mpacc <- MapLmer(newdata = long.earlyad.tplus.mpacc_with_treatment,
                                      formula.model = formula.earlyad.mpacc)
model.earlyad.neuroenriched.mpacc <- MapLmer(newdata = long.earlyad.neuroenriched.mpacc_with_treatment,
                                              formula.model = formula.earlyad.mpacc)




formula.earlyad.mpacc.simulation             <- "mPACCtrailsB ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
simulation.model.earlyad.mpacc               <- BuildSimulationModelNoPath(model.earlyad.mpacc, formula.earlyad.mpacc.simulation, long.earlyad.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.mpacc         <- BuildSimulationModelNoPath(model.earlyad.tplus.mpacc, formula.earlyad.mpacc.simulation, long.earlyad.tplus.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.mpacc <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.mpacc, formula.earlyad.mpacc.simulation, long.earlyad.neuroenriched.mpacc_with_treatment, "not-controlled")


dpm.data.earlyad.mpacc <- get_model_data(simulation.model.earlyad.mpacc, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.mpacc <- get_model_data(simulation.model.earlyad.tplus.mpacc, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.mpacc <- get_model_data(simulation.model.earlyad.neuroenriched.mpacc, terms = "new_time", type = "int") 

dpm.data.earlyad.mpacc$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.mpacc$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.mpacc$Enrichment <- "No Copathologies"

fulldpm.earlyad.mpacc <- rbind(dpm.data.earlyad.mpacc,dpm.data.earlyad.tplus.mpacc,dpm.data.earlyad.neuroenriched.mpacc)
fulldpm.earlyad.mpacc$Treatment <- fulldpm.earlyad.mpacc$group


earlyad.adas13.modeling.list <- list("data" = list("unenriched"=long.earlyad.adas13_with_treatment,
                                                  "taupos"=long.earlyad.tplus.adas13_with_treatment,
                                                  "nocopath"=long.earlyad.neuroenriched.adas13_with_treatment),
                                    "models" = list("unenriched"=model.earlyad.adas13,
                                                    "taupos"=model.earlyad.tplus.adas13,
                                                    "nocopath"=model.earlyad.neuroenriched.adas13),
                                    "formula" = replicate(3,formula.earlyad.adas13.simulation, simplify = FALSE),
                                    "fcompare" = replicate(3,formula.earlyad.adas13, simplify = FALSE),
                                    "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                    "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE))

earlyad.hipp.modeling.list <- list("data" = list("unenriched"=long.earlyad.hipp_with_treatment,
                                                   "taupos"=long.earlyad.tplus.hipp_with_treatment,
                                                   "nocopath"=long.earlyad.neuroenriched.hipp_with_treatment),
                                     "models" = list("unenriched"=model.earlyad.hipp,
                                                     "taupos"=model.earlyad.tplus.hipp,
                                                     "nocopath"=model.earlyad.neuroenriched.hipp),
                                     "formula" = replicate(3,formula.earlyad.hipp.simulation, simplify = FALSE),
                                     "fcompare" = replicate(3,formula.earlyad.hipp, simplify = FALSE),
                                     "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                     "yaxislab_dpm" = replicate(3, "Hippocampus", simplify = FALSE))


earlyad.mpacc.modeling.list <- list("data" = list("unenriched"=long.earlyad.mpacc_with_treatment,
                                                   "taupos"=long.earlyad.tplus.mpacc_with_treatment,
                                                   "nocopath"=long.earlyad.neuroenriched.mpacc_with_treatment),
                                     "models" = list("unenriched"=model.earlyad.mpacc,
                                                     "taupos"=model.earlyad.tplus.mpacc,
                                                     "nocopath"=model.earlyad.neuroenriched.mpacc),
                                     "formula" = replicate(3,formula.earlyad.mpacc.simulation, simplify = FALSE),
                                     "fcompare" = replicate(3,formula.earlyad.mpacc, simplify = FALSE),
                                     "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                     "yaxislab_dpm" = replicate(3, "mPACCtrailsB", simplify = FALSE))






