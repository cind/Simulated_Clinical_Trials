#read in early ad data
earlyadcohorts                     <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyadcohorts.rds")
earlyadcohorts.tplus               <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyadcohorts_tplus.rds")
earlyadcohorts.neuroenriched       <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyadcohorts_neuroenriched.rds")
earlyadcohorts.neuroenriched.tplus <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyadcohorts_neuroenriched_tplus.rds")



#keeping subjects with 3 or more time points
#update descriptive statistics and a few plots in powerpoint as this will remove a few subjects
earlyadcohorts                     <- Keep1YearorMore(earlyadcohorts)
earlyadcohorts.tplus               <- Keep1YearorMore(earlyadcohorts.tplus)
earlyadcohorts.neuroenriched       <- Keep1YearorMore(earlyadcohorts.neuroenriched)
earlyadcohorts.neuroenriched.tplus <- Keep1YearorMore(earlyadcohorts.neuroenriched.tplus)


######################## ADAS13 ######################## 
cs.earlyad.adas13     <- earlyadcohorts$cs$ADAS13
long.earlyad.adas13   <-  earlyadcohorts$long$ADAS13


#Tau+
cs.earlyad.tplus.adas13    <- earlyadcohorts.tplus$cs$ADAS13
long.earlyad.tplus.adas13  <-  earlyadcohorts.tplus$long$ADAS13

# Neuroenriched
cs.earlyad.neuroenriched.adas13    <- earlyadcohorts.neuroenriched$cs$ADAS13
long.earlyad.neuroenriched.adas13  <-  earlyadcohorts.neuroenriched$long$ADAS13


# Neuroenriched Tau+
cs.earlyad.neuroenriched.tplus.adas13    <- earlyadcohorts.neuroenriched.tplus$cs$ADAS13
long.earlyad.neuroenriched.tplus.adas13  <-  earlyadcohorts.neuroenriched.tplus$long$ADAS13


######################### HIPPOCAMPUS ######################## 
#LH and RH dataframes are identical and have both left/right hippocampus regions (just use LH)
cs.earlyad.hipp    <- earlyadcohorts$cs$LH
long.earlyad.hipp  <-  earlyadcohorts$long$LH

cs.earlyad.hipp$hipp_average    <- (cs.earlyad.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.hipp$hipp_average  <- (long.earlyad.hipp$ST29SV_harmonized_icv_adj + long.earlyad.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right

#Tau+
cs.earlyad.tplus.hipp           <- earlyadcohorts.tplus$cs$LH
long.earlyad.tplus.hipp         <-  earlyadcohorts.tplus$long$LH


cs.earlyad.tplus.hipp$hipp_average   <- (cs.earlyad.tplus.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.tplus.hipp$hipp_average <- (long.earlyad.tplus.hipp$ST29SV_harmonized_icv_adj + long.earlyad.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right


# Neuroenriched
cs.earlyad.neuroenriched.hipp    <- earlyadcohorts.neuroenriched$cs$LH
long.earlyad.neuroenriched.hipp  <-  earlyadcohorts.neuroenriched$long$LH


cs.earlyad.neuroenriched.hipp$hipp_average   <- (cs.earlyad.neuroenriched.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.neuroenriched.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.neuroenriched.hipp$hipp_average <- (long.earlyad.neuroenriched.hipp$ST29SV_harmonized_icv_adj + long.earlyad.neuroenriched.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right


# Neuroenriched Tau+
cs.earlyad.neuroenriched.tplus.hipp    <- earlyadcohorts.neuroenriched.tplus$cs$LH
long.earlyad.neuroenriched.tplus.hipp <-  earlyadcohorts.neuroenriched.tplus$long$LH



cs.earlyad.neuroenriched.tplus.hipp$hipp_average   <- (cs.earlyad.neuroenriched.tplus.hipp$ST29SV_harmonized_icv_adj + cs.earlyad.neuroenriched.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right
long.earlyad.neuroenriched.tplus.hipp$hipp_average <- (long.earlyad.neuroenriched.tplus.hipp$ST29SV_harmonized_icv_adj + long.earlyad.neuroenriched.tplus.hipp$ST88SV_harmonized_icv_adj) / 2 #Average of left and right


######################## MPACC######################## 
cs.earlyad.mpacc     <- earlyadcohorts$cs$MPACC
long.earlyad.mpacc   <-  earlyadcohorts$long$MPACC



#Tau+
cs.earlyad.tplus.mpacc    <- earlyadcohorts.tplus$cs$MPACC
long.earlyad.tplus.mpacc  <-  earlyadcohorts.tplus$long$MPACC






# Neuroenriched
cs.earlyad.neuroenriched.mpacc    <- earlyadcohorts.neuroenriched$cs$MPACC
long.earlyad.neuroenriched.mpacc  <-  earlyadcohorts.neuroenriched$long$MPACC



# Neuroenriched Tau+
cs.earlyad.neuroenriched.tplus.mpacc    <- earlyadcohorts.neuroenriched.tplus$cs$MPACC
long.earlyad.neuroenriched.tplus.mpacc  <-  earlyadcohorts.neuroenriched.tplus$long$MPACC





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

#Neuroenriched Tau+
earlyad.neuroenriched.tplus.adas13.desc <- cs.earlyad.neuroenriched.tplus.adas13[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                    "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                    "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                                                    "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.neuroenriched.tplus.adas13.desc$TauPos_full <- factor(earlyad.neuroenriched.tplus.adas13.desc$TauPos_full)
earlyad.neuroenriched.tplus.adas13.desc$AmyPos_full <- factor(earlyad.neuroenriched.tplus.adas13.desc$AmyPos_full)

colnames(earlyad.neuroenriched.tplus.adas13.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                       "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                       "ADAS13 (Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                       "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.tplus.adas13.desc <- table1(earlyad.neuroenriched.tplus.adas13.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.tplus.adas13.desc <- earlyad.neuroenriched.tplus.adas13.desc$Table1
earlyad.neuroenriched.tplus.adas13.desc <- earlyad.neuroenriched.tplus.adas13.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.tplus.adas13.desc)),]




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



# Neuroenriched Tau+
earlyad.neuroenriched.tplus.hipp.desc <- cs.earlyad.neuroenriched.tplus.hipp[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                "PTEDUCAT_bl", "MMSE_bl" , "CDRSB_bl",
                                                                                "hipp_average" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                                                "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.neuroenriched.tplus.hipp.desc$TauPos_full <- factor(earlyad.neuroenriched.tplus.hipp.desc$TauPos_full)
earlyad.neuroenriched.tplus.hipp.desc$AmyPos_full <- factor(earlyad.neuroenriched.tplus.hipp.desc$AmyPos_full)

colnames(earlyad.neuroenriched.tplus.hipp.desc) <- c("Diagnosis", "Age (Baseline)", "Gender", 
                                                     "Education (Baseline)","MMSE (Baseline)", "CDRSB (Baseline)", 
                                                     "Hippocampus (L.R average_ICV adjusted_Baseline)" ,"mPACCtrailsB (Baseline)", "Tau+ (Baseline)", 
                                                     "Amy+ (Baseline)", "CAA", "TDP43", "Lewy")

earlyad.neuroenriched.tplus.hipp.desc <- table1(earlyad.neuroenriched.tplus.hipp.desc, splitby = ~ Diagnosis)
earlyad.neuroenriched.tplus.hipp.desc <- earlyad.neuroenriched.tplus.hipp.desc$Table1
earlyad.neuroenriched.tplus.hipp.desc <- earlyad.neuroenriched.tplus.hipp.desc[c(1, 6:9, 11:nrow(earlyad.neuroenriched.tplus.hipp.desc)),]



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

#Neuroenriched Tau+
earlyad.neuroenriched.tplus.mpacc.desc <- cs.earlyad.neuroenriched.tplus.mpacc[,c("DX_bl", "AGE_bl", "PTGENDER", 
                                                                                  "PTEDUCAT_bl", "MMSE_bl" ,"CDRSB_bl",
                                                                                  "TOTAL13_bl" ,"mPACCtrailsB_bl", "TauPos_full", 
                                                                                  "AmyPos_full","fullcaa", "fulltdp43", "fulllewy")]

earlyad.neuroenriched.tplus.mpacc.desc$TauPos_full <- factor(earlyad.neuroenriched.tplus.mpacc.desc$TauPos_full)
earlyad.neuroenriched.tplus.mpacc.desc$AmyPos_full <- factor(earlyad.neuroenriched.tplus.mpacc.desc$AmyPos_full)

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
long.earlyad.neuroenriched.adas13 <- MMRMTime(long.earlyad.neuroenriched.adas13)
long.earlyad.neuroenriched.tplus.adas13 <- MMRMTime(long.earlyad.neuroenriched.tplus.adas13)


long.earlyad.hipp <-  MMRMTime(long.earlyad.hipp)
long.earlyad.tplus.hipp <- MMRMTime(long.earlyad.tplus.hipp)
long.earlyad.neuroenriched.hipp <- MMRMTime(long.earlyad.neuroenriched.hipp)
long.earlyad.neuroenriched.tplus.hipp <- MMRMTime(long.earlyad.neuroenriched.tplus.hipp)


long.earlyad.mpacc <- MMRMTime(long.earlyad.mpacc)
long.earlyad.tplus.mpacc <- MMRMTime(long.earlyad.tplus.mpacc)
long.earlyad.neuroenriched.mpacc <- MMRMTime(long.earlyad.neuroenriched.mpacc)
long.earlyad.neuroenriched.tplus.mpacc <- MMRMTime(long.earlyad.neuroenriched.tplus.mpacc)



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

#Neuroenriched Tau+
mmrm.earlyad.neuroenriched.tplus.adas13<- gls(TOTAL13~factor(new_time_mmrm),
                                              na.action=na.omit, data=long.earlyad.neuroenriched.tplus.adas13,
                                              correlation=nlme::corSymm(form=~new_time_mmrm|RID),
                                              weights=nlme::varIdent(form=~1|new_time_mmrm))



######################## HIPPOCAMPUS ######################## 
mmrm.earlyad.hipp   <- gls(hipp_average~factor(new_time_mmrm),
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
#Neuroenriched
mmrm.earlyad.neuroenriched.mpacc<- gls(mPACCtrailsB~factor(new_time_mmrm),
                                       na.action=na.omit, data=long.earlyad.neuroenriched.mpacc,
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

#Neuroenriched
coefs.earlyad.neuroenriched.adas13     <- unname(coef(summary(mmrm.earlyad.neuroenriched.adas13))[, "Value"])
stderror.earlyad.neuroenriched.adas13  <- unname(coef(summary(mmrm.earlyad.neuroenriched.adas13))[, "Std.Error"])

#Neuroenriched Tau+
coefs.earlyad.neuroenriched.tplus.adas13     <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.adas13))[, "Value"])
stderror.earlyad.neuroenriched.tplus.adas13  <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.adas13))[, "Std.Error"])


######################## HIPPOCAMPUS ######################## 
coefs.earlyad.hipp     <- unname(coef(summary(mmrm.earlyad.hipp))[, "Value"])
stderror.earlyad.hipp  <- unname(coef(summary(mmrm.earlyad.hipp))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.hipp     <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Value"])
stderror.earlyad.tplus.hipp  <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Std.Error"])

#Neuroenriched
coefs.earlyad.neuroenriched.hipp     <- unname(coef(summary(mmrm.earlyad.neuroenriched.hipp))[, "Value"])
stderror.earlyad.neuroenriched.hipp  <- unname(coef(summary(mmrm.earlyad.neuroenriched.hipp))[, "Std.Error"])

#Neuroenriched Tau+
coefs.earlyad.neuroenriched.tplus.hipp     <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.hipp))[, "Value"])
stderror.earlyad.neuroenriched.tplus.hipp  <- unname(coef(summary(mmrm.earlyad.neuroenriched.tplus.hipp))[, "Std.Error"])



######################## MPACC ######################## 
coefs.earlyad.mpacc     <- unname(coef(summary(mmrm.earlyad.mpacc))[, "Value"])
stderror.earlyad.mpacc  <- unname(coef(summary(mmrm.earlyad.mpacc))[, "Std.Error"])

#Enriched for TAU
coefs.earlyad.tplus.mpacc     <- unname(coef(summary(mmrm.earlyad.tplus.mpacc))[, "Value"])
stderror.earlyad.tplus.mpacc  <- unname(coef(summary(mmrm.earlyad.tplus.mpacc))[, "Std.Error"])

#Neuroenriched
coefs.earlyad.neuroenriched.mpacc     <- unname(coef(summary(mmrm.earlyad.neuroenriched.mpacc))[, "Value"])
stderror.earlyad.neuroenriched.mpacc  <- unname(coef(summary(mmrm.earlyad.neuroenriched.mpacc))[, "Std.Error"])

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

#Neuroenriched
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

#Neuroenriched
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

#Neuroenriched
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
######################################## Simulations #################################################


######################## ADAS13 ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.earlyad.adas13 <- StratifyContVar(cs.earlyad.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.adas13 <- StratifyContVar(cs.earlyad.tplus.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.adas13 <- StratifyContVar(cs.earlyad.neuroenriched.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.tplus.adas13 <- StratifyContVar(cs.earlyad.neuroenriched.tplus.adas13, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))




#Assign treatment and placebo by blocks
long.earlyad.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.adas13, long.earlyad.adas13, no.prop = TRUE)
long.earlyad.tplus.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.adas13, long.earlyad.tplus.adas13, no.prop = TRUE)
long.earlyad.neuroenriched.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.adas13, long.earlyad.neuroenriched.adas13, no.prop = TRUE)
long.earlyad.neuroenriched.tplus.adas13_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.tplus.adas13, long.earlyad.neuroenriched.tplus.adas13, no.prop = TRUE)

######################## HIPPOCAMPUS ######################## 

# No Enrichment/Overall Decline
# Categorize continuous variables for block randomization
cs.earlyad.hipp <- StratifyContVar(cs.earlyad.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.hipp <- StratifyContVar(cs.earlyad.tplus.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.hipp <- StratifyContVar(cs.earlyad.neuroenriched.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.tplus.hipp <- StratifyContVar(cs.earlyad.neuroenriched.tplus.hipp, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))





#Assign treatment and placebo by blocks
long.earlyad.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.hipp, long.earlyad.hipp, no.prop = TRUE)
long.earlyad.tplus.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.hipp, long.earlyad.tplus.hipp, no.prop = TRUE)
long.earlyad.neuroenriched.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.hipp, long.earlyad.neuroenriched.hipp, no.prop = TRUE)
long.earlyad.neuroenriched.tplus.hipp_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.tplus.hipp, long.earlyad.neuroenriched.tplus.hipp, no.prop = TRUE)


######################## MPACC ######################## 

# Categorize continuous variables for block randomization
cs.earlyad.mpacc <- StratifyContVar(cs.earlyad.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.tplus.mpacc <- StratifyContVar(cs.earlyad.tplus.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.mpacc <- StratifyContVar(cs.earlyad.neuroenriched.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))
cs.earlyad.neuroenriched.tplus.mpacc <- StratifyContVar(cs.earlyad.neuroenriched.tplus.mpacc, stratcols = c("PTEDUCAT_bl", "AGE_bl", "MMSE_bl"))





#Assign treatment and placebo by blocks
long.earlyad.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.mpacc, long.earlyad.mpacc, no.prop = TRUE)
long.earlyad.tplus.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.tplus.mpacc, long.earlyad.tplus.mpacc, no.prop = TRUE)
long.earlyad.neuroenriched.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.mpacc, long.earlyad.neuroenriched.mpacc, no.prop = TRUE)
long.earlyad.neuroenriched.tplus.mpacc_with_treatment <- RandomizeTreatment2(cs.earlyad.neuroenriched.tplus.mpacc, long.earlyad.neuroenriched.tplus.mpacc, no.prop = TRUE)


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
#random intercept
formula.earlyad.adas13 <- "TOTAL13 ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"

model.earlyad.adas13 <- MapLmer(newdata = long.earlyad.adas13_with_treatment,
                                formula.model = formula.earlyad.adas13)
model.earlyad.tplus.adas13 <- MapLmer(newdata = long.earlyad.tplus.adas13_with_treatment,
                                      formula.model = formula.earlyad.adas13)
model.earlyad.neuroenriched.adas13 <- MapLmer(newdata = long.earlyad.neuroenriched.adas13_with_treatment,
                                              formula.model = formula.earlyad.adas13)
model.earlyad.neuroenriched.tplus.adas13 <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.adas13_with_treatment,
                                                    formula.model = formula.earlyad.adas13)


formula.earlyad.adas13.simulation             <- "TOTAL13 ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
simulation.model.earlyad.adas13               <- BuildSimulationModelNoPath(model.earlyad.adas13, formula.earlyad.adas13.simulation, long.earlyad.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.adas13         <- BuildSimulationModelNoPath(model.earlyad.tplus.adas13, formula.earlyad.adas13.simulation, long.earlyad.tplus.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.adas13 <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.adas13, formula.earlyad.adas13.simulation, long.earlyad.neuroenriched.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.adas13 <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.adas13, formula.earlyad.adas13.simulation, long.earlyad.neuroenriched.tplus.adas13_with_treatment, "not-controlled")


dpm.data.earlyad.adas13 <- get_model_data(simulation.model.earlyad.adas13, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.adas13 <- get_model_data(simulation.model.earlyad.tplus.adas13, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.adas13 <- get_model_data(simulation.model.earlyad.neuroenriched.adas13, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.tplus.adas13 <- get_model_data(simulation.model.earlyad.neuroenriched.tplus.adas13, terms = "new_time", type = "int") 

dpm.data.earlyad.adas13$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.adas13$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.adas13$Enrichment <- "No Copathologies"
dpm.data.earlyad.neuroenriched.tplus.adas13$Enrichment <- "No Copathologies Tau+"


fulldpm.earlyad.adas13 <- rbind(dpm.data.earlyad.adas13,
                                dpm.data.earlyad.tplus.adas13,
                                dpm.data.earlyad.neuroenriched.adas13)
fulldpm.earlyad.adas13$Treatment <- fulldpm.earlyad.adas13$group


fulldpm.earlyad.adas13.neuro <- rbind(dpm.data.earlyad.neuroenriched.adas13,
                                      dpm.data.earlyad.neuroenriched.tplus.adas13)

fulldpm.earlyad.adas13.neuro$Treatment <- fulldpm.earlyad.adas13.neuro$group






#random intercept and slope
formula.earlyad.adas13.rs <- "TOTAL13 ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
long.earlyad.adas13_with_treatment <- StratifyContinuous(long.earlyad.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.tplus.adas13_with_treatment <- StratifyContinuous(long.earlyad.tplus.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.neuroenriched.tplus.adas13_with_treatment <- StratifyContinuous(long.earlyad.neuroenriched.tplus.adas13_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))

model.earlyad.adas13.rs <- MapLmer(newdata = long.earlyad.adas13_with_treatment,
                                   formula.model = formula.earlyad.adas13.rs)

model.earlyad.tplus.adas13.rs <- MapLmer(newdata = long.earlyad.tplus.adas13_with_treatment,
                                         formula.model = formula.earlyad.adas13.rs)
model.earlyad.neuroenriched.adas13.rs <- MapLmer(newdata = long.earlyad.neuroenriched.adas13_with_treatment,
                                                 formula.model = formula.earlyad.adas13.rs)
model.earlyad.neuroenriched.tplus.adas13.rs <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.adas13_with_treatment,
                                                       formula.model = formula.earlyad.adas13.rs)

formula.earlyad.adas13.simulation.rs             <- "TOTAL13 ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
simulation.model.earlyad.adas13.rs               <- BuildSimulationModelNoPath(model.earlyad.adas13.rs, formula.earlyad.adas13.simulation.rs, long.earlyad.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.adas13.rs         <- BuildSimulationModelNoPath(model.earlyad.tplus.adas13.rs, formula.earlyad.adas13.simulation.rs, long.earlyad.tplus.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.adas13.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.adas13.rs, formula.earlyad.adas13.simulation.rs, long.earlyad.neuroenriched.adas13_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.adas13.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.adas13.rs, formula.earlyad.adas13.simulation.rs, long.earlyad.neuroenriched.tplus.adas13_with_treatment, "not-controlled")





dpm.data.earlyad.adas13.rs <- get_model_data(simulation.model.earlyad.adas13.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.adas13.rs <- get_model_data(simulation.model.earlyad.tplus.adas13.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.adas13.rs <- get_model_data(simulation.model.earlyad.neuroenriched.adas13.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.tplus.adas13.rs <- get_model_data(simulation.model.earlyad.neuroenriched.tplus.adas13.rs, terms = "new_time", type = "int") 

dpm.data.earlyad.adas13.rs$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.adas13.rs$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.adas13.rs$Enrichment <- "No Copathologies"
dpm.data.earlyad.neuroenriched.tplus.adas13.rs$Enrichment <- "No Copathologies Tau+"


dpm.data.earlyad.adas13$Treatment <- dpm.data.earlyad.adas13$group
dpm.data.earlyad.tplus.adas13$Treatment <- dpm.data.earlyad.adas13$group

fulldpm.earlyad.adas13.rs <- rbind(dpm.data.earlyad.adas13.rs,
                                   dpm.data.earlyad.tplus.adas13.rs,
                                   dpm.data.earlyad.neuroenriched.adas13.rs)
fulldpm.earlyad.adas13.rs$Treatment <- fulldpm.earlyad.adas13.rs$group


fulldpm.earlyad.adas13.rs.neuro <- rbind(dpm.data.earlyad.neuroenriched.adas13.rs,
                                         dpm.data.earlyad.neuroenriched.tplus.adas13.rs)

fulldpm.earlyad.adas13.rs.neuro$Treatment <- fulldpm.earlyad.adas13.rs.neuro$group


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

dpms.earlyad.adas13 <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.adas13.rs,
                                     "Tau+" = simulation.model.earlyad.tplus.adas13.rs,
                                     "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.adas13.rs), ylab="ADAS13", ylim.low = 19, ylim.high = 32.5)

dpms.earlyad.tau.adas13 <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.adas13.tau.rs,
                                     "Tau+" = simulation.model.earlyad.tplus.adas13.tau.rs,
                                     "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.adas13.tau.rs), ylab="ADAS13", ylim.low = 19, ylim.high = 32.5)




dpms.earlyad.tau.adas13
######################## HIPPOCAMPUS ######################## 
formula.earlyad.hipp <- "hipp_average ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
model.earlyad.hipp <- MapLmer(newdata = long.earlyad.hipp_with_treatment,
                              formula.model = formula.earlyad.hipp)
model.earlyad.tplus.hipp <- MapLmer(newdata = long.earlyad.tplus.hipp_with_treatment,
                                    formula.model = formula.earlyad.hipp)
model.earlyad.neuroenriched.hipp <- MapLmer(newdata = long.earlyad.neuroenriched.hipp_with_treatment,
                                            formula.model = formula.earlyad.hipp)
model.earlyad.neuroenriched.tplus.hipp <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.hipp_with_treatment,
                                                  formula.model = formula.earlyad.hipp)


formula.earlyad.hipp.simulation             <- "hipp_average ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
simulation.model.earlyad.hipp               <- BuildSimulationModelNoPath(model.earlyad.hipp, formula.earlyad.hipp.simulation, long.earlyad.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.hipp         <- BuildSimulationModelNoPath(model.earlyad.tplus.hipp, formula.earlyad.hipp.simulation, long.earlyad.tplus.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.hipp <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.hipp, formula.earlyad.hipp.simulation, long.earlyad.neuroenriched.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.hipp <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.hipp, formula.earlyad.hipp.simulation, long.earlyad.neuroenriched.tplus.hipp_with_treatment, "not-controlled")


dpm.data.earlyad.hipp <- get_model_data(simulation.model.earlyad.hipp, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.hipp <- get_model_data(simulation.model.earlyad.tplus.hipp, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.hipp <- get_model_data(simulation.model.earlyad.neuroenriched.hipp, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.tplus.hipp <- get_model_data(simulation.model.earlyad.neuroenriched.tplus.hipp, terms = "new_time", type = "int") 

dpm.data.earlyad.hipp$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.hipp$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.hipp$Enrichment <- "No Copathologies"
dpm.data.earlyad.neuroenriched.tplus.hipp$Enrichment <- "No Copathologies Tau+"

fulldpm.earlyad.hipp <- rbind(dpm.data.earlyad.hipp,
                              dpm.data.earlyad.tplus.hipp,
                              dpm.data.earlyad.neuroenriched.hipp)
fulldpm.earlyad.hipp$Treatment <- fulldpm.earlyad.hipp$group


fulldpm.earlyad.hipp.neuro <- rbind(dpm.data.earlyad.neuroenriched.hipp,
                                    dpm.data.earlyad.neuroenriched.tplus.hipp)

fulldpm.earlyad.hipp.neuro$Treatment <- fulldpm.earlyad.hipp.neuro$group





#random intercept and slope
formula.earlyad.hipp.rs <- "hipp_average ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
long.earlyad.hipp_with_treatment <- StratifyContinuous(long.earlyad.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.tplus.hipp_with_treatment <- StratifyContinuous(long.earlyad.tplus.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.neuroenriched.tplus.hipp_with_treatment <- StratifyContinuous(long.earlyad.neuroenriched.tplus.hipp_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))

model.earlyad.hipp.rs <- MapLmer(newdata = long.earlyad.hipp_with_treatment,
                                 formula.model = formula.earlyad.hipp.rs)
model.earlyad.tplus.hipp.rs <- MapLmer(newdata = long.earlyad.tplus.hipp_with_treatment,
                                       formula.model = formula.earlyad.hipp.rs)
model.earlyad.neuroenriched.hipp.rs <- MapLmer(newdata = long.earlyad.neuroenriched.hipp_with_treatment,
                                               formula.model = formula.earlyad.hipp.rs)
model.earlyad.neuroenriched.tplus.hipp.rs <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.hipp_with_treatment,
                                                     formula.model = formula.earlyad.hipp.rs)


formula.earlyad.hipp.simulation.rs             <- "hipp_average ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
simulation.model.earlyad.hipp.rs               <- BuildSimulationModelNoPath(model.earlyad.hipp.rs, formula.earlyad.hipp.simulation.rs, long.earlyad.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.hipp.rs         <- BuildSimulationModelNoPath(model.earlyad.tplus.hipp.rs, formula.earlyad.hipp.simulation.rs, long.earlyad.tplus.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.hipp.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.hipp.rs, formula.earlyad.hipp.simulation.rs, long.earlyad.neuroenriched.hipp_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.hipp.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.hipp.rs, formula.earlyad.hipp.simulation.rs, long.earlyad.neuroenriched.tplus.hipp_with_treatment, "not-controlled")


dpm.data.earlyad.hipp.rs <- get_model_data(simulation.model.earlyad.hipp.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.hipp.rs <- get_model_data(simulation.model.earlyad.tplus.hipp.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.hipp.rs <- get_model_data(simulation.model.earlyad.neuroenriched.hipp.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.tplus.hipp.rs <- get_model_data(simulation.model.earlyad.neuroenriched.tplus.hipp.rs, terms = "new_time", type = "int") 

dpm.data.earlyad.hipp.rs$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.hipp.rs$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.hipp.rs$Enrichment <- "No Copathologies"
dpm.data.earlyad.neuroenriched.tplus.hipp.rs$Enrichment <- "No Copathologies Tau+"

fulldpm.earlyad.hipp.rs <- rbind(dpm.data.earlyad.hipp.rs,
                                 dpm.data.earlyad.tplus.hipp.rs,
                                 dpm.data.earlyad.neuroenriched.hipp.rs)
fulldpm.earlyad.hipp.rs$Treatment <- fulldpm.earlyad.hipp.rs$group


fulldpm.earlyad.hipp.rs.neuro <- rbind(dpm.data.earlyad.neuroenriched.hipp.rs,
                                       dpm.data.earlyad.neuroenriched.tplus.hipp.rs)

fulldpm.earlyad.hipp.rs.neuro$Treatment <- fulldpm.earlyad.hipp.rs.neuro$group


relcontr.earlyad.hipp <- GetRelContributions(model.earlyad.hipp.rs, long.earlyad.hipp_with_treatment)
relcontr.earlyad.tplus.hipp <- GetRelContributions(model.earlyad.tplus.hipp.rs, long.earlyad.tplus.hipp_with_treatment)
relcontr.earlyad.neuroenriched.tplus.hipp        <- GetRelContributions(model.earlyad.neuroenriched.tplus.hipp.rs, long.earlyad.neuroenriched.tplus.hipp_with_treatment)

relcontrcount.earlyad.hipp <- BuildNeuroCountPlot(relcontr.earlyad.hipp)
relcontrcount.earlyad.tplus.hipp <- BuildNeuroCountPlot(relcontr.earlyad.tplus.hipp)






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




dpms.earlyad.hipp <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.hipp.rs,
                                     "Tau+" = simulation.model.earlyad.tplus.hipp.rs,
                                     "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.hipp.rs), ylab="Hippocampus", ylim.low = 2800, ylim.high = 4000)

dpms.earlyad.tau.hipp <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.hipp.tau.rs,
                                         "Tau+" = simulation.model.earlyad.tplus.hipp.tau.rs,
                                         "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.hipp.tau.rs), ylab="Hippocampus", ylim.low = 2800, ylim.high = 4000)




######################## MPACC ######################## 
formula.earlyad.mpacc <- "mPACCtrailsB ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"

model.earlyad.mpacc <- MapLmer(newdata = long.earlyad.mpacc_with_treatment,
                               formula.model = formula.earlyad.mpacc)
model.earlyad.tplus.mpacc <- MapLmer(newdata = long.earlyad.tplus.mpacc_with_treatment,
                                     formula.model = formula.earlyad.mpacc)
model.earlyad.neuroenriched.mpacc <- MapLmer(newdata = long.earlyad.neuroenriched.mpacc_with_treatment,
                                             formula.model = formula.earlyad.mpacc)
model.earlyad.neuroenriched.tplus.mpacc <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.mpacc_with_treatment,
                                                   formula.model = formula.earlyad.mpacc)


formula.earlyad.mpacc.simulation             <- "mPACCtrailsB ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID)"
simulation.model.earlyad.mpacc               <- BuildSimulationModelNoPath(model.earlyad.mpacc, formula.earlyad.mpacc.simulation, long.earlyad.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.mpacc         <- BuildSimulationModelNoPath(model.earlyad.tplus.mpacc, formula.earlyad.mpacc.simulation, long.earlyad.tplus.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.mpacc <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.mpacc, formula.earlyad.mpacc.simulation, long.earlyad.neuroenriched.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.mpacc <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.mpacc, formula.earlyad.mpacc.simulation, long.earlyad.neuroenriched.tplus.mpacc_with_treatment, "not-controlled")


dpm.data.earlyad.mpacc <- get_model_data(simulation.model.earlyad.mpacc, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.mpacc <- get_model_data(simulation.model.earlyad.tplus.mpacc, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.mpacc <- get_model_data(simulation.model.earlyad.neuroenriched.mpacc, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.tplus.mpacc <- get_model_data(simulation.model.earlyad.neuroenriched.tplus.mpacc, terms = "new_time", type = "int") 

dpm.data.earlyad.mpacc$Enrichment <- "No Enrichment"
dpm.data.earlyad.tplus.mpacc$Enrichment <- "Tau+"
dpm.data.earlyad.neuroenriched.mpacc$Enrichment <- "No Copathologies"
dpm.data.earlyad.neuroenriched.tplus.mpacc$Enrichment <- "No Copathologies Tau+"

fulldpm.earlyad.mpacc <- rbind(dpm.data.earlyad.mpacc,
                               dpm.data.earlyad.tplus.mpacc,
                               dpm.data.earlyad.neuroenriched.mpacc)
fulldpm.earlyad.mpacc$Treatment <- fulldpm.earlyad.mpacc$group

fulldpm.earlyad.mpacc.neuro <- rbind(dpm.data.earlyad.neuroenriched.mpacc,
                                     dpm.data.earlyad.neuroenriched.tplus.mpacc)

fulldpm.earlyad.mpacc.neuro$Treatment <- fulldpm.earlyad.mpacc.neuro$group




#random intercept and slope
formula.earlyad.mpacc.rs <- "mPACCtrailsB ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
long.earlyad.mpacc_with_treatment <- StratifyContinuous(long.earlyad.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.tplus.mpacc_with_treatment <- StratifyContinuous(long.earlyad.tplus.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
long.earlyad.neuroenriched.tplus.mpacc_with_treatment <- StratifyContinuous(long.earlyad.neuroenriched.tplus.mpacc_with_treatment, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))

model.earlyad.mpacc.rs <- MapLmer(newdata = long.earlyad.mpacc_with_treatment,
                                  formula.model = formula.earlyad.mpacc.rs)
model.earlyad.tplus.mpacc.rs <- MapLmer(newdata = long.earlyad.tplus.mpacc_with_treatment,
                                        formula.model = formula.earlyad.mpacc.rs)
model.earlyad.neuroenriched.mpacc.rs <- MapLmer(newdata = long.earlyad.neuroenriched.mpacc_with_treatment,
                                                formula.model = formula.earlyad.mpacc.rs)
model.earlyad.neuroenriched.tplus.mpacc.rs <- MapLmer(newdata = long.earlyad.neuroenriched.tplus.mpacc_with_treatment,
                                                      formula.model = formula.earlyad.mpacc.rs)


formula.earlyad.mpacc.simulation.rs                   <- "mPACCtrailsB ~ new_time*treat + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1 + new_time|RID)"
simulation.model.earlyad.mpacc.rs                     <- BuildSimulationModelNoPath(model.earlyad.mpacc.rs, formula.earlyad.mpacc.simulation.rs, long.earlyad.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.tplus.mpacc.rs               <- BuildSimulationModelNoPath(model.earlyad.tplus.mpacc.rs, formula.earlyad.mpacc.simulation.rs, long.earlyad.tplus.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.mpacc.rs       <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.mpacc.rs, formula.earlyad.mpacc.simulation.rs, long.earlyad.neuroenriched.mpacc_with_treatment, "not-controlled")
simulation.model.earlyad.neuroenriched.tplus.mpacc.rs <- BuildSimulationModelNoPath(model.earlyad.neuroenriched.tplus.mpacc.rs, formula.earlyad.mpacc.simulation.rs, long.earlyad.neuroenriched.tplus.mpacc_with_treatment, "not-controlled")


dpm.data.earlyad.mpacc.rs                     <- get_model_data(simulation.model.earlyad.mpacc.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.tplus.mpacc.rs               <- get_model_data(simulation.model.earlyad.tplus.mpacc.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.mpacc.rs       <- get_model_data(simulation.model.earlyad.neuroenriched.mpacc.rs, terms = "new_time", type = "int") 
dpm.data.earlyad.neuroenriched.tplus.mpacc.rs <- get_model_data(simulation.model.earlyad.neuroenriched.tplus.mpacc.rs, terms = "new_time", type = "int") 

dpm.data.earlyad.mpacc.rs$Enrichment                     <- "No Enrichment"
dpm.data.earlyad.tplus.mpacc.rs$Enrichment               <- "Tau+"
dpm.data.earlyad.neuroenriched.mpacc.rs$Enrichment       <- "No Copathologies"
dpm.data.earlyad.neuroenriched.tplus.mpacc.rs$Enrichment <- "No Copathologies Tau+"

fulldpm.earlyad.mpacc.rs <- rbind(dpm.data.earlyad.mpacc.rs,
                                  dpm.data.earlyad.tplus.mpacc.rs,
                                  dpm.data.earlyad.neuroenriched.mpacc.rs)
fulldpm.earlyad.mpacc.rs$Treatment <- fulldpm.earlyad.mpacc.rs$group


fulldpm.earlyad.mpacc.rs.neuro     <- rbind(dpm.data.earlyad.neuroenriched.mpacc.rs,
                                            dpm.data.earlyad.neuroenriched.tplus.mpacc.rs)

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



dpms.earlyad.mpacc <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.mpacc.rs,
                                   "Tau+" = simulation.model.earlyad.tplus.mpacc.rs,
                                   "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.mpacc.rs), ylab="mPACCtrailsB", ylim.low = -19, ylim.high = -7)

dpms.earlyad.tau.mpacc <- DPMPlots(list("No Enrichment" = simulation.model.earlyad.mpacc.tau.rs,
                                       "Tau+" = simulation.model.earlyad.tplus.mpacc.tau.rs,
                                       "No Copathologies Tau+" = simulation.model.earlyad.neuroenriched.tplus.mpacc.tau.rs), ylab="mPACCtrailsB",  ylim.low = -19, ylim.high = -7)




############### saving RDS files ###############

if(FALSE) {
  
  
  earlyad.adas13.modeling.list <- list("sim.data" = list("unenriched"=long.earlyad.adas13_with_treatment,
                                                         "taupos"=long.earlyad.tplus.adas13_with_treatment,
                                                         "nocopath"=long.earlyad.neuroenriched.adas13_with_treatment),
                                       "model" = list("unenriched"=simulation.model.earlyad.adas13,
                                                      "taupos"=simulation.model.earlyad.tplus.adas13,
                                                      "nocopath"=simulation.model.earlyad.neuroenriched.adas13),
                                       "formula" = replicate(3,formula.earlyad.adas13.simulation, simplify = FALSE),
                                       "compare_str" = replicate(3,formula.earlyad.adas13, simplify = FALSE),
                                       "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                       "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE),
                                       "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  
  earlyad.adas13.modeling.list.rs <- list("sim.data" = list("unenriched"=long.earlyad.adas13_with_treatment,
                                                            "taupos"=long.earlyad.tplus.adas13_with_treatment,
                                                            "nocopath_tau+"=long.earlyad.neuroenriched.tplus.adas13_with_treatment),
                                          "model" = list("unenriched"    = simulation.model.earlyad.adas13.rs,
                                                         "taupos"        = simulation.model.earlyad.tplus.adas13.rs,
                                                         "nocopath_tau+" = simulation.model.earlyad.neuroenriched.tplus.adas13.rs),
                                          "formula" = replicate(3,formula.earlyad.adas13.simulation.rs, simplify = FALSE),
                                          "compare_str" = replicate(3,formula.earlyad.adas13.rs, simplify = FALSE),
                                          "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                          "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE),
                                          "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  
  earlyad.adas13.modeling.list.tau.rs <- list("sim.data" = list("unenriched"=long.earlyad.adas13_with_treatment,
                                                            "taupos"=long.earlyad.tplus.adas13_with_treatment,
                                                            "nocopath_tau+"=long.earlyad.neuroenriched.tplus.adas13_with_treatment),
                                          "model" = list("unenriched"    = simulation.model.earlyad.adas13.tau.rs,
                                                         "taupos"        = simulation.model.earlyad.tplus.adas13.tau.rs,
                                                         "nocopath_tau+" = simulation.model.earlyad.neuroenriched.tplus.adas13.tau.rs),
                                          "formula" = replicate(3,formula.earlyad.adas13.simulation.rs, simplify = FALSE),
                                          "compare_str" = replicate(3,formula.earlyad.adas13.rs, simplify = FALSE),
                                          "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                          "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE),
                                          "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  earlyad.hipp.modeling.list <- list("sim.data" = list("unenriched"=long.earlyad.hipp_with_treatment,
                                                       "taupos"=long.earlyad.tplus.hipp_with_treatment,
                                                       "nocopath"=long.earlyad.neuroenriched.hipp_with_treatment),
                                     "model" = list("unenriched"=simulation.model.earlyad.hipp,
                                                    "taupos"=simulation.model.earlyad.tplus.hipp,
                                                    "nocopath"=simulation.model.earlyad.neuroenriched.hipp),
                                     "formula" = replicate(3,formula.earlyad.hipp.simulation, simplify = FALSE),
                                     "compare_str" = replicate(3,formula.earlyad.hipp, simplify = FALSE),
                                     "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                     "yaxislab_dpm" = replicate(3, "Hippocampus", simplify = FALSE),
                                     "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  
  earlyad.hipp.modeling.list.rs <- list("sim.data" = list("unenriched"=long.earlyad.hipp_with_treatment,
                                                          "taupos"=long.earlyad.tplus.hipp_with_treatment,
                                                          "nocopath_tau+"=long.earlyad.neuroenriched.tplus.hipp_with_treatment),
                                        "model" = list("unenriched"=simulation.model.earlyad.hipp.rs,
                                                       "taupos"=simulation.model.earlyad.tplus.hipp.rs,
                                                       "nocopath_tau+"=simulation.model.earlyad.neuroenriched.tplus.hipp.rs),
                                        "formula" = replicate(3,formula.earlyad.hipp.simulation.rs, simplify = FALSE),
                                        "compare_str" = replicate(3,formula.earlyad.hipp.rs, simplify = FALSE),
                                        "breaks" = replicate(3, seq(50, 700, by=25), simplify = FALSE),
                                        "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE),
                                        "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  earlyad.hipp.modeling.list.tau.rs <- list("sim.data" = list("unenriched"=long.earlyad.hipp_with_treatment,
                                                          "taupos"=long.earlyad.tplus.hipp_with_treatment,
                                                          "nocopath_tau+"=long.earlyad.neuroenriched.tplus.hipp_with_treatment),
                                        "model" = list("unenriched"=simulation.model.earlyad.hipp.tau.rs,
                                                       "taupos"=simulation.model.earlyad.tplus.hipp.tau.rs,
                                                       "nocopath_tau+"=simulation.model.earlyad.neuroenriched.tplus.hipp.tau.rs),
                                        "formula" = replicate(3,formula.earlyad.hipp.simulation.rs, simplify = FALSE),
                                        "compare_str" = replicate(3,formula.earlyad.hipp.rs, simplify = FALSE),
                                        "breaks"= replicate(3, seq(50, 700, by=25), simplify = FALSE),
                                        "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE),
                                        "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  
  earlyad.mpacc.modeling.list <- list("sim.data" = list("unenriched"=long.earlyad.mpacc_with_treatment,
                                                        "taupos"=long.earlyad.tplus.mpacc_with_treatment,
                                                        "nocopath"=long.earlyad.neuroenriched.mpacc_with_treatment),
                                      "model" = list("unenriched"=simulation.model.earlyad.mpacc,
                                                     "taupos"=simulation.model.earlyad.tplus.mpacc,
                                                     "nocopath"=simulation.model.earlyad.neuroenriched.mpacc),
                                      "formula" = replicate(3,formula.earlyad.mpacc.simulation, simplify = FALSE),
                                      "compare_str" = replicate(3,formula.earlyad.mpacc, simplify = FALSE),
                                      "breaks" = replicate(3, seq(100, 1000, by=100), simplify = FALSE),
                                      "yaxislab_dpm" = replicate(3, "mPACCtrailsB", simplify = FALSE),
                                      "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  
  
  earlyad.mpacc.modeling.list.rs <- list("sim.data" = list("unenriched"=long.earlyad.mpacc_with_treatment,
                                                           "taupos"=long.earlyad.tplus.mpacc_with_treatment,
                                                           "nocopath_tau+"=long.earlyad.neuroenriched.tplus.mpacc_with_treatment),
                                         "model" = list("unenriched"=simulation.model.earlyad.mpacc.rs,
                                                        "taupos"=simulation.model.earlyad.tplus.mpacc.rs,
                                                        "nocopath_tau+"=simulation.model.earlyad.neuroenriched.tplus.mpacc.rs),
                                         "formula" = replicate(3,formula.earlyad.mpacc.simulation.rs, simplify = FALSE),
                                         "compare_str" = replicate(3,formula.earlyad.mpacc.rs, simplify = FALSE),
                                         "breaks"= replicate(3, seq(50, 700, by=25), simplify = FALSE),
                                         "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE),
                                         "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  earlyad.mpacc.modeling.list.tau.rs <- list("sim.data" = list("unenriched"=long.earlyad.mpacc_with_treatment,
                                                           "taupos"=long.earlyad.tplus.mpacc_with_treatment,
                                                           "nocopath_tau+"=long.earlyad.neuroenriched.tplus.mpacc_with_treatment),
                                         "model" = list("unenriched"=simulation.model.earlyad.mpacc.tau.rs,
                                                        "taupos"=simulation.model.earlyad.tplus.mpacc.tau.rs,
                                                        "nocopath_tau+"=simulation.model.earlyad.neuroenriched.tplus.mpacc.tau.rs),
                                         "formula" = replicate(3,formula.earlyad.mpacc.simulation.rs, simplify = FALSE),
                                         "compare_str" = replicate(3,formula.earlyad.mpacc.rs, simplify = FALSE),
                                         "breaks"= replicate(3, seq(50, 700, by=25), simplify = FALSE),
                                         "yaxislab_dpm" = replicate(3, "ADAS13", simplify = FALSE),
                                         "return_dpm" = replicate(3, FALSE, simplify = FALSE))
  
  
  
  
  saveRDS(list("adas13" = list("formula_largemodel" = formula.earlyad.adas13.simulation.rs , 
                          "largemodel" = simulation.model.earlyad.adas13.rs , 
                          "formula_smallmodel" = formula.earlyad.adas13 , 
                          "smallmodel" = model.earlyad.adas13.rs , 
                          "sample_sizes" = seq(100, 300, by=100), 
                          "nsim" =100),
               "hipp" = list("formula_largemodel" = formula.earlyad.hipp.simulation.rs , 
                             "largemodel" = simulation.model.earlyad.hipp.rs , 
                             "formula_smallmodel" = formula.earlyad.hipp, 
                             "smallmodel" = model.earlyad.hipp.rs , 
                             "sample_sizes" = seq(100, 200, by=50), 
                             "nsim" = 100 ),
               "mpacc" = list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs , 
                              "largemodel" = simulation.model.earlyad.mpacc.rs , 
                              "formula_smallmodel" = formula.earlyad.mpacc , 
                              "smallmodel" = model.earlyad.mpacc.rs , 
                              "sample_sizes" = seq(100, 500, by=50), 
                              "nsim" = 100),
               "adas13_tau+" = list("formula_largemodel" = formula.earlyad.adas13.simulation.rs , 
                               "largemodel" = simulation.model.earlyad.tplus.adas13.rs , 
                               "formula_smallmodel" = formula.earlyad.adas13 , 
                               "smallmodel" = model.earlyad.tplus.adas13.rs , 
                               "sample_sizes" = seq(100, 200, by=50), 
                               "nsim" = 100),
               "hipp_tau+" = list("formula_largemodel" = formula.earlyad.hipp.simulation.rs , 
                                    "largemodel" = simulation.model.earlyad.tplus.hipp.rs , 
                                    "formula_smallmodel" = formula.earlyad.hipp , 
                                    "smallmodel" = model.earlyad.tplus.hipp.rs , 
                                    "sample_sizes" = seq(100, 300, by=50), 
                                    "nsim" = 100),
               "mpacc_tau+" = list("formula_largemodel" = formula.earlyad.mpacc.simulation.rs , 
                                    "largemodel" = simulation.model.earlyad.tplus.mpacc.rs , 
                                    "formula_smallmodel" = formula.earlyad.mpacc , 
                                    "smallmodel" = model.earlyad.tplus.mpacc.rs , 
                                    "sample_sizes" = seq(100, 400, by=100), 
                                    "nsim" = 1000)), "/Users/adamgabriellang/Desktop/clinical_trial_sim/listsformanualsim.rds")
  
  
  #random intercept
  saveRDS(earlyad.adas13.modeling.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.adas13_4sim.rds")
  saveRDS(earlyad.hipp.modeling.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.hipp_4sim.rds")
  saveRDS(earlyad.mpacc.modeling.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.mpacc_4sim.rds")
  
  
  #random slope and intercept
  saveRDS(earlyad.adas13.modeling.list.rs, "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.adas13_4sim_rs.rds")
  saveRDS(earlyad.hipp.modeling.list.rs, "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.hipp_4sim_rs.rds")
  saveRDS(earlyad.mpacc.modeling.list.rs, "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.mpacc_4sim_rs.rds")
  
  #random slope and intercept
  saveRDS(list("fulldecline" = earlyad.adas13.modeling.list.rs,
               "taudecline" = earlyad.adas13.modeling.list.tau.rs), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.adas13_4sim_rs_withtau.rds")
  saveRDS(list("fulldecline" = earlyad.hipp.modeling.list.rs,
               "taudecline" = earlyad.hipp.modeling.list.tau.rs), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.hipp_4sim_rs_withtau.rds")
  saveRDS(list("fulldecline" = earlyad.mpacc.modeling.list.rs,
               "taudecline" = earlyad.mpacc.modeling.list.tau.rs), "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.mpacc_4sim_rs_withtau.rds")
  
  


#### read fitted simulation models ####

#ADAS13

earlyad.adas13.fitted.withtau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.adas13_4sim_fitted_withtau.rds")
earlyad.adas13.fitted.withtau.fulldecline <- earlyad.adas13.fitted.withtau$fulldecline
earlyad.adas13.fitted.withtau.taudecline <- earlyad.adas13.fitted.withtau$taudecline

names(earlyad.adas13.fitted.withtau.fulldecline) <- c("Unenriched", "Tau+", "NoCpthTau+")
names(earlyad.adas13.fitted.withtau.taudecline) <- c("Unenriched", "Tau+", "NoCpthTau+")

earlyad.adas13.fitted.withtau.fulldecline.plot <- CombineSimPlots(earlyad.adas13.fitted.withtau.fulldecline, limits = seq(100, 1000, by=100))
earlyad.adas13.fitted.withtau.fulldecline.plot$plot <- earlyad.adas13.fitted.withtau.fulldecline.plot$plot+ ylim(0, 100)

earlyad.adas13.fitted.withtau.taudecline.plot <- CombineSimPlots(earlyad.adas13.fitted.withtau.taudecline, limits = seq(100, 1000, by=100))
earlyad.adas13.fitted.withtau.taudecline.plot$plot <- earlyad.adas13.fitted.withtau.taudecline.plot$plot + ylim(0, 100)

earlyad.adas13.fitted.withtau.fulldecline.plot$plot
earlyad.adas13.fitted.withtau.taudecline.plot$plot

earlyad.adas13.fitted.withtau.fulldecline.dataframe <- BuildSimulationDataTable(earlyad.adas13.fitted.withtau.fulldecline.plot)
earlyad.adas13.fitted.withtau.taudecline.dataframe <- BuildSimulationDataTable(earlyad.adas13.fitted.withtau.taudecline.plot)

#Hippocampus


earlyad.hipp.fitted.withtau <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.hipp_4sim_fitted_withtau.rds")
earlyad.hipp.fitted.withtau.fulldecline <- earlyad.hipp.fitted.withtau$fulldecline
earlyad.hipp.fitted.withtau.taudecline <- earlyad.hipp.fitted.withtau$taudecline

names(earlyad.hipp.fitted.withtau.fulldecline) <- c("Unenriched", "Tau+", "NoCpthTau+")
names(earlyad.hipp.fitted.withtau.taudecline) <- c("Unenriched", "Tau+", "NoCpthTau+")

earlyad.hipp.fitted.withtau.fulldecline.plot <- CombineSimPlots(earlyad.hipp.fitted.withtau.fulldecline, limits = seq(50, 700, by=25))
earlyad.hipp.fitted.withtau.fulldecline.plot$plot <- earlyad.hipp.fitted.withtau.fulldecline.plot$plot+ ylim(0, 100)

earlyad.hipp.fitted.withtau.taudecline.plot <- CombineSimPlots(earlyad.hipp.fitted.withtau.taudecline, limits = seq(50, 700, by=25))
earlyad.hipp.fitted.withtau.taudecline.plot$plot <- earlyad.hipp.fitted.withtau.taudecline.plot$plot + ylim(0, 100)

earlyad.hipp.fitted.withtau.fulldecline.plot$plot
earlyad.hipp.fitted.withtau.taudecline.plot$plot

earlyad.hipp.fitted.withtau.fulldecline.dataframe <- BuildSimulationDataTable(earlyad.hipp.fitted.withtau.fulldecline.plot)
earlyad.hipp.fitted.withtau.taudecline.dataframe <- BuildSimulationDataTable(earlyad.hipp.fitted.withtau.taudecline.plot)




#MPACC 

earlyad.mpacc.fitted.withtau1 <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.mpacc_4sim_fitted_withtau_1.rds")
earlyad.mpacc.fitted.withtau2 <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad.mpacc_4sim_fitted_withtau_2.rds")

earlyad.mpacc.fitted.withtau.fulldecline <- earlyad.mpacc.fitted.withtau1
earlyad.mpacc.fitted.withtau.taudecline <- earlyad.mpacc.fitted.withtau2

names(earlyad.mpacc.fitted.withtau.fulldecline) <- c("Unenriched", "Tau+", "NoCpthTau+")
names(earlyad.mpacc.fitted.withtau.taudecline) <- c("Unenriched", "Tau+", "NoCpthTau+")

earlyad.mpacc.fitted.withtau.fulldecline.plot <- CombineSimPlots(earlyad.mpacc.fitted.withtau.fulldecline, limits = seq(50, 700, by=25))
earlyad.mpacc.fitted.withtau.fulldecline.plot$plot <- earlyad.mpacc.fitted.withtau.fulldecline.plot$plot+ ylim(0, 100)

earlyad.mpacc.fitted.withtau.taudecline.plot <- CombineSimPlots(earlyad.mpacc.fitted.withtau.taudecline, limits = seq(50, 700, by=25))
earlyad.mpacc.fitted.withtau.taudecline.plot$plot <- earlyad.mpacc.fitted.withtau.taudecline.plot$plot + ylim(0, 100)

earlyad.mpacc.fitted.withtau.fulldecline.plot$plot
earlyad.mpacc.fitted.withtau.taudecline.plot$plot
earlyad.mpacc.fitted.withtau.taudecline.plot$fullstats
earlyad.mpacc.fitted.withtau.fulldecline.dataframe <- BuildSimulationDataTable(earlyad.mpacc.fitted.withtau.fulldecline.plot)
earlyad.mpacc.fitted.withtau.taudecline.dataframe <- BuildSimulationDataTable(earlyad.mpacc.fitted.withtau.taudecline.plot)


saveRDS(list("adas" = list(earlyad.adas13.fitted.withtau.fulldecline.dataframe,
                           earlyad.adas13.fitted.withtau.taudecline.dataframe),
             "hipp" = list(earlyad.hipp.fitted.withtau.fulldecline.dataframe,
                           earlyad.hipp.fitted.withtau.taudecline.dataframe),
             "mpacc" = list(earlyad.mpacc.fitted.withtau.fulldecline.dataframe,
                            earlyad.mpacc.fitted.withtau.taudecline.dataframe)), "/Users/adamgabriellang/Desktop/clinical_trial_sim/simulationdatalist.rds")


}
