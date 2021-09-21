library(nlme)

#read in early ad data

earlyadcohorts       <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts.rds")
earlyadcohorts.tplus <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyagecohorts_tplus.rds")

######################## ADAS13 ######################## 
cs.earlyad.adas13    <- earlyadcohorts$cs$ADAS13
long.earlyad.adas13  <-  earlyadcohorts$long$ADAS13

#Tau+
cs.earlyad.tplus.adas13    <- earlyadcohorts.tplus$cs$ADAS13
long.earlyad.tplus.adas13  <-  earlyadcohorts.tplus$long$ADAS13


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




#fit repeated measures model on longitudinal data
# want observed scores, just use time and dont control for other variables

######################## ADAS13 ######################## 
mmrm.earlyad.adas13       <- lmer(TOTAL13 ~ factor(new_time) + (1|RID), data = long.earlyad.adas13)

#Enriched for TAU
mmrm.earlyad.tplus.adas13 <- lmer(TOTAL13 ~ factor(new_time) + (1|RID), data = long.earlyad.tplus.adas13)


######################## HIPPOCAMPUS ######################## 
mmrm.earlyad.hipp         <- lmer(hipp_average ~ factor(new_time) + (1|RID), data = long.earlyad.hipp)
#Enriched for TAU
mmrm.earlyad.tplus.hipp   <- lmer(hipp_average ~ factor(new_time) + (1|RID), data = long.earlyad.tplus.hipp)


#get coefficients and standard errors

######################## ADAS13 ######################## 
coefs.earlyad.adas13     <- unname(coef(summary(mmrm.earlyad.adas13))[, "Estimate"])
stderror.earlyad.adas13  <- unname(coef(summary(mmrm.earlyad.adas13))[, "Std. Error"])
#Enriched for TAU
coefs.earlyad.tplus.adas13     <- unname(coef(summary(mmrm.earlyad.tplus.adas13))[, "Estimate"])
stderror.earlyad.tplus.adas13  <- unname(coef(summary(mmrm.earlyad.tplus.adas13))[, "Std. Error"])

######################## HIPPOCAMPUS ######################## 
coefs.earlyad.hipp     <- unname(coef(summary(mmrm.earlyad.hipp))[, "Estimate"])
stderror.earlyad.hipp  <- unname(coef(summary(mmrm.earlyad.hipp))[, "Std. Error"])
#Enriched for TAU
coefs.earlyad.tplus.hipp     <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Estimate"])
stderror.earlyad.tplus.hipp  <- unname(coef(summary(mmrm.earlyad.tplus.hipp))[, "Std. Error"])





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



## Plots
# Combine non-enriched and tau enriched data for plotting

######################## ADAS13 ######################## 
earlyad.adas13.plotting.data <- rbind(earlyad.adas13.plotting.data,
                                      earlyad.tplus.adas13.plotting.data)
earlyad.adas13.mmrm.plot <- ggplot(earlyad.adas13.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.adas13.mmrm.plot <- earlyad.adas13.mmrm.plot + xlab("Time (Years)") + ylab("Mean ADAS-13 (95% CI)") + ylim(22, 35)
earlyad.adas13.mmrm.plot


######################## HIPPOCAMPUS ######################## 
earlyad.hipp.plotting.data <- rbind(earlyad.hipp.plotting.data,
                                    earlyad.tplus.hipp.plotting.data)
earlyad.hipp.mmrm.plot <- ggplot(earlyad.hipp.plotting.data, aes(x=time, y=means, colour=Enrichment)) + geom_point() +  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi)) + geom_line()
earlyad.hipp.mmrm.plot <- earlyad.hipp.mmrm.plot + xlab("Time (Years)") + ylab("Mean Hippocampus (95% CI)") + ylim(3000, 3500)
earlyad.hipp.mmrm.plot


