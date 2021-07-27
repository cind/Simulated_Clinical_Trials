library(simstudy)
library(ADNIMERGE)
library(dplyr)
library(furniture)
library(lme4)
adas_scores_1  <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADASSCORES.csv")
adas_scores_23 <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADAS_ADNIGO23.csv")
adas_scores_1  <- adas_scores_1[,c("RID", "VISCODE",  "TOTAL11")]
adas_scores_23 <- adas_scores_23[,c("RID", "VISCODE2", "TOTSCORE")]
colnames(adas_scores_23) <- c("RID", "VISCODE", "TOTAL11")
adas_scores              <- rbind(adas_scores_1, adas_scores_23)
adas_scores              <- adas_scores[order(adas_scores$RID, 
                                              adas_scores$VISCODE, 
                                              decreasing = FALSE), ]
column.keeps             <- c("RID", "DX", "COLPROT", "ORIGPROT", 
                               "VISCODE", "EXAMDATE", "AGE", 
                               "PTGENDER", "PTEDUCAT", 
                               "APOE4", "CDRSB", "MMSE", 
                               "M")
adas_demog       <- adnimerge[,column.keeps]
drops1           <- which(adas_scores$VISCODE == "")
drops2           <- which(adas_scores$VISCODE == "uns1")
adas_scores      <- adas_scores[-unique(union(drops1, drops2)), ]
adas_merge_demog <- merge(adas_demog, adas_scores, by=c("RID", "VISCODE"))



adas_merge_demog <- na.omit(adas_merge_demog)
av45             <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/UCBERKELEYAV45_01_14_21.csv")
av45.keeps       <- c("RID", "EXAMDATE", "VISCODE2", "SUMMARYSUVR_COMPOSITE_REFNORM")
av45             <- av45[,av45.keeps]
colnames(av45)   <- c("RID", "EXAMDATE", "VISCODE", "SUMMARYSUVR_COMPOSITE_REFNORM")
av45$AmyPos      <- rep(NA, nrow(av45))
pos              <- which(av45$SUMMARYSUVR_COMPOSITE_REFNORM >= .78)
neg              <- which(av45$SUMMARYSUVR_COMPOSITE_REFNORM < .78)
av45["AmyPos"][pos,]  <- 1
av45["AmyPos"][neg,]  <- 0
adas_merge_demog$RID <- factor(adas_merge_demog$RID)
adas_merge_demog$EXAMDATE <- as.POSIXct(adas_merge_demog$EXAMDATE)
av45$EXAMDATE <- as.POSIXct(av45$EXAMDATE)
av45$RID             <- factor(av45$RID)
adas_merge_demog_av45 <- MergeSubjectTime(adas_merge_demog, av45, "RID", "EXAMDATE", "EXAMDATE")
adas_merge_demog_av45$timediff <- abs(difftime(adas_merge_demog_av45$EXAMDATE.x, adas_merge_demog_av45$EXAMDATE.y, units = "days"))
adas_merge_demog_av45$timediff <- adas_merge_demog_av45$timediff / 365
adas_merge_demog_av45$valid <- rep(NA, nrow(adas_merge_demog_av45))
adas_merge_demog_av45["valid"][which(adas_merge_demog_av45$timediff <= .5), ] <- 1
adas_merge_demog_av45["valid"][which(adas_merge_demog_av45$timediff > .5), ] <- 0
adas_merge_demog_av45 <- adas_merge_demog_av45[order(adas_merge_demog_av45$RID,
                                                     adas_merge_demog_av45$EXAMDATE.x,
                                                     decreasing = FALSE), ]

adas_merge_demog_av45     <- na.omit(adas_merge_demog_av45)
adas_merge_demog_av45$RID <- factor(adas_merge_demog_av45$RID)
adas_merge_demog_av45$M   <- as.numeric(adas_merge_demog_av45$M)
adas_merge_demog_av45     <- TimeSinceBaseline(adas_merge_demog_av45, timecol = "M")

rows.t1 <- which(adas_merge_demog_av45$new_time==0 & adas_merge_demog_av45$valid==1)
rows.t2 <- which(adas_merge_demog_av45$new_time==6)
rows.t3 <- which(adas_merge_demog_av45$new_time==12)
rows.t4 <- which(adas_merge_demog_av45$new_time==18)
rows.t5       <- which(adas_merge_demog_av45$new_time==24)
adas_keep     <- adas_merge_demog_av45[c(rows.t1, rows.t2, rows.t3, rows.t4, rows.t5), ]
adas_keep$RID <- factor(adas_keep$RID)
adas_keep     <-  adas_keep[order(adas_keep$RID,
                                          adas_keep$new_time,
                                          decreasing = FALSE), ]


#Baseline cohort

adas_baseline           <- subset(adas_keep, new_time == 0)
adas_baseline$DX        <- factor(adas_baseline$DX)
adas_baseline$APOE4     <- factor(adas_baseline$APOE4)
adas_baseline$PTGENDER  <- factor(adas_baseline$PTGENDER)
adas_baseline$COLPROT   <- factor(adas_baseline$COLPROT)
adas_baseline$AmyPos    <- factor(adas_baseline$AmyPos)
adas_baseline$new_time  <- factor(adas_baseline$new_time)

adas_desc <- adas_baseline[,c("DX", "APOE4", "PTGENDER", "AGE", "TOTAL11", 
                              "MMSE", "CDRSB", "PTEDUCAT", "SUMMARYSUVR_COMPOSITE_REFNORM",
                              "AmyPos")]
desc.table <- table1(adas_desc, splitby = ~ DX)$Table1
desc.table <- desc.table[c(1, 6:nrow(desc.table)), ]




#### model fitting
model1.out <- lmer(TOTAL11 ~  -1 + factor(new_time) + SUMMARYSUVR_COMPOSITE_REFNORM*factor(AmyPos) + (1|RID), data = adas_keep)
summary(model1.out)

vcov(model1.out)
