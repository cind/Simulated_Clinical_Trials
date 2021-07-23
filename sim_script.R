library(simstudy)
library(ADNIMERGE)
library(dplyr)
#adas_scores <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/ADAS_ADNIGO23.csv")
#adas_scores <- adas_scores[,c("RID", "VISCODE2")]
#levels(factor(adas_scores$VISCODE))
#adas_scores$RID <- factor(adas_scores$RID)
colnames(adas_scores) <- c("RID", "VISCODE")
column.keeps          <- c("RID", "DX", "COLPROT", "ORIGPROT", "VISCODE", 
                           "EXAMDATE", "AGE", "PTGENDER", "PTEDUCAT", 
                           "APOE4", "CDRSB", "MMSE", 
                           "M", "ADAS13")
adas_demog            <- adnimerge[,column.keeps]
adas_demog$RID        <- factor(adas_demog$RID)
#adas_demog   <- merge(adas_demog, adas_scores, by=c("RID", "VISCODE"))
adas_demog            <- na.omit(adas_demog)
av45                  <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/UCBERKELEYAV45_01_14_21.csv")
av45.keeps  <- c("RID", "EXAMDATE", "SUMMARYSUVR_COMPOSITE_REFNORM")
av45        <- av45[,av45.keeps]
av45$AmyPos <- rep(NA, nrow(av45))
pos         <- which(av45$SUMMARYSUVR_COMPOSITE_REFNORM >= .78)
neg         <- which(av45$SUMMARYSUVR_COMPOSITE_REFNORM < .78)
av45["AmyPos"][pos,] <- 1
av45["AmyPos"][neg,] <- 0
av45_dict <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/UCBERKELEYAV45_DICT_01_14_21.csv")
rids.keep <- intersect(levels(factor(adas_demog$RID)), levels(factor(av45$RID)))
av45$EXAMDATE       <- as.POSIXct(av45$EXAMDATE)
adas_demog$EXAMDATE  <- as.POSIXct(adas_demog$EXAMDATE)
av45$RID      <- factor(av45$RID)
adas_demog$RID <- factor(adas_demog$RID)

### MERGING WITH AMYLOID
adas_merged   <- MergeSubjectTime(df1 = adas_demog, df2 = av45, mergecol = "RID", timecol1 = "EXAMDATE", timecol2 = "EXAMDATE")
adas_merged$years_diff <- round(abs(round(difftime(adas_merged$EXAMDATE.x, adas_merged$EXAMDATE.y, units = "days"))) / 365, 1)
adas_merged$valid      <- rep(NA, nrow(adas_merged))
adas_merged["valid"][which(adas_merged$years_diff <= 0.5), ] <- 1
adas_merged["valid"][which(adas_merged$years_diff > 0.5), ]  <- 0
adas_merged            <- adas_merged[order(adas_merged$RID, 
                                            adas_merged$EXAMDATE.x, 
                                            decreasing = FALSE), ]
adas_merged$Baseline   <- !duplicated(adas_merged$RID)

adas_merged$RID <- factor(adas_merged$RID)
adas_merged <- TimeSinceBaseline(adas_merged, timecol = "M")
adas_baseline_valid <- subset(adas_merged, valid == 1)
length(which(adas_merged$new_time==12 & adas_merged$valid == 1))
