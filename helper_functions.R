#### LIBRARIES ####

library(survey)
library(zoo)
library(ADNIMERGE)
library(plyr)
library(dplyr)
library(furniture)
library(lme4)
library(nlme)
library(simr)
library(stringr)
library(matrixStats)
library(lmerTest)
library(MASS)
library(splitstackshape)
library(purrr)
library(ggplot2)
library(reshape2)
library(sjPlot)
library(pbkrtest)
library(foreach)
library(parallel)
library(doParallel)
`%notin%` <- Negate(`%in%`)

MergeSubjectTime <- function(df1, df2, mergecol, timecol1, timecol2) {
  
  
  m_data <- lapply(intersect(df1[[mergecol]], df2[[mergecol]]),function(id) {
    d1   <- subset(df1,df1[[mergecol]] == id)
    d2   <- subset(df2,df2[[mergecol]] == id)
    
    d1$indices <- sapply(d1[[timecol1]],function(d) which.min(abs(d2[[timecol2]] - d)))
    d2$indices <- 1:nrow(d2)
    base::merge(d1, d2, by = c(mergecol, 'indices'), all.x=FALSE, all.y=FALSE)
    
  })
  mergeddata <- do.call(rbind, m_data)
  mergeddata$indices <- NULL
  return(mergeddata)
}

TimeSinceBaseline <- function(data, timecol) {
  subjlist <- split(data, data$RID)
  return.list <- list()
  for(i in 1:length(subjlist)) {
    subj <- subjlist[[i]]
    if(!all(is.na(subj[[timecol]]))) {
      val.rows   <- which(!is.na(subj[[timecol]]))
      min.index <- min(val.rows)
      time.min  <- min(subj[timecol][val.rows,])
      subj$new_time <- subj[[timecol]] - time.min
      subj <- subj[min.index:nrow(subj), ]
      return.list[[i]] <- subj
    } 
  }
  return.list <- do.call(rbind, return.list)
  return(return.list)
}


Baseline <- function(data) {
  subjlist <- split(data, data$RID)
  return.list <- list()
  
  for(i in 1:length(subjlist)) {
    subj <- subjlist[[i]]
    subj <- subj[order(subj$M_vis, decreasing = FALSE),]
    subj$new_time <- subj$M_vis - min(subj$M_vis)
    return.list[[i]] <- subj
  }
  return.list <- do.call(rbind, return.list)
  return(return.list)
}

makefactor <- function(data) {
  data$RID <- factor(data$RID)
  return(data)
}

TimeSinceBaselineValidAmy <- function(data, timecol) {
  subjlist <- split(data, data$RID)
  return.list <- list()
  for(i in 1:length(subjlist)) {
    subj <- subjlist[[i]]
    if(!all(is.na(subj$AmyPos_full))) {
      val.amy   <- which(!is.na(subj$AmyPos_full))
      min.index <- min(val.amy)
      time.min  <- min(subj[timecol][val.amy,])
      subj$new_time <- subj[[timecol]] - time.min
      subj <- subj[min.index:nrow(subj), ]
      return.list[[i]] <- subj
    } 
  }
  return.list <- do.call(rbind, return.list)
  return(return.list)
}



SetNeuroData <- function(data) {
  data$TAU_pos_path <- data$Amy_pos_path <- data$TDP_pos_path <- data$Lewy_pos_path <- data$CAA_path <- rep(NA, nrow(data))
  data["TAU_pos_path"][which(data$NPBRAAK %in% c(3:7)), ] <- 1
  data["TAU_pos_path"][which(data$NPBRAAK %in% c(1,2)), ] <- 0
  data["TAU_pos_path"][which(data$NPBRAAK %in% c(8,9)), ] <- NA
  
  data["Amy_pos_path"][which(data$NPNEUR  %in% c(2,3)), ] <- 1
  data["Amy_pos_path"][which(data$NPNEUR  %in% c(0,1)), ] <- 0
  data["Amy_pos_path"][which(data$NPNEUR  %in% c(8,9)), ] <- NA
  
  data["TDP_pos_path"][which(complete.cases(data[,c("NPTDPA","NPTDPB","NPTDPC","NPTDPD","NPTDPE")])),]     <- 0
  data["TDP_pos_path"][which(data$NPTDPB == 1 | data$NPTDPC == 1 | data$NPTDPD == 1 | data$NPTDPE == 1), ] <- 1
  data["Lewy_pos_path"][which(data$NPLBOD == 0), ] <- 0
  data["Lewy_pos_path"][which(data$NPLBOD %in% c(1,2,3,4,5)), ]  <- 1
  data["Lewy_pos_path"][which(data$NPLBOD %in% c(8,9)), ]        <- NA
  data["CAA_path"][which(!is.na(data$NPAMY)), ]      <- 0
  data["CAA_path"][which(data$NPAMY %in% c(2,3)),]   <- 1
  return(data)
}



CreateBaselineVar <- function(data, timecol, baselinecol) {
  data_split <- split(data, data$RID)
  newname <- paste(baselinecol, "_bl", sep="")
  comb.list <- list()
  for(i in 1:length(data_split)) {
    subj <- data_split[[i]]
    subj[newname] <- NA
    if(!all(is.na(subj[[baselinecol]]))) {
      val <- subj[baselinecol][min(which(!is.na(subj[[baselinecol]]))),]
      subj[newname][min(which(!is.na(subj[[baselinecol]]))): nrow(subj),] <- rep(val, 
                                                                                 length(min(which(!is.na(subj[[baselinecol]]))): nrow(subj)))
    }
    comb.list[[i]] <- subj
  }
  comb.list <- do.call(rbind, comb.list)
  return(comb.list)
}





SampleJointDistribution <- function(data, sampling.list, n) {
  data.list    <- list()
  sample.vec   <- sampling.list[[1]]
  prob.level.1 <- sampling.list[[2]]
  prob.level.2 <- sampling.list[[3]]
  split.data   <- split(data, data[[sample.vec]])
  df1      <- split.data[[1]]
  df2      <- split.data[[2]]
  df1$prob <- rep(prob.level.1, nrow(df1))
  df2$prob <- rep(prob.level.2, nrow(df2))
  df.final <- bind_rows(df1, df2)
  test.sample <- sample(rownames(df.final), size = n, replace = FALSE, prob = df.final$prob)
  df.final <- df.final[test.sample, ]
  return(df.final)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



PullLongData <- function(baseline.data, long.data) { 
  data.id <- baseline.data[["RID"]]
  longdata <- subset(long.data, RID %in% data.id)
  longdata$RID <- factor(longdata$RID)
  return(longdata)
}





PlotObsData <- function(data, formula.fixed, ylab) {
  model.fit <- lme(fixed = as.formula(formula.fixed), 
                   random = ~1|RID,
                   data=data)
  sumstat <- summarySE(data = data,
                       measurevar = word(formula.fixed, 1), 
                       groupvars = c("new_time"))
  sumstat$min.se <- sumstat[["mean"]] - (1.96*sumstat$se)
  sumstat$max.se <- sumstat[["mean"]] + (1.96*sumstat$se)
  sumstat$new_time <- as.numeric(as.character(sumstat$new_time))
  
  gplot <- ggplot(sumstat, aes_string(x="new_time", y="mean")) + 
    geom_errorbar(aes_string(ymin="min.se", ymax="max.se", width=.1)) + geom_line(group=1) + geom_point() + scale_x_discrete(name = "Weeks", limits=c(0, 6, 12, 24)) + ylab(ylab)
  return.list <- list("model" = model.fit,
                      "summary_stats" = sumstat,
                      "plot" = gplot)
  return(return.list)
  
}


SampleSizeSimulation <- function(sim.data, formula, fcompare_str, efficacy=.5, breaks, yaxislab_dpm, ptau=NULL) {
  if(!is.null(ptau)) {
    missing.vals <- unique(sim.data["RID"][which(is.na(sim.data$PTAU_bl)),])
    sim.data <- subset(sim.data, RID %notin% missing.vals)
    sim.data$RID <- factor(sim.data$RID)
    rownames(sim.data) <- 1:nrow(sim.data)
  }
  all.rids   <- levels(sim.data[["RID"]])
  nrids      <- nlevels(sim.data[["RID"]])
  n.sample   <- round(nrids/2)
  treat.rids <- sample(all.rids, n.sample, prob = rep(.5, length(all.rids)), replace = FALSE)
  base.rids  <- subset(all.rids, all.rids %notin% treat.rids)
  sim.data$treat <- rep(NA, nrow(sim.data))
  treat.rows     <- which(sim.data$RID %in% treat.rids)
  base.rows      <- which(sim.data$RID %in% base.rids)
  sim.data["treat"][treat.rows,] <- 1
  sim.data["treat"][base.rows,]  <- 0
  sim.data$treat <- factor(sim.data$treat)
  if(!is.null(ptau)) {
    sim.model       <- lmerTest::lmer(as.formula(ptau), data = sim.data)
    fixed           <- fixef(sim.model)
    fixed["treat1"]          <- 0
    fixed["new_time:treat1"] <- 0
    fixed["treat1:PTAU_bl"]  <- 0
    fixed["new_time:treat1:PTAU_bl"] <- (-1*(fixed["new_time:PTAU_bl"] /2))
    sum_model <- summary(sim.model)
    sig    <- summary(sim.model)$sigma
    varcor <- as.numeric(summary(sim.model)$varcor[[1]])
    
    constr.lme                                          <- makeLmer(as.formula(ptau), 
                                                                    fixef    = fixed, 
                                                                    VarCorr  = list(varcor), 
                                                                    sigma    = sig, 
                                                                    data     = sim.data)
    sim_ext_rid     <- extend(constr.lme, along="RID", n=max(breaks))
    return(constr.lme)
  } else {
    sim.model <- lmerTest::lmer(as.formula(formula), data = sim.data)
    fixed <- fixef(sim.model) 
    fixed["treat1"] <- 0
    fixed["treat1:new_time"] <-  (- (fixed["new_time"] /2))
    rand  <-  as.numeric(summary(sim.model)$varcor[[1]])
    sig   <- summary(sim.model)$sigma
    constr.lme <-  makeLmer(as.formula(formula), fixef = fixed, VarCorr = rand, sigma = sig, data = sim.data)
  }
  #build contrast progression plot
  treat0 <- sim.data[base.rows, ]
  treat1 <- treat0
  treat1$treat <- 1
  plot.data <- rbind(treat0, treat1)
  plot.data$prediction <- predict(constr.lme, plot.data)
  plot.disease.contr <- ggplot(plot.data, aes(x=new_time, y=prediction, colour=treat)) + geom_smooth(method = "lm")
  plot.disease.contr <- plot.disease.contr + xlab("Time (Months)") + ylab(yaxislab_dpm) + labs(colour="Treatment")
  sim_ext_rid        <- extend(constr.lme, along="RID", n=max(breaks))
  p_curve_treat_sim <- powerCurve(sim_ext_rid, test=fcompare(as.formula(fcompare_str)), along="RID", breaks=breaks)
  summ_sim<- summary(p_curve_treat_sim)
  summ_sim$mean  <- summ_sim$mean*100
  summ_sim$lower <- summ_sim$lower*100
  summ_sim$upper <- summ_sim$upper*100
  gplot.sim <- ggplot(data = summ_sim, aes(x=nlevels, y=mean)) +geom_errorbar(aes(ymin=lower, ymax=upper)) + geom_line() + geom_point() + scale_x_discrete(limits=breaks)
  gplot.sim <- gplot.sim + xlab("Sample Size per Arm") + ylab("Statistical Power (%)")
  return.list <- list("model" = constr.lme,
                      "power_curve_output" = p_curve_treat_sim,
                      "summary_stats" = summ_sim,
                      "power_curve_plot" = gplot.sim,
                      "disease_progression_plot" = plot.disease.contr)
  return(return.list)
}



GroupDiseaseTraj <- function(sim.list, yaxislab_dpm) {
  firstelement <- sim.list[[1]]
  model      <- firstelement$model
  names1     <- c(names(sim.list)[1])
  plot.data  <- firstelement$disease_progression_plot$data
  total.rids <- nlevels(factor(plot.data$RID))
  firstfail  <- paste(total.rids, " Subjects ","(Reference Group) \n", sep="")
  names1 <- paste(names1, firstfail, sep="")
  plot.data$group <- rep(names1, nrow(plot.data))
  plot.data$predicted <- predict(model, plot.data)
  screen.fail <- list()
  for(i in 2:length(sim.list)) {
    sublist <- sim.list[[i]]
    nextframe <- sublist$disease_progression_plot$data
    nextframe.model <- sublist$model
    groupname <- names(sim.list)[i]
    nrids <- nlevels(factor(nextframe$RID))
    failrate <- paste(nrids, " Subjects ", "(", 
                      round(((total.rids - nrids) / total.rids) * 100, 2), 
                      "% failure rate) \n", sep="")
    groupname <- paste(groupname, failrate, sep="")
    nextframe$group <- rep(groupname, nrow(nextframe))
    nextframe$predicted <- predict(nextframe.model, nextframe)
    plot.data <- rbind.fill(plot.data, nextframe)
    screen.fail[[i]] <- paste(nrids, " Subjects ", "(", round(((total.rids - nrids) / total.rids) * 100, 2), "% failure rate)", sep="")
    
  }
  names(screen.fail) <- cat(names(sim.list))
  plot.data$Treatment <- plot.data$treat
  plot.data$Enrichment <- plot.data$group
  
  gplot <- ggplot(plot.data, aes(y=predicted, x=new_time, colour=Enrichment, linetype = Treatment)) + geom_smooth(method = "lm", aes(fill=Enrichment), alpha=.1) 
  gplot <- gplot  + xlab("Time (Years)") + ylab(yaxislab_dpm)
  return(list("plot"=gplot, "screenfail"= screen.fail))
}





CombineSimPlots <- function(power.list, limits) {
  nonenrich <- power.list[[1]] 
  name1 <- names(power.list)[1]
  sumstats_non <- nonenrich$summary_stats
  sumstats_non$Group <- rep(name1, nrow(sumstats_non))
  plot.df <- sumstats_non
  for(i in 2:length(power.list)){
    enriched <- power.list[[i]]
    name <- names(power.list)[i]
    enriched <- enriched$summary_stats
    enriched$Group <- rep(name, nrow(enriched))
    plot.df <- rbind(plot.df, enriched)
  }
  gplot.sim <- ggplot(data = plot.df, aes(x=nlevels, y=mean, colour=Group)) +geom_errorbar(aes(ymin=lower, ymax=upper)) + geom_line() + geom_point() + scale_x_discrete(limits=limits)
  gplot.sim <- gplot.sim + xlab("Sample Size per Arm") + ylab("Statistical Power (%)")
  return(list("plot" = gplot.sim, "fullstats" = plot.df))
} 



QuickAdjust <- function(data) {
  data$RID <- factor(data$RID)
  data <- data[order(data$RID, data$M, decreasing = FALSE),]
  data <- TimeSinceBaseline(data, "M")
  data <- subset(data, new_time <= 24)
  data$new_time <- (data$new_time / 12)
  return(data)
}

Sub24 <- function(data) {
  data <- subset(data, new_time <= 24)
  data$RID <- factor(data$RID)
  return(data)
}

BuildSignificanceTable <- function(model) {
  summodel<-summary(model)
  coefdf <- summodel$coefficients
  coefdf <- coefdf[,c(1,2,4,5)]
  coefdf <- as.data.frame(coefdf)
  coefdf$Significant <- rep("", nrow(coefdf))
  for(i in 1:nrow(coefdf)) {
    if(coefdf[4][i,] > .05) {
    } else if(coefdf[4][i,] <= .05 & coefdf[4][i,] > 0.01) {
      coefdf["Significant"][i,] <- "*"
    } else if(coefdf[4][i,] <= .01 & coefdf[4][i,] > 0.001) {
      coefdf["Significant"][i,] <- "**"
    } else if(coefdf[4][i,] <= .001) {
      coefdf["Significant"][i,] <- "***"
    }
  }  
  return(coefdf)
  
}


BuildESTable <- function(data) {
  rows.int <- rownames(data[3:10,])
  rows.time <- rownames(data[11:18,])
  newframe <- data.frame(matrix(ncol = 3, nrow = 8))
  newframe$X1 <- unname(unlist(list.names[rows.int]))
  newframe$X2 <- paste(round(data["Estimate"][3:10,], 4), " (", round(data["Std. Error"][3:10,], 4), "), ", round(data["Pr(>|t|)"][3:10,], 4), data["Significant"][3:10,] , sep="")
  newframe$X3 <- paste(round(data["Estimate"][11:18,], 4), " (", round(data["Std. Error"][11:18,], 4), "), ", round(data["Pr(>|t|)"][11:18,], 4), data["Significant"][11:18,] , sep="")
  colnames(newframe) <- c("Variable", "Cognitive Level", "Cognitive Decline")
  
  return(newframe)
}



BuildNeuroPatDistr <- function(model) {
  modmat      <- as.data.frame(model.matrix(model))
  coef.modmat <- as.matrix(coef(model)[[1]])
  modmat      <- subset(modmat, new_time==0)
  modmat      <- modmat[,c("fulllewy1", "fullcaa1", 
                           "fulltdp431", "Amy_pos_path1", 
                           "TAU_pos_path1")]
  colnames(modmat) <- c("Lewy Body", "CAA", "TDP43", "Amyloid", "Tau")
  coef.modmat <- coef.modmat[,c("new_time:fulllewy1", "new_time:fullcaa1", 
                                "new_time:fulltdp431", "new_time:Amy_pos_path1", 
                                "new_time:TAU_pos_path1")]
  
  el.wise.mult <- as.matrix(coef.modmat * modmat)
  fullsum      <- rowSums2(el.wise.mult)
  checkfinal   <- el.wise.mult / fullsum
  checkfinal[is.nan(checkfinal)] <- 0 
  checkfinal   <- as.data.frame(checkfinal)
  counts       <- apply(modmat, 2, count)
  means        <- apply(checkfinal, 2, mean) * 100
  sds          <- apply(checkfinal, 2, sd) * 100
  counts <- paste(counts, " (", (counts / nrow(modmat)) * 100, " %", ")", sep="")
  means  <- paste(round(means, 3), " %", sep="")
  sds  <- paste(round(sds, 3), " %", sep="")
  returnframe <- data.frame("Neuropathology" =  colnames(modmat),
                            "Count"          =  counts,
                            "Mean"           =  means,
                            "SD"             =  sds)
  return(returnframe)
  
}


BuildNeuroPatDistrImp <- function(model) {
  modmat      <- as.data.frame(model.matrix(model))
  coef.modmat <- as.matrix(coef(model)[[1]])
  modmat      <- subset(modmat, new_time==0)
  modmat      <- modmat[,c("fulllewy1", "fullcaa1", 
                           "fulltdp431", "AmyPos_bl1", 
                           "ptau_pos_bl1")]
  colnames(modmat) <- c("Lewy Body", "CAA", "TDP43", "Amyloid", "Tau")
  coef.modmat <- coef.modmat[,c("new_time:fulllewy1", "new_time:fullcaa1", 
                                "new_time:fulltdp431", "new_time:AmyPos_bl1", 
                                "new_time:ptau_pos_bl1")]
  
  el.wise.mult <- as.matrix(coef.modmat * modmat)
  fullsum      <- rowSums2(el.wise.mult)
  checkfinal   <- el.wise.mult / fullsum
  checkfinal[is.nan(checkfinal)] <- 0 
  checkfinal   <- as.data.frame(checkfinal)
  counts       <- apply(modmat, 2, count)
  means        <- apply(checkfinal, 2, mean) * 100
  sds          <- apply(checkfinal, 2, sd) * 100
  counts <- paste(counts, " (", round((counts / nrow(modmat)), 3) * 100, " %", ")", sep="")
  means  <- paste(round(means, 3), " %", sep="")
  sds  <- paste(round(sds, 3), " %", sep="")
  returnframe <- data.frame("Neuropathology" =  colnames(modmat),
                            "Count"          =  counts,
                            "Mean"           =  means,
                            "SD"             =  sds)
  return(returnframe)
  
}



RandomizeTreatment <- function(data) {
  full.data    <- data
  sim.data     <- data
  sim.data     <- subset(sim.data, new_time==0)
  sim.data$RID <- factor(sim.data$RID)
  all.rids   <- levels(sim.data[["RID"]])
  nrids      <- nlevels(sim.data[["RID"]])
  n.sample   <- round(nrids/2)
  treat.rids <- sample(all.rids, n.sample, prob = rep(.5, length(all.rids)), replace = FALSE)
  base.rids  <- subset(all.rids, all.rids %notin% treat.rids)
  full.data$treat <- rep(NA, nrow(full.data))
  treat.rows     <- which(full.data$RID %in% treat.rids)
  base.rows      <- which(full.data$RID %in% base.rids)
  full.data["treat"][treat.rows,] <- 1
  full.data["treat"][base.rows,]  <- 0
  full.data$treat <- factor(full.data$treat)
  sim.data     <- subset(full.data, new_time==0)
  sim.data$RID <- factor(sim.data$RID)  
  data.contr <- subset(sim.data, treat==0)
  data.treat <- subset(sim.data, treat==1)
  lewy <- c(length(which(data.contr$fulllewy  ==0)), length(which(data.treat$fulllewy==0)))
  tdp  <- c(length(which(data.contr$fulltdp43 ==0)), length(which(data.treat$fulltdp43==0)))
  caa  <- c(length(which(data.contr$fullcaa   ==0)), length(which(data.treat$fullcaa==0)))
  n.compare <- rep(nrow(sim.data), 2)
  #proplist <- list("lewy" = prop.test(x=lewy, n=n.compare),
  #                 "tdp" = prop.test(x=tdp, n=n.compare),
  #                 "caa" = prop.test(x=caa, n=n.compare))
  
  return(full.data)
}



SampleSizeSimulation2 <- function(sim.data, model, formula, compare_str, breaks, yaxislab_dpm, return_dpm=FALSE) {
  placebo <- sim.data
  placebo <- placebo[order(placebo$RID, placebo$new_time, decreasing = FALSE),]
  placebo$treat <- rep(0, nrow(placebo))
  placebo$treat <- factor(placebo$treat)
  placebo$prediction <- predict(model, placebo)
  treatgroup <- sim.data
  treatgroup <- treatgroup[order(treatgroup$RID, treatgroup$new_time, decreasing = FALSE),]
  treatgroup$treat <- rep(1, nrow(treatgroup))
  treatgroup$treat <- factor(treatgroup$treat)
  treatgroup$prediction <- predict(model, treatgroup)
  plot.data <- rbind(placebo, treatgroup)
  plot.data$treat <- factor(plot.data$treat)
  plot.data$prediction <- predict(model, plot.data)
  plot.disease.contr <- ggplot(plot.data, aes(x=new_time, y=prediction, colour=treat)) + geom_smooth(method = "lm")
  plot.disease.contr <- plot.disease.contr + scale_x_discrete(limits=c(0, .5, 1, 1.5, 2)) + ylab(yaxislab_dpm) + labs(colour="Treatment") + xlim(0, 2) + xlab("Time (Years)")
  if(return_dpm) {
    return(list("model" = model,
                "disease_progression_plot"=plot.disease.contr))
  }
  
  sim_ext_rid        <- simr::extend(model, along="RID", n=max(breaks))
  p_curve_treat_sim  <- powerCurve(sim_ext_rid, test=compare(as.formula(compare_str)), along="RID", breaks=breaks)
  summ_sim           <- summary(p_curve_treat_sim)
  summ_sim$mean      <- summ_sim$mean*100
  summ_sim$lower     <- summ_sim$lower*100
  summ_sim$upper     <- summ_sim$upper*100
  gplot.sim          <- ggplot(data = summ_sim, aes(x=nlevels, y=mean)) +geom_errorbar(aes(ymin=lower, ymax=upper)) + geom_line() + geom_point() + scale_x_discrete(limits=breaks)
  gplot.sim          <- gplot.sim + xlab("Sample Size per Arm") + ylab("Statistical Power (%)")
  return.list        <- list("model"                    = model,
                             "power_curve_output"       = p_curve_treat_sim,
                             "summary_stats"            = summ_sim,
                             "power_curve_plot"         = gplot.sim
                             #"disease_progression_plot" = plot.disease.contr
  )
  return(return.list)
}


BuildSimulationModel <- function(list, formula.model, data, treatment.effect, es) {
  model            <- list[[1]]
  adjusted.decline <- list[[2]]
  fixd             <- fixef(model)
  fixd["treat1"]   <- 0
  if(treatment.effect =="controlled") {
    fixd["new_time:treat1"] <- ((adjusted.decline * .5) * -1)
  } else if(is.numeric(treatment.effect)) {
    fixd["new_time:treat1"] <- ((treatment.effect * .5) * -1)
  } else {
    fixd["new_time:treat1"] <- ((fixd["new_time"] * .5) * -1)
  }
  fixd["fulllewy1"]           <- es[1]
  fixd["fullcaa1"]            <- es[2]
  fixd["fulltdp431"]          <- es[3]
  fixd["new_time:fulllewy1"]  <- es[4]
  fixd["new_time:fullcaa1"]   <- es[5]
  fixd["new_time:fulltdp431"] <- es[6]
  
  
  fixd                        <- fixd[c( "(Intercept)", "new_time", "treat1", "PTEDUCAT_bl", "AGE_bl", "PTGENDER_blMale", 
                                         "fulllewy1", "fullcaa1", "fulltdp431", "new_time:treat1", "new_time:fulllewy1", 
                                         "new_time:fullcaa1", "new_time:fulltdp431")]
  sigma.mod        <- summary(model)$sigma
  varcor.mod       <- as.numeric(summary(model)$varcor[[1]])
  constr.lme       <- makeLmer(formula = as.formula(formula.model), fixef = fixd, VarCorr=list(varcor.mod), sigma = sigma.mod, data = data)
  return(constr.lme)
}

BuildSimulationModelNoPath <- function(list, formula.model, data, treatment.effect) {
  model            <- list
  #adjusted.decline <- list[[2]]
  fixd             <- fixef(model)
  fixd["treat1"]   <- 0
  if(treatment.effect =="controlled") {
    fixd["new_time:treat1"] <- ((adjusted.decline * .5) * -1)
  }  else if(is.numeric(treatment.effect)) {
    fixd["new_time:treat1"] <- ((treatment.effect * .5) * -1)
  } else {
    fixd["new_time:treat1"] <- ((fixd["new_time"] * .5) * -1)
  }
  if(length(unique(data$CDGLOBAL_bl)) == 1) {
  fixd                            <- fixd[c( "(Intercept)", "new_time","treat1", "PTEDUCAT_bl", "AGE_bl", "PTGENDERMale", "MMSE_bl",  "new_time:treat1")]
  } else {
     fixd                        <- fixd[c( "(Intercept)", "new_time","treat1", "PTEDUCAT_bl", "AGE_bl", "PTGENDERMale", "MMSE_bl", "CDGLOBAL_bl1",  "new_time:treat1")]
   }
  sigma.mod        <- summary(model)$sigma
  varcor.mod       <- VarCorr(model)
  constr.lme       <- makeLmer(formula = as.formula(formula.model), fixef = fixd, VarCorr=varcor.mod, sigma = sigma.mod, data = data)
  return(constr.lme)
}

BuildSimulationModelNoPathPlotting <- function(list, formula.model, data, treatment.effect) {
  model            <- list
  #adjusted.decline <- list[[2]]
  fixd             <- fixef(model)
  fixd["treat1"]   <- 0
  if(treatment.effect =="controlled") {
    fixd["new_time:treat1"] <- ((adjusted.decline * .5) * -1)
  }  else if(is.numeric(treatment.effect)) {
    fixd["new_time"] <- (treatment.effect)
  } else {
    fixd["new_time:treat1"] <- ((fixd["new_time"] * .5) * -1)
  }
  if(length(unique(data$CDGLOBAL_bl)) == 1) {
    fixd                            <- fixd[c( "(Intercept)", "new_time","treat1", "PTEDUCAT_bl", "AGE_bl", "PTGENDERMale", "MMSE_bl",  "new_time:treat1")]
  } else {
    fixd                        <- fixd[c( "(Intercept)", "new_time","treat1", "PTEDUCAT_bl", "AGE_bl", "PTGENDERMale", "MMSE_bl", "CDGLOBAL_bl1",  "new_time:treat1")]
  }
  sigma.mod        <- summary(model)$sigma
  varcor.mod       <- VarCorr(model)
  constr.lme       <- makeLmer(formula = as.formula(formula.model), fixef = fixd, VarCorr=varcor.mod, sigma = sigma.mod, data = data)
  return(constr.lme)
}

BuildSignificanceTable <- function(model) {
  list.names <- list("AGE_bl" = "Age (Baseline)", 
                     "PTGENDER_blMale" = "Gender (Male)",
                     "PTEDUCAT_bl" ="Education Years (Baseline)",
                     "fulllewy1" = "Lewy body",
                     "fullcaa1"  = "CAA",
                     "fulltdp431" = "TDP43",
                     "AmyPos_bl1" = "Amyloid",
                     "ptau_pos_bl1" = "Tau",
                     "Amy_pos_path1" = "Amyloid",
                     "TAU_pos_path1" = "Tau",
                     "new_time"      = "Mean Rate",
                     "new_time:fulllewy1" = "Lewy body (Rate)",
                     "new_time:fullcaa1"  = "CAA (Rate)",
                     "new_time:fulltdp431" = "TDP43 (Rate)",
                     "(Intercept)" = "Intercept")
  summodel<-summary(model)
  coefdf <- summodel$coefficients
  coefdf <- coefdf[,c(1,2,4,5)]
  coefdf <- as.data.frame(coefdf)
  coefdf$Significant <- rep("", nrow(coefdf))
  for(i in 1:nrow(coefdf)) {
    if(coefdf[4][i,] > .05) {
    } else if(coefdf[4][i,] <= .05 & coefdf[4][i,] > 0.01) {
      coefdf["Significant"][i,] <- "*"
    } else if(coefdf[4][i,] <= .01 & coefdf[4][i,] > 0.001) {
      coefdf["Significant"][i,] <- "**"
    } else if(coefdf[4][i,] <= .001) {
      coefdf["Significant"][i,] <- "***"
    }
  }
  rownames(coefdf) <- list.names[rownames(coefdf)]
  return(coefdf)
}



KeepCols <- function(data, cols) {
  return(data[,cols])
}

BuildDescTable <- function(data, columns, split=FALSE) {
  data <- data[,columns]
  result <- table1(data[,columns])$Table1
  return(result)
}



ZscoreAdj <- function(data, col_names, control.data) {
  colsdata <- colnames(data)
  zvarcols <- paste(col_names, "_zscore", sep="")
  for(i in col_names) {
    bline <- subset(control.data, new_time==0)
    zvar.mean <- mean(bline[[i]])
    print(zvar.mean)
    zvar.sd <- sd(bline[[i]])
    print(zvar.sd)
    zvar <- (data[[i]] - zvar.mean) / zvar.sd
    data <- cbind(data, zvar)
  }
  colnames(data) <- c(colsdata, zvarcols)
  return(data)
}


RandomizeTreatment2 <- function(data, longdata, no.prop=NULL) {
  `%notin%`  <- Negate(`%in%`)
  stratifydf <- data
  stratifydf$stratvar   <- interaction(stratifydf$CAAPos, stratifydf$LewyPos, stratifydf$TDP43Pos, stratifydf$PTGENDER, stratifydf$AGE_bl_strat,
                                     stratifydf$PTEDUCAT_bl_strat,stratifydf$MMSE_bl_strat, stratifydf$CDGLOBAL_bl, stratifydf$TauPos_bl)
  stratifydf$stratvar   <- factor(stratifydf$stratvar)
  stratifieddata        <- stratified(stratifydf, "stratvar", size = (.5), bothSets = FALSE)
  rownames(stratifieddata) <- 1:nrow(stratifieddata)
  half <- round(nrow(data) / 2)
  diff <- half - nrow(stratifieddata)
  if(diff > 0) {
    disjoint <- subset(stratifydf, RID %notin% stratifieddata$RID)
    add.rows <- sample_n(disjoint, diff)
    stratifieddata <- rbind(stratifieddata, add.rows)
    placebo        <-  subset(stratifydf, RID %notin% stratifieddata$RID)
  } else if(diff < 0)  {
    drop.rows <- sample(1:nrow(stratifieddata), abs(diff))
    stratifieddata <- stratifieddata[-drop.rows,]
    placebo        <-  subset(stratifydf, RID %notin% stratifieddata$RID)
  } else {
    placebo        <-  subset(stratifydf, RID %notin% stratifieddata$RID)
  }
  stratifieddata <- list(stratifieddata, placebo)
  names(stratifieddata) <- c("Treatment", "Placebo")
  get.props             <- Map(CalcProportionPos, stratifieddata)
  treatmentrids         <- stratifieddata[[1]]
  treatmentrids         <- unique(treatmentrids$RID)
  controlrids           <- data["RID"][which(data$RID %notin% treatmentrids),]
  treatmentrows         <- subset(longdata, RID %in% treatmentrids)
  treatmentrows$treat   <- rep(1, nrow(treatmentrows))
  controlrows           <- subset(longdata, RID %in% controlrids)
  controlrows$treat     <- rep(0, nrow(controlrows))
  returndata            <- rbind(treatmentrows, controlrows)
  returndata$treat      <- factor(returndata$treat)
  if(!is.null(no.prop)) {
    return(returndata)
  } else {
  return(list("data"= returndata, "props" = get.props))
  }
}



CalcProportionPos <- function(data, keepcols = c("CAAPos", "LewyPos", "TDP43Pos", "PTGENDER", "AGE_bl_strat",
                                                 "PTEDUCAT_bl_strat","MMSE_bl_strat", "CDGLOBAL_bl", "TauPos_bl")) {
  data <- as.data.frame(data)
  data <- data[,keepcols]
  prop <- Map(function(x) {table(x)/sum(table(x))}, data)
  prop <- lapply(prop, `[`, 1)
  prop <- unlist(prop)
  return(prop)
}



RandomizeTreatment2NoPath <- function(data, stratcolumns, longdata) {
  stratifydf          <- data[,c("RID", stratcolumns)]
  stratifydf$stratvar <- interaction(stratifydf$PTGENDER, stratifydf$AGE_bl_strat)
  stratifieddata      <- stratified(stratifydf, "stratvar", size = c(.5), bothSets = TRUE)
  treatmentrids       <- unique(stratifieddata$RID)
  controlrids         <- data["RID"][which(data$RID %notin% treatmentrids),]
  treatmentrows       <- subset(longdata, RID %in% treatmentrids)
  treatmentrows$treat <- rep(1, nrow(treatmentrows))
  controlrows         <- subset(longdata, RID %in% controlrids)
  controlrows$treat   <- rep(0, nrow(controlrows))
  returndata          <- rbind(treatmentrows, controlrows)
  returndata$treat    <- factor(returndata$treat)
  return(returndata)
}


StratifyContVar <- function(data, stratcols) {
  for(i in stratcols) {
    var <- data[,i]
    qt <- quantile(var)
    groupingvar <- cut(var, breaks = c(-Inf, unname(qt[3]), Inf), labels=c(0, 1))
    stratname <- paste(i, "_strat", sep="")
    data[stratname] <- factor(as.character(groupingvar))
  }
  return(data)
}


MapLmer <- function(newdata, formula.model) {
  newdata <- newdata
  lme.fit <- lmer(as.formula(formula.model), data = newdata, REML=TRUE, control = lmerControl(optimizer ="nmkbw"))
  return(lme.fit)
}


RemoveNormalAging <- function(control.model, enriched.model, data, formula.model) {
  control.decline  <- fixef(control.model)["new_time"]
  enriched.decline <- fixef(enriched.model)["new_time"]
  adjusted.decline <- enriched.decline - control.decline
  fixd             <- fixef(enriched.model)
  #fixd["new_time"] <- adjusted.decline
  sigma.mod        <- summary(enriched.model)$sigma
  varcor.mod       <- as.numeric(summary(enriched.model)$varcor[[1]])
  constr.lme       <- makeLmer(formula = as.formula(formula.model), fixef = fixd, VarCorr=list(varcor.mod), sigma = sigma.mod, data = data)
  return(list(constr.lme, adjusted.decline))
}

ChangeNeuroFixEf <- function(list, ES) {
  model <- list[[1]]
  adjusted.decline <- list[[2]]
  fixef(model)["fulllewy1"]  <- ES[1]
  fixef(model)["fullcaa1"]   <- ES[2]
  fixef(model)["fulltdp431"] <- ES[3]
  fixef(model)["new_time:fulllewy1"]  <- ES[4]
  fixef(model)["new_time:fullcaa1"]   <- ES[5]
  fixef(model)["new_time:fulltdp431"] <- ES[6]
  return(list(model, adjusted.decline))
}


MapNames <- function(list, nameslist) {
  names(list) <- nameslist
  return(list)
}



RegressNeuro <- function(data, es, outcome) {
  newvar <- paste(outcome, "neuradj", sep = "_")
  data$fulllewy  <- as.numeric(as.character(data$fulllewy))
  data$fulltdp43 <- as.numeric(as.character(data$fulltdp43))
  data$fullcaa   <- as.numeric(as.character(data$fullcaa))
  
  data$lewyint <- es[1]
  data$caaint  <- es[2]
  data$tdpint  <- es[3]
  
  data$lewydec <- es[4]
  data$caadec  <- es[5]
  data$tdpdec  <- es[6]
  
  full_ints <- data[,c("fulllewy", "fullcaa", "fulltdp43")] * data[,c("lewyint", "caaint", "tdpint")]
  full_decs <- data[,c("fulllewy", "fullcaa", "fulltdp43")] * data[,c("new_time")] * data[,c("lewydec", "caadec", "tdpdec")]
  
  full.rem  <- full_ints + full_decs
  full.rem  <- rowSums2(as.matrix(full.rem))
  
  data$fulllewy  <- as.factor(as.character(data$fulllewy))
  data$fulltdp43 <- as.factor(as.character(data$fulltdp43))
  data$fullcaa   <- as.factor(as.character(data$fullcaa))
  
  data[[newvar]] <- data[[outcome]] - full.rem
  return(data)
}



feature.correction <- function(training.data,  data, formula, cr.feat1, feat) {
  model <- lm(formula = as.formula(formula), data = training.data)
  mean.val1 <- mean(training.data[[cr.feat1]])
  mean.val1 <- rep(mean.val1, nrow(data))
  coef.correction1 <- model$coefficients[[cr.feat1]]
  new.feat <- data[[feat]]
  new.feat <- (new.feat - (coef.correction1*data[[cr.feat1]]))
  new.feat <- new.feat  + (coef.correction1*mean.val1)
  return(new.feat)
}

MMRMTime <- function(data) {
  #data$new_time_mmrm <- NA
  #data["new_time_mmrm"][which(data$new_time==0), ] <- 1
  #data["new_time_mmrm"][which(data$new_time==.5), ] <- 2
  #data["new_time_mmrm"][which(data$new_time==1), ] <- 3
  #data["new_time_mmrm"][which(data$new_time==1.5), ] <- 4
  #data["new_time_mmrm"][which(data$new_time==2), ] <- 5
  new_time_vec <- unique(data$new_time)
  new_time_vec <- new_time_vec[order(new_time_vec, decreasing = FALSE)]
  new_time_match <- 1:length(new_time_vec)
  data$new_time_mmrm <- NA
  for(i in 1:length(new_time_match)) {
    data["new_time_mmrm"][which(data$new_time==new_time_vec[i]), ] <- new_time_match[i]
  }
  drops <- which(duplicated(data[,c("RID", "new_time_mmrm")]) == TRUE)
  data <- data[-drops,]
  data$RID <- factor(data$RID)
  return(data)
}



CalculateSampleAtPower <- function(mean.points, conf.low.points, conf.hi.points, y=80) {
  m.mean     <- (mean.points[2] - mean.points[1]) / (mean.points[4] - mean.points[3])
  conf.low.m <- (conf.low.points[2] - conf.low.points[1]) / (conf.low.points[4] - conf.low.points[3])
  conf.hi.m <- (conf.hi.points[2] - conf.hi.points[1]) / (conf.hi.points[4] - conf.hi.points[3])
  mean.shift <- (80 - mean.points[1]) / m.mean
  new.x <- mean.points[3] + mean.shift
  conf.low.m.shift <- (80 - conf.low.points[1]) / conf.low.m
  new.conf.low <- conf.low.points[3] + conf.low.m.shift
  
  conf.hi.m.shift <- (80 - conf.hi.points[1]) / conf.hi.m
  new.conf.hi <- conf.hi.points[3] + conf.hi.m.shift
  
  return(c("mean" = new.x, "ci.low" = new.conf.low, "ci.high" = new.conf.hi))
}

CalculateSampleAtPowerModel <- function(model.list) {
  init.data <- data.frame(matrix(ncol = 3))
  modelnames <- names(model.list)
  for(i in 1:length(model.list)) {
  model_unenriched <- model.list[[i]]
  model_unenriched_powerdata <- try(model_unenriched$power.data)
  model_unenriched_powerdata$ss <- as.numeric(sapply(strsplit(rownames(model_unenriched_powerdata), "_"), "[[", 2))
  lowerval <- max(which(model_unenriched_powerdata$mean < .8))
  upperval <- min(which(model_unenriched_powerdata$mean >= .8))
  xvals <- c(model_unenriched_powerdata["ss"][lowerval,], model_unenriched_powerdata["ss"][upperval,])
  samplesize <- CalculateSampleAtPower(c(model_unenriched_powerdata["mean"][lowerval,]*100, model_unenriched_powerdata["mean"][upperval,]*100, xvals),
                                       c(model_unenriched_powerdata["ci.low"][lowerval,]*100, model_unenriched_powerdata["ci.low"][upperval,]*100, xvals),
                                       c(model_unenriched_powerdata["ci.high"][lowerval,]*100, model_unenriched_powerdata["ci.high"][upperval,]*100, xvals))
  init.data <- rbind(init.data, samplesize)
  }
  init.data <- init.data[2:nrow(init.data),]
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  modelnames <- substrRight(modelnames, 2)
  modelnames <- (as.numeric(modelnames) / 10)
  init.data$td <- modelnames
  colnames(init.data) <- c("Mean", "CI_high", "CI_low", "Trial Duration (Years)")

  init.data <- init.data[,c("Trial Duration (Years)","Mean", "CI_low", "CI_high")]
  
  return(init.data)
}






Keep1YearorMore <- function(data.list) {
  data.long <- data.list[["long"]]
  data.cs   <- data.list[["cs"]]
  newlist.cs <- list()
  newlist.long <- list()
  
  names.list <- names(data.long)
  for(i in 1:length(data.long)) {
    df <- data.long[[i]]
    df.cs <-  data.cs[[i]]
    df$RID <- factor(df$RID)
    rid.drops <- names(which(unlist(Map(nrow, (split(df, df$RID)))) <= 2))
    df <- subset(df, RID %notin% rid.drops)
    df.cs <- subset(df.cs, RID %notin% rid.drops)
    df$RID <- factor(df$RID)
    df.cs$CDGLOBAL_bl <- factor(df.cs$CDGLOBAL_bl)
    df$CDGLOBAL_bl <- factor(df$CDGLOBAL_bl)
    df.cs$RID <- factor(df.cs$RID)
    newlist.cs[[i]] <- df.cs
    newlist.long[[i]] <- df
    
  }
  names(newlist.cs)   <- names.list
  names(newlist.long) <- names.list
  newlist <- list("cs" = newlist.cs,
                  "long" = newlist.long)
  return(newlist)
}

PlotSimulationLME <- function(model, data, y, ylab, rids) {
  newdat <- expand.grid(new_time=unique(data$new_time),
                        PTEDUCAT_bl= unique(data$PTEDUCAT_bl),
                        AGE_bl = seq(min(data$AGE_bl), max(data$AGE_bl)),
                        PTGENDER_bl= unique(data$PTGENDER_bl),
                        MMSE_bl= unique(data$MMSE_bl),
                        CDGLOBAL_bl= unique(data$CDGLOBAL_bl))
  data$predict_subj_spec <- predict(model, data)
  spec.data <- subset(data, RID %in% rids)
  newdat$global_prediction <- predict(model, newdata=newdat, re.form=NA)
  plot <-ggplot(newdat, aes(x=new_time, y=global_prediction)) + geom_smooth(method="lm", formula = y ~x) + geom_line(data = spec.data, aes(x=new_time, y=predict_subj_spec, colour=RID)) + geom_point(data = spec.data, aes_string(x="new_time", y=y, colour="RID"))
  plot <- plot + scale_colour_discrete(guide = "none") + xlab("Time (Years)") + ylab(ylab)
  return(plot)
}

PlotSimulationDPM <- function(model, data, y, ylab) {
  newdattreat <- expand.grid(new_time=unique(data$new_time),
                             PTEDUCAT_bl= unique(data$PTEDUCAT_bl),
                             AGE_bl = seq(min(data$AGE_bl), max(data$AGE_bl)),
                             PTGENDER_bl= unique(data$PTGENDER_bl),
                             MMSE_bl= unique(data$MMSE_bl),
                             CDGLOBAL_bl= unique(data$CDGLOBAL_bl))
  newdatplacebo <- expand.grid(new_time=unique(data$new_time),
                               PTEDUCAT_bl= unique(data$PTEDUCAT_bl),
                               AGE_bl = seq(min(data$AGE_bl), max(data$AGE_bl)),
                               PTGENDER_bl= unique(data$PTGENDER_bl),
                               MMSE_bl= unique(data$MMSE_bl),
                               CDGLOBAL_bl= unique(data$CDGLOBAL_bl))
  newdattreat$treat<-1
  newdatplacebo$treat<-0
  newdat <- rbind(newdattreat, newdatplacebo)
  newdat$treat <- factor(newdat$treat)
  newdat$global_prediction <- predict(model, newdata=newdat, re.form=NA)
  plot <-ggplot(newdat, aes(x=new_time, y=global_prediction, linetype=treat)) + geom_smooth(method="lm", formula = y ~x) 
  plot <- plot + xlab("Time (Years)") + ylab(ylab)
  return(plot)
}


GetRelContributions <- function(model, data) {
  neuropat <- data[,c("LewyPos", "CAAPos", "TDP43Pos")]
  fullrate <- fixef(model)["new_time"]
  ranefmodel<- ranef(model)[[1]]
  ranefrate <- ranefmodel[[2]]
  ranef.ids <- rownames(ranefmodel)
  subjectdecline <- fullrate + ranefrate
  frame <- data.frame("RID" = ranef.ids,
                      "Rate_Decline" = subjectdecline)
  bline <- data[!duplicated(data$RID),]
  bline <- bline[,c("RID", "LewyPos", "CAAPos", "TDP43Pos")]
  frame <- merge(frame, bline, by="RID")
  frame$AD <- 1
  relpercentages <- data.frame("Lewy" = rep(100, nrow(frame)), 
                               "CAA" = rep(29, nrow(frame)),  
                               "TDP" = rep(52, nrow(frame)), 
                               "AD" = rep(100, nrow(frame)))
  cogstatus  <- data.frame(lapply(frame[,c("LewyPos", "CAAPos", "TDP43Pos", "AD")], function(x) as.numeric(as.character(x))))
  totalcognitivedecline <- cogstatus * relpercentages
  colnames.append <- colnames(totalcognitivedecline)
  colnames(totalcognitivedecline)<- paste(colnames.append, "_relative_AD_perc", sep="")
  totalcognitivedecline$cum_perc_decline <- rowSums2(as.matrix(totalcognitivedecline))
  releffects <- round((totalcognitivedecline[,c(1,2,3, 4)] / totalcognitivedecline[,5]) * 100, 2)
  reldeclinenames<- paste(colnames.append, "_relative_perc_decline", sep="")
  colnames(releffects) <- reldeclinenames
  totalcognitivedecline <- cbind(totalcognitivedecline, releffects)
  totalcognitivedecline$anti_tau_relative_perc <- totalcognitivedecline$AD_relative_perc / 2
  totalcognitivedecline <- cbind(frame, totalcognitivedecline)
  rel_rates <- (totalcognitivedecline[,reldeclinenames] / 100) * totalcognitivedecline[,"Rate_Decline"]
  colnames(rel_rates) <- paste(colnames.append, "_relative_contribution_decline")
  totalcognitivedecline <- cbind(totalcognitivedecline, rel_rates)
  relratesmean <- unlist(Map(function(x) {round(mean(x), 3)}, rel_rates))
  sdperc <- unlist(Map(function(x) {round(sd(x), 3)}, totalcognitivedecline[,reldeclinenames]))
  meanperc <- unlist(Map(function(x) {round(mean(x), 3)}, totalcognitivedecline[,reldeclinenames]))
  num_subjects <- unlist(Map(function(x) {length(which(x == 1))}, totalcognitivedecline[,c("LewyPos", "CAAPos", "TDP43Pos", "AD")]))
  totalratedecline <- round(fixef(model)["new_time"], 3)
  rates.data <- list("Path" = names(relratesmean),
                     "Mean_Rate" = relratesmean,
                     "SD_Perc" = sdperc,
                     "Mean_Perc" = meanperc,
                     "Num_Subjects" = num_subjects,
                     "Total_Decline" = totalratedecline)
  return(list("Table" = totalcognitivedecline,
              "Rates_Data" = rates.data))
  
}


BuildNeuroCountPlot <- function(RelContr) {
  trytable <- RelContr$Table[,c("LewyPos", "CAAPos", "TDP43Pos", "AD")]
  OnlyAD <- rep(0, nrow(trytable))
  OnlyAD[which(trytable$LewyPos==0 & trytable$CAAPos==0 & trytable$TDP43Pos==0)] <- 1
  trytable <- model.matrix(AD ~ .^2, data=trytable)
  colnames(trytable) <- c("AD", "Lewy", "CAA", "TDP43", "Lewy & CAA", "Lewy & TDP43", "CAA & TDP43")
  trytable <- as.data.frame(trytable)
  trytable$`Only AD` <- OnlyAD
  trytable
  tryframe <- unlist(Map(function(x) {length(which(x==1))}, trytable))
  tryframe <- data.frame("Group" = names(tryframe),
                         "Counts" = unname(tryframe))
  
  tryframe$Group <- factor(tryframe$Group,                                    # Factor levels in decreasing order
                           levels = tryframe$Group[order(tryframe$Counts, decreasing = TRUE)])
  plotcolors <-
    setNames( c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
              , c("AD", "Lewy", "CAA", "TDP43", "Lewy & CAA", "Lewy & TDP43", "CAA & TDP43", "Only AD"))
  return(ggplot(tryframe, aes(y=Counts, x=Group, fill=Group)) + geom_bar(stat = "identity")+  scale_fill_manual(values = plotcolors) +
           
           theme(axis.title.x=element_blank()))
}


BuildSimulationDataTable <- function(list) {
  simstats <- list$fullstats
  simstats <- simstats[,c("nlevels", "successes", "mean", "lower", "upper", "Group")]
  colnames(simstats) <- c("SampleSize", "Successes", "Mean", "Lower (95% CI)", "Upper (95% CI)", "Group")
  simstatslist <- split(simstats, simstats$Group)
  return(simstatslist)
}


DPMPlots <- function(list, ylab, ylim.low, ylim.high) {
  return.list <- list()
  for(i in 1:length(list)) { 
    mod <- list[[i]]
    get.data <-  get_model_data(model = mod, terms = "new_time", type = "int")
    get.data$Treatment <- get.data$group
    plot.dpm <- ggplot(get.data, 
                       aes(x=x, y=predicted, fill=Treatment)) + geom_point() + geom_line() + geom_ribbon(aes(ymin=conf.low, ymax=conf.high), stat = "identity",  alpha=.5) + xlab("Time From Baseline (Years)") +ylab(ylab) + ylim(ylim.low, ylim.high)
    return.list[[i]] <- plot.dpm
  }
  names(return.list) <- names(list)
  return(return.list)
}


GetConfInt <- function(data_succ) {
  names_data  <- colnames(data_succ)
  returnframe <- data.frame()
  for(i in 1:ncol(data_succ)) {
    x         <- unname(unlist(data_succ[i]))
    succ      <- length(which(x == "y"))
    total     <- length(x)
    mean      <- succ / total
    binom     <- binom.test(succ, total)
    confinter <-(stats::binom.test(succ, total))
    confinter <- as.numeric(confinter$conf.int)
    row.val   <- c(mean, confinter)
    returnframe <- rbind(returnframe, row.val)
  }
  colnames(returnframe) <- c("mean", "ci.low", "ci.high")
  rownames(returnframe) <- names_data
  return(returnframe)
}


CombineOutcomes <- function(insidelist) {
  pvals      <- which(names(insidelist) == "pval")
  treatments <- which(names(insidelist) == "Treatment")
  placebs    <- which(names(insidelist) == "Placebo")
  vals       <- insidelist[pvals]
  treatments <- insidelist[treatments]
  placebs    <- insidelist[placebs]
  treatments <- do.call(cbind, treatments)
  placebs    <- do.call(cbind, placebs)
  vals       <- unlist(vals)
  return(list("Treatment" = treatments,
              "Placebo" = placebs,
              "pval" = vals))
}


CalcSuccesses <- function(list, rows) {
  init_pval_data <- data.frame(matrix(nrow = rows))
  for(i in 1:length(list)) {
    sublist <- list[[i]]
    init_pval_data <- cbind(init_pval_data, sublist$pval)
  }
  init_pval_data[,1] <- NULL
  colnames(init_pval_data) <- names(list)
  init_pval_data[init_pval_data <= 0.05]      <- "y"
  init_pval_data[init_pval_data != "y"]       <- "n"
  return(init_pval_data)
}


CombineIters<- function(list, rows) {
  init_pval_data <- data.frame(matrix(nrow = rows))
  for(i in 1:length(list)) {
    sublist <- list[[i]]
    init_pval_data <- cbind(init_pval_data, sublist$pval)
  }
  init_pval_data[,1] <- NULL
  colnames(init_pval_data) <- names(list)
  return(init_pval_data)
}


StratifyContinuous <- function(longdata, stratcols) {
  bline      <- longdata[!duplicated(longdata$RID),]
  bline      <- StratifyContVar(bline, stratcols = stratcols)
  post.strat <- paste(stratcols, "_strat", sep="")
  bline      <- bline[,c("RID", post.strat)]
  longdata   <- merge(longdata, bline, by="RID", all.x=TRUE)
  return(longdata)
}


BalanceDiagnostics <- function(outer) {
  prop_test_frame <- Map(PropTestIter, outer)
  return(prop_test_frame)
}

PropTestIter <- function(outer_sublist) {
  tr <- outer_sublist$Treatment
  pl <- outer_sublist$Placebo
  ss <- outer_sublist$sample_size
  rnames <- rownames(tr)
  init_df <- data.frame(matrix(nrow = dim(tr)[1]))  
  if(!all(dim(tr) == dim(pl))) {
    stop("Dataframe dimensions are not identical")
  }
  for(i in 1:ncol(tr)) {
    col_tr    <- tr[,i]
    col_pl    <- pl[,i]
    col.tests <- c()
    for(j in 1:nrow(tr)) {
      cell_tr <- col_tr[[j]]
      cell_pl <- col_pl[[j]]
      test    <- prop.test(c(round(cell_tr*ss), 
                             round(cell_pl*ss)), 
                             c(ss, ss))$p.value
      col.tests <- append(col.tests, test)
    }
    init_df <- cbind(init_df, col.tests)
  }
  init_df[,1] <- NULL
  rownames(init_df) <- rnames
  return(init_df)
}


TrPlMeanCI <- function(outer_list) {
  nsim       <- ncol(outer_list$Treatment)
  rnames_tr  <- rownames(outer_list$Treatment)
  rmeans_tr  <- rowMeans2(as.matrix(outer_list$Treatment))
  rSd_tr     <- rowSds(as.matrix(outer_list$Treatment))
  rSE_tr     <- rSd_tr / sqrt(nsim)
  ci_low_tr  <- rmeans_tr - 1.96 * rSE_tr
  ci_high_tr <- rmeans_tr + 1.96 * rSE_tr
  tr.data <- data.frame("mean"    = rmeans_tr,
                        "se"      = rSE_tr,
                        "ci_low"  = ci_low_tr,
                        "ci_high" = ci_high_tr,
                        "cov"     = rnames_tr)
  tr.data$group <- "Treatment"
  
  rnames_pl  <- rownames(outer_list$Placebo)
  rmeans_pl  <- rowMeans2(as.matrix(outer_list$Placebo))
  rSd_pl     <- rowSds(as.matrix(outer_list$Placebo))
  rSE_pl     <- rSd_pl / sqrt(nsim)
  ci_low_pl  <- rmeans_pl - 1.96 * rSE_pl
  ci_high_pl <- rmeans_pl + 1.96 * rSE_pl
  pl.data    <- data.frame("mean"    = rmeans_pl,
                           "se"      = rSE_pl,
                           "ci_low"  = ci_low_pl,
                           "ci_high" = ci_high_pl,
                           "cov"     = rnames_pl)
  pl.data$group <- "Placebo"
  returndata <- rbind(tr.data, pl.data)
  return(returndata)
}


MeanCICovars <- function(outer) {
  returnlist <- Map(TrPlMeanCI, outer)
  return(returnlist)
}

ManualSimulation <- function(formula_largemodel, largemodel, formula_smallmodel, smallmodel, sample_sizes, nsim, data, trial_duration=NULL, t1errorsim) {
  # load in functions from GlobalEnv into current enviornment for clusterExport
  force(RandomizeTreatment2)
  force(pbkrtest::KRmodcomp)
  force(pbkrtest::getKR)
  force(lme4::lmer)
  force(lme4::lmerControl)
  force(dfoptim::nmkb)
  force(dplyr::bind_cols)
  force(splitstackshape::stratified)
  force(CalcProportionPos)
  force(`%notin%`)
  force(setTxtProgressBar)
  force(lmerTest::as_lmerModLmerTest)
  force(dplyr::sample_n)
  opts <- list(chunkSize=10)
  #################
  cat("Beginning simulation")
  cat("\n")
  init_iter_list                      <-  list()
  init_significance_list              <-  list()
  init_props_list_treatment           <-  list()
  init_props_list_placebo             <-  list()
  data_extended <- data
  form_lm_split <- strsplit(formula_largemodel, "~")[[1]][2]
  iter_form_lm  <- paste("large_model_response", form_lm_split, sep="~")
  levels_extended <- unique(data_extended$RID)
  if(!is.null(t1errorsim)) {
    fixef(largemodel)["new_time:treat1"] <- 0
  }
  #define inner loop for parallelization
    .nsiminnerloop <- function(j) {
    cat("\r", j, " out of ", nsim, " complete", sep = "")
    sample.levels            <-  sample(levels_extended, size = i * 2)
    data_sample              <-  subset(data_extended, RID %in% sample.levels)
    sample_baseline          <-  data_sample[!duplicated(data_sample$RID), ]
    
    #split into treatment and placebo groups while balancing subjects
    treatment.out            <-  RandomizeTreatment2(sample_baseline, data_sample)
    prop                     <-  treatment.out[["props"]]
    data_sample_treated      <-  treatment.out[["data"]]
    
    simulate_response_largemodel <- simulate(largemodel, newdata = data_sample_treated, allow.new.levels=TRUE, use.u=FALSE)
    
    refit_data_outcomes          <- data.frame("large_model_response" = simulate_response_largemodel)
    colnames(refit_data_outcomes) <- c("large_model_response")
    fit_iter_data                 <- bind_cols(refit_data_outcomes, data_sample_treated) 
  
    refit_large                   <- lme4::lmer(formula = as.formula(iter_form_lm), 
                                                data = fit_iter_data, REML = TRUE, control = lme4::lmerControl(optimizer = "nmkbw"))
    
    pval <-  as.numeric(summary(lmerTest::as_lmerModLmerTest(refit_large))[["coefficients"]][,"Pr(>|t|)"]["new_time:treat1"])

    
    
    return(list("pval"      = pval,
                "Treatment" = prop[["Treatment"]],
                "Placebo"   = prop[["Placebo"]]))
  }
  pb <- txtProgressBar(min = 1, max=nsim, style = 1)
  envlist    <- mget(ls())
  envlist    <- names(envlist)
  env.append <- c("RandomizeTreatment2","KRmodcomp",
                  "getKR","lmer","bind_cols","stratified", "sample_n",
                  ".nsiminnerloop", 
                  "CalcProportionPos", 
                  "%notin%")
  
  #create cluster
  cl <- makeCluster(1, outfile="")
  doSNOW::registerDoSNOW(cl)
  envlist <- append(envlist, env.append)
  clusterExport(cl, envlist, envir = environment())
  t1 <- Sys.time()
  outer <- foreach(i = sample_sizes, .options.nws=opts) %:%
            foreach(j = 1:nsim, .combine='c', .inorder=FALSE) %dopar% {
              setTxtProgressBar(pb, j)
              .nsiminnerloop(j)
            }
  t2     <- Sys.time()
  stopCluster(cl)
  names(outer) <- paste("SS_", sample_sizes, sep="")
  outer        <- Map(CombineOutcomes, outer)
  for(i in 1:length(outer)) {
   outer[[i]][["sample_size"]] <- sample_sizes[i]
  }
  ftestdf                 <- CombineIters(outer, nsim)
  successes               <- CalcSuccesses(outer, nsim) 
  conf.inter              <- GetConfInt(successes)
  balance.cov.diagnostics <- BalanceDiagnostics(outer)
  mean_cis_covariates     <- MeanCICovars(outer)
  timerun <- difftime(t2, t1, units = "mins")
  return(list("Successes"   = successes, 
              "power.data"  = conf.inter, 
              "cov.balance" = outer,
              "mean_cis_covariate" = mean_cis_covariates,
              "prop.tests"  = balance.cov.diagnostics,
              "time_to_run" = timerun,
              "ftestdf" = ftestdf))
}





CombineSimPlots_ManualSimulation <- function(simlist, limits) {
  nonenrich <- simlist[[1]] 
  name1 <- names(simlist)[1]
  sumstats_non <- nonenrich$power.data
  sumstats_non$Group <- name1
  sumstats_non$nlevels <- limits
  plot.df <- sumstats_non
  for(i in 2:length(simlist)){
    enriched <- simlist[[i]]
    name <- names(simlist)[i]
    enriched <- enriched$power.data
    enriched$Group <- name
    enriched$nlevels <- limits
    plot.df <- rbind(plot.df, enriched)
  }
  gplot.sim <- ggplot(data = plot.df, aes(x=nlevels, y=mean * 100, colour=Group)) +geom_errorbar(aes(ymin=ci.low * 100, ymax=ci.high * 100)) + geom_line() + geom_point() + scale_x_discrete(limits=limits)
  gplot.sim <- gplot.sim + xlab("Sample Size per Arm") + ylab("Statistical Power (%)") +ylim(0, 100) + geom_hline(yintercept = 80, linetype="dashed")
  return(list("plot" = gplot.sim, "fullstats" = plot.df))
} 


PlotCovariateBalance <- function(sim_fitted, sample_size) {
  balancedata <- sim_fitted$mean_cis_covariate[[sample_size]]
  covs_renamed <-  c("fullcaa.0"= "CAA(-)", "fulllewy.0" = "Lewy(-)",
                     "fulltdp43.0" = "TDP43(-)", "PTGENDER.Female" = "Sex(F)",
                     "Age_bl_strat" = "Age(Cat:0)", "PTEDUCAT_bl_strat.0" = "Educ(Cat:0)", 
                     "MMSE_bl_strat.0" = "MMSE(Cat:0)", "CDGLOBAL_bl.0.5"= "CDR(0.5)", "TauPos_full_bl.0" = "Tau(-)")
  balancedata$cov2 <- covs_renamed[balancedata$cov]
  plot.gg <- ggplot(balancedata, aes(x=cov2, y=mean, colour=group)) + geom_point() + geom_errorbar(aes(ymin=ci_low, ymax=ci_high))
  plot.gg <- plot.gg + scale_colour_discrete(name= "Treatment Group") + ylab("Proportion of Covariate per Group") +xlab("Covariate")
  return(plot.gg)
}


All.pair.ttests <- function(data.list, col.name) {
  c1 <- data.list[[1]][[col.name]]
  c2 <- data.list[[2]][[col.name]]
  c3 <- data.list[[3]][[col.name]]
  
  t1 <- t.test(c1, c2)
  t2 <- t.test(c1, c3)
  t3 <- t.test(c2, c3)
  return(list(t1, t2, t3))
}


All.pair.prop.test <- function(data.list, col.name, level) {
  c1 <- data.list[[1]][[col.name]]
  c2 <- data.list[[2]][[col.name]]
  c3 <- data.list[[3]][[col.name]]
  
  c1.length <- length(c1)
  c2.length <- length(c2)
  c3.length <- length(c3)
  
  c1.pos <- length(which(c1 == level))
  c2.pos <- length(which(c2 == level))
  c3.pos <- length(which(c3 == level))
  
  t1 <- prop.test(c(c1.pos, c2.pos), c(c1.length, c2.length))
  t2 <- prop.test(c(c1.pos, c3.pos), c(c1.length, c3.length))
  t3 <- prop.test(c(c2.pos, c3.pos), c(c2.length, c3.length))
  return(list(t1, t2, t3))
}





Adni_Age <- function(df) {
  bdayposit <- paste(df$PTDOBYY, df$PTDOBMM, "01", sep="-")
  bdayposit <- as.POSIXct(bdayposit, format="%Y-%m-%d")
  timediff <- round(difftime(df$EXAMDATE, bdayposit, units = "weeks") / 52, 2)
  return(timediff)
}


PlotT1Error <- function(data, sample_sizes, ycol, title) {
data_5 <- as.numeric(Map(function(x){length(which(x <= 0.05)) / length(x)}, data))
data_4 <- as.numeric(Map(function(x){length(which(x <= 0.04)) / length(x)}, data))
data_3 <- as.numeric(Map(function(x){length(which(x <= 0.03)) / length(x)}, data))
data_2 <- as.numeric(Map(function(x){length(which(x <= 0.02)) / length(x)}, data))
data_1 <- as.numeric(Map(function(x){length(which(x <= 0.01)) / length(x)}, data))
data_errors <- as.data.frame(do.call(cbind, list(data_1, data_2, data_3, data_4, data_5)))
colnames(data_errors) <-  c("c_01", "c_02","c_03","c_04","c_05")
rownames(data_errors) <- paste("SS_", sample_sizes, sep="")
data_errors <- as.data.frame(t(data_errors))
data_errors <- data_errors*100
data_errors$sig_level <- c(.01, .02, .03, .04, .05)
gg <- ggplot(data_errors, aes_string(y=ycol, x="sig_level")) + geom_point(shape="triangle") + ylab("Type I Error (%)") + xlab("Significance Level") + scale_x_reverse() +ylim(0, 6)
print(data_errors)
gg <- gg+ labs(title = title)
return(gg)
}





CombineTrialDuration <- function(power_list) {
  power_list_object <- Map(function(x){x$power.data}, power_list)
  trialdurations <- seq(1,4,by=.5)
  if(length(power_list_object) == 6) {
    trialdurations <- trialdurations[!trialdurations==1.5]
  }
  for(i in 1:length(power_list)) {
    power_list_object[[i]]$trial_duration <- trialdurations[i]
    power_list_object[[i]]$samplesize<- substr(rownames(power_list_object[[i]]), start = 4, stop = 6)
  }
  power_list_object <- do.call(rbind, power_list_object)
  min.samp <- min(as.numeric(power_list_object$samplesize))
  max.samp <- max(as.numeric(power_list_object$samplesize))

  
  my3d.adas.un.tau <- plot_ly(x=as.numeric(power_list_object$trial_duration), y=as.numeric(power_list_object$samplesize), z=as.numeric(power_list_object$mean*100), type="scatter3d", mode="markers")
  peep.tau<- expand.grid("v1" = seq(1, 4, by=.5), "v2" = seq(min.samp, max.samp, by=20))
  z.tau <- t(outer(peep.tau$v2, peep.tau$v1, function(x,y)  0*x +0*y + 80))
  my3d.adas.un.tau <- add_trace(my3d.adas.un.tau, y=~peep.tau$v2, x=~peep.tau$v1, z=z.tau, type="surface", colorscale=list(c(0, 1), c("red", "red")), showscale=FALSE) 
  my3d.adas.un.tau <- my3d.adas.un.tau %>% layout(scene=list(xaxis=list(title="Trial Duration (Years)"), yaxis=list(title="Sample Size per Arm"), zaxis=list(title="Statistical Power (%)")))
  return(list("plot" = my3d.adas.un.tau,
              "peep.tau" = peep.tau,
              "z.tau" = z.tau,
              "data" = power_list_object))
}


hybridapproach <- function(formula, model, time=c(0,.5,1,1.5,2), nsim, 
                           pct.change=.5, delta=NULL, parameter="new_time", sample.prop =.8) {
  cat("Beginning simulation")
  cat("\n")
  data_extended               <-  simr::getData(model)
  levels.id                   <- levels(factor(data_extended$RID))
  return.n <- list()
  for(i in 1:nsim) {
    print(i)
    levels.sample <- sample(levels.id, round((length(levels.id) * sample.prop)))
    data_subset <- subset(data_extended, RID %in% levels.sample)
    data_subset$RID <- factor(data_subset$RID)
    sample_baseline          <-  data_subset[!duplicated(data_subset$RID), ]
    treatment.out            <-  RandomizeTreatment2(sample_baseline, data_subset)
    data_sample_treated      <-  treatment.out[["data"]]
    refit_small               <- lme4::lmer(formula = as.formula(formula), 
                                            data = data_subset, REML = TRUE)
    fixef(refit_small)["treat1"] <- 0
    fixef(refit_small)["new_time:treat1"] <- ((fixef(refit_small)["new_time"]) * -1) * pct.change
    delta <-   fixef(refit_small)["new_time:treat1"] 
    long.p <- longpower::lmmpower(refit_small, parameter=parameter, t=time, delta=delta, power=.8, sig.level=.05)
    return.n[[i]] <- long.p
  }
  return(return.n)
  }

checkmean <- function(hybrid) {
  n <- unlist(Map(function(x){(x$N) / 2}, hybrid))
  return(c("mean" = mean(n), "sd" = sd(n)))
}

Cal <- function(mean.points) {
  m.mean     <- (mean.points[2] - mean.points[1]) / (mean.points[4] - mean.points[3])
  mean.shift <- (80 - mean.points[1]) / m.mean
  new.x <- mean.points[3] + mean.shift
  return(new.x)
}



props.covs <- function(list) {
  if(length(list)==2) {
    cutoff <- c(list[1] / (list[1] + list[2]))
  }
   else if(length(list)==3) {
    cutoff <- c((list[1] / (list[1] + list[2] + list[3])), (list[2] / (list[1] + list[2] + list[3])))
  } else {
    cutoff <- NULL
  }
  return(cutoff)
}


CalcCutoff <- function(mean, sd, pi) {
  return.vec <- c()
  for(i in 1:length(pi)) {
    pi.i <- pi[[i]]
    cutoff <- mean + (sd * qnorm(pi.i))
    cutoff <- exp(cutoff)
    return.vec <- append(return.vec,cutoff)
  }
  return(return.vec)
}



CalculateSampleAtPower_markdown <- function(data, y=80) {
  lowerval <- max(which(data$mean < .8))
  upperval <- min(which(data$mean >= .8))
  data$samplesize <- as.numeric(sapply(strsplit(rownames(data), "_"), "[[", 2))
  xvals <- c(data["samplesize"][lowerval,], data["samplesize"][upperval,])
  samplesize <- Cal(c(data["mean"][lowerval,]*100, data["mean"][upperval,]*100, xvals))
  low <- Cal(c(data["ci.low"][lowerval,]*100, data["ci.low"][upperval,]*100, xvals))
  hi<- Cal(c(data["ci.high"][lowerval,]*100, data["ci.high"][upperval,]*100, xvals))
  return.val <- paste(round(samplesize,2), " (", round(hi,2), ",", round(low,2),")", sep="")
  return(return.val)
}



DefineMVND <- function(data, n) {
  thresholds <- c()
  data <- data[,c("PTEDUCAT_bl",  "AGE_bl", "PTGENDER", "MMSE_bl", "CDGLOBAL_bl", "LewyPos", "CAAPos", "TDP43Pos", "TauPos_bl")]
  names.vec <- c("PTGENDER", "CDGLOBAL_bl", "LewyPos", "CAAPos", "TDP43Pos", "TauPos_bl")
  continuousdata <- data[,c("PTEDUCAT_bl",  "AGE_bl", "MMSE_bl")]
  for(i in 1:length(names.vec)) {
    name.i <- names.vec[i]
    levels.i <- levels(factor(data[,name.i]))
    newvec <- mapvalues(data[,name.i], from=levels.i, to = c(1:length(levels.i)))
    continuousdata[name.i] <- newvec
  }
  continuousdata <- continuousdata %>% mutate_all(as.character)
  continuousdata <- continuousdata %>% mutate_all(as.numeric)
  empir.distr    <- continuousdata
  cov.b4              <- empir.distr
  continuousdata      <- as.matrix(continuousdata)
  continuousdatalog   <- log(continuousdata)
  covmatrix           <- cov(continuousdatalog)
  means               <- colMeans(continuousdatalog)
  simcovs             <- MASS::mvrnorm(n=n, mu=means, Sigma = covmatrix)
  simcovs             <- exp(simcovs)
  cov.after           <- simcovs
  cont.thresh         <- continuousdatalog[,names.vec]
  means.log           <- colMeans2(cont.thresh)
  sds.log             <- colSds(cont.thresh)
  names(means.log)    <- names(sds.log) <- names.vec
  pi                  <- Map(table, empir.distr[,names.vec])
  pi                  <- Map(props.covs, pi)
  simcovs             <- as.data.frame(simcovs)
  simcovs1 <- simcovs
  for(i in names.vec) {
    levels.i <- levels(factor(data[,i]))
    if(length(levels.i) == 1) {
      simcovs[i] <- 1
    } else {
    threshold <- CalcCutoff(means.log[[i]], sds.log[[i]], pi[[i]])
     if(length(threshold)==1) {
    simcovs[i][which(simcovs[i]  <=  threshold[1]),] <- levels.i[1]
    simcovs[i][which(simcovs[i]  !=  levels.i[1]),]  <- levels.i[2]
    } else {
      simcovs[i][which(simcovs[i] <= threshold[1]),] <- levels.i[1]
      simcovs[i][which(simcovs[i]  > threshold[1] & simcovs[i]  <= threshold[2]),] <- levels.i[2]
      simcovs[i][which(simcovs[i]  > threshold[2]),] <- levels.i[3]
     }
    }
  }
  return(list("cov.b4" = data, "simcovs"=simcovs))
}


ExtendLongitudinal <- function(data, trial_duration) {
  new_time <- seq(0, trial_duration, by=.5)
  rows <- nrow(data)
  data <- data %>% slice(rep(1:n(), each=length(new_time)))
  rids <- rep(1:rows, length(new_time))
  rids <- rids[order(rids, decreasing = FALSE)]
  new_time_rows <- rep(new_time, rows)
  data$RID <- factor(rids)
  data$new_time <- new_time_rows
  return(data)
}


CompareDistributions <- function(mvnd.out) {
  b4    <- mvnd.out$cov.b4
  after <- mvnd.out$simcovs
  b4 <- b4[,c("PTEDUCAT_bl",  "AGE_bl", "MMSE_bl", "CDGLOBAL_bl", "LewyPos", "CAAPos", "TDP43Pos", "PTGENDER")]
  after <- after[,c("PTEDUCAT_bl",  "AGE_bl",  "MMSE_bl", "CDGLOBAL_bl", "LewyPos", "CAAPos", "TDP43Pos", "PTGENDER")]
  b4 <- b4 %>% mutate_all(as.character)
  b4["PTGENDER"][which(b4$PTGENDER=="Male"), ] <- 1
  b4["PTGENDER"][which(b4$PTGENDER=="Female"), ] <- 2
  b4 <- b4 %>% mutate_all(as.numeric)
  after <- after %>% mutate_all(as.character)
  after["PTGENDER"][which(after$PTGENDER=="Male"), ] <- 1
  after["PTGENDER"][which(after$PTGENDER=="Female"), ] <- 2
  after <- after %>% mutate_all(as.numeric)
  dist.test <- cramer::cramer.test(as.matrix(b4), as.matrix(after), sim="permutation")
  return(dist.test)
}


OverallvsTau <- function(dpm.list.overall, dpm.list.tau, ylab, ymin, ymax) {
  keeplist <- list()
  names.list <- names(dpm.list.overall)
  for(i in 1:length(dpm.list.overall)) {
    d1 <- dpm.list.overall[[i]]$data
    d2 <- dpm.list.tau[[i]]$data
    d1$`Treatment Definition` <- "Overall Decline"
    d2$`Treatment Definition` <- "Tau Associated Decline"
    d <- rbind(d1, d2)
    d <- subset(d, Treatment==0)
    gg <- ggplot(d, aes(x=x, y=predicted, fill=`Treatment Definition`)) + geom_point(aes(colour=`Treatment Definition`)) +geom_line(aes(colour=`Treatment Definition`))+ geom_ribbon(aes(ymin=conf.low, ymax=conf.high), stat = "identity",  alpha=.25)
    gg <- gg + xlab("Time From Baseline (Years)") + ylab(ylab) + ylim(ymin, ymax)
    gg <- gg +  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), stat = "identity",  alpha=.25)
    keeplist[[i]] <- gg
  }
  names(keeplist) <- names.list
  return(keeplist)
}
