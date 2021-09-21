library(survey)


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
  for(i in 1:nrow(data)) {
    if(!is.na(data["TDP_pos_path"][i,])) {
      data["fulltdp43"][i,] <- data["TDP_pos_path"][i,]
    } else if(!is.na(data["TDP43"][i,])) {
      data["fulltdp43"][i,] <- data["TDP43"][i,]
    }
    
    if(!is.na(data["Lewy_pos_path"][i,])) {
      data["fulllewy"][i,] <- data["Lewy_pos_path"][i,]
    } else if(!is.na(data["LEWY"][i,])) {
      data["fulllewy"][i,] <- data["LEWY"][i,]
    }
    
    if(!is.na(data["CAA_path"][i,])) {
      data["fullcaa"][i,] <- data["CAA_path"][i,]
    } else if(!is.na(data["CAA"][i,])) {
      data["fullcaa"][i,] <- data["CAA"][i,]
    }
  }
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
  data <- data[order(data$RID, data$M_vis, decreasing = FALSE),]
  data <- TimeSinceBaseline(data, "M_vis")
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



SampleSizeSimulation2 <- function(sim.data, formula, compare_str, breaks, yaxislab_dpm, model, return_dpm=FALSE) {
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
  model            <- list[[1]]
  adjusted.decline <- list[[2]]
  fixd             <- fixef(model)
  fixd["treat1"]   <- 0
  if(treatment.effect =="controlled") {
    fixd["new_time:treat1"] <- ((adjusted.decline * .5) * -1)
  } else {
    fixd["new_time:treat1"] <- ((fixd["new_time"] * .5) * -1)
  }
  
  
  fixd                        <- fixd[c( "(Intercept)", "new_time", "treat1", "PTEDUCAT_bl", "AGE_bl", "PTGENDER_blMale")]
  sigma.mod        <- summary(model)$sigma
  varcor.mod       <- as.numeric(summary(model)$varcor[[1]])
  constr.lme       <- makeLmer(formula = as.formula(formula.model), fixef = fixd, VarCorr=list(varcor.mod), sigma = sigma.mod, data = data)
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


RandomizeTreatment2 <- function(data, stratcolumns, longdata) {
  stratifydf          <- data[,c("RID", stratcolumns)]
  stratifydf$stratvar <- interaction(stratifydf$fullcaa, stratifydf$fulllewy, stratifydf$fulltdp, stratifydf$PTGENDER, stratifydf$AGE_bl_strat)
  stratifieddata      <- stratified(stratifydf, "stratvar", size = c(.5), bothSets = FALSE)
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


RandomizeTreatment2NoPath <- function(data, stratcolumns, longdata) {
  stratifydf          <- data[,c("RID", stratcolumns)]
  stratifydf$stratvar <- interaction(stratifydf$PTGENDER, stratifydf$AGE_bl_strat)
  stratifieddata      <- stratified(stratifydf, "stratvar", size = c(.5), bothSets = FALSE)
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
    lme.fit <- lmer(as.formula(formula.model), data = newdata)
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

