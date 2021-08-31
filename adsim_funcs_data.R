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



TimeSinceBaselineValidAmy <- function(data, timecol) {
  subjlist <- split(data, data$RID)
  return.list <- list()
  for(i in 1:length(subjlist)) {
    subj <- subjlist[[i]]
    if(!all(is.na(subj$AmyPos))) {
    val.amy   <- which(!is.na(subj$AmyPos))
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
  model<- firstelement$model
  names1 <- c(names(sim.list)[1])
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
  names(screen.fail) <- names(sim.list)
  plot.data$Treatment <- plot.data$treat
  plot.data$Enrichment <- plot.data$group
  
  gplot <- ggplot(plot.data, aes(y=predicted, x=new_time, colour=Enrichment, linetype = Treatment)) + geom_smooth(method = "lm", aes(fill=Enrichment), alpha=.1) 
  gplot <- gplot  + xlab("Time (Months)") + ylab(yaxislab_dpm)
  return(list("plot"=gplot, "screenfail"=screen.fail))
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
  data <- TimeSinceBaselineValidAmy(data, "M_vis")
  data <- subset(data, new_time <= 24)
  }


CalculateCognitiveLoss <- function(model) {
  
}
