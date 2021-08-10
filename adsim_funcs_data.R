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
    min.time <- which()
    subj$new_time <- subj[[timecol]] - min.time
    return.list[[i]] <- subj
  }
  return.list <- do.call(rbind, return.list)
  return(return.list)
}



TimeSinceBaselineValidAmy <- function(data, timecol) {
  subjlist <- split(data, data$RID)
  return.list <- list()
  for(i in 1:length(subjlist)) {
    subj <- subjlist[[i]]
    if(1 %in% subj$adas_csf_valid | 1 %in% subj$adas_pet_valid) {
    val.amy   <- which(subj$adas_csf_valid == 1 | subj$adas_pet_valid == 1)
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

 

CreateBaselineVar <- function(data, column.name) {
  splitsubj <- split(data, data$RID)
  comb.list <- list()
  for(i in 1:length(splitsubj)) {
    bline.name <- paste(column.name, "_bl", sep="")
    subj <- splitsubj[[i]]
    rows <- nrow(subj)
    bline.row <- which.min(subj$new_time)
    bline.val <- subj[column.name][bline.row,]
    bline.vec <- rep(bline.val, rows)
    subj[bline.name] <- bline.vec
    comb.list[[i]] <- subj
  }
  comb.list <- do.call(bind_rows, comb.list)
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


SampleSizeSimulation <- function(sim.data, formula, fcompare_str, efficacy=.5, breaks) {
  all.rids <- levels(sim.data[["RID"]])
  nrids <- nlevels(sim.data[["RID"]])
  n.sample <- round(nrids/2)
  treat.rids <- sample(all.rids, n.sample, prob = rep(.5, length(all.rids)), replace = FALSE)
  base.rids <- subset(all.rids, all.rids %notin% treat.rids)
  sim.data$treat <- rep(NA, nrow(sim.data))
  treat.rows <- which(sim.data$RID %in% treat.rids)
  base.rows <- which(sim.data$RID %in% base.rids)
  sim.data["treat"][treat.rows,] <- 1
  sim.data["treat"][base.rows,] <- 0
  sim.data$treat <- factor(sim.data$treat)
  sim.model <- lmer(as.formula(formula), data = sim.data)
  sim.model.large <- sim.model
  fixef(sim.model.large)["treat1"] <- 0
  treat.ef <- fixef(sim.model.large)["new_time"] * efficacy
  treat.ef <- treat.ef * -1
  fixef(sim.model.large)["treat1:new_time"] <- treat.ef
  sim_ext_class     <- extend(sim.model.large, along="RID", n=max(breaks))
  p_curve_treat_sim <- powerCurve(sim_ext_class, test=fcompare(as.formula(fcompare_str)), along="RID", breaks=breaks)
  summ_sim<- summary(p_curve_treat_sim)
  summ_sim$mean <- summ_sim$mean*100
  summ_sim$lower <- summ_sim$lower*100
  summ_sim$upper <- summ_sim$upper*100
  gplot.sim <- ggplot(data = summ_sim, aes(x=nlevels, y=mean)) +geom_errorbar(aes(ymin=lower, ymax=upper)) + geom_line() + geom_point() + scale_x_discrete(limits=breaks)
  gplot.sim <- gplot.sim + xlab("Sample Size per Arm") + ylab("Statistical Power (%)")
  return.list <- list("model" = sim.model.large,
                      "power_curve_output" = p_curve_treat_sim,
                      "summary_stats" = summ_sim,
                      "power_curve_plot" = gplot.sim)
  return(return.list)
}


CombineSimPlots <- function(nonenrich, enrich, limits) {
  sumstats_non <- nonenrich$summary_stats
  sumstats_enr <- enrich$summary_stats
  sumstats_non$Group <- rep("Non-Enriched", nrow(sumstats_non))
  sumstats_enr$Group <- rep("Enriched", nrow(sumstats_enr))
  fullsumstats <- rbind(sumstats_enr, sumstats_non)
  gplot.sim <- ggplot(data = fullsumstats, aes(x=nlevels, y=mean, colour=Group)) +geom_errorbar(aes(ymin=lower, ymax=upper)) + geom_line() + geom_point() + scale_x_discrete(limits=limits)
  gplot.sim <- gplot.sim + xlab("Sample Size per Arm") + ylab("Statistical Power (%)")
  return(list("plot"=gplot.sim, "fullstats"=fullsumstats))
} 




























