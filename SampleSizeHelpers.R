
#####################################################    LIBRARIES   #####################################################    
library(plyr)
library(dplyr)
library(lme4)
library(nlme)
library(simr)
library(stringr)
library(matrixStats)
library(lmerTest)
library(MASS) 
library(splitstackshape) 
library(cramer)

#####################################################    FUNCTIONS   #####################################################    

# ALL INTERNAL FUNCTIONS

#####################################################    

#' Extracts Covariates from Model Formula
#' 
#' @param model LME model (lme)
#' @param parameter Rate of change parameter of interest (character)
#' @return Vector of covariates
#'
GetCovariates <- function(model, parameter) {
  model.formula       <- formula(model)
  model.formula       <- as.character(model.formula)
  fixed.effects.model <- gsub("\\s*\\([^\\)]+\\)","", model.formula[3])
  fixed.effects.model <- unlist(strsplit(fixed.effects.model, "\\+"))
  fixed.effects.model <- gsub(" ", "", fixed.effects.model, fixed = TRUE)
  interaction_star    <- grep("\\*", fixed.effects.model)
  interaction_colon   <- grep("\\:", fixed.effects.model)
  if(length(interaction_star) > 0) {
    int_terms <- unlist(strsplit(fixed.effects.model[grep("\\*", fixed.effects.model)], "\\*"))
    fixed.effects.model <- append(fixed.effects.model, int_terms)
    fixed.effects.model <- fixed.effects.model[fixed.effects.model  != fixed.effects.model[grep("\\*", fixed.effects.model)]]
    fixed.effects.model <- unique(fixed.effects.model)
    
  }
  if(length(interaction_colon) > 0) {
    int_terms <- unlist(strsplit(fixed.effects.model[grep("\\:", fixed.effects.model)], "\\:"))
    fixed.effects.model <- append(fixed.effects.model, int_terms)
    fixed.effects.model <- fixed.effects.model[fixed.effects.model  != fixed.effects.model[grep("\\:", fixed.effects.model)]]
    fixed.effects.model <- unique(fixed.effects.model)
  }
  fixed.effects.model <- fixed.effects.model[fixed.effects.model  != parameter]
  return(fixed.effects.model)
}

#####################################################    

#' Creates new LME with pre-defined treatment term
#' 
#' @param model LME model (lme)
#' @param parameter Rate of change parameter of interest (character)
#' @param pct.change Percent change in parameter of interest (numeric)
#' @param delta Change in parameter of interest (numeric)
#' @param data Data from pilot model fit (data.frame)
#' @param rand.effect Column name for subject
#' @return New model identical to pilot model but with treatment term defined by @param pct.change or @param delta (lme)
#'
AddTreatmentTerm <- function(model, parameter, pct.change, delta, data, rand.effect) {
  fixd           <- nlme :: fixef(model)
  treatment_term <- paste(parameter, "treat", sep = ":")
  par.loc        <- which(names(fixd) == parameter)
  
  if(!is.null(pct.change)) { 
    treatment      <- fixd[parameter] * pct.change
  } else {
    treatment      <- delta
  }
  
  fixd                      <- append(fixd, treatment, after = par.loc)
  names(fixd)[par.loc + 1]  <- treatment_term
  fixd                      <- append(fixd, 0)
  names(fixd)[length(fixd)] <- "treat"
  rand.effect.str           <- str_extract_all(as.character(formula(model))[3], 
                                               "\\([^()]+\\)")[[1]]
  response       <- as.character(formula(model))[2]
  names.terms    <- attr(terms(model), "term.labels")
  par.loc        <- which(names.terms == parameter)
  names.terms    <- append(names.terms, treatment_term, after = par.loc)
  names.terms    <- append(names.terms, "treat")
  vc             <- VarCorr(model)
  sig            <- sigma(model)
  formula.model  <- paste(names.terms, collapse = " + ")
  formula.model  <- paste(formula.model, rand.effect.str, sep = " + ")
  formula.model  <- paste(response, formula.model, sep = " ~ ")
  data$treat     <- NA
  subs           <- unique(data[ ,rand.effect])
  treat1         <- sample(subs, length(subs) * .5)
  treat1         <- which(data[ ,rand.effect] %in% treat1)
  data["treat"][treat1,] <- 1
  data["treat"][which(is.na(data["treat"])),] <- 0
  data$treat          <- as.numeric(data$treat)
  
  buildtreatmentmodel <- simr::makeLmer(as.formula(formula.model), 
                                        fixef   = unname(fixd), 
                                        VarCorr = vc, 
                                        sigma   = sig, 
                                        data    = data)
  
  fixd.rearrange      <- fixef(buildtreatmentmodel)
  fixd                <- fixd[names(fixd.rearrange)]
  buildtreatmentmodel <- simr::makeLmer(as.formula(formula.model), 
                                        fixef   = unname(fixd), 
                                        VarCorr = vc, 
                                        sigma   = sig, 
                                        data     = data)
  
  return(list("model"          = buildtreatmentmodel, 
              "treatment_term" = treatment_term, 
              "formula.model"  = formula.model))
}

#####################################################    

#' Extends baseline simulated data to be longitudinal
#' 
#' @param data Simulated covariate data (data.frame)
#' @param parameter  Rate of change parameter of interest (character)
#' @param time Vector of time points (numeric)
#' @param rand.effect Column name for subject (character)
#' @return Longitudinal simulated covariate data
#'
ExtendLongitudinal <- function(data, parameter, time, rand.effect) {
    rows              <- nrow(data)
    data              <- data %>% slice(rep(1 : n(), each = length(time)))
    ids               <- rep(1:rows, length(time))
    ids               <- ids[order(ids, decreasing = FALSE)]
    new_time_rows     <- rep(time, rows)
    data[rand.effect] <- factor(ids)
    data[parameter]   <- time
    return(data)
  }

#####################################################      

#' Creates cutoff for discretizing step in multi-variate normal distribution estimation
#' 
#' @param mean Mean value for covariate (numeric)
#' @param sd Standard deviation for covariate (numeric)
#' @param pi Proportion of subjects in the empirical distribution with categorical value Xi(i ≤ k) (list)
#' @return Cutoff to re-categorize 
#'
CalcCutoff <- function(mean, sd, pi) {
    return.vec   <- c()
    for(i in 1:length(pi)) {
      pi.i       <- pi[[i]]
      cutoff     <- mean + (sd * qnorm(pi.i))
      cutoff     <- exp(cutoff)
      return.vec <- append(return.vec, cutoff)
    }
    return.vec <- c(-Inf, return.vec, Inf)
    return(return.vec)
  }

#####################################################      

#' Creates Multivariate Normal Distribution based on pilot data using the Continuous method
#' 
#' @param data Pilot data (data.frame)
#' @param rand.effect Column name for subject (character)
#' @param n Number of subjects to generate (numeric)
#' @param covariates Column names of covariates (character)
#' @param cols.numeric Column names of numeric covariates (character)
#' @param cols.factor Column names of categorical covariates (character)
#' @return Simulated covariate data with defined by pilot data distribution
#'
DefineMVND <- function(data, rand.effect, n, covariates, cols.numeric, cols.factor) {
  
  props.covs <- function(list) {
    list <- as.numeric(list)
    list <- list / sum(list)
    pi   <- c()
    for(i in 1 : length(list)) {
      pi <- append(pi, sum(list[1 : i]))
    }
    pi <- pi[1 : (length(pi) - 1)]
    return(pi)
  }
  
  thresholds      <- c()
  data            <- data[!duplicated(data[rand.effect]),]
  data            <- data[ ,covariates]
  columns.factor  <- cols.factor
  columns.numeric <- cols.numeric
  continuousdata  <- data.frame(matrix(nrow = nrow(data)))
  
  if(length(columns.numeric) > 0) {
    continuousdata  <- cbind(continuousdata, data[ ,columns.numeric])
  }
  
  if(length(columns.factor) > 0) { 
    for(i in 1:length(columns.factor)) {
      name.i   <- columns.factor[i]
      levels.i <- levels(factor(data[ ,name.i]))
      newvec   <- mapvalues(data[ ,name.i], from = levels.i, to = c(1 : length(levels.i)))
      continuousdata[name.i] <- newvec
    }
  }
  
  continuousdata[,1]  <- NULL
  continuousdata      <- continuousdata %>% mutate_all(as.character)
  continuousdata      <- continuousdata %>% mutate_all(as.numeric)
  empir.distr         <- continuousdata
  cov.b4              <- empir.distr
  continuousdata      <- as.matrix(continuousdata)
  continuousdatalog   <- log(continuousdata)
  covmatrix           <- cov(continuousdatalog)
  means               <- colMeans(continuousdatalog)
  simcovs             <- MASS::mvrnorm(n = n, mu = means, Sigma = covmatrix)
  simcovs             <- exp(simcovs)
  cov.after           <- simcovs
  cont.thresh         <- continuousdatalog[,columns.factor]
  means.log           <- colMeans2(cont.thresh)
  sds.log             <- colSds(cont.thresh)
  names(means.log)    <- names(sds.log) <- columns.factor
  pi                  <- Map(table, empir.distr[,columns.factor])
  pi                  <- Map(props.covs, pi)
  simcovs             <- as.data.frame(simcovs)
  if(length(columns.factor) > 0) {
    for(i in columns.factor) {
      threshold   <- CalcCutoff(means.log[[i]], sds.log[[i]], pi[[i]])
      simcovs[,i] <- cut(simcovs[ ,i], 
                         breaks  = threshold, 
                         labels  = 1 : (length(threshold) - 1))
      levels.i    <- levels(factor(data[ ,i]))
      simcovs[,i] <- mapvalues(simcovs[ ,i], from = c(1:length(levels.i)), to = levels.i)
      simcovs[,i] <- factor(as.character(simcovs[ ,i]))
    }
  }
  return(simcovs)
}

#####################################################    

#' Dichotomizes continuous variables for balancing treatment/placebo groups
#' 
#' @param longdata Longitudinal covariate data (data.frame)
#' @param rand.effect Column name for subject (character)
#' @param parameter  Rate of change parameter of interest (character)
#' @param stratcols Column names of continuous covariates to dichotomize (character)
#' @return Longitduinal covariate data with additional dichotomized columns 
#'
StratifyContinuous <- function(longdata, rand.effect, parameter, stratcols) {
  StratifyContVar  <- function(data, stratcols, rand.effect) {
    newdata <- data.frame(matrix(nrow = nrow(data)))
    for(i in stratcols) {
      var <- data[ ,i]
      qt  <- quantile(var)
      groupingvar        <- cut(var, breaks = c(-Inf, unname(qt[3]), Inf), labels = c(0, 1))
      stratname          <- paste(i, "_strat", sep = "")
      newdata[stratname] <- factor(as.character(groupingvar))
    }
    newdata[ ,1] <- NULL
    newdata[rand.effect] <- data[rand.effect]
    return(newdata)
  }
  data       <- longdata[!duplicated(longdata[rand.effect]),]
  bline      <- StratifyContVar(data, stratcols = stratcols, rand.effect = rand.effect)
  longdata   <- merge(longdata, bline, by = rand.effect, all.x = TRUE)
  return("longdata"= longdata)
}

#####################################################    


#' Assigns treatment and placebo groups by balancing covariates
#' 
#' @param longdata Longitudinal covariate data (data.frame)
#' @param rand.effect Column name for subject (character)
#' @param balance.covariates Columns to use for balancing treatment/placebo groups
#' @return A list with longitudinal data with new treatment term, and proportion of covariate values in treatment/placebo group for comparison
#'
RandomizeTreatment <- function(longdata, rand.effect, balance.covariates) {
  `%notin%`             <- Negate(`%in%`)
  data                  <- longdata[!duplicated(longdata[rand.effect]), ]
  data$stratvar         <- interaction(data[ ,balance.covariates])
  data$stratvar         <- factor(data$stratvar)
  treatment             <- as.data.frame(stratified(data, 
                                         "stratvar", size = (.5), 
                                         bothSets = FALSE))
  
  data.id               <- unique(data[ ,rand.effect])
  strat.data.id         <- unique(treatment[ ,rand.effect])
  if(nrow(treatment) == 0) {
    keeprids          <- sample(data.id, round(length(data.id) * .5))
    treatment         <- subset(data, data[ ,rand.effect] %in% keeprids)
  }
 
  rownames(treatment) <- 1:nrow(treatment)
  half                <- round(nrow(data) / 2)
  diff                <- half - nrow(treatment)
  if(diff > 0) {
    disjoint         <-  subset(data, data[ ,rand.effect] %notin% strat.data.id)
    add.rows         <-  sample_n(disjoint, diff)
    treatment        <-  rbind(treatment, add.rows)
    placebo          <-  subset(data, data[ ,rand.effect] %notin% treatment[ ,rand.effect])
  
    } else if(diff < 0)  {
    drop.rows <- sample(1:nrow(treatment), abs(diff))
    treatment           <-  treatment[-drop.rows, ]
    placebo             <-  subset(data, data[ ,rand.effect] %notin% treatment[ ,rand.effect])
  
    } else {
    placebo             <-  subset(data, data[ ,rand.effect] %notin% treatment[ ,rand.effect])
  }
 
  props.treatment       <- CalcProportionPos(treatment, balance.covariates)
  props.placebo         <- CalcProportionPos(placebo, balance.covariates)
  treat.rids            <- unique(treatment[ ,rand.effect])
  control.rids          <- unique(placebo[   ,rand.effect])
  
  treatmentrows         <- subset(longdata, longdata[ ,rand.effect] %in% treat.rids)
  treatmentrows$treat   <- rep(1, nrow(treatmentrows))
  
  controlrows           <- subset(longdata, longdata[ ,rand.effect] %in% control.rids)
  controlrows$treat     <- rep(0, nrow(controlrows))
  
  returndata            <- rbind(treatmentrows, controlrows)
  returndata$treat      <- returndata$treat
  
  return(list("data"  = returndata, 
              "props" = list("Treatment" = props.treatment,
                             "Placebo"   = props.placebo)))
}

#####################################################    

#' Test proportion of covariate value in each treatment group
#' 
#' @param covariate.props A list of proportion of covariate values by treatment group (list)
#' @return A p-value indicating of covariate balance is significantly differnt between groups
#'
PropTestIter <- function(covariate.props) {
  tr       <- covariate.props$Treatment
  pl       <- covariate.props$Placebo
  init.vec <- c()
  for(i in 1:length(tr)) {
    a        <- tr[[i]]
    b        <- pl[[i]]
    mat      <- matrix(c(a, b), ncol = 2)
    p.val    <- chisq.test(mat, simulate.p.value = TRUE)$p.value
    init.vec <- append(init.vec, p.val)
  }
  return(init.vec)
}

#####################################################    

#' Calculates confidence interval for power
#' 
#' @param data A data frame of p.values (data.frame)
#' @param sig_level Significance level of p.values (numeric)
#' @return Dataframe of power and 95% confidence interval of power
#'
GetConfInt <- function(data, sig_level) {
  names_data  <- colnames(data)
  returnframe <- data.frame()
  for(i in 1:ncol(data)) {
    x           <- unname(unlist(data[i]))
    success     <- length(which(x <= sig_level))
    total       <- length(x)
    mean        <- success / total
    binom       <- binom.test(success, total)
    confinter   <- stats :: binom.test(success, total)
    confinter   <- as.numeric(confinter$conf.int)
    row.val     <- c(mean, confinter)
    returnframe <- rbind(returnframe, row.val)
  }
  colnames(returnframe) <- c("mean", "ci.low", "ci.high")
  return(returnframe)
}

#####################################################    

#' Calculates proportion of covariate values
#' 
#' @param data Simulated covariate data
#' @param model.covariates Covariates for balancing
#' @return Proportion of covariate values
#'
CalcProportionPos <- function(data, model.covariates) {
  data <- as.data.frame(data)
  data <- data[ ,model.covariates]
  prop <- Map(function(x) {table(x)}, data)
  return(prop)
}
