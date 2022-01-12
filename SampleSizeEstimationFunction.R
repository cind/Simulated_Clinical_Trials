ManualSimulation <- function(formula_largemodel, largemodel, formula_smallmodel, smallmodel, sample_sizes, nsim, data, trial_duration=NULL, t1errorsim="SIMU") {
  # load in functions from GlobalEnv into current enviornment for clusterExport
  force(RandomizeTreatment2)
  force(DefineMVND)
  force(ExtendLongitudinal)
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
  force(plyr::mapvalues)  
  force(cramer::cramer.test)
  force(dplyr::mutate_all)
  force(MASS::mvrnorm)
  force(matrixStats::colMeans2)
  force(matrixStats::colSds)
  force(dplyr::slice)
  force(dplyr::n)
  
  opts <- list(chunkSize=10)
  #################
  cat("Beginning simulation")
  cat("\n")
  init_iter_list                      <-  list()
  init_significance_list              <-  list()
  init_props_list_treatment           <-  list()
  init_props_list_placebo             <-  list()
  form_lm_split <- strsplit(formula_largemodel, "~")[[1]][2]
  iter_form_lm  <- paste("large_model_response", form_lm_split, sep="~")
  if(t1errorsim=="T1") {
    fixef(largemodel)["new_time:treat1"] <- 0
  }
  #define inner loop for parallelization
  .nsiminnerloop <- function(j) {
    sim.covariates      <- DefineMVND(data = data,
                                      n    = i*2)
    sim.covariates.long <- ExtendLongitudinal(sim.covariates, trial_duration = trial_duration)
    sim.covariates.long <- StratifyContinuous(sim.covariates.long, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
    #split into treatment and placebo groups while balancing subjects
    treatment.out                <-  RandomizeTreatment2(sim.covariates.long)
    prop                         <-  treatment.out[["props"]]
    data_sample_treated          <-  treatment.out[["data"]]
    props.test                   <-   PropTestIter(prop)
    while(any(props.test <= .05)) {
      sim.covariates      <- DefineMVND(data = data,
                                        n    = i*2)
      sim.covariates.long <- ExtendLongitudinal(sim.covariates, trial_duration = trial_duration)
      sim.covariates.long <- StratifyContinuous(sim.covariates.long, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
      
      #split into treatment and placebo groups while balancing subjects
      treatment.out                <-  RandomizeTreatment2(sim.covariates.long)
      prop                         <-  treatment.out[["props"]]
      data_sample_treated          <-  treatment.out[["data"]]
      props.test                   <-  PropTestIter(prop)
    }
    
    simulate_response_largemodel <-  simulate(largemodel, newdata = data_sample_treated, allow.new.levels=TRUE, use.u=FALSE)
    refit_data_outcomes           <- data.frame("large_model_response" = simulate_response_largemodel)
    colnames(refit_data_outcomes) <- c("large_model_response")
    fit_iter_data                 <- bind_cols(refit_data_outcomes, data_sample_treated)
    refit_large                   <- lme4::lmer(formula = as.formula(iter_form_lm), 
                                                data = fit_iter_data, REML = TRUE, control = lme4::lmerControl(optimizer = "nmkbw"))

    pval <-  as.numeric(summary(lmerTest::as_lmerModLmerTest(refit_large))[["coefficients"]][,"Pr(>|t|)"]["new_time:treat1"])
    
    cat("\r", j, " out of ", nsim, " complete", sep = "")
    
    return(list("pval"      = pval,
                "Treatment" = prop[["Treatment"]],
                "Placebo"   = prop[["Placebo"]]))
  }
  pb <- txtProgressBar(min = 1, max=nsim, style = 1)
  envlist    <- mget(ls())
  envlist    <- names(envlist)
  env.append <- c("RandomizeTreatment2","KRmodcomp",
                  "getKR","lmer","bind_cols","stratified", "sample_n",
                  ".nsiminnerloop", "mapvalues", "%>%",
                  "CalcProportionPos", 
                  "%notin%", "ExtendLongitudinal", "StratifyContinuous", "DefineMVND", "PropTestIter",
                 "cramer.test", "mutate_all", "mvrnorm", "colMeans2", "colSds", "slice", "n")
  
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
  timerun <- difftime(t2, t1, units = "mins")
  return(list("Successes"   = successes, 
              "power.data"  = conf.inter, 
              "cov.balance" = outer,
              "prop.tests"  = balance.cov.diagnostics,
              "time_to_run" = timerun,
              "ftestdf" = ftestdf,
              "simulation_type" = t1errorsim))
}
