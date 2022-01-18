ManualSimulation <- function(formula_model, model, treatment_term, sample_sizes, nsim, data, trial_duration, t1errorsim = "power_simulation", sig_level=.05, verbose=TRUE) {
  t1 <- Sys.time()
  cat("Beginning simulation")
  cat("\n")
  init_iter_list                      <-  list()
  init_significance_list              <-  list()
  init_props_list_treatment           <-  list()
  init_props_list_placebo             <-  list()
  form_lm_split <- strsplit(formula_model, "~")[[1]][2]
  iter_form_lm  <- paste("model_response", form_lm_split, sep="~")
  if(t1errorsim=="T1") {
    fixef(largemodel)[treatment_term] <- 0
  }
  
  props <- list()  
  pvals <- list()

  #define inner loop for parallelization
  .nsiminnerloop <- function(i, j) {
    
    sim.covariates      <- DefineMVND(data = data,
                                      n    = i*2)
    sim.covariates.long <- ExtendLongitudinal(sim.covariates, trial_duration = trial_duration)
    sim.covariates.long <- StratifyContinuous(sim.covariates.long, c("AGE_bl", "PTEDUCAT_bl", "MMSE_bl"))
    #split into treatment and placebo groups while balancing subjects
    treatment.out                <-  RandomizeTreatment2(sim.covariates.long)
    prop                         <-  treatment.out[["props"]]
    data_sample_treated          <-  treatment.out[["data"]]
    props.test                   <-   PropTestIter(prop)
    while(any(props.test <= .05) | nlevels(factor(sim.covariates$PTGENDER)) < 2 | nlevels(factor(sim.covariates$CDGLOBAL_bl)) < 2) {
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
    simulate_response_largemodel  <-  simulate(model, newdata = data_sample_treated, allow.new.levels=TRUE, use.u=FALSE)
    refit_data_outcomes           <- data.frame("model_response" = simulate_response_largemodel)
    colnames(refit_data_outcomes) <- c("model_response")
    fit_iter_data                 <- bind_cols(refit_data_outcomes, data_sample_treated)
    refit_large                   <- lme4::lmer(formula = as.formula(iter_form_lm), 
                                                data = fit_iter_data, REML = TRUE, control = lme4::lmerControl(optimizer = "nmkbw"))

    pval <-  as.numeric(summary(lmerTest::as_lmerModLmerTest(refit_large))[["coefficients"]][,"Pr(>|t|)"][treatment_term])
    
    if(verbose) {
    cat("\r", j, " out of ", nsim, " complete", sep = "")
    }
    return(list("pval"      = pval,
                "Treatment" = prop[["Treatment"]],
                "Placebo"   = prop[["Placebo"]]))
  }
  
  for(i in 1:length(sample_sizes)) {
    props_ss <- list()
    pvals_ss  <- c()
    for(j in 1:nsim) {
      
      iter_j <- .nsiminnerloop(sample_sizes[i], j)
      props_ss[[j]] <- iter_j[c("Treatment", "Placebo")]
      pvals_ss[[j]] <- iter_j[[c("pval")]]
      }
    props[[i]] <- props_ss
    pvals[[i]] <- pvals_ss
  }
  names(props) <- names(pvals) <- paste("sample_size_", sample_sizes, sep="")
  pval.df                      <- as.data.frame(do.call(cbind, pvals))
  conf.inter                   <- GetConfInt(pval.df, sig_level)
  t2     <- Sys.time()
  timerun <- difftime(t2, t1, units = "mins")
  return(list("Power_Per_Sample"   = conf.inter, 
              "Treatment/Placebo_Balance" = props,
              "Run_Time" = timerun,
              "simulation_type" = t1errorsim,
              "pval.df" = pval.df))
}
