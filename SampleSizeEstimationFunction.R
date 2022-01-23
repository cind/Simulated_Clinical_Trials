ManualSimulation <- function(model, parameter, pct.change = NULL, delta = NULL, time = c(0, .5, 1, 1.5, 2), 
                             data, sample_sizes, nsim, trial_duration, sig_level = .05,
                             verbose = TRUE, balance.covariates = NULL) {
  t1 <- Sys.time()
  cat("Beginning simulation")
  cat("\n")
  init_iter_list                      <-  list()
  init_significance_list              <-  list()
  init_props_list_treatment           <-  list()
  init_props_list_placebo             <-  list()
  formula_model                       <-  as.character(formula(model))
  formula_model_join                  <-  paste(formula_model[2], formula_model[3], sep="~")                         
  iter_form_lm                        <-  paste("model_response", formula_model[3], sep="~")
  model.covariates                    <-  GetCovariates(model, parameter)
  rand.effect                         <-  names(ranef(model))
  if(!is.null(balance.covariates)) {
    model.covariates <- c(model.covariates, balance.covariates)
    model.covariates <- unique(model.covariates)
  }
  
  column.types               <- data[ ,model.covariates]
  cols.numeric               <- c(names(Filter(is.integer, column.types)), 
                                  names(Filter(is.numeric, column.types)))
  cols.numeric               <- cols.numeric[cols.numeric != parameter]
  cols.numeric.stratified    <- paste(cols.numeric, "_strat", sep="")
  cols.factor                <- c(names(Filter(is.factor, column.types)))
  cols.factor                <- cols.factor[cols.factor != rand.effect]
  cols.balance               <- c(cols.numeric.stratified, cols.factor)

  props <- list()  
  pvals <- list()
  #define inner loop for parallelization
  .nsiminnerloop <- function(i, j) {
    
    sim.covariates      <- DefineMVND(  data         = data,
                                        n            = i * 2,
                                        rand.effect  = rand.effect,
                                        covariates   = model.covariates,
                                        cols.numeric = cols.numeric,
                                        cols.factor  = cols.factor)
    sim.covariates.long    <- ExtendLongitudinal(sim.covariates, parameter, time, rand.effect)
    sim.covariates.long    <- StratifyContinuous(sim.covariates.long, rand.effect, parameter)
    treatment.out          <- RandomizeTreatment(sim.covariates.long, rand.effect, cols.balance)
    prop                         <-   treatment.out[["props"]]
    data_sample_treated          <-   treatment.out[["data"]]
    props.test                   <-   PropTestIter(prop)
    levels.factors               <-   any(Map(function(x){nlevels(x)}, 
                                        data_sample_treated[ ,cols.balance]) < 2)
    while(any(props.test <= .05) | levels.factors) {
      
      sim.covariates      <- DefineMVND(  data         = data,
                                          n            = i * 2,
                                          rand.effect  = rand.effect,
                                          covariates   = model.covariates,
                                          cols.numeric = cols.numeric,
                                          cols.factor  = cols.factor)
      sim.covariates.long    <- ExtendLongitudinal(sim.covariates, parameter, time, rand.effect)
      sim.covariates.long    <- StratifyContinuous(sim.covariates.long, rand.effect, parameter)
      treatment.out          <- RandomizeTreatment(sim.covariates.long, rand.effect, cols.balance)
      prop                         <-   treatment.out[["props"]]
      data_sample_treated          <-   treatment.out[["data"]]
      props.test                   <-   PropTestIter(prop)
      levels.factors               <-   any(Map(function(x){nlevels(x)}, 
                                                data_sample_treated[ ,cols.balance]) < 2)
    }
    
    simulate_response_largemodel  <- simulate(model, newdata = data_sample_treated, allow.new.levels = TRUE, use.u = FALSE)
    refit_data_outcomes           <- data.frame("model_response" = simulate_response_largemodel)
    colnames(refit_data_outcomes) <- c("model_response")
    fit_iter_data                 <- bind_cols(refit_data_outcomes, data_sample_treated)
    refit_large                   <- lme4::lmer(formula = as.formula(iter_form_lm), 
                                                data = fit_iter_data, REML = TRUE)

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
  t2      <- Sys.time()
  timerun <- difftime(t2, t1, units = "mins")
  return(list("Power_Per_Sample"          = conf.inter, 
              "Treatment/Placebo_Balance" = props,
              "Run_Time"                  = timerun,
              "pval.df"                   = pval.df))
}
