## simulating Type-I Error
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/helper_functions.R")
#read in simulation model
fitted_simulation_list    <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/adas13_power_data/earlyadadas13.rds")
fitted_simulation_model   <- fitted_simulation_list$largemodel
fitted_simulation_formula <- fitted_simulation_list$formula_largemodel
#create new subjects
fitted_simulation_model    <- simr::extend(fitted_simulation_model, along="RID", n=5000)
fitted_simulation_data     <- getData(fitted_simulation_model)
#change treatment effect to be zero
fixef(fitted_simulation_model)["new_time:treat1"] <- 0
levels_extended <- levels(factor(fitted_simulation_data$RID))
# random sample while balancing for covariates
# simulate outcome and refit model and test for significance of treatment effect
iters <- 1000
returnframe  <- data.frame(matrix(nrow = iters))
sample_sizes <- seq(20, 300, by=20)
for(i in sample_sizes) {
  cat("\n")
  print(paste("--------------------------------------------", i, "", sep="--------------------------------------------"))
  cat("\n")
  p_vec <- c()
  for(j in 1:iters) {
    print(j)
    levels_subset        <- sample(levels_extended, size = i*2)
    data_subset          <- subset(fitted_simulation_data, RID %in% levels_subset)
    data_subset_baseline <- data_subset[!duplicated(data_subset$RID), ]
    treatment.out        <- RandomizeTreatment2(data_subset_baseline, data_subset)
    data_subset_treated  <- treatment.out[["data"]]
    simulated_outcome                        <- simulate(fitted_simulation_model, 
                                                         newdata = data_subset_treated, 
                                                         allow.new.levels = TRUE, 
                                                         use.u = FALSE)$sim_1
    data_subset_treated$sim.out <- simulated_outcome
    model.refit <- lmerTest::lmer(formula = as.formula(fitted_simulation_formula), data = data_subset_treated, REML=TRUE)
    sum_fit     <- summary(model.refit)$coefficients[,5]
    treat_pval  <- sum_fit[length(sum_fit)]
    p_vec       <- append(p_vec, treat_pval)
  }
  returnframe <- cbind(returnframe, p_vec)
}

returnframe[,1] <- NULL
colnames(returnframe) <- paste("SS_", sample_sizes, sep="")
saveRDS(returnframe, "/Users/adamgabriellang/Desktop/clinical_trial_sim/t1_error.rds")


####################################################################################################################################################################################################################
if(FALSE) {
type1error.frame <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/t1_error.rds")
type1error_rate <- Map(function(x){(length(which(x <= 0.05)) / length(x)) *100 }, type1error.frame)
type1error_rate <- as.numeric(type1error_rate)
error.rate <- data.frame("error" = type1error_rate,
                         "sample_size" = sample_sizes)

ggplot(error.rate, aes(x=sample_size, y=error)) + geom_point(shape="triangle") + ylab("Type 1 Error (%)") + xlab("Sample Size") + scale_x_discrete(limits=seq(20, 300, by=20))
}                                                                                                                                                   
####################################################################################################################################################################################################################