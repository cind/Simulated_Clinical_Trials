source("/Users/adamgabriellang/Desktop/clinical_trial_sim/email_setup.R")

args                   <- commandArgs(trailingOnly = TRUE)
fitted_simulation_list <- as.character(args[1])
seq.a                  <- as.numeric(args[2])
seq.b                  <- as.numeric(args[3])
seq.by                 <- as.numeric(args[4])
iters                  <- as.numeric(args[5])
return_file            <- as.character(args[6])

## simulating Type-I Error
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/helper_functions.R")

#read in simulation model
fitted_simulation_list    <- readRDS(fitted_simulation_list)
fitted_simulation_model   <- fitted_simulation_list$largemodel
fitted_simulation_formula <- fitted_simulation_list$formula_largemodel

fitted_simulation_model_small   <- fitted_simulation_list$smallmodel
fitted_simulation_formula_small <- fitted_simulation_list$formula_smallmodel

sim.string.large <- strsplit(fitted_simulation_formula, "~")[[1]][2]
sim.string.large <- paste("sim.out.large", sim.string.large, sep="~")

sim.string.small <- strsplit(fitted_simulation_formula_small, "~")[[1]][2]
sim.string.small <- paste("sim.out.small", sim.string.small, sep="~")

#create new subjects
fitted_simulation_model    <- simr::extend(fitted_simulation_model, along="RID", n=5000)
fitted_simulation_data     <- getData(fitted_simulation_model)
#change treatment effect to be zero
fixef(fitted_simulation_model)["new_time:treat1"] <- 0
levels_extended <- levels(factor(fitted_simulation_data$RID))
# random sample while balancing for covariates
# simulate outcome and refit model and test for significance of treatment effect
returnframe  <- data.frame(matrix(nrow = iters))
sample_sizes <- seq(seq.a, seq.b, by=seq.by)

splitinput <- strsplit(args[1], "/")[[1]]
splitinput <- splitinput[length(splitinput)]

write.email.subject.start <- paste("T1 Error Simulation Start:", splitinput, Sys.time(), sep = " ")
write.email.text.start <- paste("Starting simulation at", Sys.time(), "for", splitinput, sep = " ")

email.start <- WriteEmail(subject = write.email.subject.start,
                          text = write.email.text.start)

gm_send_message(email.start)




for(i in sample_sizes) {
  cat("\n")
  print(paste("--------------------------------------------", i, "--------------------------------------------", sep=""))
  cat("\n")
  p_vec <- c()
  for(j in 1:iters) {
    print(j)
    levels_subset        <- sample(levels_extended, size = i*2)
    data_subset          <- subset(fitted_simulation_data, RID %in% levels_subset)
    data_subset_baseline <- data_subset[!duplicated(data_subset$RID), ]
    treatment.out        <- RandomizeTreatment2(data_subset_baseline, data_subset)
    data_subset_treated  <- treatment.out[["data"]]
    simulated_outcome_large                  <- simulate(fitted_simulation_model, 
                                                         newdata = data_subset_treated, 
                                                         allow.new.levels = TRUE, 
                                                         use.u = FALSE)$sim_1
    simulated_outcome_small                  <- simulate(fitted_simulation_model_small, 
                                                         newdata = data_subset_treated, 
                                                         allow.new.levels = TRUE, 
                                                         use.u = FALSE)$sim_1
    
    data_subset_treated$sim.out.large <- simulated_outcome_large
    data_subset_treated$sim.out.small <- simulated_outcome_small
    
    model.refit.large <- lme4::lmer(formula = as.formula(sim.string.large), data = data_subset_treated, REML = TRUE)
    model.refit.small <- lme4::lmer(formula = as.formula(sim.string.small), data = data_subset_treated, REML = TRUE)
    iter_kr_ftest                 <- pbkrtest::KRmodcomp(model.refit.large, model.refit.small) 
    pval                          <- getKR(iter_kr_ftest, "p.value")
    p_vec                         <- append(p_vec, pval)
  }
  returnframe <- cbind(returnframe, p_vec)
}

returnframe[,1] <- NULL
colnames(returnframe) <- paste("SS_", sample_sizes, sep="")
saveRDS(returnframe, return_file)

write.email.subject.end <- paste(" T1 Error Simulation Finished:", splitinput, Sys.time(), sep = " ")
write.email.text.end <- paste("Finished simulation at", Sys.time(), "for", splitinput, sep = " ")

email.end <- WriteEmail(subject = write.email.subject.end,
                        text = write.email.text.end)

gm_send_message(email.end)














####################################################################################################################################################################################################################
if(FALSE) {
type1error.frame <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/t1_error.rds")
type1error_rate <- Map(function(x){(length(which(x <= 0.05)) / length(x)) *100 }, type1error.frame)
type1error_rate <- as.numeric(type1error_rate)
error.rate <- data.frame("error" = type1error_rate,
                         "sample_size" = sample_sizes)

ggplot(error.rate, aes(x=sample_size, y=error)) + geom_point(shape="triangle") + ylab("Type 1 Error (%)") + xlab("Sample Size") + scale_x_discrete(limits=seq(20, 300, by=20))
                                                                                                                                                   
###########################


t1error <- cbind(t1errorearlyadadas_20_100, t1errorearlyadadas_120_200, t1errorearlyadadas_220_260, t1errorearlyadadas_280_300)
t1error <- Map(function(x) { length(which(x <= 0.05)) / length(x)}, t1error)
t1error <- as.numeric(t1error)
t1error <- data.frame("error" = t1error,
                      "sample_size" = seq(20, 300, by = 20))


errorgg <- ggplot(t1error, aes(x=sample_size, y = error * 100)) + geom_point(shape="triangle") + ylab("Type 1 Error (%)") 
errorgg <- errorgg + scale_x_discrete(limits = seq(20, 300, by=20), name= "Sample Size") 
}
#########################################################################################################################################################################################