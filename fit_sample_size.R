args                   <- commandArgs(trailingOnly = TRUE)
fitted_simulation_list <- as.character(args[1])
seq.a                  <- as.numeric(args[2])
seq.b                  <- as.numeric(args[3])
seq.by                 <- as.numeric(args[4])
nsim                   <- as.numeric(args[5])
trial_duration         <- as.numeric(args[6])
return_file            <- as.character(args[7])
t1                     <- as.character(args[8])

source("/Users/adamgabriellang/Desktop/clinical_trial_sim/helper_functions.R")
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/email_setup.R")
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/SampleSizeEstimationFunction.R")

models_list          <- readRDS(fitted_simulation_list)
data                 <- simr::getData(models_list$smallmodel)
data$PTGENDER        <- relevel(data$PTGENDER, "Male")
data                 <- data[!duplicated(data$RID),]
treatment_term       <- "new_time:treat1"

ss_fitted <- ManualSimulation(formula_model    = models_list$formula_largemodel,
                                 model         = models_list$largemodel,
                                 treatment_term = treatment_term,
                                 sample_sizes       = seq(seq.a, seq.b, by=seq.by),
                                 nsim               = nsim,
                                 data               = data,
                                 sig_level          = 0.05,
                                 trial_duration     = trial_duration,
                                 verbose            = TRUE)
if(t1=="T1") {
  split.return <- strsplit(return_file, "/")
  obj.name <- split.return[[1]][length(split.return[[1]])]
  obj<-strsplit(obj.name, ".rds")[[1]]
  obj <- paste(obj, "_t1", sep="")
  obj <- paste(obj, ".rds", sep="")
  split.return[[1]][length(split.return[[1]])] <- obj
  look <-unlist(split.return)
  return_file <- paste(look, collapse = "/")
}


ss_fitted[["args"]] <- args

saveRDS(ss_fitted, return_file)


