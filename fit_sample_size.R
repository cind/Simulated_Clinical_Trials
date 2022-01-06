args                   <- commandArgs(trailingOnly = TRUE)
fitted_simulation_list <- as.character(args[1])
seq.a                  <- as.numeric(args[2])
seq.b                  <- as.numeric(args[3])
seq.by                 <- as.numeric(args[4])
nsim                   <- as.numeric(args[5])
trial_duration         <- as.numeric(args[6])
return_file            <- as.character(args[7])
t1                     <- as.character(args[8])

if(t1 != "TRUE") {
  t1 <- NULL
}

source("/Users/adamgabriellang/Desktop/clinical_trial_sim/helper_functions.R")
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/email_setup.R")


#splitinput <- strsplit(fitted_simulation_list, "/")[[1]]
#splitinput <- splitinput[length(splitinput)]

#write.email.subject.start <- paste("Simulation Start:", splitinput, Sys.time(), sep = " ")
#write.email.text.start <- paste("Starting simulation at", Sys.time(), "for", splitinput, sep = " ")

#email.start <- WriteEmail(subject = write.email.subject.start,
#                          text = write.email.text.start)

#gm_send_message(email.start)
models_list          <- readRDS(fitted_simulation_list)
data                 <- simr::getData(models_list$smallmodel)
data$PTGENDER        <- relevel(data$PTGENDER, "Male")
data                 <- data[!duplicated(data$RID),]
ss_fitted <- ManualSimulation(formula_largemodel    = models_list$formula_largemodel,
                                 largemodel         = models_list$largemodel,
                                 formula_smallmodel = models_list$formula_smallmodel,
                                 smallmodel         = models_list$smallmodel,
                                 sample_sizes       = seq(seq.a, seq.b, by=seq.by),
                                 nsim               = nsim,
                                 data               = data,
                                 trial_duration     = trial_duration,
                                 t1errorsim         = t1)

ss_fitted[["distr_comparison"]] <- data_simulated
ss_fitted[["t1"]] <- t1
ss_fitted[["args"]] <- args
saveRDS(ss_fitted, return_file)


#write.email.subject.end <- paste("Simulation Finished:", splitinput, Sys.time(), sep = " ")
#write.email.text.end <- paste("Finished simulation at", Sys.time(), "for", splitinput, sep = " ")

#email.end <- WriteEmail(subject = write.email.subject.end,
#                          text = write.email.text.end)

#gm_send_message(email.end)

