#run imaging harmonization
# files saved: harmed_and_unharmed_imaging (csv)
#                saved to data_processed
#                 harmonized imaging data
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/harmonize_imaging_for_simulation.R")

#run data preprocessing
#files saved: all_outcomes_merged_data (csv)
#             saved to data_processed
#             all outcomes/demographics/necessary info for simulation
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/processing_raw_data.R")

#run enrichment processing
#files saved:

source("/Users/adamgabriellang/Desktop/clinical_trial_sim/processing_merged_data_for_simulation.R")

#run simulation pre-processing for early AD cohorts
#files saved:
source("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyadsim.R")




