modeling.list <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling.rds")
treatment.enrichment.list <- modeling.list$treatment.enrichment.list
simulation.formulas <- modeling.list$simulation.formulas
simulation.model.list <- modeling.list$simulation.model.list
simrOptions(nsim=1000)


mci.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$mci$ADAS13, 
                                    formula = simulation.formulas$ADAS13, 
                                    compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                    breaks= seq(100, 1000, by=100), 
                                    yaxislab_dpm="ADAS13",
                                    model = simulation.model.list$mci$ADAS13,
                                    return_dpm = FALSE)

mci.tplus.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$mci.tplus$ADAS13, 
                                          formula = simulation.formulas$ADAS13, 
                                          compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                          breaks= seq(100, 1000, by=100), 
                                          yaxislab_dpm="ADAS13",
                                          model = simulation.model.list$mci.tplus$ADAS13,
                                          return_dpm = FALSE)

mci.earlyage.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$mci.earlyage$ADAS13, 
                                             formula = simulation.formulas$ADAS13, 
                                             compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                             breaks= seq(100, 1000, by=100), 
                                             yaxislab_dpm="ADAS13",
                                             model = simulation.model.list$mci.earlyage$ADAS13,
                                             return_dpm = FALSE)



mci.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$mci$ADAS13, 
                                    formula = simulation.formulas$ADAS13, 
                                    compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                    breaks= seq(100, 1000, by=100), 
                                    yaxislab_dpm="ADAS13",
                                    model = simulation.model.list$mci$ADAS13,
                                    return_dpm = TRUE)

mci.tplus.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$mci.tplus$ADAS13, 
                                          formula = simulation.formulas$ADAS13, 
                                          compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                          breaks= seq(100, 1000, by=100), 
                                          yaxislab_dpm="ADAS13",
                                          model = simulation.model.list$mci.tplus$ADAS13,
                                          return_dpm = TRUE)

mci.earlyage.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$mci.earlyage$ADAS13, 
                                             formula = simulation.formulas$ADAS13, 
                                             compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                             breaks= seq(100, 1000, by=100), 
                                             yaxislab_dpm="ADAS13",
                                             model = simulation.model.list$mci.earlyage$ADAS13,
                                             return_dpm = TRUE)



mci.list <- list("MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 \n" = mci.adas13,
                       "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n" = mci.tplus.adas13,
                      "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n"  = mci.earlyage.adas13)


mci.list.dpm <- list("MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 \n" = mci.adas13.dpm,
                 "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n" = mci.tplus.adas13.dpm,
                 "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n"  = mci.earlyage.adas13.dpm)


combine.mci.dpm <- GroupDiseaseTraj(mci.list.dpm, yaxislab_dpm = "ADAS13 (Z-score)")
combine.mci.dpm$plot <- combine.mci.dpm$plot + labs(title = "MCI")
combine.mci.dpm




saveRDS(mci.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/mci_adas13_curves.rds")


mci.list <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/mci_adas13_curves.rds")
names(mci.list) <- c("MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 \n",
                     "MCI A+ CDR = .5 \nAge 65-85  MMSE 24-30 Tau+ (PTAU)\n",
                     "MCI A+ CDR = .5 \nAge 55-65  MMSE 24-30\n")


groupsimmci <- CombineSimPlots(mci.list, seq(100,1000, by=100))
groupsimmci.plt <- groupsimmci$plot
groupsimmci.plt <- groupsimmci.plt + ylim(0, 100) + labs(title = "ADAS-13  MCI") +geom_hline(yintercept = 80, linetype="dashed")
groupsimmci.plt
ggsave(filename="/Users/adamgabriellang/Desktop/clinical_trial_sim/mci_adas13_curves.jpeg", plot= groupsimmci.plt, width = 20, height = 20, units="cm", dpi = "retina", device = "tiff")
ggsave(filename="/Users/adamgabriellang/Desktop/clinical_trial_sim/mci_adas13_dpm.jpeg", plot= combine.mci.dpm$plot, width = 20, height = 20, units="cm", dpi = "retina", device = "tiff")

