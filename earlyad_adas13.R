modeling.list <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling_minus_control.rds")
treatment.enrichment.list <- modeling.list$treatment.enrichment.list
simulation.formulas <- modeling.list$simulation.formulas
simulation.model.list <- modeling.list$simulation.model.list

earlyad.adas13.try2 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$earlyad$ADAS13, 
                                        formula = simulation.formulas$ADAS13, 
                                        compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                        breaks= seq(100, 1000, by=100), 
                                        yaxislab_dpm="ADAS13",
                                        model = simulation.model.list$earlyad$ADAS13,
                                        return_dpm = FALSE)

earlyad.tplus.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$earlyad.tplus$ADAS13, 
                                              formula = simulation.formulas$ADAS13, 
                                              compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                              breaks= seq(100, 1000, by=100), 
                                              yaxislab_dpm="ADAS13",
                                              model = simulation.model.list$earlyad.tplus$ADAS13,
                                              return_dpm = FALSE)

earlyad.earlyage.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$earlyad.earlyage$ADAS13, 
                                                 formula = simulation.formulas$ADAS13, 
                                                 compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                                 breaks= seq(100, 1000, by=100), 
                                                 yaxislab_dpm="ADAS13",
                                                 model = simulation.model.list$earlyad.earlyage$ADAS13,
                                                 return_dpm = FALSE)




earlyad.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$earlyad$ADAS13, 
                                        formula = simulation.formulas$ADAS13, 
                                        compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                        breaks= seq(100, 1000, by=100), 
                                        yaxislab_dpm="ADAS13",
                                        model = simulation.model.list$earlyad$ADAS13,
                                        return_dpm = TRUE)


earlyad.tplus.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$earlyad.tplus$ADAS13, 
                                              formula = simulation.formulas$ADAS13, 
                                              compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                              breaks= seq(100, 1000, by=100), 
                                              yaxislab_dpm="ADAS13",
                                              model = simulation.model.list$earlyad.tplus$ADAS13,
                                              return_dpm = TRUE)

earlyad.earlyage.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$earlyad.earlyage$ADAS13, 
                                                 formula = simulation.formulas$ADAS13, 
                                                 compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                                 breaks= seq(100, 1000, by=100), 
                                                 yaxislab_dpm="ADAS13",
                                                 model = simulation.model.list$earlyad.earlyage$ADAS13,
                                                 return_dpm = TRUE)





early.ad.list <- list("AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n" = earlyad.adas13,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n" = earlyad.tplus.adas13,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"  = earlyad.earlyage.adas13)


early.ad.list.dpm <- list("AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30\n" = earlyad.adas13.dpm,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 65-85  MMSE 20-30 Tau+ (PTAU)\n" = earlyad.tplus.adas13.dpm,
                      "AD or MCI  \nA+ CDR = .5 or =1 \nAge 55-65  MMSE 20-30\n"  = earlyad.earlyage.adas13.dpm)


early.ad.combine.dpm <- GroupDiseaseTraj(early.ad.list.dpm, yaxislab_dpm = "ADAS13 (Z-score)")
early.ad.combine.dpm$plot <- early.ad.combine.dpm$plot + labs(title = "Early AD")
early.ad.combine.dpm

saveRDS(early.ad.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad_adas13_curves.rds")


early.ad.list <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyad_adas13_curves.rds")

combine.early.ad.adas13 <- CombineSimPlots(early.ad.list, seq(100, 1000, by=100))
combine.early.ad.adas13.plt <- combine.early.ad.adas13$plot
combine.early.ad.adas13.plt <- combine.early.ad.adas13.plt + ylim(0, 100) + labs(title = "ADAS-13  Early AD") +geom_hline(yintercept = 80, linetype="dashed")
combine.early.ad.adas13.plt
ggsave(filename="/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyaf_adas13_curves.jpeg", plot= combine.early.ad.adas13.plt, width = 20, height = 20, units="cm", dpi = "retina", device = "tiff")
ggsave(filename="/Users/adamgabriellang/Desktop/clinical_trial_sim/earlyaf_adas13_dpm.jpeg", plot= early.ad.combine.dpm$plot, width = 20, height = 20, units="cm", dpi = "retina", device = "tiff")
