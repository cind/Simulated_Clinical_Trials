modeling.list <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling.rds")
treatment.enrichment.list <- modeling.list$treatment.enrichment.list
simulation.formulas <- modeling.list$simulation.formulas
simulation.model.list <- modeling.list$simulation.model.list

ad.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad$ADAS13, 
                                   formula = simulation.formulas$ADAS13, 
                                   fcompare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                   breaks= seq(100, 1000, by=100), 
                                   yaxislab_dpm="ADAS13",
                                   model = simulation.model.list$ad$ADAS13,
                                   return_dpm = FALSE)

ad.tplus.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad.tplus$ADAS13, 
                                         formula = simulation.formulas$ADAS13, 
                                         fcompare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                         breaks= seq(100, 1000, by=100), 
                                         yaxislab_dpm="ADAS13",
                                         model = simulation.model.list$ad.tplus$ADAS13,
                                         return_dpm = FALSE)

ad.earlyage.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad.earlyage$ADAS13, 
                                            formula = simulation.formulas$ADAS13, 
                                            fcompare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                            breaks= seq(100, 1000, by=100), 
                                            yaxislab_dpm="ADAS13",
                                            model = simulation.model.list$ad.earlyage$ADAS13,
                                            return_dpm = FALSE)


