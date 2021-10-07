modeling.list <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/list_for_modeling.rds")
treatment.enrichment.list <- modeling.list$treatment.enrichment.list
simulation.formulas <- modeling.list$simulation.formulas
simulation.model.list <- modeling.list$simulation.model.list
simrOptions(nsim=1000)
if(FALSE) {
ad.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad$ADAS13, 
                                   formula = simulation.formulas$ADAS13, 
                                   compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                   breaks= seq(100, 1000, by=100), 
                                   yaxislab_dpm="ADAS13",
                                   model = simulation.model.list$ad$ADAS13,
                                   return_dpm = FALSE)
ad.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad$ADAS13, 
                                formula = simulation.formulas$ADAS13, 
                                compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                breaks= seq(100, 1000, by=100), 
                                yaxislab_dpm="ADAS13",
                                model = simulation.model.list$ad$ADAS13,
                                return_dpm = TRUE)

ad.tplus.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad.tplus$ADAS13, 
                                         formula = simulation.formulas$ADAS13, 
                                         compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                         breaks= seq(100, 1000, by=100), 
                                         yaxislab_dpm="ADAS13",
                                         model = simulation.model.list$ad.tplus$ADAS13,
                                         return_dpm = FALSE)
ad.tplus.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad.tplus$ADAS13, 
                                         formula = simulation.formulas$ADAS13, 
                                         compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                         breaks= seq(100, 1000, by=100), 
                                         yaxislab_dpm="ADAS13",
                                         model = simulation.model.list$ad.tplus$ADAS13,
                                         return_dpm = TRUE)


}
ad.earlyage.adas13 <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad.earlyage$ADAS13, 
                                            formula = simulation.formulas$ADAS13, 
                                            compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                            breaks= seq(100, 1000, by=100), 
                                            yaxislab_dpm="ADAS13",
                                            model = simulation.model.list$ad.earlyage$ADAS13,
                                            return_dpm = FALSE)

ad.earlyage.adas13.dpm <- SampleSizeSimulation2(sim.data = treatment.enrichment.list$ad.earlyage$ADAS13, 
                                            formula = simulation.formulas$ADAS13, 
                                            compare_str = "TOTAL13_zscore ~ new_time + PTEDUCAT_bl + AGE_bl + PTGENDER_bl + fulllewy*new_time + fullcaa*new_time + fulltdp43*new_time + (1|RID)", 
                                            breaks= seq(100, 1000, by=100), 
                                            yaxislab_dpm="ADAS13",
                                            model = simulation.model.list$ad.earlyage$ADAS13,
                                            return_dpm = TRUE)






ad.list <- list("AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n" = ad.adas13,
                "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n" = ad.tplus.adas13,
                "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n"  = ad.earlyage.adas13)

ad.list.dpm <- list("AD A+ CDR >= 1 \nAge 65-85  MMSE 20-28 \n" = ad.adas13.dpm,
                "AD A+ CDR >=1  \nAge 65-85  MMSE 20-28 Tau+ (PTAU)\n" = ad.tplus.adas13.dpm,
                "AD A+ CDR >= 1 \nAge 55-65  MMSE 20-28\n"  = ad.earlyage.adas13.dpm)



combine.ad.dpm <- GroupDiseaseTraj(ad.list.dpm, yaxislab_dpm = "ADAS13 (Z-score)")
combine.ad.dpm$plot <- combine.ad.dpm$plot + labs(title = "AD")
combine.ad.dpm 



saveRDS(ad.list, "/Users/adamgabriellang/Desktop/clinical_trial_sim/ad_adas13_curves.rds")

ad.list <- readRDS("/Users/adamgabriellang/Desktop/clinical_trial_sim/ad_adas13_curves.rds")

combine.ad.adas13       <- CombineSimPlots(ad.list, seq(100, 1000, by=100))
combine.ad.adas13.plt   <- combine.ad.adas13$plot
combine.ad.adas13.plt   <- combine.ad.adas13.plt + ylim(0, 100) + labs(title = "ADAS-13 AD")+geom_hline(yintercept = 80, linetype="dashed")
combine.ad.adas13.plt
ggsave(filename="/Users/adamgabriellang/Desktop/clinical_trial_sim/ad_adas13_curves.jpeg", plot= combine.ad.adas13.plt, width = 20, height = 20, units="cm", dpi = "retina", device = "tiff")
ggsave(filename="/Users/adamgabriellang/Desktop/clinical_trial_sim/ad_adas13_dpm.jpeg", plot= combine.ad.dpm$plot, width = 20, height = 20, units="cm", dpi = "retina", device = "tiff")
