new_merged_data <- read.csv("/Users/adamgabriellang/Desktop/clinical_trial_sim/Data/final_merged_data_4sim_remerged.csv")
aut.rows <- which(complete.cases(new_merged_data[,c("NPBRAAK", "NPNEUR", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE",
                                                     "NPLBOD", "NPAMY")]))

imputed.rows <- which(complete.cases(new_merged_data[,c("TDP43", "LEWY", "CAA")]))
no.non.ad <- new_merged_data[-aut.rows,]
no.non.ad <- new_merged_data[-imputed.rows,]
no.non.ad.imaging <- no.non.ad[-which(is.na(no.non.ad$ST103CV)),]

no.non.ad.imaging$hipp_average <- (no.non.ad.imaging$ST29SV + no.non.ad.imaging$ST88SV) /2
no.non.ad.imaging$hipp_average_harmonized <- (no.non.ad.imaging$ST29SV_harmonized + no.non.ad.imaging$ST88SV_harmonized) /2
no.non.ad.imaging$hipp_average_harmonized_icv_adj <- (no.non.ad.imaging$ST29SV_harmonized_icv_adj + no.non.ad.imaging$ST88SV_harmonized_icv_adj) /2

no.non.ad.imaging <-  no.non.ad.imaging[complete.cases(no.non.ad.imaging[,c("DX_bl", "PTEDUCAT_bl", "PTGENDER_bl", "CDGLOBAL_bl", 
                                                                  "AGE_bl", "MMSE_bl", "PTAU_bl")]), ]

ggplot(no.non.ad.imaging, aes(x=hipp_average)) + geom_histogram()
ggplot(new_merged_data[aut.rows,], aes(x=hipp_average)) + geom_histogram()
ggplot(new_merged_data[imputed.rows,], aes(x=hipp_average)) + geom_histogram()

lmer(hipp_average ~   PTEDUCAT_bl +  PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl + PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl + PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)

lmer(hipp_average ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)

lmer(hipp_average ~   PTEDUCAT_bl +  ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl + ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl + ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)

lmer(hipp_average ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)



lmer(hipp_average ~   PTEDUCAT_bl +  PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl + PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl + PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])

lmer(hipp_average ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), datanew_merged_data[aut.rows,]no.non.ad.imaging)

lmer(hipp_average ~   PTEDUCAT_bl +  ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl + ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl + ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])

lmer(hipp_average ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[aut.rows,])




lmer(hipp_average ~   PTEDUCAT_bl +  PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl + PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl + PTAU_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])

lmer(hipp_average ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl +  PTAU_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])

lmer(hipp_average ~   PTEDUCAT_bl +  ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = no.non.ad.imaging)
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl + ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl + ptau_pos_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])

lmer(hipp_average ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])
lmer(hipp_average_average_harmonized ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])
lmer(hipp_average_average_harmonized_icv_adj ~   PTEDUCAT_bl +  ptau_pos_unadjusted_bl + AGE_bl + PTGENDER_bl + MMSE_bl + CDGLOBAL_bl + (1|RID), data = new_merged_data[imputed.rows,])
