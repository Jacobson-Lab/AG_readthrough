# Predict readthrough of new data using the trained random forest regression models
library(data.table)
library(dplyr)
library(randomForest)
source("functions_analyses.R")

# Predict readthrough of CFTR PTCs ------------------------------------
  # Adapt code in prepare_data.R to get mRNA features for each to be predicted for readthrough
cftr_tab <- read.csv("../Figures/Data/CFTR_PTCs_mRNA_features.csv")

load("allAG_rf_reg.Rdata")
names(allAG_rf_reg)
rf_list <- allAG_rf_reg[1:2]

# Predict with native sequence
pred_native <- lapply(rf_list, predict_rt, newdata = cftr_tab, outside_features = NULL)

# Predict with reporter sequences, keep the stop codon and flanking 9 nt on each side of the stop codon native
features_3X3 <- c("stop_codon", "nt_p04", "nt_p05", "nt_p06", "nt_p07", "nt_p08", "nt_p09", "nt_p10", "nt_p11", "nt_p12",
                  "nt_m01", "nt_m02", "nt_m03", "nt_m04", "nt_m05", "nt_m06", "nt_m07", "nt_m08", "nt_m09", "codon_m01", "aa_m01", "codon_m02", "aa_m02",
                  "random_num", "random_factor")
pred_dualluc <- lapply(rf_list, predict_rt, newdata = cftr_tab, outside_features = 16, keep_native = features_3X3)

# Export data ------------------------------------
pred_native <- bind_rows(pred_native, .id = "G418_treatment")
pred_native$G418_treatment <- recode_factor(pred_native$G418_treatment, Untr = "Untreated", G418_500 = "G418-treated")
pred_dualluc <- bind_rows(pred_dualluc, .id = "G418_treatment")
pred_dualluc$G418_treatment <- recode_factor(pred_dualluc$G418_treatment, Untr = "Untreated", G418_500 = "G418-treated")

measured <- read.csv("../Figures/Data/CFTR_PTCs_DualLuc_measured_readthrough.csv")
rt <- left_join(measured[, c("allele_name", "Mutation", "stop_codon", "G418_treatment", "Average")], pred_native, 
                by = c("allele_name" = "allele", "Mutation", "G418_treatment"))
rt <- left_join(rt, pred_dualluc, by = c("allele_name" = "allele", "Mutation", "G418_treatment"), suffix = c("_native", "_reporter"))
colnames(rt)[5] <- "average_measured_rt"
write.table(rt, file = "../Figures/Data/measured_vs_predicted_readthrough.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Calculate fold change G418-treated / Untreated
measured_foldchange <- reshape2::melt(measured, id.vars = c("allele_name", "Mutation", "stop_codon", "G418_treatment"), 
                                      measure.vars = c("rep1", "rep2", "rep3", "rep4", "rep5", "rep6", "rep7"), 
                                      variable.name = "rep", value.name = "rt")
measured_foldchange <- reshape2::dcast(measured_foldchange, allele_name + Mutation + stop_codon + rep ~ G418_treatment, value.var = "rt")
measured_foldchange$foldchange <- measured_foldchange$`G418-treated` / measured_foldchange$Untreated
measured_foldchange <- reshape2::dcast(measured_foldchange, allele_name + Mutation + stop_codon ~ rep, value.var = "foldchange")
measured_foldchange$Average <- rowMeans(measured_foldchange[, 4:10], na.rm = TRUE)
measured_foldchange <- cbind(measured_foldchange, SD = apply(measured_foldchange[, 4:10], 1, sd, na.rm = TRUE))

predicted_foldchange <- reshape2::melt(rt, id.vars = c("allele_name", "Mutation", "stop_codon", "G418_treatment"), 
                                       measure.vars = c("predicted_rt_native", "predicted_rt_reporter"), 
                                       variable.name = "scheme", value.name = "rt")
predicted_foldchange <- reshape2::dcast(predicted_foldchange, allele_name + Mutation + stop_codon + scheme ~ G418_treatment, value.var = "rt")
predicted_foldchange$foldchange <- 2^predicted_foldchange$`G418-treated` / 2^predicted_foldchange$Untreated
predicted_foldchange <- reshape2::dcast(predicted_foldchange, allele_name + Mutation + stop_codon ~ scheme, value.var = "foldchange")

fc <- left_join(measured_foldchange[, c("allele_name", "Mutation", "stop_codon", "Average")], predicted_foldchange, by = c("allele_name", "Mutation", "stop_codon"))
colnames(fc)[4] <- "foldchange_measured_rt"
colnames(fc)[5:6] <- sub("^", "foldchange_", colnames(fc)[5:6])
write.table(fc, file = "../Figures/Data/measured_vs_predicted_readthrough_G418-treated_Untreated_foldchange.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

