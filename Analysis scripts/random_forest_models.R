# Create random forest models
load("allAG_merged.Rdata") # From prepare_data.R
load("feature_file_wangen_elife_2020.Rdata") # From prepare_data.R
source("functions_analyses.R")
library(data.table)
library(dplyr)
library(randomForest)
library(caret)

# 1. Prepare data ------------------------------------
# 1.1 Regression
allAG_reg <- lapply(allAG_merged, prep, cds_rpkm_cutoff = 5, utr3_rpkm_cutoff = 0.5, group_param = "none", feature_file = feature_file, utr3_region = "ext")
# 1.2 Classification
allAG_group <- lapply(allAG_merged, prep, cds_rpkm_cutoff = 5, utr3_rpkm_cutoff = 0.5, group_param = c("percentile", 15), feature_file = feature_file, utr3_region = "ext")

# 2. Random forest model -----------------------------
  # Choose columns for x variables
all_col <- names(allAG_reg[[1]])
feature_to_run <- c(all_col[grepl("tunnel", all_col)], "aa_m02", "aa_m01", all_col[grepl("nt_m", all_col)], "stop_codon", all_col[grepl("nt_p", all_col)], "l_utr3", "random_factor", "random_num")
# 2.1 Regression
set.seed(5490)
allAG_rf_reg <- lapply(allAG_reg, rfFit, y_val = "log_rte", col_feature = feature_to_run, 
                       ntree_ = 100, cv_fold = 5, met = "RMSE", dummy = FALSE, mtry_ = seq(1, length(feature_to_run), 10))
save(allAG_rf_reg, file = "allAG_rf_reg.Rdata")
# 2.2 Classification
set.seed(2423)
allAG_rf_class <- lapply(allAG_reg, rfFit, y_val = "Group", col_feature = feature_to_run, 
                         ntree_ = 100, cv_fold = 5, met = "ROC", dummy = FALSE, mtry_ = seq(1, length(feature_to_run), 10))
save(allAG_rf_class, file = "allAG_rf_class.Rdata")

# 3. Performance metric -----------------------------
# 3.1 Regression
nrmse <- lapply(allAG_rf_reg, function(x) {xx <- x$resample[, "RMSE"] / (max(x$trainingData$.outcome)-min(x$trainingData$.outcome)); return(xx)})
nrmse_meansd <- lapply(nrmse, function(x) {df <- data.frame(nrmse_mean = mean(x), nrmse_sd = sd(x))})
nrmse_meansd <- bind_rows(nrmse_meansd, .id = "Sample")
write.table(x = nrmse_meansd, file = "../Figures/Data/rf_reg_nrmse_meansd.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 3.2 Classification
AUROC <- lapply(allAG_rf_class, function(x) {xx <- x$resample[, "ROC"]; return(xx)})
AUROC_meansd <- lapply(AUROC, function(x) {df <- data.frame(AUROC_mean = mean(x), AUROC_sd = sd(x))})
AUROC_meansd <- bind_rows(AUROC_meansd, .id = "Sample")
write.table(x = AUROC_meansd, file = "../Figures/Data/rf_class_AUROC_meansd.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# 4. Feature importance ------------------------------
# 4.1 Regression
all.imp_reg <- lapply(allAG_rf_reg, function(x) {imp <- as.data.frame(importance(x$finalModel)); imp$feature <- row.names(imp); return(imp)})
all.imp_reg <- bind_rows(all.imp_reg, .id = "Sample")
write.table(x = all.imp_reg, file = "../Figures/Data/rf_reg_feature_importance.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# 4.2 Classification
all.imp_class <- lapply(allAG_rf_class, function(x) {imp <- as.data.frame(importance(x$finalModel)); imp$feature <- row.names(imp); return(imp)})
all.imp_class <- bind_rows(all.imp_class, .id = "Sample")
write.table(x = all.imp_class, file = "../Figures/Data/rf_class_feature_importance.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
