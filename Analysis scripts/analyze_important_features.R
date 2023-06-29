# Analysis of important mRNA features
load("allAG_merged.Rdata") # From prepare_data.R
load("feature_file_wangen_elife_2020.Rdata") # From prepare_data.R
source("functions_analyses.R")
library(data.table)
library(dplyr)

feature_file$stop_ntp04 <- paste0(feature_file$stop_codon, feature_file$nt_p04)
allAG_reg <- lapply(allAG_merged, prep, cds_rpkm_cutoff = 5, utr3_rpkm_cutoff = 0.5, group_param = "none", feature_file = feature_file, utr3_region = "ext")

## 1. Stop codon and surrounding nucleotides
vars_list <- c("random_factor", "nt_m03", "nt_m02", "nt_m01", "stop_codon", "nt_p04", "nt_p05", "nt_p06", "nt_p07", "nt_p08", "nt_p09")
dff <- bind_rows(lapply(allAG_reg, compare_rt_categorical, vars = vars_list, y_var = "log_rte"), .id = "Sample")
dff$Var <- gsub(pattern = "T", replacement = "U", x = dff$Var)
write.table(x = dff, file = "../Figures/Data/wilcoxon_stop_nts.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

## 2. Stop codon and nt +4 together
dff2 <- bind_rows(lapply(allAG_reg, compare_rt_categorical, vars = "stop_ntp04", y_var = "log_rte"), .id = "Sample")
dff2$Var <- gsub("T", "U", dff2$Var)
dff2$stop_codon <- gsub('.{1}$', '', dff2$Var)
dff2$ntp_04 <- gsub('^.{3}', '', dff2$Var)
write.table(x = dff2, file = "../Figures/Data/wilcoxon_stop4.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

  ## Frequency of and Chi-square test of association between stop codon and nt +4
allAG_reg_ref <- allAG_reg
allAG_reg_ref$Reference <- feature_file

stop4_chisq <- lapply(allAG_reg_ref, function(x) { 
  cc <- table(x$stop_codon, x$nt_p04)
  cc <- chisq.test(cc)
  # Chi-square test results
  pval <- cc$p.value
  resi <- reshape2::melt(cc$residuals)
  stdresi <- reshape2::melt(cc$stdres)
  chi_resi <- left_join(resi, stdresi, by = c("Var1", "Var2"), suffix = c("_raw_residuals", "_standardized_residuals"))
  chi_resi <- cbind(chi_resi, pval)
  # Observed
  ob_count <- reshape2::melt(cc$observed)
  ob_prop <- reshape2::melt(prop.table(cc$observed))
  ob <- left_join(ob_count, ob_prop, by = c("Var1", "Var2"), suffix = c("_observed_count", "_observed_fraction"))
  # Expected
  ex_count <- reshape2::melt(cc$expected)
  ex_prop <- reshape2::melt(prop.table(cc$expected))
  ex <- left_join(ex_count, ex_prop, by = c("Var1", "Var2"), suffix = c("_expected_count", "_expected_fraction"))
  # Create output table
  output <- left_join(ex, ob, by = c("Var1", "Var2"))
  output <- left_join(output, chi_resi, by = c("Var1", "Var2"))
  colnames(output) <- sub("value_", "", colnames(output))
  colnames(output)[1:2] <- c("stop_codon", "nt_p04")
  return(output)
})
stop4_chisq <- bind_rows(stop4_chisq, .id = "Sample")
stop4_chisq <- within(stop4_chisq, {
  p_sig = ifelse(pval < 0.05, "p < 0.05", "ns")
  p_hw = ifelse(pval < 0.05, 0.9, 0.6)
})
stop4_chisq$stop_codon <- gsub("T", "U", stop4_chisq$stop_codon)
stop4_chisq$nt_p04 <- gsub("T", "U", stop4_chisq$nt_p04)
write.table(x = stop4_chisq, file = "../Figures/Data/chisq_proportion_stop4.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

## 3. P-site codon
dfpc <- bind_rows(lapply(allAG_reg, compare_rt_categorical, vars = "codon_m01", y_var = "log_rte"), .id = "Sample")
dfpc$Var <- gsub("T", "U", dfpc$Var)
write.table(x = dfpc, file = "../Figures/Data/wilcoxon_psite_codon.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

## 4. 3'-UTR length
dfu <- bind_rows(lapply(allAG_reg, function(x) x[, c("log_rte", "l_utr3")]), .id = "Sample")
write.table(x = dfu, file = "../Figures/Data/readthrough_3utr_length.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# --------------------------------------------------
## Stop codon and nucleotides used in Dual-Luc assay
vars_list <- c(paste0("nt_m0", seq(1, 9, 1)), "stop_codon", paste0("nt_p0", seq(4, 9, 1)), paste0("nt_p", seq(10, 12, 1)))
dfdl <- compare_rt_categorical(dat = allAG_reg$G418_500, vars = vars_list, y_var = "log_rte")
dfdl$Var <- gsub(pattern = "T", replacement = "U", x = dfdl$Var)
write.table(x = dfdl, file = "../Figures/Data/wilcoxon_stop_nts_forDualLuc.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

