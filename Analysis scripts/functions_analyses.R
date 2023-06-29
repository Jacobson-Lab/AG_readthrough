# -------------------------------------------------------------------------------------------------------------
# Functions for random forest and comparative analyses
# -------------------------------------------------------------------------------------------------------------
# Prepare data for random forest: Filter, log-transformed readthrough efficiency, and combine with feature_file
prep <- function(df, cds_rpkm_cutoff = 5, utr3_rpkm_cutoff = 0.5, group_param = c("std", 1), feature_file, 
                 utr3_region = "ext", logbase = 2) {
  # df: data table containing readthrough efficiency (rte) and RPKM values for CDS, 3'-UTR, extension regions
  # feature_file: data table containing mRNA features for each mRNA
  # logbase: log base to transform readthrough efficiency values
  ## Filtering criteria
  # cds_rpkm_cutoff: mRNAs with RPKM of CDS region less than the indicated value are discarded.
  # utr3_rpkm_cutoff: mRNAs with RPKM of 3' region (choice of 3'-UTR or extension) less than the indicated value are discarded.
  # utr3_region: 3' region of readthrough column and for utr3_rpkm_cutoff
  ## Grouping criteria
  # group_param: choice of "none" (do not split data into groups) or a vector containing mode of grouping and integer
    # e.g. group_param = c("std", 1) indicate mode of grouping to be standard deviation of 1. 
      # "High" and "Low" readthrough groups are mRNAs with rte above and below 1 sd, respectively
    # e.g. group_param = c("percentile", 15) indicate mode of grouping to be 15th percentile
      # "High" and "Low" readthrough groups are mRNAs with rte in the top and bottom 15%, respectively
  require(data.table)
  # 1. Filter out un-analyzable rte
  df <- data.table(df)
  df <- df[df$rpkm_cds > cds_rpkm_cutoff & df[[paste0("rpkm_", utr3_region)]] > utr3_rpkm_cutoff, ]
  # 2. log-transformed readthrough efficiency
  df$log_rte <- log(x = df[[paste0("rte_", utr3_region)]], base = logbase)
  # 3. Split data into high and low readthrough efficiency groups
  if (group_param[1] == "std") { # 3.1 Using standard deviation
    message(paste0("Group mRNAs using standard deviation of ", group_param[2]))
    df$z <- data.table(scale(df$log_rte, center = TRUE, scale = TRUE))
    std <- as.integer(group_param[2])
    df <- data.table(df[z < -std | z > std, ])
    df$Group <- df$z > std
    df$Group <- gsub("FALSE", "Low", df$Group)
    df$Group <- gsub("TRUE", "High", df$Group)
  } else if (group_param[1] == "percentile") {  # 3.2 Using percentile
    message(paste0("Group mRNAs using top and bottom ", group_param[2]), "th percentile")
    df$percentile <- data.table(dplyr::percent_rank(df$log_rte))
    pt <- as.integer(group_param[2])/100
    df <- data.table(df[percentile < pt | percentile > 1-pt])
    df$Group <- df$percentile > pt
    df$Group <- gsub("FALSE", "Low", df$Group)
    df$Group <- gsub("TRUE", "High", df$Group)
  } else if (group_param[1] == "none") {  # 3.3 No grouping
    message("Data will not be grouped")
  } else {
    message("Invalid grouping parameter")
  }
  # 4. Combine data with feature_file
  df <- left_join(df[, c("transcript", "rpkm_cds", "rpkm_utr5", "rpkm_utr3", "rpkm_ext", "rte_ext", "log_rte")], feature_file[, -c("sequence")], by = "transcript")
  # 5. Final clean up
  df <- df[which(df$l_utr3 >= 10), ] # Remove mRNAs with 3'-UTR shorter than 9 nt
  df <- df[complete.cases(df), ] # Remove any rows with NA
  df <- df %>% mutate_if(is.character, as.factor) # Convert any character variables to factor
  df <- data.table(df)
  return(df)
}

# -------------------------------------------------------------------------------------------------------------
# Random forest model
rfFit <- function(df, y_val, col_feature, ntree_ = 500, cv_fold = 5, met = "Rsquared", dummy = FALSE, mtry_ = 9) {
  # df: data table output from prep()
  # y_val: column name (string) of the y variable
  # col_feature: vector of column names of the x variables
  # ntree_: number of trees for random forest (integer)
  # cv_fold: number of fold for cross-validation (integer)
  # met: performance metric. Make sure it matches the type of random forest: regression or classification
  # dummy: whether to make column with multiple levels dummy variables (TRUE or FALSE)
  # mtry_: number of x variables to subset at each node (integer or a vector of integers)
  require(randomForest)
  require(caret)
  require(data.table)
  message(paste(samp, "\nStart"))
  Grid <- expand.grid(mtry = mtry_)
  df <- data.table(df)
  x_dat <- df[, col_feature, with = FALSE]
  y_dat <- df[[y_val]]
  dat <- df[, c(y_val, col_feature), with = FALSE]
  if (is.numeric(y_dat) == TRUE) {
    train_control <- trainControl(method = "cv", number = cv_fold)
  } else {
    train_control <- trainControl(method = "cv", number = cv_fold, 
                                  classProbs = TRUE, savePredictions = TRUE, summaryFunction = twoClassSummary)
  }
  if(dummy == FALSE) {
    message("Dummy = False")
    all.rfFit[[i]] <- train(x = x_dat, y = y_dat, method = "rf", data = dat, # same as randomForest
                            trControl = train_control, importance = TRUE, metric = met, ntree = ntree_, tuneGrid = Grid)
  } else {
    message("Dummy = True")
    f <- as.formula(paste(y_val, "~."))
    all.rfFit[[i]] <- train(f, method = "rf", data = dat, # same as randomForest
                            trControl = train_control, importance = TRUE, metric = met, ntree = ntree_, tuneGrid = Grid)
  }
  message("Done\n")
  names(all.rfFit) <- samp_list
  return(all.rfFit)
}
# -------------------------------------------------------------------------------------------------------------
# Compare readthrough efficiency between mRNA groups, divided into groups by categorical X variable
compare_rt_categorical <- function(dat, vars, y_var = "log_rte") {
  # dat: data table output from prep() with group_param = "none"
  # vars: a vector containing column name of categorical variables to analyze
  # y_var: column name of y variable (continuous variable)
  compare_rt <- function(dat, var, y_var) {
    df <- dat[, c(y_var, var), with = FALSE]
    names(df) <- c("y_var", "Var")
    cm <- ggpubr::compare_means(formula = y_var~Var, data = df, method = "wilcox.test", ref.group = ".all.", p.adjust.method = "BH")
    tb <- boxplot(y_var~Var, data = df, plot = FALSE)
    tbs <- data.frame(tb$stats)
    colnames(tbs) <- tb$names
    tbs <- data.frame(t(tbs))
    tbs$group2 <- row.names(tbs)
    res <- left_join(cm[, c("group2", "p", "p.adj")], tbs[, c("group2", "X3")], by = "group2")
    res$samp_median <- median(df$y_var)
    names(res) <- c("Var", "p", "p.adj", "median", "samp_median")
    res$median_diff <- res$median - res$samp_median
    res <- within(res, {
      p_sig = ifelse(p.adj < 0.05, "p < 0.05", "ns")
      p_hw = ifelse(p.adj < 0.05, 0.9, 0.6)
    })
    res$feature <- var
    res$p_sig <- factor(res$p_sig, levels = c("p < 0.05", "ns"))
    return(res)
  }
  res_list <- list()
  for (x in vars) {
    res_list[[x]] <- compare_rt(dat, x, y_var)
  }
  res_df <- data.table(bind_rows(res_list))
  return(res_df)
}
# -------------------------------------------------------------------------------------------------------------
# Use random forest model to predict new data
predict_rt <- function(rf, newdata, outside_features = NULL, keep_native) {
  # rf: trained random forest model
  # newdata: dataframe containing sequence and features of mRNAs to be predicted; all features present in rf have to also be present here with the same names
  # outside_features: NULL if predicting with native sequence
  # In the case that some features need to be held constant (e.g. to mimic reporter assay):
    # keep_native: vector of column names of features to keep as is (features from native sequence)
    # outside_features: row number in newdata of the entry containing features to be used to replace native features for all other entries
  require(randomForest)
  features_to_use <- row.names(rf$finalModel$importance)
  if (is.null(outside_features)) {
    newdata_to_pred <- newdata[, features_to_use]
  } else if (is.numeric(outside_features) & outside_features <= nrow(newdata)) {
    # Use data of one representative entry in newdata
    keepcols <- union(c("allele_name", "Mutation", "l_tr", "l_cds", "l_utr5", "sequence"), keep_native)
    seq_cds <- substr(newdata$sequence, start = newdata$l_utr5 + 1, stop = newdata$l_tr - newdata$l_utr3)
    seq_codon_m30 <- substr(seq_cds, start = newdata$l_cds - 89, stop = newdata$l_cds)
    seq_codon_m30 <- seqinr::as.SeqFastadna(seq_codon_m30)
    seq_aa_m30 <- data.frame(matrix(nrow = length(seq_codon_m30), ncol = 30))
    for (i in 1:length(seq_codon_m30)) {
      aa <- seqinr::getTrans(seqinr::s2c(seq_codon_m30[i]))
      seq_aa_m30[i, (30-length(aa)+1):30] <- aa
    }
    source("aa_prop.R")
    seq_aa_m30[, 1:27] <- seq_aa_m30[outside_features, 1:27]
    newdata[, !(colnames(newdata) %in% keepcols)] <- newdata[outside_features, !(colnames(newdata) %in% keepcols)]
    newdata[, grepl(pattern = "upper", x = colnames(newdata))] <- NULL
    newdata <- aa_prop(seq_aa_m30 = seq_aa_m30, aa_range = 22:28, tunnel_part = "upper", feature_file = newdata)
    newdata_to_pred <- newdata[, features_to_use]
  } else {
    message("Invalid outside_features. Valid input is NULL or row number in newdata")
  }
  # Set levels to match with training data
  train <- rf$trainingData[, features_to_use]
  xtest <- rbind(train[1, ] , newdata_to_pred)
  xtest <- xtest[-1, ]
  res <- data.frame(allele = newdata$allele_name, Mutation = newdata$Mutation,
                    predicted_rt = predict(rf, newdata = xtest))
  return(res)
}
# -------------------------------------------------------------------------------------------------------------





