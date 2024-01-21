# ----------------------------------------
# Figure S7
# ----------------------------------------

library(dplyr)
library(ggplot2)
library(ggh4x)
library(rstatix)
library(ggpubr)

# Model vs. assay features
  # Run Plot_Figure1.R for imp_reg, then all.imp_df <- imp_reg, and all.imp_df preparation
feats <- unique(all.imp_df[, c("feature", "xtick", "xgroup", "xgroup2")])
assay_feat <- c(paste0("-", seq(1, 9, 1)), "E-site aa", "P-site aa", "Stop codon", paste0("+", seq(4, 12, 1)))

feat_model <- read.table("../Data/rf_all_model_features_used.txt", header = TRUE)
feat_model <- left_join(feat_model, feats, by = "feature")
feat_model$X.IncMSE <- NA
feat_model[feat_model$Sample == "Full", ]$X.IncMSE <- all.imp_df[all.imp_df$Sample == "G418 (0.5)", ]$X.IncMSE
feat_model$assay <- "No"
feat_model[feat_model$Sample == "Full" & feat_model$xtick %in% assay_feat, ]$assay <- "Yes"
feat_model[feat_model$Sample == "Full" & feat_model$xgroup == "aa 3-9\nfrom PTC", ]$assay <- "Yes_partial"
feat_model$assay <- factor(feat_model$assay, levels = c("Yes", "Yes_partial", "No"))
feat_model$Sample <- factor(feat_model$Sample, levels = c("Full", "Assay", "Reduced", "Mock1", "Mock2", "Mock3"))

n_feat <- data.frame(table(Sample = feat_model$Sample, used = feat_model$used != "FALSE"))
n_feat <- n_feat[n_feat$used == TRUE, c("Sample", "Freq")]
n_feat$Freq <- n_feat$Freq - 2 # Excluding the two negative controls from count
n_feat$Sample2 <- paste0(n_feat$Sample, " (", n_feat$Freq, ")")

feat_model <- left_join(feat_model, n_feat[, c("Sample", "Sample2")], by = "Sample")
feat_model$Sample2 <- factor(feat_model$Sample2, levels = c("Full (75)", "Assay (21)", "Reduced (22)", "Mock1 (22)", "Mock2 (22)", "Mock3 (22)"))

p_feat_model <- ggplot() +
  # %IncMSE
  geom_tile(data = feat_model, 
            aes(x = xtick, y = Sample2, fill = X.IncMSE)) +
  scale_fill_gradient2(name = "G418 (0.5)\n%IncMSE  ", low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-15, 15), oob = scales::squish, na.value = "white") +
  # Add green border around PTC context, light green around features partially from PTC context and reporter
  geom_tile(data = feat_model[feat_model$Sample == "Full" & feat_model$assay == "Yes", ], inherit.aes = FALSE, na.rm = TRUE,
            aes(x = xtick, y = Sample2), color = "forestgreen", fill = NA, linewidth = 0.5) +
  geom_tile(data = feat_model[feat_model$Sample == "Full" & feat_model$assay == "Yes_partial", ], inherit.aes = FALSE, na.rm = TRUE,
            aes(x = xtick, y = Sample2), color = "lightgreen", fill = NA, linewidth = 0.5) +
  # Features used in each model
  geom_point(data = feat_model, inherit.aes = FALSE,
             aes(x = xtick, y = Sample2, alpha = used, color = used), shape = 4, na.rm = TRUE, show.legend = FALSE)  +
  scale_alpha_manual(values = c(`TRUE` = 1, TRUE_not_exact = 1, `FALSE` = 0), na.value = NA) +
  scale_color_manual(values = c(`TRUE` = "black", TRUE_not_exact = "red", `FALSE` = "white"), na.value = NA) +
  # Facet
  facet_nested(.~xgroup2+xgroup, scales = "free", space = "free", nest_line = TRUE) +
  scale_y_discrete(limits = rev(levels(feat_model$Sample2))) +
  labs(tag = "B") +
  theme_bw(base_size = 7) +
  theme(legend.position = "top", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.25, "cm"), 
        legend.box.spacing = unit(0, "cm"), legend.title.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(face = "italic"), 
        axis.title = element_blank(),
        strip.background = element_blank(), plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft",
        panel.border = element_rect(colour = "black", fill = NA), panel.grid = element_blank(), panel.spacing = unit(0,"cm"))

# Performance metrics of different models
metrics_model <- read.table("../Data/rf_all_model_metrics_G418_500.txt", header = TRUE)
metrics_model$Model <- factor(metrics_model$Model, levels = c("Full", "Assay", "Reduced", "Mock1", "Mock2", "Mock3"))
stat.test <- metrics_model %>% group_by(metrics) %>%
  t_test(value~Model, p.adjust.method = "BH", paired = FALSE, ref.group = "Full") %>%
  add_significance() %>% add_xy_position(x = "Model", scales = "free_y")

p_nrmse <- ggplot(metrics_model[metrics_model$metrics == "NRMSE", ], 
                  aes(x = Model, y = value)) +
  stat_summary(geom = "errorbar", fun = mean, linewidth = 1.2, color = "purple", 
               aes(y = value, ymax = after_stat(y), ymin = after_stat(y))) +
  geom_point(alpha = 0.5) +
  stat_pvalue_manual(data = stat.test[stat.test$metrics == "NRMSE", ], 
                     label = "p.adj", tip.length = 0.01, size = 6/.pt, hide.ns = FALSE) +
  ylab("NRMSE") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), strip.background.x = element_rect(fill = "white"), strip.placement = "outside", 
        strip.background.y = element_rect(fill = "white", color = "white"))

p_auroc <- ggplot(metrics_model[metrics_model$metrics == "AUROC", ], 
                  aes(x = Model, y = value)) +
  stat_summary(geom = "errorbar", fun = mean, linewidth = 1.2, color = "purple", 
               aes(y = value, ymax = after_stat(y), ymin = after_stat(y))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0.5, color = "red") +
  stat_pvalue_manual(data = stat.test[stat.test$metrics == "AUROC", ], 
                     label = "p.adj", tip.length = 0.01, size = 6/.pt, hide.ns = FALSE) +
  ylab("AUROC") + 
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), strip.background.x = element_rect(fill = "white"), strip.placement = "outside", 
        strip.background.y = element_rect(fill = "white", color = "white"))

# Combine panels
library(patchwork)
p_metrics <- (p_nrmse + theme(axis.title.y = element_text(margin = margin(r = -30, unit = "pt")))) | p_auroc
p <- (p_feat_model / p_metrics) + plot_layout(heights = c(1.5, 4.2)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = "topleft")

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "../Plots/FigureS7.pdf", family = "Arial", width = 7, height = 4.7) 
p
dev.off()
