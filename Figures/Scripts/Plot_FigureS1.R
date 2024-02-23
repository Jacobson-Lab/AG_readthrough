# ----------------------------------------
# Figure S1
# ----------------------------------------

library(dplyr)
library(ggplot2)

metrics <- read.csv("../Data/rf_metrics.csv")
metrics <- reshape2::melt(metrics)
metrics$Sample <- recode_factor(metrics$Sample, Untr = "Untreated", Amik = "Amikacin",
                                G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",
                                Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")

# A
p_nrmse <- ggplot(metrics[which(metrics$metrics == "NRMSE"), ], aes(x = Sample, y = value)) +
  stat_summary(fun = mean, geom = "errorbar", color = "purple", linewidth = 1.2, aes(y = value, ymax = after_stat(y), ymin = after_stat(y))) +
  stat_summary(fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2, position = position_dodge(0.9)) +
  geom_jitter(alpha = 0.5, color = "black", height = 0, width = 0.25) +
  ylab("NRMSE") +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5, size = 6.5), axis.text.y = element_text(size = 6.5), axis.title.x = element_blank())

# B
p_auroc <- ggplot(metrics[which(metrics$metrics == "AUROC"), ], aes(x = Sample, y = value)) +
  stat_summary(fun = mean, geom = "errorbar", color = "purple", linewidth = 1.2, aes(y = value, ymax = after_stat(y), ymin = after_stat(y))) +
  stat_summary(fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "errorbar", width = 0.2, position = position_dodge(0.9)) +
  geom_jitter(alpha = 0.5, color = "black", height = 0, width = 0.25) +
  geom_hline(yintercept = 0.5, color = "red") +
  ylab("AUROC") + xlab("Sample") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5, size = 6.5), axis.text.y = element_text(size = 6.5), axis.title.x = element_blank())

# Combine panels
library(patchwork)
p <- (p_nrmse | p_auroc) + 
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
cairo_pdf(filename = "../Plots/FigureS1.pdf", family = "Arial", width = 7, height = 3) 
p
dev.off()
