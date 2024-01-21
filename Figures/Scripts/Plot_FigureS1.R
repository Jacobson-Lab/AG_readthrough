# ----------------------------------------
# Figure S1
# ----------------------------------------

library(dplyr)
library(ggplot2)

# A
nrmse_meansd <- read.table("../Data/rf_reg_nrmse_meansd.txt", header = TRUE)
nrmse_meansd$Sample <- recode_factor(nrmse_meansd$Sample, Untr = "Untreated", Amik = "Amikacin", 
                                     G418_2000 = "G418\n(2)", G418_500 = "G418\n(0.5)", G418_500_10min = "G418\n(0.5, 10 min)", 
                                     Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")

p_nrmse <- ggplot(nrmse_meansd, aes(x = Sample, y = nrmse_mean)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#999999") +
  geom_errorbar(aes(x = Sample, ymin = nrmse_mean-nrmse_sd, ymax = nrmse_mean+nrmse_sd), width = 0.2, position = position_dodge(0.9)) +
  ylab("NRMSE") +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank())

# B
AUROC_meansd <- read.table("../Data/rf_class_AUROC_meansd.txt", header = TRUE)
AUROC_meansd$Sample <- recode_factor(AUROC_meansd$Sample, Untr = "Untreated", Amik = "Amikacin", 
                                     G418_2000 = "G418\n(2)", G418_500 = "G418\n(0.5)", G418_500_10min = "G418\n(0.5, 10 min)", 
                                     Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")
p_auroc <- ggplot(data = AUROC_meansd, aes(x = Sample, y = AUROC_mean)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#999999") +
  geom_errorbar(aes(ymin = AUROC_mean-AUROC_sd, ymax = AUROC_mean+AUROC_sd), width = 0.2, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0.5, color = "red") +
  ylab("AUROC") + xlab("Sample") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank())

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
cairo_pdf(filename = "../Plots/FigureS1.pdf", family = "Arial", width = 7, height = 3.5) 
p
dev.off()
