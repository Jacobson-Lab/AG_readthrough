# ----------------------------------------
# Figure 2
# ----------------------------------------

library(dplyr)
library(ggplot2)

# A
dff <- read.table("../Data/wilcoxon_stop_nts.txt", header = TRUE)
dff$feature <- recode_factor(dff$feature, random_factor = "Random", nt_m03 = "-3", nt_m02 = "-2", nt_m01 = "-1", stop_codon = "Stop", 
                             nt_p04 = "+4", nt_p05 = "+5", nt_p06 = "+6", nt_p07 = "+7", nt_p08 = "+8", nt_p09 = "+9")
dff$p_sig <- factor(dff$p_sig, levels = c("p < 0.05", "ns"))
dff$Sample <- recode_factor(dff$Sample, Untr = "Untreated", Amik = "Amikacin", 
                            G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                            Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")

pA <- ggplot(dff) +
  geom_tile(aes(x = Var, y = Sample, fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~feature, scales = "free", space = "free_x") +
  scale_fill_gradient2(name = "Group's median\nrelative to sample median     ",
                       low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-1, 1), oob = scales::squish,
                       breaks = seq(-1, 1, length = 5), labels = c("< -1.0", "-0.5\nLower", "0.0", "0.5\nHigher", "> 1.0")) +
  scale_size_manual(values = c(`p < 0.05` = 0, ns = 3), name = "Significance  ") +
  scale_y_discrete(limits = rev(levels(dff$Sample))) +
  xlab("") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# B
dff2 <- read.table("../Data/wilcoxon_stop4.txt", header = TRUE)
dff2$p_sig <- factor(dff2$p_sig, levels = c("p < 0.05", "ns"))
dff2$Sample <- recode_factor(dff2$Sample, Untr = "Untreated", Amik = "Amikacin", 
                             G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                             Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")
pB <- ggplot(dff2) +
  geom_tile(aes(x = nt_p04, y = Sample, fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~stop_codon, scales = "free", space = "free_x") +
  scale_fill_gradient2(name = "Group's median\nrelative to sample median     ",
                       low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-1, 1), oob = scales::squish,
                       breaks = seq(-1, 1, length = 5), labels = c("< -1.0", "-0.5\nLower", "0.0", "0.5\nHigher", "> 1.0")
                       ) +
  scale_size_manual(values = c(`p < 0.05` = 0, ns = 3), name = "Significance  ") +
  scale_y_discrete(limits = rev(levels(dff2$Sample))) +
  xlab("") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "none", legend.box = "vertical", legend.direction = "horizontal",
        legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# C
chisq_df <- read.table("../Data/chisq_proportion_stop4.txt", header = TRUE)
chisq_df$p_sig <- factor(chisq_df$p_sig, levels = c("p < 0.05", "ns"))
chisq_df$Sample <- recode_factor(chisq_df$Sample, Untr = "Untreated", Amik = "Amikacin", 
                                 G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                                 Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin", 
                                 Reference = "Reference")
pC <- ggplot(chisq_df) +
  geom_tile(aes(x = nt_p04, y = Sample, fill = standardized_residuals, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~stop_codon, scales = "free", space = "free_x") +
  scale_fill_gradient2(name = expression(paste("Standardized ", chi^2, " residuals ")),
                       low = "forestgreen" , mid = "white", high = "deeppink", midpoint = 0, na.value = "grey80",
                       limits = c(-6, 6), oob = scales::squish,
                       breaks = seq(-6, 6, length = 5), labels = c("< -6", "-3\nRepulsion", "0", "3\nAttraction", "> 6")) +
  scale_size_manual(values = c(`p < 0.05` = 0, ns = 3), name = "Significance  ") +
  scale_y_discrete(limits = rev(levels(chisq_df$Sample))) +
  xlab("nt +4") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(5, "cm"), barheight = unit(0.25, "cm"), title.position = "top"), 
         size = "none")

# Combine plots
library(patchwork)
p <- pA / (pB + pC) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft")

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Figure2.pdf", family = "Arial", width = 7.5, height = 5.5) 
p
dev.off()
