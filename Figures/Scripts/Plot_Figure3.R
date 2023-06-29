# ----------------------------------------
# Figure 3
# ----------------------------------------

library(dplyr)
library(ggplot2)
library(Biostrings)

dfpc <- read.table("../Data/wilcoxon_psite_codon.txt", header = TRUE)
dfpc <- cbind(dfpc, aa = sapply(dfpc$Var, function(x) {y <- as.character(translate(RNAString(x), no.init.codon = TRUE)); return(y)}))
dfpc$aa <- factor(dfpc$aa, levels = c("F", "S", "Y", "*", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"))
dfpc$p_sig <- factor(dfpc$p_sig, levels = c("p < 0.05", "ns"))
dfpc$Sample <- recode_factor(dfpc$Sample, Untr = "Untreated", Amik = "Amikacin",
                             G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)", 
                             Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")

p <- ggplot(dfpc) +
  geom_tile(aes(x = Var, y = Sample, fill = median_diff, height = p_hw, width = p_hw, size = p_sig), color = NA) +
  facet_grid(.~aa, scales = "free", space = "free_x") +
  scale_fill_gradient2(name = "Group's median\nrelative to sample median     ",
                       low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-3, 3), oob = scales::squish,
                       breaks = seq(-3, 3, length = 5), labels = c("< -3.0", "-1.5\nLower", "0.0", "1.5\nHigher", "> 3.0")) +
  scale_size_manual(values = c(`p < 0.05` = 0, ns = 3), name = "Significance  ") +
  scale_y_discrete(limits = rev(levels(dfpc$Sample))) +
  xlab("") + ylab("") +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "Figure3.pdf", family = "Arial", width = 7.5, height = 2.5) 
p
dev.off()
