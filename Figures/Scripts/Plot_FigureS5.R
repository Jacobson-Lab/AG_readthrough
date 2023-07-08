# ----------------------------------------
# Figure S5
# ----------------------------------------

library(dplyr)
library(ggplot2)
library(ggh4x)

dfdl <- read.table("../Data/wilcoxon_stop_nts_forDualLuc.txt", header = TRUE)
cftr_feat <- read.csv("../Data/CFTR_PTCs_mRNA_features.csv", nrows = 15)
cftr_feat <- cftr_feat[, c("allele_name", "Mutation", unique(dfdl$feature))]
cftr_feat$axis_lab <- paste0(cftr_feat$Mutation, " (", cftr_feat$stop_codon, ")")

cftr <- reshape2::melt(cftr_feat, id.vars = c("allele_name", "Mutation", "axis_lab"))
cftr$value <- sub("T", "U", cftr$value)
cftr$axis_lab <- sub("\\(T", "\\(U", cftr$axis_lab)
cftr <- left_join(cftr, dfdl, by = c("variable" = "feature", "value" = "Var"))
cftr$variable <- recode_factor(cftr$variable, nt_m09 = "-9", nt_m08 = "-8", nt_m07 = "-7", nt_m06 = "-6", nt_m05 = "-5", nt_m04 = "-4",
                               nt_m03 = "-3", nt_m02 = "-2", nt_m01 = "-1", stop_codon = "Stop", nt_p04 = "+4", nt_p05 = "+5", nt_p06 = "+6", 
                               nt_p07 = "+7", nt_p08 = "+8", nt_p09 = "+9", nt_p10 = "+10", nt_p11 = "+11", nt_p12 = "+12")
cftr$facet_lab <- ""
cftr[cftr$variable %in% c("-9", "-8", "-7"), ]$facet_lab <- " "
cftr[cftr$variable %in% c("-6", "-5", "-4"), ]$facet_lab <- "E"
cftr[cftr$variable %in% c("-3", "-2", "-1"), ]$facet_lab <- "P"
cftr[cftr$variable %in% c("Stop"), ]$facet_lab <- "A"
cftr$facet_lab <- factor(cftr$facet_lab, levels = c(" ", "E", "P", "A", ""))
cftr$p_sig <- factor(cftr$p_sig, levels = c("p < 0.05", "ns"))

df <- read.table("../Data/measured_vs_predicted_readthrough.txt", header = TRUE, sep = "\t") # For rank order of measured readthrough
df$G418_treatment <- recode_factor(df$G418_treatment, No = "Untreated", Yes = "G418-treated")
df <- df[df$G418_treatment == "G418-treated", ]
df$axis_lab <- paste0(df$Mutation, " (", df$stop_codon, ")")

rank_measured <- df[order(df$average_measured_rt), ]$axis_lab
cftr$axis_lab <- factor(cftr$axis_lab, levels = rank_measured)

#rank_predicted <- df[order(df$predicted_rt_reporter), ]$axis_lab
#cftr$axis_lab <- factor(cftr$axis_lab, levels = rank_predicted)

cftr <- left_join(cftr, df[, c("allele_name", "stop_codon")], by = "allele_name")
cftr$stop_codon <- factor(cftr$stop_codon, levels = c("UGA", "UAG", "UAA"))

p_cftr2 <- ggplot(cftr, aes(x = variable, y = axis_lab)) +
  geom_tile(aes(fill = median_diff, height = p_hw, width = p_hw, size = p_sig)) +
  geom_text(aes(label = value), size = 8/.pt) +
  scale_fill_gradient2(name = "Group's median readthrough   \nrelative to sample median   \nG418 (0.5)   ",
                       low = "blue" , mid = "white", high = "red", midpoint = 0, na.value = "grey80",
                       limits = c(-1, 1), oob = scales::squish,
                       breaks = seq(-1, 1, length = 5), labels = c("< -1.0", "-0.5\nLower", "0.0", "0.5\nHigher", "> 1.0")) +
  scale_size_manual(values = c(`p < 0.05` = 0, ns = 3), name = "Significance  ") +
  #facet_grid(.~facet_lab, scales = "free", space = "free") + # p_cftr
  facet_grid(stop_codon~facet_lab, scales = "free", space = "free") + # p_cftr2
  xlab("") + ylab(expression(atop("Ordered by readthrough measured in G418-treated condition", "lowest" %->% "highest"))) +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25), strip.placement = "outside") +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.25, "cm")), 
         size = guide_legend(order = 2, override.aes = list(color = "white"), keywidth = unit(0.3, "cm"), keyheight = unit(0.3, "cm")))

# Combine panels
library(patchwork)
p <- (p_cftr / p_cftr2) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face = "bold"), plot.tag.position = "topleft", legend.position = "top")

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "FigureS5.pdf", family = "Arial", width = 7.5, height = 10) 
p
dev.off()
