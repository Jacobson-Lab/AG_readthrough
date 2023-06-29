# ----------------------------------------
# Figure 7
# ----------------------------------------

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(rstatix)

df <- read.table("../Data/measured_vs_predicted_readthrough.txt", header = TRUE, sep = "\t")
df$G418_treatment <- recode_factor(df$G418_treatment, No = "Untreated", Yes = "G418-treated")

fc <- read.table("../Data/measured_vs_predicted_readthrough_G418-treated_Untreated_foldchange.txt", header = TRUE, sep = "\t")

fc2 <- left_join(fc[, 1:4], df[df$G418_treatment == "Untreated", c("allele_name", "average_measured_rt")], by = "allele_name")
fc2$log2fc <- log2(fc2$foldchange_measured_rt)

# Fold change vs basal readthrough (untreated)
pA <- ggplot(fc2, aes(x = average_measured_rt, y = log2fc)) +
  geom_point(aes(color = stop_codon), alpha = 0.5) +
  geom_text_repel(aes(label = Mutation, color = stop_codon), 
                  size = 8/.pt, show.legend = FALSE, direction = "both",
                  box.padding = 0.1, point.padding = 0.1, max.overlaps = Inf, min.segment.length = 0) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(name = "Stop codon", values = c("red", "orange", "forestgreen")) + 
  xlab("Average measured readthrough\nUntreated") + 
  ylab(expression(atop("Average measured readthrough", "log"[2]*"(G418 / Untreated)"))) +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), legend.position = "right")

# Fold change vs stop codon
comp <- get_comparisons(fc2, variable = "stop_codon")
stat.test <- fc2 %>% 
  wilcox_test(log2fc~stop_codon, comparisons = comp, p.adjust.method = "BH") %>%
  add_significance() %>% add_xy_position()
pB <- ggplot(fc2, aes(x = stop_codon, y = log2fc, color = stop_codon)) +
  geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  stat_pvalue_manual(data = stat.test, size = 8/.pt, tip.length = 0.02) +
  scale_color_manual(name = "Stop codon", values = c("red", "orange", "forestgreen")) + 
  xlab("") + ylab(expression(atop("Average measured readthrough", "log"[2]*"(G418 / Untreated)"))) +
  theme_bw(base_size = 10) + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), legend.position = "none")

# Combine panels
library(patchwork)
p <- (pA | pB) +
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
cairo_pdf(filename = "Figure7.pdf", family = "Arial", width = 7.5, height = 3.5) 
p
dev.off()
