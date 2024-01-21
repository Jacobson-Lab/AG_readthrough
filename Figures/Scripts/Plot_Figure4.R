# ----------------------------------------
# Figure 4
# ----------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Dual-Luc measurement vs. Predicted readthrough
df <- read.table("../Data/measured_vs_predicted_readthrough.txt", header = TRUE, sep = "\t")
df$G418_treatment <- factor(df$G418_treatment, levels = c("Untreated", "G418-treated"))
df <- df[which(df$Mutation != "S434X"), ] # Remove S434X because Dual-Luc measurements were compromised by non-readthrough signals

p_rt_reporter <- ggplot(df, aes(x = average_measured_rt, y = predicted_rt_reporter)) + 
  geom_point(aes(color = stop_codon), alpha = 0.5) + 
  geom_text_repel(aes(label = Mutation, color = stop_codon), size = 6/.pt, show.legend = FALSE,
                  box.padding = 0.1, point.padding = 0.1, max.overlaps = Inf, min.segment.length = 0) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5, linewidth = 0.5, color = "grey50") + 
  stat_cor(method = "spearman", cor.coef.name = "rho", size = 6/.pt,
           label.x = Inf, label.y = -Inf, hjust = 1.1, vjust = -0.1, show.legend = FALSE) +
  facet_wrap(~G418_treatment, scales = "free", nrow = 1) +
  scale_color_manual(name = "Stop codon", values = c("red", "orange", "forestgreen")) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) + 
  xlab("Average measured readthrough") + ylab("Predicted readthrough\n(reporter sequence)") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), 
        legend.position = "top")

p_rt_native <- ggplot(df, aes(x = average_measured_rt, y = predicted_rt_native)) + 
  geom_point(aes(color = stop_codon), alpha = 0.5) + 
  geom_text_repel(aes(label = Mutation, color = stop_codon), size = 6/.pt, show.legend = FALSE,
                  box.padding = 0.1, point.padding = 0.1, max.overlaps = Inf, min.segment.length = 0) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5, linewidth = 0.5, color = "grey50") + 
  stat_cor(method = "spearman", cor.coef.name = "rho", size = 6/.pt, 
           label.x = Inf, label.y = -Inf, hjust = 1.1, vjust = -0.1, show.legend = FALSE) +
  facet_wrap(~G418_treatment, scales = "free", nrow = 1) +
  scale_color_manual(name = "Stop codon", values = c("red", "orange", "forestgreen")) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) + 
  xlab("Average measured readthrough") + ylab("Predicted readthrough\n(native sequence)") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), 
        legend.position = "top")

p_rt_compare <- ggplot(df, aes(x = predicted_rt_native, y = predicted_rt_reporter)) + 
  geom_point(aes(color = stop_codon), alpha = 0.5) + 
  geom_text_repel(aes(label = Mutation, color = stop_codon), size = 6/.pt, show.legend = FALSE,
                  box.padding = 0.1, point.padding = 0.1, max.overlaps = Inf, min.segment.length = 0) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5, linewidth = 0.5, color = "grey50") + 
  stat_cor(method = "spearman", cor.coef.name = "rho", size = 6/.pt, 
           label.x = Inf, label.y = -Inf, hjust = 1.1, vjust = -0.1, show.legend = FALSE) +
  facet_wrap(~G418_treatment, scales = "free", nrow = 1) +
  scale_color_manual(name = "Stop codon", values = c("red", "orange", "forestgreen")) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) + 
  xlab("Predicted readthrough (native sequence)") + ylab("Predicted readthrough\n(reporter sequence)") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), 
        legend.position = "right")

# Dual-Luc measurement vs. Predicted foldchange
fc <- read.table("../Data/measured_vs_predicted_readthrough_G418-treated_Untreated_foldchange.txt", header = TRUE, sep = "\t")
fc$facet_label <- "G418-treated / Untreated"
fc <- fc[which(fc$Mutation != "S434X"), ] # Remove S434X because Dual-Luc measurements were compromised by non-readthrough signals

p_fc_reporter <- ggplot(fc, aes(x = foldchange_measured_rt, y = foldchange_predicted_rt_reporter)) +
  geom_point(aes(color = stop_codon), alpha = 0.5) + 
  geom_text_repel(aes(label = Mutation, color = stop_codon), size = 6/.pt, show.legend = FALSE,
                  box.padding = 0.1, point.padding = 0.3, max.overlaps = Inf, min.segment.length = 0) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5, linewidth = 0.5, color = "grey50") +
  stat_cor(method = "spearman", cor.coef.name = "rho", size = 6/.pt, label.x.npc = "right", hjust = 1) +
  facet_wrap(~facet_label) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  scale_color_manual(name = "Stop codon", values = c("red", "orange", "forestgreen")) +
  xlab("Measured readthrough") + ylab("Predicted readthrough\n(reporter sequence)") +
  theme_bw(base_size = 7) + theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), 
                                   legend.position = "top")

p_fc_native <- ggplot(fc, aes(x = foldchange_measured_rt, y = foldchange_predicted_rt_native)) +
  geom_point(aes(color = stop_codon), alpha = 0.5) + 
  geom_text_repel(aes(label = Mutation, color = stop_codon), size = 6/.pt, show.legend = FALSE,
                  box.padding = 0.1, point.padding = 0.3, max.overlaps = Inf, min.segment.length = 0) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5, linewidth = 0.5, color = "grey50") +
  stat_cor(method = "spearman", cor.coef.name = "rho", size = 6/.pt, label.x.npc = "right", hjust = 1) +
  facet_wrap(~facet_label) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  scale_color_manual(name = "Stop codon", values = c("red", "orange", "forestgreen")) +
  xlab("Measured readthrough") + ylab("Predicted readthrough\n(native sequence)") +
  theme_bw(base_size = 7) + theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), 
                                   legend.position = "top")

# Combine panels
library(patchwork)
p_reporter <- (p_rt_reporter | p_fc_reporter) + plot_layout(widths = c(2, 1)) & theme(legend.position = "none")
p_compare <- (p_rt_compare | plot_spacer()) + plot_layout(widths = c(2.12, 0.86)) & theme(legend.position = "right")
p_native <- (p_rt_native | p_fc_native) + plot_layout(widths = c(2, 1)) & theme(legend.position = "none")
p <- (p_reporter / p_compare / p_native) + 
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
cairo_pdf(filename = "../Plots/Figure4.pdf", family = "Arial", width = 7, height = 7) 
p
dev.off()

# (half-column)
leg <- get_legend(p_rt_compare)
test <- p_rt_reporter / (p_fc_reporter | leg) / p_rt_compare / p_rt_native / (p_fc_native | plot_spacer()) + 
  #plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = list(c("a", "b", "", "c", "d", "e"))) &
  theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = "topleft", legend.position = "none")
cairo_pdf(filename = "../Plots/Figure4_test.pdf", family = "Arial", width = 3.4, height = 8.5) 
test
dev.off()
