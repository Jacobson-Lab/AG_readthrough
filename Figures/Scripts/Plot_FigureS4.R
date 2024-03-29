# ----------------------------------------
# Figure S4
# ----------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)

dfu <- read.table("../Data/readthrough_3utr_length_effect.txt", header = TRUE)
dfu$Sample <- recode_factor(dfu$Sample, Untr = "Untreated", Amik = "Amikacin", 
                            G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                            Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")
# A
dfu_filter <- dfu[which(dfu$l_utr3 >= 100 & dfu$l_utr3 <= 5000), ]
cor_df <-  dfu_filter %>% group_by(Sample) %>% cor_test(log_rte, l_utr3, method = "spearman")
p_filter <- ggplot(dfu_filter, aes(x = l_utr3, y = log_rte)) +
  geom_point(alpha = 0.3, color = "black", size = 0.25) +
  geom_smooth(method = "lm", formula = y~x) +
  #stat_cor(method = "spearman", geom = "text", label.x.npc = "middle", label.y = 5.5, hjust = 0.5, vjust = 0.5, size = 6/.pt, cor.coef.name = "rho", color = "blue") +
  geom_text(data = cor_df, aes(x = 550, y = 5.5, label = paste0("rho", "==", cor)), 
            parse = TRUE, hjust = 1, vjust = 0.5, size = 6/.pt, color = "blue") +
  geom_text(data = cor_df, aes(x = 550, y = 5.5, label = ", "), 
            parse = FALSE, hjust = 0, vjust = 0.5, size = 6/.pt, color = "blue") +
  geom_text(data = cor_df, aes(x = 550, y = 5.5, label = paste0("italic('p')==", p)), 
            parse = TRUE, hjust = -0.1, vjust = 0.5, size = 6/.pt, color = "blue") +
  facet_wrap(~Sample, nrow = 2) +
  scale_x_log10() +
  xlab("3'-UTR length (nt)") + ylab("Readthrough Efficiency") + 
  coord_cartesian(ylim = c(-10.5, 6.5)) +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing.x = unit(0.2,"cm"), panel.spacing.y = unit(0.2, "cm"),axis.text = element_text(size = 6.5),
        strip.text = element_text(face = "italic", size = 7), strip.background = element_rect(fill = "white"), aspect.ratio = 1)

# B
dfu_untr <- dfu[which(dfu$Sample == "Untreated"), ]

## Yeast data
dfy <- read.table("https://raw.githubusercontent.com/Jacobson-Lab/sup45-ts_readthrough/main/Figures/Data/Data_Figure8.txt", header = TRUE)
dfy_wt <- dfy[which(dfy$Sample == "sup45_wt_25C"), ]

## Combine data of two organisms
dfw <- data.frame(organism = c(rep("Yeast", times = nrow(dfy_wt)), rep("HEK293T", times = nrow(dfu_untr))), 
                  l_utr3 = c(dfy_wt$value, dfu_untr$l_utr3))

pl <- ggplot(dfw, aes(x = l_utr3, fill = organism)) + 
  geom_histogram(position = "dodge") + 
  scale_x_log10() + 
  scale_fill_manual(values = c(HEK293T = "coral1", Yeast = "goldenrod1"), name = "") + 
  xlab("3'-UTR length (nt)") + ylab("Count") + 
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), axis.text = element_text(size = 6.5), legend.text = element_text(size = 7), 
        legend.position = c(0.8, 0.8), legend.title = element_blank(), legend.background = element_blank())

# Combine panels
library(patchwork)
pb <- pl + plot_spacer() + plot_layout(nrow = 1, widths = c(0.5, 0.5))
p <- (p_filter / pb) + 
  plot_layout(heights = c(2, 1)) +
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
cairo_pdf(filename = "../Plots/FigureS4.pdf", family = "Arial", width = 7, height = 5.5) 
p
dev.off()
