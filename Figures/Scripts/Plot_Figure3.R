# ----------------------------------------
# Figure 3
# ----------------------------------------

library(dplyr)
library(ggplot2)
#library(ggpubr)

dfu <- read.table("../Data/readthrough_3utr_length_effect.txt", header = TRUE)
dfu$Sample <- recode_factor(dfu$Sample, Untr = "Untreated", Amik = "Amikacin", 
                            G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                            Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")

cor_df <-  dfu %>% group_by(Sample) %>% cor_test(log_rte, l_utr3, method = "spearman")
p_all <- ggplot(dfu, aes(x = l_utr3, y = log_rte)) +
  geom_point(alpha = 0.3, color = "black", size = 0.25) +
  geom_smooth(method = "lm", formula = y~x) +
  #stat_cor(method = "spearman", geom = "text", label.x.npc = "middle", label.y = 5.5, hjust = 0.5, vjust = 0.5, size = 6/.pt, cor.coef.name = "rho", color = "blue") + # Does not give exact p-value for p < 2.2e-16
  geom_text(data = cor_df, aes(x = 400, y = 5.5, label = paste0("rho", "==", cor)), 
            parse = TRUE, hjust = 1, vjust = 0.5, size = 6/.pt, color = "blue") +
  geom_text(data = cor_df, aes(x = 400, y = 5.5, label = ", "), 
            parse = FALSE, hjust = 0, vjust = 0.5, size = 6/.pt, color = "blue") +
  geom_text(data = cor_df, aes(x = 400, y = 5.5, label = paste0("italic('p')==", p)), 
            parse = TRUE, hjust = -0.1, vjust = 0.5, size = 6/.pt, color = "blue") +
  facet_wrap(~Sample, nrow = 2) +
  scale_x_log10() +
  xlab("3'-UTR length (nt)") + ylab("Readthrough Efficiency") + 
  coord_cartesian(ylim = c(-10.5, 6.5)) +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing.x = unit(0.2,"cm"), panel.spacing.y = unit(0.2, "cm"), axis.text = element_text(size = 6.5),
        strip.text = element_text(face = "italic", size = 7), strip.background = element_rect(fill = "white"), aspect.ratio = 1)

library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)

cairo_pdf(filename = "../Plots/Figure3.pdf", family = "Arial", width = 7, height = 3.5) 
p_all
dev.off()
