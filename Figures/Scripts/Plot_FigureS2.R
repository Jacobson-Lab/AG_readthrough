# ----------------------------------------
# Figure S2
# ----------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(rstatix)

load("../../Analysis scripts/allAG_merged.Rdata") # From prepare_data.R
load("../../Analysis scripts/feature_file_wangen_elife_2020.Rdata") # From prepare_data.R
source("../../Analysis scripts/functions_analyses.R")
feature_file$stop_ntp04 <- paste0(feature_file$stop_codon, feature_file$nt_p04)
allAG_reg <- lapply(allAG_merged, prep, cds_rpkm_cutoff = 5, utr3_rpkm_cutoff = 0.5, group_param = "none", feature_file = feature_file, utr3_region = "ext")

## List of readthough genes from Manjunath et al. (2022)
rt_genes <- data.frame(gene_name = c("HBB", "MPZ", "VEGFA", "AQP4", "OPRL1", "OPRK1", "MAPK10", "LDHB", "MDH1", "VDR", "AMD1", "AGO1", "ACP2", "MTCH2"),
                       stop_codon = c("UGA", "UAG", rep("UGA", times = 12)),
                       motif = c("-", "", "hnRNPA2/B1", rep("CUAG", times = 7), "", "Let-7a", "", "12nt"))

allAG_reg_rt <- bind_rows(allAG_reg, .id = "Sample")
allAG_reg_rt <- left_join(allAG_reg_rt, rt_genes[, c("gene_name", "motif")], by = "gene_name")

allAG_reg_rt$Group <- "Other"
allAG_reg_rt[which(allAG_reg_rt$stop_codon == "TGA" & allAG_reg_rt$nt_p04 == "C" & allAG_reg_rt$nt_p05 == "T" & allAG_reg_rt$nt_p06 == "A" & allAG_reg_rt$nt_p07 == "G"), ]$Group <- "Non-SCR mRNAs: UGA CUAG motif"
allAG_reg_rt[which(allAG_reg_rt$gene %in% c("AQP4", "OPRL1", "OPRK1", "MAPK10", "LDHB", "MDH1", "VDR")), ]$Group <- "SCR mRNAs: UGA CUAG motif"
allAG_reg_rt[which(allAG_reg_rt$gene %in% c("HBB", "MPZ", "VEGFA", "AMD1", "AGO1", "ACP2", "MTCH2")), ]$Group <- "SCR mRNAs: Other motif"
allAG_reg_rt$Sample <- recode_factor(allAG_reg_rt$Sample, Untr = "Untreated", Amik = "Amikacin", 
                                     G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                                     Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")
allAG_reg_rt$Group <- factor(allAG_reg_rt$Group, levels = c("SCR mRNAs: UGA CUAG motif", "SCR mRNAs: Other motif", "Non-SCR mRNAs: UGA CUAG motif", "Selenoprotein", "Other"))

## Density plot
toplot <- allAG_reg_rt
#toplot <- allAG_reg_rt[which(allAG_reg_rt$stop_codon == "TGA"), ]
toplot_summary <- toplot %>% group_by(Sample) %>% get_summary_stats(log_rte)
p_dist_all <- ggplot() + 
  geom_density(data = toplot, 
               aes(x = log_rte), alpha = 0.1, fill = "grey80") +
  geom_point(data = toplot[which(!(toplot$Group %in% c("Other"))), ], 
             aes(y = 0, x = log_rte, color = Group), alpha = 0.5) + 
  geom_label_repel(data = toplot[which(!(toplot$Group %in% c("Other"))), ],
                   aes(y = 0, x = log_rte, label = gene_name, color = Group), 
                   size = 6/.pt, show.legend = FALSE, nudge_y = 0.1,
                   #direction = "y", 
                   label.padding = 0.1, box.padding = 0.1, point.padding = 0.1, max.overlaps = Inf, min.segment.length = 0) + 
  geom_vline(data = toplot_summary, aes(xintercept = median), linetype = "dashed", color = "grey30", alpha = 0.5) +
  geom_vline(data = toplot_summary, aes(xintercept = mean), linetype = "dashed", color = "brown", alpha = 0.5) +
  facet_wrap(~Sample, nrow = 2) +
  xlab("Readthrough efficiency") + ylab("Density") + 
  labs(subtitle = "All mRNAs") +
  #labs(subtitle = "UGA mRNAs") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0.1, "cm"), legend.position = "bottom",
        strip.background = element_rect(fill = "white"), strip.text = element_text(face = "italic")) +
  guides(color = guide_legend(title = "", override.aes = list(alpha = 1)))

# Combine panels
library(patchwork)
p_dist <- (p_dist_all / p_dist_uga) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = "topleft", legend.position = "bottom")

# Export plot
library(Cairo)
CairoFonts(
  regular = "Arial:style=Regular",
  bold = "Arial:style=Black",
  italic = "Arial:style=Italic",
  bolditalic = "Arial:style=Black Italic",
  symbol = "Symbol"
)
cairo_pdf(filename = "../Plots/FigureS2.pdf", family = "Arial", width = 7, height = 8.5) 
p_dist
dev.off()
