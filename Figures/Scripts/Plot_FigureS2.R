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
allAG_reg <- lapply(allAG_merged, prep, cds_rpkm_cutoff = 5, utr3_rpkm_cutoff = 0.5, group_param = "none", feature_file = feature_file, utr3_region = "ext")

## List of readthough genes from Manjunath et al. (2022) and Loughran et al. (2018)
rt_genes <- data.frame(gene_name = c("HBB", "MPZ", "VEGFA", "AMD1", "AGO1", "ACP2", "MTCH2", # Other motif (Manjunath)
                                     "AQP4", "OPRL1", "OPRK1", "MAPK10", "LDHB", "MDH1", "VDR", # CUAG motif (Manjunath and Loughran)
                                     "ATP10D", "CDH23", "CDKN3", "CGGBP1", "DCTN3", "DDX58", "DUS4L", "GOT1L1", "HCLS1", "KCNB2", "MS4A5", "PHF19", "SIRPB1", "SPATA32", "TMEM86B", "ZFYVE20" # CUAG motif (additional from Loughran)
                                     ),
                       stop_codon = c("UGA", "UAG", rep("UGA", times = 28)),
                       motif = c("-", "", "hnRNPA2/B1", "", "Let-7a", "", "12nt", rep("CUAG", times = 23)))

allAG_reg_rt <- bind_rows(allAG_reg, .id = "Sample")
allAG_reg_rt <- left_join(allAG_reg_rt, rt_genes[, c("gene_name", "motif")], by = "gene_name")

allAG_reg_rt$Group <- "Other"
allAG_reg_rt[which(allAG_reg_rt$gene %in% rt_genes[which(rt_genes$motif == "CUAG"), ]$gene_name), ]$Group <- "SCR mRNAs: UGA CUAG motif"
allAG_reg_rt[which(allAG_reg_rt$gene %in% rt_genes[which(rt_genes$motif != "CUAG"), ]$gene_name), ]$Group <- "SCR mRNAs: Other motif"
allAG_reg_rt$Sample <- recode_factor(allAG_reg_rt$Sample, Untr = "Untreated", Amik = "Amikacin", 
                                     G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                                     Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")
allAG_reg_rt$Group <- factor(allAG_reg_rt$Group, levels = c("SCR mRNAs: UGA CUAG motif", "SCR mRNAs: Other motif", "Other"))

## Export data for Source Data
allAG_reg_rt <- allAG_reg_rt[, c("Sample", "transcript", "gene_name", "log_rte", "stop_codon", "stop_ntp04", "Group", "motif")]
write.csv(allAG_reg_rt, file = "../SourceData/SourceData_FigureS2.csv", quote = FALSE, row.names = FALSE)

## Density plot
plot_dist <- function(allAG_reg_rt, stop_subset = "all") {
  if(stop_subset == "all") {
    toplot <- allAG_reg_rt
    title <- "All mRNAs"
  } else {
    toplot <- allAG_reg_rt[which(allAG_reg_rt$stop_codon == stop_subset), ]
    stop <- gsub("T", "U", stop_subset)
    title <- paste(stop, " mRNAs")
  }
  toplot_summary <- toplot %>% group_by(Sample) %>% get_summary_stats(log_rte)
  p <- ggplot() + 
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
    scale_color_manual(values = c("#F8766D", "#00BA38")) +
    xlab("Readthrough efficiency") + ylab("Density") + 
    labs(subtitle = title) +
    theme_bw(base_size = 7) + 
    theme(panel.grid = element_blank(), panel.spacing = unit(0.1, "cm"), legend.position = "bottom", legend.text = element_text(size = 7),
          axis.text = element_text(size = 6.5),
          strip.background = element_rect(fill = "white"), strip.text = element_text(face = "italic", size = 7)) +
    guides(color = guide_legend(title = "", override.aes = list(alpha = 1)))
  return(p)
}

p_dist_all <- plot_dist(allAG_reg_rt, stop_subset = "all")
p_dist_uga <- plot_dist(allAG_reg_rt, stop_subset = "TGA")

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
