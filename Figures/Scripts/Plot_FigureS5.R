# ----------------------------------------
# Figure S5
# ----------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(rstatix)

dfu <- read.table("../Data/readthrough_3utr_length_effect.txt", header = TRUE)
dfu$Sample <- recode_factor(dfu$Sample, Untr = "Untreated", Amik = "Amikacin", 
                            G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                            Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")

## Readthrough efficiency vs. 3'-UTR length vs. MFE
dfu2 <- dfu[, c("Sample", "log_rte", "l_utr3", "MFE")]
colnames(dfu2) <- c("Sample", "Readthrough\nefficiency", "3'-UTR length\n(nt)", "MFE")
dfu2 <- split(dfu2, dfu2$Sample)
cormat <- lapply(dfu2, function(x) {
  cc <- cor(x[, 2:4], method = "spearman", use = "pairwise.complete.obs")
  return(cc)
})

  ### from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

cor_df <- bind_rows(lapply(cormat, function(x) reshape2::melt(get_upper_tri(x))), .id = "Sample")
cor_df$Sample <- recode_factor(cor_df$Sample, Untr = "Untreated", Amik = "Amikacin", 
                               G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                               Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")
p_mfe <- ggplot(cor_df) +
  geom_tile(aes(x = Var1, y = reorder(Var2, desc(Var2)), fill = value, height = 0.9, width = 0.9), color = NA) +
  facet_wrap(~Sample, nrow = 2) +
  scale_fill_gradient2(name = "Spearman's correlation coefficient   ", low = "blue", mid = "white", high = "red", midpoint = 0, na.value = NA, limits = c(-1, 1)) +
  xlab("") + ylab("") +
  geom_text(aes(x = Var1, y = reorder(Var2, desc(Var2)), label = round(value, 2)), color = "black", size = 6/.pt) +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"), strip.text = element_text(face = "italic"), 
        legend.position = "top", legend.key.height = unit(0.25, "cm"), legend.box.spacing = unit(0, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        aspect.ratio = 1)

## Readthrough efficiency vs. 3'-UTR length, grouped by stop codon and nt +4 identities
cor_all <-  dfu %>% group_by(Sample) %>% cor_test(log_rte, l_utr3, method = "spearman")
cor_stop4 <- dfu %>% group_by(Sample, stop_ntp04) %>% cor_test(log_rte, l_utr3, method = "spearman")

cor_all$stop_ntp04 <- "All"
cor_combine <- rbind(cor_all, cor_stop4)
cor_combine$psize <- ifelse(test = cor_combine$p < 0.05, yes = 0.9, no = 0.6)
cor_combine$Sample <- factor(cor_combine$Sample, levels = levels(dfu$Sample))
cor_combine$stop_ntp04 <- gsub("T", "U", cor_combine$stop_ntp04)

p_stop4 <- ggplot(cor_combine, aes(x = stop_ntp04, y = Sample, fill = cor, height = psize, width = psize)) + 
  geom_tile() + geom_text(aes(label = cor), size = 6/.pt) + 
  scale_fill_distiller(name = "Spearman's correlation coefficient   \nReadthrough efficiency vs. 3'-UTR length   ", 
                       palette = "Spectral", direction = -1, limits = c(0, 0.4)) + 
  scale_y_discrete(limits = rev(levels(dfu$Sample))) +
  xlab("Stop codon and nt +4") + ylab("") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), legend.position = "top", legend.title.align = 0.5, axis.text.y = element_text(face = "italic"), axis.title.y = element_blank()) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(5, "cm"), barheight = unit(0.25, "cm"), title.position = "left"))

# Combine panels
library(patchwork)
p <- (p_mfe / p_stop4) + 
  plot_layout(heights = c(4.5, 2)) +
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
cairo_pdf(filename = "../Plots/FigureS5.pdf", family = "Arial", width = 7, height = 6.5) 
p
dev.off()
