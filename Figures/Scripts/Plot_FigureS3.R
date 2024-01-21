# ----------------------------------------
# Figure S3
# ----------------------------------------

library(dplyr)
library(ggplot2)

stop4_chisq <- read.table("../Data/chisq_proportion_stop4.txt", header = TRUE)
stop4_chisq <- stop4_chisq %>% group_by(Sample, stop_codon) %>% mutate(total_by_stop = sum(observed_count))
stop4_chisq$expected_fraction_by_stop <- stop4_chisq$expected_count / stop4_chisq$total_by_stop
stop4_chisq$observed_fraction_by_stop <- stop4_chisq$observed_count / stop4_chisq$total_by_stop
stop4_chisq$Sample <- recode_factor(stop4_chisq$Sample, Untr = "Untreated", Amik = "Amikacin", 
                                    G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                                    Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin", 
                                    Reference = "Reference")

# Number of mRNAs in each group
dfn <- stop4_chisq %>% group_by(Sample) %>% summarise(n = sum(observed_count))
p_n <- ggplot(dfn, aes(x = n, y = Sample)) + 
  geom_col(fill = "#999999") + geom_text(aes(label = n), hjust = 0, size = 6/.pt) +
  scale_y_discrete(limits = rev(levels(stop4_chisq$Sample))) +
  coord_cartesian(xlim = c(0, 15000)) +
  xlab("Number of mRNAs") + ylab("") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic"))

# Observed frequency of stop codon and nt +4 together
p_stop4_ob <- ggplot(stop4_chisq) +
  geom_tile(aes(x = nt_p04, y = Sample, fill = observed_fraction*100), height = 0.9, width = 0.9) +
  facet_grid(.~stop_codon) +
  scale_fill_distiller(name = "Observed frequency (%) of\nstop codon and nt +4", palette = "Spectral", limit = c(2, 21)) +
  scale_y_discrete(limits = rev(levels(stop4_chisq$Sample))) +
  xlab("nt +4") + ylab("") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"),
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(5, "cm"), barheight = unit(0.25, "cm"), 
                                title.position = "top"))

# Expected frequency of stop codon and nt +4 together
p_stop4_ex <- ggplot(stop4_chisq) +
  geom_tile(aes(x = nt_p04, y = Sample, fill = expected_fraction*100), height = 0.9, width = 0.9) +
  facet_grid(.~stop_codon) +
  scale_fill_distiller(name = "Expected frequency (%) of\nstop codon and nt +4", palette = "Spectral", limit = c(2, 21)) +
  scale_y_discrete(limits = rev(levels(stop4_chisq$Sample))) +
  xlab("nt +4") + ylab("") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"),
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(5, "cm"), barheight = unit(0.25, "cm"), 
                                title.position = "top"))

# Observed frequency of stop codon and nt +4, separately
df_stop <- stop4_chisq %>% group_by(Sample, stop_codon) %>% summarise(count = sum(observed_count))
df_stop <- df_stop %>% group_by(Sample) %>% mutate(fraction = count / sum(count))
colnames(df_stop)[2] <- "position"
df_stop$feature <- "Stop codon"
df_nt4 <- stop4_chisq %>% group_by(Sample, nt_p04) %>% summarise(count = sum(observed_count))
df_nt4 <- df_nt4 %>% group_by(Sample) %>% mutate(fraction = count / sum(count))
colnames(df_nt4)[2] <- "position"
df_nt4$feature <- "nt +4"
df <- rbind(df_stop, df_nt4)
df$feature <- factor(df$feature, levels = c("Stop codon", "nt +4"))

p_sep <- ggplot(df) +
  geom_tile(aes(x = position, y = Sample, fill = fraction*100), height = 0.9, width = 0.9) +
  facet_grid(.~feature, scales = "free_x", space = "free_x") +
  scale_fill_distiller(name = "Observed frequency (%)", palette = "Spectral", limit = c(10, 60)) +
  scale_y_discrete(limits = rev(levels(df_stop$Sample))) +
  xlab("") + ylab("") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"),
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"), 
                                title.position = "top"))

# Observed frequency, grouped by stop codon identity
p_ob <- ggplot(stop4_chisq) +
  geom_tile(aes(x = nt_p04, y = Sample, fill = observed_fraction_by_stop*100), height = 0.9, width = 0.9) +
  facet_grid(.~stop_codon) +
  scale_fill_distiller(name = "Observed frequency (%) of\nnt +4 for each stop codon", palette = "Spectral", limits = c(10, 45)) +
  scale_y_discrete(limits = rev(levels(stop4_chisq$Sample))) +
  xlab("nt +4") + ylab("") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), legend.title.align = 0.5,
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(5, "cm"), barheight = unit(0.25, "cm"), 
                                title.position = "top"))

# Expected frequency, grouped by stop codon identity
p_ex <- ggplot(stop4_chisq) +
  geom_tile(aes(x = nt_p04, y = Sample, fill = expected_fraction_by_stop*100), height = 0.9, width = 0.9) +
  facet_grid(.~stop_codon) +
  scale_fill_distiller(name = "Expected frequency (%) of\nnt +4 for each stop codon", palette = "Spectral", limits = c(10, 45)) +
  scale_y_discrete(limits = rev(levels(stop4_chisq$Sample))) +
  xlab("nt +4") + ylab("") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "cm"), panel.border = element_rect(size = 0.25),
        legend.position = "top", legend.box = "horizontal", legend.title = element_text(hjust = 0.5, vjust = 1), legend.margin = margin(r = 1, unit = "cm"), 
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "white", size = 0.25)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(5, "cm"), barheight = unit(0.25, "cm"), 
                                title.position = "top"))

library(patchwork)
p <- (p_n | p_sep) / (p_stop4_ob | p_ob) / (p_stop4_ex | p_ex) +
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
cairo_pdf(filename = "../Plots/FigureS3.pdf", family = "Arial", width = 7, height = 8.5) 
p
dev.off()
