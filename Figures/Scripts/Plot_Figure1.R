# ----------------------------------------
# Figure 1
# ----------------------------------------

library(dplyr)
library(ggplot2)
library(ggh4x)

# A
imp_reg <- read.table("../Data/rf_reg_feature_importance.txt", header = TRUE)
imp_reg <- imp_reg[, c(1, 2, 4)]

# B
imp_class <- read.table("../Data/rf_class_feature_importance.txt", header = TRUE)
imp_class <- imp_class[, c(1, 4, 6)]

# Combine A and B data
all.imp_df <- full_join(imp_reg, imp_class, by = c("Sample", "feature"))

# Prepare for plot
all.imp_df$xtick <- as.factor(all.imp_df$feature) # for x-axis label
all.imp_df$xtick <- recode_factor(all.imp_df$xtick,
                                  tunnel_lower_H_count = "Helical", tunnel_lower_neg_charge = "- charged", tunnel_lower_no_charge = "No charge", tunnel_lower_pos_charge = "+ charged", tunnel_lower_aromatic = "Aromatic", tunnel_lower_polar = "Polar", tunnel_lower_nonpolar = "Nonpolar", tunnel_lower_hydrophilic = "Hydrophylic", tunnel_lower_neutral = "Neutral", tunnel_lower_hydrophobic = "Hydrophobic", tunnel_lower_v_hydrophobic = "Very hydrophobic", 
                                  tunnel_central_neg_charge = " - charged", tunnel_central_no_charge = " No charge", tunnel_central_pos_charge = " + charged", tunnel_central_aromatic = " Aromatic", tunnel_central_polar = " Polar", tunnel_central_nonpolar = " Nonpolar", tunnel_central_hydrophilic = " Hydrophylic", tunnel_central_neutral = " Neutral", tunnel_central_hydrophobic = " Hydrophobic", tunnel_central_v_hydrophobic = " Very hydrophobic",
                                  tunnel_constriction_neg_charge = "  - charged", tunnel_constriction_no_charge = "  No charge", tunnel_constriction_pos_charge = "  + charged", tunnel_constriction_aromatic = "  Aromatic", tunnel_constriction_polar = "  Polar", tunnel_constriction_nonpolar = "  Nonpolar", tunnel_constriction_hydrophilic = "  Hydrophylic", tunnel_constriction_neutral = "  Neutral", tunnel_constriction_hydrophobic = "  Hydrophobic", tunnel_constriction_v_hydrophobic = "  Very hydrophobic",
                                  tunnel_upper_neg_charge = "   - charged", tunnel_upper_no_charge = "   No charge", tunnel_upper_pos_charge = "   + charged", tunnel_upper_aromatic = "   Aromatic", tunnel_upper_polar = "   Polar", tunnel_upper_nonpolar = "   Nonpolar", tunnel_upper_hydrophilic = "   Hydrophylic", tunnel_upper_neutral = "   Neutral", tunnel_upper_hydrophobic = "   Hydrophobic", tunnel_upper_v_hydrophobic = "   Very hydrophobic",
                                  aa_m02 = "E-site aa", aa_m01 = "P-site aa", nt_m18 = "-18", nt_m17 = "-17", nt_m16 = "-16", nt_m15 = "-15", nt_m14 = "-14", nt_m13 = "-13", nt_m12 = "-12", nt_m11 = "-11", nt_m10 = "-10", nt_m09 = "-9", nt_m08 = "-8", nt_m07 = "-7", nt_m06 = "-6", nt_m05 = "-5", nt_m04 = "-4", nt_m03 = "-3", nt_m02 = "-2", nt_m01 = "-1", stop_codon = "Stop codon", nt_p04 = "+4", nt_p05 = "+5", nt_p06 = "+6", nt_p07 = "+7", nt_p08 = "+8", nt_p09 = "+9", nt_p10 = "+10", nt_p11 = "+11", nt_p12 = "+12", nt_p13 = "+13", nt_p14 = "+14", nt_p15 = "+15", nt_p16 = "+16",
                                  nis_stop = "1st 3'UTR stop", l_utr3 = "3'-UTR length", MFE = "MFE", dist_bp = "Distance from stop", random_factor = "Random factor", random_num = "Random number")

all.imp_df$xgroup <- NA # for faceting
all.imp_df$xgroup <- gsub("tunnel_lower_.*", "aa 20-30\nfrom PTC", all.imp_df$feature)
all.imp_df$xgroup <- gsub("tunnel_central_.*", "aa 13-19\nfrom PTC", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("tunnel_constriction_.*", "aa 10-12\nfrom PTC", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("tunnel_upper_.*", "aa 3-9\nfrom PTC", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("aa_m.*", "", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("nt_m.*", " nt from stop ", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("stop_codon", " ", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("nt_p.*", "nt from stop", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("l_utr3", "  ", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("nis_stop", "   ", all.imp_df$xgroup)
all.imp_df$xgroup <- gsub("random_.*", "NC", all.imp_df$xgroup)
all.imp_df$xgroup <- factor(all.imp_df$xgroup, levels = c("aa 20-30\nfrom PTC", "aa 13-19\nfrom PTC", "aa 10-12\nfrom PTC", "aa 3-9\nfrom PTC", "", " nt from stop ", " ", "nt from stop", "  ", "   ", "2°", "NC"))

all.imp_df$xgroup2 <- NA # for nested faceting
all.imp_df[which(grepl(pattern = "aa", x = all.imp_df$xgroup)), ]$xgroup2 <- "Nascent peptide in the exit tunnel"
all.imp_df[which(grepl(pattern = "aa", x = all.imp_df$xtick)), ]$xgroup2 <- "Nascent peptide in the exit tunnel"
all.imp_df[which(grepl(pattern = "nt_m", x = all.imp_df$feature)), ]$xgroup2 <- ""
all.imp_df[which(grepl(pattern = "stop_codon", x = all.imp_df$feature)), ]$xgroup2 <- " "
all.imp_df[which(grepl(pattern = "nt_p", x = all.imp_df$feature)), ]$xgroup2 <- "  "
all.imp_df[which(grepl(pattern = "l_utr3", x = all.imp_df$feature)), ]$xgroup2 <- "   "
all.imp_df[which(grepl(pattern = "NC", x = all.imp_df$xgroup)), ]$xgroup2 <- "    "
all.imp_df$xgroup2 <- factor(all.imp_df$xgroup2, levels = c("Nascent peptide in the exit tunnel", "", " ", "  ", "   ", "    "))

all.imp_df$Sample <- recode_factor(all.imp_df$Sample, Untr = "Untreated", Amik = "Amikacin", 
                                   G418_2000 = "G418 (2)", G418_500 = "G418 (0.5)", G418_500_10min = "G418 (0.5, 10 min)",  
                                   Genta = "Gentamicin", Neo = "Neomycin", Parom = "Paromomycin", Tobra = "Tobramycin")
# Plot
A <- ggplot() +
  geom_tile(data = all.imp_df, 
            mapping = aes(x = xtick, y = Sample, fill = X.IncMSE)) +
  facet_nested(.~xgroup2+xgroup, scales = "free", space = "free", switch = "y", nest_line = TRUE) +
  scale_fill_gradient2(name = "%IncMSE  ", low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-15, 15), oob = scales::squish) +
  scale_y_discrete(limits = rev(levels(all.imp_df$Sample))) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.25, "cm"), 
        legend.box.spacing = unit(0.5, "cm"), legend.title.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(face = "italic"), 
        axis.title = element_blank(),
        strip.background = element_blank(), strip.text.y = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA), 
        panel.grid = element_blank(), panel.spacing = unit(0,"cm"))
B <- ggplot() +
  geom_tile(data = all.imp_df, 
            mapping = aes(x = xtick, y = Sample, fill = MeanDecreaseAccuracy)) +
  facet_nested(.~xgroup2+xgroup, scales = "free", space = "free", switch = "y", nest_line = TRUE) +
  scale_fill_gradient2(name = "MDA  ", low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-15, 15), oob = scales::squish) +
  scale_y_discrete(limits = rev(levels(all.imp_df$Sample))) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.25, "cm"), 
        legend.box.spacing = unit(0.5, "cm"), legend.title.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(face = "italic"), 
        axis.title = element_blank(),
        strip.background = element_blank(), strip.text.y = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA), 
        panel.grid = element_blank(), panel.spacing = unit(0,"cm"))

# Combine panels
library(patchwork)
p <- (A / B) + 
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
cairo_pdf(filename = "Figure1.pdf", family = "Arial", width = 8, height = 7) 
p
dev.off()
