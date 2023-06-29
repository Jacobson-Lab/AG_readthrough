# Properties of the nascent peptide in the ribosome's exit tunnel

load("amino_acid_properties.Rdata")

# Function to sum the number of matches and calculate fractions
sm <- function(x, match_string) {
  s <- sum(x == match_string, na.rm = TRUE)
  t <- sum(!is.na(x)) 
  f <- s/t
  return(f)
}

# Function to calculate fractions of properties
aa_prop <- function(seq_aa_m30, aa_range = 1:11, tunnel_part = "lower", feature_file) {
  aa_tab <- seq_aa_m30[, aa_range]
  x <- paste0("tunnel_", tunnel_part)
  # 1. Charge of side chain
  charge <- aa_tab
  map = setNames(aa_dict$charge, aa_dict$aa)
  charge[] <- map[unlist(charge)]
  charge_m <- apply(charge, 1, sm, match_string = "-1")
  charge_n <- apply(charge, 1, sm, match_string = "0")
  charge_p <- apply(charge, 1, sm, match_string = "1")
  charge_mnp <- data.frame(neg_charge = charge_m, no_charge = charge_n, pos_charge = charge_p)
  names(charge_mnp) <- c(paste0(x, "_neg_charge"), paste0(x, "_no_charge"), paste0(x, "_pos_charge"))
  feature_file <- cbind(feature_file, charge_mnp)
  rm(map, charge, charge_m, charge_n, charge_p, charge_mnp)
  # 2. Aromaticity
  aroma <- aa_tab
  map = setNames(aa_dict$aromaticity, aa_dict$aa)
  aroma[] <- map[unlist(aroma)]
  aroma_a <- data.frame(apply(aroma, 1, sm, match_string = "Aromatic"))
  names(aroma_a) <- c(paste0(x, "_aromatic"))
  feature_file <- cbind(feature_file, aroma_a)
  rm(map, aroma, aroma_a)
  # 3. Polarity
  polarity <- aa_tab
  map = setNames(aa_dict$polarity, aa_dict$aa)
  polarity[] <- map[unlist(polarity)]
  polar <- apply(polarity, 1, sm, match_string = "Polar")
  nonpolar <- apply(polarity, 1, sm, match_string = "Nonpolar")
  polar_all <- data.frame(p = polar, n = nonpolar)
  names(polar_all) <- c(paste0(x, "_polar"), paste0(x, "_nonpolar"))
  feature_file <- cbind(feature_file, polar_all)
  rm(map, polarity, polar, nonpolar, polar_all)
  # 4. Hydrophobicity
  hb <- aa_tab
  map = setNames(aa_dict$hydrophobicity, aa_dict$aa)
  hb[] <- map[unlist(hb)]
  hb_vhb <- apply(hb, 1, sm, match_string = "Very hydrophobic")
  hb_hb <- apply(hb, 1, sm, match_string = "Hydrophobic")
  hb_n <- apply(hb, 1, sm, match_string = "Neutral")
  hb_hp <- apply(hb, 1, sm, match_string = "Hydrophilic")
  hb_all <- data.frame(v_hydrophobic = hb_vhb, hydrophobic = hb_hb, neutral = hb_n, hydrophilic = hb_hp)
  names(hb_all) <- c(paste0(x, "_v_hydrophobic"), paste0(x, "_hydrophobic"), paste0(x, "_neutral"), paste0(x, "_hydrophilic"))
  feature_file <- cbind(feature_file, hb_all)
  rm(map)
  rm(list = ls(pattern  = "^hb*"))
  return(feature_file)
}