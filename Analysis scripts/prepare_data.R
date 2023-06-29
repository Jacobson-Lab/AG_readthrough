# Obtain data and add mRNA features
## Load required libraries
library(readxl)
library(Biostrings)
library(dplyr)
library(data.table)

## Get file containing readthrough information from Wangen and Green, eLife (2020) ----------------------------------------------------
  ## Download data from https://cdn.elifesciences.org/articles/52611/elife-52611-fig2-data1-v2.xlsx
xpath <- "elife-52611-fig2-data1-v2.xlsx"
xsheet <- excel_sheets(xpath)
xsheet <- xsheet[c(5:13)] # Choose sheets to load

allAG_merged <- list()
for (i in 1:length(xsheet)) {
  dat <- as.data.frame(read_excel(path = xpath, sheet = xsheet[i]))
  dat <- dat[, c(1, 5:8, 21:23, 29, 30)]
  colnames(dat) <- c("transcript", "l_tr", "l_cds", "l_utr5", "l_utr3", "rpkm_cds", "rpkm_utr5", "rpkm_utr3", "rpkm_ext", "rte_ext")
  # Fix lengths of mRNA regions to real lengths (Wangen had adjusted it to the need of their analysis)
  dat$l_utr5 <- dat$l_utr5+6
  dat$l_cds <- dat$l_cds+30
  dat$l_utr3 <- dat$l_utr3+9
  dat$l_tr <- dat$l_tr+6+30+9
  allAG_merged[[xsheet[i]]] <- data.table(dat)
}
names(allAG_merged) <- sub(" merged", "", names(allAG_merged))
names(allAG_merged)[c(1, 2, 8, 9)] <- c("Untr", "G418_500", "G418_2000", "G418_500_10min")
save(allAG_merged, file = "allAG_merged.Rdata")

## Get cDNA sequences FASTA from ensembl ----------------------------------------------------------------------------------------------
  ## Find unique ENST among all samples
ut <- unique(bind_rows(lapply(allAG_merged, function(x) x[, 1:5])))
  ## Get FASTA sequence
library(biomaRt)
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
seqq <- getSequence(id = ut$transcript, type = "ensembl_transcript_id_version", seqType = "cdna", mart = ensembl)
colnames(seqq) <- c("sequence", "transcript")

## Record mRNA features ---------------------------------------------------------------------------------------------------------------
feature_file <- left_join(ut, seqq, by = "transcript")
feature_file <- feature_file[!is.na(feature_file$sequence), ] # Remove rows that the sequence is missing
sum(feature_file$l_tr != nchar(feature_file$sequence)) # Check that number of nt from sequence column is equal to transcript length (l_tr). Sum of unequal numbers should be zero

### Random numbers and factors (negative controls)
feature_file$random_num <- sample(100, size = nrow(feature_file), replace = TRUE)
feature_file$random_factor <- sample(c("A", "C", "G", "T"), size = nrow(feature_file), replace = TRUE)

### Identity of the stop codon
feature_file$stop_codon <- substr(feature_file$sequence, 
                                  start = feature_file$l_tr - feature_file$l_utr3 + 1, 
                                  stop = feature_file$l_tr - feature_file$l_utr3 + 3) # stop codon is considered part of 3'UTR

### Identities of nucleotides after the stop codon (stop codon is considered +1, +2, +3)
for (i in 4:16) {
  feature_file[, paste0("nt_p", sprintf("%02d", i))] <- substr(feature_file$sequence, 
                                                               start = feature_file$l_tr - feature_file$l_utr3 + i, 
                                                               stop = feature_file$l_tr - feature_file$l_utr3 + i)
}

### Identities of nucleotides priot to the stop codon (stop codon is considered +1, +2, +3)
seq_cds <- substr(feature_file$sequence, start = feature_file$l_utr5 + 1, stop = feature_file$l_tr - feature_file$l_utr3)
for (i in 1:18) {
  feature_file[, paste0("nt_m", sprintf("%02d", i))] <- substr(seq_cds, 
                                                               start = feature_file$l_cds - i + 1, 
                                                               stop = feature_file$l_cds - i + 1)
}

### Identity of P-site codon/aa (position -1)
feature_file$codon_m01 <- substr(feature_file$sequence, start = feature_file$l_tr - feature_file$l_utr3 -2, stop = feature_file$l_tr - feature_file$l_utr3) 
feature_file <- cbind(feature_file, aa_m01 = apply(t(apply((feature_file[, "codon_m01"]), 1, FUN = seqinr::s2c)), 1, FUN = seqinr::translate))

### Identity of E-site codon/aa (position -2)
feature_file$codon_m02 <- substr(feature_file$sequence, start = feature_file$l_tr - feature_file$l_utr3 -5, stop = feature_file$l_tr - feature_file$l_utr3 -3) 
feature_file <- cbind(feature_file, aa_m02 = apply(t(apply((feature_file[, "codon_m02"]), 1, FUN = seqinr::s2c)), 1, FUN = seqinr::translate))

### Properties of the nascent peptide in the ribosome's exit tunnel (excluding P- and E-site codons)
source("aa_prop.R")
seq_cds <- substr(feature_file$sequence, start = feature_file$l_utr5 + 1, stop = feature_file$l_tr - feature_file$l_utr3)
seq_codon_m30 <- substr(seq_cds, start = feature_file$l_cds - 89, stop = feature_file$l_cds) # First nt/codon is the furthest from stop
seq_codon_m30 <- seqinr::as.SeqFastadna(seq_codon_m30)
seq_aa_m30 <- data.frame(matrix(nrow = length(seq_codon_m30), ncol = 30))
for (i in 1:length(seq_codon_m30)) { # For mRNA that's shorter than 30 aa, make data "right aligned"
  aa <- seqinr::getTrans(seqinr::s2c(seq_codon_m30[i]))
  seq_aa_m30[i, (30-length(aa)+1):30] <- aa
}
feature_file <- aa_prop(seq_aa_m30 = seq_aa_m30, aa_range = 1:11, tunnel_part = "lower", feature_file = feature_file) # Near exit tunnel
feature_file <- aa_prop(seq_aa_m30 = seq_aa_m30, aa_range = 12:18, tunnel_part = "central", feature_file = feature_file)
feature_file <- aa_prop(seq_aa_m30 = seq_aa_m30, aa_range = 19:21, tunnel_part = "constriction", feature_file = feature_file)
feature_file <- aa_prop(seq_aa_m30 = seq_aa_m30, aa_range = 22:28, tunnel_part = "upper", feature_file = feature_file) # Near peptidyl transferase center (excluding excluding P- and E-site codons)

### Save file
save(feature_file, file = "feature_file_wangen_elife_2020.Rdata")
