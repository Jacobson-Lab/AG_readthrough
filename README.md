# AG_readthrough

This repository contains a collection of codes used to train random forest models to predict stop codon readthrough efficiency derived from published ribosome prilfing data of aminoglycoside-treated cells by Wangen and Green, eLife (2020) (https://doi.org/10.7554/eLife.52611) from mRNA sequences features and use the model to predict new data. The paper detailing these results is under review.

## System Requirements
### Hardware Requirements
Building random forest models requires a standard computer with enough RAM or access to high performance computing cluster (HPCC). The rest of the analyses can be done locally on a standard computer.

### Software Requirements
Analyses were performed in R version 4.2 on a laptop with macOS Monterey 12.6.3 and R version 3.6 on high performance computing cluster (HPCC) linux system. However, all codes should run on any OS that can accommodate R version 3.6 or higher.

#### R packages used on R version 3.6 (on HPCC)
Random forest:
```
caret_6.0-86
randomForest_4.6-14
```
#### R packages used on R version 4.2 (on local computer)
General data handling:
```
readxl_1.4.0
dplyr_1.0.8
data.table_1.14.2
reshape2_1.4.4
```
Biological sequence data handling:
```
biomaRt_2.52.0
seqinr_4.2-8
Biostrings_2.64.0
```
Random forest:
```
caret_6.0-92
randomForest_4.7-1
```
Statistical analysis:
```
rstatix_0.7.0
```
Data visualization (plot and export):
```
ggplot2_3.4.0
ggpubr_0.4.0
ggrepel_0.9.1
ggh4x_0.2.3
scales_1.2.1
patchwork_1.1.1
Cairo_1.5-15
```
## Demo and expected output
*  Source data from Wangen and Green, eLife (2020) was too large to upload here, but it can be downloaded at  https://cdn.elifesciences.org/articles/52611/elife-52611-fig2-data1-v2.xlsx. Move/copy the file into _**Analysis scripts**_ folder to use with the codes there.
*  _**Analysis scripts**_ folder contains codes used to prepare data, create models, and use models for prediction as well as intermediate files (i.e., expected output at different stages of analyses) in Rdata format or csv/txt format for reference. All files should be in the same folder as appeared for ease of use.
*  _**Figures**_ folder contains tab- or comma-delimited files underlying each figure and the code to plot the figure.
