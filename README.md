# AG_readthrough

This repository contains a collection of codes used to train random forest models to predict stop codon readthrough efficiency derived from published ribosome prilfing data of aminoglycoside-treated cells (https://doi.org/10.7554/eLife.52611) from mRNA sequences features and use the model to predict new data. All analyses were done in the R programming environment. Random forest model training was carried out on high performance computing cluster. Others were done on RStudio on local computer. All codes use available, established R packages (see in codes and Methods). The paper detailing these results is under review.

_**Analysis scripts**_ folder contains codes used to prepare data, create models, and use models for prediction as well as intermediate files in Rdata format for reference. All files should be in the same folder as appeared for ease of use.

_**Figures**_ folder contains tab- or comma-delimited files underlying each figure and the code to plot the figure.
