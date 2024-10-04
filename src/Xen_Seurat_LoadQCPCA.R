# Script 1 of Xenium Seurat Transcript Processing Pipeline
# Last updated: 04Oct2024
# Author: jrose

# This script covers:
#       1. Reading in Xenium data from onboard analysis output folders
#       2. Reading in additional cell data (area, MFI, etc) from additional sources
#       3. Performing cell-level QC
#       4. Integrating multiple samples from the same experiment into one SeuratObject
#       5. Count Normalization
#       6. Dimensionaltiy reduction (PCA)



#################################################


library(Seurat)
library(tidyverse)
library(future)
library(clustree)
library(spacexr)
library(scBubbletree)
library(here)
library(DESeq2)
library(renv)

