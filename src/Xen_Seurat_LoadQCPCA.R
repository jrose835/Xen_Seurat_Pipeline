# Script 1 of Santangelo lab Xenium Seurat Transcript Processing Pipeline

# This script covers:
#       1. Reading in Xenium data from onboard analysis output folders
#       2. Reading in additional cell data (area, MFI, etc) from additional sources
#       3. Performing cell-level QC
#       4. Integrating multiple samples from the same experiment into one SeuratObject
#       5. Count Normalization
#       6. Dimensionaltiy reduction (PCA)

# Description: This script performs the initial preprocessing steps for processing transcipt data from 10x Xenium experiments using the Seurat package. It is meant to be paired with Xen_Suerat_ClustAnnot.R for complete analysis

# Usage:

# Inputs:
#     - Manifest file (.csv) with 'run_name', 'experiment', 'condition', 'xen_dir', 'alt_input', and 'alt_input_dir' columns

# Outputs:

# Last updated: 04Oct2024
# Author: jrose
#################################################
#Set up

library(Seurat)
library(tidyverse)
library(future)
library(clustree)
library(spacexr)
library(scBubbletree)
library(here)
source(here("src", "Xen_Seurat_functions.R"))  #accompanying custom functions

set.seed(1984)

experiment_name = "Abbie_Lung"

man_file <- "sample_manifest.csv"
manifest <- read_csv(here(man_file)) # Get run manifest

outdir <- "output"
prepOutDir(outdir=outdir, manifest=manifest) # Prepare output dir structure. Custom func

#################################################
#Module 1 load data and generate basic qc stats

objs <- LoadMultiXenium(manifest) # Custom LoadXenium() wrapper, resturns list of SeuratObjects

cellstats <- do.call(rbind, lapply(objs, XenSeuratStats)) #Summary statistics for each run. Custom XenSeuratStats func
write.csv(cellstats, file=here(outdir, "qc", paste0(experiment_name,"_","cellstats.csv")))

for (i in 1:length(objs)){
  ImageDimPlot(objs[[i]], fov="fov", border.size=NA)
  ggsave(here(outdir,"qc","images", paste(names(objs)[i], "fullimg.png",sep=".")), dpi="retina", height=6, width=7)
}
#^Generate full image plots for each run

#################################################
#Module 2 Cell level QC