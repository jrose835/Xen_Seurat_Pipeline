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

# Last updated: 08Oct2024
# Author: jrose
#################################################
#Set up

library(Seurat)
library(tidyverse)
library(future)
library(clustree)
library(spacexr)
library(scBubbletree)
library(nanoparquet)
library(here)
source(here("src", "Xen_Seurat_functions.R"))  #accompanying custom functions

set.seed(1984)


man_file <- "sample_manifest.csv"
manifest <- read_csv(here(man_file)) # Get run manifest

outdir <- "output"
prepOutDir(outdir=outdir, manifest=manifest) # Prepare output dir structure. Custom func

n_exprs <- length(unique(manifest$experiment)) #Checking for number of distinct experiments. Only support 1 per manifest as of now

if (n_exprs==1){
  experiment_name <- unique(manifest$experiment)
} else {
  stop("Error: Please create a single manifest file per experiment and run script individually")
}

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
#Module 2 append additional cell info

objs <- map2(objs,  manifest$xen_dir, AppendAreaParquet) # Adds cell and nucleus area

for (i in 1:length(objs)){          
  if (manifest$alt_input[i]){
    objs[[i]] <- AppendIFdata(objs[[i]], manifest$alt_input_dir[i]) 
  } else {
    print(paste("No additional IF data specified for",manifest$run_name[i]))
  }
}
#^Only adds IF data if alt_input value in manifest == True

#################################################
#Module 3 cell level QC

min_nCount = 40
min_nFeature = 15
min_cellarea = 10
max_cellarea = 200

### QC plots for filters

qc_plots <- map(objs, CellQCPlots, min_nCount = min_nCount, min_nFeature = min_nFeature, min_cellarea = min_cellarea, max_cellarea = max_cellarea)

for (i in 1:length(qc_plots)){
  run_name = names(qc_plots)[i]
  imap(qc_plots[[i]], ~ ggsave(filename=here(outdir, "qc/plots", paste0(run_name, "_",.y, ".png")), plot=.x, height=6, width=6.5))
}

saveRDS(qc_plots, file=here(outdir, "qc/plots", paste0(experiment_name,"_","qc_plots.RDS")))

### Applying filters

objs <- map(objs, subset, subset=nCount_Xenium > min_nCount)
objs <- map(objs, subset, subset=nFeature_Xenium > min_nFeature)
objs <- map(objs, subset, subset=((cell_area > min_cellarea)&(cell_area<max_cellarea)))

### Creating stats and image plots of filtered cells

filtered_cellstats <- do.call(rbind, lapply(objs, XenSeuratStats)) #Summary statistics for each run. Custom XenSeuratStats func
write.csv(filtered_cellstats, file=here(outdir, "qc", paste0(experiment_name,"_","filtered_cellstats.csv")))

for (i in 1:length(objs)){
  ImageDimPlot(objs[[i]], fov="fov", border.size=NA)
  ggsave(here(outdir,"qc","images", paste(names(objs)[i], "filteredimg.png",sep=".")), dpi="retina", height=6, width=7)
}


#################################################
#Module 4 multi-sample integration




#################################################
#Module 5 normalization




#################################################
#Module 6 dimension Reduction



