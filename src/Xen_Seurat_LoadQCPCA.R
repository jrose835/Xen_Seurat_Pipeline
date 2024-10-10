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
#     - man_file: Manifest file (.csv) with 'run_name', 'experiment', 'condition', 'xen_dir', 'alt_input', and 'alt_input_dir' columns
#     - outdir: Main output directory path

# Parameters:
#     - min_nCount: Minimum number of counts per cell
#     - min_nFeature: Minimum number of unique genes per cell
#     - min_cellarea: Minimum cell area
#     - max_cellarea: Maximum cell area
#     - integraton_method: Parameter for IntegrateLayers (Seurat v5) specifying what method to be used for integration. 
#          Options include: 
#           - None: No integration method is performed, just simple concatenation
#           - Seuart_CCA: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6700744/
#           - Seurat_RPCA: Better for large datasets, and those with unequal representation of cells. https://satijalab.org/seurat/articles/integration_rpca.html
#           - harmony: https://github.com/immunogenomics/harmony

# Outputs:
#     - output/pipeline/objs Integrated & processed Seurat object .rds file
#     - output/qc/experiment_cellstats.csv
#     - output/qc/experiment_filtered_cellstats.csv
#     - output/qc/images/ full and filtered images .png
#     - output/qc/plots {cell area hist, nCount threshold, nFeat threshold, count by feature scatter plot} .png
#     - output/pipeline/objs Individual Seruat Objects, filtered for qc but no other processing
#     - output/pipeline/pca experiment_PCAelbow.png
#     - output/pipeline/pca/joint_unintegrated dim loading and pca plots for several PCs
#     - Rsession info .txt file

# Last updated: 10Oct2024
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

#plan("multisession", workers=20)
options(future.globals.maxSize = 16000 * 1024^2)

#################################################
# Inputs, Parameters, and Setup

man_file <- "sample_manifest.csv" #--INPUT--
outdir <- "output" #--INPUT--

min_nCount = 40 #--PARAM--
min_nFeature = 15 #--PARAM--
min_cellarea = 10 #--PARAM--
max_cellarea = 200 #--PARAM--

integraton_method = "Seurat_RPCA" #--PARAM-- Other options: "Seura_CCA", "harmony"

### Set up

manifest <- read_csv(here(man_file)) # Get run manifest
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

### QC plots for filters

qc_plots <- map(objs, CellQCPlots, min_nCount = min_nCount, min_nFeature = min_nFeature, min_cellarea = min_cellarea, max_cellarea = max_cellarea)

for (i in 1:length(qc_plots)){
  run_name = names(qc_plots)[i]
  imap(qc_plots[[i]], ~ ggsave(filename=here(outdir, "qc/plots", paste0(run_name, "_",.y, ".png")), plot=.x, height=6, width=6.5))
}

#saveRDS(qc_plots, file=here(outdir, "qc/plots", paste0(experiment_name,"_","qc_plots.RDS"))) # This is causing issues for some reason

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

imap(objs, ~saveRDS(object=.x,file=here(outdir,"pipeline/objs", paste0(experiment_name,"_", .y,"_", "filteredSeuratobj.rds")))) #Saving intermediate objects

#################################################
#Module 4.1 multi-sample integration (continued below)

#For quick restart of pipeline, used for testing only
# XN006A_chol1 <- readRDS(here(outdir, "pipeline/objs", "Abbie_lung_XN006A_chol1_filteredSeuratobj.rds"))
# XN006A_PBS <- readRDS(here(outdir, "pipeline/objs", "Abbie_lung_XN006A_PBS_filteredSeuratobj.rds"))
# objs <- list(XN006A_chol1,XN006A_PBS)
# names(objs) <- c("XN006A_chol1", "XN006A_PBS")

obj.full <- reduce(objs,merge)

if (integraton_method=="None"){
 obj.full<- JoinLayers(obj.full)
}

#################################################
#Module 5 normalization

obj.full <- SCTransform(obj.full, assay = "Xenium")

# I plan on expanding options here later...cell area based normalization?

#################################################
#Module 6 dimension reduction

obj.full <- RunPCA(obj.full)

PCA_outdir <- here(outdir, "pipeline/PCA", paste0("int-",integraton_method))
dir.create.check(PCA_outdir)
dir.create.check(here(PCA_outdir, "dim_loadings"))
dir.create.check(here(PCA_outdir, "pca_plots"))

elbow <- ElbowPlot(obj.full, ndims=50)
ggsave(here(PCA_outdir, paste0(experiment_name, "_","PCAelbow.png")), plot = elbow)

AllVizDimLoadings(obj.full, outdir=PCA_outdir)

PCAplots(obj.full, ndim=10, outdir=PCA_outdir)

#################################################
#Module 4.2 multi-sample integration

# Note: Will need to change if normalization method changes

if (integration_method=="Seurat_RPCA"){
  obj.full <- IntegrateLayers(object = obj.full, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method="SCT",
                  verbose = TRUE)
} elseif (integraton_method=="harmony"){
  obj.full <- IntegrateLayers(object = obj.full, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harm", normalization.method="SCT",
                              verbose = TRUE)
} elseif (integraton_method=="Seuart_CCA"){
  obj.full <- IntegrateLayers(object = obj.full, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.CCA", normalization.method="SCT",
                              verbose = TRUE)
}

#################################################
#Save output

saveRDS(obj.full, file=here(outdir, "pipeline/objs", paste0(experiment_name,"_","int-",integraton_method, "_", "scTrans_PCA.rds")))

#################################################
#Session Info

capture.output(sessionInfo(), file=here(outdir, "pipeline",paste0("Rsession.info.LoadQCPCA.",gsub("\\D", "", Sys.time()), ".txt")))
