# Script 2 of Santangelo lab Xenium Seurat Transcript Processing Pipeline

# This script covers:
# 1. Finding nearest neighbors
# 2. Clustering
# 3. UMAP 
# 4. Annotaiton of clusters

# Description: 

# Usage:

# Inputs:
#     - input_obj_path: Path to processed SeuratObject from Script 1 (Xen_Seurat_LoadQCPCA.R)
#     - outdir: Main output directory path


# Parameters:
#     - ndim: Number of PCA components to use in nearest neighbor, clustering 
#     - reduction: Type of reduction to use for NN and UMAP. Options include: pca, integrated.rpca, integrated.harm, integrated.CCA
#     - clus_res: Clustering resolution


# Outputs:
#     - output/pipeline/objs Seurat object with UMAP and clustering at resolution specified as well as multiple other cluster resolutions
#     - output/pipeline/umap UMAP plots at each resolution cluster
#     - output/pipeline/clustree.png clustree graph plotted for each resolution range specified
#     - output/visualization UMAP plots with both cluster ID and run ID to color each cell
#     - output/visualization Image plots colored by 


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
source(here("src", "Xen_Seurat_functions.R"))

set.seed(1984)

#plan("multisession", workers=20)
options(future.globals.maxSize = 16000 * 1024^2)

#################################################
#Inputs & Parameters

indir <- "output/pipeline/objs" #--INPUT-- Note: If you used LoadQCPCA.R to produce you may not need to change
infile <- "Abbie_lung_int-None_scTrans_PCA.rds" #--INPUT--
input_obj_path <- here(indir, infile) 

experiment_name <- "Abbie_lung" #--INPUT--

outdir <- "output" #--INPUT--
dir.create.check(here(outdir))

ndim <- 35 #--PARAM--
reduction <- "pca" #--PARAM--
clus_res <- 0.4

#################################################
#Load object, NearestNeighbors, UMAP

obj <- readRDS(input_obj_path)

obj <- FindNeighbors(obj, dims=1:ndim, reduction=reduction)
obj <- RunUMAP(obj, dims=1:ndim, reduction=reduction)

#################################################
#Clustering

obj <- MultiResCluster(obj, res_range = c(0,1), 
                       outdir=here(outdir, "pipeline"), 
                       perform_clustering = TRUE,
                       clustree=TRUE,
                       plot_UMAP=TRUE
                       )
Idents(obj) <- obj@meta.data[[paste0("res.",clus_res)]]

#################################################
#General Visualization

vis_out <- here(outdir, "visualization")
dir.create.check(vis_out)

###
# UMAP plots

umap1 <- DimPlot(obj, reduction="umap", raster.dpi=c(1024,1024))
umap2 <- DimPlot(obj, reduction="umap", group.by=c("run_name"), raster.dpi=c(1024,1024))
umap3 <- DimPlot(obj, reduction="umap", split.by=c("run_name"), raster.dpi=c(1024,1024))

umaps <- umap1 + umap2
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "umaps.png")), plot=umaps, height=8, width=16)
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "umaps.pdf")), plot=umaps, height=8, width=16)

ggsave(filename=here(vis_out, paste0(experiment_name, "_", "runs_umaps.png")), plot=umap3, height=8, width=16)
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "runs_umaps.pdf")), plot=umap3, height=8, width=16)

### 
# Image plots

img_plots_cmb <- ImageDimPlot(obj, fov=names(obj@images),size=0.75, border.color="black", combine=TRUE)
img_plots_list <- ImageDimPlot(obj, fov=names(obj@images),size=0.75, border.color="black", combine=FALSE)
names(img_plots_list) <- names(obj@images)

ggsave(filename=here(vis_out, paste0(experiment_name, "_", "imgclusters.png")), plot=img_plots_cmb, height=8, width=16)
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "imgclusters.pdf")), plot=img_plots_cmb, height=8, width=16)

for (i in 1:length(img_plots_list)){
  ggsave(filename=here(vis_out, paste0(experiment_name, "_", names(img_plots_list)[i],"_imgclusters.png")), plot=img_plots_list[[i]], height=8, width=9)
  ggsave(filename=here(vis_out, paste0(experiment_name, "_", names(img_plots_list)[i],"_imgclusters.pdf")), plot=img_plots_list[[i]], height=8, width=9)
}

#################################################
#Automated annotation

## UNDER CONSTRUCTION ##

#################################################
#Save output

outfile <- gsub("PCA", paste0("clusRes.", clus_res), infile)
saveRDS(obj, file=here(outdir, "pipeline/objs", outfile))

#################################################
#Session Info

capture.output(sessionInfo(), file=here(outdir, "pipeline",paste0("Rsession.info.ClusterAnnot.",gsub("\\D", "", Sys.time()), ".txt")))
