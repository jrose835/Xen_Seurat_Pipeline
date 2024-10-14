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

indir <- "output/pipeline/objs" #--INPUT--
infile <- "Abbie_lung_int-None_scTrans_PCA.rds" #--INPUT--
input_obj_path <- here(indir, "Abbie_lung_int-None_scTrans_PCA.rds") 

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
                       clustree=FALSE
                       )
Idents(obj) <- obj@meta.data[[paste0("res.",clus_res)]]

#################################################
#General Visualization

vis_out <- here(outdir, "visualization")
dir.create.check(vis_out)

umap1 <- DimPlot(obj, reduction="umap", raster.dpi=c(1024,1024))
umap2 <- DimPlot(obj, reduction="umap", group.by=c("run_name"), raster.dpi=c(1024,1024))

umaps <- umap1 + umap2
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "umaps.png")), plot=umaps, height=8, width=16)
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "umaps.pdf")), plot=umaps, height=8, width=16)

fov1 <- ImageDimPlot(obj, fov="fov", group.by="seurat_clusters", size=0.75, border.color = "black")
fov2 <- ImageDimPlot(obj, fov="fov.2", group.by="seurat_clusters",size=0.75, border.color = "black")
#Need to change this to not be hardcoded and use as many fov as exist in object

imgplot <- fov1 + fov2
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "imgclusters.png")), plot=imgplot, height=8, width=16)
ggsave(filename=here(vis_out, paste0(experiment_name, "_", "imgclusters.pdf")), plot=imgplot, height=8, width=16)


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
