# Functions for Xenium Seurat Pipleine

# Last updated: 07Oct2024
# Author: jrose

#################################################

library(Seurat)
library(tidyverse)
library(future)
library(clustree)
#library(spacexr)
library(nanoparquet)
#library(scBubbletree)
library(here)

#################################################
# LoadMultiXenium
# module 1 function. Loads multiple xenium experiments according to manifest and assigns metaadata
# input: manifest dataframe
# output: list of SeuratObjects

LoadMultiXenium <- function(manifest, assay="Xenium"){
  
  outs <- list()
  
  dirs <- manifest$xen_dir
  
  for (i in 1:length(dirs)){
    
    run <- manifest$run_name[i]
    
    print(paste0("Loading", " ", run, "--", Sys.time()))
    
    obj <- LoadXenium(dirs[i], fov=run, assay=assay)
    
    obj$run_name <- run
    obj$experiment <- manifest$experiment[i]
    obj$condition <- manifest$condition[i]
    
    outs[[run]] <- obj
    
    print(paste0("Finished", " ", run, "--", Sys.time()))
  }

  return(outs)
}

#################################################
# XenCellStats
# module 1 function. Generates descriptive statistics on transcript data from Xenium experiments
# input: SeuratObjects (with spatial fov)
# output: dataframe of cell numbers, counts per cell, 

XenSeuratStats <- function(obj){
  
  #dir.create.check(outdir)
  #Check to make sure outptu dir exisits
  
  out <- data.frame(run= obj$run_name[1], 
                    n.cell=ncol(obj), 
                    n.feature=nrow(obj), 
                    counts.per.cell.min=summary(obj$nCount_Xenium)[[1]],
                    counts.per.cell.Q1 = summary(obj$nCount_Xenium)[[2]],
                    counts.per.cell.median = summary(obj$nCount_Xenium)[[3]],
                    counts.per.cell.mean = summary(obj$nCount_Xenium)[[4]],
                    counts.per.cell.Q3 = summary(obj$nCount_Xenium)[[5]],
                    counts.per.cell.max = summary(obj$nCount_Xenium)[[6]],
                    feats.per.cell.min = summary(obj$nFeature_Xenium)[[1]],
                    feats.per.cell.Q1 = summary(obj$nFeature_Xenium)[[2]],
                    feats.per.cell.median = summary(obj$nFeature_Xenium)[[3]],
                    feats.per.cell.mean = summary(obj$nFeature_Xenium)[[4]],
                    feats.per.cell.Q3 = summary(obj$nFeature_Xenium)[[5]],
                    feats.per.cell.max = summary(obj$nFeature_Xenium)[[6]]
                    )
  rownames(out) <- NULL
  
  return(out)
  
}


#################################################
# PrepOutDir
# module 1 function. Generates output directory for pipeline
# input: manifest dataframe, output directory path
# output: A series of directories customized for the pipeline 

prepOutDir <- function(outdir="output", manifest){
  base_dir <- here(outdir)
  
  #Set up pipeline and qc directories
  dir.create.check(here(base_dir))
  dir.create.check(here(base_dir, "pipeline"))
  dir.create.check(here(base_dir, "pipeline", "PCA"))
  dir.create.check(here(base_dir, "pipeline", "objs"))
  dir.create.check(here(base_dir, "qc"))
  dir.create.check(here(base_dir, "qc", "images"))
  dir.create.check(here(base_dir, "qc", "plots"))
  
  
  # Create PCA directories for each run
  # runs <- manifest$run_name
  # 
  # for (run in runs) {
  #   run_dir <- here(base_dir, "pipeline", "PCA", run)
  #   dim_loadings_dir <- here(run_dir, "dim_loadings")
  #   pca_plots_dir <- here(run_dir, "pca_plots")
  #   
  #   # Create the directories for the current run
  #   dir.create.check(dim_loadings_dir, recursive=TRUE)
  #   dir.create.check(pca_plots_dir, recursive=TRUE)
  # }
  
}

#I'm not sure why I have to make this instead of it being a parameter in dir.create()...
dir.create.check <- function(path, ...) {
   if (dir.exists(path)) {
     print(paste("Directory already exists:", path))
   } else {
     dir.create(path, ...)
   }
}
#################################################
# AppendAreaParquet
# module 2 function. Adds cell and nucleus area data for each cell from parquet files
# input: SeuratObject, path to xenium output directory
# output: SeuratObject with new metadata appended 

AppendAreaParquet <- function(obj, xen_dir){
  
  parquet_file <- here(xen_dir, "cells.parquet")
  add_dat <- read_parquet(parquet_file) %>% select(cell_id, ends_with("_area"))
  
  add_dat <- add_dat[match(add_dat$cell_id, rownames(obj[[]])),] #making sure order of cells is the same
  
  obj$cell_area <- add_dat$cell_area
  obj$nucleus_area <- add_dat$nucleus_area
  
  return(obj)
  
}

#################################################
# AppendIFdata
# module 2 function. Add MFI data from immunofluroesence pipeline .txt file
# input: SeuratObject, path to IF data
# output: SeuratObject with MFI data

AppendIFdata <- function(obj, IF_data_path){
  
  add_dat <- read_tsv(IF_data_path)
  
  add_dat <- add_dat[match(add_dat$Name, rownames(obj[[]])),] #making sure order of cells is the same
  
  obj$MFI <- add_dat$MFI
  
  return(obj)
}

#################################################
# CellQCPlots
# module 3 function. Plot qc plots related to nCounts and nFeatures
# input: SeuratObject, thresholds for nCount and nFeature
# output: List of qc plots

CellQCPlots <- function(obj, min_nCount, min_nFeature, min_cellarea, max_cellarea){
  
  meta <- obj[[]] #extract meta data
  run_name <- unique(meta$run_name)
  
  plots <- list()
  
  plots[["nCount_threshold_bar"]] <- ThresholdBar(obj)
  
  plots[["nFeat_threshold_bar"]] <- ThresholdBar(obj,  assay="nFeature_Xenium", from=5, to=40, by=5)
  
  plots[["count_feat_scatter"]] <- ggplot(meta, aes(x=nCount_Xenium, y=nFeature_Xenium)) + 
                                        geom_point() + 
                                        #geom_density_2d() + 
                                        geom_hline(yintercept = min_nFeature, linetype="dashed",color="red") +
                                        geom_vline(xintercept = min_nCount, linetype="dashed", color="red") +
                                        labs(x="Counts per cell", y="Genes per cell") +
                                        theme_light()
  
  plots[["area_hist"]] <- ggplot(meta, aes(x=cell_area)) + 
                              geom_histogram(binwidth = 10) + 
                              geom_vline(xintercept=min_cellarea,linetype="dashed",color="red") + 
                              geom_vline(xintercept=max_cellarea,linetype="dashed",color="red") + 
                              theme_light()
  return(plots)
  # UNDER CONSTRUCTION
  
}

#################################################
# ThresholdBar
# module 3 subfunction. Plot effects of changing thresholds
# input: 
# output: 

ThresholdBar <- function(obj, from=10, to=200, by=20, assay="nCount_Xenium"){
  thresholds <- seq(from=from, to=to, by=by)
  
  filter_outputs <- data.frame(threshold=thresholds, kept=rep(NA, length(thresholds)), removed=rep(NA, length(thresholds)))
  
  for (i in 1:length(thresholds)){
    tmp <- table(obj[["nCount_Xenium"]]< thresholds[i])
    filter_outputs$kept[i] <- tmp[["FALSE"]]
    filter_outputs$removed[i] <- tmp[["TRUE"]]
  }
  
  plot <- pivot_longer(filter_outputs, cols=c('kept', "removed")) %>%
    #subset(name=="removed") %>%
    ggplot(aes(x=threshold, y=value)) + 
    labs(y="Number of Cells", x="threshold", title=assay) +
    #geom_line() +
    geom_col(aes(fill=name)) + 
    theme_light()
  
  return(plot)
}

#################################################
# AllVizDimLoadings
# module 4 subfunction. Plot a series of dim loading plots for each PC of PCA
# input: Seurat object, number of PCs (dim), output directory
# output: Dim_loaing_n.png

AllVizDimLoadings <- function(obj, dim=50, outdir,...) {
  
  dims <- 1:dim
  num_chunks <- 5
  dim_chunks <- split(dims, gl(num_chunks, ceiling(length(dims) / num_chunks), length(dims)))
  
  for (i in 1:length(dim_chunks)){
    plot <- VizDimLoadings(obj, dims=dim_chunks[[i]], ncol=5)
    ggsave(filename=here(outdir, "dim_loadings",paste0("Dim_loadings_", i, ".png")),plot=plot, height=15, width=15)
  }
    
}

#################################################
# PCAplots
# module 4 subfunction. Plot a series of PCA plots
# input: Seurat object, number of PCs, output directory
# output: 

PCAplots <- function(obj, ndim=10, outdir) {
  
  pairs = combn(1:ndim, 2)
  
  for (i in 1:ncol(pairs)){
    plot <- DimPlot(obj, reduction="pca", dims=pairs[,i]) + NoLegend()
    ggsave(filename=here(outdir, "pca_plots",paste0("pca_plot_", pairs[1,i], ".v.",pairs[2,i], ".png")),plot=plot)
  }
  
}


#################################################
# MultiResCluster
# Script 2 module 1 subfunction. Performs clustering 
# input: Seurat object, resolution range, output directory
# output: Seurat object with cluster ids for a range of resolutions, clustree graph, umaps showing clusters at each resolution

MultiResCluster <- function(obj, res_range=c(0,1), outdir, perform_clustering=TRUE, clustree=TRUE, plot_UMAP=TRUE) {
  res <- seq(from=res_range[1], to=res_range[2], by=0.1)
  
  if (perform_clustering==TRUE){
    
    obj <- FindClusters(obj, resolution=res, verbose=TRUE)
    
    prefix <- GetClusterPrefix(obj) #See helper function below
    
    if (plot_UMAP==TRUE){
      for (i in 1:length(res)){
        outdir_umap <- here(outdir, "umap")
        dir.create.check(here(outdir, "umap"))
        umap <- DimPlot(obj, group.by=paste0(prefix,res[i]))
        ggsave(filename=here(outdir_umap, paste0("umap.res.",res[i],".png")), plot=umap,height=5, width=6.5)
      }
    }
  }
  
  prefix <- GetClusterPrefix(obj)
  
  if (clustree==TRUE){
    clustree_plot <- clustree(obj, prefix = prefix, node_colour = "sc3_stability", node_label=prefix, node_text_colour="grey90")
    ggsave(filename=here(outdir, "clustree.png"), plot=clustree_plot, width=14, height=12)
  }
  
  return(obj)
}
#################################################
# GetClusterPrefix
# Script 2 helper function for getting/setting clustering metadata variable name
# input: Seurat object
# output: prefix string for clustering variables

GetClusterPrefix <- function(obj){
  
  res_names <- names(obj@meta.data)[grep("res.", names(obj@meta.data))]
  prefix <- sub("(.*_res\\.).*", "\\1", res_names) %>% unique()
  
  if (length(prefix)>1){
    warning("Multiple clustering prefixes detected")
  }
  
  return(prefix)
}