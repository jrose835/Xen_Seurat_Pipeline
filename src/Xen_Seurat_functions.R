# Functions for Xenium Seurat Pipleine

# Last updated: 07Oct2024
# Author: jrose

#################################################

library(Seurat)
library(tidyverse)
library(future)
library(clustree)
library(spacexr)
library(scBubbletree)
library(here)

#################################################
# LoadMultiXenium
# module 1 function. Loads multiple xenium experiments according to manifest and assigns metaadata
# input: manifest dataframe
# output: list of SeuratObjects

LoadMultiXenium <- function(manifest, fov="fov", assay="Xenium"){
  
  outs <- list()
  
  dirs <- manifest$xen_dir
  
  for (i in 1:length(dirs)){
    
    run <- manifest$run_name[i]
    
    print(paste0("Loading", " ", run, "--", Sys.time()))
    
    obj <- LoadXenium(dirs[i], fov=fov, assay=assay)
    
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
  dir.create.check(here(base_dir, "qc"))
  dir.create.check(here(base_dir, "qc", "images"))
  dir.create.check(here(base_dir, "qc", "plots"))
  
  
  # Create PCA directories for each run
  runs <- manifest$run_name
  
  for (run in runs) {
    run_dir <- here(base_dir, "pipeline", "PCA", run)
    dim_loadings_dir <- here(run_dir, "dim_loadings")
    pca_plots_dir <- here(run_dir, "pca_plots")
    
    # Create the directories for the current run
    dir.create.check(dim_loadings_dir, recursive=TRUE)
    dir.create.check(pca_plots_dir, recursive=TRUE)
  }
  
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
