# Functions for Xenium Seurat Pipleine

# Last updated: 04Oct2024
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
# input: manifest data frame
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

XenCellStats <- function(obj){
  
}