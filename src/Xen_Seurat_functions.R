# Functions for Xenium Seurat Pipleine

# Last updated: 13Dec2024
# Author: jrose

#################################################

library(Seurat)
library(tidyverse)
library(future)
library(clustree)
#library(spacexr)
library(nanoparquet)
#library(scBubbletree)
library(uwot)
library(Matrix)
library(here)

#################################################
# LoadMultiXenium
# module 1 function. Loads multiple xenium experiments according to manifest and assigns metaadata
# input: manifest dataframe
# output: list of SeuratObjects

LoadMultiXenium <- function(manifest, assay="Xenium", segmentations = NULL){
  
  outs <- list()
  
  dirs <- manifest$xen_dir
  
  for (i in 1:length(dirs)){
    
    run <- manifest$run_name[i]
    
    print(paste0("Loading", " ", run, "--", Sys.time()))
    
    obj <- LoadXenium(dirs[i], fov=run, assay=assay, segmentations = segmentations)
    
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

#################################################
#LogAreaNormalize
#Module 5. Very much inspired by Seurat v5 LogNormalize.deafult
# input: sparse matrix of counts data, meta data data frame containing "cell_area" values
# output: prefix string for clustering variables
#@PARAM: scale.factor : transcripts per X um
#@PARAM: margin : 2 for column, 1 for row

#Expression data operator
LogAreaNormalize.default <- function(
    data,
    metadata,
    scale.factor=100,
    margin=2L,
    ...
){
  # Check margin validity
  margin <- .CheckFmargin(fmargin = margin)
  
  # Save original row/column names
  rownames_original <- rownames(data)
  colnames_original <- colnames(data)
  
  # Ensure 'cell_area' exists in metadata
  if (!"cell_area" %in% colnames(metadata)) {
    stop("Metadata must contain a 'cell_area' column for normalization.")
  }
  
  # Retrieve and validate cell_area
  cell_area <- metadata$cell_area
  if (any(is.na(cell_area) | cell_area <= 0)) {
    stop("All 'cell_area' values must be positive numbers.")
  }
  
  # Normalize the data matrix
  if (margin == 2L) {
    # Normalize over columns
    cell_area_matrix <- Matrix::Diagonal(x = 1 / cell_area) # Efficient scaling
    data <- data %*% cell_area_matrix
  } else {
    # Normalize over rows
    cell_area_matrix <- Matrix::Diagonal(x = 1 / cell_area)
    data <- cell_area_matrix %*% data
  }
  
  # Apply scaling factor and log1p in bulk
  data <- log1p(data * scale.factor)
  
  # Restore original row/column names
  rownames(data) <- rownames_original
  colnames(data) <- colnames_original
  
  return(data)
}  

#################################################
#LogAreaNormalize
#A Seruat Object "Method" level function calling LogAreaNormalize.default on Seurat objects
#input: Seurat object with counts data in assay counts slot selected
#output: Seurat object with area normalized counts in data slot

LogAreaNormalize <- function(object, 
                             assay="Xenium",
                             scale.factor=100,
                             margin=2L,
                             ...
) {
  print(paste0("Area Normalizing Data", "--", Sys.time()))
  object[[assay]]$data <- LogAreaNormalize.default(data=object[[assay]]$counts,
                                                   metadata = object@meta.data,
  )
  object <- LogSeuratCommand(object=object)
  return(object)
}


#################################################
#UMAP functions
#input: Matrix of embedings (use Seurat::Embeddings())
#output: Umap embedding

# For calculating UMAP off of area-normalized counts
# This is literally just the function copied from Seurat v5 github
# I've run into a bug when trying to apply RunUMAP() on objects with custom normalization. So here I am avoiding it by running it myself
# I've diagnosed the problem has having somethign to do with RunUMAP.Seurat method from Seruat v5 / SeruatObjects but can't say for sure what's causing it.

RunUMAP.default <- function(
    object,
    reduction.key = 'UMAP_',
    assay = NULL,
    reduction.model = NULL,
    return.model = FALSE,
    umap.method = 'uwot',
    n.neighbors = 30L,
    n.components = 2L,
    metric = 'cosine',
    n.epochs = NULL,
    learning.rate = 1.0,
    min.dist = 0.3,
    spread = 1.0,
    set.op.mix.ratio = 1.0,
    local.connectivity = 1L,
    repulsion.strength = 1,
    negative.sample.rate = 5,
    a = NULL,
    b = NULL,
    uwot.sgd = FALSE,
    seed.use = 42,
    metric.kwds = NULL,
    angular.rp.forest = FALSE,
    densmap = FALSE,
    dens.lambda = 2,
    dens.frac = 0.3,
    dens.var.shift = 0.1,
    verbose = TRUE,
    ...
) {
  CheckDots(...)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (umap.method != 'umap-learn' && getOption('Seurat.warn.umap.uwot', TRUE)) {
    warning(
      "The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric",
      "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'",
      "\nThis message will be shown once per session",
      call. = FALSE,
      immediate. = TRUE
    )
    options(Seurat.warn.umap.uwot = FALSE)
  }
  if (umap.method == 'uwot-learn') {
    warning("'uwot-learn' is deprecated. Set umap.method = 'uwot' and return.model = TRUE")
    umap.method <- "uwot"
    return.model <- TRUE
  }
  if (densmap && umap.method != 'umap-learn'){
    warning("densmap is only supported by umap-learn method. Method is changed to 'umap-learn'")
    umap.method <- 'umap-learn'
  }
  if (return.model) {
    if (verbose) {
      message("UMAP will return its model")
    }
    umap.method = "uwot"
  }
  if (inherits(x = object, what = "Neighbor")) {
    object <- list( idx = Indices(object),
                    dist = Distances(object) )
  }
  if (!is.null(x = reduction.model)) {
    if (verbose) {
      message("Running UMAP projection")
    }
    umap.method <- "uwot-predict"
  }
  umap.output <- switch(
    EXPR = umap.method,
    'umap-learn' = {
      if (!py_module_available(module = 'umap')) {
        stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
      }
      if (!py_module_available(module = 'sklearn')) {
        stop("Cannot find sklearn, please install through pip (e.g. pip install scikit-learn).")
      }
      if (!is.null(x = seed.use)) {
        py_set_seed(seed = seed.use)
      }
      if (typeof(x = n.epochs) == "double") {
        n.epochs <- as.integer(x = n.epochs)
      }
      umap_import <- import(module = "umap", delay_load = TRUE)
      sklearn <- import("sklearn", delay_load = TRUE)
      if (densmap &&
          numeric_version(x = umap_import$pkg_resources$get_distribution("umap-learn")$version) <
          numeric_version(x = "0.5.0")) {
        stop("densmap is only supported by versions >= 0.5.0 of umap-learn. Upgrade umap-learn (e.g. pip install --upgrade umap-learn).")
      }
      random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
      umap.args <- list(
        n_neighbors = as.integer(x = n.neighbors),
        n_components = as.integer(x = n.components),
        metric = metric,
        n_epochs = n.epochs,
        learning_rate = learning.rate,
        min_dist = min.dist,
        spread = spread,
        set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity,
        repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        random_state = random.state,
        a = a,
        b = b,
        metric_kwds = metric.kwds,
        angular_rp_forest = angular.rp.forest,
        verbose = verbose
      )
      if (numeric_version(x = umap_import$pkg_resources$get_distribution("umap-learn")$version) >=
          numeric_version(x = "0.5.0")) {
        umap.args <- c(umap.args, list(
          densmap = densmap,
          dens_lambda = dens.lambda,
          dens_frac = dens.frac,
          dens_var_shift = dens.var.shift,
          output_dens = FALSE
        ))
      }
      umap <- do.call(what = umap_import$UMAP, args = umap.args)
      umap$fit_transform(as.matrix(x = object))
    },
    'uwot' = {
      if (is.list(x = object)) {
        umap(
          X = NULL,
          nn_method = object,
          n_threads = nbrOfWorkers(),
          n_components = as.integer(x = n.components),
          metric = metric,
          n_epochs = n.epochs,
          learning_rate = learning.rate,
          min_dist = min.dist,
          spread = spread,
          set_op_mix_ratio = set.op.mix.ratio,
          local_connectivity = local.connectivity,
          repulsion_strength = repulsion.strength,
          negative_sample_rate = negative.sample.rate,
          a = a,
          b = b,
          fast_sgd = uwot.sgd,
          verbose = verbose,
          ret_model = return.model
        )
      } else {
        umap(
          X = object,
          n_threads = nbrOfWorkers(),
          n_neighbors = as.integer(x = n.neighbors),
          n_components = as.integer(x = n.components),
          metric = metric,
          n_epochs = n.epochs,
          learning_rate = learning.rate,
          min_dist = min.dist,
          spread = spread,
          set_op_mix_ratio = set.op.mix.ratio,
          local_connectivity = local.connectivity,
          repulsion_strength = repulsion.strength,
          negative_sample_rate = negative.sample.rate,
          a = a,
          b = b,
          fast_sgd = uwot.sgd,
          verbose = verbose,
          ret_model = return.model
        )
      }
    },
    'uwot-predict' = {
      if (metric == 'correlation') {
        warning(
          "UWOT does not implement the correlation metric, using cosine instead",
          call. = FALSE,
          immediate. = TRUE
        )
        metric <- 'cosine'
      }
      if (is.null(x = reduction.model) || !inherits(x = reduction.model, what = 'DimReduc')) {
        stop(
          "If running projection UMAP, please pass a DimReduc object with the model stored to reduction.model.",
          call. = FALSE
        )
      }
      model <- Misc(
        object = reduction.model,
        slot = "model"
      )
      # add num_precomputed_nns to <v0.1.13 uwot models to prevent errors with newer versions of uwot
      if (!"num_precomputed_nns" %in% names(model)) {
        model$num_precomputed_nns <- 1
      }
      if (length(x = model) == 0) {
        stop(
          "The provided reduction.model does not have a model stored. Please try running umot-learn on the object first",
          call. = FALSE
        )
      }
      if (!"num_precomputed_nns" %in% names(x = model)) {
        model$num_precomputed_nns <- 0
      }
      if (is.list(x = object)) {
        if (ncol(object$idx) != model$n_neighbors) {
          warning("Number of neighbors between query and reference ",
                  "is not equal to the number of neighbors within reference")
          model$n_neighbors <- ncol(object$idx)
        }
        umap_transform(
          X = NULL,
          nn_method = object,
          model = model,
          n_threads = nbrOfWorkers(),
          n_epochs = n.epochs,
          verbose = verbose
        )
      } else {
        umap_transform(
          X = object,
          model = model,
          n_threads = nbrOfWorkers(),
          n_epochs = n.epochs,
          verbose = verbose
        )
      }
    },
    stop("Unknown umap method: ", umap.method, call. = FALSE)
  )
  if (return.model) {
    umap.output$nn_index <- NULL
    umap.model <- umap.output
    umap.output <- umap.output$embedding
  }
  colnames(x = umap.output) <- paste0(reduction.key, 1:ncol(x = umap.output))
  if (inherits(x = object, what = 'dist')) {
    rownames(x = umap.output) <- attr(x = object, "Labels")
  } else if (is.list(x = object)) {
    rownames(x = umap.output) <- rownames(x = object$idx)
  } else {
    rownames(x = umap.output) <- rownames(x = object)
  }
  umap.reduction <- CreateDimReducObject(
    embeddings = umap.output,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  if (return.model) {
    Misc(umap.reduction, slot = "model") <- umap.model
  }
  return(umap.reduction)
}


### NOTE:

# You can run the above like this to perform UMAP from area-normalized counts

# dims=1:30
# reduction="pca"
# data.use <- Embeddings(object[[reduction]])[, dims]
# 
# umap <- RunUMAP.default(object=data.use,
#                            reduction.key="umap_",
#                            assay="Xenium"
# )
#
# object[["umap"]] <- umap


#################################################
