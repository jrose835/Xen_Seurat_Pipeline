# Xenium Seruat-based Transcript Processing Pipeline

Created for the Santangelo lab at Emory University


### How to Run:

```
#!/bin/bash 

cd /path/to/your/project 

R -e "renv::restore()" 

Rscript Xen_Seurat_LoadQCPCA.R 2>&1 | tee Xen_Seurat_LoadQCPCA.R.log

#Inspect outputs here to determine inputs for next script
#Then run:

Rscript Xen_Seruat_ClusterAnnot.R | tee Xen_Seruat_ClusterAnnot.R.log

```

