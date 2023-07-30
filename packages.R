if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

lib_list = c("Seurat","knitr","ggplot2","plotly","stringr","tidyverse","reshape2","pheatmap",
             "pheatmap","dendsort","RColorBrewer","grid","igraph","MAST","circlize","viridis","ComplexHeatmap",
             "gridBase","formattable","kableExtra","reticulate","hash")


for (pack in lib_list){
  if (!require(pack, quietly = TRUE))
    BiocManager::install(pack)
}