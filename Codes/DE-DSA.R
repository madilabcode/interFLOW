rm(list = ls())
#library(dplyr)
library(Seurat)
#library(knitr)
library(ggplot2)
#library(plotly)
#library(stringr)
#library(tidyverse)  
#library(reshape2)
#library(pheatmap)
library(dendsort)
library(RColorBrewer)
#library(grid)
library(igraph)
library(MAST)
library(circlize)
#library(viridis)
library(ComplexHeatmap)
#library(gridBase)
#library(formattable)
#library(kableExtra)
library(reticulate)
library(hash)
source('Codes/signture_utils.R', echo=TRUE)
source('Codes/Ligand_Receptor_pipeline.R', echo=TRUE)
source_python("Codes/utils.py")
source_python("Codes/neo_connection_db.py")

conf = read.csv("config.csv")
pathTret = conf[conf$Var == "objTr","Value"]
pbmc_Tr = readRDS(paste("InputObjs/",as.character(pathTret),sep =""))
name_obj_Tr = as.character(conf[conf$Var == "nameObjTr","Value"])
pathControl = conf[conf$Var == "objCo","Value"]
pbmc_Co = readRDS(paste("InputObjs/",as.character(pathControl),sep =""))
name_obj_Co = as.character(conf[conf$Var == "nameObjCo","Value"])
projName = conf[conf$Var == "Proj","Value"]

graph = start_connection()
lst =  return_all_DSA_and_legRet_of_proj(graph,projName,c(name_obj_Tr,name_obj_Co))

DSA_LST = lst[[1]]
legRet_LST = lst[[2]]

result = data.frame()


fromClusters = intersect(names(DSA_LST[[name_obj_Tr]]), names(DSA_LST[[name_obj_Co]]))

for (fromCluster in fromClusters){
  print(fromCluster)
  toClusters = intersect(names(DSA_LST[[name_obj_Tr]][[fromCluster]]), names(DSA_LST[[name_obj_Co]][[fromCluster]]))
  for (toCluster in toClusters){
    print(toCluster)
    DSA_Table_Tr = DSA_LST[[name_obj_Tr]][[fromCluster]][[toCluster]]
    DSA_Table_Co = DSA_LST[[name_obj_Co]][[fromCluster]][[toCluster]]
    
    legRet_Tr = legRet_LST[[name_obj_Tr]][[fromCluster]][[toCluster]]
    legRet_Co = legRet_LST[[name_obj_Co]][[fromCluster]][[toCluster]]
    
    toTr = subset(pbmc_Tr,ident = toCluster)
    toCo = subset(pbmc_Co,ident = toCluster)
    
    fromTr = subset(pbmc_Tr,ident = fromCluster)
    fromCo = subset(pbmc_Co,ident = fromCluster)
    
    toExpTr = GetAssayData(object = toTr, slot = "counts") %>% as.data.frame()
    toExpCo = GetAssayData(object = toCo, slot = "counts") %>% as.data.frame()
    
    fromExpTr = GetAssayData(object = fromTr, slot = "counts") %>% as.data.frame()
    fromExpCo = GetAssayData(object = fromCo, slot = "counts") %>% as.data.frame()
    
    
    DSA_Tr = dsa_with_lig(toExpTr, fromExpTr, DSA_Table_Tr, legRet_Tr)
    DSA_Co = dsa_with_lig(toExpCo, fromExpCo, DSA_Table_Co, legRet_Co)
    tlst = DSA_DE(DSA_Tr,DSA_Co)
    
    temp = c(as.numeric(fromCluster),as.numeric(toCluster),as.numeric(tlst[[1]]),as.numeric(tlst[[2]]))
    result = rbind(result, temp)
  }
  
}

colnames(result) = c("From","To" ,"Stat"  ,"Pvalue" )
write.csv(result,"outputObj/DE-DSA.csv")

upDSAMean = return_all_DSA_Mean_of_proj(graph, projName, de_plot = result)
downDSAMean = return_all_DSA_Mean_of_proj(graph, projName, de_plot = result, Up = FALSE)