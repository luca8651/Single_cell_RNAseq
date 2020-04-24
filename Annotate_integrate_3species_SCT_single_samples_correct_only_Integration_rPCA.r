library(Seurat)
library(future)
options(future.globals.maxSize=40000 * 1024^2)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(monocle3)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(scran)
library(batchelor)
library(scater)

load("Annotate_integrate_3species_SCT_single_samples_correct.RData") 
rm(seu_3species)  
rm(seu_hsa)  
rm(seu_mmu)  
rm(seu_nmr)  
rm(seu3)
rm(nmr.data)

###DOES NOT WORK WITHOUT SCALE DATA

for (i in 1:length(seu_list)) {
DefaultAssay(seu_list[[1]])="RNA"
}


seu_list <- lapply(X = seu_list, FUN = function(x) {
    x <- NormalizeData(x, verbose = TRUE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seu_list)

#    x <- ScaleData(x, features = features, verbose = TRUE)
#    x <- RunPCA(x, features = features, verbose = FALSE)

seu_list <- lapply(X = seu_list, FUN = function(x) {
    x <- ScaleData(x,features = features,  verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

#for (i in 1:length(seu_list)) {
#DefaultAssay(seu_list[[1]])="RNA"
#seu_list[[i]]<- NormalizeData(seu_list[[i]], verbose = FALSE)
#seu_list[[i]]<- FindVariableFeatures(seu_list[[i]], verbose = FALSE)
#seu_list[[i]] <- ScaleData(seu_list[[i]], features = features, verbose = FALSE)
#seu_list[[i]] <- RunPCA(seu_list[[i]], features = features, verbose = FALSE)
#}

anchors <- FindIntegrationAnchors(object.list = seu_list,reference=c(1,6,15,22,31), reduction = "rpca", dims = 1:30)

saveRDS(anchors,file = "/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Annotate_integrate_3species_SCT_single_samples_correct_RpcaAnchors.rds")
#seu_3species_new <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)
seu_3species_new <- IntegrateData(anchorset = anchors,  verbose = TRUE)
saveRDS(seu_3species_new,file = "/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Annotate_integrate_3species_SCT_single_samples_correct_seuObject_bsub_rPCA.rds" )


