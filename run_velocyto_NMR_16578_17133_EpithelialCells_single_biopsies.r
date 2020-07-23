library(Seurat)
library(SeuratWrappers)
library(monocle)
library(monocle3)
library(future)
options(future.globals.maxSize=8500 * 1024^2)
library(dplyr)
library(gridExtra)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(scran)
library(batchelor)
library(scater)
library(velocyto.R)
library(gridExtra)
library(gridBase)
library(cowplot)
library(grDevices)
library(png)
library(grid)
system("module load gdal/3.0.2")
system("module load proj/6.2.1")
system("module load hdf5/1.8.18")

###################################################
#Script to visualise RNA Velocity of single cells##
##(using data for each sample/biopsy separately)###
###################################################

#Define some useful functions:

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

Seu2Dataframe = function(x) { 
  
  umap=data.frame(x@reductions[["umap"]]@cell.embeddings)
  data=data.frame(cluster=x$seurat_clusters,cell_type=x$cell_type_level2,UMAP_1=umap$UMAP_1,UMAP_2=umap$UMAP_2)
  return(data)
}


#Load the data (list of Seurat objects constructed from the loom Velocyto output)

seu_list=readRDS("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/NMR_velocyto_16578_17133_seu_list_filtered.rds" )

seu=readRDS("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Integrated_NMR_Doublet_removed.rds")


#Select a group of cells. In this case, we select Epithelial cells
seu_krt=subset(seu,subset= cell_type =="Epithelial cells"  )



vet2=c()

for ( i in 1:length(seu_list)) { 
vet2=unique(c(vet2,unique(seu_list[[i]]$name) ) )
}


vet=levels(factor(seu_krt$cell_type_level2))
ident.colors <- (scales::hue_pal())(n = length(x = vet ))
names(ident.colors) <- vet


#Choose whether to use Seurat's SC Transformation or Log Normalization
SCT=TRUE

for ( i in 1:length(seu_list)) {     
DefaultAssay(seu_list[[i]])="spliced"

cat("analysing sample ",unique(seu_list[[i]]$name),"\n" )

#Normalize and scale (if appropriate) data
if (SCT==FALSE) { 
  seu_list[[i]] <- NormalizeData(seu_list[[i]], normalization.method = "LogNormalize")
  seu_list[[i]]=FindVariableFeatures(seu_list[[i]],selection.method = "vst",nfeatures = 3000)
  seu_list[[i]] <- ScaleData(seu_list[[i]],)
method="logNorm"
} else { 
  seu_list[[i]]=FindVariableFeatures(seu_list[[i]],selection.method = "vst",nfeatures = 3000)
seu_list[[i]]  <- SCTransform(object = seu_list[[i]] , assay = "spliced",verbose = TRUE,variable.features.n = 3000)  
method="SCT"
}

#Clustering and UMAP calculation
seu_list[[i]] <- RunPCA(seu_list[[i]])    #features = VariableFeatures(seu_list[[i]])
seu_list[[i]] <- FindNeighbors(seu_list[[i]], dims = 1:30)
seu_list[[i]] <- FindClusters(seu_list[[i]], resolution = 0.5)
seu_list[[i]] <- RunUMAP(seu_list[[i]], dims = 1:30)
Idents(seu_list[[i]])="cell_type_level2"


#Define color labelling
ident.colors2 <- ident.colors[vet%in%levels(seu_list[[i]])]
levels(seu_list[[i]])=names(ident.colors2)                     

cell.colors <- ident.colors2[Idents(object = seu_list[[i]])]    
names(x = cell.colors) <- colnames(x = seu_list[[i]])
#data_cols=data.frame(colors=factor(ident.colors2),cell_type=factor(names(ident.colors2),levels=names(ident.colors2) )  )
data_cols=data.frame(colors=factor(ident.colors2),cell_type=factor(names(ident.colors2))  )

data=Seu2Dataframe(seu_list[[i]])

#Save a plot on which to run "g_legend"
plot1=print(ggplot(data,aes(x=UMAP_1,y=UMAP_2))+geom_point( aes(color=cell_type), size=0.2, alpha=0.8,shape=1  )+ 
  theme(plot.title = element_text(hjust = 0.5) , panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black") ) )

legend <- g_legend(plot1)


#Run Velocity and save plot
seu_list[[i]] <- RunVelocity(object = seu_list[[i]], deltaT = 1, kCells = 25, fit.quantile = 0.02,verbose = TRUE)

png(paste("~/Velocyto_plot_NMR_kerat_16578_17133_",unique(seu_list[[i]]$name),"_",method,".png",sep=""),res = 600,width=8,height=8,units='in')
print(show.velocity.on.embedding.cor(emb = Embeddings(object = seu_list[[i]], reduction = "umap"), vel = Tool(object = seu_list[[i]], 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1) )
dev.off()


#The default plot is without legend!
#Reload the PNG and use grid.arrange to add the legend, then save the final figure
plot2 <- readPNG(paste("~/Velocyto_plot_NMR_kerat_16578_17133_",unique(seu_list[[i]]$name),"_",method,".png",sep=""))

png(paste("~/Velocyto_plot_NMR_kerat_16578_17133_",unique(seu_list[[i]]$name),"_",method,"_Legend.png",sep=""),res = 600,width=8,height=8,units='in')
print(grid.arrange(rasterGrob(plot2),legend, ncol=2, nrow=1,  widths=c(4/6,2/6) )  )
dev.off()

#Save Seurat object
seu11=seu_list[[i]]
saveRDS(seu11,file=paste("Velocyto_seu_object_",unique(seu_list[[i]]$name),"_",method,".rds",sep=""))


}

