
library(Seurat)
library(future)
options(future.globals.maxSize=20000 * 1024^2)
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
library(cellAlign)
library(viridis)

setwd("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/")


#The value of this variable determines whether to exclude proliferative basal or not
name="AllTreat_Allclusters_noProlif"

#labels "DT" and "UN" in the variable names refer to the different datasets:
#query is labelled as UN, ref is labelled as DT

cds2_UN=readRDS("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/monocle3/AllTreat/Monocle3_NMR_Kerat_noDoublets_sameClusters_allTreatmentsUMAP_UN_sce.rds")

if (  length(grep("noProlif",name))==1  ) { 
cds2_UN <- cds2_UN[,pData(cds2_UN)$cell_type_level2!="Proliferative basal"]
}

cds2_DT=readRDS(file = "/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/monocle3/DT_UN/Monocle3_human_Kerat_noDoublets_sameClusters_UN_sce.rds" )


expGlobalUN=data.frame(logcounts(cds2_UN))
expGlobalDT=data.frame(logcounts(cds2_DT))
expLocalUN=data.frame(logcounts(cds2_UN))
expLocalDT=data.frame(logcounts(cds2_DT))



##Remove naming differences between datasets and exclude specific lines (the alignment will run on highly variable genes only)

row.names(expGlobalUN)=gsub("^MITO-", "MT-", row.names(expGlobalUN))

row.names(expGlobalUN)=chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", row.names(expGlobalUN))
row.names(expGlobalDT)=chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", row.names(expGlobalDT))

row.names(expLocalUN)=gsub("^MITO-", "MT-", row.names(expLocalUN))
row.names(expLocalUN)=chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", row.names(expLocalUN))
row.names(expLocalDT)=chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", row.names(expLocalDT))

trajDT=pseudotime(cds2_DT)
trajUN=pseudotime(cds2_UN)

trajDT = as.numeric(trajDT)
names(trajDT) = colnames(expGlobalDT)
trajUN = as.numeric(trajUN)
names(trajUN) = colnames(expGlobalUN)

seu=readRDS("Seurat_nmr_Keratinocytes_noDoublets.rds")

ngenes=2000

seu=FindVariableFeatures(seu,assay="RNA",nfeatures = ngenes)

geneset=paste(ngenes,"_VarFeatures",sep="") #"monocle_graphTest"

#genes1=read.table(file="~/Monocle3_DEgenes_graph_test_q001_NMR_Kerat_noDoublets_sameClusters_UN_ids.csv", sep="," )
#genes2=read.table(file="~/Monocle3_DEgenes_graph_test_q001_NMR_Kerat_noDoublets_sameClusters_DT_ids.csv", sep="," )

genes=VariableFeatures(seu,assay="RNA") #as.vector(unlist(unique(c(genes1,genes2))))

genes=gsub("^MITO-", "MT-", genes)
genes=chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", genes)

seu=readRDS("Seurat_human_Keratinocytes_noDoublets.rds")

seu=FindVariableFeatures(seu,assay="RNA",nfeatures = ngenes)

genes2=VariableFeatures(seu,assay="RNA")

genes2=chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", genes2)

genes=unique(c(genes,genes2))

genes=chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", genes)


expGlobalDT=expGlobalDT[ genes, ]     
expGlobalUN=expGlobalUN[ genes, ]
expLocalDT=expLocalDT[ genes, ]
expLocalUN=expLocalUN[ genes, ]

numPts = 200
interGlobalDT = cellAlign::interWeights(expDataBatch = expGlobalDT, trajCond = trajDT,
                                         winSz = 0.1, numPts = numPts)
interGlobalUN = cellAlign::interWeights(expDataBatch = expGlobalUN, trajCond = trajUN,
                                         winSz = 0.1, numPts = numPts)

require(ggplot2)
require(reshape2)
require(pheatmap)
sharedMarkers = intersect(rownames(expGlobalDT), rownames(expGlobalUN))
whichgene=sharedMarkers[1]
selectedDT<-interGlobalDT$interpolatedVals[whichgene,]
selectedUN<-interGlobalUN$interpolatedVals[whichgene,]

dfDTi = data.frame(traj = interGlobalDT$traj, value=(selectedDT), error=interGlobalDT$error[whichgene,])
dfDT = data.frame(traj = trajDT, t(expGlobalDT[whichgene,]))
dfUNi = data.frame(traj = interGlobalUN$traj, value=(selectedUN), error=interGlobalUN$error[whichgene,])
dfUN = data.frame(traj = trajUN, t(expGlobalUN[whichgene,]))
dfDTM = melt(dfDT, id.vars = 'traj')
dfUNM = melt(dfUN, id.vars = 'traj')
#plot of an example gene and its interpolation with error bars

png(paste("~/Seurat_monocle_",name,"_NMR_human_kerat_pseudotime_interpolation_",geneset,".png",sep=""),res = 600,width=8,height=7,units='in')
print(ggplot(dfDTi, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=dfDTM, aes(x=traj,y=value)) + ggtitle(whichgene) )
dev.off()

#scale the interpolated data (Recommended):
interScaledGlobalDT = cellAlign::scaleInterpolate(interGlobalDT)
interScaledGlobalUN = cellAlign::scaleInterpolate(interGlobalUN)

A=calcDistMat(interScaledGlobalUN$scaledData[,1:10],interScaledGlobalDT$scaledData[,1:10], dist.method = 'Euclidean')

png(paste("~/Seurat_monocle_",name,"_NMR_human_kerat_pseudotime_dist_matrix",geneset,".png",sep=""),res = 600,width=8,height=7,units='in')
print(pheatmap(A, cluster_cols = F, cluster_rows=F, main = "Human vs NMR distances, 1st 10 points",
         show_rownames = F, show_colnames = F,display_numbers = TRUE) )
dev.off()

#perform global alignment of all genes:
alignment = globalAlign(interScaledGlobalUN$scaledData, interScaledGlobalDT$scaledData,
                        scores = list(query = interScaledGlobalUN$traj, 
                                      ref = interScaledGlobalDT$traj),
                        sigCalc = F, numPerm = 20)

png(paste("~/Seurat_monocle_",name,"_NMR_human_kerat_pseudotime_globAlign_alignment",geneset,".png",sep=""),res = 600,width=8,height=7,units='in')
print(plotAlign(alignment))
dev.off()

mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalUN$traj, realTrajQuery = trajUN,
                            intTrajRef = interScaledGlobalDT$traj, realTrajRef = trajDT)

png(paste("~/Seurat_monocle_",name,"_NMR_human_kerat_pseudotime_globAlign_mapping",geneset,".png",sep=""),res = 600,width=8,height=7,units='in')
print(plotMapping(mapping))
dev.off()

rm(seu)

save.image(paste("Monocle3_NMR_Human_Kerat_",name,"noDoublets_sameClusters_CellAlign_global_",geneset,".RData",sep=""))


#####local alignment:

local=FALSE

if (local==TRUE)      {  

Thresh=0.2
numPts = 200
interLocalDT = interWeights(expDataBatch = expLocalDT, trajCond = trajDT, winSz = 0.1, numPts = numPts)
interLocalUN = interWeights(expDataBatch = expLocalUN, trajCond = trajUN, winSz = 0.1, numPts = numPts)
interScaledLocalDT = cellAlign::scaleInterpolate(interLocalDT)
interScaledLocalUN = cellAlign::scaleInterpolate(interLocalUN)

A=calcDistMat(interScaledLocalUN$scaledData,interScaledLocalDT$scaledData, dist.method = 'Euclidean')
A[A > 10*Thresh] <- max(A)
alignment = localAlign(interScaledLocalUN$scaledData,interScaledLocalDT$scaledData,threshPercent = Thresh)

costMat = t(apply(A,1,function(x){return(as.numeric(x))}))
linearInd = cellAlign::sub2ind(nrow(A), alignment$align[[1]]$index1, alignment$align[[1]]$index2)
costMat[linearInd] = NA
costMat = data.frame(costMat, row.names=1:numPts)
colnames(costMat) = 1:numPts
pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
         main = 'gated search region',
         show_rownames = F, show_colnames = F)

plotAlign(alignment)

BDT=colMeans(interScaledLocalDT$scaledData)
BUN=colMeans(interScaledLocalUN$scaledData)
p=unique(alignment$align[[1]]$index1)
q=unique(alignment$align[[1]]$index2)
plot(1:200,BUN,xlab = "pseudotime",ylab = "mean interpolated expression", main = "unaligned mean expression",ylim = c(0,1.1))
points(p,BUN[p],col="red")
points(1:200,BDT,col="grey60")
points(q,BDT[q],col="red")
text(90,1,"DT")
text(150,1,"UN")
text(125,.3,"red points are conserved")



}


rm(cds2_UN)
rm(cds2_DT)

data=data.frame(index1=alignment$align[[1]]$index1,index2=alignment$align[[1]]$index2)


plot_data=FALSE

#cds2_UN=readRDS(file = "/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_NMR_Kerat_noDoublets_sameClusters_removed3Clusters_UN_sce.rds" )
#cds2_DT=readRDS(file = "/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_human_Kerat_noDoublets_sameClusters_UN_sce.rds" )

#query=UN

if (plot_data==TRUE)  { 
  
  
  metaNodePt = mapping$metaNodesPt
  metaNodePt = metaNodePt[order(metaNodePt$ptQuery),]
  metaNodePt$align = 1:nrow(metaNodePt)
  metaNodePtLong = melt(metaNodePt[,c('ptQuery','ptRef','align')], id.vars = c('align'))
  metaNodePtLong = melt(metaNodePt, id.vars = c('align','metaNodeQuery','metaNodeRef'))  
  
  metaNodePtLong[,"variable"]=gsub("ptQuery","NMR_epithelial" ,metaNodePtLong[,"variable"])
  metaNodePtLong[,"variable"]=gsub("ptRef","human_epithelial" ,metaNodePtLong[,"variable"])
  
  typevet=unique(cds2_UN@colData$cell_type_level2)
  typevet2=unique(cds2_DT@colData$cell_type_level2)
  
  metaNodePtLong$cell_type_nmr=NA
  metaNodePtLong$cell_type_human=NA         #do I need another column?
  
  for (i in 1:length(typevet))  { 
    
    x <- colnames( cds2_UN[,pData(cds2_UN)$cell_type_level2==typevet[i]] )
    
    metaNodePtLong$cell_type_nmr[metaNodePtLong$metaNodeQuery%in%x]=typevet[i]
    
  }
  
  
  for (i in 1:length(typevet2))  { 
    
    x <- colnames( cds2_DT[,pData(cds2_DT)$cell_type_level2==typevet2[i]] )
    
    metaNodePtLong$cell_type_human[metaNodePtLong$metaNodeRef%in%x]=typevet2[i]
  }
  
  #In every row I have human and nmr annotation. I mix them based on column "variable":
  metaNodePtLong$cell_type=metaNodePtLong$cell_type_human
  
  metaNodePtLong$cell_type[metaNodePtLong$variable=="NMR_epithelial"]=metaNodePtLong$cell_type_nmr[metaNodePtLong$variable=="NMR_epithelial"]
  
  x=gsub("Spinous/Granular","Spinous/granular",metaNodePtLong$cell_type)
  metaNodePtLong$cell_type=factor(x)
  
  pdf(paste("~/Seurat_monocle_",name,"_NMR_human_kerat_pseudotime_globAlign_mapping",geneset,"_2.pdf",sep=""),width=12,height=5)
  ggplot(metaNodePtLong[metaNodePtLong$variable!="diff",], aes(x = variable, y = value, group = align)) + geom_line(aes(color=cell_type_nmr)) + theme_bw() + 
    geom_point(aes(color=cell_type)) +
    coord_flip() + ggtitle('Global alignment of human and nmr pseudotime')+theme(plot.title = element_text(hjust = 0.5) )
dev.off()
  ######find expression of K14
  c=as.data.frame(logcounts(cds2_DT))
  
  cont=1
  metaNodePtLong$metaNodeRef_KRT14=NA
  for (g in as.vector(metaNodePtLong$metaNodeRef)  )  {       #labels(mapping$refAssign)
    
  x=match(g, colnames( cds2_DT ) )  
    
  value=c[match("KRT14",row.names(c)),x]
  metaNodePtLong$metaNodeRef_KRT14[cont]=value
  cont=cont+1
    }
  
  c=as.data.frame(logcounts(cds2_UN))  ######??
  cont=1
  metaNodePtLong$metaNodeQuery_KRT14=NA
  for (g in as.vector(metaNodePtLong$metaNodeQuery)  )  {       #labels(mapping$refAssign)
    
    x=match(g, colnames( cds2_UN ) )  
    
    value=c[match("KRT14",row.names(c)),x]
    metaNodePtLong$metaNodeQuery_KRT14[cont]=value
    cont=cont+1
  }
  
  #create a column where I store expression from human AND NMR 
  metaNodePtLong$KRT14=metaNodePtLong$metaNodeRef_KRT14
  metaNodePtLong$KRT14[metaNodePtLong$variable=="NMR_epithelial"]=metaNodePtLong$metaNodeQuery_KRT14[metaNodePtLong$variable=="NMR_epithelial"]
  
  
  pdf(paste("~/Seurat_monocle_",name,"_NMR_human_kerat_pseudotime_globAlign_mapping",geneset,"_K14.pdf",sep=""),width=12,height=5)
  ggplot(metaNodePtLong[metaNodePtLong$variable!="diff",], aes(x = variable, y = value, group = align)) + 
  geom_line(aes(color=cell_type_nmr)) + theme_bw() + ggnewscale::new_scale("color")+
  geom_point(aes(color=KRT14,alpha=0.5),shape=1) +
  coord_flip() + ggtitle('Global alignment of human and nmr pseudotime')+theme(plot.title = element_text(hjust = 0.5) )+scale_color_gradient2(low = "green", mid = "yellow",high="red",midpoint = 3)
  dev.off()
  ####################
  
  scaled_traj=interGlobalUN$traj
  
  cell_anno=c()
  
  values=metaNodePtLong[metaNodePtLong$variable=="NMR_epithelial","value"]
  anno=factor(metaNodePtLong[metaNodePtLong$variable=="NMR_epithelial","cell_type"])
  for (i in 1:200)  { 
    
  x=abs(values-scaled_traj[i]) 
  ind=match(min(x),x) 
  cell_anno[i]=as.character(anno[ind])
  
   }
  
  
  cell_anno_nmr=cell_anno

  
  scaled_traj=interGlobalDT$traj
  
  cell_anno=c()
  
  values=metaNodePtLong[metaNodePtLong$variable=="human_epithelial","value"]
  anno=factor(metaNodePtLong[metaNodePtLong$variable=="human_epithelial","cell_type"])
  for (i in 1:200)  { 
    
    x=abs(values-scaled_traj[i]) 
    ind=match(min(x),x) 
    cell_anno[i]=as.character(anno[ind])
    
  }
  
  cell_anno_human=cell_anno  
  
  ########################
  
  plotAlign <- function(alignment){
    costMat = alignment$localCostMatrix
    costMat = t(apply(costMat,1,function(x){return(as.numeric(x))}))
    linearInd = sub2ind(nrow(costMat), alignment$align[[1]]$index1, alignment$align[[1]]$index2)
    costMat[linearInd] = NA
    costMat = data.frame(costMat, row.names=1:nrow(costMat))
    colnames(costMat) = 1:ncol(costMat)
    #for global alignment, where there is a pseudotime shift vector:
    if(!is.null(alignment$ptShift)){
      
      annotRows = data.frame(nmr_cell_type=factor(cell_anno_nmr,levels=unique(sort(cell_anno_nmr))),ptShift = abs(alignment$ptShift), sign = factor(sign(alignment$ptShift)),row.names = colnames(costMat))
      annotCols = data.frame(human_cell_type=factor(cell_anno_human,levels=unique(sort(cell_anno_human))),row.names = row.names(costMat))
      
      my_colour = list(
        #nmr_cell_type = c("Basal" = "#F8766D", "Granular" = "#A3A500", "Migrating spinous" = "#00BF7D", "Spinous" = "#00B0F6", "Spinous/granular"= "#E76BF3" ),
        nmr_cell_type = c("Basal" = "#F8766D", "Spinous" = "#00B0F6","Proliferative basal" =  "#CD9600" ,"Granular" = "#A3A500",  "Migrating spinous" = "#00BF7D", "Spinous/granular"= "#E76BF3" ),
        # human_cell_type = c("Basal" = "#F8766D", "Granular" = "#A3A500","Spinous" = "#00B0F6", "Spinous/granular"= "#E76BF3" )
        human_cell_type = c("Basal" = "#F8766D", "Spinous" = "#00B0F6", "Spinous/granular"= "#E76BF3","Sweat glands" = "blue"  ,"Undefined Epithelial cells (KRT16+)"="#00BE67" )
        
        )
      
      
      pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
               main = 'alignment plot',
               show_rownames = F, show_colnames = F, annotation_col = annotCols, annotation_row = annotRows, annotation_colors = my_colour)
    }else{
      pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
               main = 'alignment plot', show_rownames = F, show_colnames = F)
    }
    
    return(NA)
  }
  
  
  
  }

