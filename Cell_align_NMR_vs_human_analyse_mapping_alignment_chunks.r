library(Seurat)
library(future)
options(future.globals.maxSize=20000 * 1024^2)
library(dplyr)
library(gridExtra)
library(monocle3)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(scran)
library(batchelor)
library(scater)
library(cellAlign)
system("module load gdal/3.0.2")
system("module load proj/6.2.1")
system("module load hdf5/1.8.18")
system("module load proj/6.2.1")

#load("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_NMR_Human_Kerat_AllTreat_Allclusters_400ptsnoDoublets_sameClusters_CellAlign_global_2000_VarFeatures.RData")

#This script is run on the pseudotime alignment of naked molerat versus human epithelial cells

#Loop across different runs of CellAlign, with different number of points used to interpolate the data
for (numPts in c(400,1000) )  { 

load(paste("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_NMR_Human_Kerat_AllTreat_Allclusters_",numPts,"ptsnoDoublets_sameClusters_CellAlign_global_2000_VarFeatures.RData",sep=""))


#Load Seurat objects
seu_nmr=readRDS("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/shared_data/Seurat_nmr_Keratinocytes_noDoublets.rds")
seu_hsa=readRDS("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Seurat_human_Keratinocytes_noDoublets.rds")

  
#Add cell annotation to the metanode mapping data 
metaNodePt = mapping$metaNodesPt
metaNodePt = metaNodePt[order(metaNodePt$ptQuery),]
metaNodePt$align = 1:nrow(metaNodePt)
metaNodePtLong = melt(metaNodePt[,c('ptQuery','ptRef','align')], id.vars = c('align'))
metaNodePtLong = melt(metaNodePt, id.vars = c('align','metaNodeQuery','metaNodeRef'))  

metaNodePtLong[,"variable"]=gsub("ptQuery","NMR_epithelial" ,metaNodePtLong[,"variable"])
metaNodePtLong[,"variable"]=gsub("ptRef","human_epithelial" ,metaNodePtLong[,"variable"])

typevet=unique(cds2_UN@colData$cell_type_level2)               #UN-->naked molerat
typevet2=unique(cds2_DT@colData$cell_type_level2)              #DT-->human

metaNodePtLong$cell_type_nmr=NA
metaNodePtLong$cell_type_human=NA         

for (i in 1:length(typevet))  { 
  
  x <- colnames( cds2_UN[,pData(cds2_UN)$cell_type_level2==typevet[i]] )
  
  metaNodePtLong$cell_type_nmr[metaNodePtLong$metaNodeQuery%in%x]=typevet[i]
  
}


for (i in 1:length(typevet2))  { 
  
  x <- colnames( cds2_DT[,pData(cds2_DT)$cell_type_level2==typevet2[i]] )
  
  metaNodePtLong$cell_type_human[metaNodePtLong$metaNodeRef%in%x]=typevet2[i]
}


#In every row I have human and nmr annotation. I want to combine them:
metaNodePtLong$cell_type=metaNodePtLong$cell_type_human
metaNodePtLong$cell_type[metaNodePtLong$variable=="NMR_epithelial"]=metaNodePtLong$cell_type_nmr[metaNodePtLong$variable=="NMR_epithelial"]
#metaNodePtLong$cell_type=factor(metaNodePtLong$cell_type)
x=gsub("Spinous/Granular","Spinous/granular",metaNodePtLong$cell_type)
metaNodePtLong$cell_type=factor(x)



#########Find chunks of alignments, where there are changes in the pseudotime expansion-contraction in one species versus the other:

costMat = alignment$localCostMatrix
costMat = t(apply(costMat,1,function(x){return(as.numeric(x))}))
linearInd = sub2ind(nrow(costMat), alignment$align[[1]]$index1, alignment$align[[1]]$index2)


coordMat=data.frame(NMR=alignment$align[[1]]$index1,human=alignment$align[[1]]$index2)


##Define groups by identifying direction changes along the dissimilarity matrix (:alignment)
n1=coordMat[1,1]
n2=coordMat[1,2]
cont=1

coordMat$group=cont 

same_n1=coordMat[2,1]==n1

for (k in 2:length(coordMat[,1])) { 

if ( same_n1 && coordMat[k,1]!=n1  ) { 
  
  cont=cont+1
  coordMat[c(k),"group"]=cont #paste(cont,"node",sep="_")
  
  }  
  else if ( same_n1==FALSE && coordMat[k,2]!=n2  ) { 
    
    cont=cont+1
    coordMat[c(k),"group"]=cont #paste(cont,"node",sep="_")
  } 
  else {
  coordMat$group[k]=cont 
  }
  
  n1=coordMat[k,1]
  n2=coordMat[k,2]
  
  same_n1=coordMat[k-1,1]==n1
  
  }


metaNodePt = metaNodePt[order(metaNodePt$ptQuery),]       #Query is naked molerat


run_loop=TRUE

if (run_loop)  { 
  
cds2_UN@colData$CellAlign_mapping_groups=""
cds2_DT@colData$CellAlign_mapping_groups=""


for (k in 1:length(coordMat[,1])) { 

cat("working on step",k," \n")
cell_nmr=as.character(metaNodePt[k,"metaNodeQuery"])
cell_human=as.character(metaNodePt[k,"metaNodeRef"])
  
if (str_length(cds2_UN@colData$CellAlign_mapping_groups[match(cell_nmr,colnames(cds2_UN)) ])==0 ) { 
cds2_UN@colData$CellAlign_mapping_groups[match(cell_nmr,colnames(cds2_UN)) ]=coordMat[k,"group"]
} else if (sample(0:1,1)==1) { 
  cds2_UN@colData$CellAlign_mapping_groups[match(cell_nmr,colnames(cds2_UN)) ]=coordMat[k,"group"]
  }

##If a label has already been assigned, choose randomly whether to replace it or not:

if (str_length(cds2_DT@colData$CellAlign_mapping_groups[match(cell_human,colnames(cds2_DT)) ])==0 ) { 
cds2_DT@colData$CellAlign_mapping_groups[match(cell_human,colnames(cds2_DT)) ]=coordMat[k,"group"]
} else if (sample(0:1,1)==1) { 
  cds2_DT@colData$CellAlign_mapping_groups[match(cell_human,colnames(cds2_DT)) ]=coordMat[k,"group"]
  }

id_human=cds2_DT@colData$CellAlign_mapping_groups[match(cell_human,colnames(cds2_DT)) ]
id_nmr=cds2_UN@colData$CellAlign_mapping_groups[match(cell_nmr,colnames(cds2_UN)) ]

for (cell in mapping$queryAssign[[cell_nmr]] ) {      
  
  ind_cds=match(cell,colnames(cds2_UN))
    cds2_UN@colData[ind_cds,"CellAlign_mapping_groups"]=id_nmr  #coordMat[k,"group"]
 
 
                                                }

for (cell in mapping$refAssign[[cell_human]] ) {      
  
  ind_cds=match(cell,colnames(cds2_DT))
  cds2_DT@colData[ind_cds,"CellAlign_mapping_groups"]=id_human #coordMat[k,"group"]
 
                                               }

}


cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]=cds2_UN@colData$CellAlign_mapping_groups
cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]=cds2_DT@colData$CellAlign_mapping_groups

x=labels(seu_hsa$orig.ident)
y=colnames(cds2_DT)
seu_hsa@meta.data[x%in%y,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]=cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]

x=labels(seu_nmr$orig.ident)
y=colnames(cds2_UN)
seu_nmr@meta.data[x%in%y,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]=cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]

n=max(coordMat$group)
cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]=factor(cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")], levels=1:n)
cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]=factor(cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")], levels=1:n)

pdf(paste("~/CellAlign_human_nmr_epithelial_AllTreatUMAP_",numPts,"pts_2kVarFeatures_UMAP_human_byAlignGroups_randChoice.pdf",sep=""),width=11,height=7)
plot_cells(cds2_DT,color_cells_by=paste("CellAlign_mapping_groups_",numPts,"randChoice",sep=""), label_leaves=FALSE,label_branch_points=FALSE,cell_size = 1,group_label_size = 5,label_cell_groups = FALSE)+scale_color_manual(values=col[unique(sort(as.numeric((cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"randChoice",sep="")]))))])
dev.off()

            }



#########Alternative chunk assignment (no if conditions)

run_loop=TRUE

if (run_loop)  { 
  
  cds2_UN@colData$CellAlign_mapping_groups=""
  cds2_DT@colData$CellAlign_mapping_groups=""
  
  #cell CATCGGGGTCCTGCTT_2 is problematic...
  #There is an issue in the change point: the cells are labelled with multiple group names. keeping the first label which is assigned
  
  for (k in 1:length(coordMat[,1])) { 
    
    cat("working on step",k," \n")
    cell_nmr=as.character(metaNodePt[k,"metaNodeQuery"])
    cell_human=as.character(metaNodePt[k,"metaNodeRef"])
  
      cds2_UN@colData$CellAlign_mapping_groups[match(cell_nmr,colnames(cds2_UN)) ]=coordMat[k,"group"]
    
      cds2_DT@colData$CellAlign_mapping_groups[match(cell_human,colnames(cds2_DT)) ]=coordMat[k,"group"]
    
    for (cell in mapping$queryAssign[[cell_nmr]] ) {      
      
      ind_cds=match(cell,colnames(cds2_UN))
      
      cds2_UN@colData[ind_cds,"CellAlign_mapping_groups"]=coordMat[k,"group"]
      
    }
    
    for (cell in mapping$refAssign[[cell_human]] ) {      
      
      ind_cds=match(cell,colnames(cds2_DT))
        cds2_DT@colData[ind_cds,"CellAlign_mapping_groups"]=coordMat[k,"group"]
    }
    
  }
  
  cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]=cds2_UN@colData$CellAlign_mapping_groups
  cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]=cds2_DT@colData$CellAlign_mapping_groups
  
  n=max(coordMat$group)
  cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]=factor(cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")], levels=1:n)
  cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]=factor(cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")], levels=1:n)
  
  x=labels(seu_hsa$orig.ident)
  y=colnames(cds2_DT)
  seu_hsa@meta.data[x%in%y,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]=cds2_DT@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]
  
  x=labels(seu_nmr$orig.ident)
  y=colnames(cds2_UN)
  seu_nmr@meta.data[x%in%y,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]=cds2_UN@colData[,paste("CellAlign_mapping_groups_",numPts,"lastChoice",sep="")]
   
}

}

saveRDS(cds2_DT,file="/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/monocle3/DT_UN/Monocle3_human_Kerat_noDoublets_sameClusters_UN_sce.rds")
saveRDS(cds2_UN,file="/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/monocle3/AllTreat/Monocle3_NMR_Kerat_noDoublets_sameClusters_allTreatmentsUMAP_UN_sce.rds")
saveRDS(seu_nmr,file="/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/shared_data/Seurat_nmr_Keratinocytes_noDoublets.rds")
saveRDS(seu_hsa,file="/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Seurat_human_Keratinocytes_noDoublets.rds")






