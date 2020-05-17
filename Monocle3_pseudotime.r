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

system("module load gdal/3.0.2")
system("module load proj/6.2.1")
system("module load hdf5/1.8.18")
system("module load proj/6.2.1")

name="NMR_Kerat_noDoublets_sameClusters"      #UMAP calculated on DT+UN only

#load Seurat object

seu=readRDS("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Integrated_NMR_Doublet_removed.rds")

DefaultAssay(seu)="RNA"
seu=NormalizeData(seu)

#select subset of interest. Here we select treated and untreated epithelial cells:

seu_krt=subset(seu,subset= cell_type =="Epithelial cells" & (seu$group=="DT" | seu$group=="UN") )

DefaultAssay(seu_krt)="integrated"
seu_krt=RunPCA(seu_krt,dims = 1:30)
seu_krt <- FindNeighbors(seu_krt, dims = 1:30)
seu_krt=FindClusters(seu_krt,resolution = 0.5)
seu_krt=RunUMAP(seu_krt,dims = 1:30)

pca=seu_krt@reductions[["pca"]]@cell.embeddings#seu_krt@reductions$pca
umap=seu_krt@reductions[["umap"]]@cell.embeddings#seu_krt@reductions$umap
colnames(umap)=c("V1","V2")
seu_sce=as.SingleCellExperiment(seu_krt,assay = "RNA")
seu.data=seu_krt[["RNA"]]@counts

cds=new_cell_data_set(expression_data = seu.data,cell_metadata = seu_sce@colData,gene_metadata = data.frame(row.names =row.names(seu.data), gene_short_name=row.names(seu.data)) )
counts(cds)<-as(counts(cds),"dgCMatrix")
size.factors <- sizeFactors(cds) #libsizes/mean(libsizes)


reducedDim(cds, "PCA") <- as.matrix(pca)
reducedDim(cds, "UMAP") <- as.matrix(umap)


condition="DT"     #"UN" 

saveRDS(seu_krt,file = paste("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_",name,"_krt_DT_UN.rds", sep="" )  )


seu2=subset(seu_krt,subset= group==condition)
pca=seu2@reductions[["pca"]]@cell.embeddings 
umap=seu2@reductions[["umap"]]@cell.embeddings   

colnames(umap)=c("V1","V2")


#select cells of the desired treatment/group
cds2 <- cds[,pData(cds)$group==condition]

cds2=logNormCounts(cds2)
cds2 <- preprocess_cds(cds2, num_dim = 100)
cds2 <- reduce_dimension(cds2,reduction_method = "UMAP")

reducedDim(cds2, "PCA") <- as.matrix(pca)
reducedDim(cds2, "UMAP") <- as.matrix(umap)

saveRDS(cds2,file = paste("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_",name,"_",condition,"_sce.rds", sep="" )  )


cds2 = cluster_cells(cds2, reduction_method = c("UMAP"), resolution=1e-04) 


###########################################
# This analysis relies on Seurat clusters # 
###########################################

#Transfer the cluster annotation of the Seurat object

cds2@clusters$UMAP$clusters=seu2$seurat_clusters 


#find trajectory

cds2=learn_graph(cds2, use_partition = TRUE, close_loop = TRUE,learn_graph_control = NULL, verbose = TRUE)



#choose the starting point(s) for the pseudotime score:
cds2 <- order_cells(cds2)

saveRDS(cds2,file = paste("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_",name,"_",condition,"_sce.rds", sep="" )  )

cds2_pr_test_res=graph_test(cds2, neighbor_graph="principal_graph", cores=10)
write.table(cds2_pr_test_res,file=paste("~/Monocle3_DEgenes_graph_test_",name,"_",condition,".csv",sep=""), sep="," )

pr_deg_ids <- row.names(subset(cds2_pr_test_res, q_value <= 0.01 & abs(morans_test_statistic)>50 ) )
write.table(pr_deg_ids,file=paste("~/Monocle3_DEgenes_graph_test_q001_morans50_",name,"_",condition,"_ids.csv",sep="") ,sep="," )

pr_deg_ids <- row.names(subset(cds2_pr_test_res, q_value <= 0.01  ) )
write.table(pr_deg_ids,file=paste("~/Monocle3_DEgenes_graph_test_q001_",name,"_",condition,"_ids.csv",sep="") ,sep="," )


#identify gene modules:
gene_module_df <- find_gene_modules(cds2[pr_deg_ids,], resolution=0.0001 )  


saveRDS(gene_module_df,file = paste("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_",name,"_",condition,"_gene_modules.rds", sep="" )  )
#gene_module_df=readRDS(paste("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_",name,"_",condition,"_gene_modules.rds", sep="" )  )


write.table(gene_module_df,file=paste("~/Monocle3_DEgenes_graph_test_q001_modules_e04_",name,"_",condition,".csv",sep=""), sep="," )


cell_group_df <- tibble::tibble(cell=row.names(colData(cds2)), 
                                cell_group=cds2@colData$cell_type_level2  ) 


agg_mat <- aggregate_gene_expression(cds2, gene_module_df,cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))


png(paste("~/Seurat_",name,"_",condition,"_monocle_heatmap_gene_modules2.png",sep=""),res = 600,width=6,height=5.5,units='in')
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=8)
dev.off()


AFD_lineage_cds <- cds2[rowData(cds2)$gene_short_name %in% AFD_genes,]

png(paste("~/Seurat_",name,"_",condition,"_monocle_umap_clusters.png",sep=""),res = 600,width=5,height=5,units='in')
plot_cells(cds2,
           color_cells_by = "cluster",
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4,cell_size = 0.5)
dev.off()


png(paste("~/Seurat_",name,"_",condition,"_monocle_umap_pseudotime.png",sep=""),res = 600,width=6.5,height=5,units='in')
plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4,cell_size = 0.5)
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_monocle_genes_pseudotime.png",sep=""),res = 600,width=8,height=7,units='in')
monocle3::plot_genes_in_pseudotime(AFD_lineage_cds,min_expr=0.5)
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_monocle_umap_expression_byModule.png",sep=""),res = 600,width=8,height=4,units='in')
plot_cells(cds2,
           genes=gene_module_df %>% filter(module %in% c(1,2,3,4,5)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           group_label_size=4,cell_size = 0.5)
dev.off()

