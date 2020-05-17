#Script used to study epithelial mesenchymal transition in single cell RNAseq data. 
#The code is used to calculate different scores for sets of genes associated with epithelial and mesenchymal cell types. 

name="NMR_Kerat_noDoublets_sameClusters"
condition="UN"
cds2=readRDS(file = paste("/icgc/dkfzlsdf/analysis/B210/Luca/Rstudio/SCT_singleSamples/Monocle3_",name,"_",condition,"_sce.rds", sep="" ) )

markers_M=c("VIM","S100A4","ACTA2","FN1","CDH2","COL1A1","COL1A2","COL3A1")
markers_EMT=c("ACVR1", "AGER", "ALX1", "AXIN2", "BAMBI", "BCL9L", "BMP2", "BMP4", "BMP7", "COL1A1", "CRB2", "CTNNB1", "DAB2", "ENG", "EZH2", "FOXC1", "GCNT2", "GLIPR2", "HDAC2", "IL1B", "IL6", "ISL1", "JAG1", "LEF1", "LOXL2", "MDK", "MTOR", "MYOCD", "NOTCH1", "PDPN", "RGCC", "SDCBP", "SERPINB3D", "SMAD2", "SMAD3", "SMAD4", "SNAI1", "TCF7L2", "TGFB1", "TGFB1I1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TWIST1", "WWTR1", "ZFP703")
markers_E=c("EPCAM", "CDH1", "OCLN", "TJP1", "TJP2", "TJP3", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "CLDN1", "CLDN10", "CLDN11", "CLDN12", "CLDN13", "CLDN14", "CLDN15", "CLDN16", "CLDN17", "CLDN18", "CLDN19", "CLDN2", "CLDN20", "CLDN22", "CLDN23", "CLDN24", "CLDN3", "CLDN34A", "CLDN34B1", "CLDN34B2", "CLDN34B3", "CLDN34B4", "CLDN34C1", "CLDN34C4", "CLDN34D", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8", "CLDN9", "CLDND1", "CLDND2", "KRT1", "KRT10", "KRT12", "KRT13", "KRT14", "KRT15", "KRT16", "KRT17", "KRT18", "KRT19", "KRT2", "KRT20", "KRT222", "KRT23", "KRT24", "KRT25", "KRT26", "KRT27", "KRT28", "KRT31", "KRT32", "KRT33A", "KRT33B", "KRT34", "KRT35", "KRT36", "KRT39", "KRT4", "KRT40", "KRT42", "KRT5", "KRT6A", "KRT6B", "KRT7", "KRT71", "KRT72", "KRT73", "KRT74", "KRT75", "KRT76", "KRT77", "KRT78", "KRT79", "KRT8", "KRT80", "KRT81", "KRT82", "KRT83", "KRT84", "KRT85", "KRT86", "KRT9")

markers_stem=c("APC", "ARID1A", "ASCL2", "ASPM", "BCL9", "BCL9L", "BMP7", "BMPR1A", "BRAF", "CDC73", "CDH2", "CDKN2A", "CDX2", "CNOT1", "CNOT2", "CNOT3", "CREBBP", "CTNNB1", "CTR9", "CUL4A", "DDX6", "DICER1", "DIS3L2", "DLL1", "DPPA2", "EIF4E", "EIF4ENIF1", "ELAVL1", "ELF5", "EOMES", "ERDR1", "ESRRB", "FANCC", "FANCD2", "FGF10", "FGF4", "FGFR1", "FGFR3", "FOXO1", "FOXO3", "FUT10", "FZD7", "GATA2", "GNL3", "GSDMA3", "HES1", "HES5", "HMGA2", "HNF1B", "HOOK3", "IGF2BP1", "JMJD1C", "KAT2A", "KAT6A", "KDM2B", "KDM3A", "KDM4C", "KIT", "KLF10", "KLF4", "LBH", "LDB1", "LDB2", "LEO1", "LIF", "LIG4", "LIN28A", "LOXL2", "LRP5", "LSM1", "MAPK8", "MCPH1", "MED10", "MED12", "MED14", "MED15", "MED17", "MED21", "MED24", "MED27", "MED28", "MED30", "MED6", "MED7", "METTL14", "METTL3", "MIAT", "MIR294", "MMP24", "MTF2", "MYC", "NANOG", "NANOS2", "NCOA3", "NFYA", "NIPBL", "NKAP", "NODAL", "NOG", "NR2E1", "PADI4", "PAF1", "PANCT2", "PAX2", "PAX8", "PCM1", "PELO", "PHF19", "PIWIL2", "PLA2G2A", "POU5F1", "PRAMEL7", "PRDM14", "PRDM16", "PROX1", "PRRX1", "PTN", "RAF1", "RBPJ", "REST", "RIF1", "RTF1", "SALL1", "SALL4", "SAV1", "SCT", "SETD1A", "SETD6", "SFRP1", "SIX2", "SKI", "SMARCA4", "SMC1A", "SMC3", "SMO", "SOX2", "SOX4", "SOX9", "SPI1", "SRRT", "SS18", "STAG2", "STAT3", "SYMBOL", "TAF5L", "TAF6L", "TAL1", "TBX3", "TCF7L1", "TCF7L2", "TCL1", "TEAD1", "TEAD3", "TEAD4", "TET1", "TFAP2C", "TPT1", "TRIM8", "TRP63", "TUT4", "VANGL2", "VPS72", "WDR43", "WDR62", "WNT7A", "WNT9B", "YAP1", "ZC3H13", "ZFP322A", "ZFP36L2", "ZFP706", "ZHX2")

c=as.data.frame(logcounts(cds2))

m=data.frame(matrix(nrow=0,ncol=length(c[1,])))
 
  
y=c[row.names(c)%in%markers_E,]
means_E=colMeans2(as.matrix(y))
cds2@colData$E_score=means_E


y=c[row.names(c)%in%markers_M,]
means_M=colMeans2(as.matrix(y))
cds2@colData$M_score=means_M

y=c[row.names(c)%in%markers_EMT,]
means_EMT=colMeans2(as.matrix(y))
cds2@colData$EMT_score=means_EMT

library(viridis)

png(paste("~/Seurat_",name,"_",condition,"_E_score_byCluster.png",sep=""),res = 600,width=10,height=6,units='in')
plotColData(cds2, x="seurat_clusters", y="E_score", colour_by="E_score")
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_monocle_umap_E_score.png",sep=""),res = 600,width=6,height=5.5,units='in')
plot_cells(cds2,color_cells_by="E_score" , group_cells_by = "partition", label_groups_by_cluster=TRUE,label_leaves=FALSE,label_branch_points=FALSE,label_cell_groups = FALSE)
dev.off()
 
png(paste("~/Seurat_",name,"_",condition,"_monocle_umap_M_score.png",sep=""),res = 600,width=6,height=5.5,units='in')
plot_cells(cds2,color_cells_by="M_score" , group_cells_by = "partition", label_groups_by_cluster=TRUE,label_leaves=FALSE,label_branch_points=FALSE,label_cell_groups = FALSE)
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_M_score_byCluster.png",sep=""),res = 600,width=10,height=6,units='in')
plotColData(cds2, x="seurat_clusters", y="M_score", colour_by="M_score")
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_E_score_byCluster_M_score.png",sep=""),res = 600,width=10,height=6,units='in')
plotColData(cds2, x="seurat_clusters", y="E_score", colour_by="M_score")
dev.off()

data=data.frame(M_score=cds2@colData$M_score, E_score=cds2@colData$E_score,clusters=cds2@colData$seurat_clusters,cell_type=cds2@colData$cell_type_level2,pseudotime=pseudotime(cds2))
data2=melt(data,id.vars = c("clusters","pseudotime","cell_type"))
data2$variable=factor(data2$variable)

colnames(data2)=c("clusters", "pseudotime","cell_type","group","score")

png(paste("~/Seurat_",name,"_",condition,"_E_score_M_score_new_point_byPseudotime.png",sep=""),res = 600,width=10,height=6,units='in')
ggplot(data2,aes(x=pseudotime,y=score))+geom_point(alpha=0.5,aes(color=group))+geom_smooth(alpha=1,aes(color=group))
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_plot_M_E_score_new_byCluster.png",sep=""),res = 600,width=10,height=6,units='in')
ggplot(data,aes(x=M_score,y=E_score,color=clusters))+geom_point()+geom_abline(intercept = 0,slope = 1, linetype="dotted") #+geom_smooth(aes(x=M_score,y=E_score),color="black" )
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_plot_M_E_score_new_byPseudotime.png",sep=""),res = 600,width=10,height=6,units='in')
ggplot(data,aes(x=M_score,y=E_score,color=pseudotime))+geom_point()+geom_abline(intercept = 0,slope = 1, linetype="dotted")+scale_color_viridis(option = "inferno")
dev.off()


#EMT score

png(paste("~/Seurat_",name,"_",condition,"_monocle_umap_EMT_score.png",sep=""),res = 600,width=6,height=5.5,units='in')
plot_cells(cds2,color_cells_by="EMT_score" , group_cells_by = "partition", label_groups_by_cluster=TRUE,label_leaves=FALSE,label_branch_points=FALSE,label_cell_groups = FALSE)
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_EMT_score_byCluster.png",sep=""),res = 600,width=10,height=6,units='in')
plotColData(cds2, x="seurat_clusters", y="EMT_score", colour_by="EMT_score")
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_E_score_byCluster_EMT_score.png",sep=""),res = 600,width=10,height=6,units='in')
plotColData(cds2, x="seurat_clusters", y="E_score", colour_by="EMT_score")
dev.off()

data=data.frame(EMT_score=cds2@colData$EMT_score,E_score=cds2@colData$E_score,clusters=cds2@colData$seurat_clusters,cell_type=cds2@colData$cell_type_level2,pseudotime=pseudotime(cds2))
data2=melt(data,id.vars = c("clusters","pseudotime","cell_type"))
data2$variable=factor(data2$variable)

colnames(data2)=c("clusters", "pseudotime","cell_type","group","score")

png(paste("~/Seurat_",name,"_",condition,"_E_score_EMT_score_new_point_byPseudotime.png",sep=""),res = 600,width=10,height=6,units='in')
ggplot(data2,aes(x=pseudotime,y=score))+geom_point(alpha=0.5,aes(color=group))+geom_smooth(alpha=1,aes(color=group))
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_plot_M_E_score_new_byCluster.png",sep=""),res = 600,width=10,height=6,units='in')
ggplot(data,aes(x=EMT_score,y=E_score,color=clusters))+geom_point()+geom_abline(intercept = 0,slope = 1, linetype="dotted") #+geom_smooth(aes(x=EMT_score,y=E_score),color="black" )
dev.off()

png(paste("~/Seurat_",name,"_",condition,"_plot_M_E_score_new_byPseudotime.png",sep=""),res = 600,width=10,height=6,units='in')
ggplot(data,aes(x=EMT_score,y=E_score,color=pseudotime))+geom_point()+geom_abline(intercept = 0,slope = 1, linetype="dotted")+scale_color_viridis(option = "inferno")
dev.off()
