# Single_cell_RNAseq
Selection of R scripts for the analysis of 10X single cell RNAseq data

*Cell_align_NMR_vs_human.r* : perform an alignment of pseudotime trajectories using Cellalign

*Cell_align_NMR_vs_human_analyse_mapping_alignment_chunks.r* : identify alignment regions of pseudotime expansion/contraction between query and reference (for example, a small pseudotime region in the query aligned to a larger region in the reference)

*Integrate_3species_single_samples_rPCA.r*: perform integration of 3 Seurat objects

*Monocle3_pseudotime.r*: calculate Monocle3 pseudotime trajectories

*EMT_analysis.r* : calculate a score for sets of genes associated with epithelial and mesenchymal cell types. 

*run_velocyto_NMR_EpithelialCells_perBiopsy.r* : calculate and visualise RNA velocity
