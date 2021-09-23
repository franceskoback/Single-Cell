#!/usr/bin/bash
module load CBI r
Rscript Step1_Preprocessing.R /wynton/home/srivastava/franceskoback/Single-Cell/filtered_feature_bc_matrix /wynton/home/srivastava/franceskoback/Single-Cell/data/AggrSheet.csv /wynton/home/srivastava/franceskoback/Single-Cell/data/regev_lab_cell_cycle_genes.txt /wynton/home/srivastava/franceskoback/Single-Cell/data/noFilters_scoresAdded.rds
Rscript Visualize_Boundaries.R /wynton/home/srivastava/franceskoback/SRA_Analysis/data/noFilters_scoresAdded.rds 1 10 13 33 2500 23000 1000 5000
Rscript Step2_Clustering_Filtering.R /wynton/home/srivastava/franceskoback/SRA_Analysis/data/noFilters_scoresAdded.rds 1 10 13 33 2500 23000 1000 5000 /wynton/home/srivastava/franceskoback/SRA_Analysis/data/rds/FilteredAndClustered_onlyVarGenes.rds /wynton/home/srivastava/franceskoback/SRA_Analysis/data/rds/FilteredAndReClustered_onlyVarGenes.rds /wynton/home/srivastava/franceskoback/SRA_Analysis/data/csv/FindAllMarkers_harmony_res25e-2_05-24-2021.csv



echo Finished
