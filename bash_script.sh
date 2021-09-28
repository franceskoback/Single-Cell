#!/usr/bin/bash
module load CBI r

Rscript Step1_Preprocessing.R /wynton/home/srivastava/franceskoback/Single-Cell/filtered_feature_bc_matrix /wynton/home/srivastava/franceskoback/Single-Cell/data/AggrSheet.csv /wynton/home/srivastava/franceskoback/Single-Cell/data/regev_lab_cell_cycle_genes.txt /wynton/home/srivastava/franceskoback/Single-Cell/data/rds/noFilters_scoresAdded.rds 

echo Finished with Step1
