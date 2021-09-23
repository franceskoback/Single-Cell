# Single-Cell
Workflow for Single-Cell Analysis

To use this script:
1. git clone https://github.com/franceskoback/Single-Cell
2. Add your own inputs to the data and filtered_featured_bc folders as specified
3. Use the script as follows: 

Rscript Step1_Preprocessing.R /wynton/home/srivastava/franceskoback/SRA_Analysis/filtered_feature_bc_matrix /wynton/home/srivastava/franceskoback/SRA_Analysis/data/AggrSheet.csv /wynton/home/srivastava/franceskoback/SRA_Analysis/data/regev_lab_cell_cycle_genes.txt /wynton/home/srivastava/franceskoback/SRA_Analysis/data/noFilters_scoresAdded.rds
Rscript Visualize_Boundaries.R /wynton/home/srivastava/franceskoback/SRA_Analysis/data/noFilters_scoresAdded.rds 1 10 13 33 2500 23000 1000 5000
Rscript Step2_Clustering_Filtering.R /wynton/home/srivastava/franceskoback/SRA_Analysis/data/noFilters_scoresAdded.rds 1 10 13 33 2500 23000 1000 5000 /wynton/home/srivastava/franceskoback/SRA_Analysis/data/rds/FilteredAndClustered_onlyVarGenes.rds /wynton/home/srivastava/franceskoback/SRA_Analysis/data/rds/FilteredAndReClustered_onlyVarGenes.rds /wynton/home/srivastava/franceskoback/SRA_Analysis/data/csv/FindAllMarkers_harmony_res25e-2_05-24-2021.csv

Or, if you want to submit as a job, run a bash file, such as qsub -m ea -M frances.koback@gladstone.ucsf.edu bash_script.sh 
An example bash file that you could use is attached in this repository-- edit it as you wish on your computer
