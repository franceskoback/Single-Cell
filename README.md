# Single-Cell
Workflow for Single-Cell Analysis

To use this script:
1. git clone https://github.com/franceskoback/Single-Cell
2. Add your own inputs to the data and filtered_featured_bc folders as specified
3. Use the script as follows: 

Rscript Step1_Preprocessing.R /wynton/home/srivastava/franceskoback/SRA_Analysis/filtered_feature_bc_matrix /wynton/home/srivastava/franceskoback/SRA_Analysis/data/AggrSheet.csv /wynton/home/srivastava/franceskoback/SRA_Analysis/data/regev_lab_cell_cycle_genes.txt /wynton/home/srivastava/franceskoback/SRA_Analysis/data/noFilters_scoresAdded.rds

Or, if you want to submit as a job, run a bash file, such as qsub -m ea -M frances.koback@gladstone.ucsf.edu bash_script.sh 
An example bash file that you could use is attached in this repository-- edit it as you wish on your computer
