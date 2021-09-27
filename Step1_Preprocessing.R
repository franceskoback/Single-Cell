#!/usr/bin/env Rscript
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(Seurat)
library(optparse)
library(ggplot2)

### human to mouse gene conversion function (for cc genes)
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

parser = argparse::ArgumentParser(description="Script to preprocess single cell data")
parser$add_argument('tenX', help='filtered_feature_bc_matrix')
#parser$add_argument('project_name', help='project name')
parser$add_argument('aggregation_csv',help='sample_id, gem_group, condition, cellID')
parser$add_argument('genes_txt',help='list of genes of interest')
parser$add_argument('output_file',help='output noFilters_scoresAdded.rds')


args = parser$parse_args()

cat('Loading 10X library...\n')
data.tenXdata <- Read10X(data.dir =args$tenX)
#tenXdata <- CreateSeuratObject(counts = data.tenXdata , project = args$project_name, min.cells = 3, min.features = 200)
tenXdata <- CreateSeuratObject(counts = data.tenXdata, min.cells = 3, min.features = 200)


rm(data.tenXdata)

AggrSheet <- read.csv(args$aggregation_csv)

### same start as the old vector-based approach
cellID <- as.numeric(gsub(".*-","", (colnames(x = tenXdata)))) ### this variable name becomes the colname in the DF
names(cellID) <- colnames(x = tenXdata)

### coerce named list to dataframe
metadata2add <- data.frame(cellID)

### fill out metadata DF based on annotated aggr.csv
## pull rownames as column - clunky but easy
## TODO: merge(), but keep rownames and eliminate this rownames>colnames>rownames business
metadata2add$rownames_for_metadata <- rownames(metadata2add)
metadata2add <- merge(metadata2add, AggrSheet, by="cellID", all.x = TRUE, no.dups = FALSE, )
rownames(metadata2add) <- metadata2add$rownames_for_metadata
## drop columns we don't want in the metadata
metadata2add$cellID <- NULL
metadata2add$rownames_for_metadata <- NULL
head(metadata2add)

### AddMetaData
tenXdata <- AddMetaData(tenXdata, metadata = metadata2add)
head(tenXdata@meta.data)

levels(as.factor(tenXdata@meta.data$condition))

# calc % mt
tenXdata[["percent.mt"]] <- PercentageFeatureSet(tenXdata, pattern = "^mt-")
# calc % ribo
tenXdata[["percent.ribo"]] <- PercentageFeatureSet(tenXdata, pattern = "^Rp[sl]")
# cell cycle scoring
### Assign cell cycle scores
cc.genes <- readLines(args$genes_txt)
## convert Human to Mouse
cc.genes <- convertHumanGeneList(cc.genes)
## segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:45]
g2m.genes <- cc.genes[46:100]
## assign cell cycle scores
tenXdata <- CellCycleScoring(tenXdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(tenXdata@meta.data)

cat('Saving rds...\n')
saveRDS(tenXdata, file = args$output_file)
#saveRDS(tenXdata, file = "/wynton/home/srivastava/franceskoback/SRA_Analysis/data/noFilters_scoresAdded.rds")

cat('rds saved! Generating plots for user input in next step ...\n')
qc.metrics <- tenXdata[[c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")]]


### plot w/ ggplot
# % mito
qc.metrics %>%
  ggplot(aes(percent.mt)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.mt")
ggsave("../plots/initial_clustering/mito_init.png,", width = 12, height = 12, device = "png")


# % ribo
qc.metrics %>%
  ggplot(aes(percent.ribo)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.ribo")
ggsave("../plots/initial_clustering/ribo_init.png,", width = 12, height = 12, device = "png")


# nCount
qc.metrics %>%
  ggplot(aes(nCount_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nCount_RNA")
ggsave("../plots/initial_clustering/nCount_init.png,", width = 12, height = 12, device = "png")

# nFeature
qc.metrics %>%
  ggplot(aes(nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nFeature_RNA")
ggsave("../plots/initial_clustering/nFeature_init.png,", width = 12, height = 12, device = "png")


print("Done! Please progress to Notebook 2 with user inputs as specified depending on mito, ribo, nFeat, and nCount metrics")


