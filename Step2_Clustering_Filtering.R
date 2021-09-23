library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
library(Seurat)
library(optparse)
set.seed(7)

parser = argparse::ArgumentParser(description="Script to QC Filter and Cluster")
parser$add_argument('input_file', help='input noFilters_scoresAdded.RDS')
parser$add_argument('clustered_rds', help='output file name FilteredAndClustered_onlyVarGenes')
parser$add_argument('reclustered_rds', help='output file name FilteredAndReClustered_onlyVarGenes')
parser$add_argument('csv',help='output file FindAllMarkers_harmony_res25e-2_05-24-2021.csv')


args = parser$parse_args()

clustered_rds <- args$clustered_rds
reclustered_rds <- args$reclustered_rds
csv <- args$csv

#print(paste("input_file: ", clustered_rds))
print(paste("clustered_rds: ", clustered_rds))
print(paste("reclustered_rds: ", reclustered_rds))
print(paste("csv: ",csv ))

tenXdata <- readRDS(file = args$input_file)
head(tenXdata@meta.data)
qc.metrics <- tenXdata[[c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")]]


### plot w/ ggplot
# % mito
qc.metrics %>%
  ggplot(aes(percent.mt)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.mt") +
  geom_vline(xintercept = c(1, 10)) ### you can play with these values to get an idea of where the filters should be set...

# % ribo
qc.metrics %>%
  ggplot(aes(percent.ribo)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.ribo") +
  geom_vline(xintercept = c(13, 33))

# nCount
qc.metrics %>%
  ggplot(aes(nCount_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nCount_RNA") +
  geom_vline(xintercept = c(2500, 23000))

# nFeature
qc.metrics %>%
  ggplot(aes(nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nFeature_RNA") +
  geom_vline(xintercept = c(1000, 5000))

### filter by % mito and rido + nFeature (counts can be adjusted  momentarily, if needed)
## store the thresholds as variables for easy use later...
mito_cutoffs <- c(1,10)
ribo_cutoffs <- c(13, 33)
nFeat_cutoffs <- c(1000, 5000)
nCount_cutoffs <- c(2500, 23000)

qc.metrics_subset <- qc.metrics[qc.metrics$percent.mt > min(mito_cutoffs) & qc.metrics$percent.mt < max(mito_cutoffs) &
                                  qc.metrics$percent.ribo > min(ribo_cutoffs) & qc.metrics$percent.ribo <  max(ribo_cutoffs) &
                                  qc.metrics$nFeature_RNA > min(nFeat_cutoffs) & qc.metrics$nFeature_RNA <  max(nFeat_cutoffs) &
                                  qc.metrics$nCount_RNA > min(nCount_cutoffs) & qc.metrics$nCount_RNA < max(nCount_cutoffs),]

print(paste("# of cells before filters:", nrow(qc.metrics)))
print(paste("# of cells after filters:", nrow(qc.metrics_subset)))

# % mito
qc.metrics_subset %>%
  ggplot(aes(percent.mt)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.mt") +
  geom_vline(xintercept = c(min(mito_cutoffs), max(mito_cutoffs))) ### you can play with these values to get an idea of where the filters should be set...

# % ribo
qc.metrics_subset %>%
  ggplot(aes(percent.ribo)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.ribo") +
  geom_vline(xintercept = c(min(ribo_cutoffs), max(ribo_cutoffs)))

# nCount
qc.metrics_subset %>%
  ggplot(aes(nCount_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nCount_RNA") +
  geom_vline(xintercept = c(min(nCount_cutoffs), max(nCount_cutoffs)))

# nFeature
qc.metrics_subset %>%
  ggplot(aes(nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nFeature_RNA") +
  geom_vline(xintercept = c(min(nFeat_cutoffs), max(nFeat_cutoffs)))

### does the subset approach look familiar?
cached_data  <- subset(tenXdata,
                       subset = percent.mt > min(mito_cutoffs) & percent.mt < max(mito_cutoffs) &
                         percent.ribo > min(ribo_cutoffs) & percent.ribo <  max(ribo_cutoffs) &
                         nFeature_RNA > min(nFeat_cutoffs) & nFeature_RNA <  max(nFeat_cutoffs) &
                         nCount_RNA > min(nCount_cutoffs) & nCount_RNA < max(nCount_cutoffs))
Idents(tenXdata) <- "gem.group"
### violinplots
VlnPlot(tenXdata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0.1)
VlnPlot(tenXdata, features = c("percent.mt", "percent.ribo"), ncol = 2, pt.size=0.1)
### featurescatters
FeatureScatter(tenXdata, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.1, group.by = "gem.group")
FeatureScatter(tenXdata, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1, group.by = "gem.group")
FeatureScatter(tenXdata, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size=0.1, group.by = "gem.group")
FeatureScatter(tenXdata, feature1 = "nFeature_RNA", feature2 = "percent.ribo", pt.size=0.1, group.by = "gem.group")
FeatureScatter(tenXdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1, group.by = "gem.group")

# plot filtered data
# violinplots
VlnPlot(cached_data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0.1)
VlnPlot(cached_data, features = c("percent.mt", "percent.ribo"), ncol = 2, pt.size=0.1)
### featurescatters
FeatureScatter(cached_data, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.1, group.by = "gem.group")
FeatureScatter(cached_data, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1, group.by = "gem.group")
FeatureScatter(cached_data, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size=0.1, group.by = "gem.group")
FeatureScatter(cached_data, feature1 = "nFeature_RNA", feature2 = "percent.ribo", pt.size=0.1, group.by = "gem.group")
FeatureScatter(cached_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1, group.by = "gem.group")

tenXdata <- cached_data
#rm(cached_data)

#gdata::keep(tenXdata, sure = TRUE)
head(tenXdata@meta.data)

### SCTransform
tenXdata <- SCTransform(tenXdata, assay = "RNA", new.assay.name = "SCT",
                        variable.features.n = 4000,
                        variable.features.rv.th = 1.3,
                        vars.to.regress = c("percent.ribo", "percent.mt", "G2M.Score", "S.Score"),
                        return.only.var.genes = TRUE,conserve.memory= TRUE)

### PCA
tenXdata <- RunPCA(tenXdata, verbose = FALSE, npcs = 100)
### harmony
tenXdata <- harmony::RunHarmony(tenXdata, group.by.vars = "gem.group", assay.use="SCT")
### elbowplot, just to check
ElbowPlot(tenXdata, ndims = 50)

### set dims
dims <- 1:35 ### relatively high # of dims, UMAP can take it.
### NOTE: n.components = 3l is for 3d plotting *cool_emoji*)
### umap, neighbors, and clusters
tenXdata <- RunUMAP(tenXdata, reduction = "harmony", dims = dims) # , n.components = 3L)
tenXdata <- FindNeighbors(tenXdata, reduction = "harmony", dims = dims)
tenXdata <- FindClusters(tenXdata, resolution = 0.25)

DimPlot(tenXdata, split.by = "gem.group")
DimPlot(tenXdata)

#output_file1 = args$output_vargenes
cat('saving clustered rds ...\n')
#saveRDS(tenXdata,file=args$clustered_rds)
ls()
saveRDS(tenXdata,file= clustered_rds)
#saveRDS(tenXdata,file ="/wynton/home/srivastava/franceskoback/SRA_Analysis/data/FilteredAndClustered_onlyVarGenes.rds")

DimPlot(tenXdata, label = TRUE)
VlnPlot(tenXdata, features = "nFeature_RNA")
VlnPlot(tenXdata, features = "nCount_RNA")
# ggsave("../results/01_initial_clustering/01_FirstClustering_nFeature_05-21-2021.png", width = 12, height = 24, device = "png")
DimPlot(tenXdata, label = TRUE)
# ggsave("../results/01_initial_clustering/01_FirstClustering_Dimplot_05-21-2021.png,", width = 12, height = 12, device = "png")

clusters2keep <- 0:11
clusters2keep <- clusters2keep[clusters2keep != 9 &
                                 clusters2keep != 10 &
                                 clusters2keep != 11]
tenXdata <- subset(tenXdata, idents = clusters2keep)

### SCTransform - same settings as before
tenXdata <- SCTransform(tenXdata, assay = "RNA", new.assay.name = "SCT",
                        variable.features.n = 4000,
                        variable.features.rv.th = 1.3,
                        vars.to.regress = c("percent.ribo", "percent.mt", "G2M.Score", "S.Score"),
                        return.only.var.genes = TRUE,conserve.memory=TRUE)

### PCA
tenXdata <- RunPCA(tenXdata, verbose = FALSE, npcs = 100)
### harmony
tenXdata <- harmony::RunHarmony(tenXdata, group.by.vars = "gem.group", assay.use="SCT")
### elbowplot, just to check
ElbowPlot(tenXdata, ndims = 50)

### set dims
dims <- 1:35 ### relatively high # of dims, UMAP can take it.
### NOTE: n.components = 3l is for 3d plotting *cool_emoji*)
### umap, neighbors, and clusters
tenXdata <- RunUMAP(tenXdata, reduction = "harmony", dims = dims) # , n.components = 3L)
tenXdata <- FindNeighbors(tenXdata, reduction = "harmony", dims = dims)
tenXdata <- FindClusters(tenXdata, resolution = 0.25)

DimPlot(tenXdata, split.by = "gem.group")
DimPlot(tenXdata)

VlnPlot(tenXdata, features = "nFeature_RNA")
VlnPlot(tenXdata, features = "nCount_RNA")
#ggsave("../results/01_initial_clustering/02_SecondClustering_dropLowQual_nFeature_05-24-2021.png", width = 12, height = 24, device = "png")
DimPlot(tenXdata, label = TRUE)
#ggsave("../results/01_initial_clustering/02_SecondClustering_dropLowQual_Dimplot_05-24-2021.png,", width = 12, height = 12, device = "png")

### save RDS!
cat('saving reclustered rds ...\n')
#saveRDS(tenXdata,file=args$reclustered_rds)
saveRDS(tenXdata,file= reclustered_rds)
#saveRDS(tenXdata, file = "../data/rds/03_SCTandHarmonyObject_FilteredAndReclustered_onlyVarGenes_ps49v51_05-24-2021.RDS")
### find cluster markers
## parallelize
future::plan("multisession", workers = 6) ### 12-core // 128gb RAM
options(future.globals.maxSize = 3000 * 1024^2) ### global var size increase - 3gb limit
future::plan()

FoundMarkers <- FindAllMarkers(tenXdata, random.seed = 7)
future::plan("sequential") ### 12-core // 128gb RAM
#output_file2 = args$output_csv
#write.csv(FoundMarkers,output_file2)
cat('saving csv ...\n')
#write.csv(FoundMarkers, file = args$csv)
write.csv(FoundMarkers, file = csv)
#write.csv(FoundMarkers,"/wynton/home/srivastava/franceskoback/SRA_Analysis/data/findAllMarkers_harmony_res25e-2_05-24-2021.csv")



