library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)
 set.seed(7)
parser = argparse::ArgumentParser(description="Script to Visualize Boundaries on Mitochondrial, Ribosomal, nCounts, and nFeature Plots")
parser$add_argument('unfiltered_rds', help='unfiltered rds to visualize metrics plots with filter lines')
parser$add_argument('--mito_start', help='mito start')
parser$add_argument('--mito_end', help='mito end')
parser$add_argument('--ribo_start', help='ribo start')
parser$add_argument('--ribo_end', help='ribo end')
parser$add_argument('--nCount_start', help='nCount start')
parser$add_argument('--nCount_end', help='nCount end')
parser$add_argument('--nFeature_start', help='nFeature  start')
parser$add_argument('--nFeature_end', help='nFeature end')

args = parser$parse_args()
cat('Loading unfiltered rds...\n')
tenXdata <- readRDS(file =args$unfiltered_rds)
qc.metrics <- tenXdata[[c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")]]

mito_start<-args$mito_start
print(as.numeric(mito_start))




### plot w/ ggplot
# % mito
qc.metrics %>%
  ggplot(aes(percent.mt)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.mt with lines at") +
  geom_vline(xintercept = c(as.numeric(args$mito_start),as.numeric(args$mito_end)))
ggsave("plots/checking_filters/mito.png", width = 12, height = 12, device = "png")


# % ribo
qc.metrics %>%
  ggplot(aes(percent.ribo)) +
  geom_histogram(binwidth = 0.4, fill="yellow", colour="black") +
  ggtitle("percent.ribo")+
geom_vline(xintercept = c(as.numeric(args$ribo_start),as.numeric(args$ribo_end)))

ggsave("plots/checking_filters/ribo.png", width = 12, height = 12, device = "png")


# nCount
qc.metrics %>%
  ggplot(aes(nCount_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nCount_RNA")+
geom_vline(xintercept = c(as.numeric(args$nCount_start),as.numeric(args$nCount_end)))
ggsave("plots/checking_filters/nCount.png", width = 12, height = 12, device = "png")

# nFeature
qc.metrics %>%
  ggplot(aes(nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nFeature_RNA")+
geom_vline(xintercept = c(as.numeric(args$nFeature_start),as.numeric(args$nFeature_end)))
ggsave("plots/checking_filters/nFeature.png", width = 12, height = 12, device = "png")

mito_cutoffs <- c(args$mito_start,args$mito_end)
ribo_cutoffs <- c(args$ribo_start,args$ribo_end)
nFeat_cutoffs <- c(args$nFeature_start,args$nFeature_end)
nCount_cutoffs <- c(args$nCount_start,args$nCount_end)

qc.metrics_subset <- qc.metrics[qc.metrics$percent.mt > min(mito_cutoffs) & qc.metrics$percent.mt < max(mito_cutoffs) &
                                  qc.metrics$percent.ribo > min(ribo_cutoffs) & qc.metrics$percent.ribo <  max(ribo_cutoffs) &
                                  qc.metrics$nFeature_RNA > min(nFeat_cutoffs) & qc.metrics$nFeature_RNA <  max(nFeat_cutoffs) &
                                  qc.metrics$nCount_RNA > min(nCount_cutoffs) & qc.metrics$nCount_RNA < max(nCount_cutoffs),]
qc.metrics_subset %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.4, fill="blue", colour="black") +
  ggtitle("percent.mt") +
  geom_vline(xintercept = c(min(mito_cutoffs), max(mito_cutoffs))) ### you can play with these values to get an idea of where the filters should be set...
ggsave("plots/checking_filters/qcmetrics_plots/percent_mt.png", width = 12, height = 12, device = "png")

# % ribo
qc.metrics_subset %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.4, fill="pink", colour="black") +
  ggtitle("percent.ribo") +
  geom_vline(xintercept = c(min(ribo_cutoffs), max(ribo_cutoffs)))
ggsave("plots/checking_filters/qcmetrics_plots/percent_ribo.png", width = 12, height = 12, device = "png")

# nCount
qc.metrics_subset %>%
  ggplot(aes(nCount_RNA)) + 
  geom_histogram(binwidth = 100, fill="yellow", colour="black") +
  ggtitle("nCount_RNA") +
  geom_vline(xintercept = c(min(nCount_cutoffs), max(nCount_cutoffs)))
ggsave("plots/checking_filters/qcmetrics_plots/nCount.png", width = 12, height = 12, device = "png")

# nFeature
qc.metrics_subset %>%
  ggplot(aes(nFeature_RNA)) + 
  geom_histogram(binwidth = 100, fill="cyan", colour="black") +
  ggtitle("nFeature_RNA") +
  geom_vline(xintercept = c(min(nFeat_cutoffs), max(nFeat_cutoffs)))
print(paste("# of cells before filters:", nrow(qc.metrics)))
print(paste("# of cells after filters:", nrow(qc.metrics_subset)))
ggsave("plots/checking_filters/qcmetrics_plots/nFeature.png", width = 12, height = 12, device = "png")

cached_data  <- subset(tenXdata, 
                       subset = percent.mt > min(mito_cutoffs) & percent.mt < max(mito_cutoffs) &
                                  percent.ribo > min(ribo_cutoffs) & percent.ribo <  max(ribo_cutoffs) &
                                  nFeature_RNA > min(nFeat_cutoffs) & nFeature_RNA <  max(nFeat_cutoffs) &
                                  nCount_RNA > min(nCount_cutoffs) & nCount_RNA < max(nCount_cutoffs))

cached_data  <- subset(tenXdata, 
                       subset = percent.mt > min(mito_cutoffs) & percent.mt < max(mito_cutoffs) &
                                  percent.ribo > min(ribo_cutoffs) & percent.ribo <  max(ribo_cutoffs) &
                                  nFeature_RNA > min(nFeat_cutoffs) & nFeature_RNA <  max(nFeat_cutoffs) &
                                  nCount_RNA > min(nCount_cutoffs) & nCount_RNA < max(nCount_cutoffs))
Idents(tenXdata) <- "gem.group"
### violinplots
VlnPlot(tenXdata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0.1)
ggsave("plots/checking_filters/Prefiltering/nFeature_nCount_violin.png", width = 12, height = 12, device = "png")
VlnPlot(tenXdata, features = c("percent.mt", "percent.ribo"), ncol = 2, pt.size=0.1)
ggsave("plots/checking_filters/Prefiltering/mito_ribo_violin.png", width = 12, height = 12, device = "png")

### featurescatters
FeatureScatter(tenXdata, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.1, group.by = "gem.group")
ggsave("plots/checking_filters/feature_scatters/Prefiltering/nCount.png", width = 12, height = 12, device = "png")
FeatureScatter(tenXdata, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1, group.by = "gem.group")
ggsave("plots/checking_filters/feature_scatters/Prefiltering/nCount.png", width = 12, height = 12, device = "png")
FeatureScatter(tenXdata, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size=0.1, group.by = "gem.group")
ggsave("plots/checking_filters/feature_scatters/Prefiltering/nCount.png", width = 12, height = 12, device = "png")
FeatureScatter(tenXdata, feature1 = "nFeature_RNA", feature2 = "percent.ribo", pt.size=0.1, group.by = "gem.group")
ggsave("plots/checking_filters/feature_scatters/Prefiltering/nCount.png", width = 12, height = 12, device = "png")
FeatureScatter(tenXdata, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", pt.size=0.1, group.by = "gem.group")
ggsave("plots/checking_filters/feature_scatters/Prefiltering/nCount.png", width = 12, height = 12, device = "png")
FeatureScatter(tenXdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1, group.by = "gem.group")
ggsave("plots/checking_filters/feature_scatters/Prefiltering/nCount.png", width = 12, height = 12, device = "png")



print("Check these plots and adjust boundaries as needed until you are satisfied with the results, then progress to Step2_Clustering_Filtering.R")


