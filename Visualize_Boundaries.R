library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)
 set.seed(7)
parser = argparse::ArgumentParser(description="Script to Visualize Boundaries on Mitochondrial, Ribosomal, nCounts, and nFeature Plots")
parser$add_argument('unfiltered_rds', help='unfiltered rds to visualize metrics plots with filter lines')
parser$add_argument('mito_start', help='mito start')
parser$add_argument('mito_end', help='mito end')
parser$add_argument('ribo_start', help='ribo start')
parser$add_argument('ribo_end', help='ribo end')
parser$add_argument('nCount_start', help='nCount start')
parser$add_argument('nCount_end', help='nCount end')
parser$add_argument('nFeature_start', help='nFeature  start')
parser$add_argument('nFeature_end', help='nFeature end')

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


print("Done! Please progress to Notebook 2 with user inputs as specified depending on mito, ribo, nFeat, and nCount metrics")


