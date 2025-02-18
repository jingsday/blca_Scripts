# Load the required packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(sctransform)



setwd('/home/jing/resource/murine_blca')
outdir <- '/home/jing/resource/murine_blca/scripts/'

# get data location
dirs <- list.files(path = '/home/jing/resource/murine_blca/GSE174182_RAW', recursive = F, full.names = F)


names_list <- c()  # Create an empty list to store the names
for (x in dirs) {
  name <- unique(str_sub(x, end = 20))  # Extract the name
  names_list <- c(names_list, name)     # Append the name to the list
}
names_list <-unique(names_list)

names_list <- names_list[names_list != "GSM5288674_Sample-11"]

# Add "GSM5288674_Sample-11_" to the list
names_list <- c(names_list, "GSM5288674_Sample-11_")


for (name in names_list) {
  cts <- ReadMtx(mtx = paste0('GSE174182_RAW/',name,'filtered_matrix.mtx.gz'),
                 features = paste0('GSE174182_RAW/',name,'filtered_features.tsv.gz'),
                 cells = paste0('GSE174182_RAW/',name,'filtered_barcodes.tsv.gz'))
  
  # create seurat objects
  assign(str_sub(name, end = 10), CreateSeuratObject(counts = cts))
}







ls()

merged_seurat <- merge(GSM5288668, y = c(GSM5288669,GSM5288670, GSM5288671, GSM5288672, GSM5288673, GSM5288674),
                       add.cell.ids = ls()[3:9],
                       project = 'BLCA')

# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')


obj.list <- SplitObject(merged_seurat, split.by = 'Sample')


obj.list <- lapply(X = obj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                         anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

#saveRDS(seurat.integrated, file = "murine_intergration_sct.rds")


DefaultAssay(seurat.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)


head(seurat.integrated@meta.data@seurat)
p1 <- DimPlot(seurat.integrated, reduction = "umap")
p2 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Sample", 
              repel = TRUE,cols=c('red','green','blue','yellow','brown','orange','purple'))

p1 + p2
combined <- p1+p2

ggsave("UMAPsct29th.png",plot=combined,width=10,height=6,dpi=300)



p1 <- VlnPlot(seurat.integrated, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"),
              pt.size = 0.2, ncol = 4)
ggsave("uromarkers_sct.png",plot=p1,width=10,height=12,dpi=300)

p1 <- FeaturePlot(seurat.integrated, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"), pt.size = 0.2,
                  ncol = 3)
ggsave("urofeatures_sct.png",plot=p1,width=10,height=12,dpi=300)


#Barcodes and samples 
cluster_0_cells <- which(seurat.integrated@meta.data$seurat_clusters == 0)

# Subset seurat.integrated to include only cells belonging to cluster 0
seurat_integrated_cluster_0 <- seurat.integrated[cluster_0_cells, ]





#saveRDS(seurat.integrated, file = "murine_intergration_sct.rds")


sct_murine <- readRDS("/home/jing/resource/murine_blca/murine_intergration_sct.rds")


p1 <- VlnPlot(sct_murine, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"),
        pt.size = 0.2, ncol = 4)
ggsave("uromarkers_sct.png",plot=p1,width=10,height=12,dpi=300)

p1 <- FeaturePlot(sct_murine, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"), pt.size = 0.2,
            ncol = 3)
ggsave("urofeatures_sct.png",plot=p1,width=10,height=12,dpi=300)


DefaultAssay(seurat.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
#sct_murine <- FindClusters(sct_murine, resolution = 0.5)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)

##
seurat.integrated <- subset(sct_murine, subset = seurat_clusters == 0)
#saveRDS(sct_murine_cl00, file = "sct_murine_cl00_29th.rds")

sct_murine_cl00 <- subset(sct_murine, subset = Sample == GSM5288672)
Idents(object = sct_cl00) <- sct_cl00@meta.data$Sample

#subsets
onc_id <- c("GSM5288674","GSM5288668","GSM5288669","GSM5288670","GSM5288671")
onc <- subset(sct_cl00, idents = onc_id)

onc_test <- NormalizeData(onc, normalization.method = "LogNormalize", scale.factor = 10000)

onc.markers <- FindAllMarkers(onc_test, only.pos = TRUE)
onc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
onc.markers <- FindAllMarkers(onc,verbose = FALSE)
