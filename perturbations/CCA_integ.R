# Load the required packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)



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

# QC & filtering -----------------------

# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^Mt-')

#Murine cells were filtered to retain higher quality cells (>200 &<8000 uniquely identified genes
#<25% of reads mapped to mitochondrial genes),
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 &
                                   nFeature_RNA < 8000 &
                                   mitoPercent < 25)

merged_seurat_filtered
#Visualize QC metrics as a violin plot
#VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#
FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
#
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)

#

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_seurat_filtered), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(merged_seurat_filtered)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2


#Scale
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)

#PCA
merged_seurat_filtered <- RunPCA(merged_seurat_filtered, features = VariableFeatures(object = merged_seurat_filtered))

print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(merged_seurat_filtered, dims = 1:2, reduction = "pca")

#DimPlot(merged_seurat_filtered, reduction = "pca",group.by = "Sample",alpha = 0.5,
#        cols=c('red','green','blue','yellow','brown','orange')) 
#DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)

#DimHeatmap(merged_seurat_filtered, dims = 13:15, cells = 500, balanced = TRUE)

#Elbowplot
#ElbowPlot(merged_seurat_filtered)

#Clustering
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:30)

merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered,resolution=0.1)
head(Idents(merged_seurat_filtered), 5)

#UMAP
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:30)
#p1 <- DimPlot(merged_seurat_filtered, reduction = "umap")
#p2 <- DimPlot(merged_seurat_filtered, reduction = "umap",group.by="Sample"
#        ,alpha = 0.5,
#        cols=c('red','green','blue','yellow','brown','orange'))

#p1+p2

# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')

# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Sample')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)


# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:30)



p3 <- DimPlot(seurat.integrated, reduction = 'umap')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Sample',
              cols=c('red','green','blue','yellow','brown','orange','purple'))
combined<-p3+p4

markers <- FindAllMarkers(seurat.integrated, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)


p1 <-VlnPlot(seurat.integrated, features = c("Cdh1", "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b", "Ivl"),pt.size = 0.2, ncol = 4)

p2 <-VlnPlot(seurat.integrated, features = c("Cd8a", "Cd8b1", "Cd4", "Cd3e", "Cd3d", "Cd3g", "Il2ra"),pt.size = 0.2, ncol = 4)

# Save the plot
ggsave(".png", plot = p1, width = 10, height = 6, dpi = 300)

