# Load the required packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)


setwd('/home/jing/resource/murine_blca')

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
#murine_clean <- readRDS("/home/jing/resource/murine_blca/murine_intergration.rds")
sct_murine_cl00 <- readRDS("/home/jing/resource/murine_blca/sct_murine_cl00.rds")




subset_seurat_by_sample <- function(seurat_obj, metadata, sample_id_value) {
  # Extract metadata for the specified sample ID
  sample_metadata <- metadata[metadata[["Sample"]] == sample_id_value, ]
  
  # Extract barcodes from the metadata
  common_row_names <- sample_metadata$Barcode
  
  # Extract barcodes from the Seurat object
  barcodes <- colnames(seurat_obj)
  
  # Subset the Seurat object based on common row names
  subset_seurat <- seurat_obj[, barcodes %in% common_row_names]
  
  return(subset_seurat)
}

# Usage example
# Usage example
GSM5288674 <- subset_seurat_by_sample(GSM5288674, sct_murine_cl00@meta.data, "GSM5288674")
GSM5288673 <- subset_seurat_by_sample(GSM5288673, sct_murine_cl00@meta.data, "GSM5288673")
GSM5288672 <- subset_seurat_by_sample(GSM5288672, sct_murine_cl00@meta.data, "GSM5288672")
GSM5288671 <- subset_seurat_by_sample(GSM5288671, sct_murine_cl00@meta.data, "GSM5288671")
GSM5288670 <- subset_seurat_by_sample(GSM5288670, sct_murine_cl00@meta.data, "GSM5288670")
GSM5288669 <- subset_seurat_by_sample(GSM5288669, sct_murine_cl00@meta.data, "GSM5288669")
GSM5288668 <- subset_seurat_by_sample(GSM5288668, sct_murine_cl00@meta.data, "GSM5288668")



# Create the violin plot using the "RNA" assay
p1 <- VlnPlot(murine_clean, features = c("Cdh1", "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b", "Ivl"), assay = "RNA")

# Save the plot
ggsave("genes_violin_plot.png", plot = p1, width = 10, height = 6, dpi = 300)


merged_seurat <- merge(GSM5288668, y = c(GSM5288669,GSM5288670, GSM5288671, GSM5288672, GSM5288673, GSM5288674),
                       add.cell.ids = ls()[4:10],
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
#Visualize QC metrics as a violin #plot
#Vln#plot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#
#p1 <-FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#ggsave("feature#plot.png",#plot=#p1,width=10,height=6,dpi=300)
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
#
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)


# #plot variable features with and without labels
#plot1 <- VariableFeature#plot(merged_seurat_filtered)
#plot2 <- LabelPoints(#plot = #plot1, points = to#p10, repel = TRUE)
combined <-#plot1 + #plot2
###ggsave("HVG_0004.png",#plot=combined,width=10,height=6,dpi=300)


merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)

#PCA
merged_seurat_filtered <- RunPCA(merged_seurat_filtered, features = VariableFeatures(object = merged_seurat_filtered))
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)

#p1 <- VizDimLoadings(merged_seurat_filtered, dims = 1:2, reduction = "pca")
##ggsave("pca2.png",#plot=#p1,width=10,height=6,dpi=300)

#p1 <-Dim#plot(merged_seurat_filtered, reduction = "pca",group.by = "Sample",alpha = 0.5,
#        cols=c('red','green','blue','yellow','brown','orange',"purple")) 
#p2 <- DimHeatmap(merged_seurat_filtered, dims = 1:5, cells = 500, balanced = TRUE)
#ggsave("pcaprojection.png",#plot=#p1,width=10,height=6,dpi=300)
#ggsave("pca5.png",#plot=p2,width=10,height=6,dpi=300)


#Elbow#plot
#p1 <-Elbow#plot(merged_seurat_filtered)
#ggsave("elbow.png",#plot=#p1,width=10,height=6,dpi=300)


#Clustering
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:30)

merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered,resolution=0.5)
head(Idents(merged_seurat_filtered), 5)

#p1 <- Dim#plot(merged_seurat_filtered, reduction = "umap")
#p2 <- Dim#plot(merged_seurat_filtered, reduction = "umap",group.by="Sample"
              ,alpha = 0.5,
              cols=c('red','green','blue','yellow','brown','orange',"purple"))

combined <- #p1+p2

#ggsave("UMAP0004.png",#plot=combined,width=10,height=6,dpi=300)


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

seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.1)

saveRDS(seurat.integrated, file = "intergdims30_clusters0004.rds")

#Vis

p1 <- Dimplot(seurat.integrated, reduction = "umap")
p2 <- Dimplot(seurat.integrated, reduction = "umap",group.by="Sample"
              ,alpha = 0.5,
              cols=c('red','green','blue','yellow','brown','orange',"purple"))

combined <- p1+p2

ggsave("UMAPrepo1.png",plot=combined,width=10,height=6,dpi=300)

uro.markers <- FindAllMarkers(seurat.integrated, only.pos = TRUE)
uro.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)

uro.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> to#p10
#p1 <- DoHeatmap(seurat.integrated, features = to#p10$gene) + NoLegend()


#ggsave("uro_markers.png",#plot=#p1,width=10,height=6,dpi=300)
