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


sct_murine_cl00 <- readRDS("/home/jing/resource/murine_blca/sct_murine_cl00_29th.rds")



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
GSM5288674 <- subset_seurat_by_sample(GSM5288674, sct_murine_cl00@meta.data, "GSM5288674")



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

#Data ready to be retrieved here 
#saveRDS(seurat.integrated,"sct_cl00.rsd")
DefaultAssay(seurat.integrated) <- "integrated"


#Clustering and visualisation

# Run the standard workflow for visualization and clustering
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:20)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:20, verbose = FALSE)
#resolution to be specified
seurat.integrated <- FindClusters(seurat.integrated, verbose = FALSE,resolution=0.1)

p1 <- DimPlot(seurat.integrated, reduction = "umap")
p2 <- DimPlot(seurat.integrated, reduction = "umap", split.by  = "Sample", 
              repel = TRUE)


#ggsave("UMAPreclustering_split.png",plot=p2,width=15,height=6,dpi=300)


#Raw RNA find markers
DefaultAssay(seurat.integrated) <- "RNA"

cl00.markers <- FindAllMarkers(seurat.integrated)

#
cl00.markers <- FindConservedMarkers(onc, ident.1 = "", grouping.var = "stim", verbose = FALSE)
head(nk.markers)
