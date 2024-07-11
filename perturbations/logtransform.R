# Load the required packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)

sct_murine_cl00 <- readRDS("/home/jing/resource/murine_blca/sct_murine_cl00_29th.rds")


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



GSM5288668 <- NormalizeData(GSM5288668, normalization.method = "LogNormalize", scale.factor = 10000)

GSM5288674 <- NormalizeData(GSM5288674, normalization.method = "LogNormalize", scale.factor = 10000)
healthy <- as.data.frame(log1p(GSM5288674$RNA,base=2))


# Extract sample barcodes from metadata
sample_barcodes <- sct_murine_cl00@meta.data$sample_barcodes

# Convert to data frame
sample_barcodes_df <- as.data.frame(sample_barcodes)

# Save as a .csv file
write.csv(as.data.frame(sct_murine_cl00@meta.data[c("Barcodes","Sample")]), file = "sample_barcodes.csv", row.names = FALSE)


