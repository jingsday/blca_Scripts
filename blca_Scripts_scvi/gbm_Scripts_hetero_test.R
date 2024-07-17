# Load the required packages
library(Seurat)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(stringr)

#Part I: Read data, merge, create Seurat object and preprocessing(if provided raw)

setwd('/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_mouse_GSE174182_RAW')

# get data location
dirs <- list.files(path = '/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_mouse_GSE174182_RAW', recursive = F, full.names = F)


names_list <- c()  # Create an empty list to store the names
for (x in dirs) {
  name <- unique(str_sub(x, end = 20))  # Extract the name
  names_list <- c(names_list, name)     # Append the name to the list
}
names_list <-unique(names_list)
names_list
names_list <- names_list[names_list != "GSM5288674_Sample-11"]

# Add "GSM5288674_Sample-11_" to the list
names_list <- c(names_list, "GSM5288674_Sample-11_")


for (name in names_list) {
  cts <- ReadMtx(mtx = paste0(name,'filtered_matrix.mtx.gz'),
                 features = paste0(name,'filtered_features.tsv.gz'),
                 cells = paste0(name,'filtered_barcodes.tsv.gz'))
  
  # create seurat objects
  assign(str_sub(name, end = 10), CreateSeuratObject(counts = cts))
}


ls()

merged_seurat <- merge(GSM5288668, y = c(GSM5288669,GSM5288670, GSM5288671, GSM5288672, GSM5288673, GSM5288674),
                       add.cell.ids = ls()[3:9],
                       project = 'BLCA')
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^Mt-')

merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 &
                                   nFeature_RNA < 8000 &
                                   mitoPercent < 25)


# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')

#obj.list <- SplitObject(merged_seurat, split.by = 'Sample')


obj <- NormalizeData(merged_seurat)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Sample"))

library(reticulate)
library(SeuratWrappers)
library(patchwork)
obj <- IntegrateLayers(object=obj,method=scVIIntegration,
                       new.reduction ='integrated.scvi',
                       conda_env='/home/jing/miniforge3/envs/scvi',verbose=FALSE)


###
merged_seurat_filtered <- RunPCA(merged_seurat_filtered, verbose = FALSE)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:30, verbose = FALSE)

merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:30, verbose = FALSE)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, verbose = FALSE)
DimPlot(merged_seurat_filtered)

merged_seurat_filtered@reductions$

options(future.globals.maxSize = 3e+09)
obj <- SCTransform(obj)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "integrated.dr")
obj <- FindClusters(obj, resolution = 2)


obj <-


obj <- subset(merged_seurat_filtered, nFeature_RNA > 1000)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Sample"))


obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p1 <- DimPlot(
  obj,
  reduction = "harmony",
  group.by = c("Sample", "harmony_clusters"),
  combine = FALSE, label.size = 2
)
p1
obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "scvi_clusters")



obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Sample')

common_genes <- Reduce(intersect, lapply(obj.list, rownames))
length(common_genes)

obj.list_aligned <- lapply(X = obj.list_aligned, FUN = function(x) SCTransform(x, return.only.var.genes = FALSE, assay = "RNA", new.assay.name = "SCT"))
obj.list_aligned <- lapply(obj.list, function(x) {
  x <- subset(x, features = common_genes)
  return(x)
})

obj.list_aligned <- lapply(X = obj.list_aligned, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = obj.list_aligned, nfeatures = 3000)
obj.list_aligned <- PrepSCTIntegration(object.list = obj.list_aligned, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list_aligned, normalization.method = "SCT",
                                  anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")








