library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)

sct_murine_cl00 <- readRDS("/home/jing/resource/murine_blca/sct_murine_cl00_29th.rds")

corrected_UMI <- sct_murine_cl00[["SCT"]]$data
corrected_UMI <-t(corrected_UMI)

colnames_substr <- substr(colnames(corrected_UMI), 1, 10)

#input interested datasets
indices <- which(colnames_substr == "GSM5288674")
GSM5288674_corrected_UMI <- corrected_UMI[, indices]
GSM5288674_corrected_UMI<- t(GSM5288674_corrected_UMI)
writeMM(GSM5288674_corrected_UMI, "sct_full_UM/GSM5288674_corrected_UMI.mtx")
write.table(rownames(GSM5288674_corrected_UMI), file = "sct_full_UM/GSM5288674_corrected_UMI_cells.txt", col.names = FALSE, row.names = FALSE)
write.table(colnames(GSM5288674_corrected_UMI), file = "sct_full_UM/GSM5288674_corrected_UMI_genes.txt", col.names = FALSE, row.names = FALSE)

#Read data in scanpy 