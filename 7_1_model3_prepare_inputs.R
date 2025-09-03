# model3_prepare_inputs.R
library(Seurat)
library(dplyr)

car_subset <- readRDS('4.2.1.2/car_combined.rds')
car_subset$region <- sub('.*_(AL|VA|VS)$', '\\1', car_subset$orig.ident)

subset_cells <- WhichCells(car_subset, expression = cell_types %in% c('Ventricular cardiomyocytes','Fibroblasts','Myeloid','Lymphoid'))
car_subset_filtered <- subset(car_subset, cells = subset_cells)

car_subset_filtered <- NormalizeData(car_subset_filtered)
car_subset_filtered <- FindVariableFeatures(car_subset_filtered)
car_subset_filtered <- ScaleData(car_subset_filtered)
car_subset_filtered <- RunPCA(car_subset_filtered, npcs = 30)

pc_mat <- Embeddings(car_subset_filtered, 'pca')[, 1:30]
meta_df <- car_subset_filtered@meta.data

pc_df <- as.data.frame(pc_mat)
colnames(pc_df) <- paste0('PC', 1:30)
pc_df$cell_id <- rownames(meta_df)
pc_df$region <- meta_df$region
pc_df$cell_types <- meta_df$cell_types

write.csv(pc_df, './model3_fourtypes_region.csv', row.names = FALSE)
