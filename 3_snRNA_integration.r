library(Seurat)
library(dplyr)

read_patient_parts <- function(pid, parts = c("AL","VA","VS")) {
  objs <- lapply(parts, function(p) {
    fn <- file.path(paste0("H", pid, "_", p), paste0("H", pid, "_", p, "_filtered.rds"))
    readRDS(fn)
  })
  obj <- Reduce(function(a, b) merge(a, y = b), objs)
  obj@meta.data <- obj@meta.data[, 1:5]
  obj@meta.data$patient <- paste0("H", pid)
  obj
}

patients <- lapply(1:8, read_patient_parts)
names(patients) <- paste0("H", 1:8)
car_raw <- Reduce(function(a, b) merge(a, y = b), patients)

car_raw <- NormalizeData(car_raw, normalization.method = "LogNormalize", scale.factor = 10000)
car_raw <- FindVariableFeatures(car_raw, selection.method = "vst")
car_raw <- ScaleData(car_raw, vars.to.regress = c("nCount_RNA", "percent.mt"))

car.list <- SplitObject(car_raw, split.by = "patient")
car.list <- lapply(car.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x
})
features <- SelectIntegrationFeatures(object.list = car.list)
anchors  <- FindIntegrationAnchors(object.list = car.list, anchor.features = features)
car.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(car.combined) <- "integrated"
car.combined <- ScaleData(car.combined, verbose = FALSE)
car.combined <- RunPCA(car.combined, npcs = 50, verbose = FALSE)
car.combined <- FindNeighbors(car.combined, reduction = "pca", dims = 1:50)
car.combined <- FindClusters(car.combined, resolution = 0.5)
car.combined <- RunUMAP(car.combined, reduction = "pca", dims = 1:50)

new.cluster.ids <- c(
  "0"="Ventricular cardiomyocytes", "1"="Endothelial", "2"="Endothelial",
  "3"="Fibroblasts", "4"="Pericytes", "5"="Myeloid",
  "6"="Ventricular cardiomyocytes", "7"="Fibroblasts", "8"="Endothelial",
  "9"="Endothelial", "10"="Lymphoid", "11"="Smooth muscle",
  "12"="Ventricular cardiomyocytes", "13"="Endothelial",
  "16"="Smooth muscle", "17"="Neuronal", "18"="Adipocyte",
  "19"="Mast", "20"="Endothelial",
  "22"="Fibroblasts", "23"="Ventricular cardiomyocytes",
  "24"="Myeloid"
)
overlap_keys <- intersect(names(new.cluster.ids), levels(car.combined))
if (length(overlap_keys) > 0) {
  car.combined <- RenameIdents(car.combined, new.cluster.ids[overlap_keys])
}
car.combined$cell_types <- Idents(car.combined)

DimPlot(car.combined, label = TRUE, repel = TRUE, raster = FALSE)

saveRDS(car.combined, file = "4.2.1.2/car_combined.rds")
