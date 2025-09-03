library(Seurat)

car.combined <- readRDS('4.2.1.2/car_combined.rds')
CM <- subset(car.combined, subset = cell_types == "Ventricular cardiomyocytes")

DefaultAssay(CM) <- "RNA"
CM.list <- SplitObject(CM, split.by = "patient")
CM.list <- lapply(CM.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x
})

features <- SelectIntegrationFeatures(object.list = CM.list)
anchors <- FindIntegrationAnchors(object.list = CM.list, anchor.features = features)
CM.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(CM.combined) <- "integrated"
CM.combined <- ScaleData(CM.combined, verbose = FALSE)
CM.combined <- RunPCA(CM.combined, npcs = 50, verbose = FALSE)
CM.combined <- RunUMAP(CM.combined, reduction = "pca", dims = 1:50)
CM.combined <- FindNeighbors(CM.combined, reduction = "pca", dims = 1:50)
CM.combined <- FindClusters(CM.combined, resolution = 0.3)

DefaultAssay(CM.combined) <- "RNA"
check_gene <- c(
  "NRXN3","MID1","CTNND2",
  "FGF12","CORIN","ERBB4",
  "TENM3","AUTS2","NTN1",
  "ACTA1","MYL2","NPPB",
  "KCNJ3","GRID2","GJA5"
)
DotPlot(CM.combined, features = check_gene)

new.cluster.ids <- c(
  "0"="CM2","1"="CM3","2"="CM1","3"="CM4","4"="CM1",
  "5"="CM2","6"="CM2","7"="CM1","8"="CM1","9"="CM2","10"="CM5"
)
names(new.cluster.ids) <- levels(CM.combined)
CM.combined <- RenameIdents(CM.combined, new.cluster.ids)
CM.combined$cell_subtypes <- Idents(CM.combined)
