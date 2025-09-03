suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(DropletUtils)
})

run_soupx <- function(toc, tod, rho = NULL) {
  tod <- tod[rownames(toc), , drop = FALSE]
  all <- CreateSeuratObject(toc)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
  all <- ScaleData(all, features = rownames(all))
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = FALSE)
  all <- FindNeighbors(all, dims = 1:30)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:30)
  matx <- all@meta.data
  sc <- SoupChannel(tod, toc)
  sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
  if (is.null(rho)) {
    sc <- tryCatch(autoEstCont(sc), error = function(e) setContaminationFraction(sc, 0.2))
  } else {
    sc <- setContaminationFraction(sc, rho)
  }
  adjustCounts(sc, roundToInt = TRUE)
}

run_one <- function(sample_name,
                    filtered_suffix = "/filtered_feature_bc_matrix",
                    raw_suffix = "/raw_feature_bc_matrix",
                    out_root = "process/0.seurat_project/SoupX") {
  toc <- Read10X(data.dir = paste0(sample_name, filtered_suffix))
  tod <- Read10X(data.dir = paste0(sample_name, raw_suffix))
  out_mat <- run_soupx(toc, tod, rho = NULL)
  out_dir <- file.path(out_root, sample_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  DropletUtils:::write10xCounts(file.path(out_dir, "soupX_matrix"), out_mat, version = "3")
}

# 病人编号 H1–H8
samples <- paste0("H", 1:8)
# 区域
regions <- c("AL", "VA", "VS")

# 循环批量运行
for (s in samples) {
  message("---------- ", s, " ----------")
  for (r in regions) {
    sample_name <- paste0(s, "_", r)
    run_one(sample_name)
  }
}
