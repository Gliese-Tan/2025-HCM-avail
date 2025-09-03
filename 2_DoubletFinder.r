suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(DoubletFinder)
})

read_counts <- function(sample) {
  Read10X(data.dir = file.path(sample, "soupX_matrix"))
}

make_obj <- function(sample, counts) {
  CreateSeuratObject(counts, min.cells = 3, min.features = 200, project = sample)
}

add_mt <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

qc_filter <- function(obj, nmin = 200, nmax = Inf, mtmax = 5) {
  subset(obj, subset = nFeature_RNA > nmin & nFeature_RNA < nmax & percent.mt < mtmax)
}

prep_for_df <- function(obj, pcs = 50, dims_use = 1:50) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use)
  obj <- FindClusters(obj, resolution = 1)
  obj <- RunUMAP(obj, dims = dims_use)
  obj
}

df_filter <- function(obj, dims_use = 1:50, pN = 0.25, rate = 0.075) {
  sweep.res <- paramSweep_v3(obj, PCs = dims_use, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- obj@meta.data$seurat_clusters
  homo <- modelHomotypic(annotations)
  nExp <- round(rate * ncol(obj))
  nExp.adj <- round(nExp * (1 - homo))
  obj2 <- doubletFinder_v3(obj, PCs = dims_use, pN = pN, pK = mpK, nExp = nExp, reuse.pANN = FALSE, sct = FALSE)
  df_col <- grep("^DF\\.classifications_", colnames(obj2@meta.data), value = TRUE)
  if (length(df_col) == 0) stop("No DF.classifications_ column found.")
  df_col <- tail(df_col, 1)
  obj2 <- obj2[, obj2@meta.data[[df_col]] == "Singlet", drop = FALSE]
  obj2
}

set_position <- function(obj, pos) {
  obj@meta.data$position <- pos
  obj
}

save_obj <- function(obj, sample) {
  dir.create(sample, showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, file = file.path(sample, paste0(sample, "_filtered.rds")))
}

cfg <- data.frame(
  sample = c(
    "H1_AL","H1_VA","H1_VS",
    "H2_AL","H2_VA","H2_VS",
    "H3_AL","H3_VA","H3_VS",
    "H4_AL","H4_VA","H4_VS",
    "H5_AL","H5_VA","H5_VS",
    "H6_AL","H6_VA","H6_VS",
    "H7_AL","H7_VA","H7_VS",
    "H8_AL","H8_VA","H8_VS"
  ),
  nmax = c(
    Inf, Inf, 8000,
    8000, 9000, Inf,
    Inf, Inf, Inf,
    10000, 10000, 10000,
    10000, Inf, 10000,
    10000, 10000, 10000,
    Inf, Inf, Inf,
    10000, 10000, 10000
  ),
  pcs = c(
    40,40,40,
    50,50,50,
    40,40,40,
    50,40,40,
    50,40,40,
    50,50,50,
    50,50,50,
    50,50,50
  ),
  rate = c(
    0.075,0.075,0.075,
    0.075,0.075,0.075,
    0.075,0.075,0.075,
    0.075,0.075,0.075,
    0.075,0.075,0.075,
    0.075,0.075,0.075,
    0.075,0.075,0.075,
    0.085,0.085,0.085
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(cfg))) {
  s <- cfg$sample[i]
  pos <- sub(".*_(AL|VA|VS)$", "\\1", s)
  counts <- read_counts(s)
  obj <- make_obj(s, counts)
  obj <- add_mt(obj)
  obj <- qc_filter(obj, nmin = 200, nmax = cfg$nmax[i], mtmax = 5)
  dims_use <- seq_len(cfg$pcs[i])
  obj <- prep_for_df(obj, pcs = cfg$pcs[i], dims_use = dims_use)
  obj <- df_filter(obj, dims_use = dims_use, pN = 0.25, rate = cfg$rate[i])
  obj <- set_position(obj, pos)
  save_obj(obj, s)
}
