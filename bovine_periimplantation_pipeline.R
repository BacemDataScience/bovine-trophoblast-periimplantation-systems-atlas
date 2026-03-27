# bovine_periimplantation_pipeline.R
# GSE234335 bovine peri-implantation scRNA-seq analysis
# Manuscript workflow:
# Integrative single-cell systems analysis reveals coordinated IFNT signaling
# and extracellular matrix remodeling during bovine peri-implantation

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# ============================================================
# 1) USER SETTINGS
# ============================================================
DATA_DIR <- "data/raw"
OUT_DIR  <- "results"

MIN_FEATURES <- 200
MAX_FEATURES <- 8000
MAX_MT       <- 20
MITO_PATTERN <- "^MT-"

DIMS       <- 1:30
RES_MAIN   <- 0.6
RES_SWEEP  <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
SEED       <- 123

DO_INTEGRATION    <- TRUE
DO_MARKERS_MAIN   <- TRUE
DO_COMPOSITION    <- TRUE
DO_STAGE_DE       <- TRUE
DO_MODULE_SCORES  <- TRUE
DO_MONOTONIC_GENES <- TRUE
DO_PSEUDOBULK_DE  <- TRUE
DO_PSEUDOTIME     <- TRUE
DO_RES_SWEEP      <- FALSE

DOWN_N_MARKERS <- 30000
TOPN_HEATMAP   <- 10
MAX_HEAT_GENES <- 250

set.seed(SEED)

# ============================================================
# 2) OUTPUT FOLDERS
# ============================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (d in c(
  "figures/main",
  "figures/supplement",
  "figures/insight",
  "tables",
  "tables/insight",
  "tables/supp",
  "objects",
  "metrics"
)) {
  dir.create(file.path(OUT_DIR, d), showWarnings = FALSE, recursive = TRUE)
}

# ============================================================
# 3) HELPERS
# ============================================================
parse_stage_rep <- function(prefix) {
  token <- sub("^.*_(D[0-9]+[A-Z]).*$", "\\1", prefix)
  stage <- sub("([A-Z])$", "", token)
  rep   <- sub("^D[0-9]+", "", token)
  list(stage = stage, rep = rep, token = token)
}

safe_require <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) message("Package missing: ", pkg)
  ok
}

detect_mito_pattern <- function(genes) {
  patterns <- c("^MT-", "^mt-", "^[Mm][Tt]-")
  for (p in patterns) if (any(grepl(p, genes))) return(p)
  "^MT-"
}

knn_same_sample_fraction <- function(obj, k = 30, reduction = "pca") {
  if (!safe_require("RANN")) return(NA_real_)
  emb <- Embeddings(obj, reduction)
  nn <- RANN::nn2(emb, k = k + 1)$nn.idx[, -1, drop = FALSE]
  samp <- obj$stage_rep
  fracs <- vapply(seq_len(nrow(nn)), function(i) mean(samp[nn[i, ]] == samp[i]), numeric(1))
  mean(fracs)
}

pick_fc_col <- function(df) {
  if ("avg_log2FC" %in% colnames(df)) return("avg_log2FC")
  if ("avg_logFC" %in% colnames(df)) return("avg_logFC")
  stop("No fold-change column found.")
}

diet_and_save <- function(obj_in, path) {
  obj_small <- tryCatch({
    DietSeurat(
      obj_in,
      assays = intersect(c("RNA", "integrated"), Assays(obj_in)),
      dimreducs = intersect(c("pca", "umap"), names(obj_in@reductions)),
      counts = TRUE,
      data = TRUE,
      scale.data = FALSE
    )
  }, error = function(e) obj_in)

  saveRDS(obj_small, path, compress = TRUE)
  rm(obj_small)
  gc()
}

save_dotplot <- function(obj_in, genes, group_by, fname, w = 14, h = 6) {
  g <- intersect(genes, rownames(obj_in))
  if (length(g) < 3) return(invisible(NULL))
  if (length(g) > 35) g <- g[1:35]
  p <- DotPlot(obj_in, features = g, group.by = group_by) + RotatedAxis()
  ggsave(file.path(OUT_DIR, "figures", "insight", fname), p, width = w, height = h, dpi = 300)
}

save_feature_panel <- function(obj_in, genes, fname, split_by = NULL, ncol = 4) {
  g <- intersect(genes, rownames(obj_in))
  if (length(g) < 2) return(invisible(NULL))
  if (length(g) > 16) g <- g[1:16]
  p <- FeaturePlot(obj_in, features = g, split.by = split_by, ncol = ncol, order = TRUE)
  ggsave(
    file.path(OUT_DIR, "figures", "insight", fname),
    p,
    width = 4 * ncol,
    height = 3 * ceiling(length(g) / ncol),
    dpi = 300
  )
}

# ============================================================
# 4) LOAD RAW MATRICES
# ============================================================
mtx_files <- list.files(DATA_DIR, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
stopifnot(length(mtx_files) > 0)

prefixes <- sub("_matrix\\.mtx\\.gz$", "", basename(mtx_files))
message("Found samples: ", paste(prefixes, collapse = ", "))

objs <- list()
for (pref in prefixes) {
  mtx_path  <- file.path(DATA_DIR, paste0(pref, "_matrix.mtx.gz"))
  feat_path <- file.path(DATA_DIR, paste0(pref, "_features.tsv.gz"))
  bc_path   <- file.path(DATA_DIR, paste0(pref, "_barcodes.tsv.gz"))

  if (!all(file.exists(c(mtx_path, feat_path, bc_path)))) next

  counts <- ReadMtx(mtx = mtx_path, features = feat_path, cells = bc_path)
  rownames(counts) <- make.unique(rownames(counts))

  meta <- parse_stage_rep(pref)
  gsm  <- sub("^(GSM[0-9]+)_.*$", "\\1", pref)

  obj0 <- CreateSeuratObject(counts, project = meta$token, min.cells = 3, min.features = 0)
  obj0$GSM       <- gsm
  obj0$sample_id <- pref
  obj0$stage     <- meta$stage
  obj0$replicate <- meta$rep
  obj0$stage_rep <- meta$token

  objs[[pref]] <- obj0
  rm(counts, obj0)
  gc()
}
stopifnot(length(objs) >= 2)

saveRDS(objs, file.path(OUT_DIR, "objects", "seurat_list_raw.rds"), compress = TRUE)

obj_raw <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "GSE234335")
saveRDS(obj_raw, file.path(OUT_DIR, "objects", "obj_merged_raw.rds"), compress = TRUE)

# ============================================================
# 5) QC
# ============================================================
genes <- rownames(obj_raw)
if (!any(grepl(MITO_PATTERN, genes))) {
  MITO_PATTERN <- detect_mito_pattern(genes)
}
obj_raw[["percent.mt"]] <- PercentageFeatureSet(obj_raw, pattern = MITO_PATTERN)

qc_before <- obj_raw@meta.data %>%
  mutate(cell = rownames(.)) %>%
  select(cell, GSM, sample_id, stage, replicate, stage_rep, nCount_RNA, nFeature_RNA, percent.mt)
write.csv(qc_before, file.path(OUT_DIR, "tables", "QC_before.csv"), row.names = FALSE)

p_qc_sample <- VlnPlot(
  obj_raw,
  c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "stage_rep",
  ncol = 3,
  pt.size = 0.1
) + plot_annotation(title = "QC by sample")
ggsave(file.path(OUT_DIR, "figures", "main", "Fig1_QC_by_sample.png"),
       p_qc_sample, width = 16, height = 5, dpi = 300)

p_qc_stage <- VlnPlot(
  obj_raw,
  c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "stage",
  ncol = 3,
  pt.size = 0.1
) + plot_annotation(title = "QC by stage")
ggsave(file.path(OUT_DIR, "figures", "main", "Fig1_QC_by_stage.png"),
       p_qc_stage, width = 12, height = 5, dpi = 300)

obj_qc <- subset(
  obj_raw,
  subset = nFeature_RNA >= MIN_FEATURES &
    nFeature_RNA <= MAX_FEATURES &
    percent.mt <= MAX_MT
)

qc_after <- obj_qc@meta.data %>%
  mutate(cell = rownames(.)) %>%
  select(cell, GSM, sample_id, stage, replicate, stage_rep, nCount_RNA, nFeature_RNA, percent.mt)
write.csv(qc_after, file.path(OUT_DIR, "tables", "QC_after.csv"), row.names = FALSE)

qc_summary <- qc_before %>%
  count(stage_rep, name = "cells_before") %>%
  left_join(qc_after %>% count(stage_rep, name = "cells_after"), by = "stage_rep") %>%
  mutate(retained = cells_after / cells_before)
write.csv(qc_summary, file.path(OUT_DIR, "tables", "QC_summary_by_sample.csv"), row.names = FALSE)

saveRDS(obj_qc, file.path(OUT_DIR, "objects", "obj_afterQC.rds"), compress = TRUE)

# ============================================================
# 6) MERGE-ONLY BASELINE
# ============================================================
obj_m <- NormalizeData(obj_qc, verbose = FALSE)
obj_m <- FindVariableFeatures(obj_m, nfeatures = 2000, verbose = FALSE)
obj_m <- ScaleData(obj_m, features = VariableFeatures(obj_m), verbose = FALSE)
obj_m <- RunPCA(obj_m, features = VariableFeatures(obj_m), verbose = FALSE)
obj_m <- FindNeighbors(obj_m, dims = DIMS, verbose = FALSE)
obj_m <- FindClusters(obj_m, resolution = RES_MAIN, verbose = FALSE)
obj_m <- RunUMAP(obj_m, dims = DIMS, verbose = FALSE)

saveRDS(obj_m, file.path(OUT_DIR, "objects", "obj_mergeOnly.rds"), compress = TRUE)

p_m1 <- DimPlot(obj_m, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Merge-only clusters")
p_m2 <- DimPlot(obj_m, group.by = "stage") + ggtitle("Merge-only by stage")
p_m3 <- DimPlot(obj_m, group.by = "stage_rep") + ggtitle("Merge-only by sample")
ggsave(file.path(OUT_DIR, "figures", "supplement", "Sup_UMAP_mergeOnly_panels.png"),
       p_m1 + p_m2 + p_m3, width = 18, height = 6, dpi = 300)

mix_m <- knn_same_sample_fraction(obj_m, k = 30, reduction = "pca")

# ============================================================
# 7) INTEGRATION
# ============================================================
obj <- obj_m
mix_i <- NA_real_

if (DO_INTEGRATION) {
  obj_list <- SplitObject(obj_qc, split.by = "stage_rep")
  obj_list <- lapply(obj_list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    x
  })

  anchors <- FindIntegrationAnchors(object.list = obj_list, dims = DIMS)
  obj_i <- IntegrateData(anchorset = anchors, dims = DIMS)

  DefaultAssay(obj_i) <- "integrated"
  obj_i <- FindVariableFeatures(obj_i, nfeatures = 2000, verbose = FALSE)
  obj_i <- ScaleData(obj_i, features = VariableFeatures(obj_i), verbose = FALSE)
  obj_i <- RunPCA(obj_i, features = VariableFeatures(obj_i), verbose = FALSE)
  obj_i <- FindNeighbors(obj_i, dims = DIMS, verbose = FALSE)
  obj_i <- FindClusters(obj_i, resolution = RES_MAIN, verbose = FALSE)
  obj_i <- RunUMAP(obj_i, dims = DIMS, verbose = FALSE)

  saveRDS(obj_i, file.path(OUT_DIR, "objects", "obj_integrated.rds"), compress = TRUE)

  p_i1 <- DimPlot(obj_i, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Integrated clusters")
  p_i2 <- DimPlot(obj_i, group.by = "stage") + ggtitle("Integrated by stage")
  p_i3 <- DimPlot(obj_i, group.by = "stage_rep") + ggtitle("Integrated by sample")
  ggsave(file.path(OUT_DIR, "figures", "main", "Fig2_UMAP_integrated_panels.png"),
         p_i1 + p_i2 + p_i3, width = 18, height = 6, dpi = 300)

  mix_i <- knn_same_sample_fraction(obj_i, k = 30, reduction = "pca")
  obj <- obj_i
}

mix_tbl <- data.frame(
  workflow = c("merge_only", "integrated"),
  knn_same_sample_fraction = c(mix_m, mix_i)
)
write.csv(mix_tbl, file.path(OUT_DIR, "metrics", "mixing_metric_knn_same_sample.csv"), row.names = FALSE)

if ("stage" %in% colnames(obj@meta.data)) {
  obj$stage <- factor(as.character(obj$stage), levels = c("D12", "D14", "D16", "D18"))
}
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

writeLines(c(
  paste0("Time: ", Sys.time()),
  paste0("DefaultAssay(obj): ", DefaultAssay(obj)),
  paste0("Assays: ", paste(Assays(obj), collapse = ", ")),
  paste0("Cells: ", ncol(obj)),
  paste0("Genes: ", nrow(obj))
), con = file.path(OUT_DIR, "metrics", "sanity.txt"))

diet_and_save(obj, file.path(OUT_DIR, "objects", "CHECKPOINT_afterJoinLayers_dropScale_diet.rds"))

# ============================================================
# 8) OPTIONAL RESOLUTION SWEEP
# ============================================================
if (DO_RES_SWEEP && safe_require("clustree")) {
  library(clustree)
  obj_tmp <- obj
  DefaultAssay(obj_tmp) <- if ("integrated" %in% Assays(obj_tmp)) "integrated" else "RNA"
  obj_tmp <- FindNeighbors(obj_tmp, dims = DIMS, verbose = FALSE)
  for (r in RES_SWEEP) obj_tmp <- FindClusters(obj_tmp, resolution = r, verbose = FALSE)
  prefix <- paste0(DefaultAssay(obj_tmp), "_snn_res.")
  p_ct <- clustree(obj_tmp, prefix = prefix) + ggtitle("Clustree resolution sweep")
  ggsave(file.path(OUT_DIR, "figures", "supplement", "Sup_clustree_resolution.png"),
         p_ct, width = 12, height = 10, dpi = 300)
  rm(obj_tmp)
  gc()
}

# ============================================================
# 9) MARKERS
# ============================================================
if (DO_MARKERS_MAIN) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)

  obj_mark <- obj
  if (ncol(obj_mark) > DOWN_N_MARKERS) {
    set.seed(SEED)
    obj_mark <- subset(obj_mark, cells = sample(colnames(obj_mark), DOWN_N_MARKERS))
  }

  Idents(obj_mark) <- obj_mark$seurat_clusters
  markers <- FindAllMarkers(obj_mark, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)

  if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
  if (!("cluster" %in% colnames(markers))) {
    if ("group" %in% colnames(markers)) markers$cluster <- markers$group
    if ("ident" %in% colnames(markers)) markers$cluster <- markers$ident
  }

  write.csv(markers, file.path(OUT_DIR, "tables", "cluster_markers_all_DOWN.csv"), row.names = FALSE)

  fc_col <- pick_fc_col(markers)

  topN <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = TOPN_HEATMAP, with_ties = FALSE)

  heat_genes <- unique(topN$gene)
  heat_genes <- intersect(heat_genes, rownames(obj_mark))
  if (length(heat_genes) > MAX_HEAT_GENES) heat_genes <- heat_genes[1:MAX_HEAT_GENES]

  p_heat <- DoHeatmap(obj_mark, features = heat_genes, size = 3) +
    ggtitle("Top markers per cluster")
  ggsave(file.path(OUT_DIR, "figures", "main", "Fig3_TopMarkers_heatmap.png"),
         p_heat, width = 12, height = 10, dpi = 300)

  write.csv(topN, file.path(OUT_DIR, "tables", "TopMarkers_perCluster_DOWN.csv"), row.names = FALSE)
  diet_and_save(obj, file.path(OUT_DIR, "objects", "CHECKPOINT_afterStep6_diet.rds"))
}

# ============================================================
# 10) STAGE COMPOSITION
# ============================================================
if (DO_COMPOSITION) {
  comp <- obj@meta.data %>%
    count(stage, seurat_clusters) %>%
    group_by(stage) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  write.csv(comp, file.path(OUT_DIR, "tables", "composition_by_stage_cluster.csv"), row.names = FALSE)

  p_c1 <- ggplot(comp, aes(x = seurat_clusters, y = prop, fill = stage)) +
    geom_col(position = "dodge") +
    ggtitle("Cluster composition shift across D12–D18") +
    xlab("Cluster") + ylab("Proportion")
  ggsave(file.path(OUT_DIR, "figures", "main", "Fig5_Composition_by_stage.png"),
         p_c1, width = 14, height = 6, dpi = 300)
}

# ============================================================
# 11) GLOBAL DE
# ============================================================
if (DO_STAGE_DE) {
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  Idents(obj) <- obj$stage

  do_deg <- function(a, b, outname) {
    deg <- FindMarkers(obj, ident.1 = a, ident.2 = b, logfc.threshold = 0.10, min.pct = 0.10)
    deg$gene <- rownames(deg)
    write.csv(deg, file.path(OUT_DIR, "tables", outname), row.names = FALSE)
    deg
  }

  deg_D18_D12 <- do_deg("D18", "D12", "DEG_D18_vs_D12_global.csv")

  fc_deg <- pick_fc_col(deg_D18_D12)
  deg2 <- deg_D18_D12 %>%
    mutate(sig = (p_val_adj < 0.05 & abs(.data[[fc_deg]]) >= 0.10))

  p_vol <- ggplot(deg2, aes(x = .data[[fc_deg]], y = -log10(p_val_adj + 1e-300))) +
    geom_point(aes(alpha = sig)) +
    ggtitle("Global DEG: D18 vs D12") +
    xlab(fc_deg) + ylab("-log10(adj p)")
  ggsave(file.path(OUT_DIR, "figures", "main", "Fig6_Volcano_D18_vs_D12.png"),
         p_vol, width = 7, height = 6, dpi = 300)
}

# ============================================================
# 12) MONOTONIC GENES
# ============================================================
if (DO_MONOTONIC_GENES) {
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

  avg <- AverageExpression(obj, group.by = "stage", assays = "RNA", layer = "data")$RNA
  stage_order <- c("D12", "D14", "D16", "D18")
  avg2 <- avg[, stage_order, drop = FALSE]
  write.csv(avg2, file.path(OUT_DIR, "tables", "insight", "AverageExpression_by_stage.csv"))

  stage_num <- 1:4
  sp <- apply(avg2, 1, function(v) suppressWarnings(cor(v, stage_num, method = "spearman")))
  mono <- data.frame(
    gene = rownames(avg2),
    spearman_stage = sp,
    D12 = avg2[, "D12"],
    D14 = avg2[, "D14"],
    D16 = avg2[, "D16"],
    D18 = avg2[, "D18"]
  )
  mono <- mono[order(mono$spearman_stage, decreasing = TRUE), ]
  write.csv(mono, file.path(OUT_DIR, "tables", "insight", "Monotonic_genes_ranked.csv"), row.names = FALSE)
}

# ============================================================
# 13) MODULE SCORES
# ============================================================
if (DO_MODULE_SCORES) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)

  panel <- list(
    Lineage_TB  = c("KRT8","KRT18","KRT19","EPCAM","CDX2","GATA3","TFAP2C","TEAD4","EOMES"),
    Lineage_HB  = c("GATA6","SOX17","PDGFRA","FOXA2","COL4A1","LAMA1"),
    Lineage_EPI = c("POU5F1","SOX2","NANOG","OTX2"),
    IFN_ISG     = c("IFNT","ISG15","MX1","MX2","OAS1","OAS2","OASL","IFI6","IFIT1","IFIT2","IFIT3","BST2","STAT1","IRF7"),
    ECM_ADH     = c("SPP1","FN1","VCAN","ITGAV","ITGB1","ITGB3","COL1A1","COL1A2","COL3A1","COL4A1","LAMA1","LAMC1","MMP2","MMP9","TIMP1","TIMP2","LGALS3"),
    EMT         = c("VIM","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","CDH1","CDH2","EPCAM")
  )

  for (nm in names(panel)) {
    save_dotplot(obj, panel[[nm]], group_by = "stage", fname = paste0("Dot_", nm, "_byStage.png"))
    save_dotplot(obj, panel[[nm]], group_by = "seurat_clusters", fname = paste0("Dot_", nm, "_byCluster.png"), w = 16, h = 7)
    save_feature_panel(obj, panel[[nm]], fname = paste0("FeaturePanel_", nm, ".png"))
  }

  score_cols <- c()
  for (nm in names(panel)) {
    g <- intersect(panel[[nm]], rownames(obj))
    if (length(g) >= 5) {
      obj <- AddModuleScore(obj, features = list(g), name = paste0("Score_", nm))
      score_cols <- c(score_cols, paste0("Score_", nm, "1"))
    }
  }

  if (length(score_cols)) {
    write.csv(
      obj@meta.data[, c("stage_rep", "stage", "seurat_clusters", score_cols), drop = FALSE],
      file.path(OUT_DIR, "tables", "insight", "ModuleScores.csv"),
      row.names = FALSE
    )
  }
}

# ============================================================
# 14) PSEUDOBULK
# ============================================================
if (DO_PSEUDOBULK_DE) {
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

  counts <- GetAssayData(obj, layer = "counts")
  meta <- obj@meta.data
  grp <- as.character(meta$stage_rep)
  grp_levels <- unique(grp)

  pb_mat <- sapply(grp_levels, function(g) {
    cells <- rownames(meta)[grp == g]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  pb_mat <- Matrix::Matrix(pb_mat, sparse = TRUE)

  pb_meta <- data.frame(
    stage_rep = grp_levels,
    stage = sub("([A-Z])$", "", grp_levels),
    rep = sub("^D[0-9]+", "", grp_levels),
    stringsAsFactors = FALSE
  )
  pb_meta$stage <- factor(pb_meta$stage, levels = c("D12", "D14", "D16", "D18"))

  pb <- CreateSeuratObject(pb_mat, meta.data = pb_meta, project = "PSEUDOBULK")
  pb <- NormalizeData(pb, verbose = FALSE)
  pb_data <- GetAssayData(pb, layer = "data")

  cor_mat <- cor(as.matrix(pb_data), method = "spearman")
  write.csv(cor_mat, file.path(OUT_DIR, "tables", "insight", "PSEUDOBULK_sample_cor_spearman.csv"))

  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("sample1", "sample2", "spearman")
  p_cor <- ggplot(cor_df, aes(x = sample1, y = sample2, fill = spearman)) +
    geom_tile() +
    ggtitle("Pseudobulk sample Spearman correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_DIR, "figures", "insight", "Insight_PSEUDOBULK_sample_cor_heat.png"),
         p_cor, width = 10, height = 8, dpi = 300)

  saveRDS(pb, file.path(OUT_DIR, "objects", "pb_stage_rep.rds"), compress = TRUE)
}

# ============================================================
# 15) REVISED MANUSCRIPT ADDITIONAL ANALYSES
# ============================================================
# IFNT heterogeneity statistics, pseudotime inference,
# and supplementary table export

if (DO_PSEUDOTIME) {
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)

  ADDON_OUT <- file.path(OUT_DIR, "ADDON_reviewer_upgrades")
  dir.create(ADDON_OUT, showWarnings = FALSE, recursive = TRUE)
  for (d in c("figures", "tables", "tables/supp", "objects", "metrics")) {
    dir.create(file.path(ADDON_OUT, d), showWarnings = FALSE, recursive = TRUE)
  }

  MODULES <- list(
    IFN_ISG = c("IFNT","IFNT2","IFNT3","ISG15","MX1","MX2","OAS1","OAS2","OASL","IFI6","IFIT1","IFIT2","IFIT3","BST2","STAT1","IRF7"),
    ECM_ADH = c("SPP1","FN1","VCAN","ITGAV","ITGB1","ITGB3","COL1A1","COL1A2","COL3A1","COL4A1","LAMA1","LAMC1","MMP2","MMP9","TIMP1","TIMP2","LGALS3"),
    EMT1 = c("VIM","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","CDH1","CDH2","EPCAM"),
    TB_LINEAGE = c("KRT8","KRT18","KRT19","EPCAM","CDX2","GATA3","TFAP2C","TEAD4","EOMES")
  )
  MODULES <- lapply(MODULES, function(g) intersect(g, rownames(obj)))

  mod_tbl <- bind_rows(lapply(names(MODULES), function(nm) {
    data.frame(module = nm, gene = MODULES[[nm]], stringsAsFactors = FALSE)
  }))
  write.csv(mod_tbl, file.path(ADDON_OUT, "tables/supp", "Supplement_ModuleGeneLists.csv"), row.names = FALSE)

  obj <- NormalizeData(obj, verbose = FALSE)
  score_names <- character(0)
  for (nm in names(MODULES)) {
    g <- MODULES[[nm]]
    if (length(g) >= 5) {
      obj <- AddModuleScore(obj, features = list(g), name = paste0("Score_", nm))
      score_names <- c(score_names, paste0("Score_", nm, "1"))
    }
  }

  if (length(score_names) > 0) {
    ms <- obj@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      transmute(cell, stage, stage_rep, cluster = as.character(seurat_clusters)) %>%
      bind_cols(obj@meta.data %>% dplyr::select(dplyr::all_of(score_names)))
    write.csv(ms, file.path(ADDON_OUT, "tables/supp", "Supplement_ModuleScores_perCell.csv"), row.names = FALSE)
  }

  ifnt_gene <- c("IFNT", "IFNT2", "IFNT3")[c("IFNT", "IFNT2", "IFNT3") %in% rownames(obj)][1]
  if (!is.na(ifnt_gene)) {
    ifnt_expr <- FetchData(obj, vars = ifnt_gene)[, 1]
    obj$IFNT_pos <- ifnt_expr > 0

    tab_cluster <- table(cluster = obj$seurat_clusters, IFNT_pos = obj$IFNT_pos)
    write.csv(as.data.frame.matrix(tab_cluster),
              file.path(ADDON_OUT, "tables/supp", "Supplement_IFNTpos_table_byCluster.csv"))

    chisq_res <- chisq.test(tab_cluster)
    writeLines(
      paste0("Chi-square test IFNT_pos ~ cluster: p=", signif(chisq_res$p.value, 3), " | IFNT gene=", ifnt_gene),
      con = file.path(ADDON_OUT, "metrics", "IFNT_heterogeneity_stats.txt")
    )

    ifnt_by_cluster_stage <- obj@meta.data %>%
      mutate(cluster = as.character(seurat_clusters)) %>%
      group_by(stage, cluster) %>%
      summarise(n_cells = n(), pct_IFNT_pos = 100 * mean(IFNT_pos), .groups = "drop") %>%
      arrange(stage, suppressWarnings(as.integer(cluster)))

    write.csv(ifnt_by_cluster_stage,
              file.path(ADDON_OUT, "tables/supp", "Supplement_IFNTpos_percent_byStageCluster.csv"),
              row.names = FALSE)

    ifn_score_col <- grep("^Score_IFN_ISG1$", colnames(obj@meta.data), value = TRUE)
    if (length(ifn_score_col) == 1) {
      kw <- kruskal.test(obj@meta.data[[ifn_score_col]] ~ as.factor(obj$seurat_clusters))
      writeLines(
        paste0("Kruskal–Wallis IFN score across clusters: p=", signif(kw$p.value, 3)),
        con = file.path(ADDON_OUT, "metrics", "IFNscore_KW_stats.txt")
      )
    }

    p_ifnt_bar <- ggplot(ifnt_by_cluster_stage, aes(x = cluster, y = pct_IFNT_pos, fill = stage)) +
      geom_col(position = position_dodge(width = 0.8)) +
      labs(title = paste0("IFNT-positive fraction by cluster and stage (", ifnt_gene, ")"),
           x = "Cluster", y = "% IFNT+ cells") +
      theme_classic(base_size = 11)

    ggsave(file.path(ADDON_OUT, "figures", "ADDON_IFNTpos_byClusterStage.png"),
           p_ifnt_bar, width = 12, height = 5, dpi = 400)
  }

  if (safe_require("SingleCellExperiment") && safe_require("slingshot")) {
    library(SingleCellExperiment)
    library(slingshot)

    tb_score_col <- grep("^Score_TB_LINEAGE1$", colnames(obj@meta.data), value = TRUE)
    if (length(tb_score_col) == 1) {
      thr <- quantile(obj@meta.data[[tb_score_col]], 0.75, na.rm = TRUE)
      tb_cells <- rownames(obj@meta.data)[obj@meta.data[[tb_score_col]] > thr]

      if (length(tb_cells) >= 1000) {
        obj_tb <- subset(obj, cells = tb_cells)

        if (!"umap" %in% names(obj_tb@reductions)) {
          obj_tb <- RunUMAP(obj_tb, reduction = "pca", dims = DIMS)
        }

        sce <- as.SingleCellExperiment(obj_tb)
        reducedDim(sce, "UMAP") <- Embeddings(obj_tb, "umap")
        cl <- as.factor(obj_tb$seurat_clusters)
        sce <- slingshot(sce, clusterLabels = cl, reducedDim = "UMAP")

        pt <- slingPseudotime(sce)
        pt1 <- pt[, 1]
        obj_tb$pseudotime <- pt1[match(colnames(obj_tb), rownames(pt))]

        pt_tbl <- obj_tb@meta.data %>%
          tibble::rownames_to_column("cell") %>%
          transmute(cell, stage, stage_rep, cluster = as.character(seurat_clusters), pseudotime)
        write.csv(pt_tbl, file.path(ADDON_OUT, "tables/supp", "Supplement_Pseudotime_TBsubset.csv"), row.names = FALSE)

        p_pt <- FeaturePlot(obj_tb, features = "pseudotime", order = TRUE) +
          ggtitle("TB-enriched subset pseudotime (Slingshot, UMAP)")
        ggsave(file.path(ADDON_OUT, "figures", "ADDON_Pseudotime_UMAP.png"),
               p_pt, width = 7, height = 6, dpi = 400)

        ifnt_gene_tb <- c("IFNT", "IFNT2", "IFNT3")[c("IFNT", "IFNT2", "IFNT3") %in% rownames(obj_tb)][1]
        if (!is.na(ifnt_gene_tb)) {
          ifnt_expr_tb <- FetchData(obj_tb, vars = ifnt_gene_tb)[, 1]
          df_tr <- data.frame(
            pseudotime = obj_tb$pseudotime,
            IFNT_expr = ifnt_expr_tb,
            stage = as.character(obj_tb$stage),
            cluster = as.character(obj_tb$seurat_clusters),
            stringsAsFactors = FALSE
          )

          write.csv(df_tr,
                    file.path(ADDON_OUT, "tables/supp", "Supplement_IFNT_along_pseudotime_TBsubset.csv"),
                    row.names = FALSE)

          ct <- cor.test(df_tr$IFNT_expr, df_tr$pseudotime, method = "spearman", exact = FALSE)
          sink(file.path(ADDON_OUT, "metrics", "IFNT_vs_pseudotime_spearman.txt"))
          cat("Gene: ", ifnt_gene_tb, "\n")
          print(ct)
          sink()

          p_tr <- ggplot(df_tr, aes(x = pseudotime, y = IFNT_expr)) +
            geom_point(alpha = 0.08, size = 0.35) +
            geom_smooth(se = FALSE) +
            theme_classic(base_size = 11) +
            labs(title = paste0(ifnt_gene_tb, " expression along pseudotime (TB subset)"),
                 x = "Pseudotime", y = "RNA (log-normalized)")
          ggsave(file.path(ADDON_OUT, "figures", "ADDON_IFNT_along_pseudotime.png"),
                 p_tr, width = 7.5, height = 5.5, dpi = 400)
        }

        sum_stage <- pt_tbl %>%
          group_by(stage) %>%
          summarise(
            n = n(),
            mean = mean(pseudotime, na.rm = TRUE),
            median = median(pseudotime, na.rm = TRUE),
            IQR = IQR(pseudotime, na.rm = TRUE),
            min = min(pseudotime, na.rm = TRUE),
            max = max(pseudotime, na.rm = TRUE),
            .groups = "drop"
          )
        write.csv(sum_stage,
                  file.path(ADDON_OUT, "tables/supp", "Supplement_Pseudotime_summary_byStage.csv"),
                  row.names = FALSE)

        kw_pt <- kruskal.test(pseudotime ~ stage, data = pt_tbl)
        pw_pt <- pairwise.wilcox.test(pt_tbl$pseudotime, pt_tbl$stage, p.adjust.method = "BH")
        sink(file.path(ADDON_OUT, "metrics", "Pseudotime_stage_tests.txt"))
        cat("Kruskal–Wallis pseudotime ~ stage\n")
        print(kw_pt)
        cat("\nPairwise Wilcoxon (BH-adjusted)\n")
        print(pw_pt)
        sink()

        comp <- obj@meta.data %>%
          mutate(cluster = as.character(seurat_clusters)) %>%
          count(stage, cluster) %>%
          group_by(stage) %>%
          mutate(prop = n / sum(n)) %>%
          ungroup()
        write.csv(comp,
                  file.path(ADDON_OUT, "tables/supp", "Supplement_Composition_byStageCluster.csv"),
                  row.names = FALSE)

        saveRDS(obj_tb, file.path(ADDON_OUT, "objects", "obj_TBsubset_pseudotime_diet.rds"), compress = TRUE)
      }
    }
  }
}

# ============================================================
# 16) FINAL SAVE
# ============================================================
diet_and_save(obj, file.path(OUT_DIR, "objects", "obj_FINAL_diet.rds"))
writeLines(paste0("DONE: ", Sys.time()), con = file.path(OUT_DIR, "metrics", "DONE.txt"))
message("Pipeline finished. Outputs written to: ", normalizePath(OUT_DIR))
