# SECCIÓN 0: CONFIGURACIÓN

PAQUETES <- c(CRAN = "Seurat", CRAN = "Matrix", CRAN = "data.table", CRAN = "ggplot2",
              CRAN = "lmerTest", CRAN = "emmeans", GitHub = "BayesPrism")
faltan <- PAQUETES[!vapply(PAQUETES, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) {
  cmd <- c(
    if (any(names(faltan) == "CRAN"))
      sprintf('install.packages(c(%s))', paste0('"', faltan[names(faltan) == "CRAN"], '"', collapse = ", ")),
    if (any(names(faltan) == "Bioc"))
      sprintf('BiocManager::install(c(%s))', paste0('"', faltan[names(faltan) == "Bioc"], '"', collapse = ", ")),
    if (any(names(faltan) == "GitHub"))
      'devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")')
  stop("Faltan paquetes: ", paste(faltan, collapse = ", "),
       "\nInstalar con:\n  ", paste(cmd, collapse = "\n  "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(BayesPrism)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

BASE_DIR  <- "Resultados/Sc vs Sc"
REF_A_DIR <- file.path("datos", "GSE137846")
NEW_DIR   <- file.path("datos", "Data")
OUT_DIR   <- file.path(BASE_DIR, "resultados_bayesprism")

SEED                  <- 1234
N_CORES               <- 4
MAX_CELLS_PER_SAMPLE  <- 6000
MAX_CELLS_PER_TYPE    <- 1500
MIN_CELLS_PER_TYPE    <- 30
QC_MIN_FEATURES <- 200
QC_MIN_COUNTS   <- 400
QC_MAX_MT       <- 20
REF_A_CELLTYPE_COL <- NULL

SAMPLE_META <- data.frame(
  folder = c("Ldiab230111con", "LnonD230921cont",
             "YnonD231226cont", "YnonD240126con1", "YnonD240126con2",
             "Ycont20240912le",
             "Ycont241029le", "Ycont241029ri",
             "ydiab20240517le", "ydiab20240517ri",
             "Ydia20250213LE", "Ydia20250213RI",
             "YdiabR240521le", "YdiabR240521ri",
             "YdiabR240627le", "YdiabR240627ri",
             "YdiabR240702le", "YdiabR240702ri",
             "YdiabR241031le", "YdiabR241031ri"),
  donor  = c("D1_DM_live", "D2_nonDM_live",
             "Dlive_nonDM_231226", "Dlive_nonDM_240126a", "Dlive_nonDM_240126b",
             "D6_nonDM_cad",
             "D7_nonDM_cad", "D7_nonDM_cad",
             "DM_cad_20240517", "DM_cad_20240517",
             "DM_cad_20250213", "DM_cad_20250213",
             "DR_cad_240521", "DR_cad_240521",
             "DR_cad_240627", "DR_cad_240627",
             "DR_cad_240702", "DR_cad_240702",
             "DR_cad_241031", "DR_cad_241031"),
  source = c("Live", "Live", "Live", "Live", "Live",
             rep("Cadaver", 15)),
  group  = c("DM", "nonDM", "nonDM", "nonDM", "nonDM",
             "nonDM", "nonDM", "nonDM",
             "DM", "DM", "DM", "DM",
             rep("DR", 8)),
  eye    = c(NA, NA, NA, NA, NA, "L", "L", "R", "L", "R", "L", "R",
             rep(c("L", "R"), 4)),
  sex_expected = c("M", "M", NA, NA, NA, "M", "F", "F", "F", "F", "F", "F",
                   rep("M", 8)),
  stringsAsFactors = FALSE
)
SAMPLE_META$group <- factor(SAMPLE_META$group, levels = c("nonDM", "DM", "DR"))

CELLTYPE_ORDER <- c("Rod", "Cone", "BC", "AC", "HC", "RGC", "MGC", "AST",
                    "Microglia", "Vascular", "RPE", "TC", "Neutrophil")

MARKERS <- list(
  Rod        = c("RHO", "PDE6G", "GNGT1", "NRL", "SAG"),
  Cone       = c("ARR3", "GNAT2", "GNGT2", "PDE6H", "OPN1SW"),
  BC         = c("GRM6", "GRIK1", "VSX2", "PCP2", "TRPM1"),
  AC         = c("TFAP2B", "SLC32A1", "GAD1", "C1QL2", "TFAP2A"),
  HC         = c("LHX1", "ONECUT1", "ONECUT2", "TNR"),
  RGC        = c("RBPMS", "NEFL", "SNCG", "POU4F2", "THY1"),
  MGC        = c("RLBP1", "DKK3", "GLUL", "CRYAB", "APOE"),
  AST        = c("GFAP", "PAX2", "S100B", "AQP4"),
  Microglia  = c("C1QA", "C1QB", "TYROBP", "CD74", "AIF1"),
  Vascular   = c("PECAM1", "VWF", "CLDN5", "RGS5", "ACTA2"),
  RPE        = c("RPE65", "BEST1", "TTR", "PMEL"),
  TC         = c("CD3D", "CD3E", "CD52", "CD69", "IL7R"),
  Neutrophil = c("CXCR2", "CSF3R", "FCGR3B", "S100A12", "S100A8")
)

DO_DOUBLETS              <- TRUE
DO_LABEL_TRANSFER_CHECK  <- TRUE
CLUSTER_RESOLUTION       <- 0.6
RUN_Y_NONDM_LODO         <- TRUE
RUN_Y_TODAS_LODO         <- TRUE
MAX_CELLS_PER_TYPE_LODO  <- 600
MIN_CELLS_PER_STATE      <- 20
FORZAR_REF_A             <- FALSE
FORZAR_YANG              <- FALSE

set.seed(SEED)
REF_OUT <- file.path(OUT_DIR, "01_referencias")
DEC_OUT <- file.path(OUT_DIR, "02_deconvolucion")
FIG_OUT <- file.path(OUT_DIR, "figuras")
TAB_OUT <- file.path(OUT_DIR, "tablas")

# SECCIÓN 1: FUNCIONES AUXILIARES

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (d in c("01_referencias", "02_deconvolucion", "figuras", "tablas"))
  dir.create(file.path(OUT_DIR, d), showWarnings = FALSE)

msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

find_file <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) stop("No se encontro '", pattern, "' en ", dir)
  plain <- f[!grepl("\\.gz$", f)]
  if (length(plain) > 0) plain[1] else f[1]
}

collapse_rows <- function(m) {
  if (!anyDuplicated(rownames(m))) return(m)
  f <- factor(rownames(m), levels = unique(rownames(m)))
  agg <- Matrix::sparseMatrix(i = as.integer(f), j = seq_along(f), x = 1,
                              dims = c(nlevels(f), length(f)))
  out <- agg %*% m
  rownames(out) <- levels(f)
  out
}

harmonize_label <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  rules <- list(
    BC         = "bipolar|^bcs?$|^bps?$|^bcs?\\b|^db\\d|^fmb|^imb|^rbc",
    AC         = "amacrine|^ac$|^acs?\\b",
    HC         = "horizontal|^hcs?$",
    RGC        = "ganglion|^rgc|midget|parasol",
    MGC        = "m[uü]ller|muller|^mgcs?$|macroglia",
    AST        = "astro",
    Microglia  = "microglia|macrophage|myeloid",
    Vascular   = "endothel|^endo|vascular|pericyte|smooth muscle|^vsmc",
    RPE        = "^rpe|pigment",
    TC         = "t cell|t-cell|^t$|lymphocyte",
    Neutrophil = "neutrophil",
    Rod        = "rod",
    Cone       = "cone|^ml_cone|^s_cone"
  )
  for (ct in names(rules)) {
    hit <- is.na(out) & grepl(rules[[ct]], y, perl = TRUE)
    out[hit] <- ct
  }
  out
}

subsample_by_type <- function(meta, max_n) {
  idx <- unlist(lapply(split(seq_len(nrow(meta)), meta$celltype), function(ii) {
    if (length(ii) <= max_n) return(ii)
    by_d  <- split(ii, meta$donor[ii])
    avail <- lengths(by_d)
    alloc <- integer(length(by_d))
    while (sum(alloc) < max_n) {
      open <- which(alloc < avail)
      add  <- min(length(open), max_n - sum(alloc))
      alloc[open[seq_len(add)]] <- alloc[open[seq_len(add)]] + 1L
    }
    unlist(mapply(function(v, k) v[sample.int(length(v), k)], by_d, alloc, SIMPLIFY = FALSE))
  }))
  meta[sort(idx), ]
}

deconvolve <- function(R, cells, mixture, tag, max_n, state_from = c("cluster", "label")) {
  state_from <- match.arg(state_from)
  meta <- R$meta[match(cells, R$meta$cell), ]
  meta <- meta[meta$celltype %in% COMMON_TYPES, ]
  set.seed(SEED)
  meta <- subsample_by_type(meta, max_n)
  st <- if (state_from == "label") paste0(meta$celltype, "__", meta$label_orig) else
    paste0(meta$celltype, "__c", meta$cluster)
  small <- st %in% names(which(table(st) < MIN_CELLS_PER_STATE))
  st[small] <- paste0(meta$celltype[small], "__otros")

  cnt <- R$counts[genes_use, meta$cell]
  cnt <- cnt[Matrix::rowSums(cnt > 0) >= 5, ]
  sc  <- t(as.matrix(cnt))
  pc_arg <- intersect(c("psuedo.count", "pseudo.count"), names(formals(get.exp.stat)))[1]
  args <- list(sc.dat = sc, cell.type.labels = meta$celltype, cell.state.labels = st,
               cell.count.cutoff = 50, n.cores = N_CORES)
  args[[pc_arg]] <- 0.1
  stat   <- do.call(get.exp.stat, args)
  sc_sig <- select.marker(sc.dat = sc, stat = stat, pval.max = 0.01, lfc.min = 0.1)
  msg("   [", tag, "] ", nrow(sc), " celulas | ", ncol(sc_sig), " marcadores | ",
      nrow(mixture), " pseudobulk(s)")
  prism <- new.prism(reference = sc_sig, mixture = mixture, input.type = "count.matrix",
                     cell.type.labels = meta$celltype, cell.state.labels = st,
                     key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1)
  bp <- run.prism(prism = prism, n.cores = N_CORES)
  theta <- get.fraction(bp = bp, which.theta = "final", state.or.type = "type")
  theta <- theta[, intersect(COMMON_TYPES, colnames(theta)), drop = FALSE]

  cnt_all <- R$counts[genes_use, meta$cell]
  phi <- sapply(split(meta$cell, factor(meta$celltype, levels = colnames(theta))),
                function(cc) Matrix::rowSums(cnt_all[, cc, drop = FALSE]))
  phi <- t(sweep(phi, 2, colSums(phi), "/"))
  list(theta = theta, phi = phi)
}

recon_cor <- function(theta, phi, mixture) {
  g <- intersect(colnames(phi), colnames(mixture))
  obs <- mixture[rownames(theta), g, drop = FALSE]
  obs <- obs / rowSums(obs) * 1e6
  keep <- colMeans(obs) >= 1
  hat <- theta %*% phi[colnames(theta), g, drop = FALSE]
  hat <- hat / rowSums(hat) * 1e6
  sapply(rownames(theta), function(s)
    cor(log2(obs[s, keep] + 1), log2(hat[s, keep] + 1)))
}

RDS_A <- file.path(REF_OUT, "refA_GSE137846.rds")
RDS_B <- file.path(REF_OUT, "refB_Yang2026.rds")
RDS_Y <- file.path(REF_OUT, "Yang2026_seurat.rds")
HACER_REF_A <- FORZAR_REF_A || !file.exists(RDS_A)
HACER_YANG  <- FORZAR_YANG  || !file.exists(RDS_B) || !file.exists(RDS_Y)
if (!HACER_REF_A) msg("Seccion 2: referencia sana ya preparada, se salta")
if (!HACER_YANG)  msg("Secciones 3-6: Yang ya preparado, se salta")

if (HACER_REF_A) {

# SECCIÓN 2: REFERENCIA SANA GSE137846

msg("2) Leyendo GSE137846 ...")
f_mtx   <- find_file(REF_A_DIR, "counts\\.mtx")
f_genes <- find_file(REF_A_DIR, "gene_names")
f_ann   <- find_file(REF_A_DIR, "sample_annotations")
msg("   matriz:   ", basename(f_mtx))
msg("   genes:    ", basename(f_genes))
msg("   anotacion:", basename(f_ann))

refA_counts <- as(readMM(f_mtx), "CsparseMatrix")
genesA <- fread(f_genes, header = FALSE)[[1]]
annA   <- fread(f_ann)

if (length(genesA) == nrow(refA_counts) + 1 || length(genesA) == ncol(refA_counts) + 1)
  genesA <- genesA[-1]

if (nrow(refA_counts) != length(genesA) && ncol(refA_counts) == length(genesA)) {
  msg("   transponiendo matriz (venia celulas x genes)")
  refA_counts <- t(refA_counts)
}
stopifnot(nrow(refA_counts) == length(genesA))
if (nrow(annA) != ncol(refA_counts))
  stop("La anotacion tiene ", nrow(annA), " filas pero la matriz tiene ",
       ncol(refA_counts), " celulas. Revisar sample_annotations.")

cat("\nColumnas de sample_annotations (GSE137846):\n"); print(head(annA, 3))

if (is.null(REF_A_CELLTYPE_COL)) {
  cand <- names(annA)[sapply(annA, function(v)
    (is.character(v) || is.factor(v)) && data.table::uniqueN(v) %in% 3:80)]
  score <- sapply(cand, function(cl) mean(!is.na(harmonize_label(unique(annA[[cl]])))))
  if (length(score) == 0 || max(score) == 0)
    stop("No pude detectar la columna de tipo celular. Definir REF_A_CELLTYPE_COL en 00_config.R")
  REF_A_CELLTYPE_COL <- names(which.max(score))
}
msg("   columna de tipo celular usada: '", REF_A_CELLTYPE_COL, "'")

labA_orig <- as.character(annA[[REF_A_CELLTYPE_COL]])
labA      <- harmonize_label(labA_orig)
mapA <- as.data.frame(table(original = labA_orig, armonizado = labA, useNA = "ifany"))
mapA <- mapA[mapA$Freq > 0, ]
mapA <- mapA[order(mapA$armonizado, -mapA$Freq), ]
fwrite(mapA, file.path(REF_OUT, "GSE137846_mapeo_etiquetas.csv"))
cat("\nMapeo de etiquetas GSE137846 -> nomenclatura comun (REVISAR):\n"); print(mapA)

bc_col    <- grep("barcode|cell|^name$", names(annA), ignore.case = TRUE, value = TRUE)[1]
donor_col <- grep("donor|sample|patient|subject|individual", names(annA),
                  ignore.case = TRUE, value = TRUE)
donor_col <- setdiff(donor_col, c(bc_col, REF_A_CELLTYPE_COL))[1]
cellsA <- if (!is.na(bc_col)) make.unique(as.character(annA[[bc_col]])) else
  paste0("refA_cell", seq_len(ncol(refA_counts)))
colnames(refA_counts) <- paste0("A_", cellsA)
rownames(refA_counts) <- genesA
refA_counts <- collapse_rows(refA_counts)

metaA <- data.frame(
  cell      = colnames(refA_counts),
  celltype  = labA,
  label_orig = labA_orig,
  donor     = if (!is.na(donor_col)) as.character(annA[[donor_col]]) else "refA",
  stringsAsFactors = FALSE
)
keepA <- !is.na(metaA$celltype)
msg("   celulas sin tipo armonizable descartadas: ", sum(!keepA))
refA_counts <- refA_counts[, keepA]
metaA <- metaA[keepA, ]
if (any(abs(refA_counts@x - round(refA_counts@x)) > 1e-6))
  warning("GSE137846 no parece contener cuentas enteras. BayesPrism necesita cuentas crudas.")
cat("\nCelulas por tipo en GSE137846:\n"); print(table(metaA$celltype))

saveRDS(list(counts = refA_counts, meta = metaA), RDS_A)
msg("   referencia sana guardada: ", basename(RDS_A))
}

# SECCIÓN 3: YANG 2026 LECTURA, CONTROL DE CALIDAD Y DOBLETES

if (HACER_YANG) {
msg("3) Leyendo las 20 muestras de ", NEW_DIR)
folders <- list.dirs(NEW_DIR, recursive = FALSE, full.names = FALSE)
faltan  <- setdiff(SAMPLE_META$folder, folders)
sobran  <- setdiff(folders, SAMPLE_META$folder)
if (length(faltan)) stop("Carpetas esperadas que no estan: ", paste(faltan, collapse = ", "))
if (length(sobran)) warning("Carpetas no contempladas en SAMPLE_META (se ignoran): ",
                            paste(sobran, collapse = ", "))

gene_map <- NULL
obj_list <- list()
qc_tab   <- list()
for (s in SAMPLE_META$folder) {
  d <- file.path(NEW_DIR, s)
  m <- Read10X(d, gene.column = 1)
  if (is.list(m)) m <- m[["Gene Expression"]]
  if (is.null(gene_map)) {
    ff <- list.files(d, pattern = "^(features|genes)\\.tsv", full.names = TRUE)[1]
    gene_map <- fread(ff, header = FALSE)[, 1:2]
    setnames(gene_map, c("ensembl", "symbol"))
  }
  sym <- gene_map$symbol[match(rownames(m), gene_map$ensembl)]
  sym[is.na(sym)] <- rownames(m)[is.na(sym)]
  rownames(m) <- sym
  m <- collapse_rows(m)
  colnames(m) <- paste0(s, "_", colnames(m))

  o <- CreateSeuratObject(m, project = s, min.cells = 0, min.features = 0)
  o$sample <- s
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  n0 <- ncol(o)
  o <- subset(o, subset = nFeature_RNA >= QC_MIN_FEATURES &
                          nCount_RNA > QC_MIN_COUNTS &
                          percent.mt < QC_MAX_MT)
  n1 <- ncol(o)

  n_dbl <- NA
  if (DO_DOUBLETS && requireNamespace("scDblFinder", quietly = TRUE)) {
    sce <- scDblFinder::scDblFinder(
      SingleCellExperiment::SingleCellExperiment(
        list(counts = GetAssayData(o, layer = "counts"))))
    o$doublet <- sce$scDblFinder.class
    n_dbl <- sum(o$doublet == "doublet")
    o <- subset(o, subset = doublet == "singlet")
  } else if (DO_DOUBLETS && s == SAMPLE_META$folder[1]) {
    warning("scDblFinder no instalado: se omite la remocion de dobletes. ",
            "Instalar con BiocManager::install('scDblFinder').")
  }
  n2 <- ncol(o)
  if (ncol(o) > MAX_CELLS_PER_SAMPLE) o <- o[, sample(colnames(o), MAX_CELLS_PER_SAMPLE)]
  qc_tab[[s]] <- data.frame(sample = s, barcodes = n0, post_QC = n1,
                            dobletes = n_dbl, singletes = n2, usadas = ncol(o))
  obj_list[[s]] <- o
  msg("   ", s, ": ", n0, " -> ", n1, " (QC) -> ", n2, " (dobletes) -> ", ncol(o))
}
qc_tab <- do.call(rbind, qc_tab)
fwrite(merge(SAMPLE_META, qc_tab, by.x = "folder", by.y = "sample"),
       file.path(REF_OUT, "Yang_QC_por_muestra.csv"))
fwrite(gene_map, file.path(REF_OUT, "Yang_mapa_ensembl_simbolo.csv"))

yang <- merge(obj_list[[1]], obj_list[-1])
rm(obj_list); invisible(gc())
yang <- JoinLayers(yang)
md <- SAMPLE_META[match(yang$sample, SAMPLE_META$folder), ]
for (cl in c("donor", "source", "group", "eye", "sex_expected"))
  yang[[cl]] <- md[[cl]]

# SECCIÓN 4: YANG 2026 VERIFICACIÓN DE SEXO

msg("4) Verificando sexo por muestra (XIST vs genes de chrY) ...")
cnt <- GetAssayData(yang, layer = "counts")
ygenes <- intersect(c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y"), rownames(cnt))
lib <- colSums(cnt)
sex_df <- do.call(rbind, lapply(split(colnames(yang), yang$sample), function(cc) {
  data.frame(XIST_cpm = sum(cnt["XIST", cc]) / sum(lib[cc]) * 1e6,
             chrY_cpm = sum(cnt[ygenes, cc]) / sum(lib[cc]) * 1e6)
}))
sex_df$sample <- rownames(sex_df)
sex_df$sex_inferred <- ifelse(sex_df$chrY_cpm > sex_df$XIST_cpm, "M", "F")
sex_df <- merge(SAMPLE_META[, c("folder", "donor", "source", "group", "sex_expected")],
                sex_df, by.x = "folder", by.y = "sample")
sex_df$coincide <- ifelse(is.na(sex_df$sex_expected), "sin dato",
                          ifelse(sex_df$sex_expected == sex_df$sex_inferred, "OK", "NO"))
fwrite(sex_df, file.path(REF_OUT, "Yang_verificacion_sexo.csv"))
cat("\nVerificacion de sexo:\n"); print(sex_df)
if (any(sex_df$coincide == "NO"))
  warning("Hay muestras cuyo sexo inferido no coincide con la Tabla 1. ",
          "Revisar la asignacion de carpetas en SAMPLE_META antes de seguir.")
live_unknown <- sex_df$sex_inferred[is.na(sex_df$sex_expected)]
if (!identical(sort(live_unknown), c("F", "F", "M")))
  warning("Los 3 vivos YnonD deberian ser 1 M + 2 F (D3 M, D4 F, D5 F) y salen: ",
          paste(live_unknown, collapse = ","))
rm(cnt); invisible(gc())

# SECCIÓN 5: YANG 2026 NORMALIZACIÓN, INTEGRACIÓN Y CLUSTERING

msg("5) Normalizando y agrupando ...")
yang <- NormalizeData(yang, verbose = FALSE)
yang <- FindVariableFeatures(yang, nfeatures = 3000, verbose = FALSE)
yang <- ScaleData(yang, verbose = FALSE)
yang <- RunPCA(yang, npcs = 40, verbose = FALSE)
red <- "pca"
if (requireNamespace("harmony", quietly = TRUE)) {
  yang <- harmony::RunHarmony(yang, group.by.vars = "sample", verbose = FALSE)
  red <- "harmony"
} else warning("harmony no instalado: se agrupa sin correccion de batch.")
yang <- FindNeighbors(yang, reduction = red, dims = 1:30, verbose = FALSE)
yang <- FindClusters(yang, resolution = CLUSTER_RESOLUTION, verbose = FALSE)
yang <- RunUMAP(yang, reduction = red, dims = 1:30, verbose = FALSE)

# SECCIÓN 6: YANG 2026 ANOTACIÓN DE TIPOS CELULARES

msg("6) Anotando clusters con marcadores ...")
mk <- lapply(MARKERS, intersect, rownames(yang))
mk <- mk[lengths(mk) >= 2]
yang <- AddModuleScore(yang, features = mk, name = "MK_", seed = SEED)
sc_cols <- paste0("MK_", seq_along(mk))
score_cl <- aggregate(yang@meta.data[, sc_cols],
                      by = list(cluster = yang$seurat_clusters), FUN = mean)
names(score_cl)[-1] <- names(mk)
sc_mat <- as.matrix(score_cl[, -1]); rownames(sc_mat) <- score_cl$cluster
best   <- colnames(sc_mat)[apply(sc_mat, 1, which.max)]
top2   <- apply(sc_mat, 1, function(r) { s <- sort(r, decreasing = TRUE); s[1] - s[2] })
annot_cl <- data.frame(cluster = rownames(sc_mat), celltype = best,
                       score_max = apply(sc_mat, 1, max), margen_vs_2do = top2,
                       n_cells = as.integer(table(yang$seurat_clusters)[rownames(sc_mat)]))
annot_cl$dudoso <- annot_cl$margen_vs_2do < 0.1 | annot_cl$score_max < 0.1

override_f <- file.path(REF_OUT, "Yang_anotacion_manual.csv")
if (file.exists(override_f)) {
  ov <- fread(override_f, colClasses = "character")
  annot_cl$celltype[match(ov$cluster, annot_cl$cluster)] <- ov$celltype
  msg("   aplicada anotacion manual desde ", basename(override_f))
}
fwrite(annot_cl, file.path(REF_OUT, "Yang_anotacion_por_cluster.csv"))
cat("\nAnotacion por cluster (revisar los 'dudoso'):\n"); print(annot_cl)

yang$celltype <- annot_cl$celltype[match(yang$seurat_clusters, annot_cl$cluster)]
yang$celltype <- factor(yang$celltype, levels = intersect(CELLTYPE_ORDER, unique(yang$celltype)))

ggsave(bg = "white", file.path(FIG_OUT, "01_Yang_UMAP_tipos.png"),
       DimPlot(yang, group.by = "celltype", label = TRUE, raster = TRUE) + NoLegend(),
       width = 7, height = 6, dpi = 200)
ggsave(bg = "white", file.path(FIG_OUT, "01_Yang_UMAP_clusters.png"),
       DimPlot(yang, group.by = "seurat_clusters", label = TRUE, raster = TRUE) + NoLegend(),
       width = 7, height = 6, dpi = 200)
ggsave(bg = "white", file.path(FIG_OUT, "01_Yang_UMAP_origen_grupo.png"),
       DimPlot(yang, group.by = c("source", "group"), raster = TRUE),
       width = 13, height = 6, dpi = 200)
ggsave(bg = "white", file.path(FIG_OUT, "01_Yang_dotplot_marcadores_por_cluster.png"),
       DotPlot(yang, features = unique(unlist(mk)), group.by = "seurat_clusters") +
         RotatedAxis() + theme(axis.text.x = element_text(size = 6)),
       width = 18, height = 8, dpi = 200)

prop_tab <- as.data.frame(prop.table(table(yang$sample, yang$celltype), 1))
names(prop_tab) <- c("folder", "celltype", "fraction")
prop_tab <- merge(prop_tab, SAMPLE_META, by = "folder")
fwrite(prop_tab, file.path(REF_OUT, "Yang_proporciones_por_muestra.csv"))
ggsave(bg = "white", file.path(FIG_OUT, "01_Yang_proporciones_por_muestra.png"),
       ggplot(prop_tab, aes(folder, fraction, fill = celltype)) +
         geom_col() + facet_grid(~ group, scales = "free_x", space = "free_x") +
         theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
         labs(x = NULL, y = "Fraccion de celulas", fill = "Tipo"),
       width = 12, height = 6, dpi = 200)

metaB <- data.frame(cell = colnames(yang),
                    celltype = as.character(yang$celltype),
                    cluster  = as.character(yang$seurat_clusters),
                    yang@meta.data[, c("sample", "donor", "source", "group", "eye")],
                    stringsAsFactors = FALSE)
saveRDS(list(counts = GetAssayData(yang, layer = "counts"), meta = metaB), RDS_B)
saveRDS(yang, RDS_Y)
msg("   Yang guardado: ", basename(RDS_B), ", ", basename(RDS_Y))
}

# SECCIÓN 7: CONTRASTE POR TRANSFERENCIA DE ETIQUETAS

if (DO_LABEL_TRANSFER_CHECK && (HACER_REF_A || HACER_YANG)) {
  msg("7) Transferencia de etiquetas desde GSE137846 (control) ...")
  if (!exists("refA_counts")) { tmp <- readRDS(RDS_A); refA_counts <- tmp$counts; metaA <- tmp$meta; rm(tmp) }
  if (!exists("yang")) yang <- readRDS(RDS_Y)
  refA_seu <- CreateSeuratObject(refA_counts, meta.data = data.frame(metaA, row.names = metaA$cell))
  refA_seu <- NormalizeData(refA_seu, verbose = FALSE)
  refA_seu <- FindVariableFeatures(refA_seu, nfeatures = 3000, verbose = FALSE)
  refA_seu <- ScaleData(refA_seu, verbose = FALSE)
  refA_seu <- RunPCA(refA_seu, npcs = 30, verbose = FALSE)
  anch <- FindTransferAnchors(reference = refA_seu, query = yang, dims = 1:30,
                              reference.reduction = "pca", verbose = FALSE)
  pred <- TransferData(anch, refdata = refA_seu$celltype, dims = 1:30, verbose = FALSE)
  yang$celltype_transfer <- pred$predicted.id
  yang$celltype_transfer_score <- pred$prediction.score.max
  conc <- as.data.frame.matrix(table(marcadores = yang$celltype,
                                     transferencia = yang$celltype_transfer))
  fwrite(cbind(celltype = rownames(conc), conc),
         file.path(REF_OUT, "Yang_concordancia_marcadores_vs_transferencia.csv"))
  cat("\nConcordancia marcadores vs transferencia:\n"); print(conc)
  msg("   concordancia global: ",
      round(mean(as.character(yang$celltype) == yang$celltype_transfer) * 100, 1), "%")
  rm(refA_seu, anch); invisible(gc())
}

# SECCIÓN 8: CIERRE DE LA PREPARACIÓN

msg("8) Preparacion terminada. Conviene revisar en 01_referencias:")
msg("  - GSE137846_mapeo_etiquetas.csv (ninguna etiqueta sin armonizar)")
msg("  - Yang_verificacion_sexo.csv")
msg("  - Yang_anotacion_por_cluster.csv y el dotplot (clusters 'dudoso')")
if (exists("yang")) rm(yang)
if (exists("refA_counts")) rm(refA_counts, metaA)
invisible(gc())

# SECCIÓN 9: CARGA DE LAS REFERENCIAS

msg("9) Cargando referencias ...")
refA <- readRDS(file.path(REF_OUT, "refA_GSE137846.rds"))
refB <- readRDS(file.path(REF_OUT, "refB_Yang2026.rds"))
mB   <- refB$meta

# SECCIÓN 10: UNIVERSO COMÚN DE GENES

msg("10) Definiendo universo de genes ...")
all_genes <- intersect(rownames(refA$counts), rownames(refB$counts))
genes_use <- tryCatch({
  dummy <- matrix(1, nrow = 2, ncol = length(all_genes),
                  dimnames = list(c("c1", "c2"), all_genes))
  dummy <- cleanup.genes(input = dummy, input.type = "count.matrix", species = "hs",
                         gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
                         exp.cells = 1)
  dummy <- select.gene.type(dummy, gene.type = "protein_coding")
  colnames(dummy)
}, error = function(e) {
  warning("cleanup.genes fallo (", conditionMessage(e), "); se usa filtro por nombre.")
  all_genes[!grepl("^RP[SL]\\d|^MRP[SL]|^MT-|^MALAT1$|^XIST$", all_genes)]
})
STRESS_GENES <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "IER2",
                  "IER3", "ATF3", "DUSP1", "NR4A1", "NR4A2", "HSPA1A", "HSPA1B", "HSPA6",
                  "HSPB1", "DNAJB1", "HSP90AA1", "ZFP36", "KLF2", "KLF4", "KLF6", "CCN1", "CCN2")
genes_use <- setdiff(genes_use, STRESS_GENES)
msg("   genes comunes: ", length(all_genes), " -> tras filtros: ", length(genes_use))

# SECCIÓN 11: TIPOS CELULARES COMUNES

nA <- table(refA$meta$celltype)
nB <- table(mB$celltype)
nBn <- table(mB$celltype[mB$group == "nonDM"])
has <- function(tab) names(tab)[tab >= MIN_CELLS_PER_TYPE]
COMMON_TYPES <- intersect(CELLTYPE_ORDER, Reduce(intersect, list(has(nA), has(nB), has(nBn))))
tt <- data.frame(celltype = union(names(nA), names(nB)))
tt$GSE137846 <- as.integer(nA[tt$celltype]); tt$Yang <- as.integer(nB[tt$celltype])
tt$Yang_nonDM <- as.integer(nBn[tt$celltype]); tt[is.na(tt)] <- 0
tt$usado <- tt$celltype %in% COMMON_TYPES
fwrite(tt, file.path(TAB_OUT, "tipos_celulares_usados.csv"))
msg("11) Tipos comunes (", length(COMMON_TYPES), "): ", paste(COMMON_TYPES, collapse = ", "))
if (any(!tt$usado)) msg("   excluidos: ", paste(tt$celltype[!tt$usado], collapse = ", "))

# SECCIÓN 12: PSEUDOBULKS DE YANG

msg("12) Armando pseudobulks ...")
mPB <- mB[mB$celltype %in% COMMON_TYPES, ]
cntPB <- refB$counts[genes_use, mPB$cell]
samples <- SAMPLE_META$folder[SAMPLE_META$folder %in% mPB$sample]
pb <- sapply(samples, function(s) Matrix::rowSums(cntPB[, mPB$cell[mPB$sample == s], drop = FALSE]))
ct_f <- factor(mPB$celltype, levels = COMMON_TYPES)
truth_cells <- prop.table(table(factor(mPB$sample, levels = samples), ct_f), 1)
truth_umi   <- tapply(Matrix::colSums(cntPB), list(factor(mPB$sample, levels = samples), ct_f), sum)
truth_umi[is.na(truth_umi)] <- 0
truth_umi   <- prop.table(truth_umi, 1)
mix_all <- t(pb)
pb_info <- as.data.table(SAMPLE_META)[folder %in% samples]
pb_info[, n_cells := as.integer(table(mPB$sample)[folder])]
fwrite(pb_info, file.path(TAB_OUT, "pseudobulks_info.csv"))
fwrite(data.table(sample = samples, as.data.frame.matrix(truth_umi)),
       file.path(TAB_OUT, "verdad_fraccion_UMI.csv"))
fwrite(data.table(sample = samples, as.data.frame.matrix(truth_cells)),
       file.path(TAB_OUT, "verdad_fraccion_celulas.csv"))
rm(cntPB); invisible(gc())

# SECCIÓN 13: DECONVOLUCIÓN CON BAYESPRISM

est_list <- list(); rec_list <- list()

msg("13A) Referencia sana GSE137846 ...")
resA <- deconvolve(refA, refA$meta$cell, mix_all, "A_sana_GSE137846",
                   MAX_CELLS_PER_TYPE, state_from = "label")
est_list[["A_sana_GSE137846"]] <- resA$theta
rec_list[["A_sana_GSE137846"]] <- recon_cor(resA$theta, resA$phi, mix_all)

donor_of <- setNames(SAMPLE_META$donor, SAMPLE_META$folder)

if (RUN_Y_NONDM_LODO) {
  msg("13B) Yang nonDM (leave-one-donor-out) ...")
  th <- list(); rc <- list()
  nonDM_donors <- unique(SAMPLE_META$donor[SAMPLE_META$group == "nonDM"])
  s_out <- samples[!donor_of[samples] %in% nonDM_donors]
  r <- deconvolve(refB, mB$cell[mB$group == "nonDM"], mix_all[s_out, , drop = FALSE],
                  "Y_nonDM_para_DM_DR", MAX_CELLS_PER_TYPE_LODO)
  th[[1]] <- r$theta; rc[[1]] <- recon_cor(r$theta, r$phi, mix_all)
  for (d in nonDM_donors) {
    s_d <- samples[donor_of[samples] == d]
    r <- deconvolve(refB, mB$cell[mB$group == "nonDM" & mB$donor != d],
                    mix_all[s_d, , drop = FALSE], paste0("Y_nonDM_sin_", d),
                    MAX_CELLS_PER_TYPE_LODO)
    th[[d]] <- r$theta; rc[[d]] <- recon_cor(r$theta, r$phi, mix_all)
  }
  est_list[["Y_nonDM_LODO"]] <- do.call(rbind, th)
  rec_list[["Y_nonDM_LODO"]] <- unlist(unname(rc))
}

if (RUN_Y_TODAS_LODO) {
  msg("13C) Yang completa (leave-one-donor-out) ...")
  th <- list(); rc <- list()
  for (d in unique(SAMPLE_META$donor)) {
    s_d <- samples[donor_of[samples] == d]
    r <- deconvolve(refB, mB$cell[mB$donor != d], mix_all[s_d, , drop = FALSE],
                    paste0("Y_todas_sin_", d), MAX_CELLS_PER_TYPE_LODO)
    th[[d]] <- r$theta; rc[[d]] <- recon_cor(r$theta, r$phi, mix_all)
  }
  est_list[["Y_todas_LODO"]] <- do.call(rbind, th)
  rec_list[["Y_todas_LODO"]] <- unlist(unname(rc))
}
REF_LEVELS <- names(est_list)
for (nm in REF_LEVELS)
  fwrite(data.table(sample = rownames(est_list[[nm]]), est_list[[nm]]),
         file.path(DEC_OUT, paste0("proporciones_estimadas_", nm, ".csv")))

# SECCIÓN 14: ESTIMADO FRENTE A VERDADERO

msg("14) Evaluando ...")
long <- rbindlist(lapply(REF_LEVELS, function(nm) {
  th <- est_list[[nm]]
  data.table(reference = nm, sample = rep(rownames(th), ncol(th)),
             celltype = rep(colnames(th), each = nrow(th)), estimado = as.vector(th))
}))
long[, verdad_umi    := truth_umi[cbind(sample, celltype)]]
long[, verdad_celulas := truth_cells[cbind(sample, celltype)]]
long[, error := estimado - verdad_umi]
long <- merge(long, as.data.table(SAMPLE_META)[, .(sample = folder, donor, source, group)], by = "sample")
long[, reference := factor(reference, levels = REF_LEVELS)]
long[, celltype := factor(celltype, levels = COMMON_TYPES)]
fwrite(long, file.path(TAB_OUT, "estimado_vs_verdadero.csv"))

per_sample <- long[, .(RMSE  = sqrt(mean(error^2)),
                       MAE   = mean(abs(error)),
                       pearson_composicion = cor(estimado, verdad_umi),
                       RMSE_vs_celulas = sqrt(mean((estimado - verdad_celulas)^2))),
                   by = .(reference, sample, donor, source, group)]
rec_dt <- rbindlist(lapply(REF_LEVELS, function(nm)
  data.table(reference = nm, sample = names(rec_list[[nm]]), reconstruccion = rec_list[[nm]])))
per_sample <- merge(per_sample, rec_dt, by = c("reference", "sample"), all.x = TRUE)
fwrite(per_sample[order(reference, group, sample)], file.path(TAB_OUT, "metricas_por_pseudobulk.csv"))

resumen <- per_sample[, .(n_muestras = .N, n_donantes = uniqueN(donor),
                          RMSE_mediana = median(RMSE), MAE_mediana = median(MAE),
                          pearson_mediana = median(pearson_composicion),
                          reconstruccion_mediana = median(reconstruccion, na.rm = TRUE)),
                      by = .(reference, group)][order(reference, group)]
fwrite(resumen, file.path(TAB_OUT, "RESUMEN_error_por_grupo.csv"))
cat("\n==== RESUMEN: error de deconvolucion por grupo ====\n"); print(resumen)

sesgo <- long[, .(sesgo_medio = mean(error), sesgo_mediano = median(error),
                  verdad_media = mean(verdad_umi), estimado_medio = mean(estimado)),
              by = .(reference, celltype, group)][order(reference, celltype, group)]
fwrite(sesgo, file.path(TAB_OUT, "sesgo_por_tipo_y_grupo.csv"))

# SECCIÓN 15: MODELOS ESTADÍSTICOS

msg("15) Modelos ...")
per_sample[, group := factor(group, levels = c("nonDM", "DM", "DR"))]
per_sample[, reference := factor(reference, levels = REF_LEVELS)]
long[, group := factor(group, levels = c("nonDM", "DM", "DR"))]
has_lmer <- requireNamespace("lmerTest", quietly = TRUE)
has_emm  <- requireNamespace("emmeans", quietly = TRUE)

per_sample[, source := factor(source, levels = c("Cadaver", "Live"))]
diseno <- dcast(unique(per_sample[, .(donor, group, source)]), group ~ source,
                fun.aggregate = length, value.var = "donor")
contrastes <- list()

sink(file.path(TAB_OUT, "MODELOS_error_vs_grupo.txt"))
cat("Unidad: pseudobulk (muestra). Error = RMSE entre fraccion estimada y\n",
    "fraccion verdadera (UMIs) sobre todos los tipos celulares.\n",
    "Modelo principal (bloques incompletos por origen):\n",
    "   error ~ group + source + (1 | donor)\n",
    "   - source = bloque (vivo / cadaver), efecto fijo\n",
    "   - (1 | donor) = los dos ojos de un donante no son independientes\n\n",
    "Donantes por grupo y origen (bloques incompletos: DR sin vivos):\n", sep = "")
print(diseno)
cat("\n")
for (nm in REF_LEVELS) {
  d <- per_sample[reference == nm]
  cat("=====================================================\n")
  cat("Referencia:", nm, "\n")
  cat("=====================================================\n")
  for (resp in c("RMSE", "reconstruccion")) {
    if (has_lmer) {
      cat("\n--", resp, "~ group + source + (1 | donor)   [bloqueado por origen] --\n")
      fit_b <- tryCatch(lmerTest::lmer(as.formula(paste(resp, "~ group + source + (1 | donor)")),
                                       data = d), error = function(e) NULL)
      if (!is.null(fit_b)) {
        print(anova(fit_b))
        cf <- summary(fit_b)$coefficients
        if ("sourceLive" %in% rownames(cf))
          cat("Efecto del origen (vivo - cadaver):", signif(cf["sourceLive", "Estimate"], 3),
              " p =", signif(cf["sourceLive", ncol(cf)], 3), "\n")
        if (has_emm) {
          pr <- pairs(emmeans::emmeans(fit_b, ~ group), adjust = "holm")
          print(pr)
          contrastes[[paste(nm, resp)]] <- data.table(reference = nm, respuesta = resp,
                                                      as.data.frame(summary(pr)))
        }
      }
      cat("\n--", resp, "~ group + (1 | donor)   [sin bloquear, comparacion] --\n")
      fit_u <- tryCatch(lmerTest::lmer(as.formula(paste(resp, "~ group + (1 | donor)")),
                                       data = d), error = function(e) NULL)
      if (!is.null(fit_u)) {
        print(anova(fit_u))
        if (has_emm) print(pairs(emmeans::emmeans(fit_u, ~ group), adjust = "holm"))
      }
    }
    dd <- d[, .(v = mean(get(resp))), by = .(donor, group)]
    cat("Kruskal-Wallis a nivel donante (sin bloquear): p =",
        signif(kruskal.test(v ~ group, data = dd)$p.value, 3), "\n")
  }
  dc <- d[source == "Cadaver"]
  cat("\n-- Sensibilidad solo cadaveres (n donantes por grupo:",
      paste(dc[, uniqueN(donor), by = group]$V1, collapse = "/"), ") --\n")
  ddc <- dc[, .(v = mean(RMSE)), by = .(donor, group)]
  if (uniqueN(ddc$group) > 1)
    cat("Kruskal-Wallis RMSE a nivel donante: p =",
        signif(kruskal.test(v ~ group, data = ddc)$p.value, 3), "\n")
  print(dc[, .(RMSE_mediana = median(RMSE)), by = group])

  cat("\n-- Sesgo por tipo celular (estimado - verdadero) ~ group + source + (1 | donor) --\n")
  tb <- long[reference == nm, {
    dt <- data.table(error = error, group = group, donor = donor,
                     source = factor(source, levels = c("Cadaver", "Live")))
    p_b <- NA_real_
    if (has_lmer) {
      fb <- tryCatch(suppressMessages(lmerTest::lmer(error ~ group + source + (1 | donor), data = dt)),
                     error = function(e) NULL)
      if (!is.null(fb)) p_b <- anova(fb)["group", "Pr(>F)"]
    }
    dd <- dt[, .(e = mean(error)), by = .(donor, group)]
    .(sesgo_nonDM = mean(dt$error[dt$group == "nonDM"]),
      sesgo_DM    = mean(dt$error[dt$group == "DM"]),
      sesgo_DR    = mean(dt$error[dt$group == "DR"]),
      p_grupo_bloqueado = p_b,
      p_kruskal_donante = tryCatch(kruskal.test(e ~ group, data = dd)$p.value,
                                   error = function(err) NA_real_))
  }, by = celltype]
  tb[, padj_bloqueado := p.adjust(p_grupo_bloqueado, "BH")]
  tb[, padj_kruskal := p.adjust(p_kruskal_donante, "BH")]
  print(tb)
  fwrite(tb, file.path(TAB_OUT, paste0("sesgo_por_tipo_modelo_", nm, ".csv")))
  cat("\n")
}

if (has_lmer && length(REF_LEVELS) > 1) {
  cat("=====================================================\n")
  cat("Comparacion entre referencias (bloqueado por origen):\n")
  cat("RMSE ~ reference * group + reference * source + (1|donor) + (1|sample)\n")
  cat("(el efecto del origen puede ser distinto en cada referencia)\n")
  cat("=====================================================\n")
  fit_all <- tryCatch(lmerTest::lmer(RMSE ~ reference * group + reference * source +
                                       (1 | donor) + (1 | sample), data = per_sample),
                      error = function(e) NULL)
  if (!is.null(fit_all)) {
    print(anova(fit_all))
    if (has_emm) {
      cat("\nAumento del error en DM y DR respecto de nonDM, por referencia:\n")
      print(pairs(emmeans::emmeans(fit_all, ~ group | reference), reverse = TRUE, adjust = "holm"))
      cat("\nDiferencias entre referencias dentro de cada grupo:\n")
      print(pairs(emmeans::emmeans(fit_all, ~ reference | group), adjust = "holm"))
    }
  }
}
sink()
if (length(contrastes))
  fwrite(rbindlist(contrastes, fill = TRUE),
         file.path(TAB_OUT, "CONTRASTES_grupo_bloqueado_por_origen.csv"))

# SECCIÓN 16: FIGURAS

msg("16) Figuras ...")
gcol <- c(nonDM = "#4C9A6A", DM = "#D9A23B", DR = "#C2483D")

g1 <- ggplot(long[reference == "A_sana_GSE137846"],
             aes(verdad_umi, estimado, color = group, shape = source)) +
  geom_abline(linetype = 2, color = "grey50") + geom_point(size = 2.2) +
  facet_wrap(~ celltype, scales = "free") + scale_color_manual(values = gcol) +
  theme_bw() + labs(x = "Fraccion verdadera", y = "Fraccion estimada",
                    title = "Referencia sana (GSE137846) sobre pseudobulks de Yang 2026",
                    color = "Grupo", shape = "Origen")
ggsave(bg = "white", file.path(FIG_OUT, "02_A_estimado_vs_verdadero.png"), g1, width = 12, height = 9, dpi = 200)

g2 <- ggplot(per_sample, aes(group, RMSE, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(aes(shape = source), position = position_jitter(width = 0.12, seed = 1), size = 2) +
  facet_wrap(~ reference, nrow = 1) + scale_fill_manual(values = gcol) + theme_bw() +
  labs(x = NULL, y = "RMSE (estimado vs verdadero)", fill = "Grupo", shape = "Origen",
       title = "Error de deconvolucion segun el estado de la retina")
ggsave(bg = "white", file.path(FIG_OUT, "02_RMSE_por_grupo.png"), g2, width = 11, height = 5, dpi = 200)

g3 <- ggplot(sesgo, aes(group, celltype, fill = sesgo_medio)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.3f", sesgo_medio)), size = 2.6) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C2483D", midpoint = 0) +
  facet_wrap(~ reference, nrow = 1) + theme_bw() +
  labs(x = NULL, y = NULL, fill = "Sesgo medio\n(est - verd)",
       title = "Sobre (+) / subestimacion (-) de cada tipo celular")
ggsave(bg = "white", file.path(FIG_OUT, "02_sesgo_por_tipo.png"), g3, width = 11, height = 6, dpi = 200)

g4 <- ggplot(per_sample, aes(group, reconstruccion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(aes(shape = source), position = position_jitter(width = 0.12, seed = 1), size = 2) +
  facet_wrap(~ reference, nrow = 1) + scale_fill_manual(values = gcol) + theme_bw() +
  labs(x = NULL, y = "Correlacion observado vs reconstruido (log2 CPM)", fill = "Grupo",
       shape = "Origen", title = "Que tan bien explica cada referencia la expresion del pseudobulk")
ggsave(bg = "white", file.path(FIG_OUT, "02_reconstruccion_por_grupo.png"), g4, width = 11, height = 5, dpi = 200)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
msg("Listo. Empezar por tablas/RESUMEN_error_por_grupo.csv y tablas/MODELOS_error_vs_grupo.txt")
