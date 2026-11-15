# 0. CONFIGURACIÓN

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(future)
})

plan("sequential")
set.seed(42)

ruta <- "C:/Users/silvi/Desktop/Codigo_deconvolucion/01_atlas_referencia/GSE137846"

dir_base <- "C:/Users/silvi/Desktop/Codigo_deconvolucion"

dir_out <- file.path(dir_base, "01_atlas_referencia")
dir_fig <- file.path(dir_out, "figuras")
dir_tab <- file.path(dir_out, "tablas")
for (d in c(dir_out, dir_fig, dir_tab)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

col_tipo    <- "Labels"
col_region  <- "tissue"
col_donante <- NULL
region_macula    <- "MR"
region_periferia <- "PR"

umbrales <- list(min_genes = 200, max_genes = 6000, max_umis = 40000, max_mt = 10)

min_celulas_estado <- 30

log_file <- file.path(dir_out, "registro_metodos.txt")
cat("REGISTRO DE ANÁLISIS - ATLAS DE REFERENCIA GSE137846\n",
    "Fecha: ", format(Sys.time()), "\n\n", file = log_file, sep = "")

registrar <- function(...) {
  txt <- paste0(...)
  message(txt)
  cat(txt, "\n", file = log_file, append = TRUE, sep = "")
}
guardar_tabla <- function(x, nombre, row.names = TRUE) {
  write.csv(x, file.path(dir_tab, nombre), row.names = row.names)
}
guardar_fig <- function(p, nombre, w = 8, h = 6) {
  ggsave(file.path(dir_fig, nombre), p, width = w, height = h, dpi = 300, bg = "white")
}
get_counts <- function(o) {
  tryCatch(GetAssayData(o, assay = "RNA", layer = "counts"),
           error = function(e) GetAssayData(o, assay = "RNA", slot = "counts"))
}
pseudobulk <- function(cnt, grupos) {
  g <- factor(grupos)
  M <- sparse.model.matrix(~ 0 + g)
  colnames(M) <- levels(g)
  as.matrix(cnt %*% M)
}

registrar("Parámetros de QC: nFeature > ", umbrales$min_genes, ", nFeature < ", umbrales$max_genes,
          ", nCount < ", umbrales$max_umis, ", %mitocondrial < ", umbrales$max_mt)
registrar("Mínimo de células por región para definir estados regionales: ", min_celulas_estado)


# 1. CARGA DE DATOS CRUDOS E INVENTARIO DE METADATOS
counts <- readMM(file.path(ruta, "GSE137846_Seq-Well_counts.mtx"))
counts <- as(counts, "CsparseMatrix")
genes  <- read.table(file.path(ruta, "GSE137846_Seq-Well_gene_names.txt"),
                     header = FALSE, stringsAsFactors = FALSE)[, 1]
meta   <- read.table(file.path(ruta, "GSE137846_Seq-Well_sample_annotations.txt"),
                     header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

stopifnot(nrow(counts) == length(genes), ncol(counts) == nrow(meta))

registrar("\n== 1. DATOS CRUDOS ==")
registrar("Genes duplicados en la lista de genes: ", sum(duplicated(genes)),
          " (se hacen únicos con make.unique)")
rownames(counts) <- make.unique(genes)
colnames(counts) <- rownames(meta)
registrar("Matriz cruda: ", nrow(counts), " genes x ", ncol(counts), " células")
registrar("Columnas de metadatos: ", paste(colnames(meta), collapse = ", "))

if (is.null(col_donante)) {
  cand <- intersect(c("donor", "Donor", "donor_id", "patient", "individual",
                      "sample", "Sample", "orig.ident", "batch"), colnames(meta))
  col_donante <- if (length(cand) > 0) cand[1] else NA
}
registrar("Columna usada como donante/muestra: ", col_donante,
          " (revisar que efectivamente identifica donantes)")

inventario <- lapply(colnames(meta), function(cn) {
  v <- meta[[cn]]
  if (length(unique(v)) <= 50) {
    data.frame(columna = cn, nivel = names(table(v, useNA = "ifany")),
               n_celulas = as.vector(table(v, useNA = "ifany")))
  }
})
guardar_tabla(bind_rows(inventario), "01_inventario_metadatos.csv", row.names = FALSE)

sin_tipo <- sum(is.na(meta[[col_tipo]]) | meta[[col_tipo]] == "")
registrar("Células sin anotación de tipo celular en los metadatos originales: ", sin_tipo)

tab_pre <- table(tipo = meta[[col_tipo]], region = meta[[col_region]], useNA = "ifany")
guardar_tabla(as.data.frame.matrix(tab_pre), "01_composicion_preQC.csv")


# 2. CONTROL DE CALIDAD
obj <- CreateSeuratObject(counts = counts, meta.data = meta,
                          project = "GSE137846", min.cells = 3, min.features = 100)
registrar("\n== 2. CONTROL DE CALIDAD ==")
registrar("Tras CreateSeuratObject (min.cells = 3, min.features = 100): ",
          nrow(obj), " genes x ", ncol(obj), " células")
rm(counts); invisible(gc())

n_mt <- length(grep("^MT-", rownames(obj)))
registrar("Genes mitocondriales detectados (patrón ^MT-): ", n_mt,
          ifelse(n_mt == 0, "  <-- ATENCIÓN: revisar nomenclatura de genes", ""))
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^HB[AB]")

md <- obj@meta.data
obj$pasa_QC <- md$nFeature_RNA > umbrales$min_genes &
               md$nFeature_RNA < umbrales$max_genes &
               md$nCount_RNA   < umbrales$max_umis &
               md$percent.mt   < umbrales$max_mt

qc_df <- obj@meta.data
lineas <- list(nFeature_RNA = c(umbrales$min_genes, umbrales$max_genes),
               nCount_RNA = umbrales$max_umis, percent.mt = umbrales$max_mt,
               percent.ribo = NULL)
p_qc <- lapply(names(lineas), function(v) {
  p <- ggplot(qc_df, aes(x = .data[[col_tipo]], y = .data[[v]], fill = .data[[col_region]])) +
    geom_violin(scale = "width", linewidth = 0.2) +
    labs(x = NULL, y = v, fill = "Región") + theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(lineas[[v]])) p <- p + geom_hline(yintercept = lineas[[v]], linetype = "dashed", colour = "red")
  if (v == "nCount_RNA") p <- p + scale_y_log10()
  p
})
guardar_fig(wrap_plots(p_qc, ncol = 2, guides = "collect"), "02_QC_por_tipo_y_region.png", w = 12, h = 8)

retencion <- obj@meta.data %>%
  group_by(tipo = .data[[col_tipo]], region = .data[[col_region]]) %>%
  summarise(n_antes = n(), n_despues = sum(pasa_QC),
            pct_retenido = round(100 * mean(pasa_QC), 1), .groups = "drop")
guardar_tabla(retencion, "02_retencion_QC_por_tipo_region.csv", row.names = FALSE)
registrar("Células que pasan el QC: ", sum(obj$pasa_QC), " de ", ncol(obj),
          " (", round(100 * mean(obj$pasa_QC), 1), "%)")

obj <- subset(obj, subset = pasa_QC)
obj <- obj[, !is.na(obj@meta.data[[col_tipo]]) & obj@meta.data[[col_tipo]] != ""]
registrar("Referencia final tras QC y anotación: ", nrow(obj), " genes x ", ncol(obj), " células")


# 3. COMPOSICIÓN CELULAR DE LA REFERENCIA
registrar("\n== 3. COMPOSICIÓN DE LA REFERENCIA ==")
md <- obj@meta.data

tab_tr <- table(tipo = md[[col_tipo]], region = md[[col_region]])
guardar_tabla(as.data.frame.matrix(tab_tr), "03_celulas_tipo_x_region.csv")
registrar("Células por región: ", paste(names(table(md[[col_region]])),
                                        table(md[[col_region]]), sep = " = ", collapse = "; "))

if (!is.na(col_donante)) {
  guardar_tabla(as.data.frame.matrix(table(tipo = md[[col_tipo]], donante = md[[col_donante]])),
                "03_celulas_tipo_x_donante.csv")
  guardar_tabla(as.data.frame.matrix(table(donante = md[[col_donante]], region = md[[col_region]])),
                "03_celulas_donante_x_region.csv")
  registrar("Número de donantes/muestras: ", length(unique(md[[col_donante]])))
}

resumen_tipo <- md %>%
  group_by(tipo = .data[[col_tipo]]) %>%
  summarise(n_celulas      = n(),
            n_macula       = sum(.data[[col_region]] == region_macula),
            n_periferia    = sum(.data[[col_region]] == region_periferia),
            pct_macula     = round(100 * mean(.data[[col_region]] == region_macula), 1),
            mediana_UMIs   = median(nCount_RNA),
            mediana_genes  = median(nFeature_RNA),
            mediana_pct_mt = round(median(percent.mt), 2),
            .groups = "drop") %>%
  arrange(desc(n_celulas))
if (!is.na(col_donante)) {
  resumen_tipo <- resumen_tipo %>%
    left_join(md %>% group_by(tipo = .data[[col_tipo]]) %>%
                summarise(n_donantes = n_distinct(.data[[col_donante]])), by = "tipo")
}
guardar_tabla(resumen_tipo, "03_resumen_por_tipo_celular.csv", row.names = FALSE)

df_comp <- as.data.frame(tab_tr)
p_comp <- ggplot(df_comp, aes(x = reorder(tipo, -Freq), y = Freq, fill = region)) +
  geom_col(position = "stack") +
  labs(x = NULL, y = "Número de células", fill = "Región",
       title = "Composición de la referencia por tipo celular y región") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_frac <- ggplot(df_comp, aes(x = reorder(tipo, -Freq), y = Freq, fill = region)) +
  geom_col(position = "fill") + geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(x = NULL, y = "Fracción de células", fill = "Región",
       title = "Origen regional de cada tipo celular") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
guardar_fig(p_comp / p_frac, "03_composicion_tipo_region.png", w = 9, h = 9)


# 4. NORMALIZACIÓN, GENES VARIABLES, PCA Y UMAP
registrar("\n== 4. NORMALIZACIÓN Y REDUCCIÓN DE DIMENSIONES ==")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
guardar_fig(ElbowPlot(obj, ndims = 50) + ggtitle("Varianza explicada por componente"),
            "04_elbowplot.png", w = 6, h = 4)

n_pcs <- 30
obj <- RunUMAP(obj, dims = 1:n_pcs, seed.use = 42, verbose = FALSE)
registrar("LogNormalize (factor 1e4); 2000 genes variables (vst); PCA 50 componentes; ",
          "UMAP y vecinos con ", n_pcs, " componentes")

p_u1 <- DimPlot(obj, group.by = col_tipo, label = TRUE, repel = TRUE) + ggtitle("Tipo celular (autores)")
p_u2 <- DimPlot(obj, group.by = col_region) + ggtitle("Región")
p_umaps <- p_u1 + p_u2
if (!is.na(col_donante)) {
  p_umaps <- p_umaps + DimPlot(obj, group.by = col_donante) + ggtitle("Donante / muestra")
}
guardar_fig(p_umaps, "04_UMAP_tipo_region_donante.png", w = 18, h = 6)


# 5. CLUSTERING NO SUPERVISADO Y CONCORDANCIA CON LA ANOTACIÓN DE LOS AUTORES
registrar("\n== 5. CLUSTERING Y CONCORDANCIA ==")
obj <- FindNeighbors(obj, dims = 1:n_pcs, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 42, verbose = FALSE)
registrar("Clustering de Louvain (resolución 0.5): ", length(levels(obj$seurat_clusters)), " clusters")

conc <- table(cluster = obj$seurat_clusters, tipo = obj@meta.data[[col_tipo]])
guardar_tabla(as.data.frame.matrix(conc), "05_concordancia_clusters_vs_tipos.csv")

df_conc <- as.data.frame(prop.table(conc, 1))
p_conc <- ggplot(df_conc, aes(x = tipo, y = cluster, fill = Freq)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(Freq > 0.05, round(Freq, 2), "")), size = 2.5) +
  scale_fill_gradient(low = "white", high = "darkred", name = "Fracción\ndel cluster") +
  labs(x = "Tipo celular (autores)", y = "Cluster no supervisado") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
guardar_fig(p_conc, "05_concordancia_clusters_tipos.png", w = 8, h = 7)

if (requireNamespace("mclust", quietly = TRUE)) {
  ari <- mclust::adjustedRandIndex(obj$seurat_clusters, obj@meta.data[[col_tipo]])
  registrar("Índice de Rand ajustado (clusters vs anotación de autores): ", round(ari, 3))
}
pureza <- apply(prop.table(conc, 1), 1, max)
registrar("Pureza mediana de los clusters: ", round(median(pureza), 3),
          "; clusters con pureza < 0.8: ", paste(names(pureza)[pureza < 0.8], collapse = ", "))


# 6. MARCADORES CANÓNICOS POR TIPO CELULAR
canonicos <- list(
  Bastones   = c("RHO", "NRL", "NR2E3", "SAG", "GNAT1"),
  Conos      = c("ARR3", "GNAT2", "PDE6H", "OPN1SW", "OPN1MW"),
  Bipolares  = c("VSX2", "GRM6", "TRPM1", "CABP5", "PRKCA"),
  Amacrinas  = c("GAD1", "GAD2", "SLC6A9", "TFAP2A", "TFAP2B"),
  Horizont   = c("ONECUT1", "ONECUT2", "LHX1"),
  RGCs       = c("RBPMS", "SNCG", "NEFL", "THY1", "POU4F2"),
  Muller     = c("RLBP1", "GLUL", "SLC1A3", "SOX9", "AQP4"),
  Astrocitos = c("GFAP", "PAX2", "S100B"),
  Microglia  = c("C1QA", "C1QB", "P2RY12", "CX3CR1", "AIF1"),
  Endotelio  = c("PECAM1", "CLDN5", "VWF", "KDR"),
  Pericitos  = c("PDGFRB", "RGS5", "ACTA2"),
  RPE        = c("RPE65", "BEST1", "TTR")
)
canonicos <- lapply(canonicos, function(g) g[g %in% rownames(obj)])
canonicos <- canonicos[lengths(canonicos) > 0]

p_dot <- DotPlot(obj, features = canonicos, group.by = col_tipo) +
  RotatedAxis() + theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Marcadores canónicos por tipo celular (anotación de autores)")
guardar_fig(p_dot, "06_dotplot_marcadores_canonicos.png", w = 16, h = 5)


# 7. MARCADORES DE NOVO POR TIPO CELULAR
registrar("\n== 7. MARCADORES DE NOVO ==")
Idents(obj) <- col_tipo
marc <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                       test.use = "wilcox", verbose = FALSE)
marc_sig <- marc %>% filter(p_val_adj < 0.05)
guardar_tabla(marc_sig, "07_marcadores_por_tipo_celular.csv", row.names = FALSE)
registrar("Marcadores significativos (p ajustado < 0.05, log2FC > 0.25, min.pct 0.25) por tipo: ",
          paste(names(table(marc_sig$cluster)), table(marc_sig$cluster), sep = " = ", collapse = "; "))

top5 <- marc_sig %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5, with_ties = FALSE)
obj_hm <- subset(obj, downsample = 100)
obj_hm <- ScaleData(obj_hm, features = unique(top5$gene), verbose = FALSE)
p_hm <- DoHeatmap(obj_hm, features = unique(top5$gene), group.by = col_tipo, size = 3) +
  theme(axis.text.y = element_text(size = 6))
guardar_fig(p_hm, "07_heatmap_top5_marcadores.png", w = 11, h = 9)
rm(obj_hm); invisible(gc())


# 8. HETEROGENEIDAD REGIONAL DENTRO DE CADA TIPO CELULAR
registrar("\n== 8. HETEROGENEIDAD REGIONAL ==")
n_tr <- table(obj@meta.data[[col_tipo]], obj@meta.data[[col_region]])
tipos_div <- rownames(n_tr)[apply(n_tr, 1, min) >= min_celulas_estado]
registrar("Tipos con >= ", min_celulas_estado, " células en ambas regiones: ",
          paste(tipos_div, collapse = ", "))

res_reg <- list()
for (t in tipos_div) {
  sub <- subset(obj, idents = t)
  de <- FindMarkers(sub, ident.1 = region_macula, ident.2 = region_periferia,
                    group.by = col_region, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  de$gene <- rownames(de); de$tipo <- t
  res_reg[[t]] <- de
}
de_reg <- bind_rows(res_reg)
guardar_tabla(de_reg, "08_DE_macula_vs_periferia_por_tipo.csv", row.names = FALSE)

resumen_reg <- de_reg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  group_by(tipo) %>%
  summarise(n_DE = n(), mayor_en_macula = sum(avg_log2FC > 0),
            mayor_en_periferia = sum(avg_log2FC < 0), .groups = "drop")
guardar_tabla(resumen_reg, "08_resumen_DE_regional.csv", row.names = FALSE)
registrar("Genes con diferencia regional (p ajustado < 0.05, |log2FC| > 0.5) por tipo: ",
          paste(resumen_reg$tipo, resumen_reg$n_DE, sep = " = ", collapse = "; "))


# 9. SUBCLUSTERING DE LA MACROGLÍA
registrar("\n== 9. SUBCLUSTERING DE MACROGLÍA ==")
macro <- subset(obj, idents = "Macroglia")
macro <- NormalizeData(macro, verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:15, verbose = FALSE) %>%
  FindClusters(resolution = 0.3, random.seed = 42, verbose = FALSE) %>%
  RunUMAP(dims = 1:15, seed.use = 42, verbose = FALSE)

tab_macro <- table(subcluster = macro$seurat_clusters, region = macro@meta.data[[col_region]])
guardar_tabla(as.data.frame.matrix(tab_macro), "09_macroglia_subcluster_x_region.csv")

genes_glia <- c("RLBP1", "GLUL", "SLC1A3", "AQP4", "SOX9",
                "GFAP", "PAX2", "S100B",
                "VIM", "CRYAB")
genes_glia <- genes_glia[genes_glia %in% rownames(macro)]
p_macro <- (DimPlot(macro, label = TRUE) + ggtitle("Subclusters de macroglía")) +
  (DimPlot(macro, group.by = col_region) + ggtitle("Región")) +
  (DotPlot(macro, features = genes_glia) + RotatedAxis() + ggtitle("Müller vs astrocitos"))
guardar_fig(p_macro, "09_macroglia_subclustering.png", w = 18, h = 5)

Idents(macro) <- "seurat_clusters"
marc_macro <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25,
                             logfc.threshold = 0.25, verbose = FALSE) %>%
  filter(p_val_adj < 0.05)
guardar_tabla(marc_macro, "09_marcadores_subclusters_macroglia.csv", row.names = FALSE)
registrar("Subclusters de macroglía: ", length(levels(macro$seurat_clusters)),
          ". Ver tabla 09 para composición regional y marcadores.")


# 10. PERFILES DE REFERENCIA Y COLINEALIDAD
registrar("\n== 10. PERFILES DE REFERENCIA Y COLINEALIDAD ==")

obj$estado <- ifelse(obj@meta.data[[col_tipo]] %in% tipos_div,
                     paste(obj@meta.data[[col_tipo]], obj@meta.data[[col_region]], sep = "_"),
                     as.character(obj@meta.data[[col_tipo]]))
guardar_tabla(as.data.frame(table(estado = obj$estado)), "10_celulas_por_estado.csv", row.names = FALSE)

cnt <- get_counts(obj)
pb  <- pseudobulk(cnt, obj$estado)
lcpm <- log2(t(t(pb) / colSums(pb)) * 1e6 + 1)

genes_inf <- marc_sig %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50, with_ties = FALSE) %>% pull(gene) %>% unique()
cor_est <- cor(lcpm[genes_inf, ], method = "pearson")
guardar_tabla(round(cor_est, 3), "10_correlacion_perfiles_estados.csv")

orden <- hclust(as.dist(1 - cor_est))$order
df_cor <- as.data.frame(as.table(cor_est))
df_cor$Var1 <- factor(df_cor$Var1, levels = colnames(cor_est)[orden])
df_cor$Var2 <- factor(df_cor$Var2, levels = colnames(cor_est)[orden])
p_cor <- ggplot(df_cor, aes(Var1, Var2, fill = Freq)) + geom_tile() +
  geom_text(aes(label = round(Freq, 2)), size = 2) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred", midpoint = 0.5,
                       name = "Pearson") +
  labs(x = NULL, y = NULL,
       title = "Similitud entre perfiles de referencia (log-CPM, genes marcadores)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
guardar_fig(p_cor, "10_correlacion_perfiles.png", w = 10, h = 9)

tipo_base <- sub(paste0("_(", region_macula, "|", region_periferia, ")$"), "", colnames(cor_est))
pares <- which(upper.tri(cor_est), arr.ind = TRUE)
pares <- data.frame(perfil_1 = colnames(cor_est)[pares[, 1]],
                    perfil_2 = colnames(cor_est)[pares[, 2]],
                    pearson  = cor_est[pares],
                    mismo_tipo = tipo_base[pares[, 1]] == tipo_base[pares[, 2]]) %>%
  arrange(desc(pearson))
guardar_tabla(pares, "10_pares_perfiles_ordenados.csv", row.names = FALSE)
top_dist <- head(filter(pares, !mismo_tipo), 5)
registrar("Pares de tipos distintos más similares: ",
          paste0(top_dist$perfil_1, "-", top_dist$perfil_2, " (", round(top_dist$pearson, 2), ")",
                 collapse = "; "))


# 11. GENES DE miRNA (MIR*)
registrar("\n== 11. GENES MIR* ==")
mir <- grep("^MIR", rownames(obj), value = TRUE)
if (length(mir) > 0) {
  det <- Matrix::rowSums(cnt[mir, , drop = FALSE] > 0)
  df_mir <- data.frame(gene = mir, n_celulas_detectado = det,
                       es_gen_huesped = grepl("HG", mir)) %>% arrange(desc(n_celulas_detectado))
  guardar_tabla(df_mir, "11_genes_MIR_detectados.csv", row.names = FALSE)
  registrar("Genes MIR* en la matriz: ", length(mir), " (", sum(df_mir$es_gen_huesped),
            " genes huésped); detectados en >= 10 células: ", sum(det >= 10))
} else {
  registrar("No hay genes con prefijo MIR en la matriz filtrada.")
}


# 12. GUARDADO
dir_rds <- file.path(dir_out, "objetos")
dir.create(dir_rds, showWarnings = FALSE, recursive = TRUE)
saveRDS(obj,   file.path(dir_rds, "GSE137846_referencia_QC_anotada.rds"))
saveRDS(macro, file.path(dir_rds, "GSE137846_macroglia_subclusters.rds"))
writeLines(capture.output(sessionInfo()), file.path(dir_out, "sessionInfo.txt"))
registrar("\nObjeto final guardado: GSE137846_referencia_QC_anotada.rds (",
          ncol(obj), " células; columnas 'estado' y 'seurat_clusters' añadidas)")
registrar("Versiones: R ", R.version$major, ".", R.version$minor,
          "; Seurat ", as.character(packageVersion("Seurat")))
