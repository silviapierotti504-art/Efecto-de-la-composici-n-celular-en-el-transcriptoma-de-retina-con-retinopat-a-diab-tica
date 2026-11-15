# SECCIÓN 0: CONFIGURACIÓN

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(limma)
  library(edgeR)
})

BASE_DIR   <- "Resultados/Sc vs Sc"
REFB_RDS   <- file.path(BASE_DIR, "resultados_bayesprism", "01_referencias", "refB_Yang2026.rds")
GENES_FILE <- file.path(BASE_DIR, "genes_estres_oxidativo.csv")
OUT_DIR    <- file.path(BASE_DIR, "resultados_DE_por_tipo")

MIN_CELULAS_PSEUDOBULK <- 20
MIN_MUESTRAS_GRUPO     <- 2
TIPOS <- c("Rod", "Cone", "BC", "AC", "HC", "RGC", "MGC", "Microglia")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "DE_completo"), showWarnings = FALSE)
msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

# SECCIÓN 1: DATOS

msg("1) Cargando datos ...")
refB  <- readRDS(REFB_RDS)
meta  <- refB$meta
genes_ox <- unique(fread(GENES_FILE)[[1]])
msg("   genes de estres oxidativo en la lista: ", length(genes_ox),
    " | presentes en el single-cell: ", sum(genes_ox %in% rownames(refB$counts)))

# SECCIÓN 2: PSEUDOBULK POR MUESTRA Y TIPO CELULAR

msg("2) Armando pseudobulks por muestra y tipo celular ...")
pb_list <- list(); info_list <- list()
for (ct in TIPOS) {
  m <- meta[meta$celltype == ct, ]
  n_por_muestra <- table(m$sample)
  ok <- names(n_por_muestra)[n_por_muestra >= MIN_CELULAS_PSEUDOBULK]
  if (length(ok) == 0) next
  cnt <- sapply(ok, function(s) Matrix::rowSums(refB$counts[, m$cell[m$sample == s], drop = FALSE]))
  smeta <- unique(m[m$sample %in% ok, c("sample", "donor", "source", "group")])
  smeta <- smeta[match(ok, smeta$sample), ]
  smeta$n_celulas <- as.integer(n_por_muestra[ok])
  pb_list[[ct]] <- cnt
  info_list[[ct]] <- data.table(celltype = ct, smeta)
}
info <- rbindlist(info_list)
fwrite(info, file.path(OUT_DIR, "pseudobulks_por_tipo.csv"))
print(dcast(info, celltype ~ group, value.var = "sample", fun.aggregate = length))

# SECCIÓN 3: LIMMA-VOOM POR TIPO CELULAR

msg("3) Expresion diferencial por tipo celular ...")
CONTRASTES <- c(nonDM_vs_DM = "groupnonDM", DR_vs_DM = "groupDR")
res_all <- list(); fry_all <- list(); resumen <- list()

for (ct in names(pb_list)) {
  sm <- info[celltype == ct]
  sm[, group := factor(group, levels = c("DM", "nonDM", "DR"))]
  n_grupo <- table(sm$group)
  if (any(n_grupo < MIN_MUESTRAS_GRUPO)) {
    msg("   ", ct, ": se omite (muestras por grupo: ",
        paste(names(n_grupo), n_grupo, sep = "=", collapse = ", "), ")")
    next
  }
  usa_origen <- uniqueN(sm$source) > 1
  design <- if (usa_origen) model.matrix(~ group + source, data = sm) else
    model.matrix(~ group, data = sm)

  dge  <- DGEList(pb_list[[ct]])
  keep <- filterByExpr(dge, design = design)
  dge  <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])

  v  <- voom(dge, design)
  cf <- suppressWarnings(duplicateCorrelation(v, design, block = sm$donor)$consensus)
  if (!is.finite(cf)) cf <- 0
  v  <- voom(dge, design, block = sm$donor, correlation = cf)
  cf <- suppressWarnings(duplicateCorrelation(v, design, block = sm$donor)$consensus)
  if (!is.finite(cf)) cf <- 0
  fit <- eBayes(lmFit(v, design, block = sm$donor, correlation = cf))

  idx_ox <- which(rownames(v) %in% genes_ox)
  for (cn in names(CONTRASTES)) {
    coef <- CONTRASTES[[cn]]
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
    tt <- data.table(gen = rownames(tt), tt)
    fwrite(tt, file.path(OUT_DIR, "DE_completo", paste0(ct, "_", cn, ".csv")))

    ox <- tt[gen %in% genes_ox]
    ox[, FDR_en_conjunto := p.adjust(P.Value, "BH")]
    res_all[[paste(ct, cn)]] <- data.table(celltype = ct, contraste = cn,
                                           ox[, .(gen, logFC, AveExpr, t, P.Value,
                                                  adj.P.Val_global = adj.P.Val, FDR_en_conjunto)])
    fr <- fry(v, index = list(estres_oxidativo = idx_ox), design = design,
              contrast = which(colnames(design) == coef),
              block = sm$donor, correlation = cf)
    fry_all[[paste(ct, cn)]] <- data.table(celltype = ct, contraste = cn, fr)

    resumen[[paste(ct, cn)]] <- data.table(
      celltype = ct, contraste = cn, muestras = nrow(sm), donantes = uniqueN(sm$donor),
      origen_en_modelo = usa_origen, correlacion_donante = round(cf, 3),
      genes_testeados = nrow(tt), genes_FDR05 = sum(tt$adj.P.Val < 0.05),
      genes_ox_testeados = nrow(ox), genes_ox_FDR05 = sum(ox$FDR_en_conjunto < 0.05),
      ox_suben = sum(ox$FDR_en_conjunto < 0.05 & ox$logFC > 0),
      ox_bajan = sum(ox$FDR_en_conjunto < 0.05 & ox$logFC < 0),
      conjunto_direccion = fr$Direction, conjunto_p = fr$PValue,
      conjunto_p_mixto = fr$PValue.Mixed)
  }
  msg("   ", ct, ": ", nrow(sm), " muestras, ", nrow(v), " genes, correlacion donante = ", round(cf, 2))
}

# SECCIÓN 4: TABLAS DE RESULTADOS

msg("4) Guardando tablas ...")
ox_res <- rbindlist(res_all)
fry_res <- rbindlist(fry_all)
resumen <- rbindlist(resumen)
resumen[, conjunto_FDR := p.adjust(conjunto_p, "BH"), by = contraste]
fwrite(ox_res[order(contraste, celltype, P.Value)],
       file.path(OUT_DIR, "estres_oxidativo_por_gen_y_tipo.csv"))
fwrite(resumen, file.path(OUT_DIR, "RESUMEN_por_tipo_y_contraste.csv"))
cat("\n==== RESUMEN ====\n"); print(resumen[, .(celltype, contraste, muestras, genes_FDR05,
                                                genes_ox_FDR05, ox_suben, ox_bajan,
                                                conjunto_direccion, conjunto_p, conjunto_FDR)])

# SECCIÓN 5: FIGURAS

msg("5) Figuras ...")
orden <- intersect(TIPOS, unique(resumen$celltype))

resumen[, celltype := factor(celltype, levels = rev(orden))]
resumen[, signo := ifelse(conjunto_direccion == "Up", 1, -1) * -log10(conjunto_p)]
g1 <- ggplot(resumen, aes(contraste, celltype, fill = signo)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("p=%.3g", conjunto_p)), size = 3) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C2483D", midpoint = 0,
                       name = "-log10(p)\ncon signo\n(+ sube)") +
  theme_bw() + labs(x = NULL, y = NULL,
                    title = "Conjunto de estres oxidativo por tipo celular (fry)",
                    subtitle = "Grupo de referencia: DM")
ggsave(file.path(OUT_DIR, "estres_oxidativo_conjunto_por_tipo.png"), g1,
       width = 7, height = 5, dpi = 200, bg = "white")

sig_genes <- unique(ox_res[FDR_en_conjunto < 0.05, gen])
if (length(sig_genes) == 0) sig_genes <- unique(ox_res[order(P.Value)][1:30, gen])
hm <- ox_res[gen %in% sig_genes]
hm[, celltype := factor(celltype, levels = orden)]
hm[, marca := ifelse(FDR_en_conjunto < 0.05, "*", "")]
g2 <- ggplot(hm, aes(celltype, gen, fill = logFC)) +
  geom_tile(color = "white") + geom_text(aes(label = marca), size = 4) +
  facet_wrap(~ contraste) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C2483D", midpoint = 0) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     axis.text.y = element_text(size = 6)) +
  labs(x = NULL, y = NULL, title = "Genes de estres oxidativo por tipo celular",
       subtitle = "* = FDR < 0.05 dentro del conjunto")
ggsave(file.path(OUT_DIR, "estres_oxidativo_genes_por_tipo.png"), g2,
       width = 9, height = max(4, 0.15 * length(sig_genes) + 2), dpi = 200, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
msg("Listo. Empezar por RESUMEN_por_tipo_y_contraste.csv")
