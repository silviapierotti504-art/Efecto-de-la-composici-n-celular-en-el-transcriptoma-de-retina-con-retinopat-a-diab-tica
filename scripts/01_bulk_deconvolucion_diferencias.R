# SECCIÓN 1: LIBRERÍAS

library(BayesPrism)
library(Seurat)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(GEOquery)
library(limma)
library(edgeR)
library(readxl)
library(writexl)
options(timeout = 600)

# SECCIÓN 2: DIRECTORIOS Y CARGA DE ARCHIVOS

dir     <- "C:/Users/silvi/Desktop/Codigo_deconvolucion/"
dir_res <- paste0(dir, "GitHub/Resultados/")
dir_img <- paste0(dir, "GitHub/imagenes/")
dir.create(dir_res, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_img, showWarnings = FALSE, recursive = TRUE)

bulk <- as.data.frame(data.table::fread(paste0(dir, "GSE160306_human_retina_DR_totalRNA_counts.txt")))
sc   <- readRDS(paste0(dir, "01_atlas_referencia/objetos/GSE137846_referencia_QC_anotada.rds"))

# SECCIÓN 3: PREPROCESAMIENTO DEL BULK (ENSEMBL A SÍMBOLOS)

rownames(bulk) <- bulk[[1]]
bulk    <- as.matrix(bulk[, vapply(bulk, is.numeric, logical(1)), drop = FALSE])
ens     <- sub("\\.\\d+$", "", rownames(bulk))
simbolo <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL")
simbolo[is.na(simbolo)] <- mapIds(EnsDb.Hsapiens.v86, ens[is.na(simbolo)], "SYMBOL", "GENEID")
bulk <- rowsum(bulk[!is.na(simbolo), ], toupper(simbolo[!is.na(simbolo)]))

# SECCIÓN 4: PREPROCESAMIENTO DE LA REFERENCIA SINGLE-CELL

counts <- GetAssayData(sc, layer = "counts")
sc.dat <- t(as.matrix(rowsum(counts, toupper(rownames(counts)))))
ct     <- sc$Labels

sc.dat <- sc.dat[, !grepl("^(RPL|RPS|MRPL|MRPS|MT-)|^MALAT1$", colnames(sc.dat))]

crom <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86, "GENEID"),
                              keytype = "GENEID", columns = c("SYMBOL", "SEQNAME"))
sexuales <- toupper(unique(crom$SYMBOL[crom$SEQNAME %in% c("X", "Y")]))
sc.dat   <- sc.dat[, !colnames(sc.dat) %in% sexuales]

# SECCIÓN 5: DECONVOLUCIÓN CON BAYESPRISM

correr_bayesprism <- FALSE

if (correr_bayesprism) {
  prism <- new.prism(reference = sc.dat, mixture = t(bulk), input.type = "count.matrix",
                     cell.type.labels = ct, cell.state.labels = ct, key = NULL)
  res <- run.prism(prism, n.cores = 2)
  fr  <- get.fraction(res, "final", "type")
  write.table(data.frame(title = rownames(fr), fr, check.names = FALSE),
              paste0(dir, "fracciones.csv"), sep = ";", row.names = FALSE, quote = FALSE)
}

# SECCIÓN 6: METADATOS DEL BULK

archivo_meta <- paste0(dir, "metadatos_bulk.xlsx")

if (!file.exists(archivo_meta)) {
  archivo <- paste0(dir, "GSE160306_series_matrix.txt.gz")
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160306/matrix/GSE160306_series_matrix.txt.gz",
                archivo, mode = "wb")
  write_xlsx(pData(getGEO(filename = archivo)), archivo_meta)
}

# SECCIÓN 7: UNIÓN DE MUESTRAS Y VARIABLES DEL MODELO

meta <- as.data.frame(read_excel(archivo_meta))
frac <- read.csv(paste0(dir, "fracciones.csv"), sep = ";", row.names = 1, check.names = FALSE)

m      <- intersect(intersect(colnames(bulk), meta$title), rownames(frac))
bulk_m <- bulk[, m]
meta   <- meta[match(m, meta$title), ]
frac   <- frac[m, ]

grupo  <- relevel(factor(meta$`disease_group:ch1`), ref = "Diabetic")
tejido <- factor(meta$`sample_site:ch1`)
edad   <- as.numeric(meta$`age:ch1`)
sexo   <- factor(meta$`gender:ch1`)
rin    <- as.numeric(meta$`rin:ch1`)
donor  <- factor(meta$`donor:ch1`)

# SECCIÓN 8: COMPOSICIÓN CELULAR (PCA) Y MATRIZ DE DISEÑO

pcs    <- prcomp(asin(sqrt(frac)), scale. = TRUE)$x[, 1:3]
diseno <- model.matrix(~ grupo + tejido + edad + sexo + rin + pcs)

# SECCIÓN 9: MODELO A (CON COMPOSICIÓN CELULAR)

y <- DGEList(bulk_m)
y <- y[filterByExpr(y, diseno, min.count = 20), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

v    <- voom(y, diseno)
corr <- duplicateCorrelation(v, diseno, block = donor)$consensus
v    <- voom(y, diseno, block = donor, correlation = corr)
fit  <- eBayes(lmFit(v, diseno, block = donor, correlation = corr))

# SECCIÓN 10: MODELO B (SIN COMPOSICIÓN CELULAR)

diseno_B <- model.matrix(~ grupo + tejido + edad + sexo + rin)
v_B      <- voom(y, diseno_B)
corr_B   <- duplicateCorrelation(v_B, diseno_B, block = donor)$consensus
v_B      <- voom(y, diseno_B, block = donor, correlation = corr_B)
fit_B    <- eBayes(lmFit(v_B, diseno_B, block = donor, correlation = corr_B))

# SECCIÓN 11: PARTICIÓN DE VARIANZA: TEJIDO VS COMPOSICIÓN (R2)

X_base <- model.matrix(~ grupo + edad + sexo + rin)
X_tej  <- model.matrix(~ grupo + edad + sexo + rin + tejido)
X_pcs  <- model.matrix(~ grupo + edad + sexo + rin + pcs)

rss <- function(X) sapply(seq_len(nrow(v)), function(i) {
  w <- v$weights[i, ]
  sum(w * lm.wfit(X, v$E[i, ], w)$residuals^2)
})
rss_base <- rss(X_base)
rss_tej  <- rss(X_tej)
rss_pcs  <- rss(X_pcs)
rss_full <- rss(diseno)

R2_solo_tejido      <- (rss_pcs - rss_full) / rss_base
R2_solo_composicion <- (rss_tej - rss_full) / rss_base
R2_compartida       <- (rss_base - rss_full) / rss_base - R2_solo_tejido - R2_solo_composicion

# SECCIÓN 12: CONTRASTES DE INTERÉS Y EXPORTACIÓN

coefs_F  <- c("grupoControl", "grupoNPDR")
coef_tej <- grep("^tejido", colnames(diseno), value = TRUE)
coef_pcs <- grep("^pcs",    colnames(diseno), value = TRUE)

columnas_t <- function(f, cf, modelo) {
  r <- topTable(f, coef = cf, number = Inf, sort.by = "none", adjust.method = "BH")
  d <- data.frame(r$logFC, r$P.Value, r$adj.P.Val)
  names(d) <- paste0(modelo, c("_log2FC", "_p_valor", "_p_ajustado_FDR"))
  d
}
columnas_F <- function(f, modelo) {
  r <- topTable(f, coef = coefs_F, number = Inf, sort.by = "none", adjust.method = "BH")
  d <- data.frame(r$grupoControl, r$grupoNPDR, r$F, r$P.Value, r$adj.P.Val)
  names(d) <- paste0(modelo, c("_log2FC_Control_vs_Diabetico", "_log2FC_RDNP_vs_Diabetico",
                               "_F", "_p_valor", "_p_ajustado_FDR"))
  d
}
clasificar <- function(fdr_A, fdr_B) {
  ifelse(fdr_A < 0.05 & fdr_B < 0.05, "Significativo en ambos modelos",
  ifelse(fdr_B < 0.05,                "Solo en Modelo B (explicado por la composicion celular)",
  ifelse(fdr_A < 0.05,                "Solo en Modelo A",
                                      "No significativo")))
}
base <- data.frame(gen = rownames(fit), expresion_media = fit$Amean)

test_t <- function(cf) {
  d <- cbind(base, columnas_t(fit, cf, "ModeloA"), columnas_t(fit_B, cf, "ModeloB"))
  d[order(d$ModeloA_p_valor), ]
}

test_F <- cbind(base, columnas_F(fit, "ModeloA"), columnas_F(fit_B, "ModeloB"))
test_F$resultado <- clasificar(test_F$ModeloA_p_ajustado_FDR, test_F$ModeloB_p_ajustado_FDR)
test_F <- test_F[order(test_F$ModeloA_p_valor), ]

tab_tejido <- cbind(base, columnas_t(fit, coef_tej, "ModeloA"), columnas_t(fit_B, coef_tej, "ModeloB"),
                    R2_solo_tejido      = round(R2_solo_tejido, 4),
                    R2_solo_composicion = round(R2_solo_composicion, 4),
                    R2_compartida       = round(R2_compartida, 4))
tab_tejido$resultado <- clasificar(tab_tejido$ModeloA_p_ajustado_FDR, tab_tejido$ModeloB_p_ajustado_FDR)
tab_tejido <- tab_tejido[order(tab_tejido$ModeloB_p_valor), ]

r <- topTable(fit, coef = coef_pcs, number = Inf, adjust.method = "BH")
composicion <- data.frame(gen = rownames(r), r)

salida <- list(
  Test_t_Control_vs_Diabetico = test_t("grupoControl"),
  Test_t_RDNP_vs_Diabetico    = test_t("grupoNPDR"),
  Test_t_RDP_vs_Diabetico     = test_t("grupoPDR"),
  Test_F_enfermedad           = test_F,
  Tejido_macula_vs_periferia  = tab_tejido,
  Composicion_F               = composicion
)

n_sig     <- function(x) sum(x < 0.05)
explicado <- "Solo en Modelo B (explicado por la composicion celular)"

salida$Resumen <- data.frame(
  contraste = c("Control_vs_Diabetic", "NPDR_vs_Diabetic", "PDR_vs_Diabetic", "Enfermedad_F_global",
                "Enfermedad_total (modelo B)", "Enfermedad_explicada_por_composicion",
                "Tejido_total (modelo B)", "Tejido_indep_composicion", "Tejido_explicado_por_composicion",
                "Composicion_F"),
  FDR_0.05  = c(n_sig(salida$Test_t_Control_vs_Diabetico$ModeloA_p_ajustado_FDR),
                n_sig(salida$Test_t_RDNP_vs_Diabetico$ModeloA_p_ajustado_FDR),
                n_sig(salida$Test_t_RDP_vs_Diabetico$ModeloA_p_ajustado_FDR),
                n_sig(test_F$ModeloA_p_ajustado_FDR),
                n_sig(test_F$ModeloB_p_ajustado_FDR),
                sum(test_F$resultado == explicado),
                n_sig(tab_tejido$ModeloB_p_ajustado_FDR),
                n_sig(tab_tejido$ModeloA_p_ajustado_FDR),
                sum(tab_tejido$resultado == explicado),
                n_sig(composicion$adj.P.Val))
)
salida$Resumen_R2_tejido <- data.frame(
  componente = c("Solo tejido", "Solo composición", "Compartida"),
  mediana_R2 = round(c(median(R2_solo_tejido), median(R2_solo_composicion), median(R2_compartida)), 4)
)

write_xlsx(salida, paste0(dir_res, "DE_contrastes_interes.xlsx"))

# SECCIÓN 13: ESTRÉS OXIDATIVO (TEST F GEN POR GEN Y CONJUNTO)

oxi <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0006979", keytype = "GOALL", columns = "SYMBOL")$SYMBOL
oxi <- intersect(toupper(unique(oxi)), rownames(fit))

columnas_F_oxi <- function(f, modelo) {
  r <- topTable(f, coef = coefs_F, number = Inf, sort.by = "none")[oxi, ]
  d <- data.frame(r$grupoControl, r$grupoNPDR, r$F, r$P.Value, p.adjust(r$P.Value, method = "BH"))
  names(d) <- paste0(modelo, c("_log2FC_Control_vs_Diabetico", "_log2FC_RDNP_vs_Diabetico",
                               "_F", "_p_valor", "_p_ajustado_FDR"))
  d
}
genes_oxi <- cbind(data.frame(gen = oxi, expresion_media = fit$Amean[oxi]),
                   columnas_F_oxi(fit, "ModeloA"), columnas_F_oxi(fit_B, "ModeloB"))
genes_oxi$resultado <- clasificar(genes_oxi$ModeloA_p_ajustado_FDR, genes_oxi$ModeloB_p_ajustado_FDR)
genes_oxi <- genes_oxi[order(genes_oxi$ModeloA_p_valor), ]

indice    <- list(estres_oxidativo = match(oxi, rownames(v)))
direccion <- function(x) ifelse(x == "Up", "Mayor que en Diabetico", "Menor que en Diabetico")
nombre_ct <- c(grupoControl = "Control vs Diabetico", grupoNPDR = "RDNP vs Diabetico")

conjunto_oxi <- data.frame()
for (g in coefs_F) {
  a <- fry(v,   index = indice, design = diseno,   contrast = which(colnames(diseno) == g),
           block = donor, correlation = corr)
  b <- fry(v_B, index = indice, design = diseno_B, contrast = which(colnames(diseno_B) == g),
           block = donor, correlation = corr_B)
  conjunto_oxi <- rbind(conjunto_oxi, data.frame(
    comparacion       = nombre_ct[[g]],
    genes_en_conjunto = a$NGenes,
    ModeloA_direccion = direccion(a$Direction), ModeloA_p_valor = a$PValue,
    ModeloB_direccion = direccion(b$Direction), ModeloB_p_valor = b$PValue))
}

resumen_oxi <- data.frame(
  analisis = c("Genes de estres oxidativo analizados (GO:0006979)",
               "Genes con FDR < 0,05 en el test F, Modelo A",
               "Genes con FDR < 0,05 en el test F, Modelo B"),
  n_genes  = c(length(oxi), sum(genes_oxi$ModeloA_p_ajustado_FDR < 0.05),
               sum(genes_oxi$ModeloB_p_ajustado_FDR < 0.05))
)

write_xlsx(list(Resumen           = resumen_oxi,
                Test_F_genes      = genes_oxi,
                Conjunto_completo = conjunto_oxi),
           paste0(dir_res, "DE_estres_oxidativo.xlsx"))
