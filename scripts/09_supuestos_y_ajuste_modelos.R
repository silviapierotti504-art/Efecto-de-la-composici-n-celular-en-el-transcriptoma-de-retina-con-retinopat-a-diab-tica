# SECCIÓN 1: LIBRERÍAS

library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(limma)
library(edgeR)
library(readxl)
library(lmerTest)
library(writexl)

# SECCIÓN 2: DIRECTORIOS Y CARGA DE ARCHIVOS

dir     <- "C:/Users/silvi/Desktop/Codigo_deconvolucion/"
dir_res <- paste0(dir, "GitHub/Resultados/")
dir_img <- paste0(dir, "GitHub/imagenes/")

bulk <- as.data.frame(data.table::fread(paste0(dir, "GSE160306_human_retina_DR_totalRNA_counts.txt")))
meta <- as.data.frame(read_excel(paste0(dir, "metadatos_bulk.xlsx")))
frac <- read.csv(paste0(dir, "fracciones.csv"), sep = ";", row.names = 1, check.names = FALSE)

# SECCIÓN 3: PREPARACIÓN DEL BULK Y DE LAS VARIABLES

rownames(bulk) <- bulk[[1]]
bulk    <- as.matrix(bulk[, vapply(bulk, is.numeric, logical(1)), drop = FALSE])
ens     <- sub("\\.\\d+$", "", rownames(bulk))
simbolo <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL")
simbolo[is.na(simbolo)] <- mapIds(EnsDb.Hsapiens.v86, ens[is.na(simbolo)], "SYMBOL", "GENEID")
bulk <- rowsum(bulk[!is.na(simbolo), ], toupper(simbolo[!is.na(simbolo)]))

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

pcs      <- prcomp(asin(sqrt(frac)), scale. = TRUE)$x[, 1:3]
diseno   <- model.matrix(~ grupo + tejido + edad + sexo + rin + pcs)
diseno_B <- model.matrix(~ grupo + tejido + edad + sexo + rin)

# SECCIÓN 4: AJUSTE DE LOS MODELOS A Y B

y <- DGEList(bulk_m)
y <- y[filterByExpr(y, diseno, min.count = 20), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

v      <- voom(y, diseno)
corr   <- duplicateCorrelation(v, diseno, block = donor)$consensus
v      <- voom(y, diseno, block = donor, correlation = corr, save.plot = TRUE)
fit    <- eBayes(lmFit(v, diseno, block = donor, correlation = corr))

v_B    <- voom(y, diseno_B)
corr_B <- duplicateCorrelation(v_B, diseno_B, block = donor)$consensus
v_B    <- voom(y, diseno_B, block = donor, correlation = corr_B, save.plot = TRUE)
fit_B  <- eBayes(lmFit(v_B, diseno_B, block = donor, correlation = corr_B))

# SECCIÓN 5: RELACIÓN MEDIA-VARIANZA Y VARIANZA RESIDUAL

png(paste0(dir_img, "S1_media_varianza_modelos.png"), width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1))
plot(v$voom.xy, pch = 16, cex = 0.25, col = "grey40", xlab = "log2(recuentos + 0,5)", ylab = "Raíz de la desviación estándar")
lines(v$voom.line, col = "red", lwd = 2)
mtext("A", side = 3, adj = 0, font = 2)
plot(v_B$voom.xy, pch = 16, cex = 0.25, col = "grey40", xlab = "log2(recuentos + 0,5)", ylab = "Raíz de la desviación estándar")
lines(v_B$voom.line, col = "red", lwd = 2)
mtext("B", side = 3, adj = 0, font = 2)
plotSA(fit, main = "", xlab = "Expresión media (log2-CPM)", ylab = "Raíz de la desviación estándar residual")
mtext("C", side = 3, adj = 0, font = 2)
plotSA(fit_B, main = "", xlab = "Expresión media (log2-CPM)", ylab = "Raíz de la desviación estándar residual")
mtext("D", side = 3, adj = 0, font = 2)
dev.off()

# SECCIÓN 6: RESIDUOS ESTANDARIZADOS DE LOS MODELOS A Y B

res_est <- function(f, vv) {
  r <- residuals(f, vv) * sqrt(vv$weights) / f$sigma
  r
}
r_A <- res_est(fit, v)
r_B <- res_est(fit_B, v_B)
a_A <- fitted(fit)
a_B <- fitted(fit_B)

set.seed(2024)
muestra <- sample(length(r_A), 20000)

png(paste0(dir_img, "S2_residuos_modelos.png"), width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1))
qqnorm(r_A[muestra], pch = 16, cex = 0.3, main = "", xlab = "Cuantiles teóricos", ylab = "Residuos estandarizados")
qqline(r_A[muestra], col = "red", lwd = 2)
mtext("A", side = 3, adj = 0, font = 2)
qqnorm(r_B[muestra], pch = 16, cex = 0.3, main = "", xlab = "Cuantiles teóricos", ylab = "Residuos estandarizados")
qqline(r_B[muestra], col = "red", lwd = 2)
mtext("B", side = 3, adj = 0, font = 2)
smoothScatter(a_A[muestra], r_A[muestra], xlab = "Valor ajustado (log2-CPM)", ylab = "Residuos estandarizados")
abline(h = 0, col = "red", lwd = 2)
mtext("C", side = 3, adj = 0, font = 2)
smoothScatter(a_B[muestra], r_B[muestra], xlab = "Valor ajustado (log2-CPM)", ylab = "Residuos estandarizados")
abline(h = 0, col = "red", lwd = 2)
mtext("D", side = 3, adj = 0, font = 2)
dev.off()

# SECCIÓN 7: DISTRIBUCIÓN DE LOS P-VALORES

p_F_A <- topTable(fit,   coef = c("grupoControl", "grupoNPDR"), number = Inf, sort.by = "none")$P.Value
p_F_B <- topTable(fit_B, coef = c("grupoControl", "grupoNPDR"), number = Inf, sort.by = "none")$P.Value

png(paste0(dir_img, "S3_histograma_p_valores.png"), width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1))
hist(fit$p.value[, "grupoControl"], breaks = 40, col = "grey80", main = "", xlab = "p-valor, Control vs Diabético", ylab = "Número de genes")
hist(fit$p.value[, "grupoNPDR"],    breaks = 40, col = "grey80", main = "", xlab = "p-valor, RDNP vs Diabético",    ylab = "Número de genes")
hist(p_F_A, breaks = 40, col = "grey80", main = "", xlab = "p-valor, test F, Modelo A", ylab = "Número de genes")
hist(p_F_B, breaks = 40, col = "grey80", main = "", xlab = "p-valor, test F, Modelo B", ylab = "Número de genes")
dev.off()

# SECCIÓN 8: BONDAD DE AJUSTE POR GEN (R2)

r2_gen <- function(vv, X) {
  r2 <- numeric(nrow(vv))
  for (i in seq_len(nrow(vv))) {
    w   <- vv$weights[i, ]
    e   <- vv$E[i, ]
    scr <- sum(w * lm.wfit(X, e, w)$residuals^2)
    sct <- sum(w * (e - weighted.mean(e, w))^2)
    r2[i] <- 1 - scr / sct
  }
  r2
}
R2_A <- r2_gen(v, diseno)
R2_B <- r2_gen(v_B, diseno_B)
n    <- ncol(v)
R2a_A <- 1 - (1 - R2_A) * (n - 1) / (n - ncol(diseno))
R2a_B <- 1 - (1 - R2_B) * (n - 1) / (n - ncol(diseno_B))

aic <- selectModel(v, list(Modelo_A = diseno, Modelo_B = diseno_B), criterion = "aic")

por_gen <- data.frame(gen = rownames(fit), expresion_media = fit$Amean,
                      R2_ModeloA = round(R2_A, 4), R2_ModeloB = round(R2_B, 4),
                      R2_ajustado_ModeloA = round(R2a_A, 4), R2_ajustado_ModeloB = round(R2a_B, 4),
                      sigma_ModeloA = round(fit$sigma, 4), sigma_ModeloB = round(fit_B$sigma, 4),
                      modelo_preferido_AIC = aic$pref)

resumen <- data.frame(
  medida = c("Genes analizados", "Correlación intradonante (duplicateCorrelation)",
             "Grados de libertad residuales", "Grados de libertad a priori (eBayes)",
             "Varianza a priori (eBayes)", "Mediana de sigma", "Mediana de R2", "Mediana de R2 ajustado",
             "Genes con mayor R2 ajustado en este modelo", "Genes con menor AIC en este modelo"),
  Modelo_A = c(nrow(fit), round(corr, 3), fit$df.residual[1], round(median(fit$df.prior), 2),
               round(fit$s2.prior[1], 4), round(median(fit$sigma), 4), round(median(R2_A), 3),
               round(median(R2a_A), 3), sum(R2a_A > R2a_B), sum(aic$pref == "Modelo_A")),
  Modelo_B = c(nrow(fit_B), round(corr_B, 3), fit_B$df.residual[1], round(median(fit_B$df.prior), 2),
               round(fit_B$s2.prior[1], 4), round(median(fit_B$sigma), 4), round(median(R2_B), 3),
               round(median(R2a_B), 3), sum(R2a_B > R2a_A), sum(aic$pref == "Modelo_B")))

png(paste0(dir_img, "S4_R2_modelos.png"), width = 10, height = 4.5, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))
boxplot(list("Modelo A" = R2a_A, "Modelo B" = R2a_B), col = c("#C2185B", "grey70"), ylab = "R2 ajustado por gen")
mtext("A", side = 3, adj = 0, font = 2)
smoothScatter(R2a_B, R2a_A, xlab = "R2 ajustado, Modelo B", ylab = "R2 ajustado, Modelo A")
abline(0, 1, col = "red", lwd = 2)
mtext("B", side = 3, adj = 0, font = 2)
dev.off()

# SECCIÓN 9: MULTICOLINEALIDAD DE LAS COVARIABLES (VIF)

vif_tabla <- function(X) {
  X   <- X[, colnames(X) != "(Intercept)"]
  vif <- numeric(ncol(X))
  for (j in seq_len(ncol(X))) {
    r2     <- summary(lm(X[, j] ~ X[, -j]))$r.squared
    vif[j] <- 1 / (1 - r2)
  }
  data.frame(termino = colnames(X), VIF = round(vif, 2))
}
vif_A <- vif_tabla(diseno)
vif_B <- vif_tabla(diseno_B)
vif   <- merge(vif_A, vif_B, by = "termino", all.x = TRUE, suffixes = c("_ModeloA", "_ModeloB"))

# SECCIÓN 10: SUPUESTOS Y AJUSTE DE LOS MODELOS MIXTOS DE COMPOSICIÓN

mixtos <- data.frame()
res_m  <- list()
aj_m   <- list()
for (ct in colnames(frac)) {
  y_ct   <- asin(sqrt(frac[[ct]]))
  modelo <- lmer(y_ct ~ grupo + tejido + edad + sexo + rin + (1 | donor))
  vc     <- as.data.frame(VarCorr(modelo))
  var_d  <- vc$vcov[vc$grp == "donor"]
  var_e  <- vc$vcov[vc$grp == "Residual"]
  var_f  <- var(predict(modelo, re.form = NA))
  res_m[[ct]] <- residuals(modelo)
  aj_m[[ct]]  <- fitted(modelo)
  p_don  <- NA
  if (!isSingular(modelo)) p_don <- shapiro.test(ranef(modelo)$donor[, 1])$p.value
  mixtos <- rbind(mixtos, data.frame(
    tipo_celular        = ct,
    ajuste_singular     = isSingular(modelo),
    varianza_donante    = signif(var_d, 3),
    varianza_residual   = signif(var_e, 3),
    ICC_donante         = round(var_d / (var_d + var_e), 3),
    R2_marginal         = round(var_f / (var_f + var_d + var_e), 3),
    R2_condicional      = round((var_f + var_d) / (var_f + var_d + var_e), 3),
    p_Shapiro_residuos  = signif(shapiro.test(residuals(modelo))$p.value, 3),
    p_Shapiro_donante   = signif(p_don, 3)))
}

png(paste0(dir_img, "S5_residuos_modelos_mixtos.png"), width = 10, height = 10, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(4.5, 4.5, 2, 1))
for (ct in colnames(frac)) {
  qqnorm(res_m[[ct]], pch = 16, cex = 0.6, main = "", xlab = "Cuantiles teóricos", ylab = "Residuos")
  qqline(res_m[[ct]], col = "red", lwd = 2)
  mtext(ct, side = 3, adj = 0, font = 2)
}
dev.off()

png(paste0(dir_img, "S6_residuos_vs_ajustados_mixtos.png"), width = 10, height = 10, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(4.5, 4.5, 2, 1))
for (ct in colnames(frac)) {
  plot(aj_m[[ct]], res_m[[ct]], pch = 16, cex = 0.6, xlab = "Valor ajustado", ylab = "Residuos")
  abline(h = 0, col = "red", lwd = 2)
  mtext(ct, side = 3, adj = 0, font = 2)
}
dev.off()

# SECCIÓN 11: SUPUESTOS Y AJUSTE DEL MODELO DE miRNAs

mir <- as.data.frame(data.table::fread(paste0(dir, "GSE160308_human_retina_DR_smallRNA_counts.txt")))
rownames(mir) <- mir[[1]]
mir <- as.matrix(mir[, vapply(mir, is.numeric, logical(1)), drop = FALSE])

comunes <- intersect(m, colnames(mir))
idx     <- match(comunes, m)
mir     <- mir[, comunes]

dm <- data.frame(grupo = droplevels(grupo[idx]), tejido = tejido[idx],
                 edad = edad[idx], sexo = sexo[idx], rin = rin[idx])
dm$pcs <- pcs[idx, ]
diseno_mir <- model.matrix(~ grupo + tejido + edad + sexo + rin + pcs, data = dm)

ym <- DGEList(mir)
ym <- ym[filterByExpr(ym, diseno_mir, min.count = 20), , keep.lib.sizes = FALSE]
ym <- calcNormFactors(ym)

vm      <- voom(ym, diseno_mir)
corr_m  <- duplicateCorrelation(vm, diseno_mir, block = donor[idx])$consensus
vm      <- voom(ym, diseno_mir, block = donor[idx], correlation = corr_m, save.plot = TRUE)
fit_mir <- eBayes(lmFit(vm, diseno_mir, block = donor[idx], correlation = corr_m))

r_mir <- res_est(fit_mir, vm)
a_mir <- fitted(fit_mir)
p_F_mir <- topTable(fit_mir, coef = c("grupoControl", "grupoNPDR"), number = Inf, sort.by = "none")$P.Value
R2_mir  <- r2_gen(vm, diseno_mir)
R2a_mir <- 1 - (1 - R2_mir) * (ncol(vm) - 1) / (ncol(vm) - ncol(diseno_mir))

png(paste0(dir_img, "S7_supuestos_miRNAs.png"), width = 10, height = 12, units = "in", res = 300)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2, 1))
plot(vm$voom.xy, pch = 16, cex = 0.4, col = "grey40", xlab = "log2(recuentos + 0,5)", ylab = "Raíz de la desviación estándar")
lines(vm$voom.line, col = "red", lwd = 2)
mtext("A", side = 3, adj = 0, font = 2)
plotSA(fit_mir, main = "", xlab = "Expresión media (log2-CPM)", ylab = "Raíz de la desviación estándar residual")
mtext("B", side = 3, adj = 0, font = 2)
qqnorm(r_mir, pch = 16, cex = 0.3, main = "", xlab = "Cuantiles teóricos", ylab = "Residuos estandarizados")
qqline(r_mir, col = "red", lwd = 2)
mtext("C", side = 3, adj = 0, font = 2)
smoothScatter(a_mir, r_mir, xlab = "Valor ajustado (log2-CPM)", ylab = "Residuos estandarizados")
abline(h = 0, col = "red", lwd = 2)
mtext("D", side = 3, adj = 0, font = 2)
hist(fit_mir$p.value[, "grupoControl"], breaks = 40, col = "grey80", main = "", xlab = "p-valor, Control vs Diabético", ylab = "Número de miRNAs")
mtext("E", side = 3, adj = 0, font = 2)
hist(p_F_mir, breaks = 40, col = "grey80", main = "", xlab = "p-valor, test F", ylab = "Número de miRNAs")
mtext("F", side = 3, adj = 0, font = 2)
dev.off()

resumen_mir <- data.frame(
  medida = c("miRNAs analizados", "Muestras", "Correlación intradonante (duplicateCorrelation)",
             "Grados de libertad residuales", "Grados de libertad a priori (eBayes)",
             "Varianza a priori (eBayes)", "Mediana de sigma", "Mediana de R2", "Mediana de R2 ajustado"),
  modelo_miRNAs = c(nrow(fit_mir), ncol(vm), round(corr_m, 3), fit_mir$df.residual[1],
                    round(median(fit_mir$df.prior), 2), round(fit_mir$s2.prior[1], 4),
                    round(median(fit_mir$sigma), 4), round(median(R2_mir), 3), round(median(R2a_mir), 3)))

# SECCIÓN 12: EXPORTACIÓN DE RESULTADOS

write_xlsx(list(Resumen_limma = resumen, Ajuste_por_gen = por_gen, VIF = vif, Modelos_mixtos = mixtos,
                Resumen_miRNAs = resumen_mir),
           paste0(dir_res, "Supuestos_y_ajuste_modelos.xlsx"))
