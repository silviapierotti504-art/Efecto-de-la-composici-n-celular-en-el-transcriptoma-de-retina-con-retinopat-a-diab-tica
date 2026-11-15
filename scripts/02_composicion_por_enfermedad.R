# SECCIÓN 1: LIBRERÍAS

library(lmerTest)
library(writexl)

# SECCIÓN 2: MEDIAS DESCRIPTIVAS POR GRUPO Y POR TEJIDO

medias_grupo  <- aggregate(frac * 100, by = list(grupo = grupo),   FUN = mean)
medias_tejido <- aggregate(frac * 100, by = list(tejido = tejido), FUN = mean)

# SECCIÓN 3: MODELO MIXTO POR TIPO CELULAR

tabla <- data.frame()
for (ct in colnames(frac)) {
  y_ct   <- asin(sqrt(frac[[ct]]))
  modelo <- lmer(y_ct ~ grupo + tejido + edad + sexo + rin + (1 | donor))
  co     <- summary(modelo)$coefficients
  co     <- co[grep("^grupo|^tejido", rownames(co)), , drop = FALSE]
  tabla  <- rbind(tabla, data.frame(tipo_celular = ct,
                                    termino      = rownames(co),
                                    estimado     = co[, "Estimate"],
                                    p            = co[, "Pr(>|t|)"]))
}

# SECCIÓN 4: CORRECCIÓN POR COMPARACIONES MÚLTIPLES (BH)

tabla$FDR <- ave(tabla$p, tabla$termino, FUN = function(p) p.adjust(p, method = "BH"))
tabla     <- tabla[order(tabla$termino, tabla$p), ]

# SECCIÓN 5: GRÁFICO DE COMPOSICIÓN POR GRUPO

png(paste0(dir_img, "Composicion_por_grupo.png"), width = 10, height = 9, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
for (ct in colnames(frac)) {
  boxplot(frac[[ct]] * 100 ~ grupo, main = ct, ylab = "% estimado", xlab = "", las = 2)
}
dev.off()

# SECCIÓN 6: EXPORTACIÓN DE RESULTADOS

write_xlsx(list(Modelo_mixto = tabla, Medias_por_grupo = medias_grupo, Medias_por_tejido = medias_tejido),
           paste0(dir_res, "Composicion_por_grupo.xlsx"))
