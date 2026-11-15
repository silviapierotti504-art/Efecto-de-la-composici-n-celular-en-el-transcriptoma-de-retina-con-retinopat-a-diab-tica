# SECCIÓN 1: LIBRERÍAS Y PARÁMETROS

library(writexl)

alfa_BH   <- 0.00106
potencia  <- 0.80
coefs_int <- c("grupoControl", "grupoNPDR")
df        <- fit$df.total

# SECCIÓN 2: MDE POR GEN Y PACIENTES NECESARIOS

por_gen <- list()
resumen <- data.frame()

for (cf in coefs_int) {

  se  <- fit$stdev.unscaled[, cf] * sqrt(fit$s2.post)
  mde <- (qt(1 - alfa_BH / 2, df) + qt(potencia, df)) * se

  por_gen[[cf]] <- data.frame(gen = rownames(fit), AveExpr = fit$Amean,
                              logFC_observado = fit$coefficients[, cf], SE = se,
                              MDE = mde, cambio_min_pct = round(100 * (2^mde - 1), 1))

  n_actual <- min(length(unique(donor[grupo == sub("grupo", "", cf)])),
                  length(unique(donor[grupo == "Diabetic"])))

  for (objetivo in c(0.2, 0.3, 0.5, 1)) {
    resumen <- rbind(resumen, data.frame(
      comparacion          = paste0(sub("grupo", "", cf), "_vs_Diabetic"),
      alfa_BH              = alfa_BH,
      pacientes_por_grupo  = n_actual,
      MDE_mediano          = round(median(mde), 3),
      logFC_objetivo       = objetivo,
      cambio_objetivo_pct  = round(100 * (2^objetivo - 1)),
      pacientes_necesarios = ceiling(n_actual * (median(mde) / objetivo)^2)))
  }
}
print(resumen)

# SECCIÓN 3: MDE SEGÚN NIVEL DE EXPRESIÓN

por_expresion <- data.frame()
for (cf in coefs_int) {
  d     <- por_gen[[cf]]
  nivel <- cut(rank(d$AveExpr), 3, labels = c("baja", "media", "alta"))
  por_expresion <- rbind(por_expresion, data.frame(
    comparacion    = cf,
    nivel          = levels(nivel),
    MDE_mediano    = round(tapply(d$MDE, nivel, median), 3),
    cambio_min_pct = round(tapply(d$cambio_min_pct, nivel, median), 1)))
}
print(por_expresion)

# SECCIÓN 4: EXPORTACIÓN DE RESULTADOS

write_xlsx(c(list(Resumen = resumen, Por_expresion = por_expresion), por_gen),
           paste0(dir_res, "Efecto_minimo_detectable.xlsx"))
