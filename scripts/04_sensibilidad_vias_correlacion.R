# SECCIÓN 1: LIBRERÍAS Y MODELOS

library(limma)
library(writexl)

modelos <- list(ajustado = list(fit = fit,   v = v,   dis = diseno),
                total    = list(fit = fit_B, v = v_B, dis = diseno_B))

tab_cor  <- data.frame()
tab_sens <- data.frame()

# SECCIÓN 2: CORRELACIÓN ENTRE GENES Y SENSIBILIDAD DE CAMERAPR

for (mo in names(modelos)) {
  M <- modelos[[mo]]

  for (co in names(colecciones)) {

    genes  <- rownames(M$fit)[rownames(M$fit) %in% unlist(colecciones[[co]])]
    indice <- ids2indices(colecciones[[co]], genes)
    indice <- indice[lengths(indice) >= 15 & lengths(indice) <= 500]

    E   <- M$v$E[genes, ]
    rho <- sapply(indice, function(i) interGeneCorrelation(E[i, ], M$dis)$correlation)
    rho <- pmax(rho, 0)

    tab_cor <- rbind(tab_cor, data.frame(
      modelo = mo, coleccion = co, vias = length(rho),
      rho_mediana = round(median(rho), 3),
      rho_p25     = round(quantile(rho, 0.25), 3),
      rho_p75     = round(quantile(rho, 0.75), 3)))

    for (ct in names(contrastes)) {
      t_g <- M$fit$t[genes, contrastes[ct]]

      p_def <- cameraPR(t_g, indice, sort = FALSE)$PValue
      p_est <- sapply(names(indice), function(n)
                 cameraPR(t_g, indice[n], inter.gene.cor = rho[n])$PValue)
      sig_def <- p.adjust(p_def, "BH") < 0.05
      sig_est <- p.adjust(p_est, "BH") < 0.05

      tab_sens <- rbind(tab_sens, data.frame(
        modelo = mo, contraste = ct, coleccion = co, vias = length(indice),
        sig_cor_0.01     = sum(sig_def),
        sig_cor_estimada = sum(sig_est),
        se_mantienen     = sum(sig_def & sig_est)))
    }
  }
}

# SECCIÓN 3: EXPORTACIÓN

print(tab_cor)
print(tab_sens)
write_xlsx(list(Correlacion_genes = tab_cor, Sensibilidad_camera = tab_sens),
           paste0(dir_res, "Sensibilidad_vias_correlacion.xlsx"))
