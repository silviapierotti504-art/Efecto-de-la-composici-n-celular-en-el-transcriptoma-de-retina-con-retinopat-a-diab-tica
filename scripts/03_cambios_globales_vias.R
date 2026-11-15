# SECCIÓN 1: LIBRERÍAS

library(limma)
library(msigdbr)
library(writexl)

# SECCIÓN 2: COLECCIONES DE VÍAS (MSIGDB)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
gobp     <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")

colecciones <- list(Hallmark = split(toupper(hallmark$gene_symbol), hallmark$gs_name),
                    GOBP     = split(toupper(gobp$gene_symbol),     gobp$gs_name))

# SECCIÓN 3: CONTRASTES DE INTERÉS Y MODELOS

contrastes <- c(Control_vs_Diab = "grupoControl", NPDR_vs_Diab = "grupoNPDR")
nombre_ct  <- c(Control_vs_Diab = "Control_vs_Diabetico", NPDR_vs_Diab = "RDNP_vs_Diabetico")

clasificar_vias <- function(fdr_A, fdr_B) {
  ifelse(fdr_A < 0.05 & fdr_B < 0.05, "Significativa en ambos modelos",
  ifelse(fdr_B < 0.05,                "Solo en Modelo B (explicada por la composicion celular)",
  ifelse(fdr_A < 0.05,                "Solo en Modelo A",
                                      "No significativa")))
}

test_via <- function(f, co, cf, modelo) {
  t_genes <- f$t[, cf]
  t_co    <- t_genes[names(t_genes) %in% unlist(colecciones[[co]])]
  indice  <- ids2indices(colecciones[[co]], names(t_co))
  indice  <- indice[lengths(indice) >= 15 & lengths(indice) <= 500]
  r <- cameraPR(t_co, indice, sort = FALSE)
  d <- data.frame(rownames(r), r$NGenes, direccion(r$Direction), r$PValue, r$FDR)
  names(d) <- c("via", "genes_en_via", paste0(modelo, c("_direccion", "_p_valor", "_p_ajustado_FDR")))
  d
}

# SECCIÓN 4: TEST POR VÍA (CAMERAPR) EN AMBOS MODELOS

salida  <- list()
resumen <- data.frame()

for (co in names(colecciones)) {
  for (ct in names(contrastes)) {

    a   <- test_via(fit,   co, contrastes[ct], "ModeloA")
    b   <- test_via(fit_B, co, contrastes[ct], "ModeloB")
    tab <- cbind(a, b[, 3:5])
    tab$resultado <- clasificar_vias(tab$ModeloA_p_ajustado_FDR, tab$ModeloB_p_ajustado_FDR)
    tab <- tab[order(tab$ModeloA_p_valor), ]
    salida[[paste0(co, "_", nombre_ct[ct])]] <- tab

    sigA <- tab$ModeloA_p_ajustado_FDR < 0.05
    sigB <- tab$ModeloB_p_ajustado_FDR < 0.05
    resumen <- rbind(resumen, data.frame(
      coleccion          = co,
      comparacion        = gsub("_", " ", nombre_ct[ct]),
      vias_analizadas    = nrow(tab),
      ModeloA_vias_FDR_5 = sum(sigA),
      ModeloA_mayor      = sum(sigA & tab$ModeloA_direccion == "Mayor que en Diabetico"),
      ModeloA_menor      = sum(sigA & tab$ModeloA_direccion == "Menor que en Diabetico"),
      ModeloB_vias_FDR_5 = sum(sigB),
      ModeloB_mayor      = sum(sigB & tab$ModeloB_direccion == "Mayor que en Diabetico"),
      ModeloB_menor      = sum(sigB & tab$ModeloB_direccion == "Menor que en Diabetico")))
  }
}

# SECCIÓN 5: RESUMEN Y EXPORTACIÓN

print(resumen)
write_xlsx(c(list(Resumen = resumen), salida), paste0(dir_res, "Cambios_globales_vias.xlsx"))
