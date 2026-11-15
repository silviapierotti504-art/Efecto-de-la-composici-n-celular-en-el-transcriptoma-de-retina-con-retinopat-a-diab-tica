# SECCIÓN 1: LIBRERÍAS Y CONTRASTES DE INTERÉS

library(limma)
library(msigdbr)
library(writexl)

contr      <- c(Control_vs_Diab = "grupoControl", NPDR_vs_Diab = "grupoNPDR")
salida_epi <- list()

# SECCIÓN 2: MAQUINARIA EPIGENÉTICA (LISTA A PRIORI)

maquinaria <- list(
  Metilacion_ADN     = c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3"),
  Desacetilasas_HDAC = c("HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6",
                         "HDAC7", "HDAC8", "HDAC9", "HDAC10", "HDAC11"),
  Sirtuinas          = c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7"),
  Acetiltransferasas = c("EP300", "CREBBP", "KAT2A", "KAT2B"),
  PRC2_H3K27me3      = c("EZH2", "EZH1", "SUZ12", "EED", "KDM6A", "KDM6B"),
  Memoria_metabolica = c("SETD7", "SUV39H1", "KDM1A", "EHMT2", "SETDB1")
)
genes_maq <- intersect(unlist(maquinaria), rownames(fit))
familia   <- setNames(rep(names(maquinaria), lengths(maquinaria)), unlist(maquinaria))

# SECCIÓN 3: MAQUINARIA EPIGENÉTICA GEN POR GEN

for (ct in names(contr)) {
  tt <- topTable(fit, coef = contr[ct], number = Inf, sort.by = "none")[genes_maq, ]
  tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
  tt <- data.frame(gen = rownames(tt), familia = familia[rownames(tt)], tt)
  tt$miRNAs_que_lo_regulan <- sapply(tt$gen, function(g)
    paste(names(sets)[sapply(sets, function(s) g %in% s)], collapse = ", "))
  salida_epi[[paste0("Maquinaria_", ct)]] <- tt[order(tt$P.Value), ]
}

# SECCIÓN 4: MAQUINARIA EPIGENÉTICA EN BLOQUE (FRY)

indice_maq <- c(list(Toda_la_maquinaria = match(genes_maq, rownames(v))),
                lapply(maquinaria, function(g) match(intersect(g, rownames(v)), rownames(v))))
indice_maq <- indice_maq[lengths(indice_maq) >= 3]

bloque <- data.frame()
for (ct in names(contr)) {
  r <- fry(v, index = indice_maq, design = diseno, contrast = which(colnames(diseno) == contr[ct]),
           block = donor, correlation = corr)
  bloque <- rbind(bloque, data.frame(comparacion = ct, conjunto = rownames(r), r))
}
salida_epi$Maquinaria_en_bloque <- bloque

# SECCIÓN 5: CONJUNTOS DE GENES REGULADOS EPIGENÉTICAMENTE

cgp      <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
cgp      <- cgp[grepl("PRC2|H3K27ME3|SUZ12|EED_TARGETS|EZH2|HDAC|DNMT|SIRT", cgp$gs_name), ]
sets_epi <- split(toupper(cgp$gene_symbol), cgp$gs_name)

# SECCIÓN 6: TEST DE GENES REGULADOS (CAMERAPR)

resumen_epi <- data.frame()
for (ct in names(contr)) {
  t_cod  <- fit$t[rownames(fit) %in% codif, contr[ct]]
  indice <- ids2indices(sets_epi, names(t_cod))
  indice <- indice[lengths(indice) >= 15 & lengths(indice) <= 1000]
  r <- cameraPR(t_cod, indice)
  r <- data.frame(conjunto = rownames(r), r)
  salida_epi[[paste0("Blancos_epi_", ct)]] <- r
  resumen_epi <- rbind(resumen_epi, data.frame(comparacion = ct, conjuntos = nrow(r),
                                               FDR_0.05 = sum(r$FDR < 0.05),
                                               suben    = sum(r$FDR < 0.05 & r$Direction == "Up"),
                                               bajan    = sum(r$FDR < 0.05 & r$Direction == "Down")))
}
salida_epi$Resumen_blancos_epi <- resumen_epi
print(resumen_epi)

# SECCIÓN 7: EXPORTACIÓN DE RESULTADOS

write_xlsx(salida_epi, paste0(dir_res, "Epigenetica_inferida.xlsx"))
