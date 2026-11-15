# SECCIÓN 1: LIBRERÍAS

library(limma)
library(edgeR)
library(multiMiR)
library(writexl)

# SECCIÓN 2: CARGA Y ALINEACIÓN DE LOS CONTEOS DE miRNA

mir <- as.data.frame(data.table::fread(paste0(dir, "GSE160308_human_retina_DR_smallRNA_counts.txt")))
rownames(mir) <- mir[[1]]
mir <- as.matrix(mir[, vapply(mir, is.numeric, logical(1)), drop = FALSE])

comunes <- intersect(m, colnames(mir))
idx     <- match(comunes, m)
mir     <- mir[, comunes]

# SECCIÓN 3: MODELO DE EXPRESIÓN DIFERENCIAL DE miRNAs

dm <- data.frame(grupo = droplevels(grupo[idx]), tejido = tejido[idx],
                 edad = edad[idx], sexo = sexo[idx], rin = rin[idx])
dm$pcs <- pcs[idx, ]
diseno_mir <- model.matrix(~ grupo + tejido + edad + sexo + rin + pcs, data = dm)

ym <- DGEList(mir)
ym <- ym[filterByExpr(ym, diseno_mir, min.count = 20), , keep.lib.sizes = FALSE]
ym <- calcNormFactors(ym)

vm      <- voom(ym, diseno_mir)
corr_m  <- duplicateCorrelation(vm, diseno_mir, block = donor[idx])$consensus
vm      <- voom(ym, diseno_mir, block = donor[idx], correlation = corr_m)
fit_mir <- eBayes(lmFit(vm, diseno_mir, block = donor[idx], correlation = corr_m))

# SECCIÓN 4: CONTRASTES DE miRNAs

salida_mir <- list()
for (g in c("grupoControl", "grupoNPDR")) {
  r <- topTable(fit_mir, coef = g, number = Inf, adjust.method = "BH")
  salida_mir[[paste0("DE_", sub("grupo", "", g), "_vs_Diab")]] <- data.frame(miRNA = rownames(r), r)
}
r <- topTable(fit_mir, coef = c("grupoControl", "grupoNPDR"), number = Inf, adjust.method = "BH")
salida_mir$DE_F_global <- data.frame(miRNA = rownames(r), r)

# SECCIÓN 5: miRNAs CANDIDATOS Y GENES BLANCO (miRTarBase)

candidatos <- unique(c(salida_mir$DE_Control_vs_Diab$miRNA[salida_mir$DE_Control_vs_Diab$P.Value < 0.05],
                       salida_mir$DE_NPDR_vs_Diab$miRNA[salida_mir$DE_NPDR_vs_Diab$P.Value < 0.05]))

blancos <- get_multimir(org = "hsa", mirna = candidatos, table = "mirtarbase", summary = FALSE)@data
blancos <- blancos[blancos$support_type == "Functional MTI", ]
blancos <- unique(blancos[, c("mature_mirna_id", "target_symbol")])
sets    <- split(toupper(blancos$target_symbol), blancos$mature_mirna_id)

# SECCIÓN 6: GENES CODIFICANTES (FONDO DE COMPARACIÓN)

anot  <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86, "GENEID"),
                               keytype = "GENEID", columns = c("SYMBOL", "GENEBIOTYPE"))
codif <- unique(toupper(anot$SYMBOL[anot$GENEBIOTYPE == "protein_coding"]))

# SECCIÓN 7: TEST DE GENES BLANCO (CAMERAPR)

for (g in c("grupoControl", "grupoNPDR")) {
  t_genes <- fit$t[rownames(fit) %in% codif, g]
  indice  <- ids2indices(sets, names(t_genes))
  indice  <- indice[lengths(indice) >= 10 & lengths(indice) <= 500]
  r <- cameraPR(t_genes, indice)
  r <- data.frame(miRNA = rownames(r), r)
  tt_mir <- topTable(fit_mir, coef = g, number = Inf, sort.by = "none")
  r$logFC_miRNA <- tt_mir[r$miRNA, "logFC"]
  r$coherente   <- sign(r$logFC_miRNA) != ifelse(r$Direction == "Up", 1, -1)
  salida_mir[[paste0("Blancos_", sub("grupo", "", g), "_vs_Diab")]] <- r
}

# SECCIÓN 8: BLANCOS DE miRNA EN GENES DE ESTRÉS OXIDATIVO

oxi_cod  <- intersect(oxi, codif)
sets_oxi <- lapply(sets, function(genes) intersect(genes, oxi_cod))

desc <- data.frame(miRNA           = names(sets),
                   blancos_totales = lengths(sets),
                   blancos_estres  = lengths(sets_oxi),
                   genes_estres    = sapply(sets_oxi, paste, collapse = ", "))
desc$pct_estres <- round(100 * desc$blancos_estres / desc$blancos_totales, 1)
salida_mir$Blancos_estres_descriptivo <- desc[order(-desc$blancos_estres), ]

# SECCIÓN 9: TEST POR miRNA DENTRO DEL ESTRÉS OXIDATIVO

for (g in c("grupoControl", "grupoNPDR")) {
  t_oxi  <- fit$t[oxi_cod, g]
  indice <- ids2indices(sets_oxi, names(t_oxi))
  indice <- indice[lengths(indice) >= 5]
  r <- cameraPR(t_oxi, indice)
  r <- data.frame(miRNA = rownames(r), r)
  tt_mir <- topTable(fit_mir, coef = g, number = Inf, sort.by = "none")
  r$logFC_miRNA <- tt_mir[r$miRNA, "logFC"]
  r$P_miRNA     <- tt_mir[r$miRNA, "P.Value"]
  r$coherente   <- sign(r$logFC_miRNA) != ifelse(r$Direction == "Up", 1, -1)
  r$genes       <- sapply(r$miRNA, function(x) paste(sets_oxi[[x]], collapse = ", "))
  salida_mir[[paste0("Estres_", sub("grupo", "", g), "_vs_Diab")]] <- r
}

# SECCIÓN 10: EXPORTACIÓN DE RESULTADOS

write_xlsx(salida_mir, paste0(dir_res, "miRNAs_resultados.xlsx"))
