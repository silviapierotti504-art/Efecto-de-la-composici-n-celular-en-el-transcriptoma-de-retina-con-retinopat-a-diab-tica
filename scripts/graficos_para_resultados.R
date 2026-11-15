# SECCIÓN 1: LIBRERÍAS Y DIRECTORIOS

library(readxl)
library(ggplot2)
library(ggrepel)
library(scales)
library(grid)

dir_res <- "C:/Users/silvi/Desktop/Codigo_deconvolucion/GitHub/Resultados/"
dir_img <- "C:/Users/silvi/Desktop/Codigo_deconvolucion/GitHub/imagenes/"

coma <- label_number(decimal.mark = ",", big.mark = ".")
dec  <- function(x, d) format(round(x, d), nsmall = d, decimal.mark = ",")

tema <- theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 13),
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

col_ambos <- "#7B1FA2"
col_A     <- "#C2185B"
col_B     <- "#2C7BB6"
col_mayor <- "#C2185B"
col_menor <- "#2C7BB6"

# SECCIÓN 2: CARGA DE RESULTADOS

tejido <- as.data.frame(read_excel(paste0(dir_res, "DE_contrastes_interes.xlsx"), sheet = "Tejido_macula_vs_periferia"))
test_t <- as.data.frame(read_excel(paste0(dir_res, "DE_contrastes_interes.xlsx"), sheet = "Test_t_Control_vs_Diabetico"))
oxi    <- as.data.frame(read_excel(paste0(dir_res, "DE_estres_oxidativo.xlsx"), sheet = "Test_F_genes"))

leer_vias <- function(hoja, contraste) {
  d <- as.data.frame(read_excel(paste0(dir_res, "Cambios_globales_vias.xlsx"), sheet = hoja))
  d$contraste <- contraste
  d
}
gobp <- rbind(leer_vias("GOBP_Control_vs_Diabetico", "Control vs Diabético"),
              leer_vias("GOBP_RDNP_vs_Diabetico",    "RDNP vs Diabético"))
hall <- rbind(leer_vias("Hallmark_Control_vs_Diabetico", "Control vs Diabético"),
              leer_vias("Hallmark_RDNP_vs_Diabetico",    "RDNP vs Diabético"))

# SECCIÓN 3: FIGURA 06, MÁCULA VS PERIFERIA

n_genes <- data.frame(
  modelo = c("Modelo B\n(sin composición)", "Modelo A\n(ajustado por composición)"),
  genes  = c(sum(tejido$ModeloB_p_ajustado_FDR < 0.05), sum(tejido$ModeloA_p_ajustado_FDR < 0.05)))
n_genes$modelo <- factor(n_genes$modelo, levels = n_genes$modelo)

p_barras <- ggplot(n_genes, aes(x = modelo, y = genes, fill = modelo)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = coma(genes)), vjust = -0.4, size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c(col_B, col_A)) +
  scale_y_continuous(labels = coma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Genes con FDR < 5 %", title = "A") +
  tema

r2 <- data.frame(
  componente = rep(c("Solo región", "Solo composición\ncelular", "Compartida"), each = nrow(tejido)),
  valor      = 100 * c(tejido$R2_solo_tejido, tejido$R2_solo_composicion, tejido$R2_compartida))
r2$componente <- factor(r2$componente, levels = c("Solo región", "Solo composición\ncelular", "Compartida"))
medianas <- aggregate(valor ~ componente, data = r2, FUN = median)

p_r2 <- ggplot(r2, aes(x = componente, y = valor, fill = componente)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, show.legend = FALSE) +
  geom_label(data = medianas, aes(label = paste0("Mediana: ", dec(valor, 1), " %")),
             y = 75, size = 4.5, fill = "white", show.legend = FALSE) +
  scale_fill_manual(values = c("grey75", col_A, "#F4A6C0")) +
  coord_cartesian(ylim = c(0, 80)) +
  scale_y_continuous(labels = function(x) paste0(x, " %")) +
  labs(x = NULL, y = "Varianza explicada (%)", title = "B") +
  tema

png(paste0(dir_img, "06_macula_vs_periferia.png"), width = 13, height = 6.5, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_barras, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_r2,     vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# SECCIÓN 4: FIGURA 07, VÍAS GO:BP EN EL MODELO A FRENTE AL MODELO B

con_signo <- function(direccion, p) ifelse(direccion == "Mayor que en Diabetico", 1, -1) * -log10(p)
gobp$x <- con_signo(gobp$ModeloB_direccion, gobp$ModeloB_p_valor)
gobp$y <- con_signo(gobp$ModeloA_direccion, gobp$ModeloA_p_valor)

gobp$categoria <- "No significativa"
gobp$categoria[gobp$resultado == "Solo en Modelo B (explicada por la composicion celular)"] <- "Solo en el Modelo B"
gobp$categoria[gobp$resultado == "Solo en Modelo A"] <- "Solo en el Modelo A"
gobp$categoria[gobp$resultado == "Significativa en ambos modelos"] <- "En ambos modelos"
gobp$categoria <- factor(gobp$categoria, levels = c("En ambos modelos", "Solo en el Modelo A",
                                                    "Solo en el Modelo B", "No significativa"))
gobp <- gobp[order(gobp$categoria, decreasing = TRUE), ]

etiquetas_gobp <- data.frame(
  via      = c("GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS",
               "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_CYTOPLASMIC_TRANSLATION"),
  contraste = c("Control vs Diabético", "Control vs Diabético", "RDNP vs Diabético", "RDNP vs Diabético"),
  texto    = c("Fosforilación oxidativa", "Exocitosis de vesículas sinápticas",
               "Fosforilación oxidativa", "Traducción citoplasmática"))
etiquetas_gobp <- merge(etiquetas_gobp, gobp[, c("via", "contraste", "x", "y")])

lim <- max(abs(c(gobp$x, gobp$y)))

p_gobp <- ggplot(gobp, aes(x = x, y = y)) +
  annotate("rect", xmin = -lim, xmax = 0, ymin = 0, ymax = lim, fill = "grey95") +
  annotate("rect", xmin = 0, xmax = lim, ymin = -lim, ymax = 0, fill = "grey95") +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = categoria, size = categoria), alpha = 0.8) +
  geom_label_repel(data = etiquetas_gobp, aes(label = texto), size = 3.8, min.segment.length = 0,
                   box.padding = 0.8, fill = "white") +
  facet_wrap(~ contraste) +
  scale_colour_manual(values = c(col_ambos, col_A, col_B, "grey80"), drop = FALSE) +
  scale_size_manual(values = c(2.6, 2, 2, 1), drop = FALSE) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "Modelo B (sin composición): −log10(p) con signo",
       y = "Modelo A (ajustado por composición): −log10(p) con signo",
       colour = "Vía significativa (FDR < 5 %):", size = "Vía significativa (FDR < 5 %):") +
  tema

ggsave(paste0(dir_img, "07_vias_GOBP_modelo_A_vs_B.png"), p_gobp, width = 14, height = 8, dpi = 300, bg = "white")

# SECCIÓN 5: FIGURA 08, VÍAS HALLMARK SIGNIFICATIVAS

nombres_hall <- c(HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Fosforilación oxidativa",
                  HALLMARK_MTORC1_SIGNALING          = "Señalización de mTORC1",
                  HALLMARK_MYC_TARGETS_V1            = "Dianas de MYC (V1)",
                  HALLMARK_KRAS_SIGNALING_DN         = "Genes reprimidos por KRAS",
                  HALLMARK_G2M_CHECKPOINT            = "Punto de control G2/M",
                  HALLMARK_MITOTIC_SPINDLE           = "Huso mitótico",
                  HALLMARK_DNA_REPAIR                = "Reparación del ADN")
hall <- hall[hall$via %in% names(nombres_hall), ]

hall_A <- data.frame(via = hall$via, contraste = hall$contraste, modelo = "Modelo A",
                     direccion = hall$ModeloA_direccion, fdr = hall$ModeloA_p_ajustado_FDR)
hall_B <- data.frame(via = hall$via, contraste = hall$contraste, modelo = "Modelo B",
                     direccion = hall$ModeloB_direccion, fdr = hall$ModeloB_p_ajustado_FDR)
hall_l <- rbind(hall_A, hall_B)
hall_l$via <- factor(nombres_hall[hall_l$via], levels = rev(nombres_hall))
hall_l$estado <- ifelse(hall_l$fdr >= 0.05, "No significativa",
                 ifelse(hall_l$direccion == "Mayor que en Diabetico",
                        "Mayor que en diabéticos sin RD", "Menor que en diabéticos sin RD"))
hall_l$estado <- factor(hall_l$estado, levels = c("Mayor que en diabéticos sin RD",
                                                  "Menor que en diabéticos sin RD", "No significativa"))
hall_l$texto <- ifelse(hall_l$fdr < 0.001, "<0,001", dec(hall_l$fdr, 3))

p_hall <- ggplot(hall_l, aes(x = modelo, y = via)) +
  geom_point(aes(fill = estado), shape = 21, size = 15, colour = "white") +
  geom_text(aes(label = texto, colour = estado == "No significativa"), size = 3.3, fontface = "bold",
            show.legend = FALSE) +
  facet_wrap(~ contraste) +
  scale_fill_manual(values = c(col_mayor, col_menor, "grey88"), drop = FALSE) +
  scale_colour_manual(values = c("white", "grey35")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  tema + theme(panel.grid = element_blank())

ggsave(paste0(dir_img, "08_vias_hallmark.png"), p_hall, width = 11, height = 7, dpi = 300, bg = "white")

# SECCIÓN 6: FIGURA 09, PERFILES DE LOS GENES DE ESTRÉS OXIDATIVO SIGNIFICATIVOS

sig <- oxi[oxi$ModeloA_p_ajustado_FDR < 0.05, ]
efectos  <- cbind(Control = sig$ModeloA_log2FC_Control_vs_Diabetico, "Diabético" = 0,
                  RDNP = sig$ModeloA_log2FC_RDNP_vs_Diabetico)
relativa <- efectos - rowMeans(efectos)

patron <- ifelse(sig$ModeloA_log2FC_Control_vs_Diabetico > 0 & sig$ModeloA_log2FC_RDNP_vs_Diabetico < 0, "Control > Diabético > RDNP",
          ifelse(sig$ModeloA_log2FC_Control_vs_Diabetico < 0 & sig$ModeloA_log2FC_RDNP_vs_Diabetico > 0, "Control < Diabético < RDNP",
          ifelse(sig$ModeloA_log2FC_Control_vs_Diabetico > 0, "Menor en Diabético", "Mayor en Diabético")))
niveles <- c("Control > Diabético > RDNP", "Control < Diabético < RDNP", "Menor en Diabético", "Mayor en Diabético")
n_pat   <- table(factor(patron, levels = niveles))
etiq    <- setNames(paste0(niveles, "\n(n = ", n_pat, ")"), niveles)

df_pf <- data.frame(gen    = rep(sig$gen, times = 3),
                    grupo  = rep(colnames(relativa), each = nrow(sig)),
                    valor  = as.vector(relativa),
                    patron = rep(patron, times = 3))
df_pf$grupo  <- factor(df_pf$grupo, levels = c("Control", "Diabético", "RDNP"))
df_pf$patron <- factor(etiq[df_pf$patron], levels = etiq[n_pat > 0])
mediana <- aggregate(valor ~ grupo + patron, data = df_pf, FUN = median)

p_oxi <- ggplot(df_pf, aes(x = grupo, y = valor)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(aes(group = gen), colour = "#F4A6C0", linewidth = 0.6) +
  geom_line(data = mediana, aes(group = 1), colour = col_A, linewidth = 1.6) +
  geom_point(data = mediana, colour = col_A, size = 3) +
  geom_text_repel(data = df_pf[df_pf$grupo == "Control", ], aes(label = gen), size = 3, nudge_x = -0.35,
                  direction = "y", hjust = 1, segment.colour = "grey70", max.overlaps = Inf, fontface = "italic") +
  facet_wrap(~ patron, nrow = 1) +
  scale_x_discrete(expand = expansion(add = c(0.9, 0.3))) +
  scale_y_continuous(labels = coma) +
  labs(x = NULL, y = "Expresión relativa (log2, centrada por gen)") +
  tema

ggsave(paste0(dir_img, "09_perfiles_estres_oxidativo.png"), p_oxi,
       width = 4.5 * sum(n_pat > 0), height = 7, dpi = 300, bg = "white")

# SECCIÓN 7: FIGURA 10, CAMBIO DEL CONJUNTO DE ESTRÉS OXIDATIVO FRENTE AL RESTO DE LOS GENES

z_signo <- function(log2fc, p) sign(log2fc) * qnorm(p / 2, lower.tail = FALSE)
es_oxi  <- test_t$gen %in% oxi$gen

dist <- rbind(
  data.frame(modelo = "Modelo A (ajustado por composición)",
             z = z_signo(test_t$ModeloA_log2FC, test_t$ModeloA_p_valor), oxi = es_oxi),
  data.frame(modelo = "Modelo B (sin composición)",
             z = z_signo(test_t$ModeloB_log2FC, test_t$ModeloB_p_valor), oxi = es_oxi))
dist$conjunto <- factor(ifelse(dist$oxi, paste0("Estrés oxidativo (", sum(es_oxi), " genes)"),
                               "Resto de los genes"))
dist$conjunto <- relevel(dist$conjunto, ref = "Resto de los genes")
med_dist <- aggregate(z ~ modelo + conjunto, data = dist, FUN = median)


p_dist <- ggplot(dist, aes(x = z, fill = conjunto, colour = conjunto)) +
  geom_density(alpha = 0.35, linewidth = 0.8, adjust = 1.2) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(data = med_dist, aes(xintercept = z, colour = conjunto), linetype = "dashed", linewidth = 0.9,
             show.legend = FALSE) +
  facet_wrap(~ modelo, ncol = 1) +
  scale_fill_manual(values = c("grey60", col_A)) +
  scale_colour_manual(values = c("grey40", col_A)) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "Cambio de cada gen en Control vs Diabético (estadístico z con signo)", y = "Densidad",
       fill = NULL, colour = NULL) +
  tema

ggsave(paste0(dir_img, "10_conjunto_estres_oxidativo.png"), p_dist, width = 11, height = 8.5, dpi = 300, bg = "white")

# SECCIÓN 8: FIGURA 11, PERFILES DE LOS miRNAs SIGNIFICATIVOS EN EL TEST F

mir <- as.data.frame(read_excel(paste0(dir_res, "miRNAs_resultados.xlsx"), sheet = "DE_F_global"))
mir <- mir[mir$adj.P.Val < 0.05, ]
mir$nombre <- sub("^hsa-", "", mir$miRNA)

efectos_m  <- cbind(Control = mir$grupoControl, "Diabético" = 0, RDNP = mir$grupoNPDR)
relativa_m <- efectos_m - rowMeans(efectos_m)

patron_m <- ifelse(mir$grupoControl > 0 & mir$grupoNPDR < 0, "Control > Diabético > RDNP",
            ifelse(mir$grupoControl < 0 & mir$grupoNPDR > 0, "Control < Diabético < RDNP",
            ifelse(mir$grupoControl > 0, "Menor en Diabético", "Mayor en Diabético")))
n_pat_m <- table(factor(patron_m, levels = niveles))
etiq_m  <- setNames(paste0(niveles, "\n(n = ", n_pat_m, ")"), niveles)

df_m <- data.frame(mirna  = rep(mir$nombre, times = 3),
                   grupo  = rep(colnames(relativa_m), each = nrow(mir)),
                   valor  = as.vector(relativa_m),
                   patron = rep(patron_m, times = 3))
df_m$grupo  <- factor(df_m$grupo, levels = c("Control", "Diabético", "RDNP"))
df_m$patron <- factor(etiq_m[df_m$patron], levels = etiq_m[n_pat_m > 0])
mediana_m <- aggregate(valor ~ grupo + patron, data = df_m, FUN = median)

p_mir <- ggplot(df_m, aes(x = grupo, y = valor)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(aes(group = mirna), colour = "#F4A6C0", linewidth = 0.6) +
  geom_line(data = mediana_m, aes(group = 1), colour = col_A, linewidth = 1.6) +
  geom_point(data = mediana_m, colour = col_A, size = 3) +
  geom_text_repel(data = df_m[df_m$grupo == "Control", ], aes(label = mirna), size = 3, nudge_x = -0.35,
                  direction = "y", hjust = 1, segment.colour = "grey70", max.overlaps = Inf) +
  facet_wrap(~ patron, nrow = 1) +
  scale_x_discrete(expand = expansion(add = c(1.1, 0.3))) +
  scale_y_continuous(labels = coma) +
  labs(x = NULL, y = "Expresión relativa (log2, centrada por miRNA)") +
  tema

ggsave(paste0(dir_img, "11_perfiles_miRNAs.png"), p_mir, width = 18, height = 7.5, dpi = 300, bg = "white")

# SECCIÓN 9: FIGURA 12, MAQUINARIA EPIGENÉTICA (LOG2FC E INTERVALO DE CONFIANZA DEL 95 %)

leer_maq <- function(hoja, contraste) {
  d <- as.data.frame(read_excel(paste0(dir_res, "Epigenetica_inferida.xlsx"), sheet = hoja))
  d$contraste <- contraste
  d
}
maq <- rbind(leer_maq("Maquinaria_Control_vs_Diab", "Control vs Diabético"),
             leer_maq("Maquinaria_NPDR_vs_Diab",    "RDNP vs Diabético"))

nombres_fam <- c(Metilacion_ADN     = "Metilación\ndel ADN",
                 Desacetilasas_HDAC = "Desacetilasas\n(HDAC)",
                 Sirtuinas          = "Sirtuinas",
                 Acetiltransferasas = "Acetil-\ntransferasas",
                 PRC2_H3K27me3      = "PRC2 /\nH3K27me3",
                 Memoria_metabolica = "Memoria\nmetabólica")
maq$familia <- factor(nombres_fam[maq$familia], levels = nombres_fam)
maq$ee      <- maq$logFC / maq$t
maq$ic_inf  <- maq$logFC - 1.96 * maq$ee
maq$ic_sup  <- maq$logFC + 1.96 * maq$ee
maq$gen     <- factor(maq$gen, levels = rev(unique(maq$gen[order(maq$familia, maq$gen)])))

p_maq <- ggplot(maq, aes(x = logFC, y = gen)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_errorbar(aes(xmin = ic_inf, xmax = ic_sup), width = 0, colour = col_B, linewidth = 0.8) +
  geom_point(colour = col_B, size = 2.2) +
  facet_grid(familia ~ contraste, scales = "free_y", space = "free_y") +
  scale_x_continuous(labels = coma) +
  labs(x = "log2FC (intervalo de confianza del 95 %)", y = NULL) +
  tema +
  theme(axis.text.y = element_text(face = "italic", size = 10),
        strip.text.y = element_text(angle = 0, size = 11))

ggsave(paste0(dir_img, "12_maquinaria_epigenetica.png"), p_maq, width = 10, height = 11, dpi = 300, bg = "white")

# SECCIÓN 10: FIGURA 13, FRACCIÓN ESTIMADA FRENTE A OBSERVADA EN EL CONJUNTO DE YANG ET AL. (2026)

dir_yang <- "C:/Users/silvi/Desktop/TFM_Bioinformatica/Sc vs Sc/"

estimada <- read.csv(paste0(dir_yang, "resultados_bayesprism/02_deconvolucion/proporciones_estimadas_A_sana_GSE137846.csv"))
celulas  <- read.csv(paste0(dir_yang, "resultados_DE_por_tipo/pseudobulks_por_tipo.csv"))

tipos <- c("Rod", "Cone", "BC", "AC", "HC", "RGC", "MGC", "Microglia")
nombres_tipo <- c(Rod = "Bastones", Cone = "Conos", BC = "Bipolares", AC = "Amacrinas",
                  HC = "Horizontales", RGC = "CGR", MGC = "Glía de Müller", Microglia = "Microglía")

observada <- tapply(celulas$n_celulas, list(celulas$sample, celulas$celltype), sum)
observada[is.na(observada)] <- 0
observada <- observada[, tipos]
observada <- observada / rowSums(observada)

estimada  <- estimada[match(rownames(observada), estimada$sample), ]
grupo_y   <- celulas$group[match(rownames(observada), celulas$sample)]
grupo_y   <- factor(c(nonDM = "Sin diabetes", DM = "Diabetes sin RD", DR = "RD")[grupo_y],
                    levels = c("Sin diabetes", "Diabetes sin RD", "RD"))

comp <- data.frame()
for (ct in tipos) {
  comp <- rbind(comp, data.frame(tipo      = nombres_tipo[[ct]],
                                 observada = 100 * observada[, ct],
                                 estimada  = 100 * estimada[[ct]],
                                 grupo     = grupo_y))
}
comp$tipo <- factor(comp$tipo, levels = nombres_tipo)

r_tipo <- data.frame(tipo = factor(nombres_tipo, levels = nombres_tipo),
                     r    = sapply(nombres_tipo, function(t) cor(comp$observada[comp$tipo == t],
                                                                 comp$estimada[comp$tipo == t])))
r_tipo$texto <- paste0("r = ", dec(r_tipo$r, 2))

p_yang <- ggplot(comp, aes(x = observada, y = estimada)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = grupo), size = 2.6, alpha = 0.85) +
  geom_text(data = r_tipo, aes(x = -Inf, y = Inf, label = texto), hjust = -0.15, vjust = 1.5,
            size = 4.2, fontface = "bold") +
  facet_wrap(~ tipo, nrow = 2, scales = "free") +
  scale_colour_manual(values = c("#4C9A6A", "#D9A23B", "#C2185B")) +
  scale_x_continuous(labels = coma) +
  scale_y_continuous(labels = coma) +
  labs(x = "Fracción observada (% de células)", y = "Fracción estimada con GSE137846 (%)", colour = NULL) +
  tema

error <- data.frame(grupo = grupo_y,
                    rmse  = 100 * sqrt(rowMeans((as.matrix(estimada[, tipos]) - observada)^2)))

p_rmse <- ggplot(error, aes(x = grupo, y = rmse, fill = grupo)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.5, show.legend = FALSE) +
  geom_point(aes(colour = grupo), size = 2.6, position = position_jitter(width = 0.08, seed = 1), show.legend = FALSE) +
  scale_fill_manual(values = c("#4C9A6A", "#D9A23B", "#C2185B")) +
  scale_colour_manual(values = c("#4C9A6A", "#D9A23B", "#C2185B")) +
  scale_y_continuous(labels = coma, limits = c(0, NA)) +
  labs(x = NULL, y = "Error por muestra (RMSE, puntos %)", title = "A") +
  tema

p_yang <- p_yang + labs(title = "B")

png(paste0(dir_img, "13_referencia_estimada_vs_observada.png"), width = 13, height = 12, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3, heights = unit(c(1, 1.7), "null"))))
print(p_rmse, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_yang, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:3))
dev.off()
