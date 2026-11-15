# SECCIÓN 1: LIBRERÍAS

library(limma)
library(ggplot2)
library(scales)
library(readxl)
library(writexl)
set.seed(2024)

# SECCIÓN 2: PARÁMETROS Y TAMAÑOS DE EFECTO DE REFERENCIA

cambios    <- seq(0.05, 0.60, by = 0.05)
frac_genes <- c(1, 0.5)
n_rep      <- 100
q          <- 0.05
contraste  <- "grupoControl"

referencias <- data.frame(cambio   = c(0.30, 0.60),
                          etiqueta = c("GPX1 -30 %\n(ratón)", "SOD2 -60 %\n(rata)"))

nombre_modelo <- c(A = "Modelo A (ajustado por composición)",
                   B = "Modelo B (sin composición)")

idx_fit <- match(oxi, rownames(fit))
idx_v   <- match(oxi, rownames(v))
n_oxi   <- length(oxi)

pct  <- label_percent(decimal.mark = ",")
tema <- theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.box = "vertical", strip.text = element_text(face = "bold"))

# SECCIÓN 3: POTENCIA GEN POR GEN (CÁLCULO ANALÍTICO)

umbrales <- c("FDR 5 % (umbral conservador: 0,05 / n.º de genes)" = q / n_oxi,
              "p < 0,05 sin corregir (límite optimista)"          = q)

potencia_gen <- function(f, modelo) {
  se  <- f$stdev.unscaled[idx_fit, contraste] * sqrt(f$s2.post[idx_fit])
  df  <- f$df.total[idx_fit]
  res <- data.frame()
  for (u in names(umbrales)) {
    for (x in cambios) {
      d  <- abs(log2(1 - x))
      tc <- qt(1 - umbrales[[u]] / 2, df)
      pw <- pt(tc, df, ncp = d / se, lower.tail = FALSE) + pt(-tc, df, ncp = d / se)
      res <- rbind(res, data.frame(modelo = nombre_modelo[[modelo]], umbral = u, cambio = x,
                                   potencia_mediana = median(pw),
                                   potencia_p25 = unname(quantile(pw, 0.25)),
                                   potencia_p75 = unname(quantile(pw, 0.75))))
    }
  }
  res
}
gen <- rbind(potencia_gen(fit, "A"), potencia_gen(fit_B, "B"))

# SECCIÓN 4: POTENCIA DEL CONJUNTO COMPLETO (SIMULACIÓN CON FRY)

correr_simulacion <- FALSE

if (correr_simulacion) {
  donantes_cd <- unique(as.character(donor[grupo %in% c("Control", "Diabetic")]))
  n_control   <- length(unique(as.character(donor[grupo == "Control"])))
  es_cd       <- grupo %in% c("Control", "Diabetic")

  una_repeticion <- function(vv, modelo, corr_m, d, f) {
    ctrl <- sample(donantes_cd, n_control)
    g <- as.character(grupo)
    g[es_cd] <- ifelse(as.character(donor[es_cd]) %in% ctrl, "Control", "Diabetic")
    g <- relevel(factor(g), ref = "Diabetic")
    if (modelo == "A") dis <- model.matrix(~ g + tejido + edad + sexo + rin + pcs)
    if (modelo == "B") dis <- model.matrix(~ g + tejido + edad + sexo + rin)
    cambian <- idx_v
    if (f < 1) cambian <- sample(idx_v, round(f * length(idx_v)))
    E <- vv
    E$E[cambian, g == "Control"] <- E$E[cambian, g == "Control"] + d
    r <- fry(E, index = list(estres_oxidativo = idx_v), design = dis,
             contrast = which(colnames(dis) == "gControl"), block = donor, correlation = corr_m)
    r$PValue < q
  }

  conjunto <- data.frame()
  for (mo in c("A", "B")) {
    if (mo == "A") { vv <- v;   corr_m <- corr   }
    if (mo == "B") { vv <- v_B; corr_m <- corr_B }
    for (f in frac_genes) {
      for (x in cambios) {
        det <- replicate(n_rep, una_repeticion(vv, mo, corr_m, abs(log2(1 - x)), f))
        ic  <- binom.test(sum(det), n_rep)$conf.int
        conjunto <- rbind(conjunto, data.frame(
          modelo = nombre_modelo[[mo]], genes_que_cambian = paste0(f * 100, " % de los genes"),
          cambio = x, potencia = mean(det), IC95_inferior = ic[1], IC95_superior = ic[2]))
      }
    }
  }
} else {
  conjunto <- as.data.frame(read_excel(paste0(dir_res, "Potencia_estres_oxidativo.xlsx"),
                                       sheet = "Conjunto_completo"))
}

# SECCIÓN 5: FIGURAS

capa_referencias <- list(
  geom_hline(yintercept = 0.8, linetype = "dotted", colour = "grey40"),
  geom_vline(data = referencias, aes(xintercept = cambio), inherit.aes = FALSE,
             linetype = "dashed", colour = "grey50"),
  geom_text(data = referencias, aes(x = cambio, y = 1.02, label = etiqueta), inherit.aes = FALSE,
            vjust = 0, size = 3.5, colour = "grey30", lineheight = 0.9),
  scale_x_continuous(labels = pct, breaks = seq(0.1, 0.6, by = 0.1)),
  scale_y_continuous(labels = pct, breaks = seq(0, 1, by = 0.2)),
  coord_cartesian(ylim = c(0, 1.12)),
  scale_colour_manual(values = c("#C2185B", "#2C7BB6")),
  scale_fill_manual(values = c("#C2185B", "#2C7BB6"))
)

p_gen <- ggplot(gen, aes(x = cambio, y = potencia_mediana, colour = modelo, linetype = umbral)) +
  geom_ribbon(data = gen[gen$umbral == names(umbrales)[1], ],
              aes(ymin = potencia_p25, ymax = potencia_p75, fill = modelo),
              alpha = 0.15, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.1) +
  capa_referencias +
  labs(x = "Disminución de la expresión en diabéticos (por gen)",
       y = "Potencia (mediana de los genes)", colour = NULL, linetype = NULL,
       title = "Potencia gen por gen: genes de estrés oxidativo",
       subtitle = paste0("Control vs Diabético; ", n_oxi, " genes (GO:0006979). Banda: rango intercuartílico entre genes")) +
  tema

ggsave(paste0(dir_img, "04_curva_potencia_estres_oxidativo_genes.png"), p_gen,
       width = 10, height = 7.5, dpi = 300, bg = "white")

p_conj <- ggplot(conjunto, aes(x = cambio, y = potencia, colour = modelo, linetype = genes_que_cambian)) +
  geom_ribbon(aes(ymin = IC95_inferior, ymax = IC95_superior, fill = modelo,
                  group = interaction(modelo, genes_que_cambian)),
              alpha = 0.12, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.1) +
  capa_referencias +
  labs(x = "Disminución de la expresión en diabéticos (por gen)",
       y = "Potencia (proporción de simulaciones con p < 0,05)", colour = NULL, linetype = NULL,
       title = "Potencia del conjunto completo: estrés oxidativo",
       subtitle = paste0("Control vs Diabético; test fry, ", n_rep, " simulaciones por punto. Banda: IC 95 %")) +
  tema

ggsave(paste0(dir_img, "05_curva_potencia_estres_oxidativo_conjunto.png"), p_conj,
       width = 10, height = 7.5, dpi = 300, bg = "white")

# SECCIÓN 6: EXPORTACIÓN

parametros <- data.frame(
  parametro = c("Contraste", "Genes analizados", "Cambios evaluados", "Umbral gen por gen",
                "Umbral del conjunto", "Simulaciones por punto", "Fracción de genes que cambian",
                "Referencias marcadas"),
  valor = c("Control vs Diabetico", paste0(n_oxi, " genes de GO:0006979"), "5 % a 60 % de disminucion",
            "FDR 5 % conservador (0,05 / n.º de genes) y p < 0,05 sin corregir",
            "p < 0,05 (test fry)", n_rep, "100 % y 50 %",
            "GPX1 -30 % (Malaviya y Kowluru, 2024); SOD2 -60 % (Zhong y Kowluru, 2011)")
)

write_xlsx(list(Parametros = parametros, Gen_por_gen = gen, Conjunto_completo = conjunto),
           paste0(dir_res, "Potencia_estres_oxidativo.xlsx"))
