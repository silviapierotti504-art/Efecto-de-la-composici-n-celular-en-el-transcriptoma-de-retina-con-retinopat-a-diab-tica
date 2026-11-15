# Deconvolución bayesiana del transcriptoma de retina humana para el análisis del estrés oxidativo en la retinopatía diabética

Material complementario del Trabajo Fin de Máster de Silvia Pierotti (Máster Universitario en Bioinformática, UNIR).

## Resumen del trabajo

Los estudios de transcriptómica bulk de retina con retinopatía diabética suelen atribuir a la enfermedad cambios de expresión que podrían deberse a cambios en la composición celular del tejido. En este trabajo se estimó la composición celular del bulk GSE160306 por deconvolución con BayesPrism, usando como referencia el single-cell de retina sana GSE137846. Luego se compararon dos modelos de expresión diferencial, el Modelo A (ajustado por composición celular) y el Modelo B (sin ajustar), con foco en los genes de estrés oxidativo. La referencia se evaluó con el single-cell de retina sana, diabética y con retinopatía de Yang et al. (2026).

## Contenido del repositorio

| Carpeta | Contenido |
|---|---|
| `scripts/` | Códigos de R del trabajo |
| `docs/` | Explicación de los scripts por sección |
| `Resultados/` | Tablas de resultados que generan los scripts |
| `imagenes/` | Figuras que generan los scripts |
| `datos/` | No se incluye por su tamaño. Ver la sección Datos |

## Scripts y orden de ejecución

Todos los scripts se ejecutan desde la carpeta raíz del repositorio, abriendo el proyecto de RStudio, porque las rutas son relativas a esa carpeta.

### 1. Referencia single-cell

| Script | Qué hace |
|---|---|
| `00._descripcion_atlas_referencia.R` | Control de calidad y caracterización de la referencia GSE137846. Guarda la referencia anotada en `datos/objetos/`, que usa el script 01 |

### 2. Análisis del bulk

El script 01 se ejecuta primero. Los scripts 02 a 08 se ejecutan a continuación en la misma sesión de R, porque usan los objetos que crea el script 01.

| Script | Qué hace |
|---|---|
| `01_bulk_deconvolucion_diferencias.R` | Deconvolución del bulk con BayesPrism, ajuste del Modelo A y del Modelo B, partición de varianza, contrastes y análisis de los genes de estrés oxidativo |
| `02_composicion_por_enfermedad.R` | Cambios en la composición celular según el grupo, con modelos mixtos por tipo celular |
| `03_cambios_globales_vias.R` | Análisis de vías (MSigDB) con cameraPR en ambos modelos |
| `04_sensibilidad_vias_correlacion.R` | Sensibilidad de cameraPR a la correlación entre genes |
| `05_efecto_minimo_detectable.R` | Efecto mínimo detectable por gen y pacientes necesarios |
| `06_potencia_estres_oxidativo.R` | Potencia gen por gen y del conjunto de estrés oxidativo |
| `07_miRNAs.R` | Expresión diferencial de miRNAs (GSE160308) y análisis de sus genes blanco |
| `08_epigenetica_inferida.R` | Maquinaria epigenética y genes regulados epigenéticamente |
| `09_supuestos_y_ajuste_modelos.R` | Supuestos y bondad de ajuste de los modelos. Se puede ejecutar en una sesión nueva, porque lee las fracciones que guarda el script 01 |

### 3. Evaluación de la referencia con Yang et al. (2026)

Se ejecutan en este orden. La explicación detallada está en `docs/Scripts_Yang_resumen_por_seccion.docx`.

| Script | Qué hace |
|---|---|
| `Yang_deconvolucion_retina_scRNAseq.R` | Evalúa si la referencia sana deconvoluciona bien retina diabética y con retinopatía, usando pseudobulks de composición conocida |
| `Yang_DE_por_tipo_celular_estres_oxidativo.R` | Expresión diferencial de los genes de estrés oxidativo dentro de cada tipo celular |
| `Yang_DE_por_tipo_celular_estres_oxidativo_sensibilidad.R` | Análisis de sensibilidad del script anterior |

### 4. Figuras de resultados

| Script | Qué hace |
|---|---|
| `graficos_para_resultados.R` | Genera las figuras de la sección de resultados del trabajo. Se ejecuta al final, porque lee los resultados de los pasos 2 y 3 |

## Datos

Los datos no se incluyen en el repositorio por su tamaño. Se descargan de sus fuentes originales y se colocan en la carpeta `datos/`.

| Conjunto | Uso | Carpeta |
|---|---|---|
| GSE160306 | Bulk de retina (RNA total) | `datos/` |
| GSE160308 | Bulk de retina (small RNA, miRNAs) | `datos/` |
| GSE137846 | Referencia single-cell de retina sana | `datos/GSE137846/` |
| Yang et al. (2026) | Single-cell de retina sana, diabética y con retinopatía | `datos/Data/` |
