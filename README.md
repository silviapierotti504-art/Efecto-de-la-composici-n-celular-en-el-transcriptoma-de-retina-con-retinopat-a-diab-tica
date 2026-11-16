# Deconvolución bayesiana del transcriptoma de retina humana para el análisis del estrés oxidativo en la retinopatía diabética

Material complementario del Trabajo Fin de Máster. Máster Universitario en Bioinformática, UNIR.

## Resumen del trabajo

Los estudios de transcriptómica bulk de retina con retinopatía diabética suelen atribuir a la enfermedad cambios de expresión que podrían deberse a cambios en la composición celular del tejido. En este trabajo se estimó la composición celular del bulk GSE160306 por deconvolución con BayesPrism, usando como referencia el single-cell de retina sana GSE137846. Luego se compararon dos modelos de expresión diferencial, el Modelo A (ajustado por composición celular) y el Modelo B (sin ajustar), con foco en los genes de estrés oxidativo. La referencia se evaluó con el single-cell de retina sana, diabética y con retinopatía de Yang et al. (2026).

## Contenido del repositorio

| Carpeta | Contenido |
|---|---|
| `scripts/` | Códigos de R del trabajo |
| `docs/` | Un Word por script con la explicación de cada sección, sus datos de entrada y sus resultados |
| `Resultados/01_referencia/` | Tablas y figuras de la referencia single-cell GSE137846 |
| `Resultados/02_bulk/` | Tablas y figuras del análisis del bulk |
| `Resultados/03_Yang/` | Tablas y figuras de la evaluación de la referencia y del análisis por tipo celular con Yang et al. (2026) |
| `Resultados/04_supuestos/` | Tablas y figuras de los supuestos y el ajuste de los modelos |

## Scripts y orden de ejecución

Las rutas de los directorios se definen al inicio de cada script y deben adaptarse a la carpeta donde se guarden los datos.

### 1. Referencia single-cell

| Script | Qué hace | Explicación |
|---|---|---|
| `00._descripcion_atlas_referencia.R` | Control de calidad y caracterización de la referencia GSE137846. Guarda la referencia anotada que usa el script 01 | `00._Descripcion_script_atlas_referencia.docx` |

### 2. Análisis del bulk

El script 01 se ejecuta primero. Los scripts 02 a 08 se ejecutan a continuación, en orden y en la misma sesión de R, porque usan los objetos que crean los anteriores.

| Script | Qué hace | Explicación |
|---|---|---|
| `01_bulk_deconvolucion_diferencias.R` | Deconvolución del bulk con BayesPrism, ajuste del Modelo A y del Modelo B, partición de varianza, contrastes y análisis de los genes de estrés oxidativo | `01_Tabla_codigo_bulk_deconvolucion_diferencias.docx` |
| `02_composicion_por_enfermedad.R` | Cambios en la composición celular según el grupo, con modelos mixtos por tipo celular | `02_Tabla_codigo_composicion_por_enfermedad.docx` |
| `03_cambios_globales_vias.R` | Análisis de vías (MSigDB) con cameraPR en ambos modelos | `03_Tabla_codigo_cambios_globales_vias.docx` |
| `04_sensibilidad_vias_correlacion.R` | Sensibilidad de cameraPR a la correlación entre genes | `04_Tabla_codigo_sensibilidad_vias_correlacion.docx` |
| `05_efecto_minimo_detectable.R` | Efecto mínimo detectable por gen y pacientes necesarios | `05_Tabla_codigo_efecto_minimo_detectable.docx` |
| `06_potencia_estres_oxidativo.R` | Potencia gen por gen y del conjunto de estrés oxidativo | `06_Tabla_codigo_potencia_estres_oxidativo.docx` |
| `07_miRNAs.R` | Expresión diferencial de miRNAs (GSE160308) y análisis de sus genes blanco | `07_Tabla_codigo_miRNAs.docx` |
| `08_epigenetica_inferida.R` | Maquinaria epigenética y genes regulados epigenéticamente | `08_Tabla_codigo_epigenetica_inferida.docx` |
| `09_supuestos_y_ajuste_modelos.R` | Supuestos y bondad de ajuste de los modelos. Se ejecuta en una sesión nueva, porque vuelve a ajustar los modelos | `09_Tabla_codigo_supuestos_y_ajuste_modelos.docx` |

### 3. Evaluación con Yang et al. (2026)

Se ejecutan en este orden. Los tres se explican en `Yang_Tabla_codigo_scripts.docx`.

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

Los datos no se incluyen en el repositorio por su tamaño. Se descargan de sus fuentes originales.

| Conjunto | Uso |
|---|---|
| GSE160306 | Bulk de retina (RNA total) |
| GSE160308 | Bulk de retina (small RNA, miRNAs) |
| GSE137846 | Referencia single-cell de retina sana |
| Yang et al. (2026) | Single-cell de retina sana, diabética y con retinopatía |
