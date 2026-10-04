# Long_reads_data_analysis

Benchmark de transcriptomas de long reads (isoseq, masseq, ont) con mezclas de tejidos
(K, B y mezclas B20K80 / B80K20) y su procesado con SQANTI3.

## Estructura

| Carpeta | Contenido |
|---|---|
| `data/` | Matrices y salidas de SQANTI descargadas de Garnatxa (no versionado). `homogeneizar_nombres.sh` unifica los nombres de muestra. |
| `output/` | Gráficas y tablas generadas (no versionado). |
| `workflows/` | Scripts maestros que descargan de Garnatxa y lanzan los análisis en local. `master_workflow_r61.sh` es el de las pull requests de release 6.1. |
| `docs/` | Documentación: `contexto/` (contexto de cada PR), `diagramas/` (flujo del rescue con control de evidencia), `referencias/` (diseño de SIRVs), `resultados_sueltos/` (CSV y salidas de análisis antiguos). |
| `scripts/garnatxa/` | Scripts de Slurm para Garnatxa: `isocall/`, `SQ/`, `tama/` (antiguos) y `prueba_release61/` (prueba de las PRs de release 6.1). |
| `scripts/lanzadores/` | Lanzadores locales de los análisis de filtro de réplicas (`lanzar_*.sh`). |
| `scripts/scripts_r/release61/` | Scripts nuevos para las pull requests de mi repositorio (DonGazpachoWR/SQANTI3, rama `release-6.1`). |
| `scripts/scripts_r/antiguo/` | Scripts anteriores, clasificados por categoría. |

### `scripts/scripts_r/release61/`

- `correlacion_mezclas_r61.sh`: R² y RMSE de las mezclas observadas frente a las esperadas, por fase de SQANTI (raw, qc, fl, rq), para la ejecución r61 frente a kb2 y sin filtro.
- `correlacion_mezclas_individual.r`: cálculo de una combinación, fase y ejecución (lo llama el anterior).
- `stats_r2_fases.r`, `stats_delta_fases.r`: evolución por fase y comparación entre ejecuciones.
- `verificar_prs_r61.r`: comprobaciones de cada PR sobre los datos reales (min_intron_length, requantificación, prevalencias de QC, filtro frente a kb2 y control de evidencia recalculado).
- `tema_release61.r`: colores y tema comunes.

### `scripts/scripts_r/antiguo/`

| Categoría | Scripts |
|---|---|
| `01_matrices_conteos` | Fusión de matrices de conteo y TSV para TAMA. |
| `02_correlacion_mezclas_filtros_R` | Dotplots de mezclas y estadísticas de R²/RMSE con los filtros de expresión en R (`master_workflow.sh`, `master_workflow_prev2.sh`). |
| `03_filtro_replicados` | Diseño del filtro de réplicas: alpha, epsilon, umbral n, tres brazos, qué filtra, diagnóstico antisense. |
| `04_controles_sirvs` | Controles con SIRVs. |
| `05_prevalencia_sqanti_kb2` | Análisis de la ejecución kb2 (min_prevalence dentro de SQANTI, implementación anterior). |
| `06_rescue_prototipos` | Prototipos de rescue coherente con la prevalencia y línea vertical del rescue. |
| `07_diagnosticos` | Solapamiento de tejidos y violines de conteos crudos frente a QC. |
