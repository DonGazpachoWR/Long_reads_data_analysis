#!/usr/bin/env bash
# Predicted vs observed tissue mixtures (R² and RMSE) for the release-6.1 run
# (r61: PR 1 + 2 + 3 of DonGazpachoWR/SQANTI3), compared with the kb2 run
# (old min_prevalence implementation, rescue without evidence check) and with
# the run without prevalence filter (sin_filtro).
#
# Usage: ./correlacion_mezclas_r61.sh [OUT_DIR]

set -uo pipefail

BASE_DIR="/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis"
DATA_R61="$BASE_DIR/data/data_expression_matrix_r61"
DATA_KB2="$BASE_DIR/data/data_expression_matrix_prev2"      # QC prev2 + filter/rescue kb2
DATA_SIN_FILTRO="$BASE_DIR/data/data_expression_matrix_v2"  # ejecución sin prevalencia
OUT_DIR=${1:-"$BASE_DIR/output/benchmark_r61"}
CSV="$OUT_DIR/resultados_sqanti.csv"

cd "$(dirname "$0")" || exit 1
mkdir -p "$OUT_DIR"

# Parámetros: ext seq plataforma modo ruta ejecucion
PARAM_FILE=$(mktemp)
for s in isoseq masseq ont; do
  if [ "$s" == "ont" ]; then plats="bambu flair isoquant"; else plats="isoseq isocall bambu flair isoquant"; fi
  for p in $plats; do
    for m in raw qc fl rq; do
      echo "class $s $p $m $DATA_R61 r61"              >> "$PARAM_FILE"
      echo "TPM $s $p $m $DATA_R61 r61"                >> "$PARAM_FILE"
      echo "TPM $s $p $m $DATA_KB2 kb2"                >> "$PARAM_FILE"
      echo "TPM $s $p $m $DATA_SIN_FILTRO sin_filtro"  >> "$PARAM_FILE"
    done
  done
done

echo "modo_medida,ejecucion,tipo_de_seqs,tipo_de_plataformas,tipo_de_modos,comparacion,r_cuadrado,rmse,num_isoformas_modo,num_isoformas_filter_result,num_isoformas_fsm,isoform_sin_filtro_r" > "$CSV"

echo "[1/4] Calculando R² y RMSE ($(wc -l < "$PARAM_FILE") ejecuciones)..."
parallel --colsep ' ' -j 10 \
  Rscript correlacion_mezclas_individual.r {1} {2} {3} {4} {5} "$OUT_DIR" {6} \
  < "$PARAM_FILE" >> "$CSV"
rm "$PARAM_FILE"

# 21 filas por modo y bloque (isoseq 5x2 + masseq 5x1 + ont 3x2), 4 modos, 4 bloques
n_filas=$(( $(wc -l < "$CSV") - 1 ))
echo "     Filas generadas: $n_filas (esperadas 336)"
[ "$n_filas" -lt 336 ] && echo "AVISO: faltan filas. Revisa los mensajes 'ERROR en ...' de arriba."

echo "[2/4] Evolución de R² y RMSE por fase (r61)..."
Rscript stats_r2_fases.r "$OUT_DIR" r61 || exit 1

echo "[3/4] r61 frente a kb2 (TPM)..."
Rscript stats_delta_fases.r "$OUT_DIR" r61 kb2 || exit 1

echo "[4/4] r61 frente a sin filtro (TPM)..."
Rscript stats_delta_fases.r "$OUT_DIR" r61 sin_filtro || exit 1

echo "Terminado. Salidas en $OUT_DIR"
