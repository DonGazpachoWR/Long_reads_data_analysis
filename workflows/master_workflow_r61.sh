#!/bin/bash
# Release-6.1 test (DonGazpachoWR/SQANTI3 PR 1 + 2 + 3), run on Garnatxa by
# scripts/garnatxa/prueba_release61/sq_array_r61.sh:
#   1. downloads its outputs to data/data_expression_matrix_r61 and homogenises
#      the sample names,
#   2. predicted vs observed tissue mixtures (R², RMSE) against kb2 and sin_filtro,
#   3. checks of each pull request (verificar_prs_r61.r).
#
# Usage: ./master_workflow_r61.sh [garnatxa_results_dir]

BASE_DIR="/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis"
DATA_DIR="$BASE_DIR/data/data_expression_matrix_r61"
KB2_DIR="$BASE_DIR/data/data_expression_matrix_prev2"
OUT_DIR="$BASE_DIR/output/benchmark_r61"
SCRIPTS_DIR="$BASE_DIR/scripts/scripts_r/release61"

GARNATXA_HOST="adrianbe@garnatxa.uv.es"
RES=${1:-/home/adrianbe/practicas/prueba_release61/resultados}
INPUTS=/storage/gge/home_members/adrianbe/practicas/intento1
N_COMBI=13

echo "============================================================"
echo " RELEASE 6.1: descarga, correlación de mezclas y verificación"
echo " Remoto: $RES"
echo " Local:  $DATA_DIR"
echo "============================================================"

if ! ssh -o ConnectTimeout=10 "$GARNATXA_HOST" "test -d $RES"; then
    echo "ERROR: no hay conexión con Garnatxa (¿VPN de la UV?) o no existe $RES."
    exit 1
fi
mkdir -p "$DATA_DIR" "$OUT_DIR"
cd "$DATA_DIR" || exit 1

# ---------------------------------------------------------------- PASO 1
echo "[PASO 1] Descargando salidas de SQANTI (qc, filter, rescue)..."
rm -f counts_*.tsv default_*_classification.txt rules_default_* rq_*
ssh "$GARNATXA_HOST" "cd $RES && find . -maxdepth 3 -type f \\( \
      -path '*/qc/default_*_classification.txt' \
   -o -path '*/filter/rules_default_*_RulesFilter_classification.txt' \
   -o -path '*/filter/rules_default_*_filtering_reasons.txt' \
   -o -path '*/rescue/rq_*_rescued_classification.txt' \
   -o -path '*/rescue/rq_*_reassigned_counts*.tsv' \
   -o -path '*/rescue/rq_*_rescue_table.tsv' \
   -o -path '*/rescue/rq_*_rescue_inclusion_list.tsv' \
   -o -path '*/rescue/rq_*_rescue_summary.tsv' \
   -o -path '*/rescue/rq_*_rescue_artifact_outcomes.tsv' \
   -o -path '*/rescue/rq_*_SQANTI3_rescue_report.pdf' \
   -o -name 'counts_design_*.json' \\) -print0 | tar -czf - --null -T - --transform='s|.*/||'" | tar -xzf -
echo " -> Descargando counts crudos..."
ssh "$GARNATXA_HOST" "cd $INPUTS && find -L . -maxdepth 3 -name 'counts_*.tsv' -print0 | tar -chzf - --null -T - --transform='s|.*/||'" | tar -xzf -
echo " -> Descargando logs de la regresión..."
mkdir -p "$OUT_DIR/logs_garnatxa"
ssh "$GARNATXA_HOST" "cd $RES/../log && tar -czf - r61_*.out ref_*.out" | tar -xzf - -C "$OUT_DIR/logs_garnatxa"

ok=1
for patron in 'counts_*.tsv' 'default_*_classification.txt' 'rules_default_*_RulesFilter_classification.txt' \
              'rq_*_rescued_classification.txt' 'rq_*_reassigned_counts.tsv' 'rq_*_rescue_table.tsv' \
              'rq_*_rescue_summary.tsv' 'rq_*_SQANTI3_rescue_report.pdf'; do
    n=$(ls $patron 2>/dev/null | wc -l)
    printf "     %-50s %2d / %d\n" "$patron" "$n" "$N_COMBI"
    [ "$n" -ge "$N_COMBI" ] || ok=0
done
[ "$ok" -eq 1 ] || { echo "ERROR: descarga incompleta."; exit 1; }

# ---------------------------------------------------------------- PASO 2
echo "[PASO 2] Homogeneizando nombres de muestras..."
cd "$BASE_DIR/data" || exit 1
./homogeneizar_nombres.sh "$DATA_DIR" > /dev/null || exit 1

# ---------------------------------------------------------------- PASO 3
echo "[PASO 3] Correlación predicted vs observed de las mezclas..."
"$SCRIPTS_DIR/correlacion_mezclas_r61.sh" "$OUT_DIR" || exit 1

# ---------------------------------------------------------------- PASO 4
echo "[PASO 4] Verificación de las pull requests..."
Rscript "$SCRIPTS_DIR/verificar_prs_r61.r" "$DATA_DIR" "$KB2_DIR" "$OUT_DIR/verificacion" || exit 1

echo "============================================================"
echo " TERMINADO. Salidas en $OUT_DIR"
echo "============================================================"
