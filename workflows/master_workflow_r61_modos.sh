#!/bin/bash
# Release-6.1 test with the evidence check that uses the redistribution of
# requantification (DonGazpachoWR/SQANTI3 6d0ac36), in two modes run on Garnatxa
# by scripts/garnatxa/prueba_release61/sq_array_r61.sh:
#   prevalencia: resultados/         (array 3245677, rescue rerun)
#   sinprev:     resultados_sinprev/ (array 3245678, filter and rescue without prevalence)
#
#   1. checks that both arrays have finished,
#   2. prevalencia: downloads, homogenises, mixture correlations and verification
#      (master_workflow_r61.sh),
#   3. sinprev: downloads filter and rescue (QC and counts are shared) and homogenises,
#   4. randomisation test and gene bootstrap for both modes.
#
# Usage: ./master_workflow_r61_modos.sh [B] [CPUS]

B=${1:-1000}
CPUS=${2:-10}
BASE_DIR="/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis"
DATA_R61="$BASE_DIR/data/data_expression_matrix_r61"
DATA_SINPREV="$BASE_DIR/data/data_expression_matrix_r61_sinprev"
OUT_DIR="$BASE_DIR/output/benchmark_r61"
SCRIPTS_DIR="$BASE_DIR/scripts/scripts_r/release61"
GARNATXA_HOST="adrianbe@garnatxa.uv.es"
REMOTO=/home/adrianbe/practicas/prueba_release61
JOBS="3245677,3245678"

echo "[1/4] Estado de los arrays $JOBS en Garnatxa..."
estado=$(ssh -o BatchMode=yes -o ConnectTimeout=15 "$GARNATXA_HOST" "sacct -j $JOBS -X -n -o JobID%16,State%12") || {
    echo "ERROR: no hay conexión con Garnatxa."; exit 1; }
echo "$estado" | awk '{print $2}' | sort | uniq -c
if echo "$estado" | grep -qE "PENDING|RUNNING|REQUEUED"; then
    echo "Aún hay tareas en marcha. Vuelve a lanzarlo cuando terminen."; exit 1
fi
if echo "$estado" | grep -vq "COMPLETED"; then
    echo "AVISO: hay tareas que no terminaron bien (revisa $REMOTO/log/r61_<job>_<tarea>.err):"
    echo "$estado" | grep -v COMPLETED
    exit 1
fi

echo "[2/4] Modo prevalencia: descarga, correlación de mezclas y verificación..."
# Los gráficos existentes no se redibujan: se apartan los de la ejecución anterior
mkdir -p "$OUT_DIR/previo_evidencia"
for m in class TPM; do
    [ -d "$OUT_DIR/graficos/_${m}_r61" ] && mv "$OUT_DIR/graficos/_${m}_r61" "$OUT_DIR/previo_evidencia/_${m}_r61"
done
for f in resultados_sqanti.csv verificacion test_reasignacion_azar; do
    [ -e "$OUT_DIR/$f" ] && mv "$OUT_DIR/$f" "$OUT_DIR/previo_evidencia/"
done
"$BASE_DIR/workflows/master_workflow_r61.sh" "$REMOTO/resultados" || exit 1

echo "[3/4] Modo sin prevalencia: descarga..."
mkdir -p "$DATA_SINPREV"
cd "$DATA_SINPREV" || exit 1
rm -f rules_default_* rq_* counts_design_*
ssh "$GARNATXA_HOST" "cd $REMOTO/resultados_sinprev && find . -maxdepth 3 -type f \\( \
      -path '*/filter/rules_default_*_RulesFilter_classification.txt' \
   -o -path '*/filter/rules_default_*_filtering_reasons.txt' \
   -o -path '*/rescue/rq_*_rescued_classification.txt' \
   -o -path '*/rescue/rq_*_reassigned_counts*.tsv' \
   -o -path '*/rescue/rq_*_rescue_table.tsv' \
   -o -path '*/rescue/rq_*_rescue_inclusion_list.tsv' \
   -o -path '*/rescue/rq_*_rescue_summary.tsv' \
   -o -path '*/rescue/rq_*_SQANTI3_rescue_report.pdf' \
   -o -name 'counts_design_*.json' \\) -print0 | tar -czf - --null -T - --transform='s|.*/||'" | tar -xzf - || exit 1
# QC y counts son los mismos que en el modo prevalencia
cp "$DATA_R61"/counts_*.tsv "$DATA_R61"/default_*_classification.txt .
n=$(ls rq_*_rescued_classification.txt 2>/dev/null | wc -l)
echo "     rescues descargados: $n / 13"
[ "$n" -ge 13 ] || { echo "ERROR: descarga incompleta."; exit 1; }
cd "$BASE_DIR/data" || exit 1
./homogeneizar_nombres.sh "$DATA_SINPREV" > /dev/null || exit 1

echo "[4/4] Test de aleatorización y bootstrap fl vs rq (B = $B)..."
# Copia congelada del script: se puede seguir editando el original mientras corre
TMP=$(mktemp -d)
cp "$SCRIPTS_DIR"/test_reasignacion_azar_r61.r "$SCRIPTS_DIR"/requant_r61.r "$SCRIPTS_DIR"/tema_release61.r "$TMP"/
for modo in r61 sinprev; do
    datos=$DATA_R61; [ "$modo" = sinprev ] && datos=$DATA_SINPREV
    salida="$OUT_DIR/test_reasignacion_azar/$modo"
    mkdir -p "$salida"
    Rscript "$TMP/test_reasignacion_azar_r61.r" "$datos" "$salida" "$modo" "$B" "$CPUS" > "$salida/log.txt" 2>&1 \
        || { echo "ERROR en el test de $modo (ver $salida/log.txt)"; exit 1; }
    echo "     $modo terminado: $salida"
done
rm -rf "$TMP"

echo "TERMINADO. Salidas en $OUT_DIR (verificacion/, test_reasignacion_azar/r61 y sinprev)"
