#!/bin/bash
#SBATCH --job-name r61
#SBATCH --qos medium
#SBATCH --mem 60G
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --array 1-13%5
#SBATCH -o log/r61_%A_%a.out
#SBATCH -e log/r61_%A_%a.err
# ------------------------------------------------------------------------------
# Test of the three pull requests of DonGazpachoWR/SQANTI3 (release-6.1):
#   PR 1: requant by default, min_intron_length, rescue report
#   PR 2: QC prevalence_<group> columns (--counts_design), generic rules filter,
#         JSON cleaned in rescue before filtering the reference
#   PR 3: evidence check of the reintroduced reference transcripts
#
# Same inputs and filter as the kb2 run: expressed in at least 2 replicates of
# K OR of B, now written as OR rules on prevalence_K / prevalence_B. Expression
# uses the default --min_expression of QC (0: count > 0, so the fractional counts
# of bambu are kept; kb2 used count >= 1).
# Each step is skipped if its output already exists, so a relaunch resumes
# (delete the step directory to repeat it).
# At the end, regression against kb2 (filter) and evidence_kb2 (rescue).
# ------------------------------------------------------------------------------
set -euo pipefail

module load anaconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sqanti_conda

BASE=/home/adrianbe/practicas/prueba_release61
SQ=$BASE/sqanti
JSON=$BASE/filter_r61.json
params=$BASE/combinaciones.txt
fasta=/storage/gge/home_members/adrianbe/practicas/intento1/genoma_raton.fa
gtf=/storage/gge/home_members/adrianbe/practicas/intento1/copia_limpio_anotacion_raton.gtf
path_referencia=/storage/gge/home_members/adrianbe/practicas/intento1
ref_class=$BASE/resultados/referencia/reference_default_classification.txt

read -r seq plataforma < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$params")
if [ -z "${seq:-}" ] || [ -z "${plataforma:-}" ]; then
  echo "ERROR: no hay línea $SLURM_ARRAY_TASK_ID en $params" >&2
  exit 1
fi
combi=${seq}_${plataforma}
path=$path_referencia/$seq/$plataforma
counts=$path/counts_${combi}.tsv
out=$BASE/resultados/$combi
qc_dir=$out/qc
filter_dir=$out/filter
rq_dir=$out/rescue
old_filter=$path/sqanti/filter/rules_kb2
old_rq=$path/sqanti/rq/evidence_kb2
mkdir -p "$qc_dir" "$filter_dir" "$rq_dir"

echo "=== Tarea $SLURM_ARRAY_TASK_ID: $combi ==="
echo "SQANTI3: $(cat "$SQ/COMMIT")"
[ -f "$ref_class" ] || { echo "ERROR: falta el QC de referencia ($ref_class)" >&2; exit 1; }

# --- Counts design: K and B replicates; mixtures stay out ---
design=$out/counts_design_${combi}.json
awk -F'\t' 'NR == 1 {
    for (i = 2; i <= NF; i++) {
        if ($i ~ /(^|_)K3[0-9]$/)      k = k (k ? ", " : "") "\"" $i "\""
        else if ($i ~ /(^|_)B3[0-9]$/) b = b (b ? ", " : "") "\"" $i "\""
    }
    printf "{\"K\": [%s], \"B\": [%s]}\n", k, b
    exit
}' "$counts" > "$design"
if [ "$seq" = "masseq" ]; then n_rep=3; else n_rep=5; fi
n_k=$(grep -o '"[^"]*K3[0-9]"' "$design" | wc -l)
n_b=$(grep -o '"[^"]*B3[0-9]"' "$design" | wc -l)
echo ">>> Diseño: $(cat "$design")"
if [ "$n_k" -ne "$n_rep" ] || [ "$n_b" -ne "$n_rep" ]; then
  echo "ERROR: se esperaban $n_rep réplicas de K y de B, hay K=$n_k B=$n_b" >&2
  exit 1
fi

# --- QC (PR 1 min_intron_length, PR 2 prevalence_<group>) ---
qc_class=$qc_dir/default_${combi}_classification.txt
if [ -s "$qc_class" ]; then
  echo ">>> QC ya hecho, se reutiliza"
else
  echo ">>> QC"
  python "$SQ/sqanti3_qc.py" \
    --isoforms "$path/transcriptoma_${combi}.gtf" \
    --fl_count "$counts" \
    --counts_design "$design" \
    --refGTF "$gtf" \
    --refFasta "$fasta" \
    --report skip \
    --dir "$qc_dir" \
    --output "default_${combi}" \
    -t "$SLURM_CPUS_PER_TASK"
fi
header=$(head -1 "$qc_class" | tr -d '\r')   # DictWriter writes \r\n line endings
for col in min_intron_length prevalence prevalence_K prevalence_B; do
  echo "$header" | tr '\t' '\n' | grep -qx "$col" || { echo "ERROR: falta la columna $col en QC" >&2; exit 1; }
done
echo "    columnas min_intron_length, prevalence, prevalence_K, prevalence_B presentes"
echo "    $(grep MinExpression "$qc_dir/default_${combi}.qc_params.txt")"

# --- Filter (PR 2: prevalence through the generic rules) ---
filter_class=$filter_dir/rules_default_${combi}_RulesFilter_classification.txt
if [ -s "$filter_class" ]; then
  echo ">>> Filter ya hecho, se reutiliza"
else
  echo ">>> Filter"
  python "$SQ/sqanti3_filter.py" rules \
    --sqanti_class "$qc_class" \
    --filter_gtf "$qc_dir/default_${combi}_corrected.gtf" \
    --json_filter "$JSON" \
    --dir "$filter_dir" \
    --output "rules_default_${combi}"
fi

# --- Rescue (PR 2 JSON cleaning, PR 3 evidence check, PR 1 report and requant by default) ---
if [ -s "$rq_dir/rq_${combi}_reassigned_counts.tsv" ]; then
  echo ">>> Rescue ya hecho, se reutiliza"
else
  echo ">>> Rescue"
  python "$SQ/sqanti3_rescue.py" \
    -s rules \
    --filter_class "$filter_class" \
    --refGTF "$gtf" \
    --refFasta "$fasta" \
    --refClassif "$ref_class" \
    --mode full \
    --json_filter "$JSON" \
    --counts_design "$design" \
    --corrected_isoforms_fasta "$qc_dir/default_${combi}_corrected.fasta" \
    --filtered_isoforms_gtf "$filter_dir/rules_default_${combi}.filtered.gtf" \
    --counts "$counts" \
    --dir "$rq_dir" \
    --output "rq_${combi}"
fi

# --- Checks ---
echo ">>> Comprobaciones"
echo "    JSON de referencia: $(tr -d ' \n' < "$rq_dir/reference_rules_filter/reference_rules.json")"
[ -s "$rq_dir/rq_${combi}_SQANTI3_rescue_report.pdf" ] && echo "    informe de rescue: OK" || echo "    informe de rescue: FALTA"
[ -s "$rq_dir/rq_${combi}_reassigned_counts.tsv" ] && echo "    requant por defecto: OK" || echo "    requant por defecto: FALTA"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="evidence_check") c=i; next} c{n[$c]++}
            END{for(k in n) print "    evidence_check", k, n[k]}' "$rq_dir/rq_${combi}_rescue_table.tsv"

echo ">>> Regresión"
cmp_sorted() {   # $1 etiqueta, $2 antiguo, $3 nuevo
  if [ ! -e "$2" ]; then echo "SIN_REFERENCIA  $1"
  elif cmp -s <(sort "$2") <(sort "$3"); then echo "IDENTICO  $1"
  else echo "DISTINTO  $1 ($(diff <(sort "$2") <(sort "$3") | grep -c '^[<>]') líneas)"; fi
}
cut_cols() {     # isoform y filter_result de una classification
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="isoform")a=i; if($i=="filter_result")b=i}; next}{print $a"\t"$b}' "$1"
}
old_fclass=$old_filter/rules_default_${combi}_RulesFilter_classification.txt
if [ -f "$old_fclass" ]; then
  cmp_sorted "filter_result (vs kb2)" <(cut_cols "$old_fclass") <(cut_cols "$filter_class")
else
  echo "SIN_REFERENCIA  filter_result (vs kb2)"
fi
for f in rescue_table.tsv rescue_inclusion_list.tsv reassigned_counts.tsv; do
  cmp_sorted "$f (vs evidence_kb2)" "$old_rq/rq_${combi}_$f" "$rq_dir/rq_${combi}_$f"
done

# Previous r61 run with count >= 1 (moved to resultados_ge1 before relaunching bambu)
prev=$BASE/resultados_ge1/$combi
if [ -d "$prev" ]; then
  echo ">>> Comparación con la ejecución r61 anterior (count >= 1)"
  cmp_sorted "filter_result (vs r61 >= 1)" \
    <(cut_cols "$prev/filter/rules_default_${combi}_RulesFilter_classification.txt") <(cut_cols "$filter_class")
  awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++) if($i=="filter_result") c=i; next}
               {n[(FILENAME==ARGV[1] ? "antes" : "ahora") " " $c]++}
               END{for(k in n) print "    " k, n[k]}' \
    "$prev/filter/rules_default_${combi}_RulesFilter_classification.txt" "$filter_class" | sort
fi

echo "=== $combi terminado ==="
