#!/bin/bash
#SBATCH --job-name r61_ref
#SBATCH --qos medium
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH -o log/ref_%j.out
#SBATCH -e log/ref_%j.err
# ------------------------------------------------------------------------------
# QC of the reference annotation against itself with the release-6.1 code
# (adds min_intron_length). Needed by rescue; run once before the array.
# ------------------------------------------------------------------------------
set -euo pipefail

module load anaconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sqanti_conda

BASE=/home/adrianbe/practicas/prueba_release61
SQ=$BASE/sqanti
fasta=/storage/gge/home_members/adrianbe/practicas/intento1/genoma_raton.fa
gtf=/storage/gge/home_members/adrianbe/practicas/intento1/copia_limpio_anotacion_raton.gtf
out=$BASE/resultados/referencia

echo "SQANTI3: $(cat "$SQ/COMMIT")"
mkdir -p "$out"
python "$SQ/sqanti3_qc.py" \
  --isoforms "$gtf" \
  --refGTF "$gtf" \
  --refFasta "$fasta" \
  --report skip \
  --dir "$out" \
  --output reference_default \
  -t "$SLURM_CPUS_PER_TASK"

echo "QC de referencia terminado"
