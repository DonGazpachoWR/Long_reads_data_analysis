#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 50G
#SBATCH -c 30
#SBATCH -t 5:00:00
#SBATCH -o log/log.out
#SBATCH -e log/error.out

seq=ont
plataforma=bambu


path=/storage/gge/home_members/adrianbe/practicas/intento1/$seq/$plataforma
sqanti_path=/storage/gge/home_members/adrianbe/practicas/SQANTI3
path_referencia=/storage/gge/home_members/adrianbe/practicas/intento1

fasta=/storage/gge/home_members/adrianbe/practicas/intento1/genoma_raton.fa
gtf=/storage/gge/home_members/adrianbe/practicas/intento1/copia_limpio_anotacion_raton.gtf


mkdir -p $path/sqanti

path_out=$path/sqanti

mkdir -p $path_out/qc
mkdir -p $path_out/filter
mkdir -p $path_out/filter/rules
mkdir -p $path_out/rq
mkdir -p $path_out/rq/complete_rescue

mkdir -p $path_referencia/sqanti_qc





module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sqanti_conda

set -eou pipefail

# QC
python $sqanti_path/sqanti3_qc.py \
  --isoforms $path/transcriptoma_${seq}_${plataforma}.gtf \
  --fl_count $path/counts_${seq}_${plataforma}.tsv \
  --refGTF $gtf \
  --refFasta $fasta \
  --dir $path_out/qc/muestra \
  --output default_${seq}_${plataforma} \
  -t 30
echo QC con éxito
# Filter
python  $sqanti_path/sqanti3_filter.py rules \
  --sqanti_class $path_out/qc/muestra/default_${seq}_${plataforma}_classification.txt \
  --filter_gtf $path_out/qc/muestra/default_${seq}_${plataforma}_corrected.gtf \
  --json_filter $sqanti_path/src/utilities/filter/filter_default.json \
  --dir $path_out/filter/rules \
  --output rules_default_${seq}_${plataforma}
echo Filter con éxito
# Rescue, primero QC del GTF de referencia contra sí mismo. RefGTF e isoformas mismo archivo
if [ ! -f   $path_referencia/sqanti_qc/reference_default_classification.txt ]; then

    python $sqanti_path/sqanti3_qc.py \
    --isoforms $gtf \
    --refGTF $gtf \
    --refFasta $fasta \
    --dir $path_referencia/sqanti_qc/ \
    --output reference_default \
    -t 30

  echo QC para la referencia con éxito
else
  echo QC para la referencia encontrado
fi

# Rescue parte 2
python $sqanti_path/sqanti3_rescue.py \
  -s rules \
  --filter_class $path_out/filter/rules/rules_default_${seq}_${plataforma}_RulesFilter_classification.txt \
  --refGTF $gtf \
  --refFasta $fasta \
  --refClassif $path_referencia/sqanti_qc/reference_default_classification.txt \
  --mode full \
  --json_filter $sqanti_path/src/utilities/filter/filter_default.json \
  --corrected_isoforms_fasta $path_out/qc/muestra/default_${seq}_${plataforma}_corrected.fasta \
  --filtered_isoforms_gtf $path_out/filter/rules/rules_default_${seq}_${plataforma}_RulesFilter.filtered.gtf \
  --dir $path_out/rq/complete_rescue \
  --output rq_${seq}_${plataforma}

echo Rescue con éxito


/home/adrianbe/practicas/intento1/ont/bambu/transcriptoma_ont_bambu.gtf
