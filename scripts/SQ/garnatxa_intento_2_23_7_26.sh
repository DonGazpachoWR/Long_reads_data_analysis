#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 300G
#SBATCH -c 30
#SBATCH -t 24:00:00
#SBATCH -o log/log.out
#SBATCH -e log/error.out

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sqanti_conda

fasta=/storage/gge/home_members/adrianbe/practicas/intento1/genoma_raton.fa
gtf=/storage/gge/home_members/adrianbe/practicas/intento1/copia_limpio_anotacion_raton.gtf

for seq in ont isoseq;
do
  echo --------------------------------------------
  echo                    $seq
  echo --------------------------------------------
  echo
  echo Procesando $seq ...
  for plataforma in bambu flair isoquant;
  do
    echo    -----------------------------------
    echo                $plataforma
    echo    -----------------------------------
    echo
    echo      Procesando $plataforma


    path=/storage/gge/home_members/adrianbe/practicas/intento1/$seq/$plataforma
    sqanti_path=/storage/gge/home_members/adrianbe/practicas/SQANTI3
    path_referencia=/storage/gge/home_members/adrianbe/practicas/intento1

    mkdir -p $path/sqanti

    path_out=$path/sqanti

    mkdir -p $path_out/qc
    mkdir -p $path_out/filter
    mkdir -p $path_out/filter/rules
    mkdir -p $path_out/rq
    mkdir -p $path_out/rq/complete_rescue

    mkdir -p $path_referencia/sqanti_qc

    set -eou pipefail

    # QC
    echo      Ejecutando QC...
    python $sqanti_path/sqanti3_qc.py \
      --isoforms $path/transcriptoma_${seq}_${plataforma}.gtf \
      --fl_count $path/counts_${seq}_${plataforma}.tsv \
      --refGTF $gtf \
      --refFasta $fasta \
      --dir $path_out/qc/muestra \
      --output default_${seq}_${plataforma} \
      -t $SLURM_CPUS_PER_TASK

    estado_qc=$?
    if [ $estado_qc -eq 0 ]; then
      echo      QC con éxito
    else 
      echo      QC falló
    fi
    # Filter
    echo      Ejecutando Filter...
    python  $sqanti_path/sqanti3_filter.py rules \
      --sqanti_class $path_out/qc/muestra/default_${seq}_${plataforma}_classification.txt \
      --filter_gtf $path_out/qc/muestra/default_${seq}_${plataforma}_corrected.gtf \
      --json_filter $sqanti_path/src/utilities/filter/filter_default.json \
      --dir $path_out/filter/rules \
      --output rules_default_${seq}_${plataforma}

    estado_f=$?
    if [ $estado_f -eq 0 ]; then
      echo      Filter con éxito
    else 
      echo      Filter falló
    fi

    # Rescue, primero QC del GTF de referencia contra sí mismo. RefGTF e isoformas mismo archivo
    if [ ! -f   $path_referencia/sqanti_qc/reference_default_classification.txt ]; then

        python $sqanti_path/sqanti3_qc.py \
        --isoforms $gtf \
        --refGTF $gtf \
        --refFasta $fasta \
        --dir $path_referencia/sqanti_qc/ \
        --output reference_default \
        -t $SLURM_CPUS_PER_TASK

      estado_qc2=$?
      if [ $estado_qc2 -eq 0 ]; then
        echo      QC para la referencia con éxito
      fi
      
    else
      echo      QC para la referencia encontrado
    fi

    # Rescue parte 2
    echo      Ejecutando Rescue...
    python $sqanti_path/sqanti3_rescue.py \
      -s rules \
      --filter_class $path_out/filter/rules/rules_default_${seq}_${plataforma}_RulesFilter_classification.txt \
      --refGTF $gtf \
      --refFasta $fasta \
      --refClassif $path_referencia/sqanti_qc/reference_default_classification.txt \
      --mode full \
      --json_filter $sqanti_path/src/utilities/filter/filter_default.json \
      --corrected_isoforms_fasta $path_out/qc/muestra/default_${seq}_${plataforma}_corrected.fasta \
      --filtered_isoforms_gtf $path_out/filter/rules/rules_default_${seq}_${plataforma}.filtered.gtf \
      --dir $path_out/rq/complete_rescue \
      --output rq_${seq}_${plataforma}

      estado_r=$?
      if [ $estado_r -eq 0 ]; then
        echo      Rescue con éxito
      else
        echo      Rescue falló
      fi

  done

  echo
  echo
done
