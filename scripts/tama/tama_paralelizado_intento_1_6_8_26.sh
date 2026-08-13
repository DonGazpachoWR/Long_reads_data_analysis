#!/bin/bash

set -eou pipefail

tama_path=/home/adrianbe/practicas/tama
tama_go_path=/home/adrianbe/practicas/tama/tama_go/format_converter
python_path=/storage/gge/home_members/adrianbe/.conda/envs/tama_env/bin/python

for seq in isoseq masseq; do
    echo "--------------------------------------------"
    echo "            Preparando $seq                 "
    echo "--------------------------------------------"


    path=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/$seq
    out_dir=$path/tama
    mkdir -p $path/tama $path/log $path/tama/gtf $path/tama/bed $path/tama/gtf_merged

    > $out_dir/input.txt

    if [ "$seq" == "isoseq" ]; then
        lista=("brain" "kidney" "b20" "b80")
    else
        lista=("brain" "kidney" "b20")
    fi

    # Array para guardar los Job IDs de las descompresiones de esta secuencia
    job_ids=()

    # Descomprimir .gtf.gz y transformar a .bed

    for condicion in "${lista[@]}"; do

        gtf_gz=$path/merged_condicion/call/$condicion.isoforms.gtf.gz 
        gtf=$out_dir/gtf/$condicion.isoforms.gtf
        bed=$out_dir/bed/$condicion.isoforms.bed

        echo -e "$bed\tcapped\t1,1,1\t$condicion" >> $out_dir/input.txt

        if [ -f $bed ]; then
            echo Descompresión y conversión YA REALIZADA para: $seq - $condicion
        
        else 

            echo Job de descompresión y conversión para: $seq - $condicion...

            # Capturar el Job ID devuelto por Slurm (--parsable)
            JOB_ID=$(sbatch --parsable \
                --qos=medium \
                --mem=10G \
                -c 1 \
                -t 01:00:00 \
                -o "$path/log/prep_${condicion}.out" \
                -e "$path/log/prep_${condicion}.err" \
                --wrap="
                    module load samtools
                    module load anaconda
                    source \$(conda info --base)/etc/profile.d/conda.sh
                    mamba activate tama_env
                    set -eou pipefail

                    if [ ! -f $gtf ]; then
                        echo 'Descomprimiendo $gtf_gz...'
                        gunzip -c $gtf_gz > $gtf
                    else
                        echo 'El archivo $gtf ya existe.'
                    fi

                    echo 'Transformando GTF a BED...'
                    if [ ! -f $bed ]; then
                        $python_path $tama_go_path/tama_format_gtf_to_bed12_ncbi.py $gtf $bed
                    fi
                ")

            job_ids+=("$JOB_ID")
        fi


    done

    # Crear lista de dependencias separada por dos puntos (ej. 12345:12346:12347)
    DEP_LIST=$(IFS=:; echo "${job_ids[*]}")

    # Tama merge + .bed a .gtf
    prefix=$path/tama/gtf_merged/transcriptoma_${seq}_isocall

    echo "Job de fusión ($seq) con dependencia de Jobs: $DEP_LIST"

    sbatch --dependency=afterok:${DEP_LIST} \
        --qos=medium \
        --mem=20G \
        -c 1 \
        -t 04:00:00 \
        -o "$path/log/merge_${seq}.out" \
        -e "$path/log/merge_${seq}.err" \
        --wrap="
            module load samtools
            module load anaconda
            source \$(conda info --base)/etc/profile.d/conda.sh
            mamba activate tama_env
            set -eou pipefail

            echo 'Fusionando .BED de $seq...'
            $python_path $tama_path/tama_merge.py \
                -f $out_dir/input.txt \
                -p $prefix \
                -m 0 \
                -d merge_dup \
                -a 50 \
                -z 50

            echo 'Transformando BED a GTF...'
            $python_path $tama_go_path/tama_convert_bed_gtf_ensembl_no_cds.py ${prefix}.bed ${prefix}.gtf
        "

    echo "Pipeline enviado con éxito para $seq."
    echo
done