#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 50G
#SBATCH -c 10
#SBATCH -t 10:00:00
#SBATCH -o log/log.out
#SBATCH -e log/erro.out

# ESTE SCRIPT ES SIN COMBIAR MUESTRAS. SE PRUEBA A QUE FUNCIONE TODO ANTES DEL INTENTO 9 EN QUE YA SE COMBINAN LAS MUESTRAS
module load samtools

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
mamba activate isocall

set -eou pipefail

gtf2=/home/adrianbe/practicas/intento1/gtf_sincronizado.gtf
fasta=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/isoseq/arreglar_fasta/genoma_raton_reordenado.fa
bed=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/isoseq/regiones_validas.bed

for seq in isoseq masseq;
do
    echo Procesando muestras de $seq ...
    path=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/$seq

    mkdir -p $path/profile
    mkdir -p $path/call
    mkdir -p $path/log
    mkdir -p $path/filtrado

    if [[ ! -f $path/ref_seq.isoforms.gz ]]; then

        echo Ejecutando prep isoforms...
        isocall prep-isoforms \
            --gtf $gtf2 \
            --output $path/ref_seq.isoforms.gz
        echo Prep isoforms ejecutado con éxito

    fi

    echo Ejecutando bams

    set +e

    for bam in $path/aligned/*.bam;
    do
        nombre=$(basename $bam)
        echo 
        echo Procesando bam: $nombre...

        echo Cortando cromosoma X...

        samtools view -b -@ 10 -L $bed $bam >$path/filtrado/filtrado_${nombre}

        bam=$path/filtrado/filtrado_${nombre}

        samtools index -@ 10 $bam

        isocall profile \
            --reads $bam \
            --output $path/profile/filtrado_${nombre}_sample.gz \
            --io-threads $SLURM_CPUS_PER_TASK

        estado_profile=$?
        estado_call=1

        if [ $estado_profile -eq 0 ]; then

            isocall call \
                --merged-profile $path/profile/filtrado_${nombre}_sample.gz \
                --known-isoforms $path/ref_seq.isoforms.gz \
                --reference $fasta \
                --output-prefix  $path/call/${nombre}_merged \
                --threads $SLURM_CPUS_PER_TASK

            estado_call=$?

            if [ $estado_call -ne 0 ];then
                echo Fallo tras cortar cromosoma en isocall call para el archivo $nombre
            else
                echo Éxito tras cortar para el archivo $nombre
            fi
        else
            echo Fallo tras cortar cromosoma en isocall profile para el archivo $nombre
        fi

        if [ $estado_call -eq 0 ] && [ $estado_profile -eq 0 ]; then
            echo $nombre procesado con éxito
        else
            echo $nombre no procesado
        fi

    done

done


