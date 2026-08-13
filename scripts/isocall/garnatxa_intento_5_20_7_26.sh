#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 50G
#SBATCH -c 10
#SBATCH -t 5:00:00
#SBATCH -o log/log.out
#SBATCH -e log/erro.out

path=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/isoseq

mkdir -p $path/profile
mkdir -p $path/call
mkdir -p $path/log


mkdir -p $path/filtrado

set -eou pipefail

module load samtools

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
mamba activate isocall

gtf2=/home/adrianbe/practicas/intento1/gtf_sin_randomX.gtf
fasta=$path/arreglar_fasta/genoma_raton_reordenado.fa
bed=$path/regiones_validas.bed

if [[ ! -f $path/ref_seq.isoforms.gz ]]; then

echo Ejecutando prep isoforms...
isocall prep-isoforms \
    --gtf $gtf2 \
    --output ref_seq.isoforms.gz
echo Prep isoforms ejecutado con éxito

fi


echo Ejecutando bams

set +e


for bam in $path/aligned/*.bam;
do
    nombre=$(basename $bam)
    echo 
    echo Procesando bam: $nombre...

    isocall profile \
        --reads $bam \
        --output profile/${nombre}_sample.gz \
        --io-threads 8

    estado_profile=$?

    if [ $estado_profile -eq 0 ]; then

        isocall call \
            --merged-profile profile/${nombre}_sample.gz \
            --known-isoforms $path/ref_seq.isoforms.gz \
            --reference $fasta \
            --output-prefix  call/${nombre}_merged \
            --threads $SLURM_CPUS_PER_TASK

        estado_call=$?

        if [ $estado_call -ne 0 ];then

            echo Fallo en isocall call para el archivo $nombre
            echo Cortando cromosoma X...

            samtools view -b -@ 10 -L $bed $bam >$path/filtrado/filtrado_${nombre}

            bam=$path/filtrado/filtrado_${nombre}

            samtools index -@ 10 $bam

            isocall profile \
                --reads $bam \
                --output profile/filtrado_${nombre}_sample.gz \
                --io-threads 8

            estado_profile=$?

            if [ $estado_profile -eq 0 ]; then

                isocall call \
                    --merged-profile profile/filtrado_${nombre}_sample.gz \
                    --known-isoforms $path/ref_seq.isoforms.gz \
                    --reference $fasta \
                    --output-prefix  call/${nombre}_merged \
                    --threads $SLURM_CPUS_PER_TASK

                estado_call=$?

                if [ $estado_call -ne 0 ];then
                    echo Volvió a fallar tras cortar cromosoma en isocall call para el archivo $nombre
                else
                    echo Éxito tras cortar para el archivo $nombre
                fi
            else
                echo Volvió a fallar tras cortar cromosoma en isocall profile para el archivo $nombre
            fi

        fi

    else
        echo Fallo en isocall profile para el archivo $nombre
    fi

    if [ $estado_call -eq 0 ] && [ $estado_profile -eq 0 ]; then
        echo $archivo procesado con éxito
    fi

done


