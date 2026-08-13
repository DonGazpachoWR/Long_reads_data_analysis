#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 50G
#SBATCH -c 10
#SBATCH -t 10:00:00
#SBATCH -o log/log.out
#SBATCH -e log/erro.out

module load samtools

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
mamba activate isocall

set -eou pipefail

gtf2=/home/adrianbe/practicas/intento1/gtf_sincronizado.gtf
fasta=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/isoseq/arreglar_fasta/genoma_raton_reordenado.fa
bed=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/isoseq/regiones_validas.bed

espacio1="  "
espacio2="       "

for seq in isoseq masseq;
do
    echo --------------------------------------------
    echo "                   $seq"
    echo --------------------------------------------
    echo
    echo Procesando muestras de $seq ...
    path=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/$seq
    path_bams=$path

    mkdir -p $path/merged_condicion
    path=$path/merged_condicion

    mkdir -p $path/profile
    mkdir -p $path/call
    mkdir -p $path/log
    mkdir -p $path_bams/filtrado
    mkdir -p $path/merged
    mkdir -p $path/call

    if [[ ! -f $path/ref_seq.isoforms.gz ]]; then

        echo $espacio1 Ejecutando prep isoforms...
        isocall prep-isoforms \
            --gtf $gtf2 \
            --output $path/ref_seq.isoforms.gz
        echo $espacio1 Prep isoforms ejecutado con éxito

    fi

    brain=()
    kidney=()
    b20=()
    b80=()


    echo Ejecutando bams

    #set +e

    for bam in $path_bams/aligned/*.bam;
    do
        nombre=$(basename $bam)
        echo 
        echo $espacio1 Procesando bam: $nombre...

        # Cortamos todos los archivos, si no ha sido cortado lo cortamos
        if [ ! -f $path_bams/filtrado/filtrado_${nombre} ]; then

            echo $espacio2 Cortando cromosoma X...

            samtools view -b -@ 10 -L $bed $bam > $path_bams/filtrado/filtrado_${nombre}

            bam=$path_bams/filtrado/filtrado_${nombre}

            samtools index -@ 10 $bam
        else
            bam=$path_bams/filtrado/filtrado_${nombre}

        fi


        isocall profile \
            --reads $bam \
            --output $path/profile/filtrado_${nombre}_sample.gz \
            --io-threads $SLURM_CPUS_PER_TASK

        estado_profile=$?
        
        if [ $estado_profile -ne 0 ]; then
            echo $espacio2 Fallo tras cortar cromosoma en isocall profile para el archivo $nombre
        
        else
            echo $espacio2 Profile ejecutado con éxito para $nombre
            path_profile=$path/profile/filtrado_${nombre}_sample.gz
            if [[ $nombre =~ B80 ]];
            then
                b80+=("$path_profile")
            elif [[ $nombre =~ B20 ]];
            then
                b20+=("$path_profile")
            elif [[ $nombre =~ B ]];
            then
                brain+=("$path_profile")
            elif [[ $nombre =~ K ]];
            then
                kidney+=("$path_profile")
        fi

        fi
        # Combinado de muestras

    done
    echo Muestras clasificadas en las siguientes listas:
    echo
    echo brain: ${brain[@]}
    echo kidney: ${kidney[@]}
    echo B20_K80: ${b20[@]}
    echo B80_K20: ${b80[@]}
    echo
    for lista in brain kidney b20 b80;
    do
        echo
        echo $espacio1 Procesando lista $lista...
        estado_call=1

        declare -n muestras_grupo=$lista
        
        isocall merge \
            --profiles ${muestras_grupo[@]} \
            --output $path/merged/$lista.gz
        
        estado_merge=$?

        if [ $estado_merge -eq 0 ]; then

            echo $espacio1 Éxito en isocall merge para $lista

            isocall call \
                --merged-profile $path/merged/$lista.gz \
                --known-isoforms $path/ref_seq.isoforms.gz \
                --reference $fasta \
                --output-prefix  $path/call/$lista \
                --threads $SLURM_CPUS_PER_TASK
            
            estado_call=$?

            if [ $estado_call -ne 0 ];then
                    echo $espacio1 Fallo en isocall call para archivos de $lista
                else
                    echo $espacio1 Éxito para isocall para archivos de $lista
            fi
        else 
            echo $espacio1 Fallo en isocall merge para $lista
        fi

        if [ $estado_merge -eq 0 ] && [ $estado_call -eq 0 ]; then
            echo $espacio1 Éxito total al procesar $lista
        else
            echo $espacio1 Fallo total al procesar $lista
        fi

    done

        
done


