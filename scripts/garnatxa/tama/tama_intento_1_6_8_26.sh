#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 20G
#SBATCH -c 1
#SBATCH -t 10:00:00
#SBATCH -o log/log.out
#SBATCH -e log/erro.out

module load samtools

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
mamba activate tama_env

set -eou pipefail

tama_path=/home/adrianbe/practicas/tama
tama_go_path=/home/adrianbe/practicas/tama/tama_go/format_converter


for seq in isoseq masseq;
do
    echo --------------------------------------------
    echo "                   $seq"
    echo --------------------------------------------
    echo
    echo Procesando muestras de $seq ...
    path=/storage/gge/home_members/adrianbe/practicas/intento1/isocall/$seq
    

    mkdir -p $path/tama
    mkdir -p $path/log
    mkdir -p $path/tama/gtf
    mkdir -p $path/tama/bed
    mkdir -p $path/tama/gtf_merged

    out_dir=$path/tama

    > $out_dir/input.txt

    if [ $seq == "isoseq" ]; then
        lista=("brain" "kidney" "b20" "b80")
    else
        lista=("brain" "kidney" "b20")
    fi

    # generar .bed
    for condicion in ${lista[@]};
    do  
        echo
        echo Procesando $condicion...

        gtf_gz=$path/merged_condicion/call/$condicion.isoforms.gtf.gz 
        
        bed=$out_dir/bed/$condicion.isoforms.bed

        gtf=$out_dir/gtf/$condicion.isoforms.gtf

        # Descomprimir .gtf.gz
        if [ ! -f $gtf ]; then

            echo Descomprimiendo $condicion.isoforms.gtf.gz ...

            gunzip -c $gtf_gz > $gtf

            if [ $?-eq 0 ]; then

                echo Descompresión con éxito
            else
                echo Fallo en la descompresión

            fi

        fi

        

        # .gtf a .bed

        echo Transformando .gtf a .bed ...

        python $tama_go_path/tama_format_gtf_to_bed12_ncbi.py \
        $gtf \
        $bed

        estado_tama=$?

        if [ $estado_tama -eq 0 ]; then

            echo GTF a BED con éxito
        else
            echo Fallo en GTF a BED

        fi
        
        # Generar archivo input

        echo -e "$bed\tcapped\t1,1,1\t$condicion" >> $out_dir/input.txt

    done

    # tama merge
    echo
    echo Fusionado .BED de $seq ...

    prefix=$path/tama/gtf_merged/transcriptoma_${seq}_isocall

    python $tama_path/tama_merge.py \
    -f  $out_dir/input.txt \
    -p  $prefix \
    -m 0 \
    -d merge_dup \
    -a 50 \
    -z 50

    estado_tama=$?

    if [ $estado_tama -eq 0 ]; then

        echo Fusión de BEDs con éxito
    else
        echo Fallo en fusión de BEDs

    fi

    # transformar .bed a .gtf
    
    echo Transformando BED a GTF ...
    
    python $tama_go_path/tama_convert_bed_gtf_ensembl_no_cds.py \
    $prefix.bed \
    $prefix.gtf

    estado_tama=$?

    if [ $estado_tama -eq 0 ]; then

        echo BED a GTF con éxito
    else
        echo Fallo en BED a GTF

    fi

done


