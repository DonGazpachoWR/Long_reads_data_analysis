#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 50G
#SBATCH -c 10
#SBATCH -t 5:00:00
#SBATCH -o log/log.out
#SBATCH -e log/erro.out


set -eou pipefail

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
mamba activate isocall

gtf2=/home/adrianbe/practicas/intento1/gtf_sincronizado.gtf
fasta=/home/adrianbe/practicas/intento1/isocall/isoseq/arreglar_fasta/genoma_raton_reordenado.fa

echo Ejecutando prep isoforms...
isocall prep-isoforms \
    --gtf $gtf2 \
    --output ref_seq.isoforms.gz
echo Prep isoforms ejecutado con éxito


echo Ejecutando bams

set +e


for bam in /storage/gge/home_members/adrianbe/practicas/intento1/isocall/isoseq/aligned/*.bam;
do
nombre=$(basename $bam)

echo Procesando bam: $nombre...

isocall profile \
    --reads $bam \
    --output profile/${nombre}_sample.gz \
    --io-threads 8

estado_profile=$?

if [ $estado_profile -eq 0 ]; then

isocall call \
    --merged-profile profile/${nombre}_sample.gz \
    --known-isoforms ref_seq.isoforms.gz \
    --reference $fasta \
    --output-prefix  call/${nombre}_merged \
    --threads $SLURM_CPUS_PER_TASK

estado_call=$?

if [ $estado_call -ne 0 ];then

echo Fallo en isocall call para el archivo $nombre
fi

else
echo Fallo en isocall profile para el archivo $nombre
fi

done


