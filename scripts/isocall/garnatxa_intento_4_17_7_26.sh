#!/bin/bash

# submit options

#SBATCH --qos medium
#SBATCH --mem 50G
#SBATCH -c 10
#SBATCH -t 5:00:00
#SBATCH -o log/log.out
#SBATCH -e log/erro.out


set -eou pipefail

module load samtools

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
mamba activate isocall

gtf2=/home/adrianbe/practicas/intento1/gtf_sincronizado.gtf
fasta=/home/adrianbe/practicas/intento1/isocall/isoseq/arreglar_fasta/genoma_raton_reordenado.fa



samtools view -h -@ 8 /storage/gge/home_members/adrianbe/practicas/intento1/isocall/isoseq/aligned/isoseq_B152.bam | \
  awk '$3 != "chrX" || $4 < 168500000' | \
  samtools view -b -@ 10 - > isoseq_B152.filtrado.bam

samtools index isoseq_B152.filtrado.bam


isocall profile \
    --reads isoseq_B152.filtrado.bam \
    --output profile/isoseq_B152.filtrado_sample.gz \

isocall call \
    --merged-profile profile/isoseq_B152.filtrado_sample.gz \
    --known-isoforms ref_seq.isoforms.gz \
    --reference /home/adrianbe/practicas/intento1/isocall/isoseq/arreglar_fasta/genoma_raton_reordenado.fa \
    --output-prefix call/isoseq_B152.filtrado_merged \


