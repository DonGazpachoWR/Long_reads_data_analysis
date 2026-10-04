#!/usr/bin/env bash

# Crear archivo temporal con los parámetros válidos
PARAM_FILE=$(mktemp)

for ext in class TPM; do
  for fil in 0 1 2 3; do
    for s in isoseq masseq ont; do
      if [ "$s" == "ont" ]; then
        plats="bambu flair isoquant"
      else
        plats="isoseq isocall bambu flair isoquant"
      fi
      for p in $plats; do
        for m in raw qc fl rq; do
          echo "$ext $fil $s $p $m" >> "$PARAM_FILE"
        done
      done
    done
  done
done

# Ejecutar en paralelo con 8 núcleos y guardar salida CSV
HEADER="modo_medida,filtro,tipo_de_seqs,tipo_de_plataformas,tipo_de_modos,comparacion,r_cuadrado,num_isoformas_sin_filtro,num_isoformas_filtro_expresion,num_isoformas_filtro_isoformas"

path=/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/output

echo "$HEADER" > $path/resultados_globales.csv

cat "$PARAM_FILE" | parallel --colsep ' ' -j 10 Rscript paralelizado_dotplot_individual.r {1} {2} {3} {4} {5} >> $path/resultados_globales.csv

rm "$PARAM_FILE"
