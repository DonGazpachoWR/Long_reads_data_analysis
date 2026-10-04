library(dplyr)
library(readr)
library(stringr)
outdir <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/practica/datos"
# Isoseq
brain   <- read.table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/brain.count_matrix.txt", sep = ",", header = T)
kidney <- read.table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/kidney.count_matrix.txt", sep = ",", header = T)
b20 <- read.table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/b20.count_matrix.txt", sep = ",", header = T)
b80 <- read.table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/b80.count_matrix.txt", sep = ",", header = T)

df <- brain %>%
        full_join(
                kidney %>% select(id, filtrado_isoseq_K31, filtrado_isoseq_K32,
                                 filtrado_isoseq_K33, filtrado_isoseq_K34, filtrado_isoseq_K35),
                by = c("id")
        ) %>% 
        full_join( 
                b20 %>% select(id, filtrado_isoseq_B20K80_1, filtrado_isoseq_B20K80_2,
                   filtrado_isoseq_B20K80_3, filtrado_isoseq_B20K80_4, filtrado_isoseq_B20K80_5),
               by = c("id") 
        ) %>% 
        full_join(
                b80 %>% select(id, filtrado_isoseq_B80K20_1, filtrado_isoseq_B80K20_2,
               filtrado_isoseq_B80K20_3, filtrado_isoseq_B80K20_4, filtrado_isoseq_B80K20_5),
               by = c("id")
        ) %>%
        rename_with(~ str_remove(., "filtrado_isoseq_"))
# En realidad no es "superPBID" sino "associated transcript", pero lo dejamos así para el siguiente codigo de R
colnames(df)[1] <- "superPBID"

ruta_tsv <- file.path(outdir, "counts_isoseq_isocall.tsv")
write_tsv(df, file = ruta_tsv)

# Masseq

brain   <- read.table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/masseq/brain.count_matrix.txt", sep = ",", header = T)
kidney <- read.table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/masseq/kidney.count_matrix.txt", sep = ",", header = T)
b20 <- read.table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/masseq/b20.count_matrix.txt", sep = ",", header = T)

df <- brain %>%
        full_join(
                kidney %>% select(id, filtrado_masseq_K31, filtrado_masseq_K32,
                                  filtrado_masseq_K33),
                by = c("id")
        ) %>% 
        full_join( 
                b20 %>% select(id, filtrado_masseq_B20K80_1, filtrado_masseq_B20K80_2,
                               filtrado_masseq_B20K80_3, ),
                by = c("id") 
        )  %>%
        rename_with(~ str_remove(., "filtrado_masseq_"))
# En realidad no es "superPBID" sino "associated transcript", pero lo dejamos así para el siguiente codigo de R
colnames(df)[1] <- "superPBID"

ruta_tsv <- file.path(outdir, "counts_masseq_isocall.tsv")
write_tsv(df, file = ruta_tsv)
