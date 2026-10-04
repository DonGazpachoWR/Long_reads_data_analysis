library(readr)
df_isoseq <- data.frame(condicion=c("brain", "kidney", "b20", "b80"),
                 ruta=c("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/brain.count_matrix.txt",
                        "~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/kidney.count_matrix.txt",
                        "~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/b20.count_matrix.txt",
                        "~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/isoseq/b80.count_matrix.txt"))

write_tsv(df_isoseq, "/home/adrian/Documentos/Conesa_Lab/VSCODE/practica/datos_tama/isoseq.tsv", col_names = FALSE)

df_masseq <- data.frame(condicion=c("brain", "kidney", "b20"),
                        ruta=c("~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/masseq/brain.count_matrix.txt",
                               "~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/masseq/kidney.count_matrix.txt",
                               "~/Documentos/Conesa_Lab/VSCODE/practica/datos_isocall/masseq/b20.count_matrix.txt"))
write_tsv(df_masseq, "/home/adrian/Documentos/Conesa_Lab/VSCODE/practica/datos_tama/masseq.tsv", col_names = FALSE)


merged_isoseq <- read_table("~/Documentos/Conesa_Lab/VSCODE/practica/datos_tama/transcriptoma_isoseq_isocall_merge.txt")
