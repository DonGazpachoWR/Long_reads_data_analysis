#!/usr/bin/env Rscript

# Packages
suppressPackageStartupMessages({
        library(tidyverse)
        library(stringr)
        library(dplyr)
})

# Input Arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
        stop("Uso: Rscript fusionador_correcto_matriz_conteos.r <quantification_fofn> <tama_merge_file> <out_file>")
}

quantification_fofn_path <- args[1]
tama_merge_file          <- args[2]
out_file                 <- args[3]

# 1. Leer FOFN expandiendo la tilde '~' si la hubiera
quantification_fofn <- readLines(quantification_fofn_path, warn = FALSE)
quantification_fofn <- quantification_fofn[trimws(quantification_fofn) != ""]

# 2. Leer tabla TAMA deshabilitando comillas para evitar colapsos
tama_merge_table <- read.table(
        tama_merge_file, 
        header = FALSE, 
        sep = "\t", 
        quote = "", 
        comment.char = "", 
        stringsAsFactors = FALSE
)



# Obtener IDs de muestras limpios
all_sample_ids <- sapply(strsplit(quantification_fofn, "\t"), `[`, 1)
all_sample_ids <- trimws(all_sample_ids)
all_sample_ids <- all_sample_ids[order(nchar(all_sample_ids), decreasing = TRUE)]

# Parsear columna 4 de TAMA
tama_id_map <- bind_rows(lapply(as.character(tama_merge_table[, 4]), function(x) {
        parts  <- strsplit(x, ";", fixed = TRUE)[[1]]
        new_id <- parts[1]
        rest   <- parts[2]
        
        matched_sample <- NA
        old_id         <- NA
        
        for (sid in all_sample_ids) {
                prefix <- paste0(sid, "_")
                if (startsWith(rest, prefix)) {
                        matched_sample <- sid
                        old_id         <- substr(rest, nchar(prefix) + 1, nchar(rest))
                        break
                }
        }
        return(data.frame(new_id = new_id, sample_id = matched_sample, old_id = old_id, stringsAsFactors = FALSE))
}))

# Limpiar comillas parásitas en los IDs parseados de TAMA por seguridad
tama_id_map$old_id <- gsub('"', '', tama_id_map$old_id)

cat("[INFO] Head de tama_id_map:\n")
print(head(tama_id_map))

quant_mat_list <- list()

# 3. Procesar cada archivo de conteo
for (i in seq_along(quantification_fofn)) {
        
        quant_fofn_line <- strsplit(quantification_fofn[i], "\t")[[1]]
        sample_id       <- trimws(quant_fofn_line[1])
        raw_path        <- trimws(quant_fofn_line[2])
        
        # Expandir tilde '~' a ruta absoluta del sistema
        file_path       <- path.expand(raw_path)
        
        if (!file.exists(file_path)) {
                stop(sprintf("ERROR: No se encuentra el archivo de conteo en la ruta: %s", file_path))
        }
        
        # Autodetección de separador (coma o tabulador)
        first_line <- readLines(file_path, n = 1, warn = FALSE)
        sep_char   <- if (grepl(",", first_line)) "," else "\t"
        
        # Leer archivo de conteo deshabilitando el procesamiento de comillas (quote = "")
        quant_file <- read.table(
                file_path, 
                header = TRUE, 
                sep = sep_char, 
                quote = "", 
                comment.char = "", 
                check.names = FALSE, 
                stringsAsFactors = FALSE
        )
        
        # Renombrar primera columna a 'pbid' y limpiar comillas en IDs
        colnames(quant_file)[1] <- "pbid"
        quant_file$pbid <- gsub('"', '', quant_file$pbid)
        
        cat(sprintf("\n[INFO] Muestra: %s | Dimensiones leídas: %d filas x %d columnas\n", 
                    sample_id, nrow(quant_file), ncol(quant_file)))
        
        if (ncol(quant_file) < 2) {
                stop(sprintf("ERROR: El archivo '%s' se ha leído con 1 sola columna. Revisa el delimitador.", file_path))
        }
        
        # Filtrar mapeo TAMA por muestra actual
        tama_name_by_sample <- tama_id_map[!is.na(tama_id_map$sample_id) & tama_id_map$sample_id == sample_id, ]
        
        # Cruzar IDs
        match_ids <- match(quant_file$pbid, tama_name_by_sample$old_id)
        quant_file$transcript_id <- tama_name_by_sample$new_id[match_ids]
        
        n_matched <- sum(!is.na(quant_file$transcript_id))
        cat(sprintf("[INFO] IDs mapeados a TAMA: %d de %d (%.2f%%)\n", 
                    n_matched, nrow(quant_file), (n_matched / nrow(quant_file)) * 100))
        
        # Eliminar columna 'pbid' original
        quant_file <- quant_file %>% select(-pbid)
        
        # Formatear a formato ancho asignando transcript_id
        quant_file_formatted <- quant_file %>%
                filter(!is.na(transcript_id)) %>%
                pivot_longer(-transcript_id, names_to = "sample", values_to = "count") %>%
                group_by(transcript_id, sample) %>%
                summarise(tot_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
                pivot_wider(names_from = "sample", values_from = "tot_count")
        
        quant_mat_list[[i]] <- as.data.frame(quant_file_formatted)
}

# 4. Consolidar matrices de todas las muestras
quant_merged_matrix <- quant_mat_list %>%
        reduce(full_join, by = "transcript_id") %>%
        replace(is.na(.), 0) %>%
        relocate(transcript_id)

names(quant_merged_matrix)[names(quant_merged_matrix) == 'transcript_id'] <- 'superPBID'
# limpiar columnas vacías
quant_merged_matrix <- quant_merged_matrix[, colnames(quant_merged_matrix) != "" & !is.na(colnames(quant_merged_matrix))] %>%
        rename_with(~ str_remove(., "filtrado_(isoseq|masseq)_"))
# 5. Guardar resultado final
write.table(
        quant_merged_matrix, 
        file = out_file,
        sep = "\t", 
        col.names = TRUE, 
        row.names = FALSE, 
        quote = FALSE
)

cat(sprintf("\n[ÉXITO] Matriz final guardada en: %s (%d transcritos procesados)\n", out_file, nrow(quant_merged_matrix)))