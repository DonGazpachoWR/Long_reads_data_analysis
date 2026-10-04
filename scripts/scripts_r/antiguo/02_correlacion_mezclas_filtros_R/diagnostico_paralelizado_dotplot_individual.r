#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
ext  <- args[1]
fil  <- as.numeric(args[2])
s    <- args[3]
p    <- args[4]
m    <- args[5]

tag <- sprintf("[%s | fil:%s | %s | %s | %s]", ext, fil, s, p, m)
cat(sprintf("\n=== INICIO %s ===\n", tag))

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpp)
    library(ggpointdensity)
    library(viridis)
    library(dplyr)
})

ruta   <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/data/data_expression_matrix"
outdir <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/output/expression_matrix/graficos"

# --- Funciones auxiliares ---
label_10_pow <- function(x) parse(text = paste0("10^", x))

cols_sel <- function(seq){
    if (seq == "masseq"){
        return(list(
            cols1 = c("K31", "K32", "K33"),
            cols2 = c("B31", "B32", "B33"),
            cols3 = c("B20K80_1", "B20K80_2", "B20K80_3")
        ))
    } else {
        return(list(
            cols1 = c("K31", "K32", "K33", "K34", "K35"),
            cols2 = c("B31", "B32", "B33", "B34", "B35"),
            cols3 = c("B20K80_1", "B20K80_2", "B20K80_3", "B20K80_4", "B20K80_5"),
            cols4 = c("B80K20_1", "B80K20_2", "B80K20_3", "B80K20_4", "B80K20_5")
        ))
    }
}

# 1. Transformar a TPM con sustitución previa de NAs y protección contra sumas NA/0
TPM <- function(df, cols){
    cols <- unlist(cols, use.names = FALSE)
    cols_pres <- intersect(cols, colnames(df))
    cat(sprintf("%s [TPM] Columnas evaluadas: %d/%d\n", tag, length(cols_pres), length(cols)))
    
    nas_antes <- sum(is.na(df[, cols_pres]))
    cat(sprintf("%s [TPM] NAs en conteos antes de imputar a cero: %d\n", tag, nas_antes))
    
    # Imputar NAs a 0 en las columnas de expresión
    df[, cols_pres][is.na(df[, cols_pres])] <- 0
    
    # Suma por columnas en millones con na.rm = TRUE
    a <- colSums(df[, cols_pres], na.rm = TRUE) / 1000000
    cat(sprintf("%s [TPM] Suma por columnas (millones): %s\n", tag, paste(round(a, 4), collapse = ", ")))
    
    # Evitar división por cero si alguna columna no tiene lecturas
    a[a == 0] <- 1
    
    df[, cols_pres] <- sweep(df[, cols_pres], 2, a, FUN = "/")
    return(df)
}

# 2. Filtrado de expresión asegurando imputación previa de NAs
min_exprs_filt <- function(df, seq, filtro = 2, ext){
    if (!(filtro %in% c(0, 1, 2, 3))) filtro <- 2
    cols <- cols_sel(seq = seq)
    
    cols_all <- unlist(cols, use.names = FALSE)
    cols_pres <- intersect(cols_all, colnames(df))
    
    # Imputar NAs a 0 antes de cualquier filtrado o normalización
    nas_conteos <- sum(is.na(df[, cols_pres]))
    df[, cols_pres][is.na(df[, cols_pres])] <- 0
    cat(sprintf("%s [min_exprs_filt] NAs encontrados e imputados a cero: %d\n", tag, nas_conteos))
    
    if (ext == "TPM") {
        df <- TPM(df, cols)
    }
    
    filas_antes <- nrow(df)
    
    if (filtro != 0) {
        if (filtro %in% c(1, 2)){
            if (filtro == 1) { n <- 1; m <- 1 }
            else if (filtro == 2) { n <- 1; m <- 2 }
            
            df <- df[which(
                rowSums(df[, cols$cols1] >= n) >= m | 
                    rowSums(df[, cols$cols2] >= n) >= m
            ), ]
        } else if (filtro == 3) {
            if (seq == "masseq"){
                df <- df[which(
                    rowSums(df[, cols$cols1] > 0) > 0 & 
                        rowSums(df[, cols$cols2] > 0) > 0 &
                        rowSums(df[, cols$cols3] > 0) > 0
                ), ]
            } else {
                df <- df[which(
                    rowSums(df[, cols$cols1] > 0) > 0 & 
                        rowSums(df[, cols$cols2] > 0) > 0 &
                        rowSums(df[, cols$cols3] > 0) > 0 &
                        rowSums(df[, cols$cols4] > 0) > 0
                ), ]
            }
        }
    }
    cat(sprintf("%s [min_exprs_filt] Filas: %d -> %d tras filtro %d\n", tag, filas_antes, nrow(df), filtro))
    return(df)
}

isoforms_filt <- function(df, seq) {
    if (seq == "masseq"){ 
        cols_exp <- c("K100", "B100", "B20", "B20_ex") 
    } else { 
        cols_exp <- c("K100", "B100", "B20", "B80", "B20_ex", "B80_ex") 
    }
    
    filas_antes <- nrow(df)
    df_filtered <- df[which(
        df[, "associated_transcript"] != "novel" &
            df[, "structural_category"] == "full-splice_match"
    ), ]
    cat(sprintf("%s [isoforms_filt] Filas FSM y no-novel: %d -> %d\n", tag, filas_antes, nrow(df_filtered)))
    
    if (nrow(df_filtered) > 0) {
        df_filtered[, cols_exp] <- log10(df_filtered[cols_exp] + 0.01)
    }
    return(df_filtered)
}

procesar_datos <- function(modo, seq, plataforma, ruta, filtro, ext) {
    combi <- paste(seq, plataforma, sep = "_")
    extension_sq <- "_classification.txt"
    
    file_data <- switch(modo,
                        "raw" = file.path(ruta, paste0("counts_", combi, ".tsv")),
                        "qc"  = file.path(ruta, paste0("default_", combi, extension_sq)),
                        "fl"  = file.path(ruta, paste0("rules_default_", combi, "_RulesFilter", extension_sq)),
                        "rq"  = file.path(ruta, paste0("rq_", combi, "_rescued", extension_sq))
    )
    
    cat(sprintf("%s [Lectura] Archivo: %s (Existe: %s)\n", tag, basename(file_data), file.exists(file_data)))
    if (!file.exists(file_data)) return(NULL)
    
    name <- ifelse(modo == "raw", "superPBID", "isoform")
    df_modo <- read.table(file_data, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
    cat(sprintf("%s [Lectura] Leídas %d filas y %d columnas\n", tag, nrow(df_modo), ncol(df_modo)))
    
    df_modo <- min_exprs_filt(df_modo, seq, filtro, ext)
    if (is.null(df_modo) || nrow(df_modo) == 0) {
        cat(sprintf("%s [ALERTA] 0 filas tras min_exprs_filt\n", tag))
        return(NULL)
    }
    
    cols <- cols_sel(seq = seq)
    df_f <- data.frame(
        tr_id = df_modo[, name],
        K100  = rowMeans(df_modo[, cols$cols1]),
        B100  = rowMeans(df_modo[, cols$cols2]),
        B20   = rowMeans(df_modo[, cols$cols3])
    )
    df_f["B20_ex"] <- df_f[, "B100"] * 0.2 + df_f[, "K100"] * 0.8
    
    if (seq != "masseq"){
        df_f["B80"] <- rowMeans(df_modo[, cols$cols4])
        df_f["B80_ex"] <- df_f[, "B100"] * 0.8 + df_f[, "K100"] * 0.2
    }
    
    if (modo == "raw"){
        file_qc <- file.path(ruta, paste0("default_", combi, extension_sq))
        df_qc <- read.table(file_qc, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
        df_f <- df_f %>%
            left_join(
                df_qc %>% select(isoform, associated_transcript, structural_category),
                by = c("tr_id" = "isoform")
            )
    } else {
        df_f["associated_transcript"] <- df_modo[, "associated_transcript"]
        df_f["structural_category"]   <- df_modo[, "structural_category"]
    }
    
    df_f <- isoforms_filt(df_f, seq)
    return(df_f)
}

generar_y_guardar_plot <- function(df_data, var_x, var_y, titulo, filename, target_dir) {
    x_vals <- df_data[[var_x]]
    y_vals <- df_data[[var_y]]
    
    casos_finitos <- is.finite(x_vals) & is.finite(y_vals)
    n_validos <- sum(casos_finitos)
    
    cat(sprintf("%s [Diagnóstico Plot: %s]\n", tag, filename))
    cat(sprintf("  - Total filas df: %d\n", nrow(df_data)))
    cat(sprintf("  - Casos válidos para lm (%s vs %s): %d\n", var_x, var_y, n_validos))
    cat(sprintf("  - %s -> NAs: %d | Infs: %d | NaNs: %d\n", var_x, sum(is.na(x_vals)), sum(is.infinite(x_vals)), sum(is.nan(x_vals))))
    cat(sprintf("  - %s -> NAs: %d | Infs: %d | NaNs: %d\n", var_y, sum(is.na(y_vals)), sum(is.infinite(y_vals)), sum(is.nan(y_vals))))
    
    if (n_validos < 2) {
        cat(sprintf("%s [DETENIDO ANTES DE LM] Menos de 2 casos válidos. Abortando ajuste.\n", tag))
        return(NA)
    }
    
    df_validos <- df_data[casos_finitos, ]
    fit <- lm(as.formula(paste(var_y, "~", var_x)), data = df_validos)
    r_squared <- summary(fit)$r.squared
    
    filepath <- file.path(target_dir, paste0(filename, ".png"))
    if (!file.exists(filepath)) {
        p <- ggplot(df_validos, aes(x = .data[[var_x]], y = .data[[var_y]])) +
            geom_pointdensity(size = 0.8) +
            scale_color_viridis_c(
                name = "Point Density",
                guide = guide_colorbar(
                    barheight = unit(0.8, "npc"),
                    barwidth  = unit(1.2, "lines"),
                    title.position = "right",
                    title.theme    = element_text(angle = -90, hjust = 0.5)
                )
            ) +
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
            geom_text_npc(
                data = data.frame(npcx = 0.05, npcy = 0.95, label = paste0("R² = ", round(r_squared, 4))),
                aes(npcx = npcx, npcy = npcy, label = label),
                inherit.aes = FALSE,
                size = 4.5
            ) +
            scale_x_continuous(labels = label_10_pow) +
            scale_y_continuous(labels = label_10_pow) +
            labs(x = "Expected expression", y = "Observed expression", title = titulo) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        ggsave(filename = paste0(filename, ".png"), plot = p, path = target_dir, create.dir = TRUE, device = "png")
    }
    return(r_squared)
}

# --- Ejecución ---
outdir2     <- file.path(outdir, paste0("_", ext, "_filtro_", fil))
dir_destino <- file.path(outdir2, s, p, m)
dir.create(dir_destino, recursive = TRUE, showWarnings = FALSE)

df_proc <- procesar_datos(m, s, p, ruta, fil, ext)

if (!is.null(df_proc) && nrow(df_proc) > 0) {
    combi_name <- paste(m, s, p, sep = "_")
    r2_b20 <- generar_y_guardar_plot(df_proc, "B20_ex", "B20", paste(toupper(m), toupper(s), toupper(p), "- B20K80"), paste0(combi_name, "_B20K80"), dir_destino)
    cat(paste(ifelse(ext == "class", "counts", "TPM"), fil, s, p, m, "B20K80", round(r2_b20, 4), sep = ","), "\n")
    
    if (s != "masseq") {
        r2_b80 <- generar_y_guardar_plot(df_proc, "B80_ex", "B80", paste(toupper(m), toupper(s), toupper(p), "- B80K20"), paste0(combi_name, "_B80K20"), dir_destino)
        cat(paste(ifelse(ext == "class", "counts", "TPM"), fil, s, p, m, "B80K20", round(r2_b80, 4), sep = ","), "\n")
    }
} else {
    cat(sprintf("%s [FIN CON ERROR] df_proc es NULL o tiene 0 filas.\n", tag))
}
cat(sprintf("=== FIN %s ===\n\n", tag))