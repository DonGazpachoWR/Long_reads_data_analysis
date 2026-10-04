#!/usr/bin/env Rscript
# R² and RMSE of observed vs expected tissue mixture for one combination, one
# SQANTI phase and one run (release-6.1 pull requests, kb2 or the no-filter
# baseline). The prevalence filter runs inside SQANTI, so no R expression
# filter is applied here. Same computation as the kb2 analysis
# (antiguo/05_prevalencia_sqanti_kb2), kept here for the release-6.1 runs.
#
# Arguments: ext (class | TPM), seq, platform, mode (raw | qc | fl | rq),
#            data dir, output dir, run label (e.g. r61, kb2, sin_filtro)
#
# Output (stdout, one CSV line per mixture):
# modo_medida,ejecucion,tipo_de_seqs,tipo_de_plataformas,tipo_de_modos,comparacion,
# r_cuadrado,rmse,num_isoformas_modo,num_isoformas_filter_result,num_isoformas_fsm,
# isoform_sin_filtro_r

args <- commandArgs(trailingOnly = TRUE)
ext       <- args[1]
s         <- args[2]
p         <- args[3]
m         <- args[4]
ruta      <- args[5]
outdir    <- args[6]
ejecucion <- args[7]

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpp)
    library(ggpointdensity)
    library(viridis)
    library(dplyr)
})

# Etiquetas en base 10 en ggplot
label_10_pow <- function(x) {
    parse(text = paste0("10^", x))
}

# Columnas según el tipo de secuenciación.
# cols1 (K) y cols2 (B) forman la mezcla computacional predictora;
# cols3 y cols4 son las mezclas observadas.
cols_sel <- function(seq) {
    if (seq == "masseq") {
        cols <- list(
            cols1 = c("K31", "K32", "K33"),
            cols2 = c("B31", "B32", "B33"),
            cols3 = c("B20K80_1", "B20K80_2", "B20K80_3")
        )
    } else {
        cols <- list(
            cols1 = c("K31", "K32", "K33", "K34", "K35"),
            cols2 = c("B31", "B32", "B33", "B34", "B35"),
            cols3 = c("B20K80_1", "B20K80_2", "B20K80_3", "B20K80_4", "B20K80_5"),
            cols4 = c("B80K20_1", "B80K20_2", "B80K20_3", "B80K20_4", "B80K20_5")
        )
    }
    return(cols)
}

# WORKFLOW PASO 1. LIMPIEZA DE NA
NA_to_0 <- function(df, cols) {
    cols <- unlist(cols, use.names = FALSE)
    df[, cols][is.na(df[, cols])] <- 0
    return(df)
}

# Valor de m según el número de réplicas M (mismo cálculo que el antiguo filtro 2)
m_value <- function(M, epsilon = 0.05, alpha = 0.05) {
    all_possible_m <- 1:M
    probabilities <- pbinom(all_possible_m - 1, size = M, prob = epsilon, lower.tail = FALSE)
    valid_m <- all_possible_m[probabilities <= alpha]

    if (length(valid_m) == 0) {
        m_calc <- M
    } else {
        m_calc <- min(valid_m)
    }
    return(m_calc)
}

# Comprobación: filas que NO cumplen el antiguo filtro 2 de R
# (conteo >= 1 en al menos m réplicas de K o de B; mezclas excluidas).
# En fl de kb2 todas las Isoform deberían cumplirlo, así que debe dar 0.
n_sin_filtro_r <- function(df, seq) {
    cols <- cols_sel(seq = seq)
    m_rep <- m_value(M = length(cols$cols1))
    cumple <- rowSums(df[, cols$cols1] >= 1) >= m_rep |
        rowSums(df[, cols$cols2] >= 1) >= m_rep
    return(sum(!cumple))
}

# WORKFLOW PASO 2. FILTRO DE SQANTI: en fl, solo las isoformas clasificadas
# como Isoform. rq ya solo contiene Isoform + rescatadas.
sqanti_filt <- function(df, modo) {
    if (modo == "fl") {
        if (!"filter_result" %in% colnames(df)) {
            stop("La clasificación del filtro no tiene columna filter_result")
        }
        df <- df[which(df[, "filter_result"] == "Isoform"), ]
    }
    return(df)
}

# WORKFLOW PASO 3. CONVERSIÓN A TPM
TPM <- function(df, cols) {
    cols <- unlist(cols, use.names = FALSE)
    a <- colSums(df[, cols], na.rm = TRUE) / 1000000
    a[a == 0] <- 1
    df[, cols] <- sweep(df[, cols], 2, a, FUN = "/")
    return(df)
}

# WORKFLOW PASO 4. MEDIA RÉPLICAS BIOLÓGICAS + MEDIAS TEÓRICAS
bio_replicate_mean <- function(df_modo, name, cols, seq) {
    df_f <- data.frame(
        tr_id = df_modo[, name],
        K100  = rowMeans(df_modo[, cols$cols1]),
        B100  = rowMeans(df_modo[, cols$cols2]),
        B20   = rowMeans(df_modo[, cols$cols3])
    )

    df_f["B20_ex"] <- df_f[, "B100"] * 0.2 + df_f[, "K100"] * 0.8

    if (seq != "masseq") {
        df_f["B80"]    <- rowMeans(df_modo[, cols$cols4])
        df_f["B80_ex"] <- df_f[, "B100"] * 0.8 + df_f[, "K100"] * 0.2
    }
    return(df_f)
}

# Anotación estructural
anotar <- function(df_f, df_modo, modo, combi, extension, ruta) {
    if (modo == "raw") {
        file_qc <- file.path(ruta, paste0("default_", combi, extension))
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
    return(df_f)
}

# WORKFLOW PASO 5. FILTRADO DE ISOFORMAS
isoforms_filt <- function(df, seq) {
    df_filtered <- df[which(
        df[, "associated_transcript"] != "novel" &
            df[, "structural_category"] == "full-splice_match"
    ), ]
    return(df_filtered)
}

# WORKFLOW PASO 6. LOG10 + PSEUDOCOUNT
log10_pseudocount <- function(df, seq) {
    if (seq == "masseq") {
        cols_exp <- c("K100", "B100", "B20", "B20_ex")
    } else {
        cols_exp <- c("K100", "B100", "B20", "B80", "B20_ex", "B80_ex")
    }
    df[, cols_exp] <- log10(df[cols_exp] + 0.01)
    return(df)
}

# Procesamiento global
procesar_datos <- function(modo, seq, plataforma, ruta, ext) {
    combi <- paste(seq, plataforma, sep = "_")
    extension <- "_classification.txt"

    file_data <- switch(modo,
                        "raw" = file.path(ruta, paste0("counts_", combi, ".tsv")),
                        "qc"  = file.path(ruta, paste0("default_", combi, extension)),
                        "fl"  = file.path(ruta, paste0("rules_default_", combi, "_RulesFilter", extension)),
                        "rq"  = file.path(ruta, paste0("rq_", combi, "_rescued", extension))
    )

    name <- ifelse(modo == "raw", "superPBID", "isoform")
    df_modo <- read.table(file_data, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
    n0 <- nrow(df_modo)

    cols <- cols_sel(seq = seq)
    df_sin_na <- NA_to_0(df = df_modo, cols = cols)

    df_sqanti <- sqanti_filt(df = df_sin_na, modo = modo)
    n1 <- nrow(df_sqanti)
    n_r <- if (modo == "fl") n_sin_filtro_r(df = df_sqanti, seq = seq) else NA

    if (ext == "TPM") {
        df_sqanti <- TPM(df = df_sqanti, cols = cols)
    }

    df_bio_repl_mean <- bio_replicate_mean(df_modo = df_sqanti, name = name, cols = cols, seq = seq)
    df_anotado <- anotar(df_f = df_bio_repl_mean, df_modo = df_sqanti, modo = modo,
                         combi = combi, extension = extension, ruta = ruta)

    df_isoformas_filtradas <- isoforms_filt(df = df_anotado, seq = seq)
    n2 <- nrow(df_isoformas_filtradas)

    df_f <- log10_pseudocount(df = df_isoformas_filtradas, seq = seq)

    return(list(df = df_f, n0 = n0, n1 = n1, n2 = n2, n_r = n_r))
}

# Generación de gráficos y cálculo de R² y RMSE
generar_y_guardar_plot <- function(df_data, var_x, var_y, titulo, filename, target_dir) {
    fit <- lm(as.formula(paste(var_y, "~", var_x)), data = df_data)
    r_squared <- summary(fit)$r.squared

    residuos <- df_data[[var_y]] - df_data[[var_x]]
    rmse_val <- sqrt(mean(residuos^2, na.rm = TRUE))

    filepath <- file.path(target_dir, paste0(filename, ".png"))

    if (!file.exists(filepath)) {
        p <- ggplot(df_data, aes(x = .data[[var_x]], y = .data[[var_y]])) +
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
            geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95, color = "red") +
            geom_text_npc(
                data = data.frame(
                    npcx = 0.05,
                    npcy = 0.93,
                    label = paste0("R² = ", round(r_squared, 4), "\nRMSE = ", round(rmse_val, 4))
                ),
                aes(npcx = npcx, npcy = npcy, label = label),
                inherit.aes = FALSE,
                size = 4.2
            ) +
            scale_x_continuous(labels = label_10_pow) +
            scale_y_continuous(labels = label_10_pow) +
            labs(
                x = "Expected expression",
                y = "Observed expression",
                title = titulo
            ) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))

        ggsave(
            filename   = paste0(filename, ".png"),
            plot       = p,
            path       = target_dir,
            create.dir = TRUE,
            device     = "png"
        )
    }

    return(list(r2 = r_squared, rmse = rmse_val))
}

# Una línea del CSV de resultados
linea_csv <- function(mezcla, res_plot, res) {
    cat(paste(
        ifelse(ext == "class", "counts", "TPM"),
        ejecucion, s, p, m, mezcla,
        round(res_plot$r2, 4),
        round(res_plot$rmse, 4),
        res$n0, res$n1, res$n2, res$n_r,
        sep = ","
    ), "\n", sep = "")
}


# ------------------------------------------------------------------------------
# RUTAS DE SALIDA
# ------------------------------------------------------------------------------
dir_destino <- file.path(outdir, "graficos", paste0("_", ext, "_", ejecucion), s, p, m)
dir.create(dir_destino, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# EJECUCIÓN
# ------------------------------------------------------------------------------
res <- tryCatch(
    procesar_datos(m, s, p, ruta, ext),
    error = function(e) {
        message("ERROR en ", ejecucion, " ", ext, " ", s, " ", p, " ", m, ": ", conditionMessage(e))
        NULL
    }
)

if (!is.null(res) && nrow(res$df) > 0) {
    combi_name <- paste(m, s, p, sep = "_")

    # Evaluación B20K80
    res_b20 <- generar_y_guardar_plot(
        df_data    = res$df,
        var_x      = "B20_ex",
        var_y      = "B20",
        titulo     = toupper(paste(ext, ejecucion, m, s, p, "- B20K80")),
        filename   = paste0(combi_name, "_B20K80"),
        target_dir = dir_destino
    )
    linea_csv("B20K80", res_b20, res)

    # Evaluación B80K20
    if (s != "masseq") {
        res_b80 <- generar_y_guardar_plot(
            df_data    = res$df,
            var_x      = "B80_ex",
            var_y      = "B80",
            titulo     = toupper(paste(ext, ejecucion, m, s, p, "- B80K20")),
            filename   = paste0(combi_name, "_B80K20"),
            target_dir = dir_destino
        )
        linea_csv("B80K20", res_b80, res)
    }
}
