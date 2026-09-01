#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
ext  <- args[1]
fil  <- as.numeric(args[2])
s    <- args[3]
p    <- args[4]
m    <- args[5]

# Cargar librerías y definir funciones 
suppressPackageStartupMessages({
    library(ggplot2); library(ggpp); library(ggpointdensity); library(viridis); library(dplyr)
})

# etiquetas en base 10 en ggplot
label_10_pow <- function(x) {
    parse(text = paste0("10^", x))
}

# obtiene las columnas segun tipo de plataforma de secuenciacion
cols_sel <- function(seq){
    if (seq == "masseq"){
        return ( list( 
            cols1 = c("K31", "K32", "K33"),
            cols2 = c("B31", "B32", "B33"),
            cols3 = c("B20K80_1", "B20K80_2","B20K80_3") 
        ))
    }
    
    else {
        return ( list( 
            cols1 = c("K31", "K32", "K33", "K34", "K35"),
            cols2 = c("B31", "B32", "B33", "B34", "B35") ,
            cols3 = c("B20K80_1", "B20K80_2","B20K80_3", "B20K80_4", "B20K80_5") ,
            cols4 = c("B80K20_1", "B80K20_2","B80K20_3", "B80K20_4", "B80K20_5")
        ))
    }
}

# ---------------------------------------------------------------------------
# FUNCIONES DE MONITORIZACIÓN (nº de filas y categorías estructurales SQANTI)
# ---------------------------------------------------------------------------

# Resume el nº de transcritos por categoría estructural de SQANTI en un 
# string compacto "categoria:n;categoria:n;..." para poder guardarlo en 
# una única celda del csv de resultados
cat_estructural_resumen <- function(df) {
    if (!("structural_category" %in% colnames(df)) || nrow(df) == 0) return("NA")
    tab <- table(df[, "structural_category"])
    paste(paste0(names(tab), ":", as.integer(tab)), collapse = ";")
}

# Guarda en disco (carpeta nueva, separada de la de gráficos) la tabla de 
# datos ya filtrada (tras el WORKFLOW PASO 5), para poder inspeccionarla 
# o reutilizarla sin tener que rehacer todo el filtrado
guardar_tabla_filtrada <- function(df, dir_tablas, combi_name) {
    dir.create(dir_tablas, recursive = TRUE, showWarnings = FALSE)
    filepath <- file.path(dir_tablas, paste0(combi_name, "_filtrado.tsv"))
    write.table(df, file = filepath, sep = "\t", quote = FALSE, row.names = FALSE)
}

# WORKFLOW PASO 0. AÑADIR CATEGORÍA ESTRUCTURAL
# Se hace ANTES de cualquier filtrado (en vez de al final, como antes) 
# para poder monitorizar cómo cambia la composición de categorías 
# estructurales tras cada filtro independiente.
anadir_categoria_estructural <- function(df_modo, modo, combi, extension, ruta, name) {
    
    if (modo == "raw") {
        # El fichero de counts (raw) no trae structural_category ni 
        # associated_transcript: se obtienen uniendo con el fichero de 
        # clasificación por defecto (QC)
        file_qc <- file.path(ruta, paste0("default_", combi, extension))
        df_qc   <- read.table(file_qc, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
        
        df_modo <- df_modo %>%
            left_join(
                df_qc %>% select(isoform, associated_transcript, structural_category),
                by = setNames("isoform", name)
            )
    }
    # Para "qc", "fl" y "rq" ambas columnas ya vienen de forma nativa en df_modo
    
    return(df_modo)
}

# WORKFLOW PASO 1. LIMPIEZA DE NA
NA_to_0 <- function(df, cols){
    cols <- unlist(cols, use.names = FALSE)
    
    
    df[, cols][is.na(df[, cols])] <- 0
    
    return(df)
}

# ---------------------------------------------------------------------------
# WORKFLOW PASO 2. FILTRO DE EXPRESION MINIMA
# Se separa en dos bloques independientes y monitorizables:
#   Bloque 1 -> expresión en ALGUNA de las dos condiciones predictoras
#   Bloque 2 -> expresión en al menos 1 (filtro=1) o "m" (filtro=2) 
#               muestras por condición
# ---------------------------------------------------------------------------

# BLOQUE 1: expresión en alguna de las dos condiciones predictoras
# (al menos 1 lectura en al menos 1 muestra de la condición K100 O de 
# la condición B100)
filt_expr_alguna_condicion <- function(df, seq) {
    cols <- cols_sel(seq = seq)
    
    df[which(
        rowSums(df[, cols$cols1] >= 1) >= 1 |
            rowSums(df[, cols$cols2] >= 1) >= 1
    ), ]
}

# Calcula el nº mínimo de muestras "m" (de un total de n_samples por 
# condición) en las que hay que superar un umbral de "n_reads" lecturas 
# para considerar que la expresión observada no es ruido de fondo.
#
# Fundamento estadístico (Binomial Negativa + Binomial):
#   1) El ruido de fondo (transcritos no expresados / cuantificación 
#      espuria) se modela como una distribución Binomial Negativa 
#      NB(mu0, size0), cuyos parámetros se estiman por el método de los 
#      momentos a partir de las propias lecturas observadas:
#           mu0   = media(X)
#           var0  = varianza(X)
#           size0 = mu0^2 / (var0 - mu0)      [pues Var(NB) = mu + mu^2/size]
#      Si var0 <= mu0 (no hay sobredispersión) se toma size0 muy grande, 
#      lo que hace que la NB colapse a una Poisson(mu0) como caso límite.
#   2) Bajo ese modelo de ruido, la probabilidad de que UNA muestra 
#      cualquiera alcance por puro azar al menos "n_reads" lecturas es:
#           p = P(X >= n_reads) = 1 - pnbinom(n_reads - 1, mu = mu0, size = size0)
#   3) Si las n_samples muestras de la condición son independientes, el 
#      nº de muestras que superan el umbral por azar sigue una 
#      Binomial(n_samples, p).
#   4) "m" se define como el menor entero tal que 
#           P( Binomial(n_samples, p) >= m ) <= alpha
#      es decir, el mínimo nº de muestras que ya NO es plausible que se 
#      deba únicamente al ruido de fondo, con un nivel de significación alpha.
#      Esto equivale, por definición de cuantil, a:
#           m = qbinom(1 - alpha, size = n_samples, prob = p) + 1
calc_m_nbinom <- function(df, cols, n_reads = 1, n_samples, alpha = 0.05) {
    
    x <- unlist(df[, cols], use.names = FALSE)
    x <- x[!is.na(x)]
    
    mu0  <- mean(x)
    var0 <- var(x)
    
    if (is.finite(var0) && var0 > mu0 && mu0 > 0) {
        size0 <- mu0^2 / (var0 - mu0)
    } else {
        # Sin sobredispersión detectable (o mu0 = 0) -> límite Poisson
        size0 <- 1e6
    }
    
    # p = P(una muestra alcanza >= n_reads lecturas bajo el modelo de ruido)
    if (mu0 <= 0) {
        p <- 0
    } else {
        p <- 1 - pnbinom(n_reads - 1, mu = mu0, size = size0)
    }
    p <- min(max(p, 1e-10), 1 - 1e-10)  # evitar degeneraciones en qbinom
    
    # m = menor entero tal que P(Binom(n_samples, p) >= m) <= alpha
    m <- qbinom(1 - alpha, size = n_samples, prob = p) + 1
    
    # Acotar m al rango [1, n_samples]
    m <- max(1, min(m, n_samples))
    
    return(m)
}

# BLOQUE 2: expresión en al menos 1 (filtro=1) o "m" (filtro=2) muestras 
# por condición. Para filtro=2, "m" se calcula automáticamente a partir 
# del nº de muestras de cada condición mediante calc_m_nbinom().
filt_expr_m_muestras_condicion <- function(df, seq, filtro, n_reads = 1, alpha = 0.05) {
    cols <- cols_sel(seq = seq)
    
    if (filtro == 1) {
        # Filtro: al menos 1 muestra por condición
        m1 <- 1; m2 <- 1
    } else {
        # Filtro: al menos "m" muestras por condición, m ajustado según 
        # el nº de muestras de cada condición y el ruido de fondo (NB)
        m1 <- calc_m_nbinom(df = df, cols = cols$cols1, n_reads = n_reads,
                            n_samples = length(cols$cols1), alpha = alpha)
        m2 <- calc_m_nbinom(df = df, cols = cols$cols2, n_reads = n_reads,
                            n_samples = length(cols$cols2), alpha = alpha)
    }
    
    df[which(
        rowSums(df[, cols$cols1] >= n_reads) >= m1 |
            rowSums(df[, cols$cols2] >= n_reads) >= m2
    ), ]
}

# FILTRO 3: expresión en TODAS las condiciones (independiente de los 
# bloques 1 y 2 anteriores)
filt_expr_todas_condiciones <- function(df, seq) {
    cols <- cols_sel(seq = seq)
    
    if (seq == "masseq"){
        df[which(
            rowSums(df[, cols$cols1] > 0) > 0 & 
                rowSums(df[, cols$cols2] > 0) > 0 &
                rowSums(df[, cols$cols3] > 0) > 0) , ] 
    }
    else {
        df[which(
            rowSums(df[, cols$cols1] > 0) > 0 & 
                rowSums(df[, cols$cols2] > 0) > 0 &
                rowSums(df[, cols$cols3] > 0) > 0 &
                rowSums(df[, cols$cols4] > 0) > 0 ) , ] 
    }
}

# WORKFLOW PASO 3. CONVERTIR A TPM.
TPM <- function(df, cols){
    cols <- unlist(cols, use.names = FALSE)
    
    # division por filas. NO SE USA EN LONG READS   
    # df <- df_normal[,cols] / df_normal[,"length"] Esto ya no hace falta pq en logn read no hay sesgo por longitud de lectura, PERO ONT SÍ TIENE
    # suma de columas
    a <- colSums(df[, cols], na.rm = TRUE) / 1000000
    a[a == 0] <- 1  # Evitar división por cero
    
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
    
    if (seq != "masseq"){ # Añadir condicion B80
        
        df_f["B80"] <- rowMeans(df_modo[, cols$cols4])
        df_f["B80_ex"] <- df_f[, "B100"] * 0.8 + df_f[, "K100"] * 0.2
        
    }
    
    return(df_f)
    
}

# Añadir columna transcrito asociado y categoria estructural. Ya no hace 
# falta releer el fichero QC aquí: ambas columnas ya están presentes en 
# df_modo desde el WORKFLOW PASO 0 (anadir_categoria_estructural)
anotar <- function(df_f, df_modo) {
    df_f["associated_transcript"] <- df_modo[, "associated_transcript"]
    df_f["structural_category"]   <- df_modo[, "structural_category"]
    
    return(df_f)
}



# WORKFLOW PASO 5. FILTRADO DE ISOFORMAS
isoforms_filt <- function(df, seq) {
    
    if (seq == "masseq"){ cols_exp <- c("K100", "B100", "B20", "B20_ex") }
    
    else { cols_exp <- c("K100", "B100", "B20", "B80", "B20_ex", "B80_ex") }
    
    # Filtrar FSM y no noveles
    
    df_filtered <- df[which(
        df[, "associated_transcript"] != "novel" &
            df[, "structural_category"] == "full-splice_match" ) , ]
    
    
    
    return(df_filtered)
}
# WORKFLOW PASO 6. Log10 + pseudocount
log10_pseudocount <- function(df, seq){
    
    if (seq == "masseq"){ cols_exp <- c("K100", "B100", "B20", "B20_ex") }
    
    else { cols_exp <- c("K100", "B100", "B20", "B80", "B20_ex", "B80_ex") }
    
    # Transformación logarítmica base 10 con pseudocont 0.01
    df[, cols_exp] <- log10(df[cols_exp] + 0.01)
    
    return(df)
    
    
}



procesar_datos <- function(modo, seq, plataforma, ruta, filtro, ext, dir_tablas) {
    combi <- paste(seq, plataforma, sep = "_")
    
    extension <- "_classification.txt"
    # Selección de archivo de clasificación
    file_data <- switch(modo,
                        "raw" = file.path(ruta, paste0("counts_", combi, ".tsv")),
                        "qc" = file.path(ruta, paste0("default_", combi, extension)),
                        "fl" = file.path(ruta, paste0("rules_default_", combi, "_RulesFilter", extension)),
                        "rq" = file.path(ruta, paste0("rq_", combi, "_rescued", extension))
    )
    # Seleccion de columna que será el identificador
    if (modo == "raw") {name <- "superPBID"} else {name <- "isoform"}
    
    
    # Lectura del archivo (qc, fl o rq)
    df_modo <- read.table(file_data, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
    #df_modo <- read.table(file_data, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
    n0 <- nrow(df_modo)
    
    # WORKFLOW PASO 0. AÑADIR CATEGORÍA ESTRUCTURAL (para poder 
    # monitorizarla en cada paso posterior)
    df_modo <- anadir_categoria_estructural(df_modo = df_modo, modo = modo, combi = combi,
                                            extension = extension, ruta = ruta, name = name)
    cat_inicial <- cat_estructural_resumen(df_modo)
    
    # WORKFLOW PASO 1. LIMPIEZA DE NA
    cols <- cols_sel(seq = seq)
    df_sin_na <- NA_to_0(df = df_modo, cols = cols)
    
    # WORKFLOW PASO 2. FILTRO DE EXPRESION MINIMA (dos bloques independientes)
    n_bloque1 <- NA; n_bloque2 <- NA; n_todas_cond <- NA
    cat_bloque1 <- "NA"
    
    df_exprs_min <- df_sin_na
    
    if (!(filtro %in% c(0, 1, 2, 3))) { filtro <- 2 }
    
    if (filtro %in% c(1, 2)) {
        # BLOQUE 1: expresión en alguna de las dos condiciones predictoras
        df_exprs_min <- filt_expr_alguna_condicion(df = df_exprs_min, seq = seq)
        n_bloque1   <- nrow(df_exprs_min)
        cat_bloque1 <- cat_estructural_resumen(df_exprs_min)
        
        # BLOQUE 2: expresión en al menos 1 (filtro=1) o m (filtro=2) muestras por condición
        df_exprs_min <- filt_expr_m_muestras_condicion(df = df_exprs_min, seq = seq, filtro = filtro)
        n_bloque2   <- nrow(df_exprs_min)
        
    } else if (filtro == 3) {
        # Expresión en TODAS las condiciones
        df_exprs_min <- filt_expr_todas_condiciones(df = df_exprs_min, seq = seq)
        n_todas_cond <- nrow(df_exprs_min)
    }
    # filtro == 0: sin filtrado de expresión mínima (df_exprs_min = df_sin_na)
    
    n1 <- nrow(df_exprs_min)
    cat_expr_final <- cat_estructural_resumen(df_exprs_min)
    
    # WORKFLOW PASO 3. CONVERTIR A TPM.
    if (ext == "TPM"){
        df_exprs_min <- TPM(df = df_exprs_min, cols = cols)
    }
    # WORKFLOW PASO 4. MEDIA RÉPLICAS BIOLÓGICAS + MEDIAS TEÓRICAS
    
    df_bio_repl_mean <- bio_replicate_mean(df_modo = df_exprs_min, name = name, cols = cols, seq = seq) 
    
    # Anotar
    df_anotado <- anotar(df_f = df_bio_repl_mean, df_modo = df_exprs_min)
    
    # WORKFLOW PASO 5. FILTRADO DE ISOFORMAS
    
    df_isoformas_filtradas <- isoforms_filt(df = df_anotado, seq = seq)
    n2 <- nrow(df_isoformas_filtradas)
    cat_final <- cat_estructural_resumen(df_isoformas_filtradas)
    
    # Guardar la tabla de datos ya filtrada en una carpeta nueva
    guardar_tabla_filtrada(df = df_isoformas_filtradas, dir_tablas = dir_tablas,
                           combi_name = paste(modo, seq, plataforma, sep = "_"))
    
    # WORKFLOW PASO 6. Log10 + pseudocount
    df_f <- log10_pseudocount(df = df_isoformas_filtradas, seq = seq)
    
    return(list(
        df = df_f,
        n0 = n0, n_bloque1 = n_bloque1, n_bloque2 = n_bloque2, n_todas_cond = n_todas_cond,
        n1 = n1, n2 = n2,
        cat_inicial = cat_inicial, cat_bloque1 = cat_bloque1,
        cat_expr_final = cat_expr_final, cat_final = cat_final
    ))
}


generar_y_guardar_plot <- function(df_data, var_x, var_y, titulo, filename, target_dir) {
    
    # modelo lineal
    fit <- lm(as.formula(paste(var_y, "~", var_x)), data = df_data)
    r_squared <- summary(fit)$r.squared
    
    # Ruta png
    filepath <- file.path(target_dir, paste0(filename, ".png"))
    
    # ver si grafico ya existe pq es lo q más tarda en ejecutarse
    if (!file.exists(filepath)) {
        p <- ggplot(df_data, aes_string(x = var_x, y = var_y)) +
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
                data = data.frame(npcx = 0.05, npcy = 0.95, label = paste0("R² = ", round(r_squared, 4))),
                aes(npcx = npcx, npcy = npcy, label = label),
                inherit.aes = FALSE,
                size = 4.5
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
    
    return(r_squared)
}

ruta   <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/data/data_expression_matrix"
outdir <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/output/expression_matrix/graficos"
outdir_tablas <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/output/expression_matrix/tablas_filtradas"

outdir2     <- file.path(outdir, paste0("_", ext, "_filtro_", fil))
dir_destino <- file.path(outdir2, s, p, m)
dir.create(dir_destino, recursive = TRUE, showWarnings = FALSE)

outdir2_tablas     <- file.path(outdir_tablas, paste0("_", ext, "_filtro_", fil))
dir_destino_tablas <- file.path(outdir2_tablas, s, p, m)

res <- tryCatch(procesar_datos(m, s, p, ruta, fil, ext, dir_destino_tablas), error = function(e) NULL)

df_proc      <- res$df
n0           <- res$n0
n_bloque1    <- res$n_bloque1
n_bloque2    <- res$n_bloque2
n_todas_cond <- res$n_todas_cond
n1           <- res$n1
n2           <- res$n2
cat_inicial     <- res$cat_inicial
cat_bloque1     <- res$cat_bloque1
cat_expr_final  <- res$cat_expr_final
cat_final       <- res$cat_final

if (!is.null(df_proc) && nrow(df_proc) > 0) {
    combi_name <- paste(m, s, p, sep = "_")
    
    # Plot B20
    r2_b20 <- generar_y_guardar_plot(df_data = df_proc, 
                                     var_x = "B20_ex", 
                                     var_y = "B20", 
                                     titulo = toupper(paste(ext, m, s, p, "- B20K80")), 
                                     filename = paste0(combi_name, "_B20K80"), 
                                     target_dir = dir_destino
    )
    
    cat(paste(ifelse(ext == "class", "counts", "TPM"), fil, s, p, m, "B20K80", round(r2_b20, 4),
              n0, n_bloque1, n_bloque2, n_todas_cond, n1, n2,
              cat_inicial, cat_bloque1, cat_expr_final, cat_final, sep = ","), "\n")
    
    # Plot B80
    if (s != "masseq") {
        r2_b80 <- generar_y_guardar_plot(df_data = df_proc, 
                                         var_x = "B80_ex", 
                                         var_y = "B80", 
                                         titulo = toupper(paste(ext, m, s, p, "- B80K20")), 
                                         filename = paste0(combi_name, "_B80K20"), 
                                         target_dir = dir_destino
        )
        
        cat(paste(ifelse(ext == "class", "counts", "TPM"), fil, s, p, m, "B80K20", round(r2_b80, 4),
                  n0, n_bloque1, n_bloque2, n_todas_cond, n1, n2,
                  cat_inicial, cat_bloque1, cat_expr_final, cat_final, sep = ","), "\n")
    }
}