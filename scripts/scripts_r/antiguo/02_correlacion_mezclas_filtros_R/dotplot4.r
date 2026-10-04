# Librerías
library(ggplot2)
library(ggpp)
library(ggpointdensity)
library(viridis)
library(dplyr)

# Rutas base
ruta   <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/data/data_expression_matrix"
outdir <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/output/expression_matrix/graficos"

# Función auxiliar para etiquetas en base 10 en ggplot
label_10_pow <- function(x) {
    parse(text = paste0("10^", x))
}

# Transformar a TPM
TPM <- function(df, cols){
    cols <- unlist(cols, use.names = F)
    # division por filas. NO SE USA EN LONG READS   
    # df <- df_normal[,cols] / df_normal[,"length"] Esto ya no hace falta pq en logn read no hay sesgo por longitud de lectura, PERO ONT SÍ TIENE
    # suma de columas
    a <- colSums(df)/ 1000000
    # divison por columnas
    df[,cols] <- sweep(df[, cols], 2, a, FUN = "/")
    
    return(df)
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



min_exprs_filt <- function(df, seq, filtro = 2, ext){
    
    if (!(filtro %in% c(0,1,2,3) )) { filtro <- 2}
    
    cols <- cols_sel(seq= seq)
    
    if (ext == "TPM"){
        df <- TPM(df, cols)
    }
    
    if (filtro != 0) { # Filtro = 0 no filtering
        if (filtro %in% c(1,2)){
            # Filter: At least 1 read per condition
            if (filtro == 1) { n <- 1; m <- 1 }
            
            # Filter: At least n reads per m samples per condition
            else if (filtro == 2) { n <- 1; m <- 2}
            
            # All filters imply expresion in at least one of the two predictor conditions
            df <-df[which(
                # rowsums cuenta cuantas columnas cumplen el criterio
                rowSums(df[, cols$cols1] >= n) >= m | 
                rowSums(df[, cols$cols2] >= n) >= m ) , ] 
        }
        else if (filtro == 3) { # rowsums cuenta cuantas columnas cumplen el criterio    
            if (seq == "masseq"){
                
                df <-df[which(
                    rowSums(df[, cols$cols1] > 0) > 0 & 
                    rowSums(df[, cols$cols2] > 0) > 0 &
                    rowSums(df[, cols$cols3] > 0) > 0) , ] 
            }
            else {
                # rowsums cuenta cuantas columnas cumplen el criterio
                df <-df[which(
                    rowSums(df[, cols$cols1] > 0) > 0 & 
                    rowSums(df[, cols$cols2] > 0) > 0 &
                    rowSums(df[, cols$cols3] > 0) > 0 &
                    rowSums(df[, cols$cols4] > 0) > 0 ) , ] 
                
            }
        }
    }
    
    return(df)
}


isoforms_filt <- function(df, seq) {
    
    if (seq == "masseq"){ cols_exp <- c("K100", "B100", "B20", "B20_ex") }
    
    else { cols_exp <- c("K100", "B100", "B20", "B80", "B20_ex", "B80_ex") }
    
    # Filtrar FSM y no noveles
    
    df_filtered <- df[which(
        df[, "associated_transcript"] != "novel" &
        df[, "structural_category"] == "full-splice_match" ) , ]
    
    # Transformación logarítmica base 10 con pseudocont 0.01
    df_filtered[, cols_exp] <- log10(df_filtered[cols_exp] + 0.01)
    
    return(df_filtered)
}




procesar_datos <- function(modo, seq, plataforma, ruta, filtro, ext) {
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
    
    # Limpieza reads
    df_modo <- min_exprs_filt(df_modo, seq, filtro, ext)
    
    cols <- cols_sel(seq=seq)
        
    # Aplicar la media y mezclas teoricas
        
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
    # Añadir columna transcrito asociado y cateogria estructural 
    if (modo == "raw"){
        file_qc <- file.path(ruta, paste0("default_", combi, extension))  
        df_qc <- read.table(file_qc, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
        
        # Unir data frames 
        df_f <- df_f %>%
            left_join(
                df_qc %>% select(isoform, associated_transcript, structural_category),
                by = c("tr_id" = "isoform")
            )
    }
    else {
        df_f["associated_transcript"] <- df_modo[, "associated_transcript"]
        df_f["structural_category"]   <- df_modo[, "structural_category"] }
    
    # Llamar a funcion para limpiar datos
    
    df_f <- isoforms_filt(df_f, seq)
    
    return(df_f)
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

# Ejecutar con bucles para multiples datos
extension <- c("class", "TPM")
filtro <- c(0,1,2,3)
modos       <- c("raw", "qc", "fl", "rq" )
seqs        <- c("isoseq", "masseq", "ont")
plataformas <- c("isoseq", "isocall", "bambu", "flair", "isoquant")


for (ext in extension){

    for (fil in filtro){
        lista_resumen <- list()

        outdir2  <- file.path(outdir, paste0("_", ext, "_filtro_", fil))
        if (! dir.exists(outdir2)){
            dir.create(outdir2)
        }
        for (s in seqs) {
            if (s == "ont"){ plataformas <- c("bambu", "flair", "isoquant")}
            else { plataformas <- c("isoseq", "isocall", "bambu", "flair", "isoquant") }
            for (p in plataformas) {
                for (m in modos) {
                    print(c(ext, fil, s, p, m))
                    
                    dir_destino <- file.path(outdir2, s, p, m)
                    
                    df_proc <- procesar_datos(modo = m, seq = s, plataforma = p, ruta = ruta, filtro = fil, ext = ext)
                
                    
                    if (!is.null(df_proc) && nrow(df_proc) > 0) {
                        
                        combi_name <- paste(m, s, p, sep = "_")
                        
                        # Plot 1: B20K80
                        nombre1 <- paste(combi_name, "B20K80", sep = "_")
                        titulo1 <- paste(toupper(m), toupper(s), toupper(p), "- B20K80", sep = " ")
                        
                        r2_b20 <- generar_y_guardar_plot(
                            df_data    = df_proc, 
                            var_x      = "B20_ex", 
                            var_y      = "B20", 
                            titulo     = titulo1, 
                            filename   = nombre1, 
                            target_dir = dir_destino
                        )
                        # Guardar registros en la lista
                        lista_resumen[[length(lista_resumen) + 1]] <- data.frame(
                            modo_medida = ifelse(ext == "class", "counts", "TPM"),
                            filtro = fil,
                            tipo_de_seqs        = s,
                            tipo_de_plataformas = p,
                            tipo_de_modos       = m,
                            comparacion         = "B20K80",
                            r_cuadrado    = round(r2_b20, 4)
                        )
                        
                        # Plot 2: B80K20
                        if (s != "masseq"){
                            nombre2 <- paste(combi_name, "B80K20", sep = "_")
                            titulo2 <- paste(toupper(m), toupper(s), toupper(p), "- B80K20", sep = " ")
                            
                            r2_b80 <- generar_y_guardar_plot(
                                df_data    = df_proc, 
                                var_x      = "B80_ex", 
                                var_y      = "B80", 
                                titulo     = titulo2, 
                                filename   = nombre2, 
                                target_dir = dir_destino
                            )
                            
                            lista_resumen[[length(lista_resumen) + 1]] <- data.frame(
                                modo_medida = ifelse(ext == "class", "counts", "TPM"),
                                filtro = fil,
                                tipo_de_seqs        = s,
                                tipo_de_plataformas = p,
                                tipo_de_modos       = m,
                                comparacion         = "B80K20",
                                r_cuadrado    = round(r2_b80, 4)
                            )
                        }
                    }
                    else {print("error tras procesar datos")}
                }
            }
        }
        
        df_resumen <- bind_rows(lista_resumen)
        
        ruta_csv <- file.path(outdir2, "resumen_r_cuadrado.csv")
        write.csv(df_resumen, file = ruta_csv, row.names = FALSE)
    }
}