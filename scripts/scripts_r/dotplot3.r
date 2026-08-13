# Librerías
library(ggplot2)
library(ggpp)
library(ggpointdensity)
library(viridis)
library(dplyr)

# Rutas base
ruta   <- "~/Documentos/Conesa_Lab/VSCODE/practica/datos/"
outdir <- "~/Documentos/Conesa_Lab/VSCODE/practica/graficos/"

# Función auxiliar para etiquetas en base 10 en ggplot
label_10_pow <- function(x) {
        parse(text = paste0("10^", x))
}

filtrar_y_transformar <- function(df, seq, filtro=1) {
        
        if (seq == "masseq"){
                df_filtered <- df[which(
                        df[, "associated_transcript"] != "novel" &
                                df[, "structural_category"] == "full-splice_match" &
                                df[, "B100"]  > 0 & 
                                df[, "K100"] > 0 &
                                df[, "B20"]  > 0) , ]
                
                # Transformación logarítmica base 10 con pseudocont 0.01
                cols_exp <- c("K100", "B100", "B20", "B20_ex")
                df_filtered[, cols_exp] <- apply(df_filtered[, cols_exp], 2, function(x) log10(x + 0.01))        
        }
        else {
                # Filtrar FSM, no noveles y 
                df_filtered <- df[which(
                        df[, "associated_transcript"] != "novel" &
                                df[, "structural_category"] == "full-splice_match" ) , ]
                # filtro expresion mínima
                # Filtro = 0 no hace nada, coge todo
                # Filtro = 1 expresion en al menos una de las dos condiciones predictoras y en todas las observadas
                # Filtro = 2 expresión simultanea en ambas condiciones predictoras y observadas
                # Filtro = 3 expresion simultanea en condiciones predictoras y observadas
                # Debe haber al menos 1 lectura por
                if (filtro == 1){
                df_filtered <-df_filtered[which(
                                df_filtered[, "B100"]  > 0 & 
                                df_filtered[, "K100"] > 0 &
                                df_filtered[, "B20"]  > 0 &
                                df_filtered[, "B80"]  > 0 ) , ] }
                else if (filtro == 2) {
                        df_filtered <-df_filtered[which(
                                df_filtered[, "B100"]  > 0 & 
                                df_filtered[, "K100"] > 0 &
                                df_filtered[, "B20"]  > 0 &
                                df_filtered[, "B80"]  > 0 ) , ] 
                }
                else if (filtro == 3){
                        df_filtered <-df_filtered[which(
                                df_filtered[, "B100"]  > 0 & 
                                df_filtered[, "K100"] > 0 &
                                df_filtered[, "B20"]  > 0 &
                                df_filtered[, "B80"]  > 0 ) , ] 
                        
                }
                
                # Expression:  At least 1 read per condition
                n <- 1 
                # Expression:  At least n reads per n samples per condition
                n <-  ncol(df_filtered[, "B100"])
                # Expression in at least 1 predictor condition
                
                df_filtered <-df_filtered[which(
                        df_filtered[, "B100"]  >= n || 
                                df_filtered[, "K100"] >= n ) , ] 
                # Expression in all predictors and observed conditions
                df_filtered <-df_filtered[which(
                        df_filtered[, "B100"]  >= n & 
                                df_filtered[, "K100"] >= n &
                                df_filtered[, "B20"]  >= n &
                                df_filtered[, "B80"]  >= n ) , ] 
                
                ############################################################
                # Expression:  At least 1 read per sample per condition
                
                # Expression in at least 1 predictor condition
                
                df_filtered <-df_filtered[which(
                        
                        apply(df_filtered[, "B100"]  > 0, 1, all) || 
                        apply(df_filtered[, "K100"]  > 0, 1, all) ) , ] 
                # Expression in all predictors and observed conditions
                df_filtered <-df_filtered[which(
                        apply(df_filtered[, "B100"]  > 0, 1, all) & 
                        apply(df_filtered[, "K100"]  > 0, 1, all) &
                        apply(df_filtered[, "B20"]  > 0, 1, all) &
                        apply(df_filtered[, "B80"]  > 0, 1, all) ) , ] 
                
        
                
        # Transformación logarítmica base 10 con pseudocont 0.01
        cols_exp <- c("K100", "B100", "B20", "B80", "B20_ex", "B80_ex")
        df_filtered[, cols_exp] <- apply(df_filtered[, cols_exp], 2, function(x) log10(x + 0.01))
        }
        return(df_filtered)
}

procesar_datos_raw <- function(seq, plataforma, ruta) {
        combi <- paste(seq, plataforma, sep = "_")
        
        file_data  <- file.path(ruta, paste0("counts_", combi, ".tsv"))
        file_qc <- file.path(ruta, paste0("default_", combi, "_classification.txt"))
        
        # Lectura de archivos
        df_raw   <- read.table(file_data, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        df_qc <- read.table(file_qc, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
        
        # Generar data frame con medias de expresión
        if (seq == "masseq"){
                df_f <- data.frame(
                        tr_id = df_raw[, "superPBID"],
                        K100  = apply(df_raw[, c("K31", "K32", "K33")],   1, mean),
                        B100  = apply(df_raw[, c("B31", "B32", "B33")],  1, mean),
                        B20   = apply(df_raw[, c("B20K80_1", "B20K80_2","B20K80_3")], 1, mean)
                )
                # Mezclas teóricas
                df_f["B20_ex"] <- df_f[, "B100"] * 0.2 + df_f[, "K100"] * 0.8
        }
        else {
                df_f <- data.frame(
                        tr_id = df_raw[, "superPBID"],
                        K100  = apply(df_raw[, c("K31", "K32", "K33", "K34", "K35")],   1, mean),
                        B100  = apply(df_raw[, c("B31", "B32", "B33", "B34", "B35")],  1, mean),
                        B20   = apply(df_raw[, c("B20K80_1", "B20K80_2","B20K80_3", "B20K80_4", "B20K80_5")], 1, mean),
                        B80   = apply(df_raw[, c("B80K20_1", "B80K20_2","B80K20_3", "B80K20_4", "B80K20_5")], 1, mean)
                )
                # Mezclas teóricas
                df_f["B20_ex"] <- df_f[, "B100"] * 0.2 + df_f[, "K100"] * 0.8
                df_f["B80_ex"] <- df_f[, "B100"] * 0.8 + df_f[, "K100"] * 0.2
        }

        
        # Unir data frames 
 
        df_f <- df_f %>%
                left_join(
                        df_qc %>% select(isoform, associated_transcript, structural_category),
                        by = c("tr_id" = "isoform")
                )
        
        # Llamar a fucnion para limpiar datos
        return(filtrar_y_transformar(df_f, seq))
}


procesar_datos_clasificados <- function(modo, seq, plataforma, ruta) {
        combi <- paste(seq, plataforma, sep = "_")
        
        # Selección de archivo de clasificación
        file_data <- switch(modo,
                            "qc" = file.path(ruta, paste0("default_", combi, "_classification.txt")),
                            "fl" = file.path(ruta, paste0("rules_default_", combi, "_RulesFilter_classification.txt")),
                            "rq" = file.path(ruta, paste0("rq_", combi, "_rescued_classification.txt"))
        )

        
        # Lectura del archivo (qc, fl o rq)
        df_modo <- read.table(file_data, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
        
        if (seq == "masseq"){
        df_f <- data.frame(
                tr_id = df_modo[, "isoform"],
                K100  = apply(df_modo[, c("FL.K31", "FL.K32", "FL.K33")],   1, mean),
                B100  = apply(df_modo[, c("FL.B31", "FL.B32", "FL.B33")],  1, mean),
                B20   = apply(df_modo[, c("FL.B20K80_1", "FL.B20K80_2","FL.B20K80_3")], 1, mean)
        ) 
        # Mezclas teóricas
        df_f["B20_ex"] <- df_f[, "B100"] * 0.2 + df_f[, "K100"] * 0.8
        } else {
        # Medias de expresión (a partir de col 55)
        df_f <- data.frame(
                tr_id = df_modo[, "isoform"],
                K100  = apply(df_modo[, c("FL.K31", "FL.K32", "FL.K33", "FL.K34", "FL.K35")],   1, mean),
                B100  = apply(df_modo[, c("FL.B31", "FL.B32", "FL.B33", "FL.B34", "FL.B35")],  1, mean),
                B20   = apply(df_modo[, c("FL.B20K80_1", "FL.B20K80_2","FL.B20K80_3", "FL.B20K80_4", "FL.B20K80_5")], 1, mean),
                B80   = apply(df_modo[, c("FL.B80K20_1", "FL.B80K20_2","FL.B80K20_3", "FL.B80K20_4", "FL.B80K20_5")], 1, mean)
        )
        # Mezclas teóricas
        df_f["B20_ex"] <- df_f[, "B100"] * 0.2 + df_f[, "K100"] * 0.8
        df_f["B80_ex"] <- df_f[, "B100"] * 0.8 + df_f[, "K100"] * 0.2
        }

        
        # Añadir columna transcrito asociado y cateogria estructural
        df_f["associated_transcript"] <- df_modo[, "associated_transcript"]
        df_f["structural_category"]   <- df_modo[, "structural_category"]
        
        # Llamar a fucnion para limpiar datos
        return(filtrar_y_transformar(df_f, seq))
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

filtro <- 1
modos       <- c("raw", "qc", "fl", "rq" )
seqs        <- c("isoseq", "masseq", "ont")
plataformas <- c("isoseq", "isocall", "bambu", "flair", "isoquant")

lista_resumen <- list()

for (s in seqs) {
        if (s == "ont"){ plataformas <- c("bambu", "flair", "isoquant")}
        for (p in plataformas) {
                for (m in modos) {
                        print(c(s, p, m))
                        
                        dir_destino <- file.path(outdir, s, p, m)
                        
                        if (m == "raw") {
                                df_proc <- procesar_datos_raw(seq = s, plataforma = p, ruta = ruta)
                        } else {
                                df_proc <- procesar_datos_clasificados(modo = m, seq = s, plataforma = p, ruta = ruta)
                        }
                        
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
                                        tipo_de_seqs        = s,
                                        tipo_de_plataformas = p,
                                        tipo_de_modos       = m,
                                        comparacion         = "B80K20",
                                        r_cuadrado    = round(r2_b80, 4)
                                )
                                }
                        }
                }
        }
}

df_resumen <- bind_rows(lista_resumen)

ruta_csv <- file.path(outdir, "resumen_r_cuadrado.csv")
write.csv(df_resumen, file = ruta_csv, row.names = FALSE)