# Librerías
library(ggplot2)
library(ggpp)
library(ggpointdensity)
library(viridis)
library(dplyr)

# Rutas base
ruta   <- "~/Documentos/Conesa_Lab/VSCODE/practica/datos/"
outdir <- "~/Documentos/Conesa_Lab/VSCODE/practica/graficos/"

# Funcion que etiqueta exponentes en base 10 para ggplot
label_10_pow <- function(x) {
        parse(text = paste0("10^", x))
}


procesar_datos <- function(modo, seq, plataforma, ruta) {
        combi <- paste(seq, plataforma, sep = "_")
        
        file_raw <- file.path(ruta, paste0("counts_", combi, ".tsv"))
        file_qc <- file.path(ruta, paste0("default", modo, "_", combi, "_classification.txt"))
        file_fl <- file.path(ruta, paste0("rules_default_", modo, "_", combi, "_RulesFilter_classification.txt"))
        file_rq  <- file.path(ruta, paste0(modo, "_", combi, "_rescued_classification.txt"))
        
        # Verificación de que los archivos existen antes de leer
        if (!file.exists(file_rq) || !file.exists(file_raw) || !file.exists(file_fl) || !file.exists(file_rq)) {
                warning(paste("Archivo(s) no encontrado(s) para:", modo, seq, plataforma))
                return(NULL)
        }
        
        # Lectura de datos
        
        df_raw <- read.table(file_raw, sep = "\t", header = TRUE)
        df_qc  <- read.table(file_qc, sep = "\t", header = TRUE)
        df_fl  <- read.table(file_fl, sep = "\t", header = TRUE)
        df_rq  <- read.table(file_rq, sep = "\t", header = TRUE)
        
        # Generar data frame con medias de expresión
        df_f_raw <- data.frame(
                tr_id = df_raw[, "superPBID"], 
                K100  = apply(df_raw[, 2:6], 1, mean),
                B100  = apply(df_raw[, 7:11], 1, mean),
                B20   = apply(df_raw[, 12:16], 1, mean),
                B80   = apply(df_raw[, 17:21], 1, mean)
        )
        
        # Mezclas teóricas
        df["B20_ex"] <- df[, "B100"] * 0.2 + df[, "K100"] * 0.8
        df["B80_ex"] <- df[, "B100"] * 0.8 + df[, "K100"] * 0.2    
        
        # Unir data frames 
        
        df <- df %>% 
                left_join(
                        df_rq %>% select("isoform", "associated_transcript", 
                                            "structural_category"),
                        by = c("tr_id" ="isoform" )
                        
                )
        # filtrar FSM, no noveles y que tengan expresión en B100 y K100
        # Si el valor esperado es 0, la mezcla no tiene nada que validar sobre la proporción $80/20$.
        df_filtered <- df[which(
                df[, "associated_transcript"] != "novel" &
                        df[, "structural_category"] == "full-splice_match") &
                        df[, "B100"] + df[, "K100"] > 0, ]
        
        # Transformación logarítmica base 10 con pseudocont 0.01
        cols_exp <- c("K100", "B100", "B20", "B80", "B20_ex", "B80_ex")
        df_filtered[, cols_exp] <- apply(df_filtered[, cols_exp], 2, function(x) log10(x + 0.01))
        
        return(df_filtered)
}

# Modelo lineal y graficar
generar_y_guardar_plot <- function(df_data, var_x, var_y, titulo, filename, outdir) {
        
        # Modelo lineal
        fit <- lm(as.formula(paste(var_y, "~", var_x)), data = df_data)
        r_squared <- summary(fit)$r.squared
        
        # Construcción del gráfico con ggplot2
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
                        data = data.frame(npcx = 0.05, npcy = 0.95, label = paste0("r = ", round(r_squared, 3))),
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
                theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(
                filename = paste0(filename, ".png"), 
                plot     = p, 
                path     = outdir, 
                create.dir = TRUE, 
                device   = "png"
        )
}

# Ejecutar con bucles para multiples datos
modos       <- c("rq", "raw")
seqs        <- c("ont", "isoseq")
plataformas <- c("bambu", "flair", "isoquant")

for (m in modos) {
        for (s in seqs) {
                for (p in plataformas) {
                        
                        # Carga y procesado
                        df_proc <- procesar_datos(modo = m, seq = s, plataforma = p, ruta = ruta)
                        
                        # Si el archivo existía y devolvió datos, generar gráficos
                        if (!is.null(df_proc) && nrow(df_proc) > 0) {
                                
                                combi_name <- paste(m, s, p, sep = "_")
                                
                                # Plot 1: B20K80 (Esperado: B20_ex vs Observado: B20)
                                nombre1 <- paste(combi_name, "B20K80", sep = "_")
                                generar_y_guardar_plot(
                                        df_data  = df_proc, 
                                        var_x    = "B20_ex", 
                                        var_y    = "B20", 
                                        titulo   = nombre1, 
                                        filename = nombre1, 
                                        outdir   = outdir
                                )
                                
                                # Plot 2: B80K20 (Esperado: B80_ex vs Observado: B80)
                                nombre2 <- paste(combi_name, "B80K20", sep = "_")
                                generar_y_guardar_plot(
                                        df_data  = df_proc, 
                                        var_x    = "B80_ex", 
                                        var_y    = "B80", 
                                        titulo   = nombre2, 
                                        filename = nombre2, 
                                        outdir   = outdir
                                )
                                
                        }
                }
        }
}