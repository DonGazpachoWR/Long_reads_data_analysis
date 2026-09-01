library(ggplot2)
library(readr)
library(dplyr)
dir <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis/output"
datos <- read_csv(paste0(dir,"/resultados_globales.csv"))

etiquetas <- c("0" = "Sin filtro",
               "1" = "Filtro suave",
               "2" = "Filtro intermedio",
               "3" = "Filtro fuerte")


# factorizar para poner en orden los modos
datos$tipo_de_modos <- factor(datos$tipo_de_modos, 
                              levels = c("raw", "qc", "fl", "rq"))
datos$tipo_de_seqs <- factor(datos$tipo_de_seqs, 
                             levels = c("masseq", "isoseq", "ont"))
datos$modo_medida <- factor(datos$modo_medida, 
                            levels = c("counts", "TPM"))

# Extraccion de datos de los dos grupos a comparar
item_a <- list(ext = "TPM", filtro = 2)
item_b <- list(ext = "TPM", filtro = 0)
datos_a <- datos[which(datos$modo_medida == item_a$ext & datos$filtro == item_a$filtro ),]
datos_b <- datos[which(datos$modo_medida == item_b$ext & datos$filtro == item_b$filtro ),]

# Diferencia
claves <- c("tipo_de_seqs", "tipo_de_plataformas", "tipo_de_modos", "comparacion")

# Cruce y cálculo exacto de delta
datos_delta <- inner_join(
    datos_a[, c(claves, "r_cuadrado")],
    datos_b[, c(claves, "r_cuadrado")],
    by = claves,
    suffix = c("_a", "_b")
) %>%
    mutate(delta = r_cuadrado_a - r_cuadrado_b)


nombre_combinado <- toupper(paste(item_a$ext, etiquetas[as.character(item_a$filtro)], "menos", 
                                          item_b$ext, etiquetas[as.character(item_b$filtro)]))
        
p <- ggplot(datos_delta, aes(   x = tipo_de_modos, 
                                y = delta, 
                                color = tipo_de_plataformas, 
                                group = tipo_de_plataformas)) +
# Formato de linea
geom_line(linewidth = 0.6, alpha = 0.9) +
geom_point(size = 2) +
                
# Hacer 4 gráficos
facet_grid(comparacion ~ tipo_de_seqs, scales = "free_y") +

# Eje Y que ponga el
expand_limits(y = 1.0) +

# Eje Y que ponga el
expand_limits(y = 0.35) +

# Colores lineas
scale_color_manual(values = c("bambu" = "orange", 
                              "flair" = "cyan", 
                              "isoquant" = "purple",
                              "isocall" = "green",
                              "isoseq" = "firebrick")) +
                
# Tema del gráfico
theme_minimal(base_size = 11) +
labs(
        title = bquote(.(nombre_combinado) ~ " - Incremento de" ~ R^2),
        x = "Tipo de Modo",
        y = expression(R^2),
        color = "Plataforma"
) +
                
theme(
# posicion del titulo
plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 12)),

# Paneles
strip.background = element_rect(fill = "grey91", color = NA),
strip.text = element_text(face = "bold", size = 10, color = "black"),

# Fondo gráfico
panel.grid.major = element_line(color = "grey87", linewidth = 0.4),
panel.grid.minor = element_blank(),
panel.border = element_rect(color = "grey86", fill = NA, linewidth = 0.5),

# Leyenda 
legend.position = "bottom",
legend.title = element_text(face = "bold", size = 9),
legend.text = element_text(size = 9)
)

ggsave(
    filename   = paste(nombre_combinado, "comparacion.png"), 
    plot       = p, 
    path       = dir, 
    create.dir = TRUE, 
    device     = "png"
)
