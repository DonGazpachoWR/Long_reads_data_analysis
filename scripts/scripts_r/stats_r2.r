library(ggplot2)
library(readr)
datos <- read_csv("~/Documentos/Conesa_Lab/VSCODE/practica/graficos_filtro2/resumen_r_cuadrado.csv")
out_dir <- "~/Documentos/Conesa_Lab/VSCODE/practica/graficos_filtro2/"
# factorizar para poner en orden los modos
datos$tipo_de_modos <- factor(datos$tipo_de_modos, 
                              levels = c("raw", "qc", "fl", "rq"))
datos$tipo_de_seqs <- factor(datos$tipo_de_seqs, 
                              levels = c("masseq", "isoseq", "ont"))


p <- ggplot(datos, aes(x = tipo_de_modos, 
                  y = r_cuadrado, 
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
                title = expression(paste("Evolución de ", R^2, " según Modo y Plataforma")),
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
    filename   = paste0("comparacion.png"), 
    plot       = p, 
    path       = out_dir, 
    create.dir = TRUE, 
    device     = "png"
)
