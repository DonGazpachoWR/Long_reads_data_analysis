# Shared colours and theme of the release-6.1 figures (sourced by the other scripts)
colores_plat <- c("bambu" = "orange", "flair" = "cyan", "isoquant" = "purple",
                  "isocall" = "green", "isoseq" = "firebrick")

tema <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5,
                                           margin = ggplot2::margin(b = 12)),
        strip.background = ggplot2::element_rect(fill = "grey91", color = NA),
        strip.text = ggplot2::element_text(face = "bold", size = 10, color = "black"),
        panel.grid.major = ggplot2::element_line(color = "grey87", linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "grey86", fill = NA, linewidth = 0.5),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(face = "bold", size = 9),
        legend.text = ggplot2::element_text(size = 9))

ordenar_factores <- function(datos) {
    datos$tipo_de_modos <- factor(datos$tipo_de_modos, levels = c("raw", "qc", "fl", "rq"))
    datos$tipo_de_seqs  <- factor(datos$tipo_de_seqs, levels = c("masseq", "isoseq", "ont"))
    datos
}
