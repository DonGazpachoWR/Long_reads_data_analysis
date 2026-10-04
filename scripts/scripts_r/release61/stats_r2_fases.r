#!/usr/bin/env Rscript
# Evolution of R² and RMSE across SQANTI phases (raw, qc, fl, rq) for one run,
# in counts and in TPM (four figures).
#
# Usage: Rscript stats_r2_fases.r <OUT_DIR> <run>   (run as in resultados_sqanti.csv, e.g. r61)

suppressPackageStartupMessages({
    library(ggplot2)
    library(readr)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
run <- args[2]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                 "tema_release61.r"))

datos <- read_csv(file.path(dir, "resultados_sqanti.csv"), show_col_types = FALSE) %>%
    filter(ejecucion == run) %>%
    ordenar_factores()

grafico_fases <- function(df, metrica, titulo, eje_y) {
    ggplot(df, aes(x = tipo_de_modos, y = .data[[metrica]],
                   color = tipo_de_plataformas, group = tipo_de_plataformas)) +
        geom_line(linewidth = 0.6, alpha = 0.9) +
        geom_point(size = 2) +
        facet_grid(comparacion ~ tipo_de_seqs, scales = "free_y") +
        scale_color_manual(values = colores_plat) +
        labs(title = titulo, x = "Fase SQANTI", y = eje_y, color = "Plataforma") +
        tema
}

for (ext in unique(datos$modo_medida)) {
    datos_ext <- datos %>% filter(modo_medida == ext)
    nombre <- toupper(paste(ext, "SQANTI", run))
    p_r2 <- grafico_fases(datos_ext, "r_cuadrado",
                          bquote(.(nombre) ~ "- Evolución de" ~ R^2), expression(R^2)) +
        expand_limits(y = c(0.35, 1.0))
    ggsave(paste(nombre, "R2.png"), p_r2, path = dir, width = 9, height = 6, bg = "white")
    p_rmse <- grafico_fases(datos_ext, "rmse", paste(nombre, "- Evolución del RMSE"), "RMSE")
    ggsave(paste(nombre, "RMSE.png"), p_rmse, path = dir, width = 9, height = 6, bg = "white")
}
