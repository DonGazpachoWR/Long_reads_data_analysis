#!/usr/bin/env Rscript
# Phase-by-phase comparison, in TPM, of two runs (a against b): % R² increase,
# % RMSE reduction and % of FSM transcripts kept by a relative to b.
#
# Usage: Rscript stats_delta_fases.r <OUT_DIR> <run_a> <run_b>   (e.g. r61 kb2)

suppressPackageStartupMessages({
    library(ggplot2)
    library(readr)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; run_a <- args[2]; run_b <- args[3]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                 "tema_release61.r"))

datos <- read_csv(file.path(dir, "resultados_sqanti.csv"), show_col_types = FALSE) %>%
    ordenar_factores()

claves <- c("tipo_de_seqs", "tipo_de_plataformas", "tipo_de_modos", "comparacion")
metricas <- c("r_cuadrado", "rmse", "num_isoformas_fsm")

delta <- inner_join(
    datos %>% filter(modo_medida == "TPM", ejecucion == run_a) %>% select(all_of(c(claves, metricas))),
    datos %>% filter(modo_medida == "TPM", ejecucion == run_b) %>% select(all_of(c(claves, metricas))),
    by = claves, suffix = c("_a", "_b")) %>%
    mutate(delta_r2 = (r_cuadrado_a - r_cuadrado_b) / r_cuadrado_b * 100,
           pct_rmse_reduction = (rmse_b - rmse_a) / rmse_b * 100,
           pct_isoformas_retenidas = num_isoformas_fsm_a / num_isoformas_fsm_b * 100) %>%
    arrange(tipo_de_seqs, tipo_de_plataformas, comparacion, tipo_de_modos)

write_csv(delta, file.path(dir, paste0("delta_", run_a, "_vs_", run_b, ".csv")))

nombre <- toupper(paste("TPM", run_a, "vs", run_b))
grafico_delta <- function(metrica, titulo, eje_y) {
    ggplot(delta, aes(x = tipo_de_modos, y = .data[[metrica]],
                      color = tipo_de_plataformas, group = tipo_de_plataformas)) +
        geom_hline(yintercept = if (metrica == "pct_isoformas_retenidas") 100 else 0,
                   linetype = "dashed", color = "grey60") +
        geom_line(linewidth = 0.7, alpha = 0.9) +
        geom_point(size = 2.5) +
        facet_grid(comparacion ~ tipo_de_seqs, scales = "free_y") +
        scale_color_manual(values = colores_plat) +
        labs(title = titulo, x = "Fase SQANTI", y = eje_y, color = "Plataforma") +
        tema
}
ggsave(paste0(nombre, "_reduccion_rmse.png"),
       grafico_delta("pct_rmse_reduction", paste(nombre, "- Reducción RMSE"), "% reducción de error"),
       path = dir, width = 9, height = 6, bg = "white")
ggsave(paste0(nombre, "_incremento_r2.png"),
       grafico_delta("delta_r2", bquote(.(nombre) ~ "- Incremento de" ~ R^2), expression(Delta ~ R^2 ~ "(%)")),
       path = dir, width = 9, height = 6, bg = "white")
ggsave(paste0(nombre, "_transcritos_retenidos.png"),
       grafico_delta("pct_isoformas_retenidas", paste(nombre, "- Transcritos FSM retenidos"), "% de transcritos retenidos"),
       path = dir, width = 9, height = 6, bg = "white")
