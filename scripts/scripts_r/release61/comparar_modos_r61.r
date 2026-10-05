#!/usr/bin/env Rscript
# Visión global de los dos modos de la prueba de release 6.1 (con prevalencia,
# r61, y sin prevalencia, sinprev), a partir de las salidas de
# test_reasignacion_azar_r61.r de cada modo:
#   - rescue frente a reasignaciones al azar dentro del gen (aleatorización),
#   - rescue (rq) frente a los datos filtrados (fl) (bootstrap de genes).
#
# Salidas: comparacion_modos_TPM.csv (una fila por modo, combinación y mezcla),
# recuento_modos.csv (nº de casos significativos) y comparacion_modos_TPM.png.
#
# Uso: Rscript comparar_modos_r61.r [DIR_TEST] (por defecto output/benchmark_r61/test_reasignacion_azar)

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(ggplot2)
})

base_dir <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis"
args <- commandArgs(trailingOnly = TRUE)
dir_test <- if (length(args) >= 1) args[1] else file.path(base_dir, "output/benchmark_r61/test_reasignacion_azar")
MODOS <- c("r61" = "Con prevalencia", "sinprev" = "Sin prevalencia")
ALFA <- 0.05

dir_script <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(dir_script, "tema_release61.r"))

leer_modo <- function(modo, fichero) {
    f <- file.path(dir_test, modo, fichero)
    if (!file.exists(f)) stop("no existe ", f)
    read_csv(f, show_col_types = FALSE) %>% mutate(ejecucion = modo)
}

azar <- bind_rows(lapply(names(MODOS), leer_modo, "resumen_test_reasignacion.csv")) %>%
    filter(medida == "TPM", nulo == "gen") %>%
    select(ejecucion, combinacion, mezcla,
           r2_azar_media, r2_azar_ic95_inf, r2_azar_ic95_sup, p_r2_azar_bh = p_r2_bh,
           rmse_azar_media, rmse_azar_ic95_inf, rmse_azar_ic95_sup, p_rmse_azar_bh = p_rmse_bh)
boot <- bind_rows(lapply(names(MODOS), leer_modo, "bootstrap_fl_rq.csv")) %>%
    filter(medida == "TPM") %>%
    select(ejecucion, combinacion, mezcla, r2_fl, r2_rq, delta_r2, delta_r2_ic95_inf, delta_r2_ic95_sup,
           p_r2_fl_bh = p_r2_bh, rmse_fl, rmse_rq, delta_rmse, delta_rmse_ic95_inf, delta_rmse_ic95_sup,
           p_rmse_fl_bh = p_rmse_bh)

tabla <- boot %>%
    inner_join(azar, by = c("ejecucion", "combinacion", "mezcla")) %>%
    mutate(modo = MODOS[ejecucion],
           plataforma = sub("^[^_]*_", "", combinacion),
           r2_mejor_que_fl = p_r2_fl_bh < ALFA, rmse_mejor_que_fl = p_rmse_fl_bh < ALFA,
           r2_mejor_que_azar = p_r2_azar_bh < ALFA, rmse_mejor_que_azar = p_rmse_azar_bh < ALFA) %>%
    relocate(modo, .after = ejecucion) %>%
    arrange(combinacion, mezcla, ejecucion)
write_csv(tabla, file.path(dir_test, "comparacion_modos_TPM.csv"))

recuento <- tabla %>%
    group_by(modo) %>%
    summarise(casos = n(),
              r2_mejor_que_fl = sum(r2_mejor_que_fl), rmse_mejor_que_fl = sum(rmse_mejor_que_fl),
              r2_mejor_que_azar = sum(r2_mejor_que_azar), rmse_mejor_que_azar = sum(rmse_mejor_que_azar),
              .groups = "drop")
recuento_plat <- tabla %>%
    group_by(modo, plataforma) %>%
    summarise(casos = n(),
              r2_mejor_que_fl = sum(r2_mejor_que_fl), rmse_mejor_que_fl = sum(rmse_mejor_que_fl),
              r2_mejor_que_azar = sum(r2_mejor_que_azar), rmse_mejor_que_azar = sum(rmse_mejor_que_azar),
              .groups = "drop")
write_csv(bind_rows(recuento %>% mutate(plataforma = "todas"), recuento_plat) %>% relocate(plataforma, .after = modo),
          file.path(dir_test, "recuento_modos.csv"))

# Figura: cambio respecto de fl (0) en cada modo; punto e intervalo bootstrap de rq,
# barra fina con el intervalo del 95% de las reasignaciones al azar dentro del gen
largo <- bind_rows(
    tabla %>% transmute(modo, combinacion, mezcla, metrica = "Δ R² (mejor a la derecha)",
                        rq = delta_r2, rq_inf = delta_r2_ic95_inf, rq_sup = delta_r2_ic95_sup,
                        azar_inf = r2_azar_ic95_inf - r2_fl, azar_sup = r2_azar_ic95_sup - r2_fl),
    tabla %>% transmute(modo, combinacion, mezcla, metrica = "Δ RMSE (mejor a la izquierda)",
                        rq = delta_rmse, rq_inf = delta_rmse_ic95_inf, rq_sup = delta_rmse_ic95_sup,
                        azar_inf = rmse_azar_ic95_inf - rmse_fl, azar_sup = rmse_azar_ic95_sup - rmse_fl)) %>%
    mutate(panel = paste(combinacion, mezcla),
           seq = factor(sub("_.*", "", combinacion), levels = c("masseq", "isoseq", "ont")),
           modo = factor(modo, levels = MODOS),
           metrica = factor(metrica, levels = c("Δ R² (mejor a la derecha)", "Δ RMSE (mejor a la izquierda)")))

colores_modo <- c("Con prevalencia" = "#2a78d6", "Sin prevalencia" = "#eb6834")
esquiva <- position_dodge(width = 0.7)
p <- ggplot(largo, aes(y = panel, colour = modo, group = modo)) +
    geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.5) +
    geom_linerange(aes(xmin = azar_inf, xmax = azar_sup), linewidth = 2.6, alpha = 0.25, position = esquiva) +
    geom_errorbar(aes(xmin = rq_inf, xmax = rq_sup), orientation = "y", width = 0.3, linewidth = 0.6,
                  position = esquiva) +
    geom_point(aes(x = rq), size = 2, position = esquiva) +
    facet_grid(seq ~ metrica, scales = "free", space = "free_y") +
    scale_colour_manual(values = colores_modo) +
    labs(title = "Rescue con requantificación frente a los datos filtrados y al azar (TPM)",
         subtitle = paste("0 = datos filtrados (fl) de cada modo. Punto y barra: rescue (rq) con su intervalo bootstrap del 95%.",
                          "Banda clara: 95% de las reasignaciones al azar dentro del gen", sep = "\n"),
         x = NULL, y = NULL, colour = NULL) +
    tema
ggsave("comparacion_modos_TPM.png", p, path = dir_test, width = 12, height = 11, bg = "white")

print(as.data.frame(recuento), row.names = FALSE)
cat("\n")
print(as.data.frame(recuento_plat), row.names = FALSE)
cat("Salidas en", dir_test, "\n")
