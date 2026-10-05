#!/usr/bin/env Rscript
# ¿Mejora la reasignación de counts del rescue R² y RMSE de las mezclas de tejidos?
# Dos contrastes por combinación y mezcla, sobre los datos de una ejecución
# (r61 con prevalencia o sin prevalencia):
#
# 1. Rescue frente a reasignaciones al azar (test de aleatorización, Monte Carlo).
#    Se reproduce la requantificación de SQANTI3 rescue (requant_r61.r) con las
#    asignaciones reales y, B veces, con asignaciones al azar. Todo lo demás queda
#    fijo: las isoformas de la clasificación rescatada, los artefactos que se
#    reasignan, cuántos targets tiene cada uno, sus counts y los que van al
#    residual del gen. Solo cambia a qué isoforma van los counts.
#    Modelos nulos:
#      gen:    cada target se sustituye por una isoforma al azar del mismo gen que
#              el target real (contrasta la elección de la isoforma dentro del gen).
#      global: cada target se sustituye por una isoforma al azar de todo el
#              transcriptoma rescatado.
#    Un artefacto con k targets recibe k isoformas distintas al azar.
#    p-valor unilateral: (1 + nº de reasignaciones al azar con R² >= observado)/(B + 1),
#    y con RMSE <= observado para RMSE.
#
# 2. Rescue (rq) frente a los datos filtrados (fl), por bootstrap de genes.
#    Se remuestrean genes con reemplazo y se recalculan R² y RMSE de fl y de rq
#    sobre los mismos genes (las isoformas de un gen entran o salen juntas, porque
#    rescue reparte counts dentro del gen). Intervalo del 95% de la diferencia
#    rq - fl y p-valor unilateral de mejora: (1 + nº de réplicas sin mejora)/(B + 1).
#
# Cada fase se evalúa como en correlacion_mezclas_individual.r (FSM no novel,
# media de réplicas, log10 + 0.01, R² del ajuste lineal y RMSE frente a la mezcla
# esperada). Ajuste de Benjamini-Hochberg entre combinaciones y mezclas.
#
# Uso: Rscript test_reasignacion_azar_r61.r [DATOS] [SALIDA] [ETIQUETA] [B] [NUCLEOS] [COMBIS]
#   ETIQUETA: nombre de la ejecución (r61, sinprev...), se añade a las tablas.
#   COMBIS: lista separada por comas (p. ej. isoseq_bambu,ont_bambu); por defecto, las 13.
#   Si SALIDA ya tiene resultados, las demás combinaciones se conservan y solo se
#   sustituyen las de COMBIS (el resumen y las figuras incluyen todas).

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(ggplot2)
    library(parallel)
})

base_dir <- "/home/adrian/Documentos/Conesa_Lab/VSCODE/Long_reads_data_analysis"
args     <- commandArgs(trailingOnly = TRUE)
datos    <- if (length(args) >= 1) args[1] else file.path(base_dir, "data/data_expression_matrix_r61")
salida   <- if (length(args) >= 2) args[2] else file.path(base_dir, "output/benchmark_r61/test_reasignacion_azar")
etiqueta <- if (length(args) >= 3) args[3] else "r61"
B        <- if (length(args) >= 4) as.integer(args[4]) else 1000L
nucleos  <- if (length(args) >= 5) as.integer(args[5]) else 10L
combis   <- c(paste0("isoseq_", c("isocall", "isoseq", "bambu", "flair", "isoquant")),
              paste0("masseq_", c("isocall", "isoseq", "bambu", "flair", "isoquant")),
              paste0("ont_", c("bambu", "flair", "isoquant")))
if (length(args) >= 6) combis <- strsplit(args[6], ",")[[1]]
dir.create(salida, recursive = TRUE, showWarnings = FALSE)

dir_script <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(dir_script, "tema_release61.r"))
source(file.path(dir_script, "requant_r61.r"))

NULOS <- c("gen", "global")
COLS_ANOT <- c("isoform", "associated_gene", "associated_transcript", "structural_category")

# ------------------------------------------------------------------------------
# Datos de una combinación
# ------------------------------------------------------------------------------
preparar <- function(combi) {
    seq <- sub("_.*", "", combi)
    cols <- columnas(seq)
    muestras <- unlist(cols, use.names = FALSE)

    m_counts <- cargar_counts(datos, combi, muestras)
    fl <- leer_tsv(file.path(datos, paste0("rules_default_", combi, "_RulesFilter_classification.txt")),
               col_select = all_of(c(COLS_ANOT, "filter_result", muestras)))
    rq <- leer_tsv(file.path(datos, paste0("rq_", combi, "_rescued_classification.txt")),
               col_select = all_of(c(COLS_ANOT, muestras)))
    tabla <- leer_tsv(file.path(datos, paste0("rq_", combi, "_rescue_table.tsv")))
    # Filas que mueve la requantificación (rows_for_requant)
    if ("evidence_check" %in% names(tabla)) tabla <- tabla[tabla$evidence_check != "failed", ]

    d <- preparar_requant(m_counts, fl, tabla, rq$isoform)
    d$combi <- combi
    d$cols <- cols
    d$muestras <- muestras
    d$rq <- rq
    d$fl <- fl[fl$filter_result == "Isoform", ]
    evaluable <- function(df) df$structural_category == "full-splice_match" & df$associated_transcript != "novel"
    d$evaluadas_rq <- evaluable(rq)
    d$evaluadas_fl <- evaluable(d$fl)
    d
}

# ------------------------------------------------------------------------------
# 1. Reasignación al azar: para cada fila, una isoforma del mismo grupo (gen del
# target real, o todo el transcriptoma) distinta de las demás del mismo
# artefacto. Se baraja cada grupo y se toma una ventana circular de k
# posiciones desde un inicio al azar: un subconjunto de k isoformas uniforme.
# ------------------------------------------------------------------------------
preparar_nulo <- function(d, nulo) {
    grupo_iso <- if (nulo == "gen") as.integer(factor(d$rq$associated_gene)) else rep(1L, nrow(d$rq))
    grupo_fila <- grupo_iso[d$target]
    bloque <- as.integer(factor(paste(d$art, grupo_fila)))
    list(grupo_iso = grupo_iso,
         tam = tabulate(grupo_iso),
         inicio = c(0L, cumsum(tabulate(grupo_iso)))[seq_len(max(grupo_iso))],
         grupo_fila = grupo_fila,
         bloque = bloque,
         grupo_bloque = grupo_fila[!duplicated(bloque)][order(unique(bloque))],
         rango = ave(seq_along(bloque), bloque, FUN = seq_along))
}

target_azar <- function(n) {
    perm <- order(n$grupo_iso, runif(length(n$grupo_iso)))
    u <- floor(runif(length(n$grupo_bloque)) * n$tam[n$grupo_bloque])
    g <- n$grupo_fila
    perm[n$inicio[g] + (u[n$bloque] + n$rango - 1L) %% n$tam[g] + 1L]
}

# ------------------------------------------------------------------------------
# 2. Bootstrap de genes de rq frente a fl
# ------------------------------------------------------------------------------
r2_rmse_pesos <- function(e, o, w) {
    sw <- sum(w)
    me <- sum(w * e) / sw
    mo <- sum(w * o) / sw
    cov_eo <- sum(w * (e - me) * (o - mo))
    c(r2 = cov_eo^2 / (sum(w * (e - me)^2) * sum(w * (o - mo)^2)),
      rmse = sqrt(sum(w * (o - e)^2) / sw))
}

bootstrap_fl_rq <- function(d, final_rq) {
    x_fl <- as.matrix(d$fl[, d$muestras])
    x_fl[is.na(x_fl)] <- 0
    g_fl <- d$fl$associated_gene[d$evaluadas_fl]
    g_rq <- d$rq$associated_gene[d$evaluadas_rq]
    genes <- unique(c(g_fl, g_rq))
    i_fl <- match(g_fl, genes)
    i_rq <- match(g_rq, genes)
    # Mismas réplicas de genes para todas las medidas y mezclas
    pesos <- rmultinom(B, length(genes), rep(1, length(genes)))

    out <- list()
    for (medida in c("counts", "TPM")) {
        mz_fl <- mezclas(x_fl, d$cols, d$evaluadas_fl, medida)
        mz_rq <- mezclas(final_rq, d$cols, d$evaluadas_rq, medida)
        for (mezcla in names(mz_fl)) {
            f <- mz_fl[[mezcla]]
            r <- mz_rq[[mezcla]]
            obs_fl <- r2_rmse_pesos(f$esperado, f$observado, rep(1, length(i_fl)))
            obs_rq <- r2_rmse_pesos(r$esperado, r$observado, rep(1, length(i_rq)))
            delta <- vapply(seq_len(B), function(b) {
                r2_rmse_pesos(r$esperado, r$observado, pesos[i_rq, b]) -
                    r2_rmse_pesos(f$esperado, f$observado, pesos[i_fl, b])
            }, numeric(2))
            out[[length(out) + 1]] <- data.frame(
                combinacion = d$combi, medida = medida, mezcla = mezcla,
                r2_fl = obs_fl[["r2"]], r2_rq = obs_rq[["r2"]],
                delta_r2 = obs_rq[["r2"]] - obs_fl[["r2"]],
                delta_r2_ic95_inf = quantile(delta["r2", ], 0.025, names = FALSE),
                delta_r2_ic95_sup = quantile(delta["r2", ], 0.975, names = FALSE),
                p_r2 = (1 + sum(delta["r2", ] <= 0)) / (B + 1),
                rmse_fl = obs_fl[["rmse"]], rmse_rq = obs_rq[["rmse"]],
                delta_rmse = obs_rq[["rmse"]] - obs_fl[["rmse"]],
                delta_rmse_ic95_inf = quantile(delta["rmse", ], 0.025, names = FALSE),
                delta_rmse_ic95_sup = quantile(delta["rmse", ], 0.975, names = FALSE),
                p_rmse = (1 + sum(delta["rmse", ] >= 0)) / (B + 1),
                genes = length(genes), B_bootstrap = B)
        }
    }
    bind_rows(out)
}

# ------------------------------------------------------------------------------
# Una combinación
# ------------------------------------------------------------------------------
analizar <- function(combi) {
    t0 <- Sys.time()
    d <- preparar(combi)

    # Comprobación: con las asignaciones reales se reproducen los counts de rq
    final_obs <- requantificar(d)
    rq_counts <- as.matrix(d$rq[, d$muestras])
    rq_counts[is.na(rq_counts)] <- 0
    dif <- max(abs(final_obs - rq_counts))
    observado <- evaluar(final_obs, d$cols, d$evaluadas_rq)

    nulas <- list()
    for (nulo in NULOS) {
        n <- preparar_nulo(d, nulo)
        res <- vector("list", B)
        for (i in seq_len(B)) {
            res[[i]] <- cbind(iter = i, evaluar(requantificar(d, target_azar(n)), d$cols, d$evaluadas_rq))
        }
        nulas[[nulo]] <- cbind(nulo = nulo, bind_rows(res))
    }
    boot <- bootstrap_fl_rq(d, final_obs)

    n_gen <- preparar_nulo(d, "gen")
    sin_alt <- 100 * mean(n_gen$tam[n_gen$grupo_fila] <= 1)
    message(sprintf("%s %s: %d filas reasignadas, %.1f%% sin alternativa en su gen, diferencia máx. con rq = %.3g (%.0f s)",
                    etiqueta, combi, length(d$target), sin_alt, dif,
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    list(observado = cbind(combinacion = combi, observado, dif_max_rq = dif,
                           filas_reasignadas = length(d$target), pct_sin_alternativa_gen = sin_alt),
         nulas = cbind(combinacion = combi, bind_rows(nulas)),
         bootstrap = boot)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(61)
res <- mclapply(combis, function(cb) tryCatch(analizar(cb), error = function(e) {
    message("ERROR en ", cb, ": ", conditionMessage(e)); NULL
}), mc.cores = min(nucleos, length(combis)), mc.set.seed = TRUE)
res <- res[!vapply(res, is.null, logical(1))]
if (!length(res)) stop("ninguna combinación terminó")

observado <- bind_rows(lapply(res, `[[`, "observado"))
nulas     <- bind_rows(lapply(res, `[[`, "nulas"))
boot      <- bind_rows(lapply(res, `[[`, "bootstrap"))

# Resultados previos de otras combinaciones en SALIDA
f_nulas <- file.path(salida, "distribuciones_nulas.csv")
f_obs   <- file.path(salida, "observado.csv")
f_boot  <- file.path(salida, "bootstrap_fl_rq.csv")
f_res   <- file.path(salida, "resumen_test_reasignacion.csv")
if (all(file.exists(c(f_nulas, f_obs, f_boot)))) {
    nuevas <- unique(observado$combinacion)
    otras <- function(f) read_csv(f, show_col_types = FALSE) %>%
        filter(!combinacion %in% nuevas) %>%
        select(-any_of("ejecucion"))
    prev_obs <- otras(f_obs)
    if (nrow(prev_obs)) {
        message("Se conservan los resultados previos de: ", paste(unique(prev_obs$combinacion), collapse = ", "))
        observado <- bind_rows(prev_obs, observado)
        nulas <- bind_rows(otras(f_nulas), nulas)
        boot <- bind_rows(otras(f_boot), boot %>% select(-any_of("p_r2_bh"), -any_of("p_rmse_bh")))
    }
}

if (any(observado$dif_max_rq > 1e-6)) {
    warning("La requantificación reproducida no coincide con rq en: ",
            paste(unique(observado$combinacion[observado$dif_max_rq > 1e-6]), collapse = ", "))
}

resumen <- nulas %>%
    inner_join(observado %>% select(combinacion, medida, mezcla, r2_obs = r2, rmse_obs = rmse),
               by = c("combinacion", "medida", "mezcla")) %>%
    group_by(nulo, combinacion, medida, mezcla) %>%
    summarise(r2_obs = first(r2_obs), r2_azar_media = mean(r2),
              r2_azar_ic95_inf = quantile(r2, 0.025), r2_azar_ic95_sup = quantile(r2, 0.975),
              r2_z = (first(r2_obs) - mean(r2)) / sd(r2),
              p_r2 = (1 + sum(r2 >= first(r2_obs))) / (n() + 1),
              rmse_obs = first(rmse_obs), rmse_azar_media = mean(rmse),
              rmse_azar_ic95_inf = quantile(rmse, 0.025), rmse_azar_ic95_sup = quantile(rmse, 0.975),
              rmse_z = (first(rmse_obs) - mean(rmse)) / sd(rmse),
              p_rmse = (1 + sum(rmse <= first(rmse_obs))) / (n() + 1),
              B = n(), .groups = "drop") %>%
    group_by(nulo, medida) %>%
    mutate(p_r2_bh = p.adjust(p_r2, "BH"), p_rmse_bh = p.adjust(p_rmse, "BH")) %>%
    ungroup() %>%
    left_join(observado %>% distinct(combinacion, filas_reasignadas, pct_sin_alternativa_gen, dif_max_rq),
              by = "combinacion") %>%
    mutate(ejecucion = etiqueta, .before = 1) %>%
    mutate(across(where(is.double), ~ signif(.x, 4)))

boot <- boot %>%
    group_by(medida) %>%
    mutate(p_r2_bh = p.adjust(p_r2, "BH"), p_rmse_bh = p.adjust(p_rmse, "BH")) %>%
    ungroup() %>%
    mutate(ejecucion = etiqueta, .before = 1) %>%
    mutate(across(where(is.double), ~ signif(.x, 4)))

write_csv(nulas %>% mutate(ejecucion = etiqueta, .before = 1), f_nulas)
write_csv(observado %>% mutate(ejecucion = etiqueta, .before = 1), f_obs)
write_csv(resumen, f_res)
write_csv(boot, f_boot)

# ------------------------------------------------------------------------------
# Figuras (TPM)
# ------------------------------------------------------------------------------
etiquetar <- function(df) {
    df %>% mutate(panel = paste(combinacion, mezcla),
                  seq = factor(sub("_.*", "", combinacion), levels = c("masseq", "isoseq", "ont")))
}
colores_nulo <- c("gen" = "#2a78d6", "global" = "#9e9e9e")
etiquetas_nulo <- c("gen" = "Al azar dentro del gen", "global" = "Al azar en todo el transcriptoma")

for (metrica in c("r2", "rmse")) {
    nombre <- if (metrica == "r2") "R²" else "RMSE"
    obs_col <- paste0(metrica, "_obs")
    p <- ggplot(etiquetar(filter(nulas, medida == "TPM")), aes(x = .data[[metrica]], fill = nulo)) +
        geom_histogram(bins = 40, alpha = 0.75, position = "identity", colour = NA) +
        geom_vline(data = etiquetar(filter(resumen, medida == "TPM", nulo == "gen")),
                   aes(xintercept = .data[[obs_col]]), colour = "#e34948", linewidth = 0.7) +
        facet_wrap(~ panel, scales = "free", ncol = 4) +
        scale_fill_manual(values = colores_nulo, labels = etiquetas_nulo) +
        labs(title = paste0(nombre, " de la reasignación del rescue frente a reasignaciones al azar (TPM, ",
                            etiqueta, ")"),
             subtitle = "Línea roja: reasignación de SQANTI3 rescue. Histogramas: reasignaciones al azar por modelo nulo",
             x = nombre, y = "Reasignaciones al azar", fill = NULL) +
        tema + theme(axis.text = element_text(size = 7))
    ggsave(paste0("nulas_", metrica, "_TPM.png"), p, path = salida, width = 13, height = 14, bg = "white")
}

# Visión global: diferencias con los datos filtrados (fl = 0). rq con su intervalo
# bootstrap y las reasignaciones al azar dentro del gen con su intervalo del 95%
fases <- bind_rows(
    boot %>% filter(medida == "TPM") %>%
        transmute(combinacion, mezcla, metrica = "R²", fl = r2_fl, rq = delta_r2,
                  rq_inf = delta_r2_ic95_inf, rq_sup = delta_r2_ic95_sup),
    boot %>% filter(medida == "TPM") %>%
        transmute(combinacion, mezcla, metrica = "RMSE", fl = rmse_fl, rq = delta_rmse,
                  rq_inf = delta_rmse_ic95_inf, rq_sup = delta_rmse_ic95_sup)) %>%
    left_join(bind_rows(
        resumen %>% filter(medida == "TPM", nulo == "gen") %>%
            transmute(combinacion, mezcla, metrica = "R²", inf = r2_azar_ic95_inf, sup = r2_azar_ic95_sup),
        resumen %>% filter(medida == "TPM", nulo == "gen") %>%
            transmute(combinacion, mezcla, metrica = "RMSE", inf = rmse_azar_ic95_inf, sup = rmse_azar_ic95_sup)),
        by = c("combinacion", "mezcla", "metrica")) %>%
    mutate(inf = inf - fl, sup = sup - fl) %>%
    etiquetar() %>%
    mutate(metrica = factor(ifelse(metrica == "R²", "Δ R² (mejor a la derecha)", "Δ RMSE (mejor a la izquierda)"),
                            levels = c("Δ R² (mejor a la derecha)", "Δ RMSE (mejor a la izquierda)")))
p <- ggplot(fases, aes(y = panel)) +
    geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.5) +
    geom_errorbar(aes(xmin = inf, xmax = sup, colour = "Al azar dentro del gen (IC 95%)"),
                  orientation = "y", width = 0.35, linewidth = 0.8,
                  position = position_nudge(y = 0.18)) +
    geom_errorbar(aes(xmin = rq_inf, xmax = rq_sup, colour = "Rescue con requant (bootstrap IC 95%)"),
                  orientation = "y", width = 0.35, linewidth = 0.8,
                  position = position_nudge(y = -0.18)) +
    geom_point(aes(x = rq, colour = "Rescue con requant (bootstrap IC 95%)"), size = 2.2,
               position = position_nudge(y = -0.18)) +
    facet_grid(seq ~ metrica, scales = "free", space = "free_y") +
    scale_colour_manual(values = c("Al azar dentro del gen (IC 95%)" = colores_nulo[["gen"]],
                                   "Rescue con requant (bootstrap IC 95%)" = "#e34948")) +
    labs(title = paste0("Cambio respecto de los datos filtrados (TPM, ", etiqueta, ")"),
         subtitle = "0 = datos filtrados (fl). Rojo: rescue con requantificación. Azul: reasignaciones al azar dentro del gen",
         x = NULL, y = NULL, colour = NULL) +
    tema
ggsave("fases_TPM.png", p, path = salida, width = 11, height = 9, bg = "white")

print(as.data.frame(resumen %>% filter(medida == "TPM", nulo == "gen") %>%
                        select(combinacion, mezcla, r2_obs, r2_azar_media, p_r2_bh,
                               rmse_obs, rmse_azar_media, p_rmse_bh)), row.names = FALSE)
print(as.data.frame(boot %>% filter(medida == "TPM") %>%
                        select(combinacion, mezcla, r2_fl, r2_rq, p_r2_bh, rmse_fl, rmse_rq, p_rmse_bh)),
      row.names = FALSE)
cat("Salidas en", salida, "\n")
