#!/usr/bin/env Rscript
# Checks of the release-6.1 pull requests of DonGazpachoWR/SQANTI3 on the real
# data (13 combinations), from the outputs downloaded by master_workflow_r61.sh:
#
#   PR 1  min_intron_length: NA exactly for mono-exonic isoforms; distribution by
#         structural category and share of multi-exonic isoforms with an intron
#         shorter than 100 bp.
#         Requantification by default: counts are conserved per sample.
#         Rescue report: overview numbers of every combination.
#   PR 2  prevalence, prevalence_K, prevalence_B written by QC equal the values
#         recomputed from the per-sample counts; filter_result equals kb2.
#   PR 3  evidence_check of every reference target equals the one recomputed
#         independently from the artifacts' counts (OR of >= 2 samples in K or B).
#
# Usage: Rscript verificar_prs_r61.r <DATA_R61> <DATA_KB2> <OUT_DIR>

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
ruta_r61 <- args[1]; ruta_kb2 <- args[2]; outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
                 "tema_release61.r"))

combis <- c(paste0("isoseq_", c("isocall", "isoseq", "bambu", "flair", "isoquant")),
            paste0("masseq_", c("isocall", "isoseq", "bambu", "flair", "isoquant")),
            paste0("ont_", c("bambu", "flair", "isoquant")))
MIN_EXPRESSION <- 0     # --min_expression de SQANTI3 QC
# Misma regla que is_expressed() de SQANTI3: con 0, count > 0; si no, count >= MIN_EXPRESSION
expresado <- function(x) if (MIN_EXPRESSION == 0) x > 0 else x >= MIN_EXPRESSION
UMBRAL_GRUPO  <- 2      # prevalence_K >= 2 OR prevalence_B >= 2 (filter_r61.json)

leer <- function(...) {
    f <- file.path(...)
    if (!file.exists(f)) stop("no existe ", f)
    read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "",
               comment.char = "", check.names = FALSE, row.names = NULL)
}

replicas <- function(seq) {
    n <- if (seq == "masseq") 3 else 5
    list(K = paste0("K3", 1:n), B = paste0("B3", 1:n))
}

verificar <- function(combi) {
    seq <- sub("_.*", "", combi)
    rep <- replicas(seq)
    counts <- leer(ruta_r61, paste0("counts_", combi, ".tsv"))
    muestras <- colnames(counts)[-1]
    qc <- leer(ruta_r61, paste0("default_", combi, "_classification.txt"))
    fl <- leer(ruta_r61, paste0("rules_default_", combi, "_RulesFilter_classification.txt"))
    tabla <- leer(ruta_r61, paste0("rq_", combi, "_rescue_table.tsv"))
    reas <- leer(ruta_r61, paste0("rq_", combi, "_reassigned_counts.tsv"))
    resumen <- leer(ruta_r61, paste0("rq_", combi, "_rescue_summary.tsv"))

    # ---------------------------------------------------------------- PR 1
    mono <- qc$exons == 1
    intron_na_ok <- all(is.na(qc$min_intron_length) == mono)
    multi <- qc[!mono, ]
    pct_corto <- 100 * mean(multi$min_intron_length < 100)

    # Only rows that are transcripts of the QC (ont_bambu has a non-transcript row "0")
    comunes <- intersect(muestras, colnames(reas))
    en_qc <- counts[[1]] %in% qc$isoform
    antes <- colSums(counts[en_qc, comunes, drop = FALSE], na.rm = TRUE)
    despues <- colSums(reas[, comunes, drop = FALSE], na.rm = TRUE)
    conserva <- isTRUE(all.equal(unname(antes), unname(despues), tolerance = 1e-6))

    ov <- resumen %>% filter(section == "overview")
    ov_n <- setNames(ov$count, ov$category)

    # ---------------------------------------------------------------- PR 2
    # Isoforms missing from --fl_count have NA counts and NA prevalence
    sin_conteos <- rowSums(!is.na(qc[, muestras])) == 0
    q <- qc[!sin_conteos, muestras]; q[is.na(q)] <- 0
    det <- expresado(q)
    con <- qc[!sin_conteos, ]
    prev_ok <- all(is.na(qc$prevalence[sin_conteos])) &&
        all(is.na(qc$prevalence_K[sin_conteos])) && all(is.na(qc$prevalence_B[sin_conteos])) &&
        all(con$prevalence == rowSums(det)) &&
        all(con$prevalence_K == rowSums(det[, rep$K, drop = FALSE])) &&
        all(con$prevalence_B == rowSums(det[, rep$B, drop = FALSE]))
    pasa_esperado <- !is.na(qc$prevalence_K) &
        (qc$prevalence_K >= UMBRAL_GRUPO | qc$prevalence_B >= UMBRAL_GRUPO)
    fl_isoform <- fl$filter_result == "Isoform"
    # Every isoform passing the filter must meet the prevalence requisite
    prev_filtro_ok <- all(pasa_esperado[match(fl$isoform[fl_isoform], qc$isoform)])

    kb2 <- tryCatch(leer(ruta_kb2, paste0("rules_default_", combi, "_RulesFilter_classification.txt")),
                    error = function(e) NULL)
    if (is.null(kb2)) {
        n_dif_kb2 <- NA
    } else {
        m <- merge(fl[, c("isoform", "filter_result")], kb2[, c("isoform", "filter_result")],
                   by = "isoform", all = TRUE, suffixes = c("_r61", "_kb2"))
        n_dif_kb2 <- sum(is.na(m$filter_result_r61) | is.na(m$filter_result_kb2) |
                         m$filter_result_r61 != m$filter_result_kb2)
    }

    # ---------------------------------------------------------------- PR 3
    ev <- tabla %>% filter(origin == "reference", evidence_check %in% c("pass", "failed"))
    pares <- ev %>% distinct(artifact, assigned_transcript) %>%
        group_by(artifact) %>% mutate(share = 1 / n()) %>% ungroup()
    m_art <- as.matrix(fl[match(pares$artifact, fl$isoform), c(rep$K, rep$B)])
    m_art[is.na(m_art)] <- 0
    agg <- rowsum(m_art * pares$share, pares$assigned_transcript)
    pasa_rec <- rowSums(expresado(agg[, rep$K, drop = FALSE])) >= UMBRAL_GRUPO |
        rowSums(expresado(agg[, rep$B, drop = FALSE])) >= UMBRAL_GRUPO
    estado <- ev %>% distinct(assigned_transcript, evidence_check)
    pasa_sq <- estado$evidence_check[match(rownames(agg), estado$assigned_transcript)] == "pass"
    n_dif_ev <- sum(pasa_rec != pasa_sq)
    cuenta <- table(factor(tabla$evidence_check, levels = c("pass", "failed", "reassigned", "not_required")))

    list(
        resumen = data.frame(
            combinacion = combi,
            pr1_intron_na_solo_monoexonicas = intron_na_ok,
            pr1_pct_multiexonicas_intron_menor_100 = round(pct_corto, 2),
            pr1_requant_conserva_conteos = conserva,
            pr1_artefactos = ov_n[["artifacts"]],
            pr1_artefactos_rescatados = ov_n[["rescued"]],
            pr1_referencias_anadidas = ov_n[["added_reference_transcripts"]],
            pr2_isoformas_sin_conteos = sum(sin_conteos),
            pr2_prevalencias_correctas = prev_ok,
            pr2_isoformas_cumplen_regla = prev_filtro_ok,
            pr2_isoformas_isoform = sum(fl_isoform),
            pr2_diferencias_filter_vs_kb2 = n_dif_kb2,
            pr3_targets_referencia = nrow(agg),
            pr3_targets_pasan = sum(pasa_sq),
            pr3_diferencias_recalculo = n_dif_ev,
            pr3_filas_pass = cuenta[["pass"]],
            pr3_filas_failed = cuenta[["failed"]],
            pr3_filas_reassigned = cuenta[["reassigned"]],
            check.names = FALSE),
        intrones = data.frame(combinacion = combi, seq = seq,
                              structural_category = multi$structural_category,
                              min_intron_length = multi$min_intron_length),
        evidencia = data.frame(combinacion = combi, estado = names(cuenta), n = as.integer(cuenta))
    )
}

res <- lapply(combis, function(cb) tryCatch(verificar(cb), error = function(e) {
    message("ERROR en ", cb, ": ", conditionMessage(e)); NULL
}))
res <- res[!vapply(res, is.null, logical(1))]
tabla <- bind_rows(lapply(res, `[[`, "resumen"))
write_csv(tabla, file.path(outdir, "verificacion_r61.csv"))
print(as.data.frame(tabla), row.names = FALSE)

# Overall verdict
checks_logicos <- c("pr1_intron_na_solo_monoexonicas", "pr1_requant_conserva_conteos",
                    "pr2_prevalencias_correctas", "pr2_isoformas_cumplen_regla")
fallos <- sum(!as.matrix(tabla[, checks_logicos])) +
    sum(tabla$pr2_diferencias_filter_vs_kb2 != 0, na.rm = TRUE) +
    sum(tabla$pr3_diferencias_recalculo != 0)
cat(sprintf("\nCombinaciones verificadas: %d de %d. Comprobaciones fallidas: %d\n",
            nrow(tabla), length(combis), fallos))

# ---------------------------------------------------------------- figures
cat_labels <- c(`full-splice_match` = "FSM", `incomplete-splice_match` = "ISM",
                novel_in_catalog = "NIC", novel_not_in_catalog = "NNC", genic = "Genic",
                antisense = "Antisense", fusion = "Fusion", intergenic = "Intergenic",
                genic_intron = "Genic intron")
intr <- bind_rows(lapply(res, `[[`, "intrones")) %>%
    mutate(cat = factor(cat_labels[structural_category], levels = cat_labels)) %>%
    filter(!is.na(cat))

p1 <- ggplot(intr, aes(x = min_intron_length, colour = cat)) +
    stat_ecdf(linewidth = 0.6) +
    geom_vline(xintercept = 100, linetype = "dashed", colour = "grey40") +
    scale_x_log10() +
    facet_wrap(~ seq) +
    labs(title = "PR 1 - min_intron_length de las isoformas multiexónicas",
         subtitle = "Distribución acumulada por categoría estructural; línea discontinua: 100 pb",
         x = "Intrón más corto (pb, escala log)", y = "Fracción acumulada", colour = "Categoría") +
    tema
ggsave("pr1_min_intron_length_ecdf.png", p1, path = outdir, width = 10, height = 5, bg = "white")

corto <- intr %>% group_by(combinacion, cat) %>%
    summarise(pct = 100 * mean(min_intron_length < 100), .groups = "drop")
p2 <- ggplot(corto, aes(x = cat, y = combinacion, fill = pct)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.1f", pct)), size = 3) +
    scale_fill_gradient(low = "#e3eefb", high = "#0b3f86", name = "% < 100 pb") +
    labs(title = "PR 1 - % de isoformas multiexónicas con algún intrón < 100 pb",
         x = NULL, y = NULL) +
    tema + theme(legend.position = "right")
ggsave("pr1_intron_menor_100.png", p2, path = outdir, width = 9, height = 6, bg = "white")

evid <- bind_rows(lapply(res, `[[`, "evidencia")) %>%
    mutate(estado = factor(estado, levels = c("pass", "failed", "reassigned", "not_required")))
p3 <- ggplot(evid, aes(x = combinacion, y = n, fill = estado)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(pass = "#2a78d6", failed = "#e34948",
                                 reassigned = "#eb6834", not_required = "#9e9e9e")) +
    labs(title = "PR 3 - Control de evidencia por combinación",
         subtitle = "Filas de la tabla de rescue según evidence_check",
         x = NULL, y = "Filas", fill = NULL) +
    tema
ggsave("pr3_evidence_check.png", p3, path = outdir, width = 9, height = 6, bg = "white")

cat("Figuras y tabla en", outdir, "\n")
