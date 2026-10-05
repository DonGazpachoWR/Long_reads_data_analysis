# Requantification of SQANTI3 rescue reproduced in R (sq_requant.py and
# requant_helpers.py of DonGazpachoWR/SQANTI3), and the evaluation of the tissue
# mixtures of correlacion_mezclas_individual.r. Sourced by
# test_reasignacion_azar_r61.r and verificar_prs_r61.r.
#
# Requantification: counts of the --counts file restricted to the QC isoforms;
# each artifact is split among its targets in proportion to the counts of each
# target in the sample (evenly if none has counts); with integer counts, floor of
# each share and the remainder to the first target of the artifact; with
# fractional counts (detected over all the artifacts), proportional shares.

columnas <- function(seq) {
    n <- if (seq == "masseq") 3 else 5
    cols <- list(K = paste0("K3", 1:n), B = paste0("B3", 1:n),
                 B20K80 = paste0("B20K80_", 1:n))
    if (seq != "masseq") cols$B80K20 <- paste0("B80K20_", 1:n)
    cols
}

leer_tsv <- function(f, ...) {
    if (!file.exists(f)) stop("no existe ", f)
    suppressWarnings(readr::read_tsv(f, show_col_types = FALSE, progress = FALSE, ...))
}

# Counts of requantification (sq_requant.load_counts): --counts file restricted
# to the QC isoforms, as a matrix isoform x sample
cargar_counts <- function(datos, combi, muestras) {
    counts <- leer_tsv(file.path(datos, paste0("counts_", combi, ".tsv")))
    names(counts)[1] <- "isoform"
    qc <- leer_tsv(file.path(datos, paste0("default_", combi, "_classification.txt")), col_select = "isoform")
    counts <- counts[counts$isoform %in% qc$isoform, ]
    m <- as.matrix(counts[, muestras])
    m[is.na(m)] <- 0
    rownames(m) <- counts$isoform
    m
}

# Prepares the requantification of the rows of a rescue table.
#   m_counts:  cargar_counts()
#   fl:        filter classification (isoform, filter_result)
#   tabla:     rescue table rows to requantify (artifact, assigned_transcript)
#   destinos:  isoforms whose final counts are computed; must contain every target
preparar_requant <- function(m_counts, fl, tabla, destinos) {
    validas    <- fl$isoform[fl$filter_result == "Isoform"]
    artefactos <- fl$isoform[fl$filter_result == "Artifact"]
    tabla <- tabla[tabla$artifact %in% artefactos, ]
    muestras <- colnames(m_counts)

    base <- matrix(0, length(destinos), length(muestras), dimnames = list(destinos, muestras))
    en_base <- destinos %in% validas & destinos %in% rownames(m_counts)
    base[en_base, ] <- m_counts[destinos[en_base], ]

    fuente <- matrix(0, nrow(tabla), length(muestras), dimnames = list(NULL, muestras))
    tiene <- tabla$artifact %in% rownames(m_counts)
    fuente[tiene, ] <- m_counts[tabla$artifact[tiene], ]

    target <- match(tabla$assigned_transcript, destinos)
    if (anyNA(target)) stop(sum(is.na(target)), " targets de la tabla no están entre los destinos")

    # distribute_integer_counts() detects fractional counts over every artifact
    # (rescued and sent to the residual)
    todos <- m_counts[intersect(artefactos, rownames(m_counts)), , drop = FALSE]
    art <- match(tabla$artifact, unique(tabla$artifact))
    list(base = base, fuente = fuente, target = target, art = art,
         fraccional = any(todos %% 1 != 0), primero = !duplicated(art),
         tabla = tabla)
}

requantificar <- function(d, target = d$target) {
    final <- d$base
    if (!length(target)) return(final)
    w <- d$base[target, , drop = FALSE]
    total <- rowsum(w, d$art, reorder = FALSE)[as.character(d$art), , drop = FALSE]
    n_k <- tabulate(d$art)[d$art]
    frac <- ifelse(total == 0, 1 / n_k, w / ifelse(total == 0, 1, total))
    suma <- d$fuente * frac
    if (!d$fraccional) {
        suma <- floor(suma)
        resto <- d$fuente - rowsum(suma, d$art, reorder = FALSE)[as.character(d$art), , drop = FALSE]
        suma[d$primero, ] <- suma[d$primero, ] + resto[d$primero, ]
    }
    anadir <- rowsum(suma, target)
    filas <- as.integer(rownames(anadir))
    final[filas, ] <- final[filas, ] + anadir
    final
}

# Expected and observed mixture (log10 + 0.01 of the replicate means) of the
# evaluated rows (FSM not novel), as in correlacion_mezclas_individual.r.
# TPM is computed over every row of x before selecting the evaluated ones.
mezclas <- function(x, cols, evaluadas, medida) {
    if (medida == "TPM") {
        a <- colSums(x) / 1e6
        a[a == 0] <- 1
        x <- sweep(x, 2, a, "/")
    }
    x <- x[evaluadas, , drop = FALSE]
    k <- rowMeans(x[, cols$K, drop = FALSE])
    b <- rowMeans(x[, cols$B, drop = FALSE])
    out <- list()
    for (mezcla in intersect(c("B20K80", "B80K20"), names(cols))) {
        pb <- if (mezcla == "B20K80") 0.2 else 0.8
        out[[mezcla]] <- list(esperado = log10(b * pb + k * (1 - pb) + 0.01),
                              observado = log10(rowMeans(x[, cols[[mezcla]], drop = FALSE]) + 0.01))
    }
    out
}

evaluar <- function(x, cols, evaluadas) {
    out <- list()
    for (medida in c("counts", "TPM")) {
        mz <- mezclas(x, cols, evaluadas, medida)
        for (mezcla in names(mz)) {
            e <- mz[[mezcla]]$esperado
            o <- mz[[mezcla]]$observado
            out[[length(out) + 1]] <- data.frame(medida = medida, mezcla = mezcla,
                                                 r2 = cor(e, o)^2, rmse = sqrt(mean((o - e)^2)))
        }
    }
    dplyr::bind_rows(out)
}
