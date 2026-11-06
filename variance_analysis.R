#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Mean and variance effects across MMPC two-cohort comparisons
# Riordan et al.
#
# Base R only - no packages required.
#
# Usage:  Rscript variance_analysis.R [trait_index.csv] [data_dir] [out.csv]
#                                    [nperm] [seed]
#
#   nperm  permutations per statistic (default 9999)
#   seed   integer seed, or "none" to draw a fresh one each run (default 9999
#          permutations with seed 20261023, which reproduces the supplied
#          results.csv).  Permutation p-values are Monte Carlo estimates, so
#          a different seed gives slightly different values; that variation is
#          expected and is roughly +/- 0.002 for a p near 0.05 at 9999
#          permutations.
#
# Inputs
#   trait_index.csv  one row per trait, giving the source file, the assay
#                    column, the two cohorts, anything held fixed, and any
#                    excluded animals.
#   data_dir         the MMPC per-animal files named in trait_index.csv
#
# Both inputs are read as latin1, which is how the MMPC files are encoded.
# Within cohort_A / cohort_B / held_fixed, "column=value" terms are separated
# by " | " rather than ";" because one strain name (B6;129S4-...) contains a
# semicolon.
#
# Output
#   one row per trait: cohort sizes, means, MADs, the six effect metrics,
#   permutation p-values for each, and three tests (Welch's t on the values,
#   Welch's t on absolute deviations from the cohort median, and
#   Brown-Forsythe).
# ---------------------------------------------------------------------------

args      <- commandArgs(trailingOnly = TRUE)
index_f   <- if (length(args) >= 1) args[1] else "trait_index.csv"
data_dir  <- if (length(args) >= 2) args[2] else "data"
out_f     <- if (length(args) >= 3) args[3] else "results.csv"
NPERM     <- if (length(args) >= 4) as.integer(args[4]) else 9999
seed_arg  <- if (length(args) >= 5) args[5] else "20261023"

## R >= 3.6.0.  RNGkind is set explicitly because the default sampling
## algorithm changed in 3.6.0; without it a given seed yields different
## permutations on older versions.
RNGkind(sample.kind = "Rejection")
if (identical(tolower(seed_arg), "none")) {
  seed_used <- as.integer(Sys.time()) %% .Machine$integer.max
} else {
  seed_used <- as.integer(seed_arg)
}
set.seed(seed_used)
message("permutations: ", NPERM, "   seed: ", seed_used)

## --- effect metrics --------------------------------------------------------
## mad(constant = 1) is the raw median absolute deviation, unscaled.
m1 <- function(x, g) mean(x[g == "A"]);  m2 <- function(x, g) mean(x[g == "B"])
d1 <- function(x, g) mad(x[g == "A"], constant = 1)
d2 <- function(x, g) mad(x[g == "B"], constant = 1)

meanDif  <- function(x, g) abs(m1(x, g) - m2(x, g))
meanDif2 <- function(x, g) abs(m1(x, g) - m2(x, g)) / (m1(x, g) + m2(x, g))
meanRat  <- function(x, g) max(m1(x, g) / m2(x, g), m2(x, g) / m1(x, g))
madDif   <- function(x, g) abs(d1(x, g) - d2(x, g))
madDif2  <- function(x, g) abs(d1(x, g) - d2(x, g)) / (d1(x, g) + d2(x, g))
madRat   <- function(x, g) max(d1(x, g) / d2(x, g), d2(x, g) / d1(x, g))

## meanRat and madRat take the larger of the two ratios, so both are >= 1 and
## neither records which cohort is the more variable one.

## --- permutation test ------------------------------------------------------
## All six statistics are evaluated on the same shuffled labels in a single
## pass, so each permutation computes the two means and two MADs once rather
## than twelve times.  p = (#{null >= observed} + 1) / (nperm + 1), one-sided,
## which suits these statistics because each is a non-negative magnitude.
STATS <- c("meanDif", "meanDif2", "meanRat", "madDif", "madDif2", "madRat")

six <- function(a, b) {
  ma <- mean(a); mb <- mean(b)
  da <- mad(a, constant = 1); db <- mad(b, constant = 1)
  c(abs(ma - mb),
    abs(ma - mb) / (ma + mb),
    max(ma / mb, mb / ma),
    abs(da - db),
    abs(da - db) / (da + db),
    max(da / db, db / da))
}

perm_all <- function(x, g, nperm = NPERM) {
  obs <- suppressWarnings(six(x[g == "A"], x[g == "B"]))
  cnt <- numeric(length(obs))
  nA  <- sum(g == "A"); n <- length(x)
  for (k in seq_len(nperm)) {
    ia  <- sample.int(n, nA)
    v   <- suppressWarnings(six(x[ia], x[-ia]))
    cnt <- cnt + (!is.na(v) & !is.na(obs) & v >= obs)
  }
  p <- (cnt + 1) / (nperm + 1)
  p[!is.finite(obs)] <- NA_real_
  names(obs) <- STATS; names(p) <- STATS
  list(obs = obs, p = p)
}

## --- tests -----------------------------------------------------------------
safe_p <- function(expr) tryCatch(expr, error = function(e) NA_real_)
devs   <- function(x, g) abs(x - ave(x, g, FUN = median))   # |value - cohort median|

welch_p <- function(x, g) safe_p(t.test(x[g == "A"], x[g == "B"], var.equal = FALSE)$p.value)
mad_t_p <- function(x, g) { d <- devs(x, g)
  safe_p(t.test(d[g == "A"], d[g == "B"], var.equal = FALSE)$p.value) }
bf_p    <- function(x, g) { d <- devs(x, g)   # Brown-Forsythe = pooled t on the same deviations
  safe_p(t.test(d[g == "A"], d[g == "B"], var.equal = TRUE)$p.value) }

## --- cohort selection ------------------------------------------------------
parse_spec <- function(s) {
  s <- trimws(s); if (!nzchar(s)) return(list())
  parts <- strsplit(s, " | ", fixed = TRUE)[[1]]
  out <- list()
  for (p in parts) {
    p <- trimws(p); if (!nzchar(p)) next
    k <- sub("=.*$", "", p); v <- sub("^[^=]*=", "", p)
    out[[trimws(k)]] <- trimws(v)
  }
  out
}
matches <- function(df, spec) {
  keep <- rep(TRUE, nrow(df))
  for (k in names(spec)) {
    if (!k %in% names(df)) return(rep(FALSE, nrow(df)))
    keep <- keep & !is.na(df[[k]]) & trimws(as.character(df[[k]])) == spec[[k]]
  }
  keep
}

## The MMPC files are latin1 and carry a few stray leading bytes, so they are
## read as raw bytes and converted explicitly rather than via fileEncoding.
read_csv_latin1 <- function(path) {
  b  <- readBin(path, "raw", file.size(path))
  s  <- iconv(rawToChar(b), from = "latin1", to = "UTF-8")
  df <- read.csv(text = strsplit(s, "\r?\n")[[1]], check.names = FALSE,
                 stringsAsFactors = FALSE, colClasses = "character")
  names(df) <- trimws(names(df))
  df
}
read_source <- read_csv_latin1

## --- main ------------------------------------------------------------------
idx   <- read_csv_latin1(index_f)
cache <- new.env(parent = emptyenv())
res   <- vector("list", nrow(idx))

for (i in seq_len(nrow(idx))) {
  r <- idx[i, ]
  if (!exists(r$source_file, envir = cache))
    assign(r$source_file, read_source(file.path(data_dir, r$source_file)), envir = cache)
  df <- get(r$source_file, envir = cache)

  fixed <- matches(df, parse_spec(r$held_fixed))
  selA  <- fixed & matches(df, parse_spec(r$cohort_A))
  selB  <- fixed & matches(df, parse_spec(r$cohort_B))

  drop <- trimws(strsplit(r$excluded_animals, " | ", fixed = TRUE)[[1]])
  drop <- drop[nzchar(drop)]
  if (length(drop) && "Animal" %in% names(df)) {
    isdrop <- trimws(as.character(df[["Animal"]])) %in% drop
    selA <- selA & !isdrop; selB <- selB & !isdrop
  }

  num <- suppressWarnings(as.numeric(df[[r$assay_column]]))
  A <- num[selA]; A <- A[!is.na(A)]
  B <- num[selB]; B <- B[!is.na(B)]
  x <- c(A, B); g <- factor(rep(c("A", "B"), c(length(A), length(B))), levels = c("A", "B"))

  pa <- perm_all(x, g)
  pa <- perm_all(x, g)
  res[[i]] <- data.frame(
    row = r$row, study_label = r$study_label, trait = r$trait_published,
    comparison_type = r$comparison_type, sex = r$sex,
    n_A = length(A), n_B = length(B),
    mean_A = mean(A), mean_B = mean(B),
    MAD_A  = mad(A, constant = 1), MAD_B = mad(B, constant = 1),
    meanDif  = pa$obs[["meanDif"]],  p_meanDif  = pa$p[["meanDif"]],
    meanDif2 = pa$obs[["meanDif2"]], p_meanDif2 = pa$p[["meanDif2"]],
    meanRat  = pa$obs[["meanRat"]],  p_meanRat  = pa$p[["meanRat"]],
    madDif   = pa$obs[["madDif"]],   p_madDif   = pa$p[["madDif"]],
    madDif2  = pa$obs[["madDif2"]],  p_madDif2  = pa$p[["madDif2"]],
    madRat   = pa$obs[["madRat"]],   p_madRat   = pa$p[["madRat"]],
    p_welch_t = welch_p(x, g), p_mad_t = mad_t_p(x, g), p_brown_forsythe = bf_p(x, g),
    stringsAsFactors = FALSE)

  if (i %% 25 == 0) message("  ", i, " / ", nrow(idx))
}

out <- do.call(rbind, res)
write.csv(out, out_f, row.names = FALSE)
message("wrote ", out_f, ": ", nrow(out), " traits  (", NPERM,
        " permutations, seed ", seed_used, ")")
