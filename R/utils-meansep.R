## Multiple comparison of means ---------------------------------------------
##
## Every comparison in CANE is driven by an explicitly supplied error mean
## square and its degrees of freedom. That is deliberate: in a split plot or
## strip plot design the correct denominator differs between strata, and a
## procedure that simply reads the residual mean square off a fitted model
## object will test main plots against the subplot error. Passing MSerror and
## DFerror explicitly makes the stratum used for each comparison visible in
## the returned object.

#' @keywords internal
#' @noRd
.cane_crit <- function(test, ntr, MSerror, DFerror, r, alpha) {
  ## Returns a function of the span (number of means in a range) giving the
  ## critical difference. stats::qtukey() is computationally expensive and
  ## .cane_letters() calls this function once per pair of means, so each
  ## distinct critical value is computed once and then reused.
  sem <- sqrt(MSerror / r)                 # standard error of a treatment mean
  switch(as.character(test),
    "1" = {
      cd <- suppressWarnings(stats::qtukey(1 - alpha, nmeans = ntr, df = DFerror)) * sem
      function(span) cd
    },
    "2" = {
      cache <- rep(NA_real_, max(2L, ntr))
      function(span) {
        span <- max(2L, as.integer(span))
        if (span > length(cache)) cache <<- c(cache, rep(NA_real_, span - length(cache)))
        if (is.na(cache[span])) {
          ap <- 1 - (1 - alpha)^(span - 1)
          cache[span] <<- suppressWarnings(stats::qtukey(1 - ap, nmeans = span, df = DFerror)) * sem
        }
        cache[span]
      }
    },
    "3" = {
      cd <- stats::qt(1 - alpha / 2, df = DFerror) * sqrt(2 * MSerror / r)
      function(span) cd
    }
  )
}

#' @keywords internal
#' @noRd
.cane_letters <- function(m, crit) {
  ## m must already be sorted in decreasing order.
  k <- length(m)
  if (k == 1L) return(stats::setNames("a", names(m)))
  ends <- integer(k)
  for (i in seq_len(k)) {
    j <- i
    for (jj in i:k) {
      if (m[[i]] - m[[jj]] <= crit(jj - i + 1L) + 1e-10) j <- jj else break
    }
    ends[i] <- j
  }
  ## keep only intervals that are not contained in an earlier/wider one
  keep <- rep(TRUE, k)
  for (i in seq_len(k)) {
    for (i2 in seq_len(k)) {
      if (i2 != i && i2 <= i && ends[i2] >= ends[i]) {
        if (i2 < i || ends[i2] > ends[i]) { keep[i] <- FALSE; break }
      }
    }
  }
  starts <- which(keep)
  lets <- .cane_letter_seq(length(starts))
  out <- character(k)
  for (g in seq_along(starts)) {
    idx <- starts[g]:ends[starts[g]]
    out[idx] <- paste0(out[idx], lets[g])
  }
  stats::setNames(out, names(m))
}

#' @keywords internal
#' @noRd
.cane_letter_seq <- function(n) {
  if (n <= 26L) return(letters[seq_len(n)])
  reps <- ceiling(n / 26)
  unlist(lapply(seq_len(reps), function(i) strrep(letters, i)))[seq_len(n)]
}

#' @keywords internal
#' @noRd
.cane_mean_sep <- function(y, trt, MSerror, DFerror, alpha, Mult_Comp_Test,
                           stratum = "residual", scale = "original",
                           back_factor = NULL) {

  trt <- droplevels(factor(trt))
  ok  <- is.finite(y) & !is.na(trt)
  y <- y[ok]; trt <- trt[ok]
  ntr <- nlevels(trt)
  if (ntr < 2L) .cane_stop("Multiple comparison requires at least two treatment levels.")
  if (!is.finite(MSerror) || MSerror <= 0 || !is.finite(DFerror) || DFerror < 1)
    .cane_stop("Multiple comparison needs a positive error mean square and at least one ",
               "error degree of freedom; the ", stratum, " stratum supplied MSerror = ",
               signif(MSerror, 6), " on ", DFerror, " df. This normally means the design ",
               "has no replication in that stratum.")

  n   <- as.vector(table(trt))
  mu  <- as.vector(tapply(y, trt, mean))
  sdv <- as.vector(tapply(y, trt, stats::sd))
  mn  <- as.vector(tapply(y, trt, min))
  mx  <- as.vector(tapply(y, trt, max))
  q   <- t(vapply(split(y, trt), stats::quantile, numeric(3), probs = c(.25, .5, .75),
                  names = FALSE, USE.NAMES = FALSE))
  names(mu) <- levels(trt)

  balanced <- length(unique(n)) == 1L
  r <- if (balanced) n[1L] else length(n) / sum(1 / n)   # harmonic mean
  if (!balanced)
    warning("Treatment means are based on unequal numbers of observations (", 
            paste(range(n), collapse = " to "),
            "). The critical difference was computed from the harmonic mean replication (",
            signif(r, 5), "); it is approximate for unbalanced data.", call. = FALSE)

  crit <- .cane_crit(Mult_Comp_Test, ntr, MSerror, DFerror, r, alpha)
  if (!all(is.finite(vapply(seq.int(2L, max(2L, ntr)), crit, numeric(1)))))
    .cane_stop("The critical value of the ", c("Tukey", "Duncan", "LSD")[Mult_Comp_Test],
               " procedure cannot be computed for the ", stratum, " stratum, which has only ",
               DFerror, " error degree(s) of freedom (the studentized range distribution ",
               "is not available on fewer than 2). Use more environments or replications, ",
               "or Mult_Comp_Test = 3 for a pre-planned comparison.")
  o  <- order(mu, decreasing = TRUE)
  ms <- mu[o]
  grp <- .cane_letters(ms, crit)

  testname <- c("Tukey's HSD", "Duncan's multiple range test (DMRT)",
                "Fisher's LSD")[Mult_Comp_Test]
  critval  <- c(stats::qtukey(1 - alpha, ntr, DFerror),
                NA_real_,
                stats::qt(1 - alpha / 2, DFerror))[Mult_Comp_Test]
  msd <- crit(ntr)

  means <- data.frame(
    Treatment = levels(trt), Mean = mu, SD = sdv, r = n,
    Min = mn, Max = mx, Q25 = q[, 1], Q50 = q[, 2], Q75 = q[, 3],
    row.names = NULL, stringsAsFactors = FALSE)

  groups <- data.frame(
    Treatment = names(ms), Mean = as.vector(ms), Groups = as.vector(grp),
    row.names = NULL, stringsAsFactors = FALSE)

  ## back-transformed means, when the analysis was run on the Aitken scale
  if (!is.null(back_factor)) {
    bt <- back_factor[levels(trt)]
    means$Mean_original_scale   <- unname(bt)
    groups$Mean_original_scale  <- unname(back_factor[groups$Treatment])
  }

  ## all pairwise differences with the stratum-specific standard error
  pr <- utils::combn(levels(trt), 2L)
  sed <- sqrt(MSerror * (1 / n[match(pr[1, ], levels(trt))] +
                         1 / n[match(pr[2, ], levels(trt))]))
  dif <- mu[pr[1, ]] - mu[pr[2, ]]
  tst <- abs(dif) / sed
  pval <- switch(as.character(Mult_Comp_Test),
    "1" = stats::ptukey(tst * sqrt(2), nmeans = ntr, df = DFerror, lower.tail = FALSE),
    "2" = rep(NA_real_, length(dif)),
    "3" = 2 * stats::pt(tst, df = DFerror, lower.tail = FALSE))
  comparison <- data.frame(
    Comparison = paste(pr[1, ], "-", pr[2, ]),
    Difference = as.vector(dif), SE_difference = as.vector(sed),
    Critical_difference = msd,
    Significant = abs(as.vector(dif)) > msd,
    p.value = as.vector(pval),
    row.names = NULL, stringsAsFactors = FALSE)

  structure(list(
    statistics = data.frame(MSerror = MSerror, Df = DFerror, Mean = mean(y),
                            CV = 100 * sqrt(MSerror) / mean(y),
                            Critical_difference = msd, row.names = NULL),
    parameters = data.frame(test = testname, ntr = ntr, r = r,
                            critical_value = critval, alpha = alpha,
                            error_stratum = stratum, scale = scale,
                            row.names = NULL, stringsAsFactors = FALSE),
    means = means, groups = groups, comparison = comparison),
    class = "cane_comparison")
}
