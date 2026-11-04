## Shared engine for the five design-specific functions ---------------------

#' @keywords internal
#' @noRd
.cane_engine <- function(data, Response, Location, alpha, Mult_Comp_Test,
                         transform, env_effect,
                         within_terms, within_spec,
                         pooled_terms, pooled_spec,
                         compare, gei_label, design) {

  dat <- data
  L   <- Location
  envs <- levels(dat[[L]])

  ## ---- 1. environment-wise analysis ------------------------------------
  ind <- vector("list", length(envs)); names(ind) <- envs
  res_list <- vector("list", length(envs)); names(res_list) <- envs
  ems <- stats::setNames(numeric(length(envs)), envs)
  dfe <- stats::setNames(integer(length(envs)), envs)
  strata_tab <- list()

  for (e in envs) {
    sub <- droplevels(dat[dat[[L]] == e, , drop = FALSE])
    fit <- .cane_fit(Response, within_terms, sub)
    tab <- .cane_strata_anova(fit, within_spec, residual_label = design$within_resid)
    ind[[e]] <- tab
    res_list[[e]] <- stats::residuals(fit)
    ems[e] <- stats::deviance(fit) / stats::df.residual(fit)
    dfe[e] <- stats::df.residual(fit)
    strata_tab[[e]] <- tab
  }
  if (any(dfe < 1L))
    .cane_stop("Environment(s) ", paste(envs[dfe < 1L], collapse = ", "),
               " leave no residual degrees of freedom. The design is saturated there, so no ",
               "error mean square can be estimated and no pooled analysis is possible.")

  ## ---- 2. homogeneity of error variances -------------------------------
  bart <- .cane_bartlett(ems, dfe, alpha)
  ## designs with more than one error stratum: test the other strata as well
  bart$other_strata <- .cane_bartlett_strata(ind, design$other_strata, alpha)
  lev  <- .cane_levene(res_list, dfe)
  norm <- .cane_normality(res_list)
  dec  <- .cane_decide_transform(bart, lev, transform, alpha)
  if (length(design$other_strata)) {
    dec$reason <- paste0(dec$reason, " The Bartlett test, and the Aitken scaling when applied, use the ",
                         tolower(design$within_resid), " mean square of each environment.")
    os <- bart$other_strata
    het <- os$Stratum[!is.na(os$heterogeneous) & os$heterogeneous]
    if (length(het))
      dec$reason <- paste0(dec$reason, " Bartlett's test also rejects homogeneity of the ",
                           paste(tolower(het), collapse = " and "),
                           " across environments (see Bartlett_Test$other_strata); tests ",
                           "in that stratum pool heterogeneous variances and should be ",
                           "interpreted with caution, or the series analysed with a mixed ",
                           "model that allows a separate variance for each environment.")
  }

  ## ---- 3. conditional Aitken transformation ----------------------------
  ## Each observation is divided by the square root of the error mean square
  ## of ITS OWN environment. The lookup is by environment label, which is
  ## order-independent and cannot be misaligned by the row order of 'data'.
  resp_used <- Response
  if (dec$apply) {
    dat[[".cane_transformed"]] <- dat[[Response]] / sqrt(ems[as.character(dat[[L]])])
    resp_used <- ".cane_transformed"
  }

  ## ---- 4. pooled ANOVA with the correct error stratum ------------------
  pfit <- .cane_fit(resp_used, pooled_terms, dat)
  if (identical(env_effect, "random")) {
    pooled_spec <- design$randomise(pooled_spec)
    ## With random environments the pooled comparison of means must use the
    ## same denominator as the F-test of the factor compared (for example the
    ## treatment-by-environment mean square), not the pooled error.
    compare <- lapply(compare, function(cmp) {
      if (!is.null(cmp$spec) && !is.null(pooled_spec[[cmp$spec]]))
        cmp$pooled_error <- pooled_spec[[cmp$spec]]$denom
      cmp
    })
  }
  ptab <- .cane_strata_anova(pfit, pooled_spec, residual_label = design$pooled_resid)

  ## diagnostic identity: on the Aitken scale the pooled error mean square
  ## of a single-stratum design is 1 by construction. Report it as a check.
  check <- NULL
  if (dec$apply) {
    i <- which(ptab$Source == design$pooled_resid)
    check <- list(pooled_error_MS = ptab$`Mean Sq`[i],
                  expected = 1,
                  note = paste0("On the Aitken scale the pooled error mean square of a ",
                                "single-error design equals 1 by construction, because each ",
                                "environment's error sum of squares is divided by its own ",
                                "error mean square. Observed: ",
                                signif(ptab$`Mean Sq`[i], 6), "."))
  }

  ## ---- 5. multiple comparisons -----------------------------------------
  .cane_lsd_warning(max(vapply(compare, function(z) nlevels(droplevels(factor(dat[[z$factor]]))), integer(1))), Mult_Comp_Test)

  loc_wise <- vector("list", length(envs)); names(loc_wise) <- envs
  for (e in envs) {
    sub <- droplevels(dat[dat[[L]] == e, , drop = FALSE])
    out_e <- list()
    for (cmp in compare) {
      st <- .cane_stratum(ind[[e]], cmp$within_error)
      if (is.null(st) || is.na(st$MS) || st$Df < 1L) next
      out_e[[cmp$name]] <- tryCatch(
        .cane_mean_sep(sub[[Response]], sub[[cmp$factor]], st$MS, st$Df, alpha,
                       Mult_Comp_Test, stratum = st$label, scale = "original"),
        error = function(ee) { warning("Environment ", e, ", ", cmp$name, ": ",
                                       conditionMessage(ee), call. = FALSE); NULL })
    }
    loc_wise[[e]] <- if (length(out_e) == 1L) out_e[[1L]] else out_e
  }

  pooled_cmp <- list()
  for (cmp in compare) {
    st <- .cane_stratum(ptab, cmp$pooled_error)
    if (is.null(st) || is.na(st$MS) || st$Df < 1L) {
      warning("No pooled comparison for ", cmp$name,
              ": the ", paste(cmp$pooled_error, collapse = ":"),
              " stratum has no estimable error. This design needs more replication.",
              call. = FALSE)
      next
    }
    bf <- tapply(dat[[Response]], droplevels(factor(dat[[cmp$factor]])), mean)
    res_cmp <- tryCatch(
      .cane_mean_sep(
        dat[[resp_used]], dat[[cmp$factor]], st$MS, st$Df, alpha, Mult_Comp_Test,
        stratum = st$label,
        scale = if (dec$apply) "Aitken-transformed" else "original",
        back_factor = if (dec$apply) bf else NULL),
      error = function(ee) { warning("No pooled comparison for ", cmp$name, ": ",
                                     conditionMessage(ee), call. = FALSE); NULL })
    if (!is.null(res_cmp)) pooled_cmp[[cmp$name]] <- res_cmp
  }

  out <- list(
    Design = design$name,
    Call = NULL,
    Individual_ANOVA = ind,
    Location_wise = loc_wise,
    Bartlett_Test = bart,
    Levene_Test = lev,
    Normality_Test = norm,
    Transformation = list(applied = dec$apply, reason = dec$reason,
                          error_mean_squares = ems, error_df = dfe,
                          variance_ratio = bart$variance_ratio,
                          scale = if (dec$apply) "Aitken-transformed" else "original",
                          diagnostic = check),
    Pooled_ANOVA = ptab,
    Treatments_Comparison = if (length(pooled_cmp) == 1L) pooled_cmp[[1L]] else pooled_cmp,
    GEI_Note = .cane_gei_note(ptab, gei_label, alpha),
    Environment_Key = data.frame(Level = envs,
                                 n = as.vector(table(dat[[L]])[envs]),
                                 Error_MS = as.vector(ems),
                                 Error_df = as.vector(dfe),
                                 row.names = NULL),
    Settings = list(Response = Response, Location = Location, alpha = alpha,
                    Mult_Comp_Test = Mult_Comp_Test, transform = transform,
                    env_effect = env_effect,
                    CANE_version = as.character(utils::packageVersion("CANE")),
                    R_version = paste(R.version$major, R.version$minor, sep = ".")))
  class(out) <- c("cane", "list")
  out
}
