rtrunc_gamma <- function(n, lower = -Inf, upper = Inf, Tg_shape, Tg_rate) {
  pfun <- function(x) pgamma(x, shape = Tg_shape, rate = Tg_rate)
  qfun <- function(p) qgamma(p, shape = Tg_shape, rate = Tg_rate)
  a <- if (is.finite(lower)) pfun(lower) else 0
  b <- if (is.finite(upper)) pfun(upper) else 1
  if (b <= a) {
    ## No mass in the requested interval — return boundary value to avoid errors
    return(rep(if (is.finite(lower)) lower else upper, n))
  }
  u <- runif(n, a, b)
  qfun(u)
}

## Build a sampling grid for the upfront sanity check on time-varying
## probability parameters. For any inputs that came from make_time_varying(),
## we include the breakpoints (and midpoints) so the grid covers every
## changepoint. Falls back to an evenly spaced 0–365 grid when no
## time-varying inputs are supplied.
build_sanity_grid <- function(params) {
  breakpoints <- numeric(0)
  for (p in params) {
    if (inherits(p, "time_varying_fn")) {
      bp <- attr(p, "times", exact = TRUE)
      if (is.null(bp)) {
        ## Fallback for time_varying_fns built before make_time_varying
        ## exposed breakpoints as attributes.
        bp <- tryCatch(environment(p)$x, error = function(e) NULL)
      }
      if (is.numeric(bp)) {
        breakpoints <- c(breakpoints, bp)
      }
    }
  }
  breakpoints <- sort(unique(breakpoints))

  if (length(breakpoints) >= 2L) {
    mids <- (breakpoints[-length(breakpoints)] + breakpoints[-1]) / 2
    sort(unique(c(breakpoints, mids)))
  } else if (length(breakpoints) == 1L) {
    bp <- breakpoints[1]
    sort(unique(c(0, bp - 1, bp, bp + 1, max(bp + 10, 100))))
  } else {
    seq(0, 365, length.out = 50)
  }
}

## Sanity-check a probability parameter (scalar or function) by sampling it
## on `grid` and confirming every resolved value is in [0, 1]. NULL inputs
## are skipped so this can be called uniformly for optional parameters.
check_probability_on_grid <- function(param, grid, param_name) {
  if (is.null(param)) return(invisible(NULL))
  values <- resolve_time_varying(param, grid, param_name)
  if (any(values < 0 | values > 1)) {
    bad <- which(values < 0 | values > 1)
    show <- bad[seq_len(min(3L, length(bad)))]
    stop(sprintf(
      "`%s` resolves to values outside [0, 1] at t = %s (got %s).",
      param_name,
      paste(round(grid[show], 3), collapse = ", "),
      paste(round(values[show], 4), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

## Sanity-check a strictly positive parameter (scalar or function) by
## sampling it on `grid` and confirming every resolved value is > 0.
check_positive_on_grid <- function(param, grid, param_name) {
  if (is.null(param)) return(invisible(NULL))
  values <- resolve_time_varying(param, grid, param_name)
  if (any(values <= 0)) {
    bad <- which(values <= 0)
    show <- bad[seq_len(min(3L, length(bad)))]
    stop(sprintf(
      "`%s` must resolve to positive value(s); failed at t = %s (got %s).",
      param_name,
      paste(round(grid[show], 3), collapse = ", "),
      paste(round(values[show], 4), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

## Resolve a strictly-positive parameter (scalar or function(t)) at the supplied
## time(s), asserting every resolved value is finite and > 0. Used for
## transmissibility-style parameters (e.g. the mn_offspring_* means) that may be
## supplied either as a single positive scalar or as a function of absolute
## calendar time. The interpolation itself is delegated to resolve_time_varying()
## so callers stay scenario-agnostic.
resolve_positive_time_varying <- function(param, t, param_name = "parameter") {
  value <- resolve_time_varying(param = param, t = t, param_name = param_name)
  if (any(!is.finite(value)) || any(value <= 0)) {
    stop(sprintf("`%s` must resolve to positive value(s).", param_name), call. = FALSE)
  }
  value
}

## Sanity-check a non-negative parameter (scalar or function) by sampling it on
## `grid` and confirming every resolved value is finite and >= 0. Unlike
## check_positive_on_grid this allows the boundary value zero (relevant for
## obv_pep_dpc, where 0 means same-day treatment).
check_nonneg_on_grid <- function(param, grid, param_name) {
  if (is.null(param)) return(invisible(NULL))
  values <- resolve_time_varying(param, grid, param_name)
  if (any(!is.finite(values)) || any(values < 0)) {
    bad <- which(!is.finite(values) | values < 0)
    show <- bad[seq_len(min(3L, length(bad)))]
    stop(sprintf(
      "`%s` must resolve to finite, non-negative value(s); failed at t = %s (got %s).",
      param_name,
      paste(round(grid[show], 3), collapse = ", "),
      paste(round(values[show], 4), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}


#' Obeldesivir PEP efficacy as a function of days post challenge/exposure
#'
#' Returns OBV PEP efficacy for each DPC value. By default (\code{shape = "flat"})
#' efficacy is \emph{constant} at \code{E0} for every DPC -- a deliberately
#' model-neutral default, since the DPC-to-efficacy relationship is not yet
#' pinned down and the previous logistic constants were placeholders. Set
#' \code{shape = "logistic"} to use the NHP-derived scaled-logistic decay
#' (efficacy \code{E0} at DPC 0, falling to 0 by \code{dpc_zero} and forced to 0
#' beyond \code{max_dpc}); its default \code{d50}/\code{k}/\code{dpc_zero}/
#' \code{max_dpc} are the maximum-likelihood estimates from
#' \code{ODV_PEP_plot1_only_with_ribbon_GOF.R}. The curve is separable in
#' \code{E0}: \code{eff(dpc) = E0 * shape(dpc)} with \code{shape(0) = 1}, so
#' \code{E0} is a pure vertical scale equal to the peak (DPC 0) efficacy.
#'
#' @param dpc Numeric vector. Days post challenge/exposure at which OBV is first
#'   received.
#' @param E0 Efficacy at DPC 0 (the peak). For \code{shape = "flat"} this is the
#'   efficacy at \emph{every} DPC.
#' @param shape Character: \code{"flat"} (default; constant \code{E0}) or
#'   \code{"logistic"} (NHP-style decay using \code{d50}, \code{k},
#'   \code{dpc_zero}, \code{max_dpc}).
#' @param d50 Fitted logistic midpoint (used when \code{shape = "logistic"}).
#' @param k Fixed logistic steepness (used when \code{shape = "logistic"}).
#' @param dpc_zero DPC at which efficacy is forced to zero in the underlying fit
#'   (used when \code{shape = "logistic"}).
#' @param max_dpc Maximum DPC retained in this model branch; later dosing is
#'   assigned zero efficacy (used when \code{shape = "logistic"}).
#'
#' @return Numeric vector of efficacy values between 0 and 1.
#' @export
obv_pep_efficacy_from_dpc <- function(dpc,
                                      E0 = 0.82342697,
                                      shape = c("flat", "logistic"),
                                      d50 = 5.59722823,
                                      k = 1,
                                      dpc_zero = 15,
                                      max_dpc = 10) {
  if (!is.numeric(dpc)) {
    stop("`dpc` must be numeric.", call. = FALSE)
  }
  if (length(dpc) == 0L) {
    return(numeric(0))
  }
  if (any(!is.finite(dpc))) {
    stop("`dpc` must contain only finite values.", call. = FALSE)
  }
  if (any(dpc < 0)) {
    stop("`dpc` must be non-negative.", call. = FALSE)
  }
  shape <- match.arg(shape)

  if (!is.numeric(E0) || length(E0) != 1L || is.na(E0) || !is.finite(E0)) {
    stop("`E0` must be a single finite numeric value.", call. = FALSE)
  }
  if (E0 < 0 || E0 > 1) {
    stop("`E0` must be in [0, 1].", call. = FALSE)
  }

  ## Default model: flat efficacy -- E0 at every DPC, no decay. The NHP-style
  ## logistic decay is opt-in (shape = "logistic"); the normalised logistic
  ## cannot be made flat by any choice of d50/k/dpc_zero (as k -> 0 it tends to a
  ## LINEAR decay to 0 at dpc_zero, not a flat line), so flat is a distinct mode.
  if (shape == "flat") {
    return(rep(E0, length(dpc)))
  }

  for (nm in c("d50", "k", "dpc_zero", "max_dpc")) {
    val <- get(nm, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || !is.finite(val)) {
      stop(sprintf("`%s` must be a single finite numeric value.", nm), call. = FALSE)
    }
  }
  if (dpc_zero <= 0 || max_dpc < 0) {
    stop("`dpc_zero` must be positive and `max_dpc` must be non-negative.", call. = FALSE)
  }

  g_logistic <- function(d) 1 / (1 + exp(k * (d - d50)))

  g0 <- g_logistic(0)
  gz <- g_logistic(dpc_zero)
  gd <- g_logistic(dpc)
  denom <- g0 - gz
  if (!is.finite(denom) || abs(denom) < 1e-10) {
    stop("Invalid OBV efficacy curve: denominator is effectively zero.", call. = FALSE)
  }

  eff <- E0 * (gd - gz) / denom
  eff[dpc >= dpc_zero] <- 0
  eff[dpc > max_dpc] <- 0
  pmin(1, pmax(0, eff))
}

## OBV PEP accumulator. The first seven aggregate counters surface both
## pre-thinning (the "treat-all-contacts" policy denominator) and post-thinning
## (the "treat-only-PPE-failures" policy denominator) views of treatment, plus
## the count of infections actually prevented by adherent OBV. These seven are
## produced per offspring-function call by apply_obv_pep_gate() and summed.
##
## `prevented_deaths` is different in kind: it is NOT a gate counter. It is the
## subset of `prevented` infections that WOULD have died had they occurred,
## resolved once after the simulation loop by replaying them through the same
## outcome model as realised cases (see branching_process_main). It stays 0L
## here and is overwritten at the end of the run. Drawing it post-loop keeps the
## simulated trajectory's RNG stream untouched (see that deferral for why).
empty_obv_pep_num_treated <- function() {
  list(
    pre_eligible     = 0L,   # # candidates eligible BEFORE PPE/quarantine thinning
    pre_treated      = 0L,   # # of pre_eligible who received OBV (Bernoulli(coverage))
    pre_adherent     = 0L,   # # of pre_treated who adhered (Bernoulli(adherence))
    post_eligible    = 0L,   # # pre_eligible candidates that survived thinning (i.e. would-be cases)
    post_treated     = 0L,   # # pre_treated candidates that survived thinning
    post_adherent    = 0L,   # # pre_adherent candidates that survived thinning
    prevented        = 0L,   # # post_adherent whose infection was prevented by efficacy(dpc)
    prevented_deaths = 0L    # # of `prevented` infections that would have died (deferred counterfactual)
  )
}

## Empty (0-row) container for the per-call snapshot of infections the OBV gate
## prevented. Carries just enough to resolve each one's counterfactual death
## later: class (drives hospitalisation/CFR), infection location, and the
## absolute infection time (the calendar clock for time-varying parameters).
empty_obv_prevented_info <- function() {
  data.frame(
    class                   = character(0),
    infection_location      = character(0),
    time_infection_absolute = numeric(0),
    stringsAsFactors        = FALSE
  )
}

## Extract the candidate infections the OBV gate prevented (those that survived
## PPE/quarantine thinning but were blocked by the efficacy draw) from an
## offspring function's local state. `keep_mask` is apply_obv_pep_gate()'s
## `keep` vector, aligned to `which(keep_infection)`; its FALSE entries are the
## prevented kept-candidates. Reads already-drawn values only -- consumes no RNG
## -- so callers can stash the result and resolve would-be deaths after the loop.
extract_obv_prevented_info <- function(pre_thinning, keep_infection, keep_mask) {
  prevented_idx <- which(keep_infection)[!keep_mask]
  if (length(prevented_idx) == 0L) {
    return(empty_obv_prevented_info())
  }
  data.frame(
    class                   = pre_thinning$offspring_class[prevented_idx],
    infection_location      = pre_thinning$infection_location[prevented_idx],
    time_infection_absolute = pre_thinning$infection_time_absolute[prevented_idx],
    stringsAsFactors        = FALSE
  )
}

## Cache for the (constant) empty offspring data frame. Building a 0-row,
## 7-column data.frame per call is material in profiles: it is hit on every
## transmission route that yields no offspring (e.g. the funeral route for every
## parent who survives). The frame and its zero-valued obv_pep counter attribute
## never vary, and data frames are copy-on-modify, so handing back a shared
## cached object is safe -- callers only rbind it or read its attribute.
.empty_offspring_df_cache <- new.env(parent = emptyenv())

empty_offspring_dataframe <- function() {
  cached <- .empty_offspring_df_cache$value
  if (!is.null(cached)) {
    return(cached)
  }
  out <- data.frame(
    infection_location      = character(0),
    time_infection_relative = numeric(0),
    class                   = character(0),
    contact_risk_level      = integer(0),
    contact_risk_category   = character(0),
    traced                  = logical(0),
    obv_pep_eligible        = logical(0),
    obv_pep_received        = logical(0),
    obv_pep_adherent        = logical(0),
    obv_pep_dpc             = numeric(0),
    stringsAsFactors = FALSE
  )
  attr(out, "obv_pep_num_treated")    <- empty_obv_pep_num_treated()
  attr(out, "obv_pep_prevented_info") <- empty_obv_prevented_info()
  attr(out, "contact_log")            <- empty_contact_log()
  .empty_offspring_df_cache$value <- out
  out
}

## Validate the optional override list for the built-in OBV efficacy curve
## (obv_pep_efficacy_from_dpc). NULL or list() means "use the curve defaults".
## Names must be a subset of the curve's tunable arguments; `dpc` is supplied by
## the simulation and cannot be overridden. The overrides feed the built-in
## curve, which is selected when obv_pep_efficacy is either NULL (E0 from the
## args/default) or a numeric scalar (the scalar IS the curve's E0 -- its peak
## efficacy at DPC 0 -- with the shape from the args/default). Hence:
##   * a custom function obv_pep_efficacy takes no curve overrides (error); and
##   * when obv_pep_efficacy is a scalar it already supplies E0, so passing E0
##     again in the args is a double-specification (error).
validate_obv_efficacy_args <- function(obv_pep_efficacy, obv_pep_efficacy_args) {
  if (is.null(obv_pep_efficacy_args)) return(invisible(NULL))
  if (!is.list(obv_pep_efficacy_args)) {
    stop("`obv_pep_efficacy_args` must be NULL or a named list.", call. = FALSE)
  }
  if (length(obv_pep_efficacy_args) == 0L) return(invisible(NULL))
  nm <- names(obv_pep_efficacy_args)
  if (is.null(nm) || any(nm == "")) {
    stop("`obv_pep_efficacy_args` must be a fully named list.", call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop("`obv_pep_efficacy_args` has duplicate names.", call. = FALSE)
  }
  allowed <- setdiff(names(formals(obv_pep_efficacy_from_dpc)), "dpc")
  bad <- setdiff(nm, allowed)
  if (length(bad)) {
    stop(sprintf("`obv_pep_efficacy_args` has unknown argument(s): %s. Allowed: %s.",
                 paste(bad, collapse = ", "), paste(allowed, collapse = ", ")),
         call. = FALSE)
  }
  if (is.function(obv_pep_efficacy)) {
    stop("`obv_pep_efficacy_args` does not apply when `obv_pep_efficacy` is a function; parameterise that function directly.",
         call. = FALSE)
  }
  if (!is.null(obv_pep_efficacy) && "E0" %in% nm) {
    stop("`E0` is supplied via `obv_pep_efficacy` (used as the curve's E0); do not also pass it in `obv_pep_efficacy_args`.",
         call. = FALSE)
  }
  invisible(NULL)
}

resolve_obv_efficacy <- function(obv_pep_efficacy, dpc, obv_pep_efficacy_args = NULL) {
  if (length(dpc) == 0L) {
    return(numeric(0))
  }

  if (is.function(obv_pep_efficacy)) {
    value <- obv_pep_efficacy(dpc)
  } else {
    ## Built-in NHP-derived curve. A scalar obv_pep_efficacy supplies the curve's
    ## E0 (peak efficacy at DPC 0); the shape parameters (d50, k, dpc_zero,
    ## max_dpc) come from obv_pep_efficacy_args or the curve defaults. NULL
    ## obv_pep_efficacy takes E0 from the args/default too. obv_pep_efficacy_from_dpc
    ## validates the resulting parameters (e.g. E0 in [0, 1]) and errors on bad ones.
    curve_args <- obv_pep_efficacy_args
    if (!is.null(obv_pep_efficacy)) {
      curve_args <- c(list(E0 = obv_pep_efficacy), curve_args)
    }
    value <- if (length(curve_args)) {
      do.call(obv_pep_efficacy_from_dpc, c(list(dpc), curve_args))
    } else {
      obv_pep_efficacy_from_dpc(dpc)
    }
  }

  if (!is.numeric(value)) {
    stop("`obv_pep_efficacy` must be NULL, a function(dpc), or a numeric value in [0, 1].", call. = FALSE)
  }
  if (length(value) == 1L && length(dpc) > 1L) {
    value <- rep(value, length(dpc))
  }
  if (length(value) != length(dpc)) {
    stop("`obv_pep_efficacy` must return one value per DPC, or a single value to recycle.", call. = FALSE)
  }
  if (any(!is.finite(value)) || any(value < 0 | value > 1)) {
    stop("`obv_pep_efficacy` must resolve to value(s) in [0, 1].", call. = FALSE)
  }

  as.numeric(value)
}

#' Apply the obeldesivir PEP gate to candidate offspring (two-phase)
#'
#' This is an infection-prevention gate split into two phases that bracket the
#' offspring functions' Swiss-cheese PPE/ETU thinning step:
#'
#' \enumerate{
#'   \item \emph{Pre-thinning phase.} On the full set of candidate transmissions
#'     produced by an offspring function (before PPE / hospital quarantine has
#'     removed any), the gate identifies candidates whose offspring class is in
#'     \code{obv_pep_target_class} AND whose infection location is in
#'     \code{obv_pep_target_locations}. For each pre-thinning eligible
#'     candidate, the gate draws a treatment status (received OBV, adherent to
#'     the course) and a days-post-challenge value -- the latter either the
#'     deterministic \code{obv_pep_dpc} mean or, when \code{obv_pep_dpc_sd} is
#'     set, an independent Gamma draw with that mean (so individuals vary in how
#'     quickly they are dosed). The pre-thinning counters record the size of the
#'     "treat-all-contacts" cohort.
#'   \item \emph{Post-thinning phase.} The gate is told which pre-thinning
#'     candidates survived PPE/ETU thinning (\code{kept_indices}). Treatment
#'     status is carried through consistently: a candidate that was assigned
#'     received/adherent pre-thinning retains that status. The post-thinning
#'     counters record the size of the "treat-only-PPE-failures" cohort.
#'     Finally, for each \emph{kept and adherent} candidate, the gate draws
#'     Bernoulli(\code{efficacy(dpc)}) to decide whether the infection is
#'     prevented; prevented infections are removed from the output via the
#'     returned \code{keep} mask.
#' }
#'
#' Each candidate has a single set of treatment-status decisions; the pre- and
#' post-thinning counts are nested as sets per individual, so doses delivered
#' under Policy A ("OBV to all contacts") and Policy B ("OBV only when PPE
#' failed and infection would have occurred") are directly comparable from the
#' same simulation.
#'
#' Per-call multi-exposure caveat: the branching process tracks exposures, not
#' individuals. A single HCW exposed by two different parents during a
#' simulation contributes to the treatment counters twice. For exposures
#' separated by less than a typical OBV course length, this overstates real-world
#' dose requirements.
#'
#' @param pre_thinning Named list of equal-length vectors describing the
#'   pre-thinning candidate set: \code{infection_location}, \code{offspring_class},
#'   \code{infection_time_absolute}.
#' @param kept_indices Integer vector of positions into \code{pre_thinning}'s
#'   vectors indicating which candidates survived PPE/ETU thinning.
#' @param obv_pep_enabled Logical scalar. If FALSE, the gate is a no-op and all
#'   returned counters are zero.
#' @param obv_pep_coverage Numeric in \code{[0,1]} or function(t). Probability
#'   a pre-thinning eligible candidate receives OBV.
#' @param obv_pep_adherence Numeric in \code{[0,1]} or function(t). Probability a
#'   recipient adheres sufficiently for efficacy to apply.
#' @param obv_pep_dpc Non-negative numeric or function(t). Days post
#'   challenge / exposure to first dose. When \code{obv_pep_dpc_sd} is NULL
#'   this value is used directly as each recipient's DPC; when an sd is
#'   supplied it is the \emph{mean} of the per-recipient Gamma draw. Evaluated at
#'   the candidate's own absolute infection time, so the \emph{mean} delay may
#'   itself vary over calendar time while the efficacy(DPC) relationship stays
#'   fixed.
#' @param obv_pep_dpc_sd NULL or a single finite positive numeric. If NULL
#'   (default), DPC is deterministic and equal to \code{obv_pep_dpc} (identical
#'   RNG stream to pre-feature runs). If supplied, it is the standard deviation
#'   of the per-recipient DPC: each recipient's DPC is drawn independently from a
#'   Gamma reparameterised from (mean, sd) as
#'   \code{Gamma(shape = obv_pep_dpc(t)^2 / obv_pep_dpc_sd^2,
#'   scale = obv_pep_dpc_sd^2 / obv_pep_dpc(t))}, which has mean
#'   \code{obv_pep_dpc(t)} and variance \code{obv_pep_dpc_sd^2}
#'   (so smaller sd = tighter spread; \code{CV = obv_pep_dpc_sd / obv_pep_dpc(t)}).
#'   The spread is therefore a \emph{fixed absolute} sd even when the mean varies
#'   over calendar time (unlike a fixed shape, which would hold the CV constant).
#'   This models individual variation in how quickly the drug is received
#'   post-exposure. A mean of 0 yields a point mass at DPC 0.
#' @param obv_pep_efficacy NULL, numeric in \code{[0,1]}, or function(dpc).
#'   Selects the efficacy model. NULL or a scalar use the built-in
#'   \code{obv_pep_efficacy_from_dpc()} curve, which is \emph{flat} by default
#'   (constant efficacy at every DPC); a \emph{scalar} is taken as that constant
#'   (the curve's \code{E0}). DPC-dependent decay is opt-in via
#'   \code{obv_pep_efficacy_args} (\code{shape = "logistic"}). A function(dpc) is
#'   used as-is.
#' @param obv_pep_efficacy_args NULL or a named list of overrides for the
#'   built-in \code{obv_pep_efficacy_from_dpc()} curve -- any of \code{shape},
#'   \code{E0}, \code{d50}, \code{k}, \code{dpc_zero}, \code{max_dpc}. Applies to
#'   the built-in curve only (when \code{obv_pep_efficacy} is NULL or a scalar);
#'   e.g. \code{list(shape = "logistic", d50 = 4)} switches on the logistic decay.
#'   Errors if combined with a function \code{obv_pep_efficacy}, if it names
#'   \code{E0} while \code{obv_pep_efficacy} is a scalar (E0 then comes from
#'   \code{obv_pep_efficacy}), or given an unknown name.
#' @param obv_pep_target_class Character vector of offspring classes eligible
#'   for OBV PEP. Defaults to \code{"HCW"}.
#' @param obv_pep_target_locations Character vector of exposure locations
#'   eligible for OBV PEP. Defaults to \code{"hospital"} for HCW occupational
#'   exposures; set e.g. to \code{c("hospital", "community", "funeral")} to
#'   target HCWs in any setting.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{keep}}{Logical vector of length \code{length(kept_indices)};
#'       FALSE entries correspond to kept candidates whose infection was
#'       prevented by OBV efficacy. Use this mask on the offspring function's
#'       kept candidates to obtain the final realised offspring set.}
#'     \item{\code{metadata}}{Data frame with one row per element of
#'       \code{kept_indices}, with columns \code{obv_pep_eligible},
#'       \code{obv_pep_received}, \code{obv_pep_adherent}, and \code{obv_pep_dpc}.
#'       Non-recipients have \code{obv_pep_dpc = NA_real_} and the three booleans FALSE.
#'       Filter by \code{obv_pep_received == TRUE} before using \code{obv_pep_dpc}.}
#'     \item{\code{num_treated}}{Named list of integer counters
#'       (see \code{empty_obv_pep_num_treated()}). The gate populates the seven
#'       treatment/prevention counters; the \code{prevented_deaths} slot is left
#'       at 0L for the caller to fill in after the run.}
#'   }
#' @noRd
apply_obv_pep_gate <- function(pre_thinning,
                               kept_indices,
                               obv_pep_enabled = FALSE,
                               obv_pep_coverage = 0,
                               obv_pep_adherence = 1,
                               obv_pep_dpc = 1,
                               obv_pep_dpc_sd = NULL,
                               obv_pep_efficacy = NULL,
                               obv_pep_efficacy_args = NULL,
                               obv_pep_target_class = "HCW",
                               obv_pep_target_locations = "hospital") {

  if (!is.list(pre_thinning) ||
      !all(c("infection_location", "offspring_class", "infection_time_absolute")
           %in% names(pre_thinning))) {
    stop("`pre_thinning` must be a list with `infection_location`, `offspring_class`, and `infection_time_absolute` components.",
         call. = FALSE)
  }

  n_pre <- length(pre_thinning$infection_location)
  if (length(pre_thinning$offspring_class) != n_pre ||
      length(pre_thinning$infection_time_absolute) != n_pre) {
    stop("`pre_thinning` components must all have the same length.", call. = FALSE)
  }
  if (!is.numeric(kept_indices) ||
      (length(kept_indices) > 0L && (any(kept_indices < 1L) || any(kept_indices > n_pre)))) {
    stop("`kept_indices` must be integer indices into the pre-thinning candidate set.",
         call. = FALSE)
  }
  n_kept <- length(kept_indices)

  ## Defaults for non-recipients: booleans FALSE, dpc NA. NA disambiguates
  ## "did not receive OBV" from "received OBV at day 0" (same-day treatment).
  ## In any analysis, filter by `obv_pep_received == TRUE` before using dpc.
  metadata <- data.frame(
    obv_pep_eligible = rep(FALSE, n_kept),
    obv_pep_received = rep(FALSE, n_kept),
    obv_pep_adherent = rep(FALSE, n_kept),
    obv_pep_dpc      = rep(NA_real_, n_kept),
    stringsAsFactors = FALSE
  )
  keep <- rep(TRUE, n_kept)
  num_treated <- empty_obv_pep_num_treated()

  if (n_pre == 0L || !isTRUE(obv_pep_enabled)) {
    return(list(keep = keep, metadata = metadata, num_treated = num_treated))
  }

  if (!is.character(obv_pep_target_class) || length(obv_pep_target_class) < 1L) {
    stop("`obv_pep_target_class` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(obv_pep_target_locations) || length(obv_pep_target_locations) < 1L) {
    stop("`obv_pep_target_locations` must be a non-empty character vector.", call. = FALSE)
  }

  validate_probability_or_time_varying <- function(param, param_name) {
    if (is.function(param)) return(invisible(NULL))
    if (!is.numeric(param) || length(param) != 1L || is.na(param) || param < 0 || param > 1) {
      stop(sprintf("`%s` must be a function(t) or a single numeric in [0, 1].", param_name),
           call. = FALSE)
    }
    invisible(NULL)
  }

  validate_nonnegative_or_time_varying <- function(param, param_name) {
    if (is.function(param)) return(invisible(NULL))
    if (!is.numeric(param) || length(param) != 1L || is.na(param) || param < 0) {
      stop(sprintf("`%s` must be a function(t) or a single non-negative numeric.", param_name),
           call. = FALSE)
    }
    invisible(NULL)
  }

  resolve_probability <- function(param, t, param_name) {
    value <- resolve_time_varying(param = param, t = t, param_name = param_name)
    if (any(value < 0 | value > 1)) {
      stop(sprintf("`%s` must resolve to value(s) in [0, 1].", param_name), call. = FALSE)
    }
    value
  }

  resolve_nonnegative <- function(param, t, param_name) {
    value <- resolve_time_varying(param = param, t = t, param_name = param_name)
    if (any(!is.finite(value)) || any(value < 0)) {
      stop(sprintf("`%s` must resolve to finite, non-negative value(s).", param_name), call. = FALSE)
    }
    value
  }

  validate_probability_or_time_varying(obv_pep_coverage, "obv_pep_coverage")
  validate_probability_or_time_varying(obv_pep_adherence, "obv_pep_adherence")
  validate_nonnegative_or_time_varying(obv_pep_dpc, "obv_pep_dpc")
  if (!is.null(obv_pep_dpc_sd) &&
      (!is.numeric(obv_pep_dpc_sd) || length(obv_pep_dpc_sd) != 1L ||
       !is.finite(obv_pep_dpc_sd) || obv_pep_dpc_sd <= 0)) {
    stop("`obv_pep_dpc_sd` must be NULL or a single finite positive numeric.", call. = FALSE)
  }
  validate_obv_efficacy_args(obv_pep_efficacy, obv_pep_efficacy_args)

  ## --- Phase 1: pre-thinning eligibility, treatment status, DPC. ---
  ## Status vectors span the full pre-thinning set so they can be indexed by
  ## kept_indices in Phase 2 without bookkeeping gymnastics.
  ##
  ## Note: Phase 1 RNG draws (coverage, adherence, dpc) run for every
  ## pre-thinning eligible candidate, even ones whose `kept_indices` is empty
  ## (e.g. a safe funeral where all candidates were thinned upstream). This is
  ## by design -- the pre_* counters represent "Policy A doses delivered" and
  ## include candidates whose infections would have been blocked by PPE or
  ## safe-burial anyway. The minor compute overhead vs. an "only run draws if
  ## anything is kept" path is a deliberate trade for cleaner per-policy
  ## semantics. LOW PRIORITY: could short-circuit if `length(kept_indices) == 0`
  ## AND callers don't care about pre_* counters in that case.
  pre_eligible   <- pre_thinning$offspring_class %in% obv_pep_target_class &
                    pre_thinning$infection_location %in% obv_pep_target_locations
  num_treated$pre_eligible <- sum(pre_eligible)

  status_received <- rep(FALSE, n_pre)
  status_adherent <- rep(FALSE, n_pre)
  status_dpc      <- rep(NA_real_, n_pre)

  if (any(pre_eligible)) {
    pre_eligible_idx <- which(pre_eligible)
    coverage_t <- resolve_probability(
      obv_pep_coverage,
      pre_thinning$infection_time_absolute[pre_eligible_idx],
      "obv_pep_coverage"
    )
    pre_received <- as.logical(rbinom(n = length(pre_eligible_idx), size = 1, prob = coverage_t))
    status_received[pre_eligible_idx] <- pre_received
    num_treated$pre_treated <- sum(pre_received)

    pre_received_idx <- pre_eligible_idx[pre_received]
    if (length(pre_received_idx) > 0L) {
      adherence_t <- resolve_probability(
        obv_pep_adherence,
        pre_thinning$infection_time_absolute[pre_received_idx],
        "obv_pep_adherence"
      )
      pre_adherent <- as.logical(rbinom(n = length(pre_received_idx), size = 1, prob = adherence_t))
      status_adherent[pre_received_idx] <- pre_adherent
      num_treated$pre_adherent <- sum(pre_adherent)

      ## obv_pep_dpc(t) gives the MEAN days-post-challenge for a candidate
      ## infected at calendar time t. With obv_pep_dpc_sd = NULL the DPC is
      ## that mean exactly (deterministic; same RNG stream as pre-feature runs).
      ## Otherwise each recipient draws its own DPC from a Gamma with that mean
      ## and the supplied fixed standard deviation. The Gamma is reparameterised
      ## from (mean, sd) to (shape, scale) via shape = mean^2 / sd^2 and
      ## scale = sd^2 / mean, so the draw has mean `mean` and a constant absolute
      ## variance sd^2 (i.e. CV = sd / mean) -- individuals vary in how quickly
      ## they receive the drug post-exposure. This per-individual draw IS the
      ## realised (receipt - infection) delay. A mean of 0 has no positive-sd
      ## Gamma (a non-negative variable with mean 0 is 0 a.s.), so it degenerates
      ## to a point mass at DPC 0; guard it to avoid the shape = 0 / scale = Inf
      ## edge in rgamma() and keep same-day dosing exact.
      dpc_mean <- resolve_nonnegative(
        obv_pep_dpc,
        pre_thinning$infection_time_absolute[pre_received_idx],
        "obv_pep_dpc"
      )
      if (is.null(obv_pep_dpc_sd)) {
        status_dpc[pre_received_idx] <- dpc_mean
      } else {
        drawn_dpc <- numeric(length(dpc_mean))
        positive  <- dpc_mean > 0
        if (any(positive)) {
          drawn_dpc[positive] <- rgamma(
            n     = sum(positive),
            shape = dpc_mean[positive]^2 / obv_pep_dpc_sd^2,
            scale = obv_pep_dpc_sd^2 / dpc_mean[positive]
          )
        }
        status_dpc[pre_received_idx] <- drawn_dpc
      }
    }
  }

  if (n_kept == 0L) {
    return(list(keep = keep, metadata = metadata, num_treated = num_treated))
  }

  ## --- Phase 2: intersect with kept; apply efficacy to (kept & adherent). ---
  kept_eligible <- pre_eligible[kept_indices]
  kept_received <- status_received[kept_indices]
  kept_adherent <- status_adherent[kept_indices]
  kept_dpc      <- status_dpc[kept_indices]

  num_treated$post_eligible <- sum(kept_eligible)
  num_treated$post_treated  <- sum(kept_received)
  num_treated$post_adherent <- sum(kept_adherent)

  if (any(kept_adherent)) {
    adh_local_idx <- which(kept_adherent)
    efficacy_vals <- resolve_obv_efficacy(obv_pep_efficacy, kept_dpc[adh_local_idx],
                                          obv_pep_efficacy_args = obv_pep_efficacy_args)
    prevented <- as.logical(rbinom(n = length(adh_local_idx), size = 1, prob = efficacy_vals))
    keep[adh_local_idx[prevented]] <- FALSE
    num_treated$prevented <- sum(prevented)
  }

  metadata$obv_pep_eligible <- kept_eligible
  metadata$obv_pep_received <- kept_received
  metadata$obv_pep_adherent <- kept_adherent
  metadata$obv_pep_dpc      <- kept_dpc

  list(keep = keep, metadata = metadata, num_treated = num_treated)
}
