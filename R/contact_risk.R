## Contact-first transmission: risk-tier structure and its runtime helpers.
##
## fiber generates transmission contact-first. The Negative Binomial draw in each
## offspring function is the number of CONTACTS (exposure events) a case has.
## Every contact is assigned a risk tier; a contact in tier l transmits with
## probability
##
##     baseline_risk(t) * relative_risk[l]
##
## resolved at that contact's own calendar time. Contacts that transmit then pass
## through the intervention layers (PPE, hospital quarantine, safe burial)
## exactly as before.
##
## By convention the REFERENCE tier -- the one `baseline_risk` describes -- is the
## highest-risk tier, so relative risks sit in (0, 1] and `baseline_risk` is the
## transmission probability of the riskiest contact (a household or caregiving
## exposure, the thing secondary attack rates measure). That keeps the feasibility
## constraint `baseline_risk * max_relative_risk <= 1` from ever binding, and gives
## a calibration a rectangular [0, 1] x [0, 1] parameter space instead of one whose
## bounds move as the spread is fitted. Anchoring on a different tier is allowed and
## changes nothing -- only the product baseline_risk * relative_risk[l] is ever used.
##
## Because each contact's tier is drawn independently, a contact's marginal
## transmission probability is `baseline_risk * mean_relative_risk` whatever tier
## it lands in. Two consequences worth knowing:
##
##   * The expected offspring count for a route is
##         mn_contacts * baseline_risk * mean_relative_risk
##     before intervention thinning, so the contact overdispersion never enters
##     R0 (thinning a Negative Binomial leaves its mean unchanged) and can be
##     used to dial superspreading independently of the calibration.
##
##   * The tier mix among realised INFECTIONS is not the tier mix among contacts:
##     higher-risk tiers are over-represented, with weights proportional to
##     fractions * relative_risk. `case_weights` below carries that distribution,
##     and it is what the R0 approximation must average over when a downstream
##     process (contact tracing) depends on the tier a case was infected in.
##
## Risk tiers also carry a per-tier contact-tracing probability, which is the
## route by which they drive the non-pharmaceutical interventions: a traced case
## is admitted sooner, and optionally more often.

#' Define a contact risk-tier structure
#'
#' Builds the risk-tier structure used by fiber's contact-first transmission
#' model. Contacts are partitioned into tiers by \code{fractions}; each tier
#' carries a \code{relative_risk} multiplier applied to a single baseline
#' per-contact transmission probability, and a \code{trace_prob} giving the
#' probability a contact in that tier is reached by contact tracing at full
#' programme coverage.
#'
#' The tier named by \code{reference} is the one the baseline risk refers to: all
#' relative risks are divided through by that tier's value, so the baseline risk
#' parameter is literally the per-contact transmission probability of the
#' reference tier.
#'
#' \strong{By default the reference is the highest-risk tier}, so relative risks
#' come out in \code{(0, 1]} and the baseline risk is the transmission probability
#' of the riskiest contact -- a household or caregiving exposure, say, which is the
#' quantity secondary attack rates actually measure. Three things follow:
#' \itemize{
#'   \item The feasibility constraint \code{baseline_risk * max_relative_risk <= 1}
#'         reduces to \code{baseline_risk <= 1}, so it can never bind. Anchoring on
#'         the lowest tier instead leaves an upper bound on the baseline risk that
#'         moves with the relative-risk spread.
#'   \item Both \code{baseline_risk} and the relative risks live in \code{[0, 1]}
#'         independently, which gives a calibration a rectangular parameter space
#'         rather than one whose bounds shift as the spread is fitted.
#'   \item \code{mean_relative_risk} is then \code{<= 1} and attenuates rather than
#'         amplifies, so R0 sweeps from 0 to \code{mn_contacts * mean_relative_risk}
#'         as the baseline risk sweeps from 0 to 1.
#' }
#' Naming a different \code{reference} is allowed and changes nothing about the
#' model -- it is a pure reparameterisation, since only the product
#' \code{baseline_risk * relative_risk[l]} enters anywhere -- but it does bring the
#' moving feasibility bound back.
#'
#' The default structure is five equally sized tiers with equal relative risk and
#' no tracing -- deliberately flat, adding the contact bookkeeping without imposing
#' any heterogeneity or intervention effect. Use
#' \code{\link{contact_risk_gradient}} to build a graded structure from a single
#' spread parameter.
#'
#' @param fractions Numeric vector of tier shares. Must be non-negative and sum
#'   to 1. Its length sets the number of tiers. Defaults to five equal shares.
#' @param relative_risk Numeric vector, same length as \code{fractions}. Risk of
#'   each tier relative to the reference tier. Must be positive and finite. Only
#'   the ratios matter -- the vector is rescaled so the reference tier is exactly
#'   1 -- so \code{c(1, 2, 5)} and \code{c(0.2, 0.4, 1)} describe the same structure.
#' @param trace_prob Numeric vector in \code{[0,1]}, same length as
#'   \code{fractions} (or a single value recycled across tiers). Probability that
#'   a contact in each tier is successfully traced, at full programme coverage.
#'   Multiplied by the time-varying \code{trace_coverage(t)} programme lever at
#'   run time. Defaults to 0 (no tracing).
#' @param reference Integer index (1-based) of the reference tier -- the tier the
#'   baseline risk parameter describes. Defaults to the \emph{highest-risk} tier,
#'   which keeps every relative risk in \code{(0, 1]} and the baseline risk
#'   unconstrained; see the details above.
#' @param labels Optional character vector of tier names, same length as
#'   \code{fractions}. Used in the contact log and summaries. Defaults to
#'   \code{"risk_1"}, \code{"risk_2"}, ...
#'
#' @return An object of class \code{"fiber_contact_risk"}: a list with the
#'   validated \code{fractions}, normalised \code{relative_risk},
#'   \code{trace_prob}, \code{reference}, \code{labels}, \code{n_levels}, and the
#'   derived quantities \code{mean_relative_risk} (the fraction-weighted mean,
#'   which scales R0), \code{max_relative_risk} (which caps the feasible baseline
#'   risk at \code{1 / max_relative_risk}), \code{case_weights} (the tier
#'   distribution among realised infections, proportional to
#'   \code{fractions * relative_risk}) and \code{mean_trace_prob_cases} (the
#'   case-weighted mean trace probability, which is what the R0 approximation
#'   needs).
#'
#' @examples
#' ## Default: five flat tiers, no tracing
#' make_contact_risk()
#'
#' ## Most contacts low risk, a small high-risk tail that is easier to trace.
#' ## Written as risks relative to the lowest tier, but rescaled on the way in so
#' ## the top tier is 1 -- the baseline risk is then that tier's own probability.
#' make_contact_risk(fractions     = c(0.6, 0.25, 0.1, 0.04, 0.01),
#'                   relative_risk = c(1, 2, 5, 10, 25),
#'                   trace_prob    = c(0.1, 0.2, 0.4, 0.7, 0.9))
#' @export
make_contact_risk <- function(fractions     = rep(0.2, 5),
                              relative_risk = rep(1, 5),
                              trace_prob    = 0,
                              reference     = which.max(relative_risk),
                              labels        = NULL) {

  if (!is.numeric(fractions) || length(fractions) < 1L ||
      any(!is.finite(fractions)) || any(fractions < 0)) {
    stop("`fractions` must be a non-empty numeric vector of finite, non-negative values.",
         call. = FALSE)
  }
  if (abs(sum(fractions) - 1) > 1e-8) {
    stop(sprintf("`fractions` must sum to 1 (got %s).", format(sum(fractions))),
         call. = FALSE)
  }

  n_levels <- length(fractions)

  if (!is.numeric(relative_risk) || length(relative_risk) != n_levels ||
      any(!is.finite(relative_risk)) || any(relative_risk <= 0)) {
    stop(sprintf("`relative_risk` must be a numeric vector of %d finite, strictly positive values (same length as `fractions`).",
                 n_levels), call. = FALSE)
  }

  if (length(trace_prob) == 1L) trace_prob <- rep(trace_prob, n_levels)
  if (!is.numeric(trace_prob) || length(trace_prob) != n_levels ||
      any(!is.finite(trace_prob)) || any(trace_prob < 0) || any(trace_prob > 1)) {
    stop(sprintf("`trace_prob` must be a single value or a numeric vector of %d values in [0, 1].",
                 n_levels), call. = FALSE)
  }

  if (!is.numeric(reference) || length(reference) != 1L || is.na(reference) ||
      reference != round(reference) || reference < 1 || reference > n_levels) {
    stop(sprintf("`reference` must be a single integer index between 1 and %d.", n_levels),
         call. = FALSE)
  }
  reference <- as.integer(reference)

  ## Normalise so the reference tier has relative risk exactly 1. The baseline
  ## risk parameter is then that tier's own per-contact transmission probability,
  ## whether or not the user happened to write its relative risk as 1. Only the
  ## product baseline_risk * relative_risk[l] is ever used, so this rescaling is a
  ## pure reparameterisation -- it cannot change a simulated outcome.
  relative_risk <- relative_risk / relative_risk[reference]

  if (is.null(labels)) {
    labels <- paste0("risk_", seq_len(n_levels))
  }
  if (!is.character(labels) || length(labels) != n_levels || anyNA(labels)) {
    stop(sprintf("`labels` must be NULL or a character vector of %d non-missing names.",
                 n_levels), call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop("`labels` must be unique.", call. = FALSE)
  }

  mean_relative_risk <- sum(fractions * relative_risk)
  ## Tier distribution among realised infections. A contact in a high-risk tier
  ## is more likely to become a case, so cases over-represent those tiers with
  ## weights proportional to fractions * relative_risk.
  case_weights <- fractions * relative_risk / mean_relative_risk

  structure(
    list(
      fractions             = as.numeric(fractions),
      relative_risk         = as.numeric(relative_risk),
      trace_prob            = as.numeric(trace_prob),
      reference             = reference,
      labels                = labels,
      n_levels              = n_levels,
      mean_relative_risk    = mean_relative_risk,
      max_relative_risk     = max(relative_risk),
      case_weights          = case_weights,
      mean_trace_prob_cases = sum(case_weights * trace_prob)
    ),
    class = "fiber_contact_risk"
  )
}

#' @export
print.fiber_contact_risk <- function(x, ...) {
  cat("<fiber contact risk structure>\n")
  cat(sprintf("  %d tiers, reference tier: %s\n", x$n_levels, x$labels[x$reference]))
  print(data.frame(
    level         = seq_len(x$n_levels),
    label         = x$labels,
    fraction      = x$fractions,
    relative_risk = x$relative_risk,
    trace_prob    = x$trace_prob,
    case_share    = x$case_weights,
    stringsAsFactors = FALSE
  ), row.names = FALSE)
  cat(sprintf("  mean relative risk: %.4f  (R0 = mean contacts x baseline risk x this)\n",
              x$mean_relative_risk))
  cat(sprintf("  max relative risk:  %.4f  (baseline risk must be <= %.4f)\n",
              x$max_relative_risk, 1 / x$max_relative_risk))
  cat(sprintf("  case-weighted trace probability: %.4f\n", x$mean_trace_prob_cases))
  invisible(x)
}

#' Build a graded contact risk structure from a single spread parameter
#'
#' Convenience constructor for a log-spaced gradient of relative risks, so the
#' amount of contact heterogeneity can be swept with one number instead of
#' hand-writing a vector of weights. Trace probabilities can be given the same
#' treatment via \code{trace_prob_range}, since higher-risk contacts (household,
#' caregiving) are typically both more infectious and easier to identify.
#'
#' \code{ratio} is the spread between the lowest and highest tier. Following the
#' package convention (see \code{\link{make_contact_risk}}), the \emph{highest}
#' tier is the reference, so the relative risks run from \code{1/ratio} up to 1
#' and the baseline risk describes the riskiest contact.
#'
#' @param n_levels Integer, number of tiers. Defaults to 5.
#' @param ratio Positive numeric. Ratio of the top tier's risk to the bottom
#'   tier's. \code{ratio = 1} reproduces a flat structure.
#' @param trace_prob_range Either a single probability applied to every tier, or
#'   a length-2 vector \code{c(lowest, highest)} linearly interpolated across
#'   tiers. Defaults to 0 (no tracing).
#' @param fractions Optional numeric vector of tier shares (must sum to 1).
#'   Defaults to equal shares.
#' @param reference Integer index of the reference tier. Defaults to the highest
#'   tier, so relative risks come out in \code{(0, 1]}.
#' @param labels Optional character vector of tier names.
#'
#' @return An object of class \code{"fiber_contact_risk"}.
#'
#' @examples
#' ## Relative risks 0.1, 0.18, 0.32, 0.56, 1 -- a ten-fold spread anchored on the top tier
#' contact_risk_gradient(n_levels = 5, ratio = 10,
#'                       trace_prob_range = c(0.1, 0.8))
#' @export
contact_risk_gradient <- function(n_levels         = 5,
                                  ratio            = 10,
                                  trace_prob_range = 0,
                                  fractions        = NULL,
                                  reference        = n_levels,
                                  labels           = NULL) {

  if (!is.numeric(n_levels) || length(n_levels) != 1L || is.na(n_levels) ||
      n_levels != round(n_levels) || n_levels < 1) {
    stop("`n_levels` must be a single positive integer.", call. = FALSE)
  }
  n_levels <- as.integer(n_levels)

  if (!is.numeric(ratio) || length(ratio) != 1L || !is.finite(ratio) || ratio <= 0) {
    stop("`ratio` must be a single finite positive numeric value.", call. = FALSE)
  }

  if (is.null(fractions)) {
    fractions <- rep(1 / n_levels, n_levels)
  }

  relative_risk <- if (n_levels == 1L) 1 else exp(seq(0, log(ratio), length.out = n_levels))

  if (!is.numeric(trace_prob_range) || !length(trace_prob_range) %in% c(1L, 2L)) {
    stop("`trace_prob_range` must be a single probability or a length-2 vector c(lowest, highest).",
         call. = FALSE)
  }
  trace_prob <- if (length(trace_prob_range) == 1L || n_levels == 1L) {
    rep(trace_prob_range[1], n_levels)
  } else {
    seq(trace_prob_range[1], trace_prob_range[2], length.out = n_levels)
  }

  make_contact_risk(fractions     = fractions,
                    relative_risk = relative_risk,
                    trace_prob    = trace_prob,
                    reference     = reference,
                    labels        = labels)
}

## Internal: accept either a fiber_contact_risk object or a plain list of
## make_contact_risk() arguments, and return a validated structure. NULL falls
## back to `default`, which lets each transmission route inherit a shared
## structure unless it overrides it.
as_contact_risk <- function(x, default = NULL, param_name = "contact_risk") {
  if (is.null(x)) {
    if (is.null(default)) {
      return(make_contact_risk())
    }
    return(as_contact_risk(default, NULL, param_name))
  }
  if (inherits(x, "fiber_contact_risk")) {
    return(x)
  }
  if (is.list(x)) {
    return(do.call(make_contact_risk, x))
  }
  stop(sprintf("`%s` must be NULL, a `fiber_contact_risk` object (see make_contact_risk()), or a named list of its arguments.",
               param_name), call. = FALSE)
}

## ---------------------------------------------------------------------------
## Runtime pieces used by the offspring functions
## ---------------------------------------------------------------------------

## Assign a risk tier to each of `n` contacts, drawn independently from the
## structure's tier fractions. Consumes exactly one RNG call.
draw_contact_risk_levels <- function(n, contact_risk) {
  if (n == 0L) return(integer(0))
  sample.int(contact_risk$n_levels, size = n, replace = TRUE,
             prob = contact_risk$fractions)
}

## Per-contact transmission probability = baseline_risk(t) * relative_risk[tier].
## `baseline_risk` may be a scalar or a function of absolute calendar time, and
## is resolved at each contact's own contact time (matching how PPE coverage and
## quarantine efficacy are resolved per event rather than per parent).
contact_transmission_prob <- function(levels,
                                      contact_times_absolute,
                                      baseline_risk,
                                      contact_risk,
                                      param_name = "baseline_risk") {
  if (length(levels) == 0L) return(numeric(0))
  p0 <- resolve_time_varying(baseline_risk, contact_times_absolute, param_name)
  if (any(p0 < 0 | p0 > 1)) {
    stop(sprintf("`%s` must resolve to value(s) in [0, 1].", param_name), call. = FALSE)
  }
  p <- p0 * contact_risk$relative_risk[levels]
  if (any(p > 1 + 1e-12)) {
    i <- which(p > 1 + 1e-12)[1]
    stop(sprintf(
      "`%s` = %s at t = %s gives tier '%s' (relative risk %s) a transmission probability of %s, which exceeds 1. Lower the baseline risk or narrow the relative-risk spread.",
      param_name,
      format(round(p0[i], 6)),
      format(round(contact_times_absolute[i], 3)),
      contact_risk$labels[levels[i]],
      format(round(contact_risk$relative_risk[levels[i]], 4)),
      format(round(p[i], 4))
    ), call. = FALSE)
  }
  pmin(p, 1)
}

## Per-contact tracing probability = trace_coverage(t) * trace_prob[tier],
## capped at 1. Coverage is the time-varying programme lever (how much tracing
## capacity exists); trace_prob is the fixed per-tier traceability (how findable
## that kind of contact is). Resolved at each contact's own calendar time.
contact_trace_prob <- function(levels,
                               contact_times_absolute,
                               trace_coverage,
                               contact_risk,
                               param_name = "trace_coverage") {
  if (length(levels) == 0L) return(numeric(0))
  cov_t <- resolve_time_varying(trace_coverage, contact_times_absolute, param_name)
  if (any(cov_t < 0 | cov_t > 1)) {
    stop(sprintf("`%s` must resolve to value(s) in [0, 1].", param_name), call. = FALSE)
  }
  pmin(cov_t * contact_risk$trace_prob[levels], 1)
}

## Empty (0-row) contact log, matching the columns built by new_contact_log().
empty_contact_log <- function() {
  data.frame(
    parent                = integer(0),
    case_id               = integer(0),
    record_type           = character(0),
    class                 = character(0),
    infection_location    = character(0),
    time_contact_relative = numeric(0),
    time_contact_absolute = numeric(0),
    contact_risk_level    = integer(0),
    contact_risk_category = character(0),
    relative_risk         = numeric(0),
    transmission_prob     = numeric(0),
    traced                = logical(0),
    blocked_by            = character(0),
    stringsAsFactors      = FALSE
  )
}

## Build the per-parent contact log. One row per contact generated, in the order
## the contacts were drawn.
##
##   record_type  "infection" if the contact ended up as a realised case,
##                "contact" otherwise.
##   blocked_by   NA for realised infections; otherwise the layer that stopped
##                it -- "no_transmission" (the risk draw), the route's
##                intervention layer, or "obv_pep".
##
## `case_id` is left NA here because ids are assigned by the simulation loop;
## branching_process_main fills it in for the "infection" rows, which appear in
## the same order as the offspring rows this call returns.
new_contact_log <- function(parent_id,
                            offspring_class,
                            infection_location,
                            time_contact_relative,
                            time_contact_absolute,
                            risk_levels,
                            transmission_prob,
                            traced,
                            contact_risk,
                            transmitted,
                            kept,
                            realised,
                            intervention_label) {

  n <- length(offspring_class)
  if (n == 0L) return(empty_contact_log())

  blocked_by <- rep(NA_character_, n)
  blocked_by[!transmitted]        <- "no_transmission"
  blocked_by[transmitted & !kept] <- intervention_label
  blocked_by[kept & !realised]    <- "obv_pep"

  data.frame(
    parent                = rep(parent_id, n),
    case_id               = rep(NA_integer_, n),
    record_type           = ifelse(realised, "infection", "contact"),
    class                 = offspring_class,
    infection_location    = infection_location,
    time_contact_relative = time_contact_relative,
    time_contact_absolute = time_contact_absolute,
    contact_risk_level    = as.integer(risk_levels),
    contact_risk_category = contact_risk$labels[risk_levels],
    relative_risk         = contact_risk$relative_risk[risk_levels],
    transmission_prob     = transmission_prob,
    traced                = traced,
    blocked_by            = blocked_by,
    stringsAsFactors      = FALSE
  )
}
