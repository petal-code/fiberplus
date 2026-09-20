## =============================================================================
## Single-type (genPop-dominant) R0 approximation at t = 0, and the inversion
## that turns a target (R0, prop_funeral) into baseline per-contact risks.
##
## This is the in-package version of the approximation already used for ABC
## calibration, updated for the contact-first parameterisation and extended to
## cover contact tracing. The structure is unchanged:
##
##     R0_direct  = mn_contacts_genPop  * baseline_risk_genPop  * rr_bar_genPop  * D
##     R0_funeral = mn_contacts_funeral * baseline_risk_funeral * rr_bar_funeral * F
##
## where rr_bar is a risk structure's mean relative risk, D is the direct
## multiplier (what survives hospital quarantine) and F is the
## funeral multiplier (chance of dying, times what survives safe burial).
##
## As before, the expensive Monte-Carlo quantities are EFFICACY-INDEPENDENT and
## are computed once by compute_r0_invariants(); D and F are then cheap
## closed-form functions of the invariants plus a particle's efficacies, so an
## ABC loop varying etu_efficacy / safe_funeral_efficacy
## re-evaluates only the closed forms.
##
## Conditioning convention (unchanged): prob_death_comm and prob_hospitalised_*
## are used DIRECTLY as P(event | symptomatic), matching complete_offspring_info().
## =============================================================================

## Resolve a parameter that may be a scalar OR a function(t), at t = 0.
.at_t0 <- function(x) {
  if (is.null(x))     return(NULL)
  if (is.function(x)) return(x(0))
  x
}

## hospital_quarantine_efficacy(0): the prop_etu(0)-weighted mixture of the fixed
## ETU and general-hospital scalars.
.hq_eff0_from_parts <- function(prop_etu_0, etu_efficacy,
                                general_hospital_quarantine_efficacy) {
  prop_etu_0 * etu_efficacy +
    (1 - prop_etu_0) * general_hospital_quarantine_efficacy
}

#' Efficacy-independent invariants for the R0 approximation
#'
#' Monte-Carlo simulation of a genPop case's natural history at \code{t = 0},
#' mirroring \code{\link{complete_offspring_info}}, to obtain the generation-time
#' mass fractions that the intervention efficacies act on. None of the returned
#' quantities depend on the conditional efficacies (ETU, general hospital,
#' safe burial), so they can be cached per scenario and reused across
#' every particle of a calibration.
#'
#' \code{Q_g} is the expected fraction of a case's generation-time mass falling
#' \emph{after} hospital admission, which is what hospital quarantine acts on.
#'
#' Contact tracing enters through admission timing: a traced case is admitted a flat
#' \code{onset_to_hospitalisation_traced} days after onset (capping its own delay), and
#' may be admitted more often. Both raise \code{Q_g}, since more of the case's
#' generation-time mass falls after admission where quarantine can act on it.
#'
#' Whether a case is traced depends on the risk tier of the contact that infected it.
#' Cases over-represent high-risk tiers, so the tier is drawn from the structure's
#' \code{case_weights} (which are proportional to \code{fractions * relative_risk})
#' rather than from the raw contact fractions.
#'
#' @param args Named list of arguments as passed to
#'   \code{\link{branching_process_main}}.
#' @param n Integer, number of Monte-Carlo draws. Defaults to 50000.
#' @param seed Optional integer seed.
#'
#' @return A named list of invariants: \code{Q_g}, \code{p_realised_hosp},
#'   \code{p_traced}, \code{p_die_comm}, \code{p_die_hosp}, and the resolved \code{t = 0} inputs
#'   \code{prop_etu_0}, \code{p_uf_cgp_0}, \code{p_uf_hgp_0}.
#' @export
compute_r0_invariants <- function(args, n = 50000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  prob_hosp_g <- .at_t0(args$prob_hospitalised_genPop)
  if (is.null(prob_hosp_g)) {
    stop("`args$prob_hospitalised_genPop` is missing.", call. = FALSE)
  }
  hdf <- .at_t0(args$hospitalisation_delay_factor)
  if (is.null(hdf)) hdf <- 1.0

  ## Contact tracing inputs. All default to "off" so a scenario with no tracing
  ## reproduces the untraced approximation exactly.
  risk <- as_contact_risk(args$contact_risk_genPop, args$contact_risk, "contact_risk_genPop")
  trace_cov_0 <- .at_t0(args$trace_coverage)
  if (is.null(trace_cov_0)) trace_cov_0 <- 0
  hosp_mult_traced <- .at_t0(args$prob_hospitalised_multiplier_traced)
  if (is.null(hosp_mult_traced)) hosp_mult_traced <- 1
  hosp_abs_traced <- .at_t0(args$prob_hospitalised_traced)
  traced_delay <- .at_t0(args$onset_to_hospitalisation_traced)

  ## 1. Incubation period
  T_incub <- args$incubation_period(n)

  ## 2. Symptomatic?
  symptomatic <- as.logical(rbinom(n, 1, args$prob_symptomatic))

  ## 3. Tier this case was infected in, and whether it was traced. Tiers are
  ## drawn from the CASE-weighted distribution, not the contact fractions.
  tier <- sample.int(risk$n_levels, size = n, replace = TRUE, prob = risk$case_weights)
  traced <- as.logical(rbinom(n, 1, pmin(trace_cov_0 * risk$trace_prob[tier], 1)))

  ## 4. Would-be community outcome (death/recovery, ignoring hospitalisation).
  would_die_comm <- symptomatic & as.logical(rbinom(n, 1, args$prob_death_comm))

  T_comm_out <- T_incub
  if (any(would_die_comm)) {
    T_comm_out[would_die_comm] <- T_comm_out[would_die_comm] +
      args$onset_to_death(sum(would_die_comm))
  }
  if (any(!would_die_comm)) {
    T_comm_out[!would_die_comm] <- T_comm_out[!would_die_comm] +
      args$onset_to_recovery(sum(!would_die_comm))
  }

  ## 5. Potential hospitalisation. Traced cases may be admitted more often and
  ## sooner (both multipliers default to 1, i.e. no effect).
  p_hosp_i <- rep(prob_hosp_g, n)
  p_hosp_i[traced] <- if (!is.null(hosp_abs_traced)) {
    hosp_abs_traced
  } else {
    pmin(prob_hosp_g * hosp_mult_traced, 1)
  }
  potentially_hosp <- symptomatic & as.logical(rbinom(n, 1, p_hosp_i))

  T_hosp <- rep(NA_real_, n)
  if (any(potentially_hosp)) {
    own_delay <- args$onset_to_hospitalisation(sum(potentially_hosp)) * hdf
    ## Tracing caps the onset-to-admission delay rather than replacing it.
    if (!is.null(traced_delay)) {
      is_traced <- traced[potentially_hosp]
      own_delay[is_traced] <- pmin(own_delay[is_traced], traced_delay)
    }
    T_hosp[potentially_hosp] <- T_incub[potentially_hosp] + own_delay
  }

  ## 6. Realised hospitalisation: admission must beat the community outcome
  realised_hosp <- potentially_hosp & !is.na(T_hosp) & (T_hosp < T_comm_out)

  ## 7. Outcome status & time
  second_chance_death <- if (args$prob_death_comm > 0) {
    args$prob_death_hosp / args$prob_death_comm
  } else 0
  dies_in_hosp <- logical(n)
  idx <- which(realised_hosp & would_die_comm)
  if (length(idx) > 0) {
    dies_in_hosp[idx] <- as.logical(rbinom(length(idx), 1, second_chance_death))
  }
  outcome_death <- would_die_comm
  outcome_death[realised_hosp] <- dies_in_hosp[realised_hosp]

  T_out <- T_comm_out
  idx <- which(realised_hosp & dies_in_hosp)
  if (length(idx) > 0) {
    T_out[idx] <- T_hosp[idx] + args$hospitalisation_to_death(length(idx))
  }
  idx <- which(realised_hosp & !dies_in_hosp)
  if (length(idx) > 0) {
    T_out[idx] <- T_hosp[idx] + args$hospitalisation_to_recovery(length(idx))
  }

  ## 8. Generation-time mass fractions. Contact times are drawn from a Gamma truncated
  ## to [lower, T_out], so the share of a case's transmission in a window [a, b] is
  ## (F(b) - F(a)) / (F(T_out) - F(lower)). The lower bound is the parent's infection
  ## (0) normally, and their symptom onset when presymptomatic transmission is off --
  ## which concentrates the remaining mass later and so raises Q_g.
  presympt <- args$presymptomatic_transmission
  if (is.null(presympt)) presympt <- TRUE
  gt_cdf <- function(x) pgamma(x, shape = args$Tg_shape_genPop, rate = args$Tg_rate_genPop)
  F_out   <- gt_cdf(T_out)
  F_lower <- if (isTRUE(presympt)) 0 else gt_cdf(T_incub)
  F_total <- F_out - F_lower
  valid   <- F_total > .Machine$double.eps

  ## Post-admission mass (hospital quarantine acts here)
  q_h <- numeric(n)
  sel <- realised_hosp & valid
  if (any(sel)) {
    q_h[sel] <- (F_out[sel] - gt_cdf(T_hosp[sel])) / F_total[sel]
  }

  prop_etu_0 <- .at_t0(args$prop_etu)
  p_uf_cgp_0 <- .at_t0(args$p_unsafe_funeral_comm_genPop)
  p_uf_hgp_0 <- .at_t0(args$p_unsafe_funeral_hosp_genPop)
  if (is.null(prop_etu_0)) stop("`args$prop_etu` is required.", call. = FALSE)
  if (is.null(p_uf_cgp_0)) stop("`args$p_unsafe_funeral_comm_genPop` is required.", call. = FALSE)
  if (is.null(p_uf_hgp_0)) stop("`args$p_unsafe_funeral_hosp_genPop` is required.", call. = FALSE)

  list(
    Q_g             = mean(q_h),
    p_realised_hosp = mean(realised_hosp),
    p_traced        = mean(traced),
    p_die_comm      = mean(outcome_death & !realised_hosp),
    p_die_hosp      = mean(outcome_death &  realised_hosp),
    prop_etu_0      = prop_etu_0,
    p_uf_cgp_0      = p_uf_cgp_0,
    p_uf_hgp_0      = p_uf_hgp_0
  )
}

#' Direct-transmission multiplier D
#'
#' The fraction of a genPop case's would-be transmission that survives the
#' post-admission quarantine layer:
#' \deqn{D = 1 - hq(0) \, Q_g}
#' PPE does not enter:
#' it only thins HCW recipients, and this is the genPop-dominant single-type
#' approximation.
#'
#' @param inv Invariants from \code{\link{compute_r0_invariants}}.
#' @param etu_efficacy,general_hospital_quarantine_efficacy Fixed scalar
#'   post-admission quarantine efficacies.
#' @return A single numeric multiplier.
#' @export
r0_direct_multiplier <- function(inv,
                                 etu_efficacy,
                                 general_hospital_quarantine_efficacy) {
  hq0 <- .hq_eff0_from_parts(inv$prop_etu_0, etu_efficacy,
                             general_hospital_quarantine_efficacy)
  1 - hq0 * inv$Q_g
}

#' Funeral-transmission multiplier F
#'
#' The probability a genPop case dies, weighted by what survives safe burial in
#' each death location.
#'
#' @param inv Invariants from \code{\link{compute_r0_invariants}}.
#' @param safe_funeral_efficacy Fixed scalar efficacy of a safe burial.
#' @return A single numeric multiplier.
#' @export
r0_funeral_multiplier <- function(inv, safe_funeral_efficacy) {
  inv$p_die_comm * (1 - safe_funeral_efficacy * (1 - inv$p_uf_cgp_0)) +
    inv$p_die_hosp * (1 - safe_funeral_efficacy * (1 - inv$p_uf_hgp_0))
}

#' Forward R0 approximation at t = 0
#'
#' Diagnostic: the expected secondary infections from a genPop index case under
#' the contact-first parameterisation, split into the direct (community plus
#' hospital) and funeral routes.
#'
#' @param args Named list of arguments as passed to
#'   \code{\link{branching_process_main}}.
#' @param n,seed Monte-Carlo settings passed to \code{\link{compute_r0_invariants}}.
#' @param invariants Optional pre-computed invariants, to avoid re-running the
#'   Monte Carlo.
#'
#' @return A named list with \code{R0}, \code{R0_direct}, \code{R0_funeral}, the
#'   multipliers \code{D} and \code{F}, and the invariants used.
#' @export
approx_r0 <- function(args, n = 50000, seed = NULL, invariants = NULL) {

  inv <- if (is.null(invariants)) compute_r0_invariants(args, n = n, seed = seed) else invariants

  D <- r0_direct_multiplier(inv, args$etu_efficacy,
                            args$general_hospital_quarantine_efficacy)
  F_fun <- r0_funeral_multiplier(inv, args$safe_funeral_efficacy)

  risk_g <- as_contact_risk(args$contact_risk_genPop,  args$contact_risk, "contact_risk_genPop")
  risk_f <- as_contact_risk(args$contact_risk_funeral, args$contact_risk, "contact_risk_funeral")

  R0_direct <- .at_t0(args$mn_contacts_genPop) * .at_t0(args$baseline_risk_genPop) *
    risk_g$mean_relative_risk * D

  R0_funeral <- .at_t0(args$mn_contacts_funeral) * .at_t0(args$baseline_risk_funeral) *
    risk_f$mean_relative_risk * F_fun

  list(
    R0         = R0_direct + R0_funeral,
    R0_direct  = R0_direct,
    R0_funeral = R0_funeral,
    D          = D,
    F          = F_fun,
    invariants = inv
  )
}

## Internal: turn a required expected-infections figure for one route into the
## baseline per-contact risk that delivers it, checking the top tier stays a
## valid probability.
.baseline_risk_for_route <- function(target, mn_contacts, risk, multiplier, route) {
  if (target == 0) return(0)
  if (!is.finite(multiplier) || multiplier <= 0) {
    stop(sprintf("The %s multiplier is %s (non-positive); cannot invert for `baseline_risk_%s`.",
                 route, signif(multiplier, 4), route), call. = FALSE)
  }
  if (is.null(mn_contacts)) {
    stop(sprintf("`mn_contacts_%s` is required to solve for `baseline_risk_%s`.", route, route),
         call. = FALSE)
  }
  p0 <- target / (mn_contacts * risk$mean_relative_risk * multiplier)
  if (p0 * risk$max_relative_risk > 1 + 1e-12) {
    stop(sprintf(
      paste0("The requested %s contribution (%s infections per case) is not achievable with ",
             "mn_contacts_%s = %s and this risk structure: it needs a baseline risk of %s, ",
             "which puts the highest-risk tier at %s (must be <= 1). The most this route can ",
             "deliver is %s -- raise the mean contact number or narrow the relative-risk spread."),
      route, format(round(target, 4)), route, format(round(mn_contacts, 4)),
      format(round(p0, 6)), format(round(p0 * risk$max_relative_risk, 4)),
      format(round(mn_contacts * risk$mean_relative_risk * multiplier / risk$max_relative_risk, 4))
    ), call. = FALSE)
  }
  p0
}

#' Solve baseline per-contact risks from a target R0
#'
#' Inversion of \code{\link{approx_r0}}: given a target R0, the share of
#' transmission that should come from funerals, and the contact distributions,
#' return the baseline per-contact transmission probabilities for the genPop and
#' funeral routes.
#'
#' \deqn{p_0^{genPop} = (1 - \pi_f) R_0 / (\bar{c}_{genPop} \, \overline{rr}_{genPop} \, D)}
#' \deqn{p_0^{funeral} = \pi_f R_0 / (\bar{c}_{funeral} \, \overline{rr}_{funeral} \, F)}
#'
#' Every tier's transmission probability must remain in \code{[0, 1]}, so the
#' highest-risk tier caps what each route can deliver; an unachievable target is
#' an error naming the ceiling rather than a silently clipped top tier.
#'
#' @param R0 Positive numeric. Target basic reproduction number at \code{t = 0}.
#' @param args Named list of arguments as passed to
#'   \code{\link{branching_process_main}}. Must contain the contact means, the
#'   natural-history distributions, and the efficacies.
#' @param proportion_transmission_from_funerals Numeric in \code{[0, 1]}. Share of
#'   R0 attributed to the funeral route.
#' @param n,seed Monte-Carlo settings passed to \code{\link{compute_r0_invariants}}.
#' @param invariants Optional pre-computed invariants, to avoid re-running the
#'   Monte Carlo across calibration particles.
#'
#' @return A named list with \code{baseline_risk_genPop_required},
#'   \code{baseline_risk_funeral_required}, the multipliers \code{D} and \code{F},
#'   the target inputs, and the invariants used.
#'
#' @examples
#' \dontrun{
#' fit <- solve_baseline_risk_for_r0(R0 = 2, args = args,
#'                                   proportion_transmission_from_funerals = 0.2)
#' args$baseline_risk_genPop  <- fit$baseline_risk_genPop_required
#' args$baseline_risk_funeral <- fit$baseline_risk_funeral_required
#' }
#' @export
solve_baseline_risk_for_r0 <- function(R0,
                                       args,
                                       proportion_transmission_from_funerals,
                                       n          = 50000,
                                       seed       = NULL,
                                       invariants = NULL) {

  if (!is.numeric(R0) || length(R0) != 1L || is.na(R0) || R0 <= 0) {
    stop("`R0` must be a single positive number.", call. = FALSE)
  }
  pi_f <- proportion_transmission_from_funerals
  if (!is.numeric(pi_f) || length(pi_f) != 1L || is.na(pi_f) || pi_f < 0 || pi_f > 1) {
    stop("`proportion_transmission_from_funerals` must be a single number in [0, 1].",
         call. = FALSE)
  }

  inv <- if (is.null(invariants)) compute_r0_invariants(args, n = n, seed = seed) else invariants

  D <- r0_direct_multiplier(inv, args$etu_efficacy,
                            args$general_hospital_quarantine_efficacy)
  F_fun <- r0_funeral_multiplier(inv, args$safe_funeral_efficacy)

  risk_g <- as_contact_risk(args$contact_risk_genPop,  args$contact_risk, "contact_risk_genPop")
  risk_f <- as_contact_risk(args$contact_risk_funeral, args$contact_risk, "contact_risk_funeral")

  p0_g <- .baseline_risk_for_route((1 - pi_f) * R0, .at_t0(args$mn_contacts_genPop),
                                   risk_g, D, "genPop")
  p0_f <- .baseline_risk_for_route(pi_f * R0, .at_t0(args$mn_contacts_funeral),
                                   risk_f, F_fun, "funeral")

  list(
    baseline_risk_genPop_required  = p0_g,
    baseline_risk_funeral_required = p0_f,
    R0_target                      = R0,
    proportion_transmission_from_funerals = pi_f,
    D_direct_multiplier            = D,
    F_funeral_multiplier           = F_fun,
    invariants                     = inv
  )
}
