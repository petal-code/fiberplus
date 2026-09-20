## =============================================================================
## How much transmission happens before the infector shows symptoms?
##
## fiber draws a case's contact times from the generation-time distribution and
## its incubation period from a separate, independent distribution. Nothing ties
## the two together, so contacts routinely fall before the infector's symptom
## onset. The size of that overlap is not a parameter anyone sets -- it is an
## emergent consequence of the two distributions, and with plausible filovirus
## parameters it can be anywhere from a quarter to a half of all transmission.
##
## This matters for any intervention keyed to symptom onset (contact tracing
## leading to faster admission, isolation, care-seeking): such interventions can
## only ever act on the post-onset share, so the presymptomatic share caps what
## they can achieve. Quantify it before interpreting a tracing scenario.
##
## `presymptomatic_transmission = FALSE` removes it, by truncating contact times
## to start at the end of the parent's incubation period rather than at their
## infection (see the offspring functions).
## =============================================================================

#' Estimate the share of transmission that happens before symptom onset
#'
#' Monte-Carlo estimate of the fraction of a case's transmission that occurs before
#' that case develops symptoms, for a given parameter set. Contact times are drawn
#' from a Gamma generation-time distribution truncated to \eqn{[0, t_{outcome}]}, so
#' for one case the presymptomatic share is
#' \deqn{F(t_{incubation}) / F(t_{outcome})}
#' where \eqn{F} is the generation-time CDF. This function averages that over
#' simulated natural histories.
#'
#' The natural history is replayed the same way \code{\link{complete_offspring_info}}
#' resolves it: symptomatic status, a would-be community outcome, potential
#' hospitalisation, and the hospital second-chance on death. It is evaluated at
#' \eqn{t = 0}, so time-varying inputs are taken at the start of the simulation.
#'
#' Reported separately for genPop and HCW parents, since they have their own
#' generation-time distributions.
#'
#' @param args Named list of arguments as passed to
#'   \code{\link{branching_process_main}}. Needs the natural-history distributions,
#'   \code{prob_symptomatic}, the death and hospitalisation probabilities, and the
#'   \code{Tg_*} generation-time parameters.
#' @param n Integer, number of Monte-Carlo draws. Defaults to 50000.
#' @param seed Optional integer seed. Supplying one makes the estimate reproducible;
#'   note it sets the global seed, so pass \code{NULL} inside a simulation.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{genPop}, \code{hcw}}{Expected share of a case's transmission falling
#'       before its own symptom onset, per route.}
#'     \item{\code{genPop_symptomatic_only}, \code{hcw_symptomatic_only}}{The same,
#'       restricted to cases that actually develop symptoms. For asymptomatic cases the
#'       incubation period still marks the end of the latent period, but calling that
#'       transmission "presymptomatic" is a stretch, so both framings are reported.}
#'     \item{\code{mean_incubation}, \code{mean_time_to_outcome}}{Diagnostics that explain
#'       the result: the overlap is driven by the incubation period relative to the
#'       generation time.}
#'   }
#'
#' @examples
#' \dontrun{
#' ps <- approx_presymptomatic_transmission(args)
#' ps$genPop   # e.g. 0.37 -- over a third of transmission precedes symptoms
#' }
#' @export
approx_presymptomatic_transmission <- function(args, n = 50000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  at_t0 <- function(x, default = NULL) {
    if (is.null(x)) return(default)
    if (is.function(x)) return(x(0))
    x
  }

  needed <- c("incubation_period", "onset_to_death", "onset_to_recovery",
              "onset_to_hospitalisation", "hospitalisation_to_death",
              "hospitalisation_to_recovery")
  missing <- needed[!vapply(needed, function(nm) is.function(args[[nm]]), logical(1))]
  if (length(missing)) {
    stop(sprintf("`args` is missing the delay distribution function(s): %s.",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  prob_symptomatic <- at_t0(args$prob_symptomatic, 1)
  prob_death_comm  <- at_t0(args$prob_death_comm, 0)
  prob_death_hosp  <- at_t0(args$prob_death_hosp, 0)
  hdf              <- at_t0(args$hospitalisation_delay_factor, 1)

  ## --- Natural history, mirroring complete_offspring_info() -------------------
  T_incub     <- args$incubation_period(n)
  symptomatic <- as.logical(rbinom(n, 1, prob_symptomatic))

  ## Would-be community outcome, ignoring hospitalisation for now.
  would_die_comm <- symptomatic & as.logical(rbinom(n, 1, prob_death_comm))
  T_comm_out <- T_incub
  if (any(would_die_comm)) {
    T_comm_out[would_die_comm] <- T_comm_out[would_die_comm] +
      args$onset_to_death(sum(would_die_comm))
  }
  if (any(!would_die_comm)) {
    T_comm_out[!would_die_comm] <- T_comm_out[!would_die_comm] +
      args$onset_to_recovery(sum(!would_die_comm))
  }

  ## Outcome time per route. Hospitalisation probabilities differ between genPop and
  ## HCW, and admission changes the time to outcome, so resolve each route separately.
  outcome_time <- function(prob_hosp) {
    potentially_hosp <- symptomatic & as.logical(rbinom(n, 1, prob_hosp))
    T_hosp <- rep(NA_real_, n)
    if (any(potentially_hosp)) {
      T_hosp[potentially_hosp] <- T_incub[potentially_hosp] +
        args$onset_to_hospitalisation(sum(potentially_hosp)) * hdf
    }
    realised_hosp <- potentially_hosp & !is.na(T_hosp) & (T_hosp < T_comm_out)

    second_chance <- if (prob_death_comm > 0) prob_death_hosp / prob_death_comm else 0
    dies_in_hosp <- logical(n)
    idx <- which(realised_hosp & would_die_comm)
    if (length(idx) > 0) {
      dies_in_hosp[idx] <- as.logical(rbinom(length(idx), 1, second_chance))
    }

    T_out <- T_comm_out
    idx <- which(realised_hosp & dies_in_hosp)
    if (length(idx) > 0) {
      T_out[idx] <- T_hosp[idx] + args$hospitalisation_to_death(length(idx))
    }
    idx <- which(realised_hosp & !dies_in_hosp)
    if (length(idx) > 0) {
      T_out[idx] <- T_hosp[idx] + args$hospitalisation_to_recovery(length(idx))
    }
    T_out
  }

  ## --- Presymptomatic share per route ----------------------------------------
  ## Contact times are Gamma truncated to [0, T_out], so the mass falling before
  ## symptom onset is F(T_incub) / F(T_out). Cases whose total mass is numerically
  ## zero (an outcome essentially at time zero) contribute nothing and are dropped.
  share <- function(T_out, shape, rate) {
    if (is.null(shape) || is.null(rate)) return(c(all = NA_real_, symp = NA_real_))
    F_out   <- pgamma(T_out,   shape = shape, rate = rate)
    F_onset <- pgamma(T_incub, shape = shape, rate = rate)
    ok <- F_out > .Machine$double.eps
    frac <- rep(NA_real_, n)
    frac[ok] <- F_onset[ok] / F_out[ok]
    c(all  = mean(frac[ok]),
      symp = if (any(ok & symptomatic)) mean(frac[ok & symptomatic]) else NA_real_)
  }

  T_out_genPop <- outcome_time(at_t0(args$prob_hospitalised_genPop, 0))
  T_out_hcw    <- outcome_time(at_t0(args$prob_hospitalised_hcw, 0))

  s_genPop <- share(T_out_genPop, args$Tg_shape_genPop, args$Tg_rate_genPop)
  s_hcw    <- share(T_out_hcw,    args$Tg_shape_hcw,    args$Tg_rate_hcw)

  list(
    genPop                  = unname(s_genPop["all"]),
    hcw                     = unname(s_hcw["all"]),
    genPop_symptomatic_only = unname(s_genPop["symp"]),
    hcw_symptomatic_only    = unname(s_hcw["symp"]),
    mean_incubation         = mean(T_incub),
    mean_time_to_outcome    = mean(T_out_genPop),
    n                       = n
  )
}
