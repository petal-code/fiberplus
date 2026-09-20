#' Generate offspring for a deceased parent via an unsafe funeral event
#'
#' Simulates the contacts made at a funeral for a deceased parent, and the subset of those
#' contacts that become secondary infections. Funeral contacts are generated only if the
#' parent died; whether the funeral was safe or unsafe is resolved upstream in
#' \code{\link{complete_offspring_info}} and stored on the parent as
#' \code{funeral_safety}.
#'
#' The number of **contacts** at the funeral is drawn from a Negative Binomial distribution.
#' Contact times are a **Gamma-distributed delay** after the parent's outcome time, with
#' small variance to represent a concentrated exposure event with modest jitter. All funeral
#' contacts occur in the \code{"funeral"} setting.
#'
#' Each contact is assigned a class (\code{"HCW"} or \code{"genPop"}) and a **risk tier**
#' from \code{contact_risk_funeral}, and transmits with probability
#' \code{baseline_risk_funeral(t) * relative_risk[tier]}. Contacts that transmit are then
#' thinned by \code{safe_funeral_efficacy} if the funeral was a safe one. Isolation does not
#' apply on this route -- the parent is dead.
#'
#' Every contact is independently marked as traced, with probability
#' \code{trace_coverage(t) * trace_prob[tier]}.
#'
#' @param parent_info One-row data.frame/list of parent attributes. Funeral transmission reads
#'   \code{hospitalisation}, \code{time_hospitalisation_relative}, \code{time_outcome_relative},
#'   \code{outcome} (did the parent die?), \code{class} ("genPop"/"HCW"), \code{funeral_safety}
#'   ("safe"/"unsafe") and \code{time_infection_absolute}.
#' @param mn_contacts_funeral Positive numeric or function(t). Mean of the Negative Binomial
#'   distribution for the number of contacts at the funeral. Resolved at the parent's absolute
#'   death (outcome) time, since the funeral occurs then rather than at the parent's infection
#'   time.
#' @param overdisp_contacts_funeral Positive numeric. Negative Binomial size (overdispersion)
#'   of the funeral contact distribution.
#' @param baseline_risk_funeral Numeric in \code{[0,1]} or function(t). Per-contact
#'   transmission probability for the reference risk tier at a funeral, resolved at each
#'   contact's own calendar time. Typically higher than the community route -- body handling
#'   is a more intense exposure.
#' @param contact_risk_funeral A \code{\link{make_contact_risk}} structure (or a named list of
#'   its arguments) for the funeral route. Usually distinct from the community structure.
#' @param Tg_shape_funeral Positive numeric. Shape of the Gamma delay distribution from
#'   outcome to funeral contacts.
#' @param Tg_rate_funeral Positive numeric. Rate of the Gamma delay distribution
#'   (mean delay = shape / rate).
#' @param trace_coverage Numeric in \code{[0,1]} or function(t). Programme-level contact
#'   tracing coverage, multiplying each tier's \code{trace_prob}.
#' @param return_contact_log Logical scalar. \code{FALSE} skips building the contact log for
#'   this parent and returns an empty one, saving the per-parent data.frame construction. See
#'   \code{\link{branching_process_main}} for why that matters at scale.
#' @param safe_funeral_efficacy Numeric in \code{[0,1]}. Efficacy of a safe burial in
#'   preventing transmission (1 = fully blocking, 0 = equivalent to an unsafe funeral). The
#'   time-varying probability that a funeral is safe or unsafe is resolved upstream.
#' @param obv_pep_enabled,obv_pep_coverage,obv_pep_adherence,obv_pep_dpc,obv_pep_dpc_sd Obeldesivir
#'   PEP gate settings; see \code{\link{offspring_function_genPop}}. By default
#'   \code{obv_pep_target_locations = "hospital"}, so funeral exposures are not eligible.
#' @param obv_pep_efficacy,obv_pep_efficacy_args,obv_pep_target_class,obv_pep_target_locations Further
#'   obeldesivir PEP gate settings; see \code{\link{offspring_function_genPop}}.
#' @param prob_hcw_cond_funeral_hcw Numeric in \code{[0,1]}. Probability a funeral contact of
#'   an HCW parent is an HCW.
#' @param prob_hcw_cond_funeral_genPop Numeric in \code{[0,1]}. Probability a funeral contact
#'   of a genPop parent is an HCW.
#'
#' @return A data.frame with one row per realised funeral infection, plus the full contact log
#'   attached as the \code{"contact_log"} attribute. Returns a 0-row data.frame if the parent
#'   survived or no contacts occurred.
#' @export
offspring_function_funeral <- function(

  ## *Parent* characteristics and properties
  parent_info,

  ## Contact distribution and per-contact risk at the funeral
  mn_contacts_funeral = NULL,        # scalar or function(t): mean CONTACT distribution (resolved at parent death time)
  overdisp_contacts_funeral = NULL,  # overdispersion of the contact distribution
  baseline_risk_funeral = NULL,      # scalar or function(t): per-contact transmission prob for the reference risk tier
  contact_risk_funeral = NULL,       # fiber_contact_risk structure

  ## Timing of funeral contacts (delay after outcome)
  Tg_shape_funeral = NULL, # gamma shape parameter for the delay distribution ### have high shape, high rate to get low variance ##
  Tg_rate_funeral = NULL,  # gamma rate parameter for the delay distribution

  ## Contact tracing
  trace_coverage = 0,
  return_contact_log = TRUE,  # FALSE skips building the per-parent contact log

  ### efficacy of a safe funeral (thinning funeral offspring)
  safe_funeral_efficacy = NULL, ## efficacy of a safe burial in reducing transmission in a funeral setting

  ## Obeldesivir PEP for exposed HCWs
  obv_pep_enabled = FALSE,
  obv_pep_coverage = 0,
  obv_pep_adherence = 1,
  obv_pep_dpc = 1,
  obv_pep_dpc_sd = NULL,
  obv_pep_efficacy = NULL,
  obv_pep_efficacy_args = NULL,
  obv_pep_target_class = "HCW",
  obv_pep_target_locations = "hospital",

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral_hcw = NULL, ### probability that the unsafe funeral infector infects a HCW
  prob_hcw_cond_funeral_genPop = NULL
) {

  ### note: high probability that this is genPop to genPop (which is what we want) but should allow for parent to be HCW

  ## Step 0: Extract parent info
  parent_hospitalised = parent_info$hospitalisation                          # whether the parent (infector) is hospitalised or not
  parent_time_to_hospitalisation = parent_info$time_hospitalisation_relative # if parent is hospitalised, the time of hospitalisation (relative to infection)
  parent_time_to_outcome = parent_info$time_outcome_relative                 # the time when the parent dies/recovers (relative to time of infection)
  parent_died = parent_info$outcome                                          # was the outcome at parent_time_to_outcome death? (yes/no)
  parent_class = parent_info$class                                           # "genPop" or "HCW" this matters now because unlike other two functions (although improbable) parent could be either
  parent_funeral = parent_info$funeral_safety                                # whether the funeral was safe or not (for those dying)

  #########################################################################################
  ## Input checks
  #########################################################################################

  if (is.null(parent_died) || !is.logical(parent_died) || length(parent_died) != 1L || is.na(parent_died))
    stop("`parent_died` must be a single logical value.", call. = FALSE)

  if (is.null(parent_hospitalised) || !is.logical(parent_hospitalised) || length(parent_hospitalised) != 1L || is.na(parent_hospitalised))
    stop("`parent_hospitalised` must be TRUE or FALSE.", call. = FALSE)

  if (is.null(parent_class) ||
      !parent_class %in% c("genPop", "HCW") ||
      length(parent_class) != 1L) {
    stop("`parent_class` must be either 'genPop' or 'HCW'.", call. = FALSE)
  }

  if (is.null(parent_time_to_outcome) || !is.numeric(parent_time_to_outcome) ||
      length(parent_time_to_outcome) != 1L || is.na(parent_time_to_outcome) ||
      parent_time_to_outcome < 0)
    stop("`parent_time_to_outcome` must be a single non-negative numeric value.", call. = FALSE)

  if (isTRUE(parent_hospitalised)) {
    if (is.null(parent_time_to_hospitalisation) || !is.numeric(parent_time_to_hospitalisation) ||
        length(parent_time_to_hospitalisation) != 1L || is.na(parent_time_to_hospitalisation) ||
        parent_time_to_hospitalisation < 0)
      stop("If parent is hospitalised, `parent_time_to_hospitalisation` must be non-negative.", call. = FALSE)
  } else {
    if (!is.na(parent_time_to_hospitalisation))
      stop("If parent is NOT hospitalised, `parent_time_to_hospitalisation` must be NA.", call. = FALSE)
  }

  # Other positive parameters (scalars). mn_contacts_funeral is handled separately
  # below because it may be a scalar OR a function(t).
  for (nm in c("overdisp_contacts_funeral",
               "Tg_shape_funeral", "Tg_rate_funeral")) {
    val <- get(nm, inherits = FALSE)
    if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop(sprintf("`%s` must be a single positive numeric value.", nm), call. = FALSE)
  }

  # mn_contacts_funeral may be a single positive scalar or a function(t) of
  # absolute calendar time (resolved at the parent's death time, see Step 3).
  if (!is.function(mn_contacts_funeral)) {
    if (is.null(mn_contacts_funeral) || !is.numeric(mn_contacts_funeral) ||
        length(mn_contacts_funeral) != 1L || is.na(mn_contacts_funeral) ||
        mn_contacts_funeral <= 0)
      stop("`mn_contacts_funeral` must be a function(t) or a single positive numeric value.",
           call. = FALSE)
  }

  if (is.null(baseline_risk_funeral)) {
    stop("`baseline_risk_funeral` is required (solve it from a target R0 with solve_baseline_risk_for_r0()).",
         call. = FALSE)
  }
  if (!is.function(baseline_risk_funeral) &&
      (!is.numeric(baseline_risk_funeral) || length(baseline_risk_funeral) != 1L ||
       is.na(baseline_risk_funeral) || baseline_risk_funeral < 0 || baseline_risk_funeral > 1)) {
    stop("`baseline_risk_funeral` must be a function(t) or a single numeric in [0, 1].",
         call. = FALSE)
  }
  contact_risk <- as_contact_risk(contact_risk_funeral, NULL, "contact_risk_funeral")

  if (!is.function(trace_coverage) &&
      (!is.numeric(trace_coverage) || length(trace_coverage) != 1L || is.na(trace_coverage) ||
       trace_coverage < 0 || trace_coverage > 1)) {
    stop("`trace_coverage` must be a function(t) or a single numeric in [0, 1].", call. = FALSE)
  }

  ## Keep safe funeral efficacy as a scalar here. The probability that a funeral
  ## is safe/unsafe is already resolved upstream in complete_offspring_info();
  ## this parameter only controls residual transmission if a funeral is labelled
  ## safe. Set safe_funeral_efficacy = 1 for fully transmission-blocking safe funerals.
  if (is.null(safe_funeral_efficacy) || !is.numeric(safe_funeral_efficacy) ||
      length(safe_funeral_efficacy) != 1L || is.na(safe_funeral_efficacy) ||
      safe_funeral_efficacy < 0 || safe_funeral_efficacy > 1) {
    stop("`safe_funeral_efficacy` must be a single numeric value in [0, 1].", call. = FALSE)
  }

  if (is.null(prob_hcw_cond_funeral_hcw) || !is.numeric(prob_hcw_cond_funeral_hcw) ||
      length(prob_hcw_cond_funeral_hcw) != 1L || is.na(prob_hcw_cond_funeral_hcw) ||
      prob_hcw_cond_funeral_hcw < 0 || prob_hcw_cond_funeral_hcw > 1)
    stop("`prob_hcw_cond_funeral_hcw` must be in [0, 1].", call. = FALSE)

  if (is.null(prob_hcw_cond_funeral_genPop) || !is.numeric(prob_hcw_cond_funeral_genPop) ||
      length(prob_hcw_cond_funeral_genPop) != 1L || is.na(prob_hcw_cond_funeral_genPop) ||
      prob_hcw_cond_funeral_genPop < 0 || prob_hcw_cond_funeral_genPop > 1)
    stop("`prob_hcw_cond_funeral_genPop` must be in [0, 1].", call. = FALSE)


  #########################################################################################
  ## Logic for whether a safe or unsafe funeral occurs
  #########################################################################################

  # Step 1: Ensure parent is dead. If parent survived, no funeral transmission
  if (!isTRUE(parent_died)) {
    return(empty_offspring_dataframe())
  }

  # Step 2: Information on whether the parent had an unsafe or safe funeral
  has_unsafe_funeral <- parent_funeral
  has_unsafe_funeral <- ifelse(has_unsafe_funeral == "unsafe", TRUE, FALSE)

  #########################################################################################
  ## Produce funeral contacts
  #########################################################################################

  # Step 3: Draw the raw number of CONTACTS at the funeral. mn_contacts_funeral may be
  #         time-varying. The funeral is a discrete event at the parent's death, so resolve
  #         it at the parent's absolute death (outcome) time rather than their infection time.
  parent_death_time_absolute <- parent_info$time_infection_absolute + parent_time_to_outcome
  mn_contacts_t <- resolve_positive_time_varying(
    mn_contacts_funeral, parent_death_time_absolute, "mn_contacts_funeral")
  num_contacts <- rnbinom(
    n    = 1,
    mu   = mn_contacts_t,
    size = overdisp_contacts_funeral
  )

  if (num_contacts == 0L) {
    return(empty_offspring_dataframe())
  }

  # Step 4: Generate contact times = outcome time + Gamma distributed 'delay', typically with
  #         little variance to represent a singular point from which infections arose
  delay_funeral <- rgamma(n = num_contacts, shape = Tg_shape_funeral, rate = Tg_rate_funeral)
  contact_times <- parent_time_to_outcome + delay_funeral
  contact_times_absolute <- parent_info$time_infection_absolute + contact_times

  # Step 5: Assign setting ("funeral") for all contacts
  contact_settings <- rep("funeral", num_contacts)

  # Step 6: Assign class (HCW or genPop) to each contact.
  contact_class <- rep("genPop", num_contacts)
  if (parent_class == "genPop") {
    flip_hcw <- as.logical(rbinom(n = num_contacts, size = 1, prob = prob_hcw_cond_funeral_genPop))
  } else if (parent_class == "HCW") {
    flip_hcw <- as.logical(rbinom(n = num_contacts, size = 1, prob = prob_hcw_cond_funeral_hcw))
  } else {
    stop("Step 6 of funeral offspring function is broken")
  }
  contact_class[flip_hcw] <- "HCW"

  # Step 7: Assign a risk tier to each contact, then decide which contacts transmit.
  risk_levels <- draw_contact_risk_levels(num_contacts, contact_risk)
  p_transmit <- contact_transmission_prob(risk_levels, contact_times_absolute,
                                          baseline_risk_funeral, contact_risk,
                                          "baseline_risk_funeral")
  transmitted <- as.logical(rbinom(n = num_contacts, size = 1, prob = p_transmit))

  # Step 8: Contact tracing, drawn for every funeral contact.
  p_trace <- contact_trace_prob(risk_levels, contact_times_absolute,
                                trace_coverage, contact_risk, "trace_coverage")
  traced <- as.logical(rbinom(n = num_contacts, size = 1, prob = p_trace))

  trans_idx <- which(transmitted)
  if (length(trans_idx) == 0L) {
    out <- empty_offspring_dataframe()
    attr(out, "contact_log") <- if (!return_contact_log) empty_contact_log() else new_contact_log(
      parent_id             = parent_info$id,
      offspring_class       = contact_class,
      infection_location    = contact_settings,
      time_contact_relative = contact_times,
      time_contact_absolute = contact_times_absolute,
      risk_levels           = risk_levels,
      transmission_prob     = p_transmit,
      traced                = traced,
      contact_risk          = contact_risk,
      transmitted           = transmitted,
      kept                  = rep(FALSE, num_contacts),
      realised              = rep(FALSE, num_contacts),
      intervention_label    = "safe_funeral"
    )
    return(out)
  }

  infection_times          <- contact_times[trans_idx]
  infection_times_absolute <- contact_times_absolute[trans_idx]
  infection_settings       <- contact_settings[trans_idx]
  offspring_class          <- contact_class[trans_idx]

  # Step 9: Snapshot the pre-thinning candidate set for the two-phase OBV PEP gate.
  pre_thinning <- list(
    infection_location      = infection_settings,
    offspring_class         = offspring_class,
    infection_time_absolute = infection_times_absolute
  )

  # Step 10: Thin if the funeral is a safe one (class-independent thinning by
  #          safe_funeral_efficacy). OBV efficacy is applied separately below.
  n_trans <- length(trans_idx)
  if (!has_unsafe_funeral) {
    keep_infection <- as.logical(rbinom(n = n_trans, size = 1, prob = 1 - safe_funeral_efficacy))
  } else {
    keep_infection <- rep(TRUE, n_trans)
  }

  # Step 11: OBV PEP gate. By default this does nothing for funeral exposures because
  #          obv_pep_target_locations = "hospital", but enabling e.g.
  #          obv_pep_target_locations = c("hospital", "funeral") turns it on here too.
  obv_gate <- apply_obv_pep_gate(
    pre_thinning             = pre_thinning,
    kept_indices             = which(keep_infection),
    obv_pep_enabled          = obv_pep_enabled,
    obv_pep_coverage         = obv_pep_coverage,
    obv_pep_adherence        = obv_pep_adherence,
    obv_pep_dpc              = obv_pep_dpc,
    obv_pep_dpc_sd           = obv_pep_dpc_sd,
    obv_pep_efficacy         = obv_pep_efficacy,
    obv_pep_efficacy_args    = obv_pep_efficacy_args,
    obv_pep_target_class     = obv_pep_target_class,
    obv_pep_target_locations = obv_pep_target_locations
  )

  ## Final realised set = transmitted AND kept AND not prevented by OBV.
  final_local  <- obv_gate$keep
  final_idx    <- which(keep_infection)[final_local]
  obv_metadata <- obv_gate$metadata[final_local, , drop = FALSE]

  kept_full     <- rep(FALSE, num_contacts)
  kept_full[trans_idx] <- keep_infection
  realised_full <- rep(FALSE, num_contacts)
  realised_full[trans_idx[final_idx]] <- TRUE

  offspring_df <- data.frame(infection_location      = infection_settings[final_idx],
                             time_infection_relative = infection_times[final_idx],
                             class                   = offspring_class[final_idx],
                             contact_risk_level      = risk_levels[trans_idx[final_idx]],
                             contact_risk_category   = contact_risk$labels[risk_levels[trans_idx[final_idx]]],
                             traced                  = traced[trans_idx[final_idx]],
                             obv_metadata,
                             stringsAsFactors = FALSE)
  attr(offspring_df, "obv_pep_num_treated") <- obv_gate$num_treated
  attr(offspring_df, "obv_pep_prevented_info") <-
    extract_obv_prevented_info(pre_thinning, keep_infection, final_local)
  attr(offspring_df, "contact_log") <- if (!return_contact_log) empty_contact_log() else new_contact_log(
    parent_id             = parent_info$id,
    offspring_class       = contact_class,
    infection_location    = contact_settings,
    time_contact_relative = contact_times,
    time_contact_absolute = contact_times_absolute,
    risk_levels           = risk_levels,
    transmission_prob     = p_transmit,
    traced                = traced,
    contact_risk          = contact_risk,
    transmitted           = transmitted,
    kept                  = kept_full,
    realised              = realised_full,
    intervention_label    = "safe_funeral"
  )
  return(offspring_df)
}
