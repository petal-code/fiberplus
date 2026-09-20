#' Generate offspring for a healthcare-worker (HCW) parent
#'
#' Simulates the contacts made by a **healthcare-worker** parent and the subset of those
#' contacts that become secondary infections.
#'
#' The number of **contacts** is drawn from a Negative Binomial distribution and contact
#' times from a **truncated Gamma** generation-time distribution on
#' \eqn{[0, t_\mathrm{outcome}]}. For HCW parents, **pre-admission** contacts may occur in
#' the **hospital** (while working, with probability
#' \code{prob_hospital_cond_hcw_preAdm}) or in the **community**; contacts at or after the
#' parent's own admission are always in the hospital. Each contact is assigned a class
#' (\code{"HCW"} or \code{"genPop"}) and a **risk tier** from \code{contact_risk_hcw}, and
#' transmits with probability \code{baseline_risk_hcw(t) * relative_risk[tier]} resolved at
#' that contact's own calendar time.
#'
#' Contacts that transmit are thinned by three multiplicatively-combined protective layers,
#' each applied only where relevant:
#' \itemize{
#'   \item **Source PPE** worn by the still-working HCW parent: all pre-admission hospital
#'         events, any recipient class.
#'   \item **Receiver PPE** worn by HCW recipients: all hospital events with an HCW
#'         recipient, pre- or post-admission.
#'   \item **Hospital quarantine**: all post-admission hospital events.
#' }
#'
#' Setting \code{presymptomatic_transmission = FALSE} truncates contact times to start at the
#' end of the parent's incubation period; see \code{\link{offspring_function_genPop}}.
#'
#' Every contact is independently marked as traced, with probability
#' \code{trace_coverage(t) * trace_prob[tier]}; traced status travels with contacts that
#' become cases.
#'
#' @param parent_info One-row data.frame/list containing parent infection, incubation,
#'   hospitalisation and outcome times.
#' @param mn_contacts_hcw Positive numeric or function(t). Mean of the Negative Binomial
#'   contact distribution for HCW parents, resolved at the parent's absolute infection time.
#' @param overdisp_contacts_hcw Positive numeric. Negative Binomial \code{size} of the
#'   contact distribution.
#' @param baseline_risk_hcw Numeric in \code{[0,1]} or function(t). Per-contact transmission
#'   probability for the reference risk tier, resolved at each contact's own calendar time.
#' @param contact_risk_hcw A \code{\link{make_contact_risk}} structure (or a named list of
#'   its arguments) for this route.
#' @param Tg_shape_hcw Positive numeric. Shape of the Gamma generation-time distribution.
#' @param Tg_rate_hcw Positive numeric. Rate of the Gamma generation-time distribution.
#' @param prob_hospital_cond_hcw_preAdm Numeric in \code{[0,1]}. Probability that a
#'   pre-admission contact occurs in the hospital setting while the HCW is working.
#' @param trace_coverage Numeric in \code{[0,1]} or function(t). Programme-level contact
#'   tracing coverage, multiplying each tier's \code{trace_prob}.
#' @param return_contact_log Logical scalar. \code{FALSE} skips building the contact log for
#'   this parent and returns an empty one, saving the per-parent data.frame construction. See
#'   \code{\link{branching_process_main}} for why that matters at scale.
#' @param presymptomatic_transmission Logical scalar. \code{TRUE} (default) lets contacts occur
#'   before the parent's symptom onset; \code{FALSE} truncates contact times to start at the end
#'   of the parent's incubation period.
#' @param ppe_coverage_hcw Numeric in \code{[0,1]} or function(t). Coverage: probability that
#'   an HCW (source or receiver) has PPE. Each PPE layer thins by
#'   \code{ppe_coverage_hcw(t) * ppe_efficacy}.
#' @param ppe_efficacy Numeric in \code{[0,1]} (fixed scalar). Efficacy of PPE conditional on
#'   the HCW having it. Shared by the source and receiver layers.
#' @param prop_etu Numeric in \code{[0,1]} or function(t). Proportion of hospitalised cases
#'   managed in ETU/ETC care.
#' @param etu_efficacy Numeric in \code{[0,1]} (fixed scalar). Post-admission quarantine
#'   efficacy in ETU/ETC care.
#' @param general_hospital_quarantine_efficacy Numeric in \code{[0,1]} (fixed scalar).
#'   Post-admission quarantine efficacy in general hospital care.
#' @param obv_pep_enabled,obv_pep_coverage,obv_pep_adherence,obv_pep_dpc,obv_pep_dpc_sd Obeldesivir
#'   PEP gate settings; see \code{\link{offspring_function_genPop}}.
#' @param obv_pep_efficacy,obv_pep_efficacy_args,obv_pep_target_class,obv_pep_target_locations Further
#'   obeldesivir PEP gate settings; see \code{\link{offspring_function_genPop}}.
#' @param prob_hcw_cond_hcw_comm Numeric in \code{[0,1]}. Probability a community-located
#'   contact of an HCW parent is an HCW.
#' @param prob_hcw_cond_hcw_hospital Numeric in \code{[0,1]}. Probability a hospital-located
#'   contact of an HCW parent is an HCW.
#'
#' @return A data.frame with one row per realised infection, plus the full contact log
#'   attached as the \code{"contact_log"} attribute. See
#'   \code{\link{offspring_function_genPop}} for the column list.
#' @export
offspring_function_hcw <- function(

  ## Characteristics and properties of the parent (who we are generating the contacts for)
  parent_info,

  ## Contact distribution and per-contact risk for HCWs
  mn_contacts_hcw = NULL,                   # scalar or function(t): mean CONTACT distribution (resolved at parent infection time)
  overdisp_contacts_hcw = NULL,             # overdispersion of the contact distribution
  baseline_risk_hcw = NULL,                 # scalar or function(t): per-contact transmission prob for the reference risk tier
  contact_risk_hcw = NULL,                  # fiber_contact_risk structure

  ## Generation time distribution
  Tg_shape_hcw = NULL,                      # gamma shape parameter for Tg distribution for HCWs
  Tg_rate_hcw = NULL,                       # gamma rate parameter for Tg distribution for HCWs

  ## Setting model for HCWs
  prob_hospital_cond_hcw_preAdm = NULL,     # probability that a pre-admission contact occurs in the hospital (whilst HCW is working)

  ## Contact tracing
  trace_coverage = 0,                       # scalar/function(t): programme-level tracing coverage
  return_contact_log = TRUE,                # FALSE skips building the per-parent contact log

  ## Whether contacts can occur before the parent develops symptoms
  presymptomatic_transmission = TRUE,       # FALSE truncates contact times to start at symptom onset

  ## PPE and post-admission quarantine
  ppe_coverage_hcw = NULL,                  # scalar or function(t): coverage/probability that an HCW (source or receiver) has PPE
  ppe_efficacy = NULL,                      # scalar: efficacy of PPE conditional on the HCW having it
  prop_etu = NULL,                          # scalar/function(t): proportion of hospitalised cases in ETU/ETC care
  etu_efficacy = NULL,                      # scalar: post-admission quarantine efficacy for ETU/ETC care
  general_hospital_quarantine_efficacy = NULL,  # scalar: post-admission quarantine efficacy for general hospital care

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

  ## Probabilities for HCW contacting either genPop or HCWs, depending on the setting
  prob_hcw_cond_hcw_comm = NULL,            # prob that a community-located contact of an HCW parent is a HCW
  prob_hcw_cond_hcw_hospital = NULL         # prob that a hospital-located contact of an HCW parent is a HCW
) {

  ## Local helpers for parameters that can be either scalars or functions of
  ## absolute calendar time. The actual interpolation is handled by
  ## resolve_time_varying(), so this function stays scenario-agnostic.
  validate_probability_or_time_varying <- function(param, param_name) {
    if (is.function(param)) {
      return(invisible(NULL))
    }
    if (!is.numeric(param) || length(param) != 1L || is.na(param) || param < 0 || param > 1) {
      stop(sprintf("`%s` must be a function(t) or a single numeric in [0, 1].", param_name),
           call. = FALSE)
    }
    invisible(NULL)
  }

  validate_probability_scalar <- function(param, param_name) {
    if (!is.numeric(param) || length(param) != 1L || is.na(param) || param < 0 || param > 1) {
      stop(sprintf("`%s` must be a single numeric in [0, 1].", param_name),
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

  resolve_hospital_quarantine_efficacy <- function(t) {
    prop_etu_t <- resolve_probability(prop_etu, t, "prop_etu")

    ## Post-admission quarantine efficacy is a coverage-style mixture over where a
    ## hospitalised case is managed: a proportion prop_etu(t) receive the (higher)
    ## ETU/ETC quarantine efficacy and the remainder receive the general-hospital
    ## quarantine efficacy. Both efficacies are fixed scalars; the only time
    ## variation enters through prop_etu(t).
    hospital_quarantine_efficacy_t <-
      prop_etu_t * etu_efficacy +
      (1 - prop_etu_t) * general_hospital_quarantine_efficacy

    if (any(hospital_quarantine_efficacy_t < 0 |
            hospital_quarantine_efficacy_t > 1)) {
      stop(
        "`hospital_quarantine_efficacy(t)` must resolve to values between 0 and 1.",
        call. = FALSE
      )
    }

    hospital_quarantine_efficacy_t
  }

  # Step 0: Extract relevant parent information from parent_info
  parent_hospitalised = parent_info$hospitalisation                          # whether the parent (infector) is hospitalised or not
  parent_time_to_hospitalisation = parent_info$time_hospitalisation_relative # if parent is hospitalised, the time of hospitalisation (relative to infection)
  parent_time_to_outcome = parent_info$time_outcome_relative                 # the time when the parent dies/recovers (relative to time of infection)
  parent_time_infection_absolute = parent_info$time_infection_absolute       # absolute calendar time of parent infection
  parent_incubation_period = parent_info$incubation_period                   # infection -> symptom onset; start of infectiousness when presymptomatic transmission is off

  #########################################################################################
  ## Checks to make sure function inputs are correctly specified
  #########################################################################################
  if (is.null(parent_hospitalised) || length(parent_hospitalised) != 1L || !is.logical(parent_hospitalised) || is.na(parent_hospitalised)) {
    stop("`parent_hospitalised` must be a single logical value: TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(parent_hospitalised)) {
    if (is.null(parent_time_to_hospitalisation) || length(parent_time_to_hospitalisation) != 1L ||
        !is.numeric(parent_time_to_hospitalisation) || is.na(parent_time_to_hospitalisation) || parent_time_to_hospitalisation < 0) {
      stop("When `parent_hospitalised` is TRUE, `parent_time_to_hospitalisation` must be a single non-negative numeric value (not NA).", call. = FALSE)
    }
  }
  if (is.null(parent_time_to_outcome) || length(parent_time_to_outcome) != 1L ||
      !is.numeric(parent_time_to_outcome) || is.na(parent_time_to_outcome) || parent_time_to_outcome < 0) {
    stop("`parent_time_to_outcome` must be a single non-negative numeric value (not NA).", call. = FALSE)
  }
  if (is.null(parent_time_infection_absolute) || length(parent_time_infection_absolute) != 1L ||
      !is.numeric(parent_time_infection_absolute) || is.na(parent_time_infection_absolute)) {
    stop("`parent_info$time_infection_absolute` must be a single non-missing numeric value.", call. = FALSE)
  }
  if (identical(parent_hospitalised, FALSE)) {
    if (!is.na(parent_time_to_hospitalisation)) {
      stop("When `parent_hospitalised` is FALSE, `parent_time_to_hospitalisation` must be NA.", call. = FALSE)
    }
  }

  ## Generation-time parameters and the contact distribution.
  for (nm in c("overdisp_contacts_hcw", "Tg_shape_hcw", "Tg_rate_hcw")) {
    val <- get(nm, inherits = FALSE)
    if (is.null(val) || length(val) != 1L || !is.numeric(val) || is.na(val) || val <= 0) {
      stop(sprintf("`%s` must be a single positive numeric value.", nm), call. = FALSE)
    }
  }
  if (is.null(mn_contacts_hcw)) {
    stop("`mn_contacts_hcw` is required.", call. = FALSE)
  }
  if (!is.function(mn_contacts_hcw) &&
      (!is.numeric(mn_contacts_hcw) || length(mn_contacts_hcw) != 1L ||
       is.na(mn_contacts_hcw) || mn_contacts_hcw <= 0)) {
    stop("`mn_contacts_hcw` must be a function(t) or a single positive numeric value.",
         call. = FALSE)
  }
  if (is.null(baseline_risk_hcw)) {
    stop("`baseline_risk_hcw` is required (solve it from a target R0 with solve_baseline_risk_for_r0()).",
         call. = FALSE)
  }
  validate_probability_or_time_varying(baseline_risk_hcw, "baseline_risk_hcw")
  contact_risk <- as_contact_risk(contact_risk_hcw, NULL, "contact_risk_hcw")

  for (nm in c("prob_hospital_cond_hcw_preAdm", "prob_hcw_cond_hcw_comm", "prob_hcw_cond_hcw_hospital")) {
    val <- get(nm, inherits = FALSE)
    if (is.null(val) || length(val) != 1L || !is.numeric(val) || is.na(val) || val < 0 || val > 1) {
      stop(sprintf("`%s` must be a single numeric value in [0, 1].", nm), call. = FALSE)
    }
  }
  validate_probability_or_time_varying(ppe_coverage_hcw, "ppe_coverage_hcw")
  validate_probability_scalar(ppe_efficacy, "ppe_efficacy")
  validate_probability_or_time_varying(prop_etu, "prop_etu")
  validate_probability_scalar(etu_efficacy, "etu_efficacy")
  validate_probability_scalar(general_hospital_quarantine_efficacy, "general_hospital_quarantine_efficacy")
  validate_probability_or_time_varying(trace_coverage, "trace_coverage")
  if (!is.logical(presymptomatic_transmission) || length(presymptomatic_transmission) != 1L ||
      is.na(presymptomatic_transmission)) {
    stop("`presymptomatic_transmission` must be a single logical value.", call. = FALSE)
  }
  if (!presymptomatic_transmission &&
      (is.null(parent_incubation_period) || length(parent_incubation_period) != 1L ||
       !is.numeric(parent_incubation_period) || is.na(parent_incubation_period) ||
       parent_incubation_period < 0)) {
    stop("`parent_info$incubation_period` must be a single non-negative numeric value when `presymptomatic_transmission = FALSE`.",
         call. = FALSE)
  }

  ########################################################################################################
  ## Generating contacts, contact times, settings, classes and risk tiers
  ########################################################################################################

  # Step 1: Draw the raw number of CONTACTS from the Negative Binomial.
  mn_contacts_t <- resolve_positive_time_varying(
    mn_contacts_hcw, parent_time_infection_absolute, "mn_contacts_hcw")
  num_contacts <- rnbinom(n = 1, mu = mn_contacts_t, size = overdisp_contacts_hcw)

  if (num_contacts == 0L) {
    return(empty_offspring_dataframe())
  }

  # Step 2: Generate the time of each contact from the generation time (truncated at outcome).
  #         With presymptomatic transmission switched off the lower bound moves from the
  #         parent's infection to their symptom onset; see offspring_function_genPop.
  gt_lower <- if (isTRUE(presymptomatic_transmission)) 0 else parent_incubation_period
  contact_times <- rtrunc_gamma(n = num_contacts,
                                lower = gt_lower,
                                upper = parent_time_to_outcome,
                                Tg_shape = Tg_shape_hcw,
                                Tg_rate = Tg_rate_hcw)
  contact_times_absolute <- parent_time_infection_absolute + contact_times

  # Step 3: Generate the setting of each contact. For HCWs, pre-admission contacts can be
  #         community *or* hospital (while working), controlled by prob_hospital_cond_hcw_preAdm.
  #         Post-admission contacts are necessarily in the hospital.
  t_hosp <- if (isTRUE(parent_hospitalised)) parent_time_to_hospitalisation else Inf
  if (is.finite(t_hosp) && parent_time_to_outcome <= t_hosp) t_hosp <- Inf
  pre_admission <- contact_times < t_hosp
  is_hosp_pre_admission <- as.logical(rbinom(n = num_contacts, size = 1, prob = prob_hospital_cond_hcw_preAdm))
  contact_settings <- ifelse(pre_admission, ifelse(is_hosp_pre_admission, "hospital", "community"), "hospital")

  # Step 4: Assign class (HCW or genPop) to each contact based on the setting.
  contact_class <- rep("genPop", num_contacts)
  prob_hcw <- ifelse(contact_settings == "community", prob_hcw_cond_hcw_comm, prob_hcw_cond_hcw_hospital)
  flip_hcw <- as.logical(rbinom(n = num_contacts, size = 1, prob = prob_hcw))
  contact_class[flip_hcw] <- "HCW"

  # Step 5: Assign a risk tier to each contact, then decide which contacts transmit.
  risk_levels <- draw_contact_risk_levels(num_contacts, contact_risk)
  p_transmit <- contact_transmission_prob(risk_levels, contact_times_absolute,
                                          baseline_risk_hcw, contact_risk,
                                          "baseline_risk_hcw")
  transmitted <- as.logical(rbinom(n = num_contacts, size = 1, prob = p_transmit))

  # Step 6: Contact tracing, drawn for every contact (see offspring_function_genPop).
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
      intervention_label    = "ppe_quarantine"
    )
    return(out)
  }

  infection_times          <- contact_times[trans_idx]
  infection_times_absolute <- contact_times_absolute[trans_idx]
  infection_settings       <- contact_settings[trans_idx]
  offspring_class          <- contact_class[trans_idx]
  infection_pre_admission  <- pre_admission[trans_idx]

  # Step 7: Compute per-event keep probability under the protective layers
  #         (Swiss-cheese multiplicative), each applied only where it is relevant.
  #   - Source PPE (worn by the still-working HCW parent): all PRE-admission hospital events.
  #   - Receiver PPE (worn by HCW recipients): all hospital events with an HCW recipient.
  #   - Hospital quarantine: all POST-admission hospital events.
  p_keep_infection <- rep(1, length(infection_times))

  hospital_idx <- which(infection_settings == "hospital")
  if (length(hospital_idx) > 0) {
    pre_adm_in_hospital  <- infection_pre_admission[hospital_idx]
    post_adm_in_hospital <- !pre_adm_in_hospital
    hcw_recv             <- offspring_class[hospital_idx] == "HCW"

    ppe_source_cov_t <- numeric(length(hospital_idx))
    if (any(pre_adm_in_hospital)) {
      ppe_source_cov_t[pre_adm_in_hospital] <- resolve_probability(
        ppe_coverage_hcw,
        infection_times_absolute[hospital_idx[pre_adm_in_hospital]],
        "ppe_coverage_hcw"
      )
    }

    ppe_recv_cov_t <- numeric(length(hospital_idx))
    if (any(hcw_recv)) {
      ppe_recv_cov_t[hcw_recv] <- resolve_probability(
        ppe_coverage_hcw,
        infection_times_absolute[hospital_idx[hcw_recv]],
        "ppe_coverage_hcw"
      )
    }

    hq_t <- numeric(length(hospital_idx))
    if (any(post_adm_in_hospital)) {
      hq_t[post_adm_in_hospital] <- resolve_hospital_quarantine_efficacy(
        infection_times_absolute[hospital_idx[post_adm_in_hospital]]
      )
    }

    p_keep_infection[hospital_idx] <- p_keep_infection[hospital_idx] *
      (1 - ppe_source_cov_t * ppe_efficacy) * (1 - ppe_recv_cov_t * ppe_efficacy) * (1 - hq_t)
  }

  # Step 7b: Snapshot the pre-thinning candidate set for the two-phase OBV PEP gate.
  pre_thinning <- list(
    infection_location      = infection_settings,
    offspring_class         = offspring_class,
    infection_time_absolute = infection_times_absolute
  )

  # Step 8: Thin (Swiss-cheese; OBV efficacy is applied separately, post-thinning, below)
  keep_infection <- as.logical(rbinom(n = length(infection_times), size = 1, prob = p_keep_infection))

  # Step 9: OBV PEP gate (two phases internally, see apply_obv_pep_gate).
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

  # Step 10: Define and output the realised-infection dataframe
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
    intervention_label    = "ppe_quarantine"
  )

  return(offspring_df)
}
