summarise_output <- function(
    tdf,
    subset      = "realised_subset",
    sim_info    = NULL,  # optional list from branching_process_main
    contact_log = NULL   # optional contact log from branching_process_main
) {

  ##--------------------------------------------------------------
  ## 0. Subset choice and basic checks
  ##--------------------------------------------------------------
  if (!(subset %in% c("total_tdf", "realised_subset"))) {
    stop('subset must be either "total_tdf" or "realised_subset"')
  }

  if (subset == "realised_subset") {
    if (is.null(tdf$offspring_generated)) {
      stop('tdf must contain column "offspring_generated" for realised_subset')
    }
    subset_vector <- tdf$offspring_generated == TRUE
  } else {
    subset_vector <- rep(TRUE, nrow(tdf))
  }

  if (!any(subset_vector)) {
    stop("No rows selected by subset; nothing to summarise.")
  }

  ##--------------------------------------------------------------
  ## 1. Outbreak timing (start / end / duration)
  ##--------------------------------------------------------------
  ## Continuous times
  first_infection_time <- min(tdf$time_infection_absolute[subset_vector],
                              na.rm = TRUE)
  last_outcome_time <- max(tdf$time_outcome_absolute[subset_vector],
                           na.rm = TRUE)

  if (!is.finite(first_infection_time) || !is.finite(last_outcome_time)) {
    stop("Could not determine outbreak start/end times (all times NA?).")
  }

  outbreak_duration_cont <- last_outcome_time - first_infection_time
  outbreak_duration_days <- round(outbreak_duration_cont)

  ##--------------------------------------------------------------
  ## 2. Basic counts by class (genPop vs HCW)
  ##--------------------------------------------------------------
  class_vec <- tdf$class

  n_cases_total <- sum(subset_vector)
  n_cases_genPop <- sum(class_vec == "genPop" & subset_vector, na.rm = TRUE)
  n_cases_HCW   <- sum(class_vec == "HCW"    & subset_vector, na.rm = TRUE)

  ## outcome == TRUE corresponds to death in the current model
  outcome_vec <- tdf$outcome

  n_deaths_total  <- sum(outcome_vec & subset_vector, na.rm = TRUE)
  n_deaths_genPop <- sum(outcome_vec & class_vec == "genPop" &
                           subset_vector, na.rm = TRUE)
  n_deaths_HCW    <- sum(outcome_vec & class_vec == "HCW" &
                           subset_vector, na.rm = TRUE)

  ## CFRs (case fatality risks)
  cfr_overall <- if (n_cases_total > 0)  n_deaths_total  / n_cases_total  else NA_real_
  cfr_genPop  <- if (n_cases_genPop > 0) n_deaths_genPop / n_cases_genPop else NA_real_
  cfr_HCW     <- if (n_cases_HCW > 0)    n_deaths_HCW    / n_cases_HCW    else NA_real_

  ##--------------------------------------------------------------
  ## 3. Attack rates (require sim_info with population & hcw_total)
  ##--------------------------------------------------------------
  population <- NA_real_
  hcw_total  <- NA_real_

  if (!is.null(sim_info)) {
    if (!is.null(sim_info$population)) {
      population <- sim_info$population
    }
    if (!is.null(sim_info$hcw_total)) {
      hcw_total <- sim_info$hcw_total
    }
  }

  ## Overall attack rate
  attack_rate_overall <- if (is.finite(population) && population > 0) {
    n_cases_total / population
  } else NA_real_

  ## Approximate genPop population as N - HCWs if both are supplied
  genpop_pop <- if (is.finite(population) && is.finite(hcw_total)) {
    population - hcw_total
  } else NA_real_

  attack_rate_genPop <- if (is.finite(genpop_pop) && genpop_pop > 0) {
    n_cases_genPop / genpop_pop
  } else NA_real_

  ## HCW attack rate = proportion of HCW stock infected
  hcw_attack_rate <- if (is.finite(hcw_total) && hcw_total > 0) {
    n_cases_HCW / hcw_total
  } else NA_real_

  ## deaths per 1,000 population (overall)
  deaths_per_1000_pop <- if (is.finite(population) && population > 0) {
    n_deaths_total / population * 1000
  } else NA_real_

  ##--------------------------------------------------------------
  ## 4. Transmission setting breakdown (community / hospital / funeral)
  ##--------------------------------------------------------------
  if (!is.null(tdf$infection_location)) {
    setting_vec <- tdf$infection_location

    n_comm    <- sum(setting_vec == "community" & subset_vector, na.rm = TRUE)
    n_hosp    <- sum(setting_vec == "hospital"  & subset_vector, na.rm = TRUE)
    n_funeral <- sum(setting_vec == "funeral"   & subset_vector, na.rm = TRUE)

    n_with_setting <- n_comm + n_hosp + n_funeral

    prop_comm    <- if (n_with_setting > 0) n_comm    / n_with_setting else NA_real_
    prop_hosp    <- if (n_with_setting > 0) n_hosp    / n_with_setting else NA_real_
    prop_funeral <- if (n_with_setting > 0) n_funeral / n_with_setting else NA_real_

  } else {

    n_comm <- n_hosp <- n_funeral <- NA_real_
    prop_comm <- prop_hosp <- prop_funeral <- NA_real_
  }

  ##--------------------------------------------------------------
  ## 5. OBV PEP summary
  ##
  ## The gate counters are nested as sets per individual:
  ##   pre_eligible >= pre_treated >= pre_adherent
  ##                ^                ^
  ##                |                |--- "Policy A: treat all contacts" denominator + treated
  ##                |
  ##   post_eligible >= post_treated >= post_adherent >= prevented
  ##                 ^                 ^
  ##                 |                 |--- "Policy B: treat only PPE failures" denominator + treated
  ##                 |--- subset of pre_eligible that also survived PPE/quarantine thinning
  ##
  ## `n_obv_pep_prevented_deaths` <= `n_obv_pep_prevented`: the subset of prevented
  ## infections that would have died had they occurred. It is a deferred within-run
  ## counterfactual (resolved post-loop via the same outcome model as realised cases;
  ## see branching_process_main), so deaths averted by OBV can be read off a single run
  ## rather than differencing a separate no-OBV run (which diverges stochastically).
  ## Interpretation: this is a *direct* count -- the would-be death of each prevented
  ## index infection only -- and so a LOWER BOUND on the total deaths the programme
  ## averts, excluding deaths in the onward transmission chains those infections would
  ## have seeded.
  ##
  ## tdf-based cohort counters: realised HCW cases in the linelist who were
  ## eligible / treated / adherent (i.e. OBV did not prevent their infection).
  ## Only `n_obv_pep_breakthroughs` (= realised cases who received AND adhered
  ## to OBV) is a clinical "breakthrough" in the vaccine-failure sense; the
  ## other two are descriptive counts of the eligible/treated cohort that still
  ## became cases and let you decompose the failure modes:
  ##   eligible_cases - treated_cases  = HCW cases the coverage gap let through
  ##   treated_cases  - breakthroughs  = HCW cases the adherence gap let through
  ##   breakthroughs                   = HCW cases adequate OBV failed to prevent
  ##--------------------------------------------------------------
  obv_treated <- attr(tdf, "obv_pep_num_treated", exact = TRUE)
  if (is.null(obv_treated) && !is.null(sim_info$obv_pep_num_treated)) {
    obv_treated <- sim_info$obv_pep_num_treated
  }

  n_obv_pep_pre_eligible  <- if (!is.null(obv_treated)) obv_treated$pre_eligible  else NA_real_
  n_obv_pep_pre_treated   <- if (!is.null(obv_treated)) obv_treated$pre_treated   else NA_real_
  n_obv_pep_pre_adherent  <- if (!is.null(obv_treated)) obv_treated$pre_adherent  else NA_real_
  n_obv_pep_post_eligible <- if (!is.null(obv_treated)) obv_treated$post_eligible else NA_real_
  n_obv_pep_post_treated  <- if (!is.null(obv_treated)) obv_treated$post_treated  else NA_real_
  n_obv_pep_post_adherent <- if (!is.null(obv_treated)) obv_treated$post_adherent else NA_real_
  n_obv_pep_prevented     <- if (!is.null(obv_treated)) obv_treated$prevented     else NA_real_
  ## `prevented_deaths` is the deferred counterfactual: the subset of `prevented`
  ## infections that would have died had they occurred (see branching_process_main).
  ## Guard the field lookup so tdfs whose obv_pep_num_treated predates it stay NA.
  n_obv_pep_prevented_deaths <- if (!is.null(obv_treated) && !is.null(obv_treated$prevented_deaths)) {
    obv_treated$prevented_deaths
  } else NA_real_

  if (!is.null(tdf$obv_pep_eligible)) {
    n_obv_pep_eligible_cases <- sum(tdf$obv_pep_eligible & subset_vector, na.rm = TRUE)
    n_obv_pep_treated_cases  <- sum(tdf$obv_pep_received & subset_vector, na.rm = TRUE)
    n_obv_pep_breakthroughs  <- sum(tdf$obv_pep_adherent & subset_vector, na.rm = TRUE)
  } else {
    n_obv_pep_eligible_cases <- NA_real_
    n_obv_pep_treated_cases  <- NA_real_
    n_obv_pep_breakthroughs  <- NA_real_
  }

  prop_obv_pep_prevented_among_adherent <- if (is.finite(n_obv_pep_post_adherent) &&
                                               n_obv_pep_post_adherent > 0) {
    n_obv_pep_prevented / n_obv_pep_post_adherent
  } else NA_real_

  ##--------------------------------------------------------------
  ## 6. Contact tracing among realised cases
  ##
  ## `traced` records whether the CONTACT that produced this case was reached by
  ## tracing. Seed cases are never traced (no index case to trace them from).
  ##--------------------------------------------------------------
  if (!is.null(tdf$traced)) {
    n_cases_traced    <- sum(tdf$traced & subset_vector, na.rm = TRUE)
    prop_cases_traced <- if (n_cases_total > 0) n_cases_traced / n_cases_total else NA_real_
  } else {
    n_cases_traced <- NA_real_
    prop_cases_traced <- NA_real_
  }

  ## Risk-tier breakdown of realised cases. The tier mix among cases is risk-weighted
  ## relative to the tier mix among contacts, since higher-risk contacts are more likely
  ## to become cases.
  ## table() on a character column sorts alphabetically, which scrambles the tiers
  ## whenever they carry real names ("casual", "household", ...) rather than
  ## numbered ones. Since the whole point of these breakdowns is that the numbers
  ## should rise with risk, order the levels by the tier index instead.
  tier_levels <- function(level, category) {
    ok <- !is.na(level) & !is.na(category)
    if (!any(ok)) return(character(0))
    m <- unique(data.frame(l = level[ok], c = as.character(category[ok])))
    m$c[order(m$l)]
  }
  tier_table <- function(level, category) {
    lv <- tier_levels(level, category)
    if (length(lv) == 0L) return(table(factor(character(0))))
    table(factor(as.character(category), levels = lv), useNA = "no")
  }

  cases_by_risk_tier <- if (!is.null(tdf$contact_risk_category)) {
    tier_table(tdf$contact_risk_level[subset_vector],
               tdf$contact_risk_category[subset_vector])
  } else NULL

  ##--------------------------------------------------------------
  ## 7. Contact-level summary (requires the run's contact_log)
  ##
  ## This is the denominator a contact-tracing or prophylaxis programme actually
  ## faces: every exposure generated, not just the ones that became infections.
  ##--------------------------------------------------------------
  if (!is.null(contact_log) && nrow(contact_log) > 0) {
    n_contacts_total    <- nrow(contact_log)
    n_contacts_traced   <- sum(contact_log$traced, na.rm = TRUE)
    n_contacts_infected <- sum(contact_log$record_type == "infection", na.rm = TRUE)

    contacts_by_risk_tier <- tier_table(contact_log$contact_risk_level,
                                        contact_log$contact_risk_category)
    contacts_by_location  <- table(contact_log$infection_location, useNA = "no")

    ## Realised per-tier attack rate: of the contacts in each tier, what share became
    ## cases. Reflects the tier's relative risk after intervention thinning.
    infected_by_tier <- table(
      factor(contact_log$contact_risk_category[contact_log$record_type == "infection"],
             levels = names(contacts_by_risk_tier))
    )
    attack_rate_by_risk_tier <- as.numeric(infected_by_tier) / as.numeric(contacts_by_risk_tier)
    names(attack_rate_by_risk_tier) <- names(contacts_by_risk_tier)

    blocked_by_reason <- table(contact_log$blocked_by, useNA = "no")

    prop_contacts_traced <- n_contacts_traced / n_contacts_total
    contacts_per_case    <- if (n_cases_total > 0) n_contacts_total / n_cases_total else NA_real_
  } else {
    n_contacts_total <- n_contacts_traced <- n_contacts_infected <- NA_real_
    prop_contacts_traced <- contacts_per_case <- NA_real_
    contacts_by_risk_tier <- contacts_by_location <- NULL
    attack_rate_by_risk_tier <- blocked_by_reason <- NULL
  }

  ##--------------------------------------------------------------
  ## 8. Return a named list
  ##--------------------------------------------------------------
  out <- list(
    ## Outbreak timing
    outbreak_start_time      = first_infection_time,
    outbreak_end_time        = last_outcome_time,
    outbreak_duration_cont   = outbreak_duration_cont,
    outbreak_duration_days   = outbreak_duration_days,

    ##  case and death counts by class
    n_cases_total            = n_cases_total,
    n_cases_genPop           = n_cases_genPop,
    n_cases_HCW              = n_cases_HCW,
    n_deaths_total           = n_deaths_total,
    n_deaths_genPop          = n_deaths_genPop,
    n_deaths_HCW             = n_deaths_HCW,

    ## Setting breakdown
    n_comm                   = n_comm,
    n_hosp                   = n_hosp,
    n_funeral                = n_funeral,
    prop_comm                = prop_comm,
    prop_hosp                = prop_hosp,
    prop_funeral             = prop_funeral,

    ## CFRs
    cfr_overall              = cfr_overall,
    cfr_genPop               = cfr_genPop,
    cfr_HCW                  = cfr_HCW,

    ## Population + attack rates
    population               = population,
    hcw_total                = hcw_total,
    attack_rate_overall      = attack_rate_overall,
    attack_rate_genPop       = attack_rate_genPop,
    hcw_attack_rate          = hcw_attack_rate,
    deaths_per_1000_pop      = deaths_per_1000_pop,

    ## OBV PEP num-treated counters (see Step 5 above for set-nesting semantics)
    n_obv_pep_pre_eligible            = n_obv_pep_pre_eligible,
    n_obv_pep_pre_treated             = n_obv_pep_pre_treated,
    n_obv_pep_pre_adherent            = n_obv_pep_pre_adherent,
    n_obv_pep_post_eligible           = n_obv_pep_post_eligible,
    n_obv_pep_post_treated            = n_obv_pep_post_treated,
    n_obv_pep_post_adherent           = n_obv_pep_post_adherent,
    n_obv_pep_prevented               = n_obv_pep_prevented,
    ## Deferred counterfactual: would-be deaths among the `prevented` infections.
    n_obv_pep_prevented_deaths        = n_obv_pep_prevented_deaths,

    ## OBV PEP tdf-based cohort counters (HCW cases who became cases despite
    ## being in the eligible / treated / adherent cohort). Only
    ## `n_obv_pep_breakthroughs` (adherent recipients who still got infected)
    ## is a clinical "breakthrough" in the vaccine-failure sense.
    n_obv_pep_eligible_cases          = n_obv_pep_eligible_cases,
    n_obv_pep_treated_cases           = n_obv_pep_treated_cases,
    n_obv_pep_breakthroughs           = n_obv_pep_breakthroughs,
    prop_obv_pep_prevented_among_adherent = prop_obv_pep_prevented_among_adherent,

    ## Was this run censored by check_final_size? A censored final size measures the
    ## cap rather than transmission, so anything comparing final sizes must check this.
    hit_final_size_cap       = if (!is.null(sim_info$hit_final_size_cap)) sim_info$hit_final_size_cap else NA,
    stop_reason              = if (!is.null(sim_info$stop_reason)) sim_info$stop_reason else NA_character_,

    ## Contact tracing among realised cases
    n_cases_traced           = n_cases_traced,
    prop_cases_traced        = prop_cases_traced,
    cases_by_risk_tier       = cases_by_risk_tier,

    ## Contact-level counts (NA unless the run's contact_log was supplied)
    n_contacts_total         = n_contacts_total,
    n_contacts_traced        = n_contacts_traced,
    n_contacts_infected      = n_contacts_infected,
    prop_contacts_traced     = prop_contacts_traced,
    contacts_per_case        = contacts_per_case,
    contacts_by_risk_tier    = contacts_by_risk_tier,
    contacts_by_location     = contacts_by_location,
    attack_rate_by_risk_tier = attack_rate_by_risk_tier,
    contacts_blocked_by      = blocked_by_reason
  )

  return(out)
}

