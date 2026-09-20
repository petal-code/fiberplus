## Add in details on parameter arguments here
## Add in details about outputs etc here
## Basically all the required documentation for this to be
## functional inside a package

## Note: need to add descriptions of each of the function inputs, which are currently missing
## Note: currently the code has a mixture of distribution parameter inputs (e.g. generation time) and others
##       where there are actual parameters of the distribution inputs (e.g. Tg_shape_funeral and Tg_rate_funeral)
##       we should harmonise this at some point
## Note: we might make a function called "generate_seeding_case_attributes" or something that does everything that we do in step 2 currently
## Note: we should change the function so that 1s/0s of parent outcome are characters i.e. "death" / "recovery" explicitly
## Note: general note that we should be actively thinking about how to ensure we don't end up in weird edge cases where all
##       of our infections end up dying or recovering before they even need healthcare. that'll require us to be careful with how
##       we approach parameterising this (maybe we put a check in place??)
## Note: prob_hospitalised_hcw and prob_hospitalised_genPop - do they need to be made specific to the location of the infection as well?

#' Run a stochastic branching-process outbreak simulation
#'
#' Top-level entry point for fiber's filovirus branching-process model.
#' Iteratively generates offspring (community, hospital, and funeral
#' transmission) from active cases until the outbreak ends or the configured
#' final-size cap is reached. Time-varying scenario inputs (probabilities,
#' delay factors, IPC / ETU coverage, and mean-offspring transmissibility) can
#' be passed as scalars or functions of time produced by [make_time_varying()].
#'
#' @param mn_contacts_genPop Positive numeric or function(t). Mean of the Negative Binomial
#'   *contact* distribution for general-population (genPop) parents; resolved at the parent's
#'   infection time. Scalar or a function of absolute calendar time (e.g. [make_time_varying()]).
#' @param overdisp_contacts_genPop Positive numeric. Negative Binomial size (overdispersion) of the
#'   genPop contact distribution. Does not affect the implied R0, so it dials superspreading
#'   independently of the calibration.
#' @param baseline_risk_genPop Numeric in `[0, 1]` or function(t). Per-contact transmission
#'   probability for the *reference* risk tier on the genPop route; other tiers scale it by their
#'   relative risk. Supply this or `r0_target`, not both.
#' @param Tg_shape_genPop Positive numeric. Shape of the Gamma generation-time distribution for
#'   genPop parents.
#' @param Tg_rate_genPop Positive numeric. Rate of the Gamma generation-time distribution for genPop
#'   parents (mean generation time = shape / rate).
#' @param mn_contacts_hcw Positive numeric or function(t). Mean Negative Binomial contact
#'   distribution for healthcare-worker (HCW) parents; resolved at the parent's infection time.
#' @param overdisp_contacts_hcw Positive numeric. Negative Binomial size (overdispersion) of the HCW
#'   contact distribution.
#' @param baseline_risk_hcw Numeric in `[0, 1]` or function(t). Per-contact transmission probability
#'   for the reference risk tier on the HCW route. When `r0_target` is used this defaults to the
#'   solved genPop baseline risk.
#' @param Tg_shape_hcw Positive numeric. Shape of the Gamma generation-time distribution for HCW
#'   parents.
#' @param Tg_rate_hcw Positive numeric. Rate of the Gamma generation-time distribution for HCW
#'   parents.
#' @param mn_contacts_funeral Positive numeric or function(t). Mean Negative Binomial number of
#'   contacts at a funeral; resolved at the parent's death (outcome) time.
#' @param overdisp_contacts_funeral Positive numeric. Negative Binomial size (overdispersion) of the
#'   funeral contact distribution.
#' @param baseline_risk_funeral Numeric in `[0, 1]` or function(t). Per-contact transmission
#'   probability for the reference risk tier at a funeral. Supply this or `r0_target`, not both.
#' @param contact_risk A [make_contact_risk()] structure (or a named list of its arguments) giving
#'   the risk tier fractions, relative risks and per-tier contact-tracing probabilities. Used by
#'   every route unless that route overrides it. Defaults to five flat tiers with no tracing.
#' @param contact_risk_genPop,contact_risk_hcw,contact_risk_funeral Optional per-route risk
#'   structures. `NULL` (the default) inherits `contact_risk`. Typically genPop and HCW share the
#'   inherited structure while funerals get their own.
#' @param r0_target Optional positive numeric. If supplied, the genPop and funeral baseline risks
#'   are solved from this target R0 via [solve_baseline_risk_for_r0()] instead of being given
#'   directly, using `r0_prop_funeral` to split transmission between the two routes. The solved
#'   values are returned in `sim_info`.
#' @param r0_prop_funeral Numeric in `[0, 1]`. Share of `r0_target` attributed to the funeral route.
#'   Only used when `r0_target` is supplied. Defaults to 0.
#' @param r0_solve_n,r0_solve_seed Monte-Carlo settings for the R0 inversion.
#' @param trace_coverage Numeric in `[0, 1]` or function(t). Programme-level contact tracing
#'   coverage, multiplying each risk tier's `trace_prob`. Defaults to 0 (no tracing).
#' @param onset_to_hospitalisation_traced Non-negative numeric or function(t), or NULL. Flat
#'   onset-to-admission delay, in days, for cases that were traced. It *caps* rather than replaces
#'   each traced case's own drawn delay, so tracing can only bring an admission forward, never
#'   push it back. `NULL` (default) means tracing does not change admission timing. A warning is
#'   raised if this is not comfortably below the untraced delay distribution, since a value above
#'   it would silently do nothing.
#' @param prob_hospitalised_traced Numeric in `[0, 1]` or function(t), or NULL. The absolute
#'   P(hospitalised | symptomatic) for traced cases, replacing the untraced value outright. Use this
#'   when a scenario says "traced cases are hospitalised with probability 0.9" — a multiplier cannot
#'   pin that down when the untraced probability is itself time-varying. `NULL` (default) means no
#'   effect. Takes precedence over `prob_hospitalised_multiplier_traced`; supplying both is an error.
#' @param prob_hospitalised_multiplier_traced Positive numeric or function(t). Multiplier on
#'   P(hospitalised | symptomatic) for traced cases, capped at 1. Defaults to 1 (no effect). Note
#'   that faster admission raises the *realised* hospitalisation rate on its own, independently of
#'   this multiplier, because admission is more likely to beat the community outcome.
#' @param presymptomatic_transmission Logical scalar. `TRUE` (default) allows contacts to occur
#'   before the infector develops symptoms — the model's natural behaviour, since contact times and
#'   incubation periods are drawn independently. `FALSE` truncates each parent's contact times to
#'   start at the end of their incubation period, removing presymptomatic transmission entirely.
#'   Applies to every parent, symptomatic or not. Note this lengthens the realised generation time,
#'   because the generation-time distribution is being conditioned rather than reshaped.
#' @param check_presymptomatic Logical scalar. If `TRUE` (default) and presymptomatic transmission
#'   is allowed, estimate its share once at the start of the run and warn when it exceeds
#'   `presymptomatic_warn_threshold`. The estimate uses its own RNG draws and restores the random
#'   seed afterwards, so it never perturbs the simulated trajectory. Set `FALSE` to skip the cost
#'   in large calibration runs. See [approx_presymptomatic_transmission()].
#' @param presymptomatic_warn_threshold Numeric in `[0, 1]`. Presymptomatic share above which
#'   `check_presymptomatic` warns. Defaults to 0.1.
#' @param return_contact_log Logical scalar. `TRUE` (default) returns the full contact log. `FALSE`
#'   skips building it per parent and accumulating it, and returns a 0-row frame. The log is roughly
#'   one row per contact, so at `check_final_size = 30000` with ~15 contacts per case it is around
#'   450,000 rows (~47 MB) **per run** — set this `FALSE` in calibration loops, where thousands of
#'   runs would otherwise hold that much each.
#' @param quiet Logical scalar. `FALSE` (default) emits a `message()` at the end of each run saying
#'   why the simulation stopped, and warns when it was censored by `check_final_size`. Set `TRUE` to
#'   silence it. A censored run's final size measures the cap, not transmission, so any analysis of
#'   final size must check `sim_info$hit_final_size_cap` first.
#' @param Tg_shape_funeral Positive numeric. Shape of the Gamma outcome-to-funeral-infection delay
#'   distribution.
#' @param Tg_rate_funeral Positive numeric. Rate of the Gamma funeral-delay distribution (mean delay
#'   = shape / rate).
#' @param incubation_period Function(n). Random generator returning n incubation-period draws
#'   (infection to symptom onset).
#' @param onset_to_hospitalisation Function(n). Random generator returning n
#'   onset-to-hospitalisation delay draws.
#' @param hospitalisation_delay_factor Positive numeric or function(t). Multiplier applied to
#'   `onset_to_hospitalisation` draws; may vary with absolute calendar time. Defaults to 1.
#' @param onset_to_death Function(n). Random generator returning n onset-to-death delay draws.
#' @param onset_to_recovery Function(n). Random generator returning n onset-to-recovery delay draws.
#' @param hospitalisation_to_death Function(n). Random generator returning n
#'   hospitalisation-to-death delay draws.
#' @param hospitalisation_to_recovery Function(n). Random generator returning n
#'   hospitalisation-to-recovery delay draws.
#' @param prob_symptomatic Numeric in `[0, 1]`. Probability an infection is symptomatic.
#' @param prob_hospitalised_hcw Numeric in `[0, 1]` or function(t). P(hospitalised | symptomatic) for
#'   HCWs; resolved at symptom-onset time.
#' @param prob_hospitalised_genPop Numeric in `[0, 1]` or function(t). P(hospitalised | symptomatic)
#'   for genPop; resolved at symptom-onset time.
#' @param prob_death_comm Numeric in `[0, 1]`. P(death | symptomatic) for cases managed in the
#'   community.
#' @param prob_death_hosp Numeric in `[0, 1]`. P(death | symptomatic) for hospitalised cases; must be
#'   no greater than `prob_death_comm`.
#' @param prob_hcw_cond_genPop_comm Numeric in `[0, 1]`. Probability a community-located infection from
#'   a genPop parent is an HCW.
#' @param prob_hcw_cond_genPop_hospital Numeric in `[0, 1]`. Probability a hospital-located infection
#'   from a genPop parent is an HCW.
#' @param prob_hcw_cond_hcw_comm Numeric in `[0, 1]`. Probability a community-located infection from an
#'   HCW parent is an HCW.
#' @param prob_hcw_cond_hcw_hospital Numeric in `[0, 1]`. Probability a hospital-located infection from
#'   an HCW parent is an HCW.
#' @param prob_hospital_cond_hcw_preAdm Numeric in `[0, 1]`. Probability that an infection generated by
#'   an HCW parent before their own admission occurs in the hospital (while still working).
#' @param ppe_efficacy_hcw Numeric in `[0, 1]` or function(t). Per-layer efficacy of PPE/IPC at
#'   reducing hospital transmission; resolved at each candidate hospital transmission time.
#' @param hospital_quarantine_efficacy Optional numeric in `[0, 1]` or function(t). Direct
#'   post-admission hospital quarantine/ETU efficacy. If NULL, derived from `prop_etu`, `ipc_helper`
#'   and `etu_efficacy_baseline`.
#' @param prop_etu Numeric in `[0, 1]` or function(t). Proportion of hospitalised cases in ETU/ETC
#'   care; used when `hospital_quarantine_efficacy` is NULL.
#' @param ipc_helper Numeric in `[0, 1]` or function(t). IPC/response-maturity proxy; used when
#'   `hospital_quarantine_efficacy` is NULL.
#' @param etu_efficacy_baseline Numeric in `[0, 1]`. Baseline ETU/ETC transmission-blocking efficacy
#'   before IPC-maturity adjustment; used when `hospital_quarantine_efficacy` is NULL.
#' @param obv_pep_enabled Logical scalar. If TRUE, apply the obeldesivir (OBV) PEP
#'   infection-prevention gate.
#' @param obv_pep_coverage Numeric in `[0, 1]` or function(t). Probability an eligible candidate
#'   receives OBV PEP.
#' @param obv_pep_adherence Numeric in `[0, 1]` or function(t). Probability a received OBV course is
#'   adhered to.
#' @param obv_pep_dpc Non-negative numeric or function(t). Days post challenge/exposure to first
#'   dose, evaluated at each candidate's own absolute infection time. With
#'   `obv_pep_dpc_sd = NULL` this is each recipient's DPC exactly; with an sd it is the
#'   *mean* DPC, so the average treatment delay can vary over calendar time while the
#'   efficacy(DPC) relationship stays fixed.
#' @param obv_pep_dpc_sd NULL or a single positive numeric. NULL (default) keeps DPC
#'   deterministic at `obv_pep_dpc` (bit-for-bit identical to pre-feature runs). If supplied,
#'   it is the standard deviation of the per-recipient DPC: each recipient draws an independent
#'   DPC from a Gamma reparameterised from (mean, sd) as
#'   `Gamma(shape = obv_pep_dpc(t)^2 / obv_pep_dpc_sd^2, scale = obv_pep_dpc_sd^2 / obv_pep_dpc(t))`
#'   -- mean `obv_pep_dpc(t)`, variance `obv_pep_dpc_sd^2` (`CV = obv_pep_dpc_sd / obv_pep_dpc(t)`)
#'   -- giving individual variation in how quickly the drug is received post-exposure. The spread
#'   is a fixed absolute sd even when the mean varies over time (a fixed shape would instead hold
#'   the CV constant).
#' @param obv_pep_efficacy NULL, numeric in `[0, 1]`, or function(dpc). Selects the efficacy model.
#'   NULL or a scalar use the built-in [obv_pep_efficacy_from_dpc()] curve, which is **flat** by
#'   default (constant efficacy at every DPC); a *scalar* is taken as that constant (the curve's
#'   `E0`). DPC-dependent decay is opt-in via `obv_pep_efficacy_args` (`shape = "logistic"`). A
#'   function(dpc) is used as-is (e.g. `function(dpc) rep(0.5, length(dpc))`).
#' @param obv_pep_efficacy_args NULL or a named list of overrides for the built-in efficacy curve
#'   [obv_pep_efficacy_from_dpc()] -- any of `shape`, `E0`, `d50`, `k`, `dpc_zero`, `max_dpc`. Applies
#'   to the built-in curve only (when `obv_pep_efficacy` is NULL or a scalar); e.g.
#'   `obv_pep_efficacy_args = list(shape = "logistic", d50 = 4)` switches on the logistic decay.
#'   Errors if combined with a function `obv_pep_efficacy`, if it names `E0` while `obv_pep_efficacy`
#'   is a scalar (E0 then comes from `obv_pep_efficacy`), or given an unknown name.
#' @param obv_pep_target_class Character vector. Offspring classes eligible for OBV PEP (default
#'   "HCW").
#' @param obv_pep_target_locations Character vector. Exposure settings eligible for OBV PEP (default
#'   "hospital").
#' @param p_unsafe_funeral_comm_hcw Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a community death of an HCW.
#' @param p_unsafe_funeral_hosp_hcw Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a hospital death of an HCW.
#' @param p_unsafe_funeral_comm_genPop Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a community death of a genPop case.
#' @param p_unsafe_funeral_hosp_genPop Numeric in `[0, 1]` or function(t). Probability of an unsafe
#'   funeral after a hospital death of a genPop case.
#' @param safe_funeral_efficacy Numeric in `[0, 1]`. Efficacy of a safe burial at preventing funeral
#'   transmission (1 = fully blocking).
#' @param prob_hcw_cond_funeral_hcw Numeric in `[0, 1]`. Probability a funeral infection from an HCW
#'   parent is an HCW.
#' @param prob_hcw_cond_funeral_genPop Numeric in `[0, 1]`. Probability a funeral infection from a
#'   genPop parent is an HCW.
#' @param population Positive integer. Total population size.
#' @param hcw_per_capita Positive numeric. HCWs per capita; total HCWs is
#'   `round(hcw_per_capita * population)`.
#' @param check_final_size Positive integer. Final-size cap; the simulation stops once this many
#'   cases have been generated.
#' @param initial_immune Non-negative integer. Number initially immune (removed from the susceptible
#'   pool). Defaults to 0.
#' @param seeding_cases Positive integer. Number of initial genPop community seeding infections.
#' @param susceptible_deplete Logical scalar. Placeholder for susceptible-depletion behaviour (not
#'   yet implemented). Defaults to FALSE.
#' @param seed Optional integer. RNG seed for reproducibility.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{`tdf`}{The simulated transmission tree: one row per realised
#'       infection, ordered by absolute infection time. Carries attributes
#'       `hcw_total`, `hcw_infected`, `hcw_remaining`, and `obv_pep_num_treated`.
#'       Each case also records the risk tier of the contact that produced it
#'       (`contact_risk_level`, `contact_risk_category`), whether that contact was
#'       traced (`traced`).}
#'     \item{`contact_log`}{Every contact generated over the run, one row each --
#'       including the contacts that never became infections. Columns: `parent`,
#'       `case_id` (the `tdf` id for contacts that became cases, `NA` otherwise),
#'       `record_type` (`"contact"` or `"infection"`), `class`,
#'       `infection_location`, `time_contact_relative`, `time_contact_absolute`,
#'       `contact_risk_level`, `contact_risk_category`, `relative_risk`,
#'       `transmission_prob` (the realised per-contact probability), `traced`, and
#'       `blocked_by` -- `NA` for realised infections, else `"no_transmission"`
#'       (the contact simply did not transmit), the route's intervention layer, or
#'       `"obv_pep"`. This is the denominator for anything contact-tracing or
#'       prophylaxis related.}
#'     \item{`prevented_completed`}{A data frame of the infections the OBV PEP gate
#'       prevented -- the averted index infections only, not their averted onward
#'       chains -- each replayed through the same outcome model as realised cases to
#'       give its counterfactual natural history. `time_infection_absolute` is the
#'       would-be calendar infection time of each averted infection, and `outcome`
#'       its would-be death/recovery status (summing `outcome` reproduces
#'       `obv_pep_num_treated$prevented_deaths`). `NULL` when the gate prevented
#'       nothing (including when `obv_pep_enabled = FALSE`). The replay uses a
#'       zero-time dummy parent, so `parent` and `generation` are `NA` and
#'       `time_infection_relative` equals `time_infection_absolute`.}
#'     \item{`sim_info`}{Scalar run metadata: `population`, `hcw_per_capita`,
#'       `hcw_total`, `seed`, `obv_pep_enabled`, and the `obv_pep_num_treated`
#'       counters.}
#'   }
#'
#' @export
branching_process_main <- function(

  ## Transmission: contacts first, then a per-contact risk-tier transmission draw.
  ## Each route has its own contact distribution, baseline per-contact risk and risk
  ## structure. Baseline risks can be supplied directly or solved from `r0_target`.
  mn_contacts_genPop = NULL,                # scalar or function(t): mean CONTACT distribution for genPop (resolved at parent infection time)
  overdisp_contacts_genPop = NULL,          # overdispersion of the contact distribution for genPop
  baseline_risk_genPop = NULL,              # scalar or function(t): per-contact transmission prob for the reference risk tier
  Tg_shape_genPop = NULL,                   # gamma shape parameter for Tg distribution for general population
  Tg_rate_genPop = NULL,                    # gamma rate parameter for Tg distribution for general population
  mn_contacts_hcw = NULL,                   # scalar or function(t): mean CONTACT distribution for HCWs (resolved at parent infection time)
  overdisp_contacts_hcw = NULL,             # overdispersion of the contact distribution for HCWs
  baseline_risk_hcw = NULL,                 # scalar or function(t): per-contact transmission prob for the reference risk tier
  Tg_shape_hcw = NULL,                      # gamma shape parameter for Tg distribution for HCWs
  Tg_rate_hcw = NULL,                       # gamma rate parameter for Tg distribution for HCWs
  mn_contacts_funeral = NULL,               # scalar or function(t): mean CONTACT distribution at a funeral (resolved at parent death time)
  overdisp_contacts_funeral = NULL,         # overdispersion of the funeral contact distribution
  baseline_risk_funeral = NULL,             # scalar or function(t): per-contact transmission prob for the reference risk tier
  Tg_shape_funeral = NULL,                  # gamma shape parameter for Tg distribution at funerals ### have high shape, high rate to get low variance ##
  Tg_rate_funeral = NULL,                   # gamma rate parameter for Tg distribution at funerals

  ## Contact risk structure: tier fractions, relative risks and per-tier trace probabilities.
  ## Per-route arguments default to the shared `contact_risk`, so genPop and HCW can share one
  ## structure while funerals use a different one.
  contact_risk = NULL,                      # fiber_contact_risk structure shared across routes (default: 5 flat tiers, no tracing)
  contact_risk_genPop = NULL,               # NULL = inherit contact_risk
  contact_risk_hcw = NULL,                  # NULL = inherit contact_risk
  contact_risk_funeral = NULL,              # NULL = inherit contact_risk

  ## Optional R0-anchored calibration: solve the genPop and funeral baseline risks from a
  ## target R0 rather than supplying them directly. See solve_baseline_risk_for_r0().
  r0_target = NULL,                         # positive scalar: target R0 at t = 0
  r0_prop_funeral = 0,                      # share of r0_target coming from the funeral route
  r0_solve_n = 50000,                       # Monte-Carlo draws for the R0 inversion
  r0_solve_seed = NULL,                     # optional seed for the R0 inversion

  ## Contact tracing. The probability of being traced is tier-specific (it lives in the
  ## contact risk structure), so the risk tiers drive these NPIs. A traced case is admitted
  ## sooner, and optionally more often.
  trace_coverage = 0,                       # scalar/function(t): programme-level tracing coverage
  onset_to_hospitalisation_traced = NULL,   # scalar/function(t): flat onset-to-admission delay for traced cases (caps their own); NULL = no effect
  prob_hospitalised_traced = NULL,          # scalar/function(t): ABSOLUTE P(hospitalised | symptomatic) for traced cases; NULL = no effect
  prob_hospitalised_multiplier_traced = 1,  # scalar/function(t): multiplier on P(hospitalised | symptomatic) for traced cases (ignored if the absolute is set)

  ## Presymptomatic transmission. TRUE (the model's natural behaviour) lets contacts happen
  ## before the infector's symptom onset; FALSE truncates contact times to start at onset.
  presymptomatic_transmission = TRUE,
  check_presymptomatic = TRUE,              # estimate the presymptomatic share once and warn if large
  presymptomatic_warn_threshold = 0.1,      # share above which check_presymptomatic warns

  ## Output control
  return_contact_log = TRUE,                # FALSE skips building and accumulating the contact log
  quiet = FALSE,                            # TRUE suppresses the end-of-run stop-reason message

  ## Natural history
  incubation_period,              # DESCRIPTION HERE
  onset_to_hospitalisation,    # DESCRIPTION HERE
  hospitalisation_delay_factor = 1.0,   # scalar or function(t): multiplier on onset_to_hospitalisation draws
  onset_to_death,
  onset_to_recovery,
  hospitalisation_to_death,              # Note: Jacob to look up whether the time -> death is the same typically as time -> recovery (or do they need to be different)
  hospitalisation_to_recovery,           # Note: Jacob to look up whether the time -> death is the same typically as time -> recovery (or do they need to be different)

  # Disease severity and healthcare seeking. Hospitalisation and death
  # probabilities are CONDITIONAL on the offspring being symptomatic.
  prob_symptomatic = NULL,           # P(symptomatic | infected)
  prob_hospitalised_hcw = NULL,      # scalar or function(t): P(hospitalised | symptomatic) for HCWs
  prob_hospitalised_genPop = NULL,   # scalar or function(t): P(hospitalised | symptomatic) for genPop
  prob_death_comm = NULL,            # P(die | symptomatic, community)
  prob_death_hosp = NULL,            # P(die | symptomatic, hospitalised); must be <= prob_death_comm

  ## Probabilities for genPop infecting either genPop or HCWs, depending on the setting
  prob_hcw_cond_genPop_comm = NULL,         # prob that a community-located infection generated by genPop is a HCW
  prob_hcw_cond_genPop_hospital = NULL,     # prob that a hospital-located infection generated by genPop is a HCW

  ## Probabilities for HCW infecting either genPop or HCWs, depending on the setting
  prob_hcw_cond_hcw_comm = NULL,            # prob that a community-located infection generated by HCW is a HCW
  prob_hcw_cond_hcw_hospital = NULL,        # prob that a hospital-located infection generated by HCW is a HCW

  ## Setting model for HCWs
  prob_hospital_cond_hcw_preAdm = NULL,     # probability that an infection generated prior to parent hospitaliation occurs in the hospital (whilst HCW is working)
  ppe_coverage_hcw = NULL,                  # scalar or function(t): coverage/probability that a relevant HCW has PPE (time-varying coverage lever)
  ppe_efficacy = NULL,                      # scalar: efficacy of PPE at preventing infection conditional on having it
  prop_etu = NULL,                          # scalar/function(t): proportion of hospitalised cases in ETU/ETC care (time-varying coverage lever)
  etu_efficacy = NULL,                      # scalar: post-admission quarantine efficacy for ETU/ETC care
  general_hospital_quarantine_efficacy = NULL,  # scalar: post-admission quarantine efficacy for general (non-ETU) hospital care

  ## Obeldesivir PEP. The gate is applied around the Swiss-cheese thinning step
  ## in each offspring function: treatment status (received, adherent, DPC) is
  ## assigned to all pre-thinning eligible candidates, and efficacy is applied
  ## only to candidates that also survive PPE/quarantine. See apply_obv_pep_gate().
  ## Per-individual treatment-status decisions compose multiplicatively with
  ## the existing PPE source / PPE receiver / hospital quarantine layers in the
  ## offspring functions; the NB overdispersion of HCW-only offspring counts is
  ## not preserved exactly under this thinning.
  ##
  ## Reproducibility caveat: when obv_pep_enabled = TRUE, the gate consumes
  ## additional rbinom() draws (coverage, adherence, prevention) per offspring
  ## call. This means OBV-enabled and OBV-disabled simulations with the same
  ## `seed` will NOT produce identical genPop / hospital / funeral trajectories
  ## -- they diverge as soon as the first OBV draw fires. For paired
  ## counterfactual comparisons use many replicates with different seeds and
  ## compare distributions (see Scenario 2 in the OBV verification script).
  obv_pep_enabled = FALSE,                   # logical: apply OBV infection-prevention gate
  obv_pep_coverage = 0,                  # scalar/function(t): probability eligible candidate receives OBV
  obv_pep_adherence = 1,                     # scalar/function(t): probability received course is effectively adhered to
  obv_pep_dpc = 1,                           # scalar/function(t): days post challenge/exposure to first dose (mean DPC when an sd is set)
  obv_pep_dpc_sd = NULL,                     # NULL = deterministic DPC; positive scalar = per-recipient Gamma(mean = obv_pep_dpc(t), sd) draw
  obv_pep_efficacy = NULL,                   # NULL/function(dpc)/scalar: efficacy; NULL uses obv_pep_efficacy_from_dpc()
  obv_pep_efficacy_args = NULL,              # NULL or named list of overrides for obv_pep_efficacy_from_dpc() (shape/E0/d50/k/dpc_zero/max_dpc); used when obv_pep_efficacy is NULL or a scalar (errors if a function)
  obv_pep_target_class = "HCW",              # character vector: offspring classes eligible for OBV PEP
  obv_pep_target_locations = "hospital",     # character vector: exposure settings eligible for OBV PEP

  ## Funeral occurrence
  p_unsafe_funeral_comm_hcw = NULL, ## scalar or function(t): probability of unsafe funeral after a community death, HCW
  p_unsafe_funeral_hosp_hcw = NULL, ## scalar or function(t): probability of unsafe funeral after a hospital death, HCW
  p_unsafe_funeral_comm_genPop = NULL, ## scalar or function(t): probability of unsafe funeral after a community death, genPop
  p_unsafe_funeral_hosp_genPop = NULL, ## scalar or function(t): probability of unsafe funeral after a hospital death, genPop
  safe_funeral_efficacy = NULL, ## efficacy of a safe burial in reducing transmission in a funeral setting

  ## HCW vs genPop at funeral
  prob_hcw_cond_funeral_hcw = NULL, ### probability that the unsafe funeral infector infects a HCW
  prob_hcw_cond_funeral_genPop = NULL, ## DESCRIPTION NEEDED HERE

  ## Misc
  population,
  hcw_per_capita = 10,
  check_final_size,
  initial_immune = 0,
  seeding_cases,
  susceptible_deplete = FALSE,  ## note - still need to add code around this as functionality
                                ## envisaging this will adapt mn_offspring to reflect susceptible depletion
  seed = NULL

) {

  ##################################################################
  ### Step 1: Set up everything we need for the simulation
  ##################################################################
  # Set seed for reproducibility
  set.seed(seed)

  ## Local helpers for resolving parameters that can be either scalars or
  ## functions of calendar time. These mirror the logic in complete_offspring_info
  ## because some time-varying parameters are needed before offspring are created.
  resolve_probability <- function(param, t, param_name) {
    value <- resolve_time_varying(param = param, t = t, param_name = param_name)
    if (any(value < 0 | value > 1)) {
      stop(sprintf("`%s` must resolve to value(s) in [0, 1].", param_name), call. = FALSE)
    }
    value
  }

  ## Post-admission hospital quarantine efficacy is always derived inside the
  ## offspring functions as a prop_etu(t) mixture of the ETU and general-hospital
  ## quarantine efficacies. Require the full set up front so missing inputs fail
  ## early with a clear message rather than failing downstream.
  if (is.null(prop_etu) ||
      is.null(etu_efficacy) ||
      is.null(general_hospital_quarantine_efficacy)) {
    stop(
      "Supply all of `prop_etu`, `etu_efficacy`, and `general_hospital_quarantine_efficacy`.",
      call. = FALSE
    )
  }

  ## Fixed scalar efficacies in [0, 1]. The time-varying response enters through the
  ## coverage levers (`ppe_coverage_hcw`, `prop_etu`), not the efficacies, which are
  ## deliberately scalar-only. The ETU and general-hospital efficacies are
  ## independently togglable; no ordering between them is enforced.
  validate_scalar_probability <- function(param, param_name) {
    if (!is.numeric(param) ||
        length(param) != 1L ||
        is.na(param) ||
        param < 0 ||
        param > 1) {
      stop(sprintf("`%s` must be a single numeric value between 0 and 1.", param_name),
           call. = FALSE)
    }
  }

  validate_scalar_probability(ppe_efficacy, "ppe_efficacy")
  validate_scalar_probability(etu_efficacy, "etu_efficacy")
  validate_scalar_probability(general_hospital_quarantine_efficacy,
                              "general_hospital_quarantine_efficacy")

  ##################################################################
  ### Step 1a: Resolve the contact risk structures and, optionally,
  ### solve the baseline per-contact risks from a target R0.
  ###
  ### Each route gets its own structure, falling back to the shared
  ### `contact_risk` when not overridden -- so genPop and HCW can share one
  ### structure while funerals use a steeper one.
  ##################################################################
  contact_risk_shared  <- as_contact_risk(contact_risk, NULL, "contact_risk")
  risk_genPop  <- as_contact_risk(contact_risk_genPop,  contact_risk_shared, "contact_risk_genPop")
  risk_hcw     <- as_contact_risk(contact_risk_hcw,     contact_risk_shared, "contact_risk_hcw")
  risk_funeral <- as_contact_risk(contact_risk_funeral, contact_risk_shared, "contact_risk_funeral")

  r0_solution <- NULL
  if (!is.null(r0_target)) {
    if (!is.null(baseline_risk_genPop) || !is.null(baseline_risk_funeral)) {
      stop("Supply either `r0_target` (to solve the baseline risks) or the `baseline_risk_*` arguments directly, not both.",
           call. = FALSE)
    }
    ## The inversion needs the same inputs the simulation uses, so hand it the
    ## already-resolved structures alongside the natural-history parameters.
    r0_args <- list(
      mn_contacts_genPop           = mn_contacts_genPop,
      mn_contacts_funeral          = mn_contacts_funeral,
      contact_risk_genPop          = risk_genPop,
      contact_risk_funeral         = risk_funeral,
      Tg_shape_genPop              = Tg_shape_genPop,
      Tg_rate_genPop               = Tg_rate_genPop,
      incubation_period            = incubation_period,
      onset_to_hospitalisation     = onset_to_hospitalisation,
      hospitalisation_delay_factor = hospitalisation_delay_factor,
      onset_to_death               = onset_to_death,
      onset_to_recovery            = onset_to_recovery,
      hospitalisation_to_death     = hospitalisation_to_death,
      hospitalisation_to_recovery  = hospitalisation_to_recovery,
      prob_symptomatic             = prob_symptomatic,
      prob_hospitalised_genPop     = prob_hospitalised_genPop,
      prob_death_comm              = prob_death_comm,
      prob_death_hosp              = prob_death_hosp,
      prop_etu                     = prop_etu,
      etu_efficacy                 = etu_efficacy,
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy        = safe_funeral_efficacy,
      p_unsafe_funeral_comm_genPop = p_unsafe_funeral_comm_genPop,
      p_unsafe_funeral_hosp_genPop = p_unsafe_funeral_hosp_genPop,
      trace_coverage               = trace_coverage,
      prob_hospitalised_multiplier_traced = prob_hospitalised_multiplier_traced,
      prob_hospitalised_traced            = prob_hospitalised_traced,
      onset_to_hospitalisation_traced     = onset_to_hospitalisation_traced,
      presymptomatic_transmission         = presymptomatic_transmission
    )
    r0_solution <- solve_baseline_risk_for_r0(
      R0   = r0_target,
      args = r0_args,
      proportion_transmission_from_funerals = r0_prop_funeral,
      n    = r0_solve_n,
      seed = r0_solve_seed
    )
    baseline_risk_genPop  <- r0_solution$baseline_risk_genPop_required
    baseline_risk_funeral <- r0_solution$baseline_risk_funeral_required
  }

  ## HCW parents share the genPop per-contact risk unless given their own. Per-contact
  ## transmission risk is a property of the exposure rather than the infector's
  ## occupation, and the single-type R0 inversion has no separate HCW term to solve
  ## against. The HCW contact distribution stays explicit: how many contacts a
  ## healthcare worker has is genuinely a different question.
  if (is.null(baseline_risk_hcw)) baseline_risk_hcw <- baseline_risk_genPop

  for (nm in c("mn_contacts_genPop", "mn_contacts_hcw", "mn_contacts_funeral",
               "overdisp_contacts_genPop", "overdisp_contacts_hcw", "overdisp_contacts_funeral",
               "baseline_risk_genPop", "baseline_risk_hcw", "baseline_risk_funeral")) {
    if (is.null(get(nm, inherits = FALSE))) {
      stop(sprintf("`%s` is required. Supply the baseline risks directly, or set `r0_target` to solve them.", nm),
           call. = FALSE)
    }
  }

  if (!is.logical(obv_pep_enabled) || length(obv_pep_enabled) != 1L || is.na(obv_pep_enabled)) {
    stop("`obv_pep_enabled` must be a single logical value.", call. = FALSE)
  }
  if (!is.character(obv_pep_target_class) || length(obv_pep_target_class) < 1L) {
    stop("`obv_pep_target_class` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(obv_pep_target_locations) || length(obv_pep_target_locations) < 1L) {
    stop("`obv_pep_target_locations` must be a non-empty character vector.", call. = FALSE)
  }
  if (isTRUE(obv_pep_enabled) &&
      is.numeric(obv_pep_coverage) && length(obv_pep_coverage) == 1L &&
      !is.na(obv_pep_coverage) && obv_pep_coverage == 0) {
    warning("`obv_pep_enabled = TRUE` but `obv_pep_coverage = 0`; the OBV PEP gate will be a no-op.",
            call. = FALSE)
  }

  ## OBV PEP per-call accumulator: 7 gate counters, see empty_obv_pep_num_treated().
  obv_num_treated <- empty_obv_pep_num_treated()
  ## Collected per-call snapshots of infections OBV prevented (no RNG drawn in
  ## the loop). Their counterfactual would-be deaths are resolved once, after the
  ## loop, to populate obv_num_treated$prevented_deaths without perturbing the
  ## simulated trajectory's RNG stream.
  obv_prevented_info_list <- list()
  ## Per-parent contact logs (every contact generated, whether or not it transmitted).
  ## Concatenated once after the loop into the returned `contact_log`.
  contact_log_list <- list()
  ##################################################################
  ### Step 1b: Upfront sanity checks on time-varying parameters
  ###
  ### Sample each probability parameter on a grid built from the
  ### make_time_varying breakpoints (where supplied) plus midpoints, or
  ### a default 0..365 grid otherwise. This catches curves that resolve
  ### outside [0, 1] somewhere in the simulation horizon before the
  ### simulation starts, rather than mid-run with a cryptic error.
  ###
  ### Also enforces the cross-parameter constraint that
  ###     prob_death_hosp <= prob_death_comm
  ### so the second-chance ratio prob_death_hosp / prob_death_comm
  ### stays in [0, 1] (see complete_offspring_info Step 3.2).
  ##################################################################
  sanity_params <- list(
    prob_hospitalised_hcw         = prob_hospitalised_hcw,
    prob_hospitalised_genPop      = prob_hospitalised_genPop,
    p_unsafe_funeral_comm_hcw     = p_unsafe_funeral_comm_hcw,
    p_unsafe_funeral_hosp_hcw     = p_unsafe_funeral_hosp_hcw,
    p_unsafe_funeral_comm_genPop  = p_unsafe_funeral_comm_genPop,
    p_unsafe_funeral_hosp_genPop  = p_unsafe_funeral_hosp_genPop,
    ppe_coverage_hcw              = ppe_coverage_hcw,
    prop_etu                      = prop_etu,
    obv_pep_coverage          = obv_pep_coverage,
    obv_pep_adherence             = obv_pep_adherence,
    trace_coverage                = trace_coverage,
    prob_hospitalised_traced      = prob_hospitalised_traced,
    baseline_risk_genPop          = baseline_risk_genPop,
    baseline_risk_hcw             = baseline_risk_hcw,
    baseline_risk_funeral         = baseline_risk_funeral
  )
  ## Build the sampling grid from ALL time-varying inputs -- the probabilities
  ## above plus the positive-valued curves below -- so the upfront check lands on
  ## every changepoint, including those of hospitalisation_delay_factor,
  ## obv_pep_dpc, and the time-varying mean-offspring parameters.
  grid_inputs <- c(
    sanity_params,
    list(
      hospitalisation_delay_factor = hospitalisation_delay_factor,
      obv_pep_dpc                  = obv_pep_dpc,
      mn_contacts_genPop           = mn_contacts_genPop,
      mn_contacts_hcw              = mn_contacts_hcw,
      mn_contacts_funeral          = mn_contacts_funeral,
      prob_hospitalised_multiplier_traced = prob_hospitalised_multiplier_traced,
      prob_hospitalised_traced            = prob_hospitalised_traced,
      onset_to_hospitalisation_traced     = onset_to_hospitalisation_traced,
      presymptomatic_transmission         = presymptomatic_transmission
    )
  )
  sanity_grid <- build_sanity_grid(grid_inputs)

  for (nm in names(sanity_params)) {
    check_probability_on_grid(sanity_params[[nm]], sanity_grid, nm)
  }

  ## hospitalisation_delay_factor is strictly positive (a multiplier), not a probability.
  check_positive_on_grid(hospitalisation_delay_factor, sanity_grid,
                         "hospitalisation_delay_factor")

  ## mn_contacts_* are strictly positive NB means and may be scalars or functions
  ## of absolute calendar time. They are resolved inside the offspring functions
  ## (genPop/HCW at the parent's infection time, funeral at the parent's death
  ## time); here we only sanity-check positivity across the simulation horizon.
  check_positive_on_grid(mn_contacts_genPop,  sanity_grid, "mn_contacts_genPop")
  check_positive_on_grid(mn_contacts_hcw,     sanity_grid, "mn_contacts_hcw")
  check_positive_on_grid(mn_contacts_funeral, sanity_grid, "mn_contacts_funeral")

  ## The traced-case hospitalisation multiplier is strictly positive (it scales a probability).
  check_positive_on_grid(prob_hospitalised_multiplier_traced, sanity_grid,
                         "prob_hospitalised_multiplier_traced")
  ## The traced admission delay may legitimately be zero (same-day admission).
  check_nonneg_on_grid(onset_to_hospitalisation_traced, sanity_grid,
                       "onset_to_hospitalisation_traced")

  ## The highest-risk tier's per-contact transmission probability is
  ## baseline_risk(t) * max_relative_risk and must stay a valid probability across the
  ## whole horizon. Catch an infeasible combination here rather than mid-run.
  check_top_tier_probability <- function(baseline_risk, risk, param_name) {
    values <- resolve_time_varying(baseline_risk, sanity_grid, param_name)
    top <- values * risk$max_relative_risk
    if (any(top > 1 + 1e-12)) {
      bad <- which(top > 1 + 1e-12)
      i <- bad[1]
      stop(sprintf(
        "`%s` = %s at t = %s gives the highest-risk tier a transmission probability of %s (must be <= 1). Lower the baseline risk, raise the mean contact number, or narrow the relative-risk spread.",
        param_name, format(round(values[i], 6)), format(round(sanity_grid[i], 3)),
        format(round(top[i], 4))
      ), call. = FALSE)
    }
    invisible(NULL)
  }
  check_top_tier_probability(baseline_risk_genPop,  risk_genPop,  "baseline_risk_genPop")
  check_top_tier_probability(baseline_risk_hcw,     risk_hcw,     "baseline_risk_hcw")
  check_top_tier_probability(baseline_risk_funeral, risk_funeral, "baseline_risk_funeral")

  if (!is.logical(presymptomatic_transmission) || length(presymptomatic_transmission) != 1L ||
      is.na(presymptomatic_transmission)) {
    stop("`presymptomatic_transmission` must be a single logical value.", call. = FALSE)
  }

  ## The absolute traced hospitalisation probability and the multiplier are two ways of
  ## saying the same thing; accepting both would silently ignore one of them.
  if (!is.null(prob_hospitalised_traced) &&
      !(is.numeric(prob_hospitalised_multiplier_traced) &&
        length(prob_hospitalised_multiplier_traced) == 1L &&
        isTRUE(all.equal(prob_hospitalised_multiplier_traced, 1)))) {
    stop("Supply either `prob_hospitalised_traced` (an absolute probability) or `prob_hospitalised_multiplier_traced` (a multiplier), not both.",
         call. = FALSE)
  }

  ####################################################################################
  ### Step 1c: Pre-flight checks that need their own random draws
  ###
  ### Both of these sample from the user's delay distributions, which would consume
  ### RNG and shift every subsequent draw in the simulation. The whole block therefore
  ### saves the random seed on entry and restores it on exit, so the checks are
  ### invisible to the simulated trajectory.
  ###
  ###  (a) `onset_to_hospitalisation_traced` is only meaningful if it is actually faster
  ###      than the untraced pathway. A value at or above the untraced delay distribution
  ###      would silently do nothing (the delay is applied as a cap), so warn rather than
  ###      let a scenario quietly have no tracing effect.
  ###  (b) Estimate how much transmission happens before symptom onset, and warn if it is
  ###      substantial. This matters because fast admission of traced cases can only act on
  ###      post-onset transmission, so a large presymptomatic share caps what tracing can
  ###      ever achieve.
  ####################################################################################
  presymptomatic_share <- NULL
  run_preflight <- (!is.null(onset_to_hospitalisation_traced)) ||
    (isTRUE(check_presymptomatic) && isTRUE(presymptomatic_transmission))

  if (run_preflight) {
    seed_existed <- exists(".Random.seed", envir = globalenv())
    saved_seed <- if (seed_existed) get(".Random.seed", envir = globalenv()) else NULL

    if (!is.null(onset_to_hospitalisation_traced)) {
      untraced_delays <- onset_to_hospitalisation(n = 2000) *
        resolve_positive_time_varying(hospitalisation_delay_factor, 0,
                                      "hospitalisation_delay_factor")
      traced_delay_0 <- resolve_time_varying(onset_to_hospitalisation_traced,
                                             sanity_grid, "onset_to_hospitalisation_traced")
      q25 <- stats::quantile(untraced_delays, 0.25, names = FALSE)
      if (max(traced_delay_0) >= q25) {
        warning(sprintf(
          paste0("`onset_to_hospitalisation_traced` (max %.2f days) is not clearly below the untraced ",
                 "onset-to-admission delay (25th percentile %.2f days, median %.2f). Because the traced ",
                 "delay caps rather than replaces each case's own delay, tracing will have little or no ",
                 "effect on admission timing. Set a smaller value."),
          max(traced_delay_0), q25, stats::median(untraced_delays)
        ), call. = FALSE)
      }
    }

    if (isTRUE(check_presymptomatic) && isTRUE(presymptomatic_transmission)) {
      ps <- approx_presymptomatic_transmission(
        list(
          incubation_period        = incubation_period,
          prob_symptomatic         = prob_symptomatic,
          prob_death_comm          = prob_death_comm,
          prob_death_hosp          = prob_death_hosp,
          prob_hospitalised_genPop = prob_hospitalised_genPop,
          prob_hospitalised_hcw    = prob_hospitalised_hcw,
          onset_to_death           = onset_to_death,
          onset_to_recovery        = onset_to_recovery,
          onset_to_hospitalisation = onset_to_hospitalisation,
          hospitalisation_delay_factor = hospitalisation_delay_factor,
          hospitalisation_to_death     = hospitalisation_to_death,
          hospitalisation_to_recovery  = hospitalisation_to_recovery,
          Tg_shape_genPop = Tg_shape_genPop, Tg_rate_genPop = Tg_rate_genPop,
          Tg_shape_hcw    = Tg_shape_hcw,    Tg_rate_hcw    = Tg_rate_hcw
        ),
        n = 10000
      )
      presymptomatic_share <- ps
      if (ps$genPop > presymptomatic_warn_threshold) {
        warning(sprintf(
          paste0("About %.0f%% of genPop transmission in this parameter set happens before the ",
                 "infector develops symptoms (HCW: %.0f%%). Contact tracing and admission-based ",
                 "interventions can only act on the remainder. Set `presymptomatic_transmission = FALSE` ",
                 "to remove it, or `check_presymptomatic = FALSE` to silence this."),
          100 * ps$genPop, 100 * ps$hcw
        ), call. = FALSE)
      }
    }

    if (seed_existed) {
      assign(".Random.seed", saved_seed, envir = globalenv())
    }
  }

  ## obv_pep_dpc is non-negative (0 = same-day treatment is a meaningful boundary value).
  check_nonneg_on_grid(obv_pep_dpc, sanity_grid, "obv_pep_dpc")

  ## obv_pep_dpc_sd (optional): when set, DPC is drawn per-recipient from a Gamma with
  ## mean obv_pep_dpc(t) and this fixed standard deviation; NULL keeps DPC deterministic.
  ## Fail fast here rather than per-parent inside the gate.
  if (!is.null(obv_pep_dpc_sd) &&
      (!is.numeric(obv_pep_dpc_sd) || length(obv_pep_dpc_sd) != 1L ||
       !is.finite(obv_pep_dpc_sd) || obv_pep_dpc_sd <= 0)) {
    stop("`obv_pep_dpc_sd` must be NULL or a single finite positive numeric.", call. = FALSE)
  }

  ## obv_pep_efficacy_args (optional): named overrides for the built-in efficacy curve. The curve
  ## is used when obv_pep_efficacy is NULL or a scalar (the scalar is the curve's E0). Validate
  ## structure/names/conflict, then force one curve evaluation so bad values (e.g. E0 outside
  ## [0, 1]) fail fast before the simulation loop. Skip a custom function (calling it could consume
  ## RNG) and the trivial NULL-efficacy/no-overrides case (curve defaults are known-valid).
  validate_obv_efficacy_args(obv_pep_efficacy, obv_pep_efficacy_args)
  if (!is.function(obv_pep_efficacy) &&
      (!is.null(obv_pep_efficacy) || length(obv_pep_efficacy_args))) {
    invisible(resolve_obv_efficacy(obv_pep_efficacy, 0, obv_pep_efficacy_args = obv_pep_efficacy_args))
  }

  ## prob_death_hosp must not exceed prob_death_comm (so second_chance_death_prob <= 1).
  ## Currently both are scalars; if they become time-varying in future, this still
  ## works because resolve_time_varying recycles scalars across the grid.
  if (!is.null(prob_death_comm) && !is.null(prob_death_hosp)) {
    pdc <- resolve_time_varying(prob_death_comm, sanity_grid, "prob_death_comm")
    pdh <- resolve_time_varying(prob_death_hosp, sanity_grid, "prob_death_hosp")
    if (any(pdh > pdc)) {
      bad <- which(pdh > pdc)
      show <- bad[seq_len(min(3L, length(bad)))]
      stop(sprintf(
        "`prob_death_hosp` must be <= `prob_death_comm` at all times (so the second-chance ratio is a valid probability); violation at t = %s (prob_death_hosp = %s, prob_death_comm = %s).",
        paste(round(sanity_grid[show], 3), collapse = ", "),
        paste(round(pdh[show], 4), collapse = ", "),
        paste(round(pdc[show], 4), collapse = ", ")
      ), call. = FALSE)
    }
  }


  ## Initialise the susceptible population
  susc <- population - initial_immune

  ## Initialise the HCW population
  hcw_total <- round(hcw_per_capita * population)
  if (hcw_total <= 0) {
    stop("number of hcws is <= 0 as currently specified by hcw_per_capita and population")
  }
  hcw_available <- hcw_total

  ## Preallocate data frame -
  max_cases <- check_final_size
  tdf <- data.frame(
    id                             = integer(max_cases),   # id of the infected individual
    class                          = NA_character_,
    infection_location             = NA_character_,
    parent                         = integer(max_cases),   # ancestor of the infected individual i.e. the parent
    generation                     = integer(max_cases),   # generation of the infected individual i.e. how many infections precede them in the transmission chain
    time_infection_relative        = NA_real_,             # time of the infection relative to the parent
    time_infection_absolute        = NA_real_,             # time of the infection in absolute calendar time (i.e. since start of outbreak)
    incubation_period              = NA_real_,
    symptomatic                    = NA,                   # are they symptomatic?
    time_symptom_onset_relative    = NA_real_,             # time of symptom onset relative to the parent
    time_symptom_onset_absolute    = NA_real_,             # time of symptom onset in absolute calendar time (i.e. since start of outbreak)
    hospitalisation                = FALSE,
    time_hospitalisation_relative  = NA_real_,
    time_hospitalisation_absolute  = NA_real_,
    outcome                        = FALSE,         # what the outcome is for that individual
    outcome_location               = NA_character_,
    time_outcome_relative          = NA_real_,
    time_outcome_absolute          = NA_real_,
    funeral_safety                 = NA_character_,
    contact_risk_level             = NA_integer_,          # risk tier of the contact that produced this case
    contact_risk_category          = NA_character_,        # its label
    traced                         = rep(FALSE, max_cases),# was that contact reached by contact tracing?
    obv_pep_eligible               = rep(FALSE, max_cases),
    obv_pep_received               = rep(FALSE, max_cases),
    obv_pep_adherent               = rep(FALSE, max_cases),
    obv_pep_dpc                    = rep(NA_real_, max_cases),
    n_offspring                    = integer(max_cases),
    offspring_generated            = FALSE,
    stringsAsFactors = FALSE
  )

  #########################################################################
  ### Step 2: Initialise conditions and features of the seeding cases
  #########################################################################

  ## Deciding whether the seeding cases are symptomatic, and if so, when they develop symptoms
  seeding_cases_time_infection <- seq(from = 0, to = 0.01, length.out = seeding_cases)
  seeding_cases_incubation <- incubation_period(n = seeding_cases)
  seeding_cases_symptomatic <- as.logical(rbinom(n = seeding_cases, size = 1, prob = prob_symptomatic))
  seeding_cases_symptom_onset <- rep(NA_real_, seeding_cases)
  seeding_cases_symptom_onset[seeding_cases_symptomatic] <- seeding_cases_incubation[seeding_cases_symptomatic]
  seeding_cases_symptom_onset_absolute <- seeding_cases_time_infection + seeding_cases_symptom_onset

  ## Deciding on the outcome for the seeding cases, and if so, when that outcome occurs
  seeding_cases_outcome <- rep(FALSE, seeding_cases)
  seeding_cases_outcome[seeding_cases_symptomatic] <- as.logical(rbinom(n = sum(seeding_cases_symptomatic), size = 1, prob = prob_death_comm))
  seeding_cases_outcome_time <- rep(NA_real_, seeding_cases)
  seeding_cases_outcome_time[seeding_cases_outcome] <- seeding_cases_incubation[seeding_cases_outcome] + onset_to_death(n = sum(seeding_cases_outcome))
  seeding_cases_outcome_time[!seeding_cases_outcome] <- seeding_cases_incubation[!seeding_cases_outcome] + onset_to_recovery(n = sum(!seeding_cases_outcome))
  seeding_cases_outcome_time_absolute <- seeding_cases_time_infection + seeding_cases_outcome_time

  ## Deciding funeral safety for seed cases who die. Seed cases are genPop and
  ## their deaths occur in the community, so use p_unsafe_funeral_comm_genPop at
  ## each seed case's death time rather than forcing all seed funerals unsafe.
  seeding_cases_funeral_safety <- rep(NA_character_, seeding_cases)
  if (any(seeding_cases_outcome)) {
    seeding_cases_p_unsafe_funeral <- resolve_probability(p_unsafe_funeral_comm_genPop,
                                                          seeding_cases_outcome_time_absolute[seeding_cases_outcome],
                                                          "p_unsafe_funeral_comm_genPop")
    seeding_cases_funeral_safety[seeding_cases_outcome] <- ifelse(
      rbinom(n = sum(seeding_cases_outcome), size = 1, prob = seeding_cases_p_unsafe_funeral) == 1,
      "unsafe", "safe"
    )
  }

  ## Initialising the dataframe with the seed cases and their attributes
  tdf[1:seeding_cases, ] <- data.frame(
    id                             = seq_len(seeding_cases),
    class                          = rep("genPop", seeding_cases),
    infection_location             = rep("community", seeding_cases),
    parent                         = NA_character_,
    generation                     = 1,
    time_infection_relative        = seeding_cases_time_infection,
    time_infection_absolute        = seeding_cases_time_infection,
    incubation_period              = seeding_cases_incubation,
    symptomatic                    = seeding_cases_symptomatic,
    time_symptom_onset_relative    = seeding_cases_symptom_onset,
    time_symptom_onset_absolute    = seeding_cases_symptom_onset_absolute,
    hospitalisation                = rep(FALSE, seeding_cases),
    time_hospitalisation_relative  = NA_real_,
    time_hospitalisation_absolute  = NA_real_,
    outcome                        = seeding_cases_outcome,
    outcome_location               = rep("community", seeding_cases),
    time_outcome_relative          = seeding_cases_outcome_time,
    time_outcome_absolute          = seeding_cases_outcome_time_absolute,
    funeral_safety                 = seeding_cases_funeral_safety,
    ## Seed cases were not produced by a contact, so they carry no risk tier and are
    ## never traced (there is no index case to trace them from).
    contact_risk_level             = NA_integer_,
    contact_risk_category          = NA_character_,
    traced                         = rep(FALSE, seeding_cases),
    obv_pep_eligible               = rep(FALSE, seeding_cases),
    obv_pep_received               = rep(FALSE, seeding_cases),
    obv_pep_adherent               = rep(FALSE, seeding_cases),
    obv_pep_dpc                    = rep(NA_real_, seeding_cases),
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors = FALSE
  )

  ## --- Columnar working store -------------------------------------------------
  ## The hot loop appends offspring by writing into pre-allocated *standalone*
  ## column vectors (refcount 1 => in-place `[<-`) rather than the row-block
  ## assignment `tdf[rows, ] <- df`, which copied the whole frame on every
  ## iteration (the dominant cost + GC churn at scale; see dev/profile.html).
  ## Columns are pulled out *after* seeding so they inherit the exact
  ## post-seeding column types (e.g. `parent` is character because the seed
  ## block assigns NA_character_ into the integer column; `generation` is double
  ## because the seed assigns `1`). Reassembled into a data.frame once after the
  ## loop, the result is identical -- values, types and order -- to the old path.
  v_id                            <- tdf$id
  v_class                         <- tdf$class
  v_infection_location            <- tdf$infection_location
  v_parent                        <- tdf$parent
  v_generation                    <- tdf$generation
  v_time_infection_relative       <- tdf$time_infection_relative
  v_time_infection_absolute       <- tdf$time_infection_absolute
  v_incubation_period             <- tdf$incubation_period
  v_symptomatic                   <- tdf$symptomatic
  v_time_symptom_onset_relative   <- tdf$time_symptom_onset_relative
  v_time_symptom_onset_absolute   <- tdf$time_symptom_onset_absolute
  v_hospitalisation               <- tdf$hospitalisation
  v_time_hospitalisation_relative <- tdf$time_hospitalisation_relative
  v_time_hospitalisation_absolute <- tdf$time_hospitalisation_absolute
  v_outcome                       <- tdf$outcome
  v_outcome_location              <- tdf$outcome_location
  v_time_outcome_relative         <- tdf$time_outcome_relative
  v_time_outcome_absolute         <- tdf$time_outcome_absolute
  v_funeral_safety                <- tdf$funeral_safety
  v_contact_risk_level            <- tdf$contact_risk_level
  v_contact_risk_category         <- tdf$contact_risk_category
  v_traced                        <- tdf$traced
  v_obv_pep_eligible              <- tdf$obv_pep_eligible
  v_obv_pep_received              <- tdf$obv_pep_received
  v_obv_pep_adherent              <- tdf$obv_pep_adherent
  v_obv_pep_dpc                   <- tdf$obv_pep_dpc
  v_n_offspring                   <- tdf$n_offspring
  v_offspring_generated           <- tdf$offspring_generated
  rm(tdf)
  ## Number of filled rows. Rows fill densely from the top (seed block, then one
  ## contiguous block per parent) and id == row index throughout, so this single
  ## counter equals both `max(which(!is.na(time_infection_absolute)))` (next free
  ## row) and `max(id)` -- replacing those per-iteration O(N) scans.
  n_filled <- seeding_cases

  #################################################################################
  ### Step 3: Loop through infections and generate offspring for each of them
  #################################################################################
  ## While we haven't hit the simulation cap size (check_final_size) and any infections exist where we have not yet generated the requisite offspring,
  ## continue to generate infections.
  ## `n_filled <= check_final_size` is exactly the old `nrow(tdf) <= check_final_size`:
  ## tdf was pre-allocated to max_cases == check_final_size, so nrow only exceeded
  ## the cap once an append extended it past max_cases, i.e. once n_filled did.
  while (any(is.na(v_n_offspring)) && susc > 0 && n_filled <= check_final_size) {

    #############################################################################################
    ## Step 1: Get earliest infection not yet expanded to act as a parent, and their attributes
    #############################################################################################
    parent_time_infection <- min(v_time_infection_absolute[!v_offspring_generated & !is.na(v_time_infection_absolute)])
    idx <- which(v_time_infection_absolute == parent_time_infection & !v_offspring_generated)[1]
    ## Rebuild the single parent row as a 1-row data.frame (the offspring/
    ## completion functions read it via `parent_info$<field>`). Same columns,
    ## same types as the old `tdf[idx, ]`; the row label differs (1 vs idx) but
    ## nothing downstream reads rownames(parent_info).
    parent_info <- data.frame(
      id                            = v_id[idx],
      class                         = v_class[idx],
      infection_location            = v_infection_location[idx],
      parent                        = v_parent[idx],
      generation                    = v_generation[idx],
      time_infection_relative       = v_time_infection_relative[idx],
      time_infection_absolute       = v_time_infection_absolute[idx],
      incubation_period             = v_incubation_period[idx],
      symptomatic                   = v_symptomatic[idx],
      time_symptom_onset_relative   = v_time_symptom_onset_relative[idx],
      time_symptom_onset_absolute   = v_time_symptom_onset_absolute[idx],
      hospitalisation               = v_hospitalisation[idx],
      time_hospitalisation_relative = v_time_hospitalisation_relative[idx],
      time_hospitalisation_absolute = v_time_hospitalisation_absolute[idx],
      outcome                       = v_outcome[idx],
      outcome_location              = v_outcome_location[idx],
      time_outcome_relative         = v_time_outcome_relative[idx],
      time_outcome_absolute         = v_time_outcome_absolute[idx],
      funeral_safety                = v_funeral_safety[idx],
      contact_risk_level            = v_contact_risk_level[idx],
      contact_risk_category         = v_contact_risk_category[idx],
      traced                        = v_traced[idx],
      obv_pep_eligible              = v_obv_pep_eligible[idx],
      obv_pep_received              = v_obv_pep_received[idx],
      obv_pep_adherent              = v_obv_pep_adherent[idx],
      obv_pep_dpc                   = v_obv_pep_dpc[idx],
      n_offspring                   = v_n_offspring[idx],
      offspring_generated           = v_offspring_generated[idx],
      stringsAsFactors              = FALSE
    )
    if (!(parent_info$class %in% c("genPop", "HCW"))) {
      stop("error with parent class")
    }

    ###################################################################################################################
    ### Step 2: Generate offspring associated with community and (if hospitalised) healthcare associated transmission
    ###################################################################################################################
    ## Pass scalar-or-time-varying response parameters into the offspring functions
    ## directly. Those functions know the candidate transmission times, so they can
    ## resolve PPE coverage and post-admission hospital quarantine/ETU efficacy at the
    ## actual absolute calendar time of each candidate hospital exposure. PPE thinning
    ## is ppe_coverage_hcw(t) * ppe_efficacy; hospital quarantine efficacy is
    ## calculated inside the offspring functions as a prop_etu(t)-weighted mixture
    ## of etu_efficacy and general_hospital_quarantine_efficacy.

    if (parent_info$class == "genPop") {
      offspring_community_healthcare_df <- offspring_function_genPop(parent_info = parent_info,
                                                                     mn_contacts_genPop = mn_contacts_genPop,
                                                                     overdisp_contacts_genPop = overdisp_contacts_genPop,
                                                                     baseline_risk_genPop = baseline_risk_genPop,
                                                                     contact_risk_genPop = risk_genPop,
                                                                     Tg_shape_genPop = Tg_shape_genPop,
                                                                     Tg_rate_genPop = Tg_rate_genPop,
                                                                     trace_coverage = trace_coverage,
                                                                     return_contact_log = return_contact_log,
                                                                     presymptomatic_transmission = presymptomatic_transmission,
                                                                     prop_etu = prop_etu,
                                                                     etu_efficacy = etu_efficacy,
                                                                     general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
                                                                     obv_pep_enabled = obv_pep_enabled,
                                                                     obv_pep_coverage = obv_pep_coverage,
                                                                     obv_pep_adherence = obv_pep_adherence,
                                                                     obv_pep_dpc = obv_pep_dpc,
                                                                     obv_pep_dpc_sd = obv_pep_dpc_sd,
                                                                     obv_pep_efficacy = obv_pep_efficacy,
                                                                     obv_pep_efficacy_args = obv_pep_efficacy_args,
                                                                     obv_pep_target_class = obv_pep_target_class,
                                                                     obv_pep_target_locations = obv_pep_target_locations,
                                                                     ppe_coverage_hcw = ppe_coverage_hcw,
                                                                     ppe_efficacy = ppe_efficacy,
                                                                     prob_hcw_cond_genPop_comm = prob_hcw_cond_genPop_comm,
                                                                     prob_hcw_cond_genPop_hospital = prob_hcw_cond_genPop_hospital)
    } else if (parent_info$class == "HCW") {
      ## PPE/IPC efficacy is passed through unresolved. The HCW offspring
      ## function resolves it separately for each candidate pre-admission
      ## hospital transmission event.

      offspring_community_healthcare_df <- offspring_function_hcw(parent_info = parent_info,
                                                                  mn_contacts_hcw = mn_contacts_hcw,
                                                                  overdisp_contacts_hcw = overdisp_contacts_hcw,
                                                                  baseline_risk_hcw = baseline_risk_hcw,
                                                                  contact_risk_hcw = risk_hcw,
                                                                  Tg_shape_hcw = Tg_shape_hcw,
                                                                  Tg_rate_hcw = Tg_rate_hcw,
                                                                  prob_hospital_cond_hcw_preAdm = prob_hospital_cond_hcw_preAdm,
                                                                  trace_coverage = trace_coverage,
                                                                  return_contact_log = return_contact_log,
                                                                  presymptomatic_transmission = presymptomatic_transmission,
                                                                  ppe_coverage_hcw = ppe_coverage_hcw,
                                                                  ppe_efficacy = ppe_efficacy,
                                                                  prop_etu = prop_etu,
                                                                  etu_efficacy = etu_efficacy,
                                                                  general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
                                                                  obv_pep_enabled = obv_pep_enabled,
                                                                  obv_pep_coverage = obv_pep_coverage,
                                                                  obv_pep_adherence = obv_pep_adherence,
                                                                  obv_pep_dpc = obv_pep_dpc,
                                                                  obv_pep_dpc_sd = obv_pep_dpc_sd,
                                                                  obv_pep_efficacy = obv_pep_efficacy,
                                                                  obv_pep_efficacy_args = obv_pep_efficacy_args,
                                                                  obv_pep_target_class = obv_pep_target_class,
                                                                  obv_pep_target_locations = obv_pep_target_locations,
                                                                  prob_hcw_cond_hcw_comm = prob_hcw_cond_hcw_comm,
                                                                  prob_hcw_cond_hcw_hospital = prob_hcw_cond_hcw_hospital)
    }

    ## Accumulate OBV PEP counters from this offspring call (community/healthcare).
    step <- attr(offspring_community_healthcare_df, "obv_pep_num_treated", exact = TRUE)
    obv_num_treated$pre_eligible  <- obv_num_treated$pre_eligible  + step$pre_eligible
    obv_num_treated$pre_treated   <- obv_num_treated$pre_treated   + step$pre_treated
    obv_num_treated$pre_adherent  <- obv_num_treated$pre_adherent  + step$pre_adherent
    obv_num_treated$post_eligible <- obv_num_treated$post_eligible + step$post_eligible
    obv_num_treated$post_treated  <- obv_num_treated$post_treated  + step$post_treated
    obv_num_treated$post_adherent <- obv_num_treated$post_adherent + step$post_adherent
    obv_num_treated$prevented     <- obv_num_treated$prevented     + step$prevented
    pinfo <- attr(offspring_community_healthcare_df, "obv_pep_prevented_info", exact = TRUE)
    if (!is.null(pinfo) && nrow(pinfo) > 0) {
      obv_prevented_info_list[[length(obv_prevented_info_list) + 1L]] <- pinfo
    }

    ## Count HCWs generated and subtract from available pool
    n_hcw_community_healthcare <- sum(offspring_community_healthcare_df$class == "HCW")
    hcw_available <- hcw_available - n_hcw_community_healthcare

    #############################################################################################
    ### Step 3: Generate offspring associated with funeral transmission
    #############################################################################################
    offspring_funeral_df <- offspring_function_funeral(parent_info = parent_info,
                                                       mn_contacts_funeral = mn_contacts_funeral,
                                                       overdisp_contacts_funeral = overdisp_contacts_funeral,
                                                       baseline_risk_funeral = baseline_risk_funeral,
                                                       contact_risk_funeral = risk_funeral,
                                                       Tg_shape_funeral = Tg_shape_funeral,
                                                       Tg_rate_funeral = Tg_rate_funeral,
                                                       trace_coverage = trace_coverage,
                                                       return_contact_log = return_contact_log,
                                                       safe_funeral_efficacy = safe_funeral_efficacy,
                                                       obv_pep_enabled = obv_pep_enabled,
                                                       obv_pep_coverage = obv_pep_coverage,
                                                       obv_pep_adherence = obv_pep_adherence,
                                                       obv_pep_dpc = obv_pep_dpc,
                                                       obv_pep_dpc_sd = obv_pep_dpc_sd,
                                                       obv_pep_efficacy = obv_pep_efficacy,
                                                       obv_pep_efficacy_args = obv_pep_efficacy_args,
                                                       obv_pep_target_class = obv_pep_target_class,
                                                       obv_pep_target_locations = obv_pep_target_locations,
                                                       prob_hcw_cond_funeral_hcw = prob_hcw_cond_funeral_hcw,
                                                       prob_hcw_cond_funeral_genPop = prob_hcw_cond_funeral_genPop)

    ## Accumulate OBV PEP counters from this offspring call (funeral).
    step <- attr(offspring_funeral_df, "obv_pep_num_treated", exact = TRUE)
    obv_num_treated$pre_eligible  <- obv_num_treated$pre_eligible  + step$pre_eligible
    obv_num_treated$pre_treated   <- obv_num_treated$pre_treated   + step$pre_treated
    obv_num_treated$pre_adherent  <- obv_num_treated$pre_adherent  + step$pre_adherent
    obv_num_treated$post_eligible <- obv_num_treated$post_eligible + step$post_eligible
    obv_num_treated$post_treated  <- obv_num_treated$post_treated  + step$post_treated
    obv_num_treated$post_adherent <- obv_num_treated$post_adherent + step$post_adherent
    obv_num_treated$prevented     <- obv_num_treated$prevented     + step$prevented
    pinfo <- attr(offspring_funeral_df, "obv_pep_prevented_info", exact = TRUE)
    if (!is.null(pinfo) && nrow(pinfo) > 0) {
      obv_prevented_info_list[[length(obv_prevented_info_list) + 1L]] <- pinfo
    }

    ## Update hcw_available after funeral transmission
    n_hcw_funeral <- sum(offspring_funeral_df$class == "HCW")
    hcw_available <- hcw_available - n_hcw_funeral

    ## Collect this parent's full contact log across both routes. The rows flagged
    ## "infection" appear in the same order as the offspring rows, which is what lets
    ## case ids be filled in once those rows are appended below.
    combined_contact_log <- rbind(
      attr(offspring_community_healthcare_df, "contact_log", exact = TRUE),
      attr(offspring_funeral_df, "contact_log", exact = TRUE)
    )

    #################################################################################################################
    ### Step 4: Complete offspring information based on parent attributes and timings; and update parent information
    ##          (e.g. num_offspring, offspring_generated == TRUE etc)
    #################################################################################################################
    ## Completing offspring information if there are any. Combine the two
    ## offspring sources once (previously rbind'd twice -- guard + call).
    combined_offspring_df <- rbind(offspring_community_healthcare_df, offspring_funeral_df)
    if (nrow(combined_offspring_df) > 0) {
      complete_offspring_df <- complete_offspring_info(parent_info = parent_info,
                                                       offspring_dataframe = combined_offspring_df,
                                                       prob_symptomatic = prob_symptomatic,
                                                       prob_hospitalised_hcw = prob_hospitalised_hcw,
                                                       prob_hospitalised_genPop = prob_hospitalised_genPop,
                                                       prob_death_comm = prob_death_comm,
                                                       prob_death_hosp = prob_death_hosp,
                                                       p_unsafe_funeral_comm_hcw = p_unsafe_funeral_comm_hcw,
                                                       p_unsafe_funeral_hosp_hcw = p_unsafe_funeral_hosp_hcw,
                                                       p_unsafe_funeral_comm_genPop = p_unsafe_funeral_comm_genPop,
                                                       p_unsafe_funeral_hosp_genPop = p_unsafe_funeral_hosp_genPop,
                                                       onset_to_hospitalisation_traced = onset_to_hospitalisation_traced,
                                                       prob_hospitalised_traced = prob_hospitalised_traced,
                                                       prob_hospitalised_multiplier_traced = prob_hospitalised_multiplier_traced,
                                                       incubation_period = incubation_period,
                                                       onset_to_hospitalisation = onset_to_hospitalisation,
                                                       hospitalisation_delay_factor = hospitalisation_delay_factor,
                                                       hospitalisation_to_death = hospitalisation_to_death,
                                                       hospitalisation_to_recovery = hospitalisation_to_recovery,
                                                       onset_to_death = onset_to_death,
                                                       onset_to_recovery = onset_to_recovery)
      n_new <- nrow(complete_offspring_df)
      v_n_offspring[idx] <- n_new
    } else {
      n_new <- 0L
      v_n_offspring[idx] <- 0   # double literal preserved (old `tdf$n_offspring[idx] <- 0`
                                # promotes the integer column to double on first 0-offspring parent)
    }
    v_offspring_generated[idx] <- TRUE

    #################################################################################################################
    ### Step 5: Adding the complete offspring dataframe (complete_offspring_df) to the main column vectors
    #################################################################################################################
    ## If offspring exist, append them by writing into the pre-allocated standalone
    ## column vectors in place. The next free block is (n_filled + 1):(n_filled + n_new)
    ## and id == row index, so this replaces both the row-block data.frame copy and the
    ## per-iteration max(which(!is.na(...))) / max(id) scans.
    if (n_new > 0) {
      rows <- (n_filled + 1L):(n_filled + n_new)
      ## Link the contact log's realised infections to the case ids they became. The
      ## offspring functions emit their infection rows in candidate order and rbind
      ## preserves it, so the two line up one-for-one; check rather than assume.
      if (nrow(combined_contact_log) > 0) {
        inf_rows <- which(combined_contact_log$record_type == "infection")
        if (length(inf_rows) != n_new) {
          stop("Internal error: contact log infections do not align with the offspring rows.",
               call. = FALSE)
        }
        combined_contact_log$case_id[inf_rows] <- rows
      }
      v_id[rows]                            <- rows
      v_class[rows]                         <- complete_offspring_df$class
      v_infection_location[rows]            <- complete_offspring_df$infection_location
      v_parent[rows]                        <- complete_offspring_df$parent
      v_generation[rows]                    <- complete_offspring_df$generation
      v_time_infection_relative[rows]       <- complete_offspring_df$time_infection_relative
      v_time_infection_absolute[rows]       <- complete_offspring_df$time_infection_absolute
      v_incubation_period[rows]             <- complete_offspring_df$incubation_period
      v_symptomatic[rows]                   <- complete_offspring_df$symptomatic
      v_time_symptom_onset_relative[rows]   <- complete_offspring_df$time_symptom_onset_relative
      v_time_symptom_onset_absolute[rows]   <- complete_offspring_df$time_symptom_onset_absolute
      v_hospitalisation[rows]               <- complete_offspring_df$hospitalisation
      v_time_hospitalisation_relative[rows] <- complete_offspring_df$time_hospitalisation_relative
      v_time_hospitalisation_absolute[rows] <- complete_offspring_df$time_hospitalisation_absolute
      v_outcome[rows]                       <- complete_offspring_df$outcome
      v_outcome_location[rows]              <- complete_offspring_df$outcome_location
      v_time_outcome_relative[rows]         <- complete_offspring_df$time_outcome_relative
      v_time_outcome_absolute[rows]         <- complete_offspring_df$time_outcome_absolute
      v_funeral_safety[rows]                <- complete_offspring_df$funeral_safety
      v_contact_risk_level[rows]            <- complete_offspring_df$contact_risk_level
      v_contact_risk_category[rows]         <- complete_offspring_df$contact_risk_category
      v_traced[rows]                        <- complete_offspring_df$traced
      v_obv_pep_eligible[rows]              <- complete_offspring_df$obv_pep_eligible
      v_obv_pep_received[rows]              <- complete_offspring_df$obv_pep_received
      v_obv_pep_adherent[rows]              <- complete_offspring_df$obv_pep_adherent
      v_obv_pep_dpc[rows]                   <- complete_offspring_df$obv_pep_dpc
      v_n_offspring[rows]                   <- complete_offspring_df$n_offspring
      v_offspring_generated[rows]           <- complete_offspring_df$offspring_generated
      n_filled <- n_filled + n_new
    }
    if (return_contact_log && nrow(combined_contact_log) > 0) {
      contact_log_list[[length(contact_log_list) + 1L]] <- combined_contact_log
    }

    ## Deplete susceptibles
    susc <- susc - v_n_offspring[idx]
  }

  ############################################################################################
  ### Why did the loop stop?
  ###
  ### Three exits, and they mean very different things for an analysis:
  ###   "outbreak_ended"  every case was expanded -- the natural end, final size is real
  ###   "final_size_cap"  we ran out of budget with cases still waiting to be expanded, so the
  ###                     final size is CENSORED and measures the cap rather than transmission
  ###   "susceptibles"    the susceptible pool was exhausted
  ### A censored run silently looks like a controlled one if you only read the final size, so
  ### the reason is surfaced in sim_info and (unless quiet) announced at the end of the run.
  ############################################################################################
  n_unexpanded <- sum(is.na(v_n_offspring))
  stop_reason <- if (n_unexpanded == 0L) {
    "outbreak_ended"
  } else if (susc <= 0) {
    "susceptibles"
  } else {
    "final_size_cap"
  }
  hit_final_size_cap <- identical(stop_reason, "final_size_cap")
  if (!quiet) {
    if (identical(stop_reason, "final_size_cap")) {
      warning(sprintf(
        paste0("Simulation stopped at the `check_final_size` cap (%d) with %d case(s) still ",
               "unexpanded. The final size is CENSORED -- it measures the cap, not transmission. ",
               "Raise `check_final_size`, or filter on `sim_info$hit_final_size_cap` before ",
               "comparing final sizes across scenarios."),
        check_final_size, n_unexpanded), call. = FALSE)
    } else {
      message(sprintf("Simulation ended: %s (%d cases).",
                      if (identical(stop_reason, "susceptibles")) "susceptible pool exhausted"
                      else "outbreak died out", n_filled))
    }
  }

  ############################################################################################
  ### Final tidy of dataframe and then outputting it
  ############################################################################################
  ## Reassemble the frame once from the column vectors (column order matches the
  ## Step-1 pre-allocation exactly). Includes the unfilled pre-allocated tail
  ## (NA time_infection_absolute), which order() sends to the bottom -- identical
  ## to the old behaviour, which never truncated the pre-allocated frame.
  tdf <- data.frame(
    id                            = v_id,
    class                         = v_class,
    infection_location            = v_infection_location,
    parent                        = v_parent,
    generation                    = v_generation,
    time_infection_relative       = v_time_infection_relative,
    time_infection_absolute       = v_time_infection_absolute,
    incubation_period             = v_incubation_period,
    symptomatic                   = v_symptomatic,
    time_symptom_onset_relative   = v_time_symptom_onset_relative,
    time_symptom_onset_absolute   = v_time_symptom_onset_absolute,
    hospitalisation               = v_hospitalisation,
    time_hospitalisation_relative = v_time_hospitalisation_relative,
    time_hospitalisation_absolute = v_time_hospitalisation_absolute,
    outcome                       = v_outcome,
    outcome_location              = v_outcome_location,
    time_outcome_relative         = v_time_outcome_relative,
    time_outcome_absolute         = v_time_outcome_absolute,
    funeral_safety                = v_funeral_safety,
    contact_risk_level            = v_contact_risk_level,
    contact_risk_category         = v_contact_risk_category,
    traced                        = v_traced,
    obv_pep_eligible              = v_obv_pep_eligible,
    obv_pep_received              = v_obv_pep_received,
    obv_pep_adherent              = v_obv_pep_adherent,
    obv_pep_dpc                   = v_obv_pep_dpc,
    n_offspring                   = v_n_offspring,
    offspring_generated           = v_offspring_generated,
    stringsAsFactors              = FALSE
  )
  tdf <- tdf[order(tdf$time_infection_absolute, tdf$id), ]
  rownames(tdf) <- NULL

  #########################################################################################
  ### Deferred OBV PEP "prevented deaths" counterfactual
  ###
  ### For each infection the OBV gate prevented, decide whether it WOULD have died had
  ### it occurred, by replaying it through the SAME outcome model as realised cases
  ### (complete_offspring_info: symptomatic -> potential hospitalisation -> community CFR
  ### -> hospital second-chance). This is a "direct" count -- the would-be death of each
  ### prevented index infection only -- mirroring `prevented`, which likewise excludes the
  ### averted onward transmission chains.
  ###
  ### The completed frame (`prevented_completed`) is also returned in `out` so callers
  ### can read each averted infection's counterfactual natural history directly --
  ### notably `time_infection_absolute`, the would-be calendar infection time. It stays
  ### NULL when nothing was prevented (so a caller can `rbind` it across replicates and
  ### the empty runs simply drop out).
  ###
  ### Why after the loop: these draws consume RNG, so doing them inline would shift every
  ### subsequent draw and change the simulated trajectory. Run once the tree is finalised,
  ### nothing downstream in this call depends on them, and each branching_process_main()
  ### call re-seeds up front -- so the trajectory and every existing output are byte-for-byte
  ### identical to a run without this counter. Skipped entirely when nothing was prevented,
  ### so zero-prevention (incl. obv-disabled) runs draw nothing extra.
  #########################################################################################
  prevented_completed <- NULL
  if (length(obv_prevented_info_list) > 0) {
    obv_prevented_info <- do.call(rbind, obv_prevented_info_list)
    if (nrow(obv_prevented_info) > 0) {
      ## Carry each prevented infection's absolute infection time as a "relative" time
      ## against a zero-time dummy parent, so complete_offspring_info's clock
      ## (parent_abs + relative) resolves time-varying parameters at the true time.
      prevented_offspring <- data.frame(
        infection_location      = obv_prevented_info$infection_location,
        time_infection_relative = obv_prevented_info$time_infection_absolute,
        class                   = obv_prevented_info$class,
        stringsAsFactors        = FALSE
      )
      dummy_parent <- data.frame(
        id                      = NA_integer_,
        generation              = NA_integer_,
        time_infection_absolute = 0,
        stringsAsFactors        = FALSE
      )
      prevented_completed <- complete_offspring_info(
        parent_info                  = dummy_parent,
        offspring_dataframe          = prevented_offspring,
        prob_symptomatic             = prob_symptomatic,
        prob_hospitalised_hcw        = prob_hospitalised_hcw,
        prob_hospitalised_genPop     = prob_hospitalised_genPop,
        prob_death_comm              = prob_death_comm,
        prob_death_hosp              = prob_death_hosp,
        p_unsafe_funeral_comm_hcw    = p_unsafe_funeral_comm_hcw,
        p_unsafe_funeral_hosp_hcw    = p_unsafe_funeral_hosp_hcw,
        p_unsafe_funeral_comm_genPop = p_unsafe_funeral_comm_genPop,
        p_unsafe_funeral_hosp_genPop = p_unsafe_funeral_hosp_genPop,
        incubation_period            = incubation_period,
        onset_to_hospitalisation     = onset_to_hospitalisation,
        hospitalisation_delay_factor = hospitalisation_delay_factor,
        hospitalisation_to_death     = hospitalisation_to_death,
        hospitalisation_to_recovery  = hospitalisation_to_recovery,
        onset_to_death               = onset_to_death,
        onset_to_recovery            = onset_to_recovery
      )
      obv_num_treated$prevented_deaths <- sum(prevented_completed$outcome, na.rm = TRUE)
    }
  }

  attr(tdf, "hcw_total") <- hcw_total
  attr(tdf, "hcw_infected") <- hcw_total - hcw_available
  attr(tdf, "hcw_remaining") <- hcw_available
  attr(tdf, "obv_pep_num_treated") <- obv_num_treated

  ## Assemble the full contact log: one row per contact generated over the whole run,
  ## ordered by parent then by the order contacts were drawn. `record_type` says whether
  ## a contact became an infection, and `case_id` joins those rows to `tdf`.
  contact_log <- if (length(contact_log_list) > 0) {
    cl <- do.call(rbind, contact_log_list)
    rownames(cl) <- NULL
    cl
  } else {
    empty_contact_log()
  }

  out <- list(
    tdf = tdf,
    ## Full contact log (see the @return docs).
    contact_log = contact_log,
    ## Counterfactual completed offspring info for the infections OBV prevented
    ## (averted index infections only; NULL when nothing was prevented). See the
    ## deferred-counterfactual block above and the @return docs for column notes.
    prevented_completed = prevented_completed,
    sim_info = list(
      population          = population,
      hcw_per_capita      = hcw_per_capita,
      hcw_total           = hcw_total,
      seed                = seed,
      obv_pep_enabled     = obv_pep_enabled,
      obv_pep_num_treated = obv_num_treated,
      ## Contact-first calibration: the risk structures actually used, the baseline
      ## per-contact risks (solved from r0_target when that was supplied), and the
      ## R0 inversion diagnostics.
      contact_risk_genPop   = risk_genPop,
      contact_risk_hcw      = risk_hcw,
      contact_risk_funeral  = risk_funeral,
      baseline_risk_genPop  = baseline_risk_genPop,
      baseline_risk_hcw     = baseline_risk_hcw,
      baseline_risk_funeral = baseline_risk_funeral,
      r0_target             = r0_target,
      r0_solution           = r0_solution,
      ## Why the run stopped. `hit_final_size_cap = TRUE` means the final size is censored
      ## by `check_final_size` and must not be compared across scenarios at face value.
      stop_reason           = stop_reason,
      hit_final_size_cap    = hit_final_size_cap,
      n_unexpanded          = n_unexpanded
    )
  )

  return(out)
}

