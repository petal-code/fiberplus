## Tests for the OBV PEP gate.
## Covers:
##   - The two-phase semantics of apply_obv_pep_gate()
##   - Disabled / zero-coverage no-op behaviour
##   - Full-prevention regression (coverage = adherence = efficacy = 1 zeroes
##     out HCW infections at hospital)
##   - Boundary behaviour of obv_pep_efficacy_from_dpc()
##   - Set-nesting invariants between the 7 counters

## --- helpers -----------------------------------------------------------

make_parent_info_genPop_hospitalised <- function(time_infection_absolute = 0,
                                                  time_to_hospitalisation = 0.001,
                                                  time_to_outcome = 30) {
  data.frame(
    id                             = 1L,
    class                          = "genPop",
    infection_location             = "community",
    parent                         = NA_integer_,
    generation                     = 1L,
    time_infection_relative        = 0,
    time_infection_absolute        = time_infection_absolute,
    incubation_period              = 5,
    symptomatic                    = TRUE,
    time_symptom_onset_relative    = 5,
    time_symptom_onset_absolute    = time_infection_absolute + 5,
    hospitalisation                = TRUE,
    time_hospitalisation_relative  = time_to_hospitalisation,
    time_hospitalisation_absolute  = time_infection_absolute + time_to_hospitalisation,
    outcome                        = TRUE,
    outcome_location               = "hospital",
    time_outcome_relative          = time_to_outcome,
    time_outcome_absolute          = time_infection_absolute + time_to_outcome,
    funeral_safety                 = "unsafe",
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors               = FALSE
  )
}

obv_off_args <- function() {
  list(
    mn_contacts_genPop          = 500,
    baseline_risk_genPop          = 1,
    overdisp_contacts_genPop    = 200,
    Tg_shape_genPop              = 4,
    Tg_rate_genPop               = 1,
    prop_etu                     = 1,
    etu_efficacy                 = 0,
    general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw             = 0,
    ppe_efficacy                 = 1,
    prob_hcw_cond_genPop_comm    = 0,
    prob_hcw_cond_genPop_hospital = 1
  )
}

## A full, valid branching_process_main() argument set for end-to-end OBV tests.
## OBV parameters fall back to the function defaults (disabled) unless overridden
## via `...`.
obv_bpm_args <- function(...) {
  fixed <- function(v) function(n) rep(v, n)
  defaults <- list(
    mn_contacts_genPop           = 1.5,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 0.5,
    Tg_shape_genPop               = 2,
    Tg_rate_genPop                = 0.15,
    mn_contacts_hcw              = 1.5,
    baseline_risk_hcw              = 1,
    overdisp_contacts_hcw        = 0.5,
    Tg_shape_hcw                  = 2,
    Tg_rate_hcw                   = 0.15,
    mn_contacts_funeral          = 1.5,
    baseline_risk_funeral          = 1,
    overdisp_contacts_funeral    = 0.5,
    Tg_shape_funeral              = 10,
    Tg_rate_funeral               = 5,
    incubation_period             = fixed(5),
    onset_to_hospitalisation      = fixed(3),
    onset_to_death                = fixed(7),
    onset_to_recovery             = fixed(10),
    hospitalisation_to_death      = fixed(7),
    hospitalisation_to_recovery   = fixed(10),
    prob_symptomatic              = 0.9,
    prob_hospitalised_hcw         = 0.5,
    prob_hospitalised_genPop      = 0.5,
    prob_death_comm               = 0.5,
    prob_death_hosp               = 0.3,
    prob_hcw_cond_genPop_comm     = 0.01,
    prob_hcw_cond_genPop_hospital = 0.3,
    prob_hcw_cond_hcw_comm        = 0.01,
    prob_hcw_cond_hcw_hospital    = 0.3,
    prob_hospital_cond_hcw_preAdm = 0.5,
    ppe_coverage_hcw              = 0.5,
    ppe_efficacy                  = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 0.5,
    general_hospital_quarantine_efficacy = 0.5,
    p_unsafe_funeral_comm_hcw     = 0.5,
    p_unsafe_funeral_hosp_hcw     = 0.2,
    p_unsafe_funeral_comm_genPop  = 0.5,
    p_unsafe_funeral_hosp_genPop  = 0.2,
    safe_funeral_efficacy         = 1.0,
    prob_hcw_cond_funeral_hcw     = 0.05,
    prob_hcw_cond_funeral_genPop  = 0.05,
    population                    = 5000,
    hcw_per_capita                = 0.02,
    check_presymptomatic = FALSE,
    quiet                = TRUE,   # these runs sit at a tiny cap; the censoring warning is expected
    check_final_size              = 200,
    seeding_cases                 = 3,
    seed                          = 1L
  )
  utils::modifyList(defaults, list(...))
}

## A high-transmission, low-thinning, HCW-at-hospital-dominant scenario with full
## OBV prevention -- engineered so that the (default) HCW-at-hospital target
## yields prevented infections. Outbreak size is still stochastic, so callers that
## need prevented > 0 should drive it via run_with_prevention() below rather than
## trusting any single seed.
obv_full_prevention_args <- function(...) {
  obv_bpm_args(
    obv_pep_enabled               = TRUE,
    obv_pep_coverage              = 1,
    obv_pep_adherence             = 1,
    obv_pep_dpc                   = 0,
    obv_pep_efficacy              = 1,
    ## Remove hospital thinning so eligible exposures reach the gate, and make
    ## most hospital infections HCWs so they are eligible.
    ppe_coverage_hcw              = 0,
    etu_efficacy                  = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hospitalised_hcw         = 0.9,
    prob_hospitalised_genPop      = 0.9,
    prob_hcw_cond_genPop_hospital = 0.9,
    prob_hcw_cond_hcw_hospital    = 0.9,
    ## Seed the outbreak heavily so it reliably establishes (a 3-case seed can
    ## fizzle before any HCW-at-hospital transmission occurs).
    seeding_cases                 = 20,
    ...
  )
}

## --- extract_obv_prevented_info (deterministic, no RNG) ----------------

test_that("extract_obv_prevented_info selects kept-but-blocked candidates", {
  pre_thinning <- list(
    infection_location      = c("hospital", "hospital", "community", "hospital", "funeral"),
    offspring_class         = c("HCW", "HCW", "genPop", "HCW", "genPop"),
    infection_time_absolute = c(10, 11, 12, 13, 14)
  )
  ## Candidate 3 was thinned out upstream; the gate then blocked the 2nd and 4th
  ## *kept* candidates (i.e. original candidates 2 and 5).
  keep_infection <- c(TRUE, TRUE, FALSE, TRUE, TRUE)   # which() -> 1, 2, 4, 5
  keep_mask      <- c(TRUE, FALSE, TRUE, FALSE)        # aligned to kept candidates

  res <- extract_obv_prevented_info(pre_thinning, keep_infection, keep_mask)
  expect_equal(nrow(res), 2L)
  expect_equal(res$class, c("HCW", "genPop"))
  expect_equal(res$infection_location, c("hospital", "funeral"))
  expect_equal(res$time_infection_absolute, c(11, 14))
})

test_that("extract_obv_prevented_info returns a typed empty frame when nothing is blocked", {
  pre_thinning <- list(
    infection_location      = rep("hospital", 3),
    offspring_class         = rep("HCW", 3),
    infection_time_absolute = rep(5, 3)
  )
  res <- extract_obv_prevented_info(pre_thinning,
                                    keep_infection = c(TRUE, TRUE, TRUE),
                                    keep_mask      = c(TRUE, TRUE, TRUE))
  expect_equal(nrow(res), 0L)
  expect_named(res, c("class", "infection_location", "time_infection_absolute"))
})

## --- branching_process_main prevented-death counter --------------------

## Outbreak size -- hence the number of prevented infections -- is highly
## seed-dependent (a chain can fizzle before any HCW-at-hospital transmission).
## To keep the CFR invariants below non-trivial (prevented > 0) without pinning a
## fragile single seed, scan seeds and use the first run that actually prevents
## something. Productive seeds are common for this scenario, so this terminates
## almost immediately; if none prevents, fail loudly rather than asserting a
## vacuous 0 == 0.
run_with_prevention <- function(arg_fun, seeds = 1:40) {
  for (sd in seeds) {
    res <- do.call(branching_process_main, arg_fun(sd))
    s <- summarise_output(res$tdf, sim_info = res$sim_info)
    if (s$n_obv_pep_prevented > 0) return(s)
  }
  stop("no seed in 1:40 produced any prevented infections; check the scenario")
}

test_that("prevented_deaths equals prevented when everyone is symptomatic and CFR = 1", {
  ## prob_symptomatic = 1 and both CFRs = 1 => every infection that would have
  ## occurred dies, so every *prevented* infection is a prevented death, whatever
  ## the (stochastic) outbreak realisation.
  s <- run_with_prevention(function(sd)
    obv_full_prevention_args(prob_symptomatic = 1, prob_death_comm = 1,
                             prob_death_hosp = 1, seed = sd))
  expect_gt(s$n_obv_pep_prevented, 0)
  expect_equal(s$n_obv_pep_prevented_deaths, s$n_obv_pep_prevented)
})

test_that("prevented_deaths is 0 when CFR = 0 even though infections are prevented", {
  s <- run_with_prevention(function(sd)
    obv_full_prevention_args(prob_death_comm = 0, prob_death_hosp = 0, seed = sd))
  expect_gt(s$n_obv_pep_prevented, 0)
  expect_equal(s$n_obv_pep_prevented_deaths, 0)
})

test_that("prevented_deaths is bounded by prevented and reproducible under a fixed seed", {
  args <- obv_full_prevention_args(seed = 21L)
  r1 <- do.call(branching_process_main, args)
  r2 <- do.call(branching_process_main, args)
  s1 <- summarise_output(r1$tdf, sim_info = r1$sim_info)
  s2 <- summarise_output(r2$tdf, sim_info = r2$sim_info)

  expect_true(is.finite(s1$n_obv_pep_prevented_deaths))
  expect_gte(s1$n_obv_pep_prevented_deaths, 0)
  expect_lte(s1$n_obv_pep_prevented_deaths, s1$n_obv_pep_prevented)
  ## Deferred draw is reproducible (and, being post-loop, does not perturb the tree).
  expect_equal(s1$n_obv_pep_prevented_deaths, s2$n_obv_pep_prevented_deaths)
  expect_equal(s1$n_obv_pep_prevented, s2$n_obv_pep_prevented)
  expect_equal(r1$tdf, r2$tdf)
})

test_that("prevented_deaths is 0 when OBV is disabled", {
  res <- do.call(branching_process_main, obv_bpm_args(obv_pep_enabled = FALSE, seed = 21L))
  s <- summarise_output(res$tdf, sim_info = res$sim_info)
  expect_equal(s$n_obv_pep_prevented, 0)
  expect_equal(s$n_obv_pep_prevented_deaths, 0)
})

## --- branching_process_main prevented_completed output ------------------

## Like run_with_prevention(), but returns the full branching_process_main()
## result (not just the summary) so the prevented_completed frame can be
## inspected. Scans seeds for a run that actually prevents something.
run_res_with_prevention <- function(arg_fun, seeds = 1:40) {
  for (sd in seeds) {
    res <- do.call(branching_process_main, arg_fun(sd))
    if (summarise_output(res$tdf, sim_info = res$sim_info)$n_obv_pep_prevented > 0) {
      return(res)
    }
  }
  stop("no seed in 1:40 produced any prevented infections; check the scenario")
}

test_that("prevented_completed surfaces one counterfactual row per prevented infection", {
  res <- run_res_with_prevention(function(sd) obv_full_prevention_args(seed = sd))
  pc  <- res$prevented_completed
  s   <- summarise_output(res$tdf, sim_info = res$sim_info)

  expect_s3_class(pc, "data.frame")
  ## One row per averted index infection (the direct count, not onward chains).
  expect_equal(nrow(pc), s$n_obv_pep_prevented)
  ## The counterfactual infection-time column this output exists to expose.
  expect_true("time_infection_absolute" %in% names(pc))
  expect_true(all(is.finite(pc$time_infection_absolute)))
  ## The full counterfactual natural history is present (replayed outcome model).
  expect_true(all(c("class", "infection_location", "symptomatic", "hospitalisation",
                    "outcome", "outcome_location", "time_outcome_absolute") %in% names(pc)))
  ## Summing the would-be deaths reproduces the deferred prevented_deaths counter.
  expect_equal(sum(pc$outcome, na.rm = TRUE), s$n_obv_pep_prevented_deaths)
  ## Zero-time dummy-parent replay artifacts: parent/generation are NA and the
  ## "relative" infection time equals the absolute one.
  expect_true(all(is.na(pc$parent)))
  expect_true(all(is.na(pc$generation)))
  expect_equal(pc$time_infection_relative, pc$time_infection_absolute)

  ## Surfacing the frame draws no RNG, so it is reproducible under a fixed seed.
  res2 <- do.call(branching_process_main, obv_full_prevention_args(seed = res$sim_info$seed))
  expect_equal(res2$prevented_completed, pc)
})

test_that("prevented_completed is NULL when OBV prevents nothing", {
  res <- do.call(branching_process_main, obv_bpm_args(obv_pep_enabled = FALSE, seed = 21L))
  expect_null(res$prevented_completed)
})

## --- obv_pep_efficacy_from_dpc: flat default ---------------------------

test_that("obv_pep_efficacy_from_dpc is flat (constant E0) by default", {
  dpc <- c(0, 2, 5, 10, 15, 30, 100)
  expect_equal(obv_pep_efficacy_from_dpc(dpc, E0 = 0.7), rep(0.7, length(dpc)))
  ## shape = "flat" is the explicit default; cutoffs do not apply in flat mode.
  expect_equal(obv_pep_efficacy_from_dpc(dpc, E0 = 0.7, shape = "flat"), rep(0.7, length(dpc)))
  ## Default E0 at any DPC equals the dpc = 0 value (no decay).
  expect_equal(obv_pep_efficacy_from_dpc(c(0, 50)), rep(obv_pep_efficacy_from_dpc(0), 2))
})

test_that("obv_pep_efficacy_from_dpc rejects an unknown shape", {
  expect_error(obv_pep_efficacy_from_dpc(0, shape = "sigmoid"))
})

## --- obv_pep_efficacy_from_dpc: logistic boundaries (shape = "logistic") ---

test_that("obv_pep_efficacy_from_dpc returns E0 at dpc = 0", {
  ## True for both modes; the dpc = 0 value is the peak / the flat constant.
  expect_equal(obv_pep_efficacy_from_dpc(0), 0.82342697, tolerance = 1e-6)
  expect_equal(obv_pep_efficacy_from_dpc(0, shape = "logistic"), 0.82342697, tolerance = 1e-6)
})

test_that("obv_pep_efficacy_from_dpc (logistic) is 0 at dpc_zero = 15", {
  expect_equal(obv_pep_efficacy_from_dpc(15, shape = "logistic"), 0)
})

test_that("obv_pep_efficacy_from_dpc (logistic) is 0 just past max_dpc = 10", {
  expect_equal(obv_pep_efficacy_from_dpc(10.01, shape = "logistic"), 0)
})

test_that("obv_pep_efficacy_from_dpc (logistic) is in [0, 1] across a sensible range", {
  vals <- obv_pep_efficacy_from_dpc(seq(0, 20, by = 0.5), shape = "logistic")
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("obv_pep_efficacy_from_dpc (logistic) is monotone non-increasing on [0, 15]", {
  dpc <- seq(0, 15, by = 0.5)
  vals <- obv_pep_efficacy_from_dpc(dpc, shape = "logistic")
  expect_true(all(diff(vals) <= 1e-10))
})

## --- apply_obv_pep_gate no-op paths ------------------------------------

test_that("gate is a no-op when obv_pep_enabled = FALSE", {
  pre_thinning <- list(
    infection_location      = rep("hospital", 5),
    offspring_class         = rep("HCW", 5),
    infection_time_absolute = rep(10, 5)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:5,
    obv_pep_enabled = FALSE,
    obv_pep_coverage = 1,
    obv_pep_adherence = 1,
    obv_pep_efficacy = 1
  )
  expect_true(all(res$keep))
  expect_equal(res$num_treated$pre_eligible, 0L)
  expect_equal(res$num_treated$prevented, 0L)
  expect_true(all(!res$metadata$obv_pep_eligible))
  expect_true(all(is.na(res$metadata$obv_pep_dpc)))
})

test_that("gate is a no-op when obv_pep_coverage = 0", {
  pre_thinning <- list(
    infection_location      = rep("hospital", 5),
    offspring_class         = rep("HCW", 5),
    infection_time_absolute = rep(10, 5)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:5,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 0,
    obv_pep_adherence = 1,
    obv_pep_efficacy = 1
  )
  expect_true(all(res$keep))
  expect_equal(res$num_treated$pre_eligible, 5L)   ## eligibility is still recorded
  expect_equal(res$num_treated$pre_treated, 0L)
  expect_equal(res$num_treated$post_treated, 0L)
  expect_equal(res$num_treated$prevented, 0L)
})

## --- full-prevention regression ----------------------------------------

test_that("gate prevents all HCW-hospital candidates when cov = adh = eff = 1", {
  pre_thinning <- list(
    infection_location      = c("hospital", "hospital", "hospital", "community", "community"),
    offspring_class         = c("HCW", "HCW", "HCW", "HCW", "genPop"),
    infection_time_absolute = rep(10, 5)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:5,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 1,
    obv_pep_adherence = 1,
    obv_pep_dpc = 0,
    obv_pep_efficacy = 1,
    obv_pep_target_class = "HCW",
    obv_pep_target_locations = "hospital"
  )
  ## All three HCW-hospital rows must be flipped to FALSE (prevented).
  expect_equal(res$keep, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(res$num_treated$pre_eligible, 3L)
  expect_equal(res$num_treated$pre_treated,  3L)
  expect_equal(res$num_treated$pre_adherent, 3L)
  expect_equal(res$num_treated$post_eligible, 3L)
  expect_equal(res$num_treated$post_treated,  3L)
  expect_equal(res$num_treated$post_adherent, 3L)
  expect_equal(res$num_treated$prevented, 3L)
})

test_that("gate respects obv_pep_target_class (genPop community example)", {
  pre_thinning <- list(
    infection_location      = c("community", "community", "hospital", "hospital"),
    offspring_class         = c("genPop", "HCW", "genPop", "HCW"),
    infection_time_absolute = rep(10, 4)
  )
  set.seed(1)
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = 1:4,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 1,
    obv_pep_adherence = 1,
    obv_pep_dpc = 0,
    obv_pep_efficacy = 1,
    obv_pep_target_class = "genPop",
    obv_pep_target_locations = "community"
  )
  ## Only the genPop-community row (index 1) is eligible.
  expect_equal(res$num_treated$pre_eligible, 1L)
  expect_equal(res$keep, c(FALSE, TRUE, TRUE, TRUE))
})

## --- pre vs post nesting invariants ------------------------------------

test_that("pre/post counts are nested for any random draw", {
  set.seed(123)
  n_pre <- 200
  pre_thinning <- list(
    infection_location      = rep("hospital", n_pre),
    offspring_class         = sample(c("HCW", "genPop"), n_pre, TRUE),
    infection_time_absolute = rep(20, n_pre)
  )
  kept <- sort(sample.int(n_pre, n_pre %/% 2))
  res <- apply_obv_pep_gate(
    pre_thinning = pre_thinning,
    kept_indices = kept,
    obv_pep_enabled = TRUE,
    obv_pep_coverage = 0.6,
    obv_pep_adherence = 0.7,
    obv_pep_dpc = 2,
    obv_pep_efficacy = 0.5
  )
  nt <- res$num_treated
  ## Pre-thinning nesting
  expect_true(nt$pre_treated  <= nt$pre_eligible)
  expect_true(nt$pre_adherent <= nt$pre_treated)
  ## Post-thinning is a subset of pre-thinning
  expect_true(nt$post_eligible <= nt$pre_eligible)
  expect_true(nt$post_treated  <= nt$pre_treated)
  expect_true(nt$post_adherent <= nt$pre_adherent)
  ## Post-thinning is also nested with itself
  expect_true(nt$post_treated  <= nt$post_eligible)
  expect_true(nt$post_adherent <= nt$post_treated)
  ## prevented is bounded above by post_adherent (only adherent kept can be prevented)
  expect_true(nt$prevented <= nt$post_adherent)
  ## metadata length matches kept_indices length
  expect_equal(nrow(res$metadata), length(kept))
  expect_equal(length(res$keep), length(kept))
})

## --- integration through offspring_function_genPop ---------------------

test_that("offspring_function_genPop returns obv_pep_num_treated attribute (disabled = zero)", {
  parent <- make_parent_info_genPop_hospitalised()
  set.seed(7)
  out <- do.call(offspring_function_genPop,
                 c(list(parent_info = parent), obv_off_args()))
  nt <- attr(out, "obv_pep_num_treated", exact = TRUE)
  expect_false(is.null(nt))
  expect_true(all(unlist(nt) == 0))
})

test_that("offspring_function_genPop with full prevention zeroes out HCW hospital offspring", {
  parent <- make_parent_info_genPop_hospitalised()
  set.seed(7)
  args <- obv_off_args()
  args$obv_pep_enabled        <- TRUE
  args$obv_pep_coverage   <- 1
  args$obv_pep_adherence      <- 1
  args$obv_pep_dpc            <- 0
  args$obv_pep_efficacy       <- 1
  out <- do.call(offspring_function_genPop, c(list(parent_info = parent), args))
  ## No HCWs should remain at hospital in the realised offspring set.
  expect_equal(sum(out$class == "HCW" & out$infection_location == "hospital"), 0)
  ## genPop offspring at community / hospital are untouched.
  nt <- attr(out, "obv_pep_num_treated", exact = TRUE)
  expect_true(nt$pre_eligible >= nt$prevented)
  expect_true(nt$prevented == nt$post_adherent)
})

## --- summarise_output OBV field contract -------------------------------
## Regression guard: pin the set of OBV-related field names surfaced by
## summarise_output(). If anyone renames or drops these in the future,
## downstream analysis scripts will break -- this test catches it.

test_that("summarise_output() exposes the expected OBV PEP field names", {
  ## Build a tiny tdf that has just enough structure for summarise_output to run.
  ## We don't care about the simulation contents -- we're just checking output
  ## field names.
  tdf <- data.frame(
    id                            = 1L,
    class                         = "genPop",
    infection_location            = "community",
    parent                        = NA_integer_,
    generation                    = 1L,
    time_infection_relative       = 0,
    time_infection_absolute       = 0,
    incubation_period             = 5,
    symptomatic                   = TRUE,
    time_symptom_onset_relative   = 5,
    time_symptom_onset_absolute   = 5,
    hospitalisation               = FALSE,
    time_hospitalisation_relative = NA_real_,
    time_hospitalisation_absolute = NA_real_,
    outcome                       = TRUE,
    outcome_location              = "community",
    time_outcome_relative         = 15,
    time_outcome_absolute         = 15,
    funeral_safety                = "safe",
    obv_pep_eligible              = FALSE,
    obv_pep_received              = FALSE,
    obv_pep_adherent              = FALSE,
    obv_pep_dpc                   = NA_real_,
    n_offspring                   = 0L,
    offspring_generated           = TRUE,
    stringsAsFactors              = FALSE
  )
  attr(tdf, "obv_pep_num_treated") <- empty_obv_pep_num_treated()
  sim_info <- list(population = 100, hcw_per_capita = 0.05, hcw_total = 5,
                   obv_pep_enabled = FALSE,
                   obv_pep_num_treated = empty_obv_pep_num_treated())

  out <- summarise_output(tdf, sim_info = sim_info)

  expected_obv_fields <- c(
    ## seven gate counters
    "n_obv_pep_pre_eligible",
    "n_obv_pep_pre_treated",
    "n_obv_pep_pre_adherent",
    "n_obv_pep_post_eligible",
    "n_obv_pep_post_treated",
    "n_obv_pep_post_adherent",
    "n_obv_pep_prevented",
    ## deferred counterfactual: would-be deaths among prevented infections
    "n_obv_pep_prevented_deaths",
    ## three tdf-based cohort counters
    "n_obv_pep_eligible_cases",
    "n_obv_pep_treated_cases",
    "n_obv_pep_breakthroughs",
    ## derived ratio
    "prop_obv_pep_prevented_among_adherent"
  )
  missing_fields <- setdiff(expected_obv_fields, names(out))
  expect_equal(missing_fields, character(0),
               info = sprintf("summarise_output() missing fields: %s",
                              paste(missing_fields, collapse = ", ")))
})

## --- per-individual stochastic DPC (obv_pep_dpc_sd) ------------------
## obv_pep_dpc(t) is the MEAN days-post-challenge; obv_pep_dpc_sd turns the
## (previously deterministic) DPC into a per-recipient Gamma draw with that mean
## and a fixed standard deviation (a fixed absolute spread, so variance ~ sd^2
## regardless of the mean), modelling individual variation in how quickly the
## drug is received. The draw IS the realised (receipt - infection) delay. NULL
## (default) preserves the deterministic behaviour bit-for-bit.

## All-eligible HCW-hospital candidate set, fully covered & adherent, with
## efficacy = 0 so nothing is prevented -- every candidate is then a recipient
## whose drawn DPC surfaces in metadata$obv_pep_dpc, isolating the DPC draw.
dpc_gate_all_treated <- function(n, dpc, sd, seed = 1L,
                                 infection_time_absolute = 10) {
  pre_thinning <- list(
    infection_location      = rep("hospital", n),
    offspring_class         = rep("HCW", n),
    infection_time_absolute = rep(infection_time_absolute, n)
  )
  set.seed(seed)
  apply_obv_pep_gate(
    pre_thinning             = pre_thinning,
    kept_indices             = seq_len(n),
    obv_pep_enabled          = TRUE,
    obv_pep_coverage         = 1,
    obv_pep_adherence        = 1,
    obv_pep_dpc              = dpc,
    obv_pep_dpc_sd           = sd,
    obv_pep_efficacy         = 0,
    obv_pep_target_class     = "HCW",
    obv_pep_target_locations = "hospital"
  )
}

test_that("obv_pep_dpc_sd = NULL keeps DPC deterministic at the mean", {
  res <- dpc_gate_all_treated(n = 200, dpc = 3, sd = NULL)
  expect_true(all(res$metadata$obv_pep_received))
  ## Every recipient's DPC is exactly the mean -- no draw, no spread.
  expect_equal(res$metadata$obv_pep_dpc, rep(3, 200))
})

test_that("obv_pep_dpc_sd draws per-recipient DPC with the requested mean and variance", {
  res <- dpc_gate_all_treated(n = 20000, dpc = 4, sd = 2, seed = 42L)
  dpc <- res$metadata$obv_pep_dpc
  expect_length(dpc, 20000)
  expect_true(all(is.finite(dpc)))
  expect_true(all(dpc >= 0))                            # Gamma support is non-negative
  expect_false(isTRUE(all.equal(dpc, rep(4, 20000))))  # genuine individual variation
  ## Sample mean recovers obv_pep_dpc; sample variance ~ sd^2 = 4.
  expect_equal(mean(dpc), 4, tolerance = 0.05)
  expect_equal(var(dpc),  4, tolerance = 0.2)
})

test_that("smaller obv_pep_dpc_sd gives a tighter DPC spread (same mean)", {
  tight <- dpc_gate_all_treated(n = 20000, dpc = 5, sd = 1, seed = 7L)
  loose <- dpc_gate_all_treated(n = 20000, dpc = 5, sd = 4, seed = 7L)
  expect_lt(var(tight$metadata$obv_pep_dpc), var(loose$metadata$obv_pep_dpc))
  expect_equal(mean(tight$metadata$obv_pep_dpc), 5, tolerance = 0.05)
  expect_equal(mean(loose$metadata$obv_pep_dpc), 5, tolerance = 0.1)
  ## Variance tracks the requested sd^2 (1 vs 16), not the mean.
  expect_equal(var(tight$metadata$obv_pep_dpc), 1,  tolerance = 0.15)
  expect_equal(var(loose$metadata$obv_pep_dpc), 16, tolerance = 0.15)
})

test_that("stochastic DPC draw is reproducible under a fixed seed", {
  r1 <- dpc_gate_all_treated(n = 500, dpc = 3, sd = 2, seed = 99L)
  r2 <- dpc_gate_all_treated(n = 500, dpc = 3, sd = 2, seed = 99L)
  expect_equal(r1$metadata$obv_pep_dpc, r2$metadata$obv_pep_dpc)
})

test_that("a zero mean yields DPC 0 even with an sd set (degenerate point mass)", {
  res <- dpc_gate_all_treated(n = 100, dpc = 0, sd = 2)
  expect_equal(res$metadata$obv_pep_dpc, rep(0, 100))
})

test_that("a time-varying mean is honoured per-recipient under stochastic DPC", {
  ## Mean DPC ramps from 2 (early) to 12 (late) in calendar time; the drawn DPCs
  ## should track the mean at each recipient's own absolute infection time.
  mean_fn <- make_time_varying(times = c(0, 100), values = c(2, 12))
  pre_thinning <- list(
    infection_location      = rep("hospital", 8000),
    offspring_class         = rep("HCW", 8000),
    infection_time_absolute = c(rep(0, 4000), rep(100, 4000))
  )
  set.seed(3L)
  res <- apply_obv_pep_gate(
    pre_thinning             = pre_thinning,
    kept_indices             = seq_len(8000),
    obv_pep_enabled          = TRUE,
    obv_pep_coverage         = 1,
    obv_pep_adherence        = 1,
    obv_pep_dpc              = mean_fn,
    obv_pep_dpc_sd           = 2,
    obv_pep_efficacy         = 0
  )
  dpc <- res$metadata$obv_pep_dpc
  expect_equal(mean(dpc[1:4000]),    2,  tolerance = 0.1)
  expect_equal(mean(dpc[4001:8000]), 12, tolerance = 0.2)
})

test_that("obv_pep_dpc_sd validation rejects non-positive / non-scalar values", {
  base <- list(
    pre_thinning = list(
      infection_location      = rep("hospital", 3),
      offspring_class         = rep("HCW", 3),
      infection_time_absolute = rep(10, 3)
    ),
    kept_indices     = 1:3,
    obv_pep_enabled  = TRUE,
    obv_pep_coverage = 1,
    obv_pep_efficacy = 0
  )
  expect_error(do.call(apply_obv_pep_gate, c(base, list(obv_pep_dpc_sd = 0))),
               "must be NULL or a single finite positive numeric")
  expect_error(do.call(apply_obv_pep_gate, c(base, list(obv_pep_dpc_sd = -1))),
               "positive")
  expect_error(do.call(apply_obv_pep_gate, c(base, list(obv_pep_dpc_sd = c(1, 2)))),
               "single")
  expect_error(do.call(apply_obv_pep_gate, c(base, list(obv_pep_dpc_sd = Inf))),
               "finite")
  expect_error(do.call(apply_obv_pep_gate, c(base, list(obv_pep_dpc_sd = "a"))),
               "numeric")
  ## NULL (the default) is accepted.
  expect_silent(do.call(apply_obv_pep_gate, c(base, list(obv_pep_dpc_sd = NULL))))
})

test_that("stochastic DPC flows through branching_process_main: varying, finite, reproducible", {
  ## Treat all eligible HCW-hospital exposures (coverage 1) but prevent nothing
  ## (efficacy 0), so the outbreak proceeds and the recorded DPCs are pure draws.
  args <- obv_bpm_args(
    obv_pep_enabled                      = TRUE,
    obv_pep_coverage                     = 1,
    obv_pep_adherence                    = 1,
    obv_pep_dpc                          = 3,
    obv_pep_dpc_sd                       = 2,
    obv_pep_efficacy                     = 0,
    ppe_coverage_hcw                     = 0,
    etu_efficacy                         = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hospitalised_hcw                = 0.9,
    prob_hospitalised_genPop             = 0.9,
    prob_hcw_cond_genPop_hospital        = 0.9,
    prob_hcw_cond_hcw_hospital           = 0.9,
    seeding_cases                        = 20,
    seed                                 = 5L
  )
  r1  <- do.call(branching_process_main, args)
  r2  <- do.call(branching_process_main, args)
  rec <- r1$tdf$obv_pep_dpc[r1$tdf$obv_pep_received]

  expect_gt(length(rec), 5)            # the scenario does treat people
  expect_true(all(is.finite(rec)))
  expect_true(all(rec >= 0))
  expect_gt(var(rec), 0)               # genuine per-individual variation in DPC
  expect_equal(r1$tdf, r2$tdf)         # reproducible under a fixed seed

  ## Deterministic counterpart: same scenario, sd = NULL => DPC pinned at the mean.
  det <- do.call(branching_process_main, utils::modifyList(args, list(obv_pep_dpc_sd = NULL)))
  rec_det <- det$tdf$obv_pep_dpc[det$tdf$obv_pep_received]
  expect_gt(length(rec_det), 5)
  expect_equal(rec_det, rep(3, length(rec_det)))
})

## --- configurable efficacy curve (obv_pep_efficacy_args) ----------------
## obv_pep_efficacy_args is a named list of shape overrides for the built-in
## obv_pep_efficacy_from_dpc() curve (E0, d50, k, dpc_zero, max_dpc), applied
## when obv_pep_efficacy is NULL or a scalar (a scalar IS the curve's E0), so the
## curve's shape can be swept without writing a closure.

test_that("resolve_obv_efficacy applies obv_pep_efficacy_args to the built-in curve", {
  dpc <- c(0, 2, 5, 8, 12)
  ## Overrides are forwarded verbatim, including the shape switch that turns on decay.
  expect_equal(resolve_obv_efficacy(NULL, dpc, list(shape = "logistic", E0 = 0.6, d50 = 4)),
               obv_pep_efficacy_from_dpc(dpc, shape = "logistic", E0 = 0.6, d50 = 4))
  ## NULL or empty list reproduces the (flat) curve defaults exactly.
  expect_equal(resolve_obv_efficacy(NULL, dpc, NULL),   obv_pep_efficacy_from_dpc(dpc))
  expect_equal(resolve_obv_efficacy(NULL, dpc, list()), obv_pep_efficacy_from_dpc(dpc))
  ## Default (no shape) is flat, so the shape override genuinely changes the result.
  expect_false(isTRUE(all.equal(
    resolve_obv_efficacy(NULL, dpc, list(E0 = 0.6)),
    resolve_obv_efficacy(NULL, dpc, list(shape = "logistic", E0 = 0.6)))))
})

test_that("obv_pep_efficacy_args validation catches misuse", {
  expect_error(validate_obv_efficacy_args(NULL, list(E0 = 0.9, foo = 1)),
               "unknown argument")
  expect_error(validate_obv_efficacy_args(NULL, list(0.9)),            # unnamed element
               "fully named")
  expect_error(validate_obv_efficacy_args(NULL, list(E0 = 1, E0 = 2)), # duplicate names
               "duplicate")
  expect_error(validate_obv_efficacy_args(NULL, 5),                    # not a list
               "named list")
  ## Curve overrides apply to the built-in curve only, not a custom function.
  expect_error(validate_obv_efficacy_args(function(d) d, list(d50 = 4)),
               "function")
  ## A scalar obv_pep_efficacy already supplies E0, so naming E0 in args conflicts.
  expect_error(validate_obv_efficacy_args(0.5, list(E0 = 0.9)),
               "E0")
  ## Valid combinations are silent.
  expect_silent(validate_obv_efficacy_args(NULL, NULL))
  expect_silent(validate_obv_efficacy_args(NULL, list()))
  expect_silent(validate_obv_efficacy_args(NULL, list(E0 = 0.9, max_dpc = 14)))
  expect_silent(validate_obv_efficacy_args(0.5, list(d50 = 4, k = 2)))   # scalar E0 + shape OK
  expect_silent(validate_obv_efficacy_args(function(d) rep(0.5, length(d)), NULL))
})

test_that("bad curve values surface through resolve_obv_efficacy", {
  ## Value validation is delegated to obv_pep_efficacy_from_dpc (E0 in [0, 1]).
  expect_error(resolve_obv_efficacy(NULL, 0, list(E0 = 1.5)), "E0")
})

test_that("scalar obv_pep_efficacy is used as the curve's E0 (flat by default)", {
  dpc <- c(0, 2, 5, 8, 12)
  ## A bare scalar is E0 of the built-in curve, which is flat by default:
  ## efficacy = the scalar at every DPC.
  expect_equal(resolve_obv_efficacy(0.6, dpc), rep(0.6, length(dpc)))
  expect_equal(resolve_obv_efficacy(0.6, dpc), obv_pep_efficacy_from_dpc(dpc, E0 = 0.6))
  ## With shape = "logistic" the scalar E0 feeds the decaying curve.
  expect_equal(resolve_obv_efficacy(0.6, dpc, list(shape = "logistic", d50 = 4, k = 2)),
               obv_pep_efficacy_from_dpc(dpc, E0 = 0.6, shape = "logistic", d50 = 4, k = 2))
  ## E0 is a pure vertical scale: halving E0 halves the whole curve (logistic too).
  expect_equal(resolve_obv_efficacy(0.5, dpc, list(shape = "logistic")),
               0.5 * resolve_obv_efficacy(1, dpc, list(shape = "logistic")))
  ## A bad scalar E0 is rejected by the curve's own validation.
  expect_error(resolve_obv_efficacy(1.5, 0), "E0")
})

test_that("scalar obv_pep_efficacy as E0 drives prevention through the gate", {
  pre_thinning <- list(
    infection_location      = rep("hospital", 8),
    offspring_class         = rep("HCW", 8),
    infection_time_absolute = rep(10, 8)
  )
  ## dpc = 0 so efficacy = E0 exactly: scalar 1 prevents all, scalar 0 none.
  run_gate <- function(eff, curve_args = NULL) {
    set.seed(1)
    apply_obv_pep_gate(
      pre_thinning = pre_thinning, kept_indices = 1:8,
      obv_pep_enabled = TRUE, obv_pep_coverage = 1, obv_pep_adherence = 1,
      obv_pep_dpc = 0, obv_pep_efficacy = eff, obv_pep_efficacy_args = curve_args
    )
  }
  expect_equal(run_gate(1)$num_treated$prevented, 8L)
  expect_equal(run_gate(0)$num_treated$prevented, 0L)
  ## Scalar E0 plus a shape override is accepted and reaches the curve.
  expect_silent(run_gate(0.5, list(d50 = 4)))
})

test_that("obv_pep_efficacy_args overrides the curve in the gate (E0 = 1 prevents all, E0 = 0 none)", {
  pre_thinning <- list(
    infection_location      = rep("hospital", 6),
    offspring_class         = rep("HCW", 6),
    infection_time_absolute = rep(10, 6)
  )
  run_gate <- function(curve_args) {
    set.seed(1)
    apply_obv_pep_gate(
      pre_thinning = pre_thinning, kept_indices = 1:6,
      obv_pep_enabled = TRUE, obv_pep_coverage = 1, obv_pep_adherence = 1,
      obv_pep_dpc = 0, obv_pep_efficacy = NULL, obv_pep_efficacy_args = curve_args
    )
  }
  ## At dpc = 0, efficacy = E0, so E0 = 1 prevents every adherent candidate and
  ## E0 = 0 prevents none -- a clean, deterministic check that overrides land.
  all_prevented  <- run_gate(list(E0 = 1))
  none_prevented <- run_gate(list(E0 = 0))
  expect_true(all(!all_prevented$keep))
  expect_equal(all_prevented$num_treated$prevented, 6L)
  expect_true(all(none_prevented$keep))
  expect_equal(none_prevented$num_treated$prevented, 0L)
})

test_that("branching_process_main fails fast on an out-of-range efficacy curve override", {
  expect_error(
    do.call(branching_process_main,
            obv_bpm_args(obv_pep_enabled = TRUE, obv_pep_efficacy_args = list(E0 = 1.5), seed = 1L)),
    "E0"
  )
})

test_that("obv_pep_efficacy_args flows through branching_process_main (E0 = 0 treats but prevents nothing)", {
  args <- obv_bpm_args(
    obv_pep_enabled                      = TRUE,
    obv_pep_coverage                     = 1,
    obv_pep_adherence                    = 1,
    obv_pep_efficacy_args                = list(E0 = 0),
    ppe_coverage_hcw                     = 0,
    etu_efficacy                         = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hospitalised_hcw                = 0.9,
    prob_hospitalised_genPop             = 0.9,
    prob_hcw_cond_genPop_hospital        = 0.9,
    prob_hcw_cond_hcw_hospital           = 0.9,
    seeding_cases                        = 20,
    seed                                 = 5L
  )
  res <- do.call(branching_process_main, args)
  s   <- summarise_output(res$tdf, sim_info = res$sim_info)
  expect_gt(s$n_obv_pep_pre_treated, 5)   # people are treated (coverage 1)
  expect_equal(s$n_obv_pep_prevented, 0)  # but E0 = 0 => zero efficacy => nothing prevented
})
