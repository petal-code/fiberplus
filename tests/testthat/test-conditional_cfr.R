## Tests for the Option B (conditional CFR) interpretation of
## prob_death_comm / prob_death_hosp / prob_hospitalised_X, the upfront
## sanity check on time-varying probability parameters, and the
## recycling warning in resolve_time_varying.

## --- helpers --------------------------------------------------------

make_parent_info <- function(time_infection_absolute = 0) {
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
    hospitalisation                = FALSE,
    time_hospitalisation_relative  = NA_real_,
    time_hospitalisation_absolute  = NA_real_,
    outcome                        = TRUE,
    outcome_location               = "community",
    time_outcome_relative          = 12,
    time_outcome_absolute          = time_infection_absolute + 12,
    funeral_safety                 = "unsafe",
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors               = FALSE
  )
}

make_offspring_df <- function(n, classes = rep("genPop", n)) {
  data.frame(
    id                             = rep(NA_integer_, n),
    class                          = classes,
    infection_location             = rep("community", n),
    parent                         = rep(NA_integer_, n),
    generation                     = rep(NA_integer_, n),
    time_infection_relative        = rep(1, n),
    time_infection_absolute        = rep(NA_real_, n),
    incubation_period              = rep(NA_real_, n),
    symptomatic                    = rep(NA, n),
    time_symptom_onset_relative    = rep(NA_real_, n),
    time_symptom_onset_absolute    = rep(NA_real_, n),
    hospitalisation                = rep(FALSE, n),
    time_hospitalisation_relative  = rep(NA_real_, n),
    time_hospitalisation_absolute  = rep(NA_real_, n),
    outcome                        = rep(FALSE, n),
    outcome_location               = rep(NA_character_, n),
    time_outcome_relative          = rep(NA_real_, n),
    time_outcome_absolute          = rep(NA_real_, n),
    funeral_safety                 = rep(NA_character_, n),
    n_offspring                    = rep(NA_integer_, n),
    offspring_generated            = rep(FALSE, n),
    stringsAsFactors               = FALSE
  )
}

fixed_inc   <- function(n) rep(5, n)
fixed_hosp  <- function(n) rep(2, n)
fixed_death <- function(n) rep(7, n)
fixed_recov <- function(n) rep(10, n)

base_args <- function(...) {
  defaults <- list(
    prob_symptomatic             = 0.5,
    prob_hospitalised_hcw        = 0,
    prob_hospitalised_genPop     = 0,
    prob_death_comm              = 0.4,
    prob_death_hosp              = 0.2,
    p_unsafe_funeral_comm_hcw    = 0.5,
    p_unsafe_funeral_hosp_hcw    = 0.2,
    p_unsafe_funeral_comm_genPop = 0.5,
    p_unsafe_funeral_hosp_genPop = 0.2,
    incubation_period            = fixed_inc,
    onset_to_hospitalisation     = fixed_hosp,
    hospitalisation_to_death     = fixed_death,
    hospitalisation_to_recovery  = fixed_recov,
    onset_to_death               = fixed_death,
    onset_to_recovery            = fixed_recov
  )
  utils::modifyList(defaults, list(...))
}

## --- Option B: symptomatic-death rate matches prob_death_comm ----------

## With prob_symptomatic = 0.5 and prob_death_comm = 0.4, under Option B
## the fraction of symptomatic non-hospitalised offspring who die should be
## ~0.4 (not 0.8 as under the old Option-A implementation).
test_that("Option B: symptomatic offspring in community die at prob_death_comm rate", {
  set.seed(2026)
  parent    <- make_parent_info()
  offspring <- make_offspring_df(2000)
  args <- base_args(prob_symptomatic = 0.5, prob_death_comm = 0.4,
                    prob_death_hosp  = 0.2,
                    prob_hospitalised_hcw = 0, prob_hospitalised_genPop = 0)

  result <- do.call(complete_offspring_info,
                    c(list(parent_info = parent, offspring_dataframe = offspring),
                      args))

  symp <- result[result$symptomatic, ]
  expect_true(nrow(symp) > 800,
              info = "need a reasonable symptomatic sample")
  death_rate <- mean(symp$outcome)
  expect_true(abs(death_rate - 0.4) < 0.04,
              info = sprintf("P(die | symp, comm) = %.3f; expected ~0.4 under Option B", death_rate))
})

## --- Option B: realised hospital CFR equals prob_death_hosp ------------

## With prob_symptomatic = 0.5, prob_death_comm = 0.8, prob_death_hosp = 0.3,
## hospitalised symptomatic offspring should die at rate ~0.3 (= prob_death_hosp).
## Under the old Option-A implementation this would resolve to 0.3 / 0.5 = 0.6.
test_that("Option B: realised hospital CFR equals prob_death_hosp", {
  set.seed(2027)
  parent    <- make_parent_info()
  offspring <- make_offspring_df(3000)
  args <- base_args(prob_symptomatic = 0.5,
                    prob_death_comm  = 0.8, prob_death_hosp = 0.3,
                    prob_hospitalised_hcw = 0.9, prob_hospitalised_genPop = 0.9)

  result <- do.call(complete_offspring_info,
                    c(list(parent_info = parent, offspring_dataframe = offspring),
                      args))

  hosp <- result[result$hospitalisation, ]
  expect_true(nrow(hosp) > 500,
              info = "need a reasonable hospitalised sample")
  hosp_cfr <- mean(hosp$outcome)
  expect_true(abs(hosp_cfr - 0.3) < 0.04,
              info = sprintf("P(die | hosp) = %.3f; expected ~0.3 under Option B", hosp_cfr))
})

## --- Option B: hospitalisation rate matches prob_hospitalised_X --------

## With prob_symptomatic = 0.5 and prob_hospitalised_genPop = 0.6, the fraction
## of symptomatic offspring who get hospitalised should be ~0.6 directly.
test_that("Option B: prob_hospitalised_X is P(hosp | symptomatic) directly", {
  set.seed(2028)
  parent    <- make_parent_info()
  offspring <- make_offspring_df(2000)
  args <- base_args(prob_symptomatic = 0.5, prob_hospitalised_genPop = 0.6,
                    prob_death_comm = 0.1, prob_death_hosp = 0.05)

  result <- do.call(complete_offspring_info,
                    c(list(parent_info = parent, offspring_dataframe = offspring),
                      args))

  symp <- result[result$symptomatic, ]
  hosp_rate <- mean(symp$hospitalisation)
  expect_true(abs(hosp_rate - 0.6) < 0.04,
              info = sprintf("P(hosp | symp) = %.3f; expected ~0.6 under Option B", hosp_rate))
})

## --- Asymptomatic offspring still cannot die ---------------------------

test_that("Asymptomatic offspring never die regardless of parameters", {
  set.seed(2029)
  parent    <- make_parent_info()
  offspring <- make_offspring_df(500)
  args <- base_args(prob_symptomatic = 0.3,
                    prob_death_comm  = 0.9, prob_death_hosp = 0.5,
                    prob_hospitalised_hcw = 0.5, prob_hospitalised_genPop = 0.5)

  result <- do.call(complete_offspring_info,
                    c(list(parent_info = parent, offspring_dataframe = offspring),
                      args))

  asymp <- result[!result$symptomatic, ]
  expect_true(nrow(asymp) > 0)
  expect_equal(sum(asymp$outcome), 0,
               info = "no asymptomatic offspring should die")
  expect_equal(sum(asymp$hospitalisation), 0,
               info = "no asymptomatic offspring should be hospitalised")
})

## --- resolve_time_varying recycling warning ---------------------------

test_that("resolve_time_varying warns when a function returns length 1 for length>1 t", {
  expect_warning(
    resolve_time_varying(function(t) 0.5, c(0, 10, 20), "test_param"),
    "test_param.*length-1 value"
  )
})

test_that("resolve_time_varying does NOT warn for scalar input recycled over vector t", {
  expect_silent(
    resolve_time_varying(0.5, c(0, 10, 20), "test_param")
  )
})

test_that("resolve_time_varying does NOT warn for vectorised function returning length(t)", {
  expect_silent(
    resolve_time_varying(function(t) rep(0.5, length(t)), c(0, 10, 20), "test_param")
  )
})

## --- Upfront sanity check in branching_process_main -------------------
##
## We need to actually call branching_process_main with otherwise-valid
## arguments to exercise the upfront sanity block.

bpm_args <- function(...) {
  defaults <- list(
    mn_contacts_genPop           = 1.0,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 0.5,
    Tg_shape_genPop               = 2,
    Tg_rate_genPop                = 0.15,
    mn_contacts_hcw              = 1.0,
    baseline_risk_hcw              = 1,
    overdisp_contacts_hcw        = 0.5,
    Tg_shape_hcw                  = 2,
    Tg_rate_hcw                   = 0.15,
    mn_contacts_funeral          = 1.0,
    baseline_risk_funeral          = 1,
    overdisp_contacts_funeral    = 0.5,
    Tg_shape_funeral              = 10,
    Tg_rate_funeral               = 5,
    incubation_period             = fixed_inc,
    onset_to_hospitalisation      = fixed_hosp,
    onset_to_death                = fixed_death,
    onset_to_recovery             = fixed_recov,
    hospitalisation_to_death      = fixed_death,
    hospitalisation_to_recovery   = fixed_recov,
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
    check_final_size              = 50,
    seeding_cases                 = 2,
    seed                          = 1L
  )
  utils::modifyList(defaults, list(...))
}

test_that("Upfront sanity check rejects prob_death_hosp > prob_death_comm", {
  args <- bpm_args(prob_death_comm = 0.3, prob_death_hosp = 0.6)
  expect_error(
    do.call(branching_process_main, args),
    "prob_death_hosp.*prob_death_comm"
  )
})

test_that("Upfront sanity check rejects time-varying prob exceeding 1", {
  bad <- make_time_varying(c(0, 30, 60), c(0.5, 0.9, 1.2))
  args <- bpm_args(prob_hospitalised_genPop = bad)
  expect_error(
    do.call(branching_process_main, args),
    "prob_hospitalised_genPop.*\\[0, 1\\]"
  )
})

test_that("Upfront sanity check rejects time-varying prob below 0", {
  bad <- make_time_varying(c(0, 30, 60), c(0.5, 0.2, -0.1))
  args <- bpm_args(p_unsafe_funeral_comm_genPop = bad)
  expect_error(
    do.call(branching_process_main, args),
    "p_unsafe_funeral_comm_genPop.*\\[0, 1\\]"
  )
})

test_that("Upfront sanity check rejects non-positive hospitalisation_delay_factor", {
  bad <- make_time_varying(c(0, 30, 60), c(1.0, 0.5, -0.1))
  args <- bpm_args(hospitalisation_delay_factor = bad)
  expect_error(
    do.call(branching_process_main, args),
    "hospitalisation_delay_factor.*positive"
  )
})

test_that("Upfront sanity check accepts a valid time-varying curve", {
  good <- make_time_varying(c(0, 30, 60, 90), c(0.1, 0.3, 0.5, 0.7))
  args <- bpm_args(prob_hospitalised_genPop = good)
  ## Should run to completion without error. The simulation itself may produce
  ## a small outbreak; we just check it didn't error in the upfront block.
  expect_silent({
    out <- suppressWarnings(do.call(branching_process_main, args))
  })
  expect_true(is.list(out))
  expect_true("tdf" %in% names(out))
})
