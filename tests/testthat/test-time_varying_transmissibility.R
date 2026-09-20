## Tests for time-varying mean-offspring transmissibility:
##   mn_contacts_genPop, mn_contacts_hcw, mn_contacts_funeral.
##
## These may be supplied as a single positive scalar (back-compatible) or as a
## function(t) of absolute calendar time. genPop/HCW are resolved at the
## parent's infection time; funeral is resolved at the parent's death time.

## --- helpers --------------------------------------------------------------

## A parent_info row with the columns the offspring functions read. Defaults to
## a non-hospitalised genPop case that dies (so funeral transmission can fire and
## genPop community offspring are not subject to hospital thinning).
make_parent_info <- function(time_infection_absolute = 0,
                             time_outcome_relative   = 12,
                             class                   = "genPop",
                             hospitalised            = FALSE,
                             outcome                 = TRUE,
                             funeral_safety          = "unsafe") {
  data.frame(
    id                             = 1L,
    class                          = class,
    infection_location             = "community",
    parent                         = NA_integer_,
    generation                     = 1L,
    time_infection_relative        = 0,
    time_infection_absolute        = time_infection_absolute,
    incubation_period              = 5,
    symptomatic                    = TRUE,
    time_symptom_onset_relative    = 5,
    time_symptom_onset_absolute    = time_infection_absolute + 5,
    hospitalisation                = hospitalised,
    time_hospitalisation_relative  = if (hospitalised) 3 else NA_real_,
    time_hospitalisation_absolute  = if (hospitalised) time_infection_absolute + 3 else NA_real_,
    outcome                        = outcome,
    outcome_location               = "community",
    time_outcome_relative          = time_outcome_relative,
    time_outcome_absolute          = time_infection_absolute + time_outcome_relative,
    funeral_safety                 = funeral_safety,
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors               = FALSE
  )
}

## Mean realised offspring count over `reps` independent calls (different seeds).
mean_offspring_genPop <- function(reps, parent, mn, base_seed = 0) {
  counts <- vapply(seq_len(reps), function(i) {
    set.seed(base_seed + i)
    nrow(offspring_function_genPop(
      parent_info                   = parent,
      mn_contacts_genPop           = mn,
      baseline_risk_genPop           = 1,
      overdisp_contacts_genPop     = 5,
      Tg_shape_genPop               = 2,
      Tg_rate_genPop                = 0.15,
      prop_etu                      = 1,
      etu_efficacy                  = 0.5,
      general_hospital_quarantine_efficacy = 0.5,
      ppe_coverage_hcw              = 0.5,
      ppe_efficacy                  = 1,
      prob_hcw_cond_genPop_comm     = 0,
      prob_hcw_cond_genPop_hospital = 0.3
    ))
  }, numeric(1))
  mean(counts)
}

mean_offspring_funeral <- function(reps, parent, mn, base_seed = 0) {
  counts <- vapply(seq_len(reps), function(i) {
    set.seed(base_seed + i)
    nrow(offspring_function_funeral(
      parent_info                  = parent,
      mn_contacts_funeral         = mn,
      baseline_risk_funeral         = 1,
      overdisp_contacts_funeral   = 5,
      Tg_shape_funeral             = 10,
      Tg_rate_funeral              = 5,
      safe_funeral_efficacy        = 1.0,
      prob_hcw_cond_funeral_hcw    = 0,
      prob_hcw_cond_funeral_genPop = 0
    ))
  }, numeric(1))
  mean(counts)
}

## A full, valid branching_process_main() argument set (scalars by default).
fixed_inc   <- function(n) rep(5, n)
fixed_hosp  <- function(n) rep(3, n)
fixed_death <- function(n) rep(7, n)
fixed_recov <- function(n) rep(10, n)

bpm_args <- function(...) {
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
    check_final_size              = 200,
    seeding_cases                 = 3,
    seed                          = 1L
  )
  utils::modifyList(defaults, list(...))
}

## --- resolve_positive_time_varying unit tests -----------------------------

test_that("resolve_positive_time_varying handles scalars and functions", {
  expect_equal(resolve_positive_time_varying(5, 10, "x"), 5)
  expect_equal(resolve_positive_time_varying(function(t) t + 1, 4, "x"), 5)
  expect_equal(resolve_positive_time_varying(function(t) rep(2, length(t)), c(0, 5, 9), "x"),
               c(2, 2, 2))
})

test_that("resolve_positive_time_varying rejects non-positive values", {
  expect_error(resolve_positive_time_varying(0, 0, "x"), "positive")
  expect_error(resolve_positive_time_varying(-1, 0, "x"), "positive")
  expect_error(resolve_positive_time_varying(function(t) -1, 5, "x"), "positive")
})

## --- Scalar back-compat & constant-function equivalence -------------------

test_that("offspring_function_genPop accepts a scalar mn (back-compatible)", {
  set.seed(1)
  res <- offspring_function_genPop(
    parent_info                   = make_parent_info(),
    mn_contacts_genPop           = 8,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 5,
    Tg_shape_genPop               = 2,
    Tg_rate_genPop                = 0.15,
    prop_etu                      = 1,
    etu_efficacy                  = 0.5,
    general_hospital_quarantine_efficacy = 0.5,
    ppe_coverage_hcw              = 0.5,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.3
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("infection_location", "time_infection_relative", "class") %in% names(res)))
})

test_that("a constant function(t) reproduces the scalar exactly (genPop)", {
  parent <- make_parent_info()

  set.seed(123)
  res_scalar <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = 8,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 5,
    Tg_shape_genPop               = 2,
    Tg_rate_genPop                = 0.15,
    prop_etu                      = 1,
    etu_efficacy                  = 0.5,
    general_hospital_quarantine_efficacy = 0.5,
    ppe_coverage_hcw              = 0.5,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.3
  )

  set.seed(123)
  res_fn <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = function(t) 8,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 5,
    Tg_shape_genPop               = 2,
    Tg_rate_genPop                = 0.15,
    prop_etu                      = 1,
    etu_efficacy                  = 0.5,
    general_hospital_quarantine_efficacy = 0.5,
    ppe_coverage_hcw              = 0.5,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.3
  )

  expect_identical(res_scalar$time_infection_relative, res_fn$time_infection_relative)
  expect_identical(res_scalar$class, res_fn$class)
})

## --- Time-variation actually changes offspring counts ---------------------

test_that("time-varying mn_contacts_genPop is keyed to the parent's infection time", {
  ## High transmissibility early (t < 100), low late.
  mn_fn <- function(t) ifelse(t < 100, 40, 2)

  parent_early <- make_parent_info(time_infection_absolute = 0)    # mn(0)   = 40
  parent_late  <- make_parent_info(time_infection_absolute = 200)  # mn(200) = 2

  mean_early <- mean_offspring_genPop(reps = 60, parent = parent_early, mn = mn_fn)
  mean_late  <- mean_offspring_genPop(reps = 60, parent = parent_late,  mn = mn_fn)

  expect_gt(mean_early, 25)   # near 40
  expect_lt(mean_late,  10)   # near 2
  expect_gt(mean_early, mean_late)
})

test_that("time-varying mn_contacts_funeral is keyed to the parent's DEATH time", {
  ## Low early (t < 100), high late. Both parents are infected at t = 0, so if
  ## the funeral mean were (wrongly) keyed to infection time both would be low.
  ## Keyed to death time, the late-dying parent should generate many more.
  mn_fn <- function(t) ifelse(t < 100, 2, 40)

  parent_dies_late  <- make_parent_info(time_infection_absolute = 0,
                                        time_outcome_relative = 150)  # death @150 -> 40
  parent_dies_early <- make_parent_info(time_infection_absolute = 0,
                                        time_outcome_relative = 10)   # death @10  -> 2

  mean_late  <- mean_offspring_funeral(reps = 60, parent = parent_dies_late,  mn = mn_fn)
  mean_early <- mean_offspring_funeral(reps = 60, parent = parent_dies_early, mn = mn_fn)

  expect_gt(mean_late, 25)   # near 40, proving death-time anchoring
  expect_lt(mean_early, 10)  # near 2
  expect_gt(mean_late, mean_early)
})

## --- Validation -----------------------------------------------------------

test_that("offspring functions reject a non-positive scalar mn", {
  expect_error(
    offspring_function_genPop(
      parent_info                   = make_parent_info(),
      mn_contacts_genPop           = -1,
      baseline_risk_genPop           = 1,
      overdisp_contacts_genPop     = 5,
      Tg_shape_genPop               = 2,
      Tg_rate_genPop                = 0.15,
      prop_etu                      = 1,
      etu_efficacy                  = 0.5,
      general_hospital_quarantine_efficacy = 0.5,
      ppe_coverage_hcw              = 0.5,
      ppe_efficacy                  = 1,
      prob_hcw_cond_genPop_comm     = 0,
      prob_hcw_cond_genPop_hospital = 0.3
    ),
    "mn_contacts_genPop"
  )
})

test_that("branching_process_main rejects a mn curve that dips to <= 0", {
  bad <- make_time_varying(c(0, 30, 60), c(2, 1, -0.5))
  expect_error(
    do.call(branching_process_main, bpm_args(mn_contacts_genPop = bad)),
    "mn_contacts_genPop.*positive"
  )
})

## --- End-to-end integration ----------------------------------------------

test_that("branching_process_main runs with time-varying mn_contacts_*", {
  args <- bpm_args(
    mn_contacts_genPop  = make_time_varying(c(0, 60, 120), c(2.0, 1.2, 0.6)),
    baseline_risk_genPop  = 1,
    mn_contacts_hcw     = function(t) 1.5,
    baseline_risk_hcw     = 1,
    mn_contacts_funeral = make_time_varying(c(0, 60), c(2.0, 0.5))
  )
  out <- do.call(branching_process_main, args)
  expect_true(is.list(out))
  expect_true(nrow(out$tdf) > 0)
})

test_that("a constant function reproduces the scalar run exactly (full sim)", {
  out_scalar <- do.call(branching_process_main, bpm_args(mn_contacts_genPop = 1.5))
  out_fn     <- do.call(branching_process_main, bpm_args(mn_contacts_genPop = function(t) 1.5))
  ## Same seed + same resolved mean + no extra RNG draws => identical tree.
  expect_identical(out_scalar$tdf$time_infection_absolute,
                   out_fn$tdf$time_infection_absolute)
  expect_identical(out_scalar$tdf$class, out_fn$tdf$class)
  expect_identical(nrow(out_scalar$tdf), nrow(out_fn$tdf))
})
