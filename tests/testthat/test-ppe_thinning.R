## Tests for PPE / hospital-quarantine thinning in the offspring functions.
## Covers the fixes to Issues A (genPop -> HCW at hospital) and B (HCW
## post-admission -> HCW receiver) and pins the multiplicative (Option B)
## PPE semantics in offspring_function_hcw.

## --- helpers -----------------------------------------------------------

make_parent_info_genPop_hospitalised <- function(time_infection_absolute = 0,
                                                  time_to_hospitalisation = 3,
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

make_parent_info_HCW <- function(time_infection_absolute = 0,
                                 hospitalised = TRUE,
                                 time_to_hospitalisation = 3,
                                 time_to_outcome = 30) {
  data.frame(
    id                             = 1L,
    class                          = "HCW",
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
    time_hospitalisation_relative  = if (hospitalised) time_to_hospitalisation else NA_real_,
    time_hospitalisation_absolute  = if (hospitalised) time_infection_absolute + time_to_hospitalisation else NA_real_,
    outcome                        = TRUE,
    outcome_location               = if (hospitalised) "hospital" else "community",
    time_outcome_relative          = time_to_outcome,
    time_outcome_absolute          = time_infection_absolute + time_to_outcome,
    funeral_safety                 = "unsafe",
    n_offspring                    = NA_integer_,
    offspring_generated            = FALSE,
    stringsAsFactors               = FALSE
  )
}

## --- Issue A regression ------------------------------------------------
## genPop parent hospitalised, ppe = 1, hq = 0: HCW recipients in hospital must
## all be thinned by PPE, while genPop recipients in hospital pass through.
test_that("genPop parent: PPE thins HCW recipients at hospital (Issue A)", {
  parent <- make_parent_info_genPop_hospitalised(time_to_hospitalisation = 0.001)
  set.seed(42)
  result <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = 500,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 200,
    Tg_shape_genPop               = 4,
    Tg_rate_genPop                = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 0,
    general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw              = 1,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.5
  )

  hospital_offspring <- result[result$infection_location == "hospital", ]
  expect_true(nrow(hospital_offspring) > 0,
              info = "expected hospital offspring with hq = 0")
  expect_equal(sum(hospital_offspring$class == "HCW"), 0,
               info = "ppe = 1 should remove every HCW hospital recipient")
  expect_true(sum(hospital_offspring$class == "genPop") > 0,
              info = "genPop hospital recipients should not be PPE-thinned")
})

## hq = 1, ppe = 0: all hospital recipients regardless of class are thinned.
test_that("genPop parent: full hospital quarantine (hq = 1) zeroes all hospital events", {
  parent <- make_parent_info_genPop_hospitalised(time_to_hospitalisation = 0.001)
  set.seed(7)
  result <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = 500,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 200,
    Tg_shape_genPop               = 4,
    Tg_rate_genPop                = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 1,
    general_hospital_quarantine_efficacy = 1,
    ppe_coverage_hcw              = 0,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.5
  )
  expect_equal(sum(result$infection_location == "hospital"), 0)
})

## --- Issue B regression ------------------------------------------------
## HCW parent post-admission, ppe = 1, hq = 0: HCW recipients in hospital must
## all be thinned by receiver PPE.
test_that("HCW post-admission parent: receiver PPE thins HCW recipients (Issue B)", {
  parent <- make_parent_info_HCW(hospitalised = TRUE,
                                 time_to_hospitalisation = 0.001,
                                 time_to_outcome = 30)
  set.seed(99)
  result <- offspring_function_hcw(
    parent_info                     = parent,
    mn_contacts_hcw                = 500,
    baseline_risk_hcw                = 1,
    overdisp_contacts_hcw          = 200,
    Tg_shape_hcw                    = 4,
    Tg_rate_hcw                     = 1,
    prob_hospital_cond_hcw_preAdm   = 0.5,
    ppe_coverage_hcw                = 1,
    ppe_efficacy                    = 1,
    prop_etu                        = 1,
    etu_efficacy                    = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hcw_cond_hcw_comm          = 0,
    prob_hcw_cond_hcw_hospital      = 0.5
  )
  hospital_offspring <- result[result$infection_location == "hospital", ]
  expect_true(nrow(hospital_offspring) > 0,
              info = "expected genPop post-admission hospital offspring (hq = 0)")
  expect_equal(sum(hospital_offspring$class == "HCW"), 0,
               info = "ppe = 1 should remove every HCW hospital recipient post-admission")
})

## --- Pre-admission HCW source: source PPE still applies to all -------
## ppe = 1 on a never-hospitalised HCW parent removes all hospital events
## (source PPE) regardless of recipient class.
test_that("HCW pre-admission only: source PPE blocks all hospital events", {
  parent <- make_parent_info_HCW(hospitalised = FALSE, time_to_outcome = 20)
  set.seed(77)
  result <- offspring_function_hcw(
    parent_info                     = parent,
    mn_contacts_hcw                = 500,
    baseline_risk_hcw                = 1,
    overdisp_contacts_hcw          = 200,
    Tg_shape_hcw                    = 4,
    Tg_rate_hcw                     = 1,
    prob_hospital_cond_hcw_preAdm   = 0.5,
    ppe_coverage_hcw                = 1,
    ppe_efficacy                    = 1,
    prop_etu                        = 1,
    etu_efficacy                    = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hcw_cond_hcw_comm          = 0.2,
    prob_hcw_cond_hcw_hospital      = 0.5
  )
  expect_equal(sum(result$infection_location == "hospital"), 0,
               info = "source PPE = 1 should block all pre-admission hospital events")
  expect_true(sum(result$infection_location == "community") > 0,
              info = "community events should still pass through")
})

## --- Option B multiplicative PPE -------------------------------------
## The per-layer PPE thinning is ppe = ppe_coverage_hcw * ppe_efficacy (= 0.5 * 1
## below). HCW pre-admission with HCW recipient: keep = (1-ppe)^2 (source +
## receiver), while HCW pre-admission with genPop recipient: keep = (1-ppe)
## (source only). With prob_hcw_cond_hcw_hospital = 0.5, the expected HCW fraction
## of the *surviving* pre-admission hospital cohort is
## 0.5*(1-ppe)^2 / [0.5*(1-ppe)^2 + 0.5*(1-ppe)] = (1-ppe) / (1 + (1-ppe)).
## With ppe = 0.5 that's 0.5/1.5 = 1/3.
test_that("HCW pre-admission HCW receiver: multiplicative PPE produces (1-ppe)^2 keep prob", {
  parent <- make_parent_info_HCW(hospitalised = FALSE, time_to_outcome = 30)
  set.seed(123)
  result <- offspring_function_hcw(
    parent_info                     = parent,
    mn_contacts_hcw                = 8000,
    baseline_risk_hcw                = 1,
    overdisp_contacts_hcw          = 5000,
    Tg_shape_hcw                    = 4,
    Tg_rate_hcw                     = 1,
    prob_hospital_cond_hcw_preAdm   = 0.5,
    ppe_coverage_hcw                = 0.5,
    ppe_efficacy                    = 1,
    prop_etu                        = 1,
    etu_efficacy                    = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hcw_cond_hcw_comm          = 0,
    prob_hcw_cond_hcw_hospital      = 0.5
  )
  hospital_offspring <- result[result$infection_location == "hospital", ]
  expect_true(nrow(hospital_offspring) > 200,
              info = "need a reasonable sample to test the ratio")
  hcw_frac <- mean(hospital_offspring$class == "HCW")
  ## Expected ≈ 1/3 under multiplicative PPE. Single-layer PPE would give 0.5.
  expect_true(abs(hcw_frac - 1/3) < 0.05,
              info = paste0("expected HCW fraction ~1/3 under Option B, got ",
                            round(hcw_frac, 3)))
})

## --- Funeral structural reorder: same distribution ---------------------
## With safe_funeral_efficacy = 1 and the parent's funeral safe, no funeral
## offspring should be produced.
test_that("Funeral: safe funeral with efficacy = 1 produces zero offspring", {
  parent <- make_parent_info_HCW(hospitalised = FALSE, time_to_outcome = 20)
  parent$funeral_safety <- "safe"

  set.seed(55)
  result <- offspring_function_funeral(
    parent_info                  = parent,
    mn_contacts_funeral         = 200,
    baseline_risk_funeral         = 1,
    overdisp_contacts_funeral   = 100,
    Tg_shape_funeral             = 10,
    Tg_rate_funeral              = 5,
    safe_funeral_efficacy        = 1,
    prob_hcw_cond_funeral_hcw    = 0.3,
    prob_hcw_cond_funeral_genPop = 0.1
  )
  expect_equal(nrow(result), 0)
})

## With safe_funeral_efficacy = 0 (i.e. a "safe" funeral with no efficacy),
## the funeral should produce the same number of offspring as an unsafe one
## given the same NB draw.
test_that("Funeral: safe funeral with efficacy = 0 behaves like unsafe", {
  parent_unsafe <- make_parent_info_HCW(hospitalised = FALSE, time_to_outcome = 20)
  parent_unsafe$funeral_safety <- "unsafe"
  parent_safe   <- parent_unsafe
  parent_safe$funeral_safety <- "safe"

  set.seed(31)
  res_unsafe <- offspring_function_funeral(
    parent_info                  = parent_unsafe,
    mn_contacts_funeral         = 200,
    baseline_risk_funeral         = 1,
    overdisp_contacts_funeral   = 100,
    Tg_shape_funeral             = 10,
    Tg_rate_funeral              = 5,
    safe_funeral_efficacy        = 0,
    prob_hcw_cond_funeral_hcw    = 0.3,
    prob_hcw_cond_funeral_genPop = 0.1
  )

  set.seed(31)
  res_safe <- offspring_function_funeral(
    parent_info                  = parent_safe,
    mn_contacts_funeral         = 200,
    baseline_risk_funeral         = 1,
    overdisp_contacts_funeral   = 100,
    Tg_shape_funeral             = 10,
    Tg_rate_funeral              = 5,
    safe_funeral_efficacy        = 0,
    prob_hcw_cond_funeral_hcw    = 0.3,
    prob_hcw_cond_funeral_genPop = 0.1
  )
  expect_equal(nrow(res_safe), nrow(res_unsafe))
})

## --- PPE coverage x efficacy decomposition --------------------------------
## The PPE layer thins by coverage * efficacy. Two parameterisations with the
## same product must give bit-identical results under a fixed seed: the keep
## probability collapses to (1 - coverage*efficacy) and the split consumes no
## extra RNG.
test_that("PPE: equal coverage*efficacy products give identical draws (genPop)", {
  parent <- make_parent_info_genPop_hospitalised(time_to_hospitalisation = 0.001)

  set.seed(2024)
  res_a <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = 500,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 200,
    Tg_shape_genPop               = 4,
    Tg_rate_genPop                = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 0,
    general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw              = 1,
    ppe_efficacy                  = 0.5,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.5
  )

  set.seed(2024)
  res_b <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = 500,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 200,
    Tg_shape_genPop               = 4,
    Tg_rate_genPop                = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 0,
    general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw              = 0.5,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.5
  )

  expect_identical(res_a$class, res_b$class)
  expect_identical(res_a$infection_location, res_b$infection_location)
  expect_identical(res_a$time_infection_relative, res_b$time_infection_relative)
})

test_that("PPE: equal coverage*efficacy products give identical draws (HCW)", {
  parent <- make_parent_info_HCW(hospitalised = FALSE, time_to_outcome = 30)

  set.seed(321)
  res_a <- offspring_function_hcw(
    parent_info                     = parent,
    mn_contacts_hcw                = 4000,
    baseline_risk_hcw                = 1,
    overdisp_contacts_hcw          = 2000,
    Tg_shape_hcw                    = 4,
    Tg_rate_hcw                     = 1,
    prob_hospital_cond_hcw_preAdm   = 0.5,
    ppe_coverage_hcw                = 1,
    ppe_efficacy                    = 0.5,
    prop_etu                        = 1,
    etu_efficacy                    = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hcw_cond_hcw_comm          = 0,
    prob_hcw_cond_hcw_hospital      = 0.5
  )

  set.seed(321)
  res_b <- offspring_function_hcw(
    parent_info                     = parent,
    mn_contacts_hcw                = 4000,
    baseline_risk_hcw                = 1,
    overdisp_contacts_hcw          = 2000,
    Tg_shape_hcw                    = 4,
    Tg_rate_hcw                     = 1,
    prob_hospital_cond_hcw_preAdm   = 0.5,
    ppe_coverage_hcw                = 0.5,
    ppe_efficacy                    = 1,
    prop_etu                        = 1,
    etu_efficacy                    = 0,
    general_hospital_quarantine_efficacy = 0,
    prob_hcw_cond_hcw_comm          = 0,
    prob_hcw_cond_hcw_hospital      = 0.5
  )

  expect_identical(res_a$class, res_b$class)
  expect_identical(res_a$infection_location, res_b$infection_location)
  expect_identical(res_a$time_infection_relative, res_b$time_infection_relative)
})

## At full coverage, the efficacy knob alone controls thinning: efficacy = 0
## leaves HCW hospital recipients un-thinned, efficacy = 1 removes them all.
test_that("PPE: efficacy modulates thinning at fixed (full) coverage (genPop)", {
  parent <- make_parent_info_genPop_hospitalised(time_to_hospitalisation = 0.001)

  set.seed(11)
  res_eff0 <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = 500,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 200,
    Tg_shape_genPop               = 4,
    Tg_rate_genPop                = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 0,
    general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw              = 1,
    ppe_efficacy                  = 0,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.5
  )
  hosp0 <- res_eff0[res_eff0$infection_location == "hospital", ]
  expect_true(sum(hosp0$class == "HCW") > 0,
              info = "efficacy = 0 should leave HCW hospital recipients un-thinned")

  set.seed(11)
  res_eff1 <- offspring_function_genPop(
    parent_info                   = parent,
    mn_contacts_genPop           = 500,
    baseline_risk_genPop           = 1,
    overdisp_contacts_genPop     = 200,
    Tg_shape_genPop               = 4,
    Tg_rate_genPop                = 1,
    prop_etu                      = 1,
    etu_efficacy                  = 0,
    general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw              = 1,
    ppe_efficacy                  = 1,
    prob_hcw_cond_genPop_comm     = 0,
    prob_hcw_cond_genPop_hospital = 0.5
  )
  hosp1 <- res_eff1[res_eff1$infection_location == "hospital", ]
  expect_equal(sum(hosp1$class == "HCW"), 0,
               info = "efficacy = 1 at full coverage should remove all HCW hospital recipients")
})

## ppe_efficacy is a required scalar in [0, 1]; an out-of-range value is rejected.
test_that("PPE: ppe_efficacy must be a scalar in [0, 1]", {
  parent <- make_parent_info_genPop_hospitalised(time_to_hospitalisation = 0.001)
  expect_error(
    offspring_function_genPop(
      parent_info                   = parent,
      mn_contacts_genPop           = 10,
      baseline_risk_genPop           = 1,
      overdisp_contacts_genPop     = 5,
      Tg_shape_genPop               = 4,
      Tg_rate_genPop                = 1,
      prop_etu                      = 1,
      etu_efficacy                  = 0,
      general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw              = 1,
      ppe_efficacy                  = 1.5,
      prob_hcw_cond_genPop_comm     = 0,
      prob_hcw_cond_genPop_hospital = 0.5
    ),
    "ppe_efficacy"
  )
})
