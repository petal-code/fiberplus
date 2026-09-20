## Tests for the contact-first transmission model: the risk-tier structure, the
## R0 approximation and its inversion, the contact log, and the contact-tracing /
## tracing pathway that lets the risk tiers drive the NPIs, and the
## presymptomatic-transmission diagnostic and switch.

## --- helpers -----------------------------------------------------------

bpm_args <- function(...) {
  args <- list(
    mn_contacts_genPop        = 10,
    overdisp_contacts_genPop  = 0.5,
    baseline_risk_genPop      = 0.1,
    Tg_shape_genPop           = 4, Tg_rate_genPop = 0.5,
    mn_contacts_hcw           = 10,
    overdisp_contacts_hcw     = 0.5,
    Tg_shape_hcw              = 4, Tg_rate_hcw = 0.5,
    mn_contacts_funeral       = 15,
    overdisp_contacts_funeral = 1,
    baseline_risk_funeral     = 0.1,
    Tg_shape_funeral          = 20, Tg_rate_funeral = 10,
    incubation_period           = function(n) rgamma(n, shape = 4, rate = 0.5),
    onset_to_hospitalisation    = function(n) rgamma(n, shape = 2, rate = 0.5),
    onset_to_death              = function(n) rgamma(n, shape = 4, rate = 0.5),
    onset_to_recovery           = function(n) rgamma(n, shape = 6, rate = 0.5),
    hospitalisation_to_death    = function(n) rgamma(n, shape = 3, rate = 0.5),
    hospitalisation_to_recovery = function(n) rgamma(n, shape = 5, rate = 0.5),
    prob_symptomatic         = 0.9,
    prob_hospitalised_hcw    = 0.5,
    prob_hospitalised_genPop = 0.4,
    prob_death_comm = 0.7, prob_death_hosp = 0.5,
    prob_hcw_cond_genPop_comm = 0.02, prob_hcw_cond_genPop_hospital = 0.3,
    prob_hcw_cond_hcw_comm = 0.05,    prob_hcw_cond_hcw_hospital = 0.3,
    prob_hospital_cond_hcw_preAdm = 0.5,
    ppe_coverage_hcw = 0.5, ppe_efficacy = 0.7,
    prop_etu = 0.5, etu_efficacy = 0.9,
    general_hospital_quarantine_efficacy = 0.3,
    p_unsafe_funeral_comm_hcw = 0.5, p_unsafe_funeral_hosp_hcw = 0.2,
    p_unsafe_funeral_comm_genPop = 0.6, p_unsafe_funeral_hosp_genPop = 0.2,
    safe_funeral_efficacy = 0.8,
    prob_hcw_cond_funeral_hcw = 0.1, prob_hcw_cond_funeral_genPop = 0.05,
    population = 1e6, hcw_per_capita = 0.01,
    check_presymptomatic = FALSE,
    ## These runs sit at a deliberately small cap, so nearly all of them end
    ## censored. That is fine here -- none of these tests reads the final size --
    ## but the cap warning would otherwise drown the suite and break every
    ## expect_no_warning() that is really about something else.
    quiet = TRUE,
    check_final_size = 300, seeding_cases = 5, seed = 1
  )
  utils::modifyList(args, list(...))
}

parent_genPop <- function(time_infection_absolute = 0,
                          time_to_hospitalisation = NA_real_,
                          time_to_outcome = 30,
                          incubation_period = 5) {
  data.frame(
    id                            = 1L,
    class                         = "genPop",
    infection_location            = "community",
    parent                        = NA_integer_,
    generation                    = 1L,
    time_infection_relative       = 0,
    time_infection_absolute       = time_infection_absolute,
    incubation_period             = incubation_period,
    symptomatic                   = TRUE,
    time_symptom_onset_relative   = incubation_period,
    time_symptom_onset_absolute   = time_infection_absolute + incubation_period,
    hospitalisation               = !is.na(time_to_hospitalisation),
    time_hospitalisation_relative = time_to_hospitalisation,
    time_hospitalisation_absolute = time_infection_absolute + time_to_hospitalisation,
    outcome                       = TRUE,
    outcome_location              = "community",
    time_outcome_relative         = time_to_outcome,
    time_outcome_absolute         = time_infection_absolute + time_to_outcome,
    funeral_safety                = "unsafe",
    n_offspring                   = NA_integer_,
    offspring_generated           = FALSE,
    stringsAsFactors              = FALSE
  )
}

## --- risk structure ----------------------------------------------------

test_that("make_contact_risk validates and derives its summary quantities", {
  r <- make_contact_risk()
  expect_s3_class(r, "fiber_contact_risk")
  expect_equal(r$n_levels, 5)
  expect_equal(r$mean_relative_risk, 1)
  expect_equal(r$max_relative_risk, 1)
  expect_equal(sum(r$case_weights), 1)
  expect_equal(r$trace_prob, rep(0, 5))

  ## Fractions must sum to 1, lengths must match, risks must be positive.
  expect_error(make_contact_risk(fractions = c(0.5, 0.4)), "sum to 1")
  expect_error(make_contact_risk(fractions = c(0.5, 0.5), relative_risk = c(1, 2, 3)),
               "same length")
  expect_error(make_contact_risk(fractions = c(0.5, 0.5), relative_risk = c(1, -2)),
               "strictly positive")
  expect_error(make_contact_risk(fractions = c(0.5, 0.5), relative_risk = c(1, 2),
                                 trace_prob = c(0.5, 1.5)), "\\[0, 1\\]")
  expect_error(make_contact_risk(reference = 9), "between 1 and 5")
})

test_that("relative risks are normalised so the reference tier is the baseline", {
  ## Writing the reference tier's risk as 4 rather than 1 must not change anything:
  ## the baseline risk parameter always describes the reference tier.
  a <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(1, 2, 4, 8))
  b <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(4, 8, 16, 32))
  expect_equal(a$relative_risk, b$relative_risk)
  expect_equal(a$mean_relative_risk, b$mean_relative_risk)

  ## Naming a different reference rescales relative to that tier.
  c3 <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(1, 2, 4, 8),
                          reference = 3)
  expect_equal(c3$relative_risk[3], 1)
  expect_equal(c3$relative_risk, a$relative_risk / a$relative_risk[3])
})

test_that("case weights over-represent high-risk tiers", {
  r <- make_contact_risk(fractions = rep(0.25, 4), relative_risk = c(1, 2, 4, 8))
  ## Cases are drawn from contacts in proportion to fractions * relative risk.
  expect_equal(r$case_weights,
               r$fractions * r$relative_risk / r$mean_relative_risk)
  expect_true(all(diff(r$case_weights) > 0))
  ## The lowest tier is under-represented among cases relative to contacts.
  expect_lt(r$case_weights[1], r$fractions[1])
  expect_gt(r$case_weights[4], r$fractions[4])
})

test_that("contact_risk_gradient builds a log-spaced gradient anchored on the top tier", {
  g <- contact_risk_gradient(n_levels = 5, ratio = 16, trace_prob_range = c(0.1, 0.9))
  ## The highest tier is the reference, so the risks run 1/ratio up to 1.
  expect_equal(g$reference, 5L)
  expect_equal(g$relative_risk[5], 1)
  expect_equal(g$relative_risk[1], 1 / 16)
  expect_equal(g$max_relative_risk, 1)
  ## Log-spaced: successive ratios are constant.
  expect_equal(diff(log(g$relative_risk)), rep(log(16) / 4, 4))
  expect_equal(g$trace_prob, seq(0.1, 0.9, length.out = 5))
  ## ratio = 1 degenerates to a flat structure.
  expect_equal(contact_risk_gradient(ratio = 1)$relative_risk, rep(1, 5))
})

test_that("the highest-risk tier is the reference by default", {
  ## However the relative risks are written, they are rescaled so the top tier is 1
  ## and the baseline risk describes the riskiest contact.
  r <- make_contact_risk(fractions = c(0.6, 0.25, 0.1, 0.04, 0.01),
                         relative_risk = c(1, 2, 5, 10, 25))
  expect_equal(r$reference, 5L)
  expect_equal(r$relative_risk, c(0.04, 0.08, 0.2, 0.4, 1))
  expect_equal(r$max_relative_risk, 1)
  ## Every relative risk is then a valid attenuation factor, and the mean
  ## attenuates rather than amplifies.
  expect_true(all(r$relative_risk > 0 & r$relative_risk <= 1))
  expect_lte(r$mean_relative_risk, 1)
})

test_that("changing the reference tier is a pure reparameterisation", {
  ## Only the product baseline_risk * relative_risk[l] is ever used, so anchoring on
  ## a different tier and rescaling the baseline risk to match must leave every
  ## per-contact probability -- and hence the whole model -- untouched.
  frac <- c(0.6, 0.25, 0.1, 0.04, 0.01)
  rr   <- c(1, 2, 5, 10, 25)
  hi <- make_contact_risk(frac, rr)                  # default: top tier
  lo <- make_contact_risk(frac, rr, reference = 1)   # bottom tier

  p_lo <- 0.02
  p_hi <- p_lo * max(rr)
  expect_equal(p_hi * hi$relative_risk, p_lo * lo$relative_risk)

  ## R0 and the case-weight distribution are likewise invariant.
  expect_equal(20 * p_hi * hi$mean_relative_risk, 20 * p_lo * lo$mean_relative_risk)
  expect_equal(hi$case_weights, lo$case_weights)
})

## --- R0 approximation and its inversion --------------------------------

test_that("solve_baseline_risk_for_r0 round-trips through approx_r0", {
  args <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 6),
                   contact_risk_funeral = contact_risk_gradient(5, ratio = 3))
  args$baseline_risk_genPop <- NULL
  args$baseline_risk_funeral <- NULL

  set.seed(11)
  fit <- solve_baseline_risk_for_r0(R0 = 2.2, args = args,
                                    proportion_transmission_from_funerals = 0.25,
                                    n = 20000, seed = 3)
  args$baseline_risk_genPop  <- fit$baseline_risk_genPop_required
  args$baseline_risk_funeral <- fit$baseline_risk_funeral_required

  fwd <- approx_r0(args, invariants = fit$invariants)
  expect_equal(fwd$R0, 2.2)
  ## The funeral share must land where it was asked to.
  expect_equal(fwd$R0_funeral / fwd$R0, 0.25)
})

test_that("R0 is linear in the mean contact number and the mean relative risk", {
  args <- bpm_args()
  inv <- compute_r0_invariants(args, n = 20000, seed = 5)

  r1 <- approx_r0(args, invariants = inv)
  r2 <- approx_r0(bpm_args(mn_contacts_genPop = 20), invariants = inv)
  ## Doubling contacts doubles the direct contribution.
  expect_equal(r2$R0_direct, 2 * r1$R0_direct)

  ## Doubling every relative risk is absorbed by the reference-tier normalisation,
  ## so it must NOT change R0; changing the SPREAD does.
  flat <- make_contact_risk(fractions = rep(0.2, 5), relative_risk = rep(2, 5))
  r3 <- approx_r0(bpm_args(contact_risk = flat), invariants = inv)
  expect_equal(r3$R0_direct, r1$R0_direct)

  graded <- contact_risk_gradient(5, ratio = 10)
  r4 <- approx_r0(bpm_args(contact_risk = graded), invariants = inv)
  expect_equal(r4$R0_direct, r1$R0_direct * graded$mean_relative_risk)
})

test_that("the contact overdispersion does not enter R0", {
  args <- bpm_args()
  inv <- compute_r0_invariants(args, n = 10000, seed = 7)
  expect_equal(approx_r0(bpm_args(overdisp_contacts_genPop = 0.1), invariants = inv)$R0,
               approx_r0(bpm_args(overdisp_contacts_genPop = 50),  invariants = inv)$R0)
})

test_that("an infeasible R0 target errors and names the achievable ceiling", {
  ## 5 contacts cannot deliver R0 = 20, whatever the baseline risk.
  args <- bpm_args(mn_contacts_genPop = 5)
  args$baseline_risk_genPop <- NULL
  args$baseline_risk_funeral <- NULL
  expect_error(
    solve_baseline_risk_for_r0(R0 = 20, args = args,
                               proportion_transmission_from_funerals = 0,
                               n = 5000, seed = 2),
    "not achievable"
  )
})

test_that("faster admission for traced cases raises Q_g and lowers D", {
  ## Admitting traced cases sooner pushes more of their generation-time mass past
  ## admission, where hospital quarantine can act on it. So Q_g rises and the direct
  ## multiplier D falls.
  base <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.9),
                   trace_coverage = 1)
  inv_slow <- compute_r0_invariants(base, n = 30000, seed = 9)
  inv_fast <- compute_r0_invariants(
    utils::modifyList(base, list(onset_to_hospitalisation_traced = 0.5)),
    n = 30000, seed = 9)

  expect_gt(inv_fast$Q_g, inv_slow$Q_g)
  expect_lt(r0_direct_multiplier(inv_fast, 0.9, 0.3),
            r0_direct_multiplier(inv_slow, 0.9, 0.3))

  ## With no tracing the fast-admission delay has nothing to act on.
  no_trace <- bpm_args(trace_coverage = 0, onset_to_hospitalisation_traced = 0.5)
  expect_equal(compute_r0_invariants(no_trace, n = 30000, seed = 9)$p_traced, 0)
})

test_that("switching off presymptomatic transmission raises Q_g", {
  ## Truncating contact times to start at symptom onset removes the earliest mass,
  ## so a larger share of what remains falls after admission.
  args <- bpm_args()
  inv_with <- compute_r0_invariants(args, n = 30000, seed = 13)
  inv_without <- compute_r0_invariants(
    utils::modifyList(args, list(presymptomatic_transmission = FALSE)),
    n = 30000, seed = 13)
  expect_gt(inv_without$Q_g, inv_with$Q_g)
})

## --- offspring-function behaviour --------------------------------------

test_that("baseline_risk scales the number of infections, not the contacts", {
  ## With a flat structure, halving the baseline risk should roughly halve infections
  ## while leaving the contact count untouched.
  run <- function(p, seed) {
    set.seed(seed)
    o <- offspring_function_genPop(
      parent_info = parent_genPop(),
      mn_contacts_genPop = 400, overdisp_contacts_genPop = 200,
      baseline_risk_genPop = p, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    )
    c(contacts = nrow(attr(o, "contact_log")), infections = nrow(o))
  }
  a <- run(0.8, 100); b <- run(0.4, 100)
  expect_equal(a[["contacts"]], b[["contacts"]])   # same seed, same contact draw
  expect_gt(a[["infections"]], b[["infections"]])
  expect_lt(abs(b[["infections"]] / a[["infections"]] - 0.5), 0.15)
})

test_that("baseline_risk = 1 with flat tiers makes every contact an infection", {
  set.seed(21)
  o <- offspring_function_genPop(
    parent_info = parent_genPop(),
    mn_contacts_genPop = 200, overdisp_contacts_genPop = 100,
    baseline_risk_genPop = 1, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
    prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw = 0, ppe_efficacy = 0,
    prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
  )
  log <- attr(o, "contact_log")
  expect_equal(nrow(o), nrow(log))
  expect_true(all(log$record_type == "infection"))
})

test_that("high-risk tiers transmit more often than low-risk tiers", {
  ## Two tiers, the top five times riskier. The top tier is the reference, so the
  ## baseline risk IS the top tier's transmission probability and the bottom tier
  ## gets a fifth of it.
  risk <- make_contact_risk(fractions = rep(0.5, 2), relative_risk = c(1, 5))
  expect_equal(risk$relative_risk, c(0.2, 1))
  set.seed(31)
  o <- offspring_function_genPop(
    parent_info = parent_genPop(),
    mn_contacts_genPop = 4000, overdisp_contacts_genPop = 2000,
    baseline_risk_genPop = 0.5, contact_risk_genPop = risk,
    Tg_shape_genPop = 4, Tg_rate_genPop = 1,
    prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
    ppe_coverage_hcw = 0, ppe_efficacy = 0,
    prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
  )
  log <- attr(o, "contact_log")
  rate <- tapply(log$record_type == "infection", log$contact_risk_level, mean)
  expect_lt(abs(rate[["1"]] - 0.1), 0.03)
  expect_lt(abs(rate[["2"]] - 0.5), 0.05)
})

test_that("the top tier's probability cannot exceed 1", {
  run <- function(risk, p0) {
    set.seed(71)
    offspring_function_genPop(
      parent_info = parent_genPop(),
      mn_contacts_genPop = 50, overdisp_contacts_genPop = 25,
      baseline_risk_genPop = p0, contact_risk_genPop = risk,
      Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    )
  }

  ## Default convention: the top tier IS the reference, so its probability equals
  ## the baseline risk. The constraint can never bind, even at a baseline risk of 1.
  default_risk <- make_contact_risk(fractions = rep(0.5, 2), relative_risk = c(1, 5))
  expect_s3_class(run(default_risk, 1), "data.frame")

  ## Anchoring on the LOWEST tier brings the moving bound back -- and it is enforced.
  low_ref <- make_contact_risk(fractions = rep(0.5, 2), relative_risk = c(1, 5),
                               reference = 1)
  expect_equal(low_ref$max_relative_risk, 5)
  expect_error(run(low_ref, 0.5), "exceeds 1")
  ## Below the bound it is fine: 0.2 * 5 = 1 exactly.
  expect_s3_class(run(low_ref, 0.2), "data.frame")
})

test_that("presymptomatic_transmission = FALSE pushes all contacts past symptom onset", {
  run <- function(presympt, seed = 41) {
    set.seed(seed)
    o <- offspring_function_genPop(
      parent_info = parent_genPop(incubation_period = 6, time_to_outcome = 30),
      mn_contacts_genPop = 3000, overdisp_contacts_genPop = 1500,
      baseline_risk_genPop = 1, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      presymptomatic_transmission = presympt,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    )
    attr(o, "contact_log")$time_contact_relative
  }
  with_pre <- run(TRUE)
  without  <- run(FALSE)

  ## With it on, a substantial share of contacts land before onset at day 6.
  expect_gt(mean(with_pre < 6), 0.3)
  ## With it off, none do -- and the realised generation time is longer.
  expect_equal(sum(without < 6), 0)
  expect_gt(min(without), 6 - 1e-8)
  expect_gt(mean(without), mean(with_pre))

  ## The contact COUNT is unchanged: only the timing is conditioned, not the number.
  expect_equal(length(with_pre), length(without))
})

test_that("approx_presymptomatic_transmission recovers a known share", {
  ## Incubation fixed at 5 days and outcome far away, so the presymptomatic share is
  ## almost exactly the generation-time CDF at day 5.
  args <- bpm_args(
    incubation_period        = function(n) rep(5, n),
    onset_to_death           = function(n) rep(300, n),
    onset_to_recovery        = function(n) rep(300, n),
    prob_hospitalised_genPop = 0,
    prob_hospitalised_hcw    = 0,
    prob_symptomatic         = 1,
    Tg_shape_genPop = 4, Tg_rate_genPop = 0.5
  )
  ps <- approx_presymptomatic_transmission(args, n = 20000, seed = 2)
  expected <- pgamma(5, shape = 4, rate = 0.5) / pgamma(305, shape = 4, rate = 0.5)
  expect_lt(abs(ps$genPop - expected), 0.01)
  expect_equal(ps$mean_incubation, 5)

  ## A longer incubation period leaves less transmission before onset.
  later <- approx_presymptomatic_transmission(
    utils::modifyList(args, list(incubation_period = function(n) rep(12, n))),
    n = 20000, seed = 2)
  expect_gt(later$genPop, ps$genPop)
})

test_that("tracing probability follows the tier and scales with coverage", {
  risk <- make_contact_risk(fractions = rep(0.5, 2), relative_risk = c(1, 1),
                            trace_prob = c(0.2, 0.8))
  traced_rates <- function(cov, seed = 61) {
    set.seed(seed)
    o <- offspring_function_genPop(
      parent_info = parent_genPop(),
      mn_contacts_genPop = 4000, overdisp_contacts_genPop = 2000,
      baseline_risk_genPop = 0.1, contact_risk_genPop = risk,
      trace_coverage = cov, Tg_shape_genPop = 4, Tg_rate_genPop = 1,
      prop_etu = 1, etu_efficacy = 0, general_hospital_quarantine_efficacy = 0,
      ppe_coverage_hcw = 0, ppe_efficacy = 0,
      prob_hcw_cond_genPop_comm = 0, prob_hcw_cond_genPop_hospital = 0
    )
    log <- attr(o, "contact_log")
    tapply(log$traced, log$contact_risk_level, mean)
  }
  full <- traced_rates(1)
  expect_lt(abs(full[["1"]] - 0.2), 0.03)
  expect_lt(abs(full[["2"]] - 0.8), 0.03)

  half <- traced_rates(0.5)
  expect_lt(abs(half[["1"]] - 0.1), 0.03)
  expect_lt(abs(half[["2"]] - 0.4), 0.03)

  ## Tracing is drawn for every contact, not only those that transmit.
  expect_equal(length(full), 2L)
})

## --- contact log integrity ---------------------------------------------

test_that("the contact log accounts for every contact exactly once", {
  out <- do.call(branching_process_main,
                 bpm_args(contact_risk = contact_risk_gradient(5, ratio = 5,
                                                               trace_prob_range = c(0.1, 0.9)),
                          trace_coverage = 0.7))
  log <- out$contact_log
  expect_gt(nrow(log), 0)

  ## Every row is either an infection (with a case id and no block reason) or a
  ## non-infection (with a reason and no case id). No third state.
  inf <- log$record_type == "infection"
  expect_true(all(log$record_type %in% c("contact", "infection")))
  expect_false(anyNA(log$case_id[inf]))
  expect_true(all(is.na(log$blocked_by[inf])))
  expect_true(all(is.na(log$case_id[!inf])))
  expect_false(anyNA(log$blocked_by[!inf]))
  expect_true(all(log$blocked_by[!inf] %in%
                    c("no_transmission", "ppe_quarantine",
                      "safe_funeral", "obv_pep")))

  ## Case ids are unique and every one resolves to a row in the tree.
  expect_false(anyDuplicated(log$case_id[inf]) > 0)
  expect_true(all(log$case_id[inf] %in% out$tdf$id))

  ## The tier, class and traced status recorded on the tree must agree with the log.
  m <- match(log$case_id[inf], out$tdf$id)
  expect_identical(out$tdf$contact_risk_level[m], log$contact_risk_level[inf])
  expect_identical(out$tdf$contact_risk_category[m], log$contact_risk_category[inf])
  expect_identical(out$tdf$traced[m], log$traced[inf])
  expect_identical(out$tdf$class[m], log$class[inf])

  ## Contacts always outnumber infections when the baseline risk is below 1.
  expect_gt(nrow(log), sum(inf))
})

test_that("seed cases carry no risk tier and are never traced", {
  out <- do.call(branching_process_main, bpm_args(trace_coverage = 1))
  seeds <- out$tdf[out$tdf$generation == 1 & !is.na(out$tdf$time_infection_absolute), ]
  expect_true(all(is.na(seeds$contact_risk_level)))
  expect_true(all(!seeds$traced))
})

## --- tracing drives the NPIs -------------------------------------------

test_that("tracing plus fast admission reduces onward transmission", {
  ## Measured as realised offspring per expanded case rather than final size: these
  ## parameters are supercritical, so both arms run into `check_final_size` and the
  ## final size says more about the cap than about transmission.
  traced_risk <- contact_risk_gradient(5, ratio = 5, trace_prob_range = c(0.5, 0.95))
  mean_offspring <- function(extra, seeds = 1:8) {
    vapply(seeds, function(s) {
      o <- do.call(branching_process_main,
                   utils::modifyList(bpm_args(contact_risk = traced_risk, seed = s), extra))
      done <- o$tdf[!is.na(o$tdf$time_infection_absolute) & o$tdf$offspring_generated, ]
      mean(done$n_offspring)
    }, numeric(1))
  }
  off <- mean_offspring(list(trace_coverage = 0))
  on  <- mean_offspring(list(trace_coverage = 1, onset_to_hospitalisation_traced = 0.5,
                             etu_efficacy = 1, general_hospital_quarantine_efficacy = 1))
  expect_lt(mean(on), mean(off))
})

test_that("the traced admission delay caps rather than replaces a case's own delay", {
  args <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.9),
                   trace_coverage = 1, onset_to_hospitalisation_traced = 1,
                   check_final_size = 600)
  out <- do.call(branching_process_main, args)
  real <- out$tdf[!is.na(out$tdf$time_infection_absolute) & out$tdf$hospitalisation, ]
  delay <- real$time_hospitalisation_relative - real$incubation_period

  ## No traced case waits longer than the cap ...
  expect_lte(max(delay[real$traced]), 1 + 1e-8)
  ## ... and some are admitted sooner than it, because the cap never slows anyone down.
  expect_true(any(delay[real$traced] < 1 - 1e-8))
})

test_that("a traced delay above the untraced distribution warns and does nothing", {
  ## onset_to_hospitalisation here has mean 4 days, so a 20-day traced delay can never
  ## bind. That is almost certainly a mistake, so it must warn rather than silently
  ## produce a no-op tracing scenario.
  expect_warning(
    do.call(branching_process_main,
            bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.9),
                     trace_coverage = 1, onset_to_hospitalisation_traced = 20)),
    "not clearly below"
  )
})

test_that("tracing can accelerate admission for traced cases", {
  args <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.9),
                   trace_coverage = 1, onset_to_hospitalisation_traced = 0.5,
                   check_final_size = 600)
  out <- do.call(branching_process_main, args)
  real <- out$tdf[!is.na(out$tdf$time_infection_absolute) & out$tdf$hospitalisation, ]
  ## Admission delay is measured from symptom onset.
  delay <- real$time_hospitalisation_relative - real$incubation_period
  expect_lt(mean(delay[real$traced]), mean(delay[!real$traced]))
})

test_that("tracing can raise the probability of admission for traced cases", {
  args <- bpm_args(contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.5),
                   trace_coverage = 1, prob_hospitalised_multiplier_traced = 2,
                   check_final_size = 800)
  out <- do.call(branching_process_main, args)
  real <- out$tdf[!is.na(out$tdf$time_infection_absolute) & out$tdf$symptomatic, ]
  expect_gt(mean(real$hospitalisation[real$traced]),
            mean(real$hospitalisation[!real$traced]))
})

## --- main-function wiring ----------------------------------------------

test_that("branching_process_main solves baseline risks from r0_target", {
  out <- do.call(branching_process_main,
                 utils::modifyList(
                   bpm_args(r0_target = 1.6, r0_prop_funeral = 0.3,
                            r0_solve_n = 10000, r0_solve_seed = 4),
                   list(baseline_risk_genPop = NULL, baseline_risk_funeral = NULL)))
  expect_equal(out$sim_info$r0_target, 1.6)
  expect_true(is.finite(out$sim_info$baseline_risk_genPop))
  expect_true(is.finite(out$sim_info$baseline_risk_funeral))
  ## HCW parents inherit the genPop per-contact risk when not given their own.
  expect_equal(out$sim_info$baseline_risk_hcw, out$sim_info$baseline_risk_genPop)
})

test_that("supplying both r0_target and a baseline risk is an error", {
  expect_error(
    do.call(branching_process_main, bpm_args(r0_target = 2)),
    "not both"
  )
})

test_that("per-route risk structures are independent and reported back", {
  rg <- contact_risk_gradient(5, ratio = 3)
  rf <- contact_risk_gradient(4, ratio = 9)
  out <- do.call(branching_process_main,
                 bpm_args(contact_risk = rg, contact_risk_funeral = rf))
  ## genPop and HCW inherit the shared structure; funerals use their own.
  expect_equal(out$sim_info$contact_risk_genPop$relative_risk, rg$relative_risk)
  expect_equal(out$sim_info$contact_risk_hcw$relative_risk, rg$relative_risk)
  expect_equal(out$sim_info$contact_risk_funeral$relative_risk, rf$relative_risk)
  expect_equal(out$sim_info$contact_risk_funeral$n_levels, 4)

  ## Funeral contacts must be labelled from the funeral structure.
  log <- out$contact_log
  fun_cats <- unique(log$contact_risk_category[log$infection_location == "funeral"])
  expect_true(all(fun_cats %in% rf$labels))
})

test_that("an infeasible baseline risk fails before the simulation starts", {
  ## Only reachable by anchoring on a tier that is not the riskiest -- under the
  ## default convention the baseline risk is bounded by 1 and nothing else.
  low_ref <- make_contact_risk(fractions = rep(0.2, 5),
                               relative_risk = c(1, 2, 4, 8, 16),
                               reference = 1)
  expect_error(
    do.call(branching_process_main,
            bpm_args(contact_risk = low_ref, baseline_risk_genPop = 0.5)),
    "highest-risk tier"
  )
  ## The same structure under the default convention runs without complaint.
  expect_s3_class(
    do.call(branching_process_main,
            bpm_args(contact_risk = contact_risk_gradient(5, ratio = 16),
                     baseline_risk_genPop = 0.5))$tdf,
    "data.frame"
  )
})

test_that("summarise_output reports contact and tracing counts", {
  out <- do.call(branching_process_main,
                 bpm_args(contact_risk = contact_risk_gradient(5, ratio = 5,
                                                               trace_prob_range = c(0.2, 0.9)),
                          trace_coverage = 0.8, onset_to_hospitalisation_traced = 1))
  s <- summarise_output(out$tdf, sim_info = out$sim_info, contact_log = out$contact_log)

  expect_equal(s$n_contacts_total, nrow(out$contact_log))
  expect_equal(s$n_contacts_infected, sum(out$contact_log$record_type == "infection"))
  expect_gt(s$contacts_per_case, 1)
  expect_gt(s$n_cases_traced, 0)
  expect_lte(s$n_cases_traced, s$n_cases_total)

  ## The per-tier attack rate must increase with the tier's relative risk.
  ar <- s$attack_rate_by_risk_tier
  expect_equal(length(ar), 5L)
  expect_gt(ar[[5]], ar[[1]])

  ## ... and the tiers must come back in TIER order, not alphabetical order.
  ## table() on a character column sorts by name, which silently scrambles the
  ## breakdown as soon as the tiers carry real labels, making a monotone attack
  ## rate look non-monotone.
  named <- make_contact_risk(
    fractions     = rep(0.2, 5),
    relative_risk = c(0.1, 0.2, 0.4, 0.7, 1.0),
    trace_prob    = rep(0.5, 5),
    ## Deliberately chosen so alphabetical order differs from risk order.
    labels        = c("casual", "repeated", "community", "close", "household"))
  out2 <- do.call(branching_process_main,
                  bpm_args(contact_risk = named, trace_coverage = 0.8))
  s3 <- summarise_output(out2$tdf, sim_info = out2$sim_info,
                         contact_log = out2$contact_log)
  expect_identical(names(s3$contacts_by_risk_tier), named$labels)
  expect_identical(names(s3$attack_rate_by_risk_tier), named$labels)
  expect_identical(names(s3$cases_by_risk_tier), named$labels)
  ## Monotone in risk once the order is right.
  expect_gt(s3$attack_rate_by_risk_tier[["household"]],
            s3$attack_rate_by_risk_tier[["casual"]])

  ## Without the log the contact fields are NA but the call still works.
  s2 <- summarise_output(out$tdf, sim_info = out$sim_info)
  expect_true(is.na(s2$n_contacts_total))
  expect_equal(s2$n_cases_traced, s$n_cases_traced)
})

## --- stop reason and the contact-log switch ----------------------------

test_that("a censored run is reported as such, and a finished one is not", {
  ## A run that stops at the cap must be distinguishable from one that ended on
  ## its own. Reading nrow(tdf) alone cannot tell them apart, which is the whole
  ## reason these fields exist.
  capped <- do.call(branching_process_main, bpm_args())
  expect_identical(capped$sim_info$stop_reason, "final_size_cap")
  expect_true(capped$sim_info$hit_final_size_cap)
  expect_gt(capped$sim_info$n_unexpanded, 0)

  ## Subcritical: baseline risk near zero, so the chains die out well short of
  ## the cap and every case gets expanded.
  ended <- do.call(branching_process_main,
                   bpm_args(baseline_risk_genPop = 0.001,
                            baseline_risk_funeral = 0.001, seeding_cases = 2))
  expect_identical(ended$sim_info$stop_reason, "outbreak_ended")
  expect_false(ended$sim_info$hit_final_size_cap)
  expect_identical(ended$sim_info$n_unexpanded, 0L)
})

test_that("the censoring warning fires unless quiet, and quiet keeps the fields", {
  expect_warning(do.call(branching_process_main, bpm_args(quiet = FALSE)),
                 "CENSORED")
  expect_no_warning(do.call(branching_process_main, bpm_args(quiet = TRUE)))
  ## quiet silences the announcement, it does not stop the bookkeeping.
  q <- do.call(branching_process_main, bpm_args(quiet = TRUE))
  expect_true(q$sim_info$hit_final_size_cap)
})

test_that("return_contact_log = FALSE drops the log without touching the tree", {
  ## The log is built from draws that have already happened, so switching it off
  ## must not consume any randomness. If it did, every calibration run would
  ## silently diverge from the equivalent diagnostic run.
  with_log    <- do.call(branching_process_main, bpm_args(return_contact_log = TRUE))
  without_log <- do.call(branching_process_main, bpm_args(return_contact_log = FALSE))

  expect_identical(with_log$tdf, without_log$tdf)
  expect_gt(nrow(with_log$contact_log), 0)
  expect_equal(nrow(without_log$contact_log), 0)
  ## Still a well-formed frame, so downstream code does not need a special case.
  expect_identical(names(without_log$contact_log), names(with_log$contact_log))
})

## --- presymptomatic transmission wiring --------------------------------

test_that("branching_process_main warns when presymptomatic transmission is substantial", {
  ## The default test parameters have a mean incubation of 8 days and a mean generation
  ## time of 8 days, so around half of transmission precedes symptoms. That is exactly
  ## the case the check exists to surface.
  expect_warning(
    do.call(branching_process_main, bpm_args(check_presymptomatic = TRUE)),
    "before the infector develops symptoms"
  )
  ## Silent when the share is below the threshold ...
  expect_no_warning(
    do.call(branching_process_main,
            bpm_args(check_presymptomatic = TRUE, presymptomatic_warn_threshold = 0.99))
  )
  ## ... when the check is switched off ...
  expect_no_warning(do.call(branching_process_main, bpm_args(check_presymptomatic = FALSE)))
  ## ... and when presymptomatic transmission has been removed outright.
  expect_no_warning(
    do.call(branching_process_main,
            bpm_args(check_presymptomatic = TRUE, presymptomatic_transmission = FALSE))
  )
})

test_that("the presymptomatic check does not perturb the simulated trajectory", {
  ## The check draws from the user's delay distributions, so it must save and restore
  ## the random seed. Otherwise every downstream draw would shift and a scenario would
  ## silently change depending on whether diagnostics were switched on.
  quiet <- do.call(branching_process_main, bpm_args(check_presymptomatic = FALSE))
  loud  <- suppressWarnings(
    do.call(branching_process_main, bpm_args(check_presymptomatic = TRUE)))
  expect_identical(quiet$tdf, loud$tdf)
})

test_that("removing presymptomatic transmission reduces final size", {
  ## Transmission that would have happened before onset is not redistributed -- it is
  ## pushed later, into the window where admission and quarantine can act on it. With
  ## quarantine switched on, that means fewer secondary cases.
  final <- function(presympt, seeds = 1:8) {
    vapply(seeds, function(s) {
      o <- do.call(branching_process_main,
                   bpm_args(seed = s, presymptomatic_transmission = presympt,
                            prob_hospitalised_genPop = 0.9, prob_hospitalised_hcw = 0.9,
                            etu_efficacy = 1, general_hospital_quarantine_efficacy = 1,
                            prop_etu = 1))
      done <- o$tdf[!is.na(o$tdf$time_infection_absolute) & o$tdf$offspring_generated, ]
      mean(done$n_offspring)
    }, numeric(1))
  }
  expect_lt(mean(final(FALSE)), mean(final(TRUE)))
})

test_that("presymptomatic_transmission must be a single logical", {
  expect_error(do.call(branching_process_main, bpm_args(presymptomatic_transmission = "no")),
               "single logical")
})

test_that("prob_hospitalised_traced sets an absolute admission probability", {
  ## The multiplier cannot pin the traced probability at a fixed value when the
  ## untraced one is time-varying, so an absolute override exists alongside it.
  args <- bpm_args(
    contact_risk = contact_risk_gradient(5, ratio = 4, trace_prob_range = 0.9),
    trace_coverage = 1,
    prob_hospitalised_genPop = 0.1, prob_hospitalised_hcw = 0.1,
    prob_hospitalised_traced = 0.95,
    check_final_size = 800
  )
  out <- do.call(branching_process_main, args)
  real <- out$tdf[!is.na(out$tdf$time_infection_absolute) & out$tdf$symptomatic, ]
  ## Traced cases are admitted far more often than untraced ones, and close to the
  ## absolute value asked for (below it, since admission must also beat the outcome).
  expect_gt(mean(real$hospitalisation[real$traced]),
            mean(real$hospitalisation[!real$traced]) + 0.3)
  expect_gt(mean(real$hospitalisation[real$traced]), 0.7)
})

test_that("the absolute and multiplier forms of traced hospitalisation are exclusive", {
  expect_error(
    do.call(branching_process_main,
            bpm_args(prob_hospitalised_traced = 0.9,
                     prob_hospitalised_multiplier_traced = 2)),
    "not both"
  )
})
