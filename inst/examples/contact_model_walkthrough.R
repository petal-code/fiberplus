## =============================================================================
## contact_model_walkthrough.R
##
## A three-run tour of the contact-first transmission model in fiberplus.
##
##   RUN 1  Defaults, one simulation. What the new outputs look like, and how
##          much presymptomatic transmission this parameter set carries.
##   RUN 2  West Africa Q curves driving the time-varying NPIs, with the contact
##          layer switched OFF -- flat risk tiers, no tracing. This is an exact
##          recapitulation of the pre-contact model.
##   RUN 3  The same scenario with the contact layer switched ON -- five risk
##          tiers, tracing probability rising with tier, and traced cases
##          hospitalised with probability 0.9 one day after symptom onset.
##
## Runs 2 and 3 are made comparable by holding the expected infections per case
## fixed (see `contacts_for_target_offspring()` below), so the only thing that
## differs between them is the contact structure and what tracing does with it.
##
## Requires:
##   * fiberplus (this package)
##   * a local checkout of petal-code/antiviral_pep_hcw_paper, for the scalar
##     parameter defaults and the West Africa Q curve. Point PAPER_REPO at it.
##
## Usage:
##   Rscript inst/examples/contact_model_walkthrough.R
##   # or, from an R session at the package root:
##   #   devtools::load_all(); source("inst/examples/contact_model_walkthrough.R")
## =============================================================================

PAPER_REPO <- Sys.getenv("PAPER_REPO", unset = "../antiviral_pep_hcw_paper")

## --- Load the package -------------------------------------------------------
if (!exists("branching_process_main", mode = "function")) {
  if (requireNamespace("devtools", quietly = TRUE) && dir.exists("R")) {
    devtools::load_all(".", quiet = TRUE)
  } else if (requireNamespace("fiberplus", quietly = TRUE)) {
    library(fiberplus)
  } else {
    for (f in list.files("R", full.names = TRUE, pattern = "[.]R$")) source(f)
  }
}

setup_file <- file.path(PAPER_REPO, "functions", "setup_model_parameters.R")
qcurve_file <- file.path(PAPER_REPO, "data-processed", "WestAfrica_QCurve",
                         "WestAfrica_QCurve_Fit.rds")
if (!file.exists(setup_file) || !file.exists(qcurve_file)) {
  stop("Cannot find the paper repo at '", PAPER_REPO, "'.\n",
       "Clone petal-code/antiviral_pep_hcw_paper and set PAPER_REPO to its path, e.g.\n",
       "  PAPER_REPO=../antiviral_pep_hcw_paper Rscript inst/examples/contact_model_walkthrough.R",
       call. = FALSE)
}
## Provides DEFAULT_SCALAR_INPUTS, make_base_args(), build_time_varying_args().
suppressWarnings(source(setup_file))

rule <- function(title) cat("\n", strrep("=", 78), "\n", title, "\n",
                            strrep("=", 78), "\n", sep = "")


## =============================================================================
## Helpers
## =============================================================================

## --- Translate the paper's parameter set to the contact parameterisation -----
##
## make_base_args() emits the pre-contact transmission arguments. Exactly six of
## them no longer exist:
##
##   mn_offspring_{genPop,hcw,funeral}  ->  mn_contacts_*  + baseline_risk_*
##   overdisp_offspring_*               ->  overdisp_contacts_*
##
## With a FLAT risk structure and baseline_risk = 1 every contact transmits, so
## the contact draw *is* the offspring draw and the translated run is exactly
## equivalent to the old model -- same distribution, same mean, same dispersion.
## That is the identity the whole contact layer is built on, and it is what makes
## RUN 2 a faithful recapitulation rather than an approximation.
as_contact_args <- function(args,
                            baseline_risk = 1,
                            contact_risk  = make_contact_risk()) {
  rr_bar <- contact_risk$mean_relative_risk
  for (route in c("genPop", "hcw", "funeral")) {
    mn_old <- args[[paste0("mn_offspring_", route)]]
    if (is.null(mn_old)) next
    ## Preserve the expected infections per case: mn_contacts * p0 * rr_bar must
    ## equal the old mean offspring count. Handles a time-varying mean too.
    args[[paste0("mn_contacts_", route)]] <- if (is.function(mn_old)) {
      local({ f <- mn_old; function(t) f(t) / (baseline_risk * rr_bar) })
    } else {
      mn_old / (baseline_risk * rr_bar)
    }
    args[[paste0("baseline_risk_", route)]]    <- baseline_risk
    args[[paste0("overdisp_contacts_", route)]] <- args[[paste0("overdisp_offspring_", route)]]
    args[[paste0("mn_offspring_", route)]]      <- NULL
    args[[paste0("overdisp_offspring_", route)]] <- NULL
  }
  args$contact_risk <- contact_risk
  args
}

## Mean contacts needed to deliver a target number of infections per case, given
## a per-contact baseline risk and a risk structure. This is the contact-layer
## identity rearranged:  R0 = mn_contacts * baseline_risk * mean_relative_risk.
contacts_for_target_offspring <- function(target_offspring, baseline_risk, contact_risk) {
  target_offspring / (baseline_risk * contact_risk$mean_relative_risk)
}

## --- West Africa Q curve -> the paper's scenario matrix ----------------------
##
## The fitted curve carries the six latent response parameters on a day grid.
## Renaming them to the published column names (as
## analyses/01_latent_response_parameter_estimation/03_Combine_QCurves.R does)
## gives a one-scenario matrix that the paper's own build_time_varying_args()
## can consume, so the time-varying NPI curves here are built by the paper's
## code, not re-implemented.
west_africa_scenario_matrix <- function(path = qcurve_file) {
  fit <- readRDS(path)
  cs  <- as.data.frame(fit$curve_summ)[, c("parameter", "relative_day", "mean")]
  wide <- stats::reshape(cs, idvar = "relative_day", timevar = "parameter",
                         direction = "wide")
  names(wide) <- sub("^mean[.]", "", names(wide))
  wide <- wide[order(wide$relative_day), ]

  data.frame(
    scenario                 = "west_africa",
    scenario_label           = "West Africa (fitted Q curve)",
    relative_day             = wide$relative_day,
    prob_hosp                = wide$p_hosp,
    delay_hosp               = wide$delay_hosp,
    prob_unsafe_funeral_comm = wide$p_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp = wide$p_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu  = 0,          # always 0 in this analysis
    prop_etu                 = wide$p_ETU,
    ipc_helper               = wide$latent_IPC,   # drives ppe_coverage_hcw(t)
    stringsAsFactors         = FALSE
  )
}

## --- Compact summary of one run ---------------------------------------------
report <- function(out, label) {
  s <- summarise_output(out$tdf, sim_info = out$sim_info,
                        contact_log = out$contact_log)
  cat(sprintf("%-34s %s\n", "scenario", label))
  cat(sprintf("%-34s %d\n",   "cases (total)",        s$n_cases_total))
  cat(sprintf("%-34s %d / %d\n", "  genPop / HCW",    s$n_cases_genPop, s$n_cases_HCW))
  cat(sprintf("%-34s %d\n",   "deaths",               s$n_deaths_total))
  cat(sprintf("%-34s %.0f days\n", "outbreak duration", s$outbreak_duration_days))
  cat(sprintf("%-34s %.3f / %.3f / %.3f\n", "transmission comm/hosp/funeral",
              s$prop_comm, s$prop_hosp, s$prop_funeral))
  cat(sprintf("%-34s %.1f%%\n", "cases hospitalised",
              100 * mean(out$tdf$hospitalisation[out$tdf$offspring_generated], na.rm = TRUE)))
  if (!is.na(s$n_contacts_total)) {
    cat(sprintf("%-34s %d\n",     "contacts generated",  s$n_contacts_total))
    cat(sprintf("%-34s %.2f\n",   "  contacts per case", s$contacts_per_case))
    cat(sprintf("%-34s %d (%.1f%%)\n", "  traced",
                s$n_contacts_traced, 100 * s$prop_contacts_traced))
    cat(sprintf("%-34s %d (%.1f%%)\n", "cases arising from traced contacts",
                s$n_cases_traced, 100 * s$prop_cases_traced))
  }
  invisible(s)
}


## =============================================================================
## RUN 1 -- defaults, one simulation, and what comes out
## =============================================================================
rule("RUN 1  Default parameters, single simulation")

args1 <- as_contact_args(suppressWarnings(make_base_args(overrides = list(
  check_final_size = 3000,
  ## The scenario curves supply these in runs 2 and 3; fix them here so run 1
  ## stands alone.
  prob_hospitalised_genPop = 0.35, prob_hospitalised_hcw = 0.35,
  p_unsafe_funeral_comm_genPop = 0.9, p_unsafe_funeral_comm_hcw = 0.9,
  p_unsafe_funeral_hosp_genPop = 0.1, p_unsafe_funeral_hosp_hcw = 0.1,
  prop_etu = 0.3, ppe_coverage_hcw = 0.3, hospitalisation_delay_factor = 5
))))
args1$seed <- 1

## How much transmission happens before the infector shows symptoms? This is not
## a parameter -- it falls out of the incubation period (mean 8.5 d) sitting
## inside the generation time (mean 15.4 d). It caps what any onset-triggered
## intervention (tracing, isolation, care-seeking) can ever achieve, so it is
## worth knowing before interpreting runs 2 and 3.
ps <- approx_presymptomatic_transmission(args1, n = 50000, seed = 1)
cat(sprintf("presymptomatic share of transmission: genPop %.1f%%, HCW %.1f%%\n",
            100 * ps$genPop, 100 * ps$hcw))
cat(sprintf("  (mean incubation %.1f d, mean time to outcome %.1f d)\n\n",
            ps$mean_incubation, ps$mean_time_to_outcome))

## Unmitigated R0 implied by this parameter set, via the single-type approximation.
r0_1 <- approx_r0(args1, n = 50000, seed = 1)
cat(sprintf("approximate R0: %.3f  (direct %.3f, funeral %.3f; D = %.3f, F = %.3f)\n\n",
            r0_1$R0, r0_1$R0_direct, r0_1$R0_funeral, r0_1$D, r0_1$F))

out1 <- suppressWarnings(do.call(branching_process_main, args1))
report(out1, "defaults")

cat("\n-- transmission tree, new columns ------------------------------------\n")
tree_cols <- c("id", "class", "infection_location", "contact_risk_category",
               "traced", "hospitalisation", "outcome")
print(utils::head(out1$tdf[out1$tdf$offspring_generated, tree_cols], 6))

cat("\n-- contact log -------------------------------------------------------\n")
cat("One row per contact generated, whether or not it became an infection.\n\n")
print(utils::head(out1$contact_log[, c("parent", "case_id", "record_type", "class",
                                       "infection_location", "contact_risk_category",
                                       "transmission_prob", "traced", "blocked_by")], 6))
cat("\nwhy contacts did not become infections:\n")
print(table(out1$contact_log$blocked_by, useNA = "ifany"))


## =============================================================================
## RUN 2 -- West Africa Q curves, contact layer OFF (the old model)
## =============================================================================
rule("RUN 2  West Africa Q curves, contact layer OFF (recapitulates the old model)")

wa_matrix <- west_africa_scenario_matrix()
cat(sprintf("Q curve grid: %d days (%.0f to %.0f)\n",
            nrow(wa_matrix), min(wa_matrix$relative_day), max(wa_matrix$relative_day)))
cat("NPI curves at day 0 -> day 357:\n")
for (nm in c("prob_hosp", "delay_hosp", "prop_etu", "ipc_helper",
             "prob_unsafe_funeral_comm")) {
  cat(sprintf("  %-26s %.3f -> %.3f\n", nm,
              wa_matrix[[nm]][1], wa_matrix[[nm]][nrow(wa_matrix)]))
}

tv <- build_time_varying_args(scenario_id = "west_africa", matrix = wa_matrix)
tv <- tv[setdiff(names(tv), c("scenario_label", "scenario_matrix"))]

## Flat tiers, baseline risk 1: every contact transmits, so this is exactly the
## pre-contact model. No tracing, no admission speed-up.
base2 <- make_base_args(overrides = list(check_final_size = 3000))
args2 <- as_contact_args(c(base2, tv))   # flat structure, baseline_risk = 1
args2$seed <- 2
args2$check_presymptomatic <- FALSE

cat(sprintf("\nrisk structure: %d flat tiers, mean relative risk %.3f, no tracing\n",
            args2$contact_risk$n_levels, args2$contact_risk$mean_relative_risk))
cat(sprintf("mean contacts genPop = %.3f, baseline risk = %.3f  -> %.3f infections/case\n",
            args2$mn_contacts_genPop, args2$baseline_risk_genPop,
            args2$mn_contacts_genPop * args2$baseline_risk_genPop *
              args2$contact_risk$mean_relative_risk))
cat("   (identical to the old mn_offspring_genPop =",
    DEFAULT_SCALAR_INPUTS$mn_offspring_genPop, ")\n\n")

out2 <- suppressWarnings(do.call(branching_process_main, args2))
s2 <- report(out2, "West Africa, contact layer off")


## =============================================================================
## RUN 3 -- same scenario, contact layer ON
## =============================================================================
rule("RUN 3  West Africa Q curves, contact layer ON (5 tiers + tracing)")

## Five tiers spanning a ten-fold risk gradient. The top tier is the reference,
## so relative risks run 0.1 -> 1 and `baseline_risk` is the transmission
## probability of the riskiest contact -- a household / caregiving exposure, the
## thing secondary attack rates actually measure.
##
## Traceability rises with risk: a casual market contact is hard to find, a
## household member is not.
risk5 <- make_contact_risk(
  fractions     = c(0.40, 0.25, 0.20, 0.10, 0.05),
  relative_risk = c(0.10, 0.18, 0.32, 0.56, 1.00),
  trace_prob    = c(0.10, 0.25, 0.45, 0.65, 0.85),
  labels        = c("casual", "community", "repeated", "close", "household")
)
print(risk5)

## Hold expected infections per case fixed at the old value, so RUN 2 vs RUN 3
## isolates the effect of tracing rather than confounding it with a change in
## transmissibility.
p0 <- 0.40   # top-tier (household) per-contact transmission probability
base3 <- make_base_args(overrides = list(check_final_size = 3000))
args3 <- c(base3, tv)
for (route in c("genPop", "hcw", "funeral")) {
  args3[[paste0("mn_contacts_", route)]] <- contacts_for_target_offspring(
    args3[[paste0("mn_offspring_", route)]], p0, risk5)
  args3[[paste0("baseline_risk_", route)]]     <- p0
  args3[[paste0("overdisp_contacts_", route)]] <- args3[[paste0("overdisp_offspring_", route)]]
  args3[[paste0("mn_offspring_", route)]]      <- NULL
  args3[[paste0("overdisp_offspring_", route)]] <- NULL
}
args3$contact_risk <- risk5

## Contact tracing, and what it does to a traced case.
args3$trace_coverage                  <- 1.0   # full programme reach; tier trace_prob does the rest
args3$prob_hospitalised_traced        <- 0.9   # absolute, replaces the Q-curve prob_hosp(t)
args3$onset_to_hospitalisation_traced <- 1.0   # admitted 1 day after onset (caps their own delay)
args3$seed                            <- 2     # same seed as RUN 2
args3$check_presymptomatic            <- FALSE

cat(sprintf("\nmean contacts genPop = %.3f, baseline risk = %.3f, mean rr = %.3f\n",
            args3$mn_contacts_genPop, p0, risk5$mean_relative_risk))
cat(sprintf("  -> %.3f infections/case, matching RUN 2\n",
            args3$mn_contacts_genPop * p0 * risk5$mean_relative_risk))
cat(sprintf("case-weighted trace probability: %.3f",  risk5$mean_trace_prob_cases))
cat(sprintf("  (vs %.3f contact-weighted)\n",
            sum(risk5$fractions * risk5$trace_prob)))
cat("  Cases over-represent high-risk tiers, which are also the traceable ones,\n")
cat("  so tracing reaches more CASES than a contact-weighted average suggests.\n\n")

out3 <- suppressWarnings(do.call(branching_process_main, args3))
s3 <- report(out3, "West Africa, contact layer on")


## =============================================================================
## Comparison
## =============================================================================
rule("RUN 2 vs RUN 3")

cmp <- data.frame(
  metric = c("cases", "deaths", "HCW cases", "% hospitalised",
             "% transmission in community", "% transmission in hospital"),
  contact_layer_off = c(
    s2$n_cases_total, s2$n_deaths_total, s2$n_cases_HCW,
    round(100 * mean(out2$tdf$hospitalisation[out2$tdf$offspring_generated], na.rm = TRUE), 1),
    round(100 * s2$prop_comm, 1), round(100 * s2$prop_hosp, 1)),
  contact_layer_on = c(
    s3$n_cases_total, s3$n_deaths_total, s3$n_cases_HCW,
    round(100 * mean(out3$tdf$hospitalisation[out3$tdf$offspring_generated], na.rm = TRUE), 1),
    round(100 * s3$prop_comm, 1), round(100 * s3$prop_hosp, 1)),
  stringsAsFactors = FALSE
)
print(cmp, row.names = FALSE)

cat("\nHCW cases RISE under tracing, even though total cases fall. That is not a bug:\n")
cat("admitting more cases, sooner, moves transmission out of the community and into\n")
cat("the hospital, where health workers are the ones exposed. Whether that trade is\n")
cat("worth it depends on PPE coverage -- which is exactly the kind of question the\n")
cat("contact layer exists to let you ask.\n")

cat("\nper-tier realised attack rate (RUN 3): share of contacts in each tier that\n")
cat("became infections, after every intervention layer.\n")
print(round(s3$attack_rate_by_risk_tier, 4))

cat("\nAdmission delay from symptom onset (RUN 3):\n")
adm <- out3$tdf[out3$tdf$offspring_generated & out3$tdf$hospitalisation, ]
if (nrow(adm) > 0) {
  d <- adm$time_hospitalisation_relative - adm$incubation_period
  cat(sprintf("  traced   n = %4d, mean %.2f days\n", sum(adm$traced), mean(d[adm$traced])))
  cat(sprintf("  untraced n = %4d, mean %.2f days\n", sum(!adm$traced), mean(d[!adm$traced])))
}

cat("\nNote: single seeds. For anything quantitative, run many replicates and\n")
cat("compare distributions -- these runs differ in their RNG consumption, so a\n")
cat("shared seed does not give a paired comparison.\n")
