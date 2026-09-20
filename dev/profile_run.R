## dev/profile_run.R
## ---------------------------------------------------------------------------
## Dev-only convenience wrapper + profiling harness for branching_process_main().
## NOT part of the package. Reconstructed from the known-good `bpm_args()`
## fixture in tests/testthat/test-time_varying_transmissibility.R, so every
## argument name/value is guaranteed valid (the passing tests use this set).
##
## Usage (from the package root, in R):
##   devtools::load_all()
##   source("dev/profile_run.R")
##   out <- run_sim(seed = 1)        # one simulation, check_final_size = 30000
##   profile_run()                   # profvis the hot loop at scale -> dev/profile.html
## ---------------------------------------------------------------------------

## --- natural-history delay distributions (fixed -> reproducible) -----------
fixed_inc   <- function(n) rep(5,  n)
fixed_hosp  <- function(n) rep(3,  n)
fixed_death <- function(n) rep(7,  n)
fixed_recov <- function(n) rep(10, n)

## run_sim(): a single branching_process_main() call with a full, valid
## parameter set. Override any argument via `...`, e.g.
##   run_sim(mn_contacts_genPop = 3, seed = 7).
##
## Two defaults differ from the bpm_args() test fixture, both profiling-relevant:
##   * check_final_size = 30000  (fixture: 200)   <- the cap you asked for.
##   * population       = 1e6    (fixture: 5000)  <- so the final-size cap is the
##     BINDING constraint. susc <- population - initial_immune is decremented by
##     every offspring (branching_process_main.R:382,650); with population = 5000
##     the loop halts (`susc > 0` fails) at ~5000 cases and never approaches
##     30000. 1e6 keeps susceptibles non-limiting.
run_sim <- function(seed = 1L, check_final_size = 30000L, population = 1e6, ...) {
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
    population                    = population,
    hcw_per_capita                = 0.02,
    check_presymptomatic = FALSE,
    check_final_size              = check_final_size,
    seeding_cases                 = 3,
    seed                          = seed
  )
  do.call(branching_process_main, utils::modifyList(defaults, list(...)))
}

## profile_run(): profvis the hot loop at scale, save an interactive HTML.
##
## The offspring means are bumped to 3 so the epidemic reliably takes off and
## reaches check_final_size. This is PURELY to make the profile representative
## of the loop at N ~ 30000 -- it is not an epidemiological statement. With the
## faithful means (1.5) a single seed can fizzle early, leaving an uninformative
## profile dominated by start-up costs.
##
## NB: at check_final_size = 30000 the per-iteration linear scan that selects the
## next case is O(N^2) overall (branching_process_main.R:~496), so this single
## run may take a few minutes. That cost showing up at the top of the profile is
## itself the headline finding for the scheduler discussion.
profile_run <- function(seed = 1L, out_html = "dev/profile.html") {
  if (!requireNamespace("profvis", quietly = TRUE)) {
    stop("profvis not installed. Run install.packages('profvis') first, ",
         "or use the base-R Rprof fallback at the bottom of this file.")
  }
  p <- profvis::profvis({
    run_sim(seed = seed,
            mn_contacts_genPop  = 3,
            mn_contacts_hcw     = 3,
            mn_contacts_funeral = 3)
  })
  if (!is.null(out_html) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    htmlwidgets::saveWidget(p, normalizePath(out_html, mustWork = FALSE),
                            selfcontained = TRUE)
    message("Saved interactive profile to ", out_html)
  }
  p  # also returned; auto-opens in RStudio's Viewer
}

## --- base-R fallback (no extra packages) -----------------------------------
## If you'd rather not install profvis, this prints the top self-time lines:
##
##   Rprof("dev/Rprof.out", line.profiling = TRUE)
##   invisible(run_sim(seed = 1, mn_contacts_genPop = 3,
##                     mn_contacts_hcw = 3, mn_contacts_funeral = 3))
##   Rprof(NULL)
##   s <- summaryRprof("dev/Rprof.out", lines = "show")
##   head(s$by.self, 20)
