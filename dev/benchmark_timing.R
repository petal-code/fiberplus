## dev/benchmark_timing.R
## ---------------------------------------------------------------------------
## Wall-clock benchmark of branching_process_main() at the 30,000-case cap with
## ALL control measures switched on. Run the SAME script on the new model (this
## branch) and the old model; same seeds => same outbreak trajectories => same
## work => a fair A/B timing comparison.
##
## Usage (from the package root, in R):
##   devtools::load_all()
##   source("dev/benchmark_timing.R")
##
## To time the OLD model: check it out (cleanest reference is commit e384a31,
## the commit immediately before the optimisation), copy this file in if needed,
## then run the same two lines. Compare the reported "mean elapsed".
##
## NOTE ON RUNTIME: at the 30k cap the per-iteration parent-selection scan is
## still O(N^2) (left unchanged on purpose -- changing it would alter results),
## so each rep takes on the order of a minute or more on the OLD model. 5 reps +
## a warm-up can therefore take ~5-20 min depending on the machine and version.
## ---------------------------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/profile_run.R")   # provides run_sim(): a full, valid parameter set
                              # (identical on the old and new model -- it lives in
                              # dev/ and was not touched by the optimisation).

## --- knobs ------------------------------------------------------------------
## All overridable from the shell without editing this file, e.g.
##   BENCH_OFFSPRING_MEAN=2 bash dev/benchmark_compare.sh
OFFSPRING_MEAN <- as.numeric(Sys.getenv("BENCH_OFFSPRING_MEAN", "3"))
                        # baseline NB mean for all three routes. The community
                        # genPop route is NOT thinned by any control measure, so
                        # this is what drives takeoff. 3 => effective R ~2-2.5 =>
                        # the cap is hit reliably. Dial toward 1.5 for a more
                        # realistic R0, but expect occasional fizzles (watch the
                        # "hit cap?" column) and more iterations (=> slower).
CAP   <- as.integer(Sys.getenv("BENCH_CAP",  "30000"))
REPS  <- as.integer(Sys.getenv("BENCH_REPS", "5"))
SEEDS <- seq_len(REPS)   # 1..REPS; reproducible and identical across versions.

## All control measures ON / at full strength. (PPE/ETU/OBV target the hospital
## & HCW routes; safe burial targets funerals; none thin community genPop.)
control_overrides <- list(
  ppe_coverage_hcw                     = 1.0,   # PPE/IPC coverage
  ppe_efficacy                         = 1.0,
  prop_etu                             = 1.0,   # all hospitalised cases in an ETU
  etu_efficacy                         = 1.0,   # full post-admission quarantine
  general_hospital_quarantine_efficacy = 1.0,
  safe_funeral_efficacy                = 1.0,   # safe burial fully effective
  p_unsafe_funeral_comm_hcw            = 0.1,   # high safe-burial coverage
  p_unsafe_funeral_hosp_hcw            = 0.1,   #   (= 1 - p_unsafe)
  p_unsafe_funeral_comm_genPop         = 0.1,
  p_unsafe_funeral_hosp_genPop         = 0.1,
  obv_pep_enabled                      = TRUE,  # OBV PEP on
  obv_pep_coverage                     = 1.0,
  obv_pep_adherence                    = 1.0,
  seeding_cases                        = 10     # several seeds => negligible
                                                # chance of an early fizzle
)

run_bench <- function(seed, check_final_size = CAP) {
  do.call(run_sim, c(
    list(seed = seed,
         check_presymptomatic = FALSE,
         check_final_size     = check_final_size,
         mn_contacts_genPop  = OFFSPRING_MEAN,
         mn_contacts_hcw     = OFFSPRING_MEAN,
         mn_contacts_funeral = OFFSPRING_MEAN),
    control_overrides
  ))
}

## --- warm-up (untimed): trigger JIT/bytecode compilation of the hot path ----
cat("warming up (untimed, small cap) ...\n")
invisible(run_bench(seed = 999L, check_final_size = 2000L))

## --- timed reps -------------------------------------------------------------
times <- numeric(REPS)
sizes <- integer(REPS)
cat(sprintf("\nbranching_process_main()  |  cap = %d  |  offspring mean = %g  |  reps = %d\n",
            CAP, OFFSPRING_MEAN, REPS))
cat("rep | seed | elapsed (s) | realised cases | hit cap?\n")
cat("----+------+-------------+----------------+----------\n")
for (i in seq_len(REPS)) {
  s <- SEEDS[i]
  el <- system.time(res <- run_bench(seed = s))[["elapsed"]]
  ## realised cases = rows that are actual cases (non-NA infection time); the
  ## unused pre-allocated tail is NA. >= cap means the size cap is what stopped it.
  n_real <- sum(!is.na(res$tdf$time_infection_absolute))
  times[i] <- el
  sizes[i] <- n_real
  cat(sprintf("%3d | %4d | %11.2f | %14d | %s\n",
              i, s, el, n_real, if (n_real >= CAP) "yes" else "NO (fizzled)"))
}
cat("----+------+-------------+----------------+----------\n")
cat(sprintf("\nmean elapsed : %.2f s   (sd %.2f, min %.2f, max %.2f)  over %d reps\n",
            mean(times), stats::sd(times), min(times), max(times), REPS))
cat(sprintf("mean realised cases : %.0f\n", mean(sizes)))
if (any(sizes < CAP)) {
  cat("\n!! Some reps did not reach the cap -- their timings are not comparable.\n",
      "   Increase OFFSPRING_MEAN at the top of this script and re-run.\n", sep = "")
}
