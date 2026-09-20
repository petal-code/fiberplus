## =============================================================================
## tracing_tier_sweep.R
##
## How much does contact tracing buy, as you widen it from the highest-risk
## contacts outward?
##
## West Africa scenario (fitted Q curves driving the time-varying NPIs), five
## contact risk tiers, R0 calibrated to 1.35. Contact tracing is then switched on
## progressively:
##
##   none        nobody traced                      (baseline)
##   top 1       only the highest-risk tier traced
##   top 2       highest two tiers traced
##   ...
##   all 5       every contact traced
##
## Traced cases are hospitalised with probability 0.9, one day after symptom
## onset. Everything else is held fixed, so the sweep isolates *who you trace*.
##
## Run twice: with presymptomatic transmission allowed (the model's natural
## behaviour) and removed. R0 is re-solved to 1.35 in each, so the two panels are
## on equal footing at t = 0.
##
## Output: a PNG boxplot (x = tracing setup, y = final size) plus an RDS of the
## raw results for re-plotting.
##
## Requires a local checkout of petal-code/antiviral_pep_hcw_paper; point
## PAPER_REPO at it.
##
## Usage:
##   PAPER_REPO=../antiviral_pep_hcw_paper Rscript inst/examples/tracing_tier_sweep.R
## =============================================================================

PAPER_REPO <- Sys.getenv("PAPER_REPO", unset = "../antiviral_pep_hcw_paper")
OUT_DIR    <- Sys.getenv("OUT_DIR",    unset = ".")
N_REPS     <- as.integer(Sys.getenv("N_REPS",   unset = "20"))
CAP        <- as.integer(Sys.getenv("CAP",      unset = "10000"))
R0_TARGET  <- as.numeric(Sys.getenv("R0_TARGET", unset = "1.35"))
N_CORES    <- max(1L, as.integer(Sys.getenv("N_CORES", unset = "4")))

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

setup_file  <- file.path(PAPER_REPO, "functions", "setup_model_parameters.R")
qcurve_file <- file.path(PAPER_REPO, "data-processed", "WestAfrica_QCurve",
                         "WestAfrica_QCurve_Fit.rds")
if (!file.exists(setup_file) || !file.exists(qcurve_file)) {
  stop("Cannot find the paper repo at '", PAPER_REPO, "'. Set PAPER_REPO to its path.",
       call. = FALSE)
}
suppressWarnings(source(setup_file))


## =============================================================================
## Scenario setup
## =============================================================================

## --- West Africa Q curve -> the paper's time-varying NPI curves --------------
## The fitted curve carries the six latent response parameters on a day grid.
## Renaming them to the published column names gives a one-scenario matrix that
## the paper's own build_time_varying_args() consumes, so the NPI mapping has a
## single source of truth rather than being re-implemented here.
fit  <- readRDS(qcurve_file)
cs   <- as.data.frame(fit$curve_summ)[, c("parameter", "relative_day", "mean")]
wide <- stats::reshape(cs, idvar = "relative_day", timevar = "parameter", direction = "wide")
names(wide) <- sub("^mean[.]", "", names(wide))
wide <- wide[order(wide$relative_day), ]

wa_matrix <- data.frame(
  scenario                 = "west_africa",
  scenario_label           = "West Africa (fitted Q curve)",
  relative_day             = wide$relative_day,
  prob_hosp                = wide$p_hosp,
  delay_hosp               = wide$delay_hosp,
  prob_unsafe_funeral_comm = wide$p_unsafe_funeral_comm,
  prob_unsafe_funeral_hosp = wide$p_unsafe_funeral_hosp,
  prob_unsafe_funeral_etu  = 0,
  prop_etu                 = wide$p_ETU,
  ipc_helper               = wide$latent_IPC,
  stringsAsFactors         = FALSE
)
tv <- build_time_varying_args("west_africa", wa_matrix)
tv <- tv[setdiff(names(tv), c("scenario_label", "scenario_matrix"))]

## --- Contact risk tiers ------------------------------------------------------
## Five tiers with relative risks 0.2 .. 1.0. The lowest-risk tier is five times
## as common as the highest, so tier frequency falls as risk rises -- many casual
## contacts, few household ones. Weights 5:4:3:2:1 give exactly that ratio.
##
## The top tier is the reference (relative risk 1), so `baseline_risk` is the
## per-contact transmission probability of the riskiest contact.
TIER_LABELS <- c("rr_0.2", "rr_0.4", "rr_0.6", "rr_0.8", "rr_1.0")
tier_weights <- c(5, 4, 3, 2, 1)

make_risk <- function(trace_prob) {
  make_contact_risk(
    fractions     = tier_weights / sum(tier_weights),
    relative_risk = c(0.2, 0.4, 0.6, 0.8, 1.0),
    trace_prob    = trace_prob,
    labels        = TIER_LABELS
  )
}
risk_no_tracing <- make_risk(rep(0, 5))

## --- Tracing setups ----------------------------------------------------------
## Widen the net from the top tier outward. Traced with probability 1 if the tier
## is included, 0 if not.
TRACING_SETUPS <- list(
  `none`  = c(0, 0, 0, 0, 0),
  `top 1` = c(0, 0, 0, 0, 1),
  `top 2` = c(0, 0, 0, 1, 1),
  `top 3` = c(0, 0, 1, 1, 1),
  `top 4` = c(0, 1, 1, 1, 1),
  `all 5` = c(1, 1, 1, 1, 1)
)

## --- Base argument list ------------------------------------------------------
## Mean contacts is fixed; the R0 solver then finds the baseline per-contact risk
## that delivers the target R0 at the specified funeral share (see FUNERAL_SHARE
## below).
MN_CONTACTS <- c(genPop = 15, hcw = 15, funeral = 20)

base_args <- function() {
  a <- c(suppressWarnings(make_base_args(overrides = list(check_final_size = CAP))), tv)
  for (r in c("genPop", "hcw", "funeral")) {
    a[[paste0("mn_contacts_", r)]]       <- MN_CONTACTS[[r]]
    a[[paste0("overdisp_contacts_", r)]] <- a[[paste0("overdisp_offspring_", r)]]
  }
  ## The paper's helper still emits the pre-contact-model offspring arguments.
  ## They have no counterpart in branching_process_main() any more, so drop them.
  a <- a[!grepl("^(mn|overdisp)_offspring_", names(a))]
  a$contact_risk         <- risk_no_tracing
  a$check_presymptomatic <- FALSE
  ## The sweep only needs final sizes, so skip the per-contact log (it is the
  ## bulk of the memory and time in a 10,000-case run) and stay quiet about
  ## hitting the cap -- we count that ourselves.
  a$return_contact_log   <- FALSE
  a$quiet                <- TRUE
  a
}

## Share of t = 0 transmission going through unsafe funerals. This is an INPUT,
## the same way the paper treats it -- `run_single_simulation.R` sets
## FUNERAL_FRAC <- 0.25 and solves the offspring means from it, and the ABC fits
## `prop_funeral` as a free parameter. Do not try to back it out of the package's
## default offspring means: those are placeholders the paper overwrites, and the
## share they imply (~0.09, because F carries the CFR and safe-burial thinning)
## is not the calibrated quantity.
FUNERAL_SHARE <- as.numeric(Sys.getenv("FUNERAL_SHARE", unset = "0.25"))
funeral_share <- FUNERAL_SHARE

## --- Calibrate R0 separately for each presymptomatic setting -----------------
## Removing presymptomatic transmission shifts Q_g and hence D, so the same
## baseline risk gives a slightly different R0. Re-solving keeps both panels at
## R0 = R0_TARGET at t = 0, which is what makes them comparable.
calibrate <- function(presympt) {
  a <- base_args()
  a$presymptomatic_transmission <- presympt
  s <- solve_baseline_risk_for_r0(R0_TARGET, a,
                                  proportion_transmission_from_funerals = funeral_share,
                                  n = 50000, seed = 1)
  list(genPop = s$baseline_risk_genPop_required,
       funeral = s$baseline_risk_funeral_required,
       D = s$D_direct_multiplier, F = s$F_funeral_multiplier)
}
cal <- list(`presymptomatic on` = calibrate(TRUE),
            `presymptomatic off` = calibrate(FALSE))

cat("== Scenario ==================================================================\n")
cat(sprintf("West Africa Q curves | R0 target %.2f | cap %d | %d reps per point\n",
            R0_TARGET, CAP, N_REPS))
cat(sprintf("mean contacts: genPop %g, HCW %g, funeral %g\n",
            MN_CONTACTS[["genPop"]], MN_CONTACTS[["hcw"]], MN_CONTACTS[["funeral"]]))
cat(sprintf("funeral share of R0 at t = 0 (specified, as in the paper): %.2f\n", funeral_share))
print(data.frame(
  tier          = TIER_LABELS,
  fraction      = round(risk_no_tracing$fractions, 4),
  relative_risk = risk_no_tracing$relative_risk,
  case_share    = round(risk_no_tracing$case_weights, 4)
), row.names = FALSE)
cat(sprintf("mean relative risk: %.4f\n\n", risk_no_tracing$mean_relative_risk))
for (nm in names(cal)) {
  cat(sprintf("%-19s baseline risk genPop %.4f, funeral %.4f  (D %.3f, F %.3f)\n",
              nm, cal[[nm]]$genPop, cal[[nm]]$funeral, cal[[nm]]$D, cal[[nm]]$F))
}


## =============================================================================
## Sweep
## =============================================================================

run_one <- function(setup_name, presympt, seed) {
  a <- base_args()
  cl <- cal[[if (presympt) "presymptomatic on" else "presymptomatic off"]]
  a$baseline_risk_genPop  <- cl$genPop
  a$baseline_risk_hcw     <- cl$genPop
  a$baseline_risk_funeral <- cl$funeral
  a$contact_risk          <- make_risk(TRACING_SETUPS[[setup_name]])
  a$presymptomatic_transmission <- presympt
  if (setup_name != "none") {
    a$trace_coverage                  <- 1
    a$prob_hospitalised_traced        <- 0.9
    a$onset_to_hospitalisation_traced <- 1
  }
  a$seed <- seed
  o <- suppressWarnings(do.call(branching_process_main, a))
  tdf <- o$tdf[!is.na(o$tdf$time_infection_absolute), ]
  data.frame(setup = setup_name, presymptomatic = presympt, seed = seed,
             final_size = nrow(tdf), deaths = sum(tdf$outcome),
             hcw_cases = sum(tdf$class == "HCW"),
             hit_cap = isTRUE(o$sim_info$hit_final_size_cap),
             stringsAsFactors = FALSE)
}

jobs <- expand.grid(setup = names(TRACING_SETUPS), presympt = c(TRUE, FALSE),
                    seed = seq_len(N_REPS), stringsAsFactors = FALSE)
cat(sprintf("\nRunning %d simulations on %d core(s)...\n", nrow(jobs), N_CORES))
t0 <- Sys.time()
res_list <- parallel::mclapply(seq_len(nrow(jobs)), function(i)
  run_one(jobs$setup[i], jobs$presympt[i], jobs$seed[i]), mc.cores = N_CORES)
bad <- !vapply(res_list, is.data.frame, logical(1))
if (any(bad)) stop("simulation failed: ", conditionMessage(res_list[[which(bad)[1]]]))
res <- do.call(rbind, res_list)
res$setup <- factor(res$setup, levels = names(TRACING_SETUPS))
cat(sprintf("done in %.0f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))


## =============================================================================
## Results
## =============================================================================

cat("== Final size by tracing setup ===============================================\n")
cat("median [IQR], and the share of runs censored by the final-size cap.\n\n")
cat(sprintf("%-8s | %-34s | %-34s\n", "", "presymptomatic ON", "presymptomatic OFF"))
cat(sprintf("%-8s | %10s %18s %4s | %10s %18s %4s\n",
            "setup", "median", "IQR", "cap", "median", "IQR", "cap"))
for (s in levels(res$setup)) {
  line <- sprintf("%-8s |", s)
  for (p in c(TRUE, FALSE)) {
    x <- res$final_size[res$setup == s & res$presymptomatic == p]
    cp <- mean(res$hit_cap[res$setup == s & res$presymptomatic == p])
    line <- paste0(line, sprintf(" %10.0f %8.0f-%-9.0f %3.0f%% |",
                                 median(x), quantile(x, .25), quantile(x, .75), 100 * cp))
  }
  cat(line, "\n")
}

cat("\n== HCW cases by tracing setup (median) =======================================\n")
for (s in levels(res$setup)) {
  cat(sprintf("%-8s  presympt ON %5.0f   presympt OFF %5.0f\n", s,
              median(res$hcw_cases[res$setup == s & res$presymptomatic]),
              median(res$hcw_cases[res$setup == s & !res$presymptomatic])))
}

saveRDS(list(results = res, calibration = cal, funeral_share = funeral_share,
             risk = risk_no_tracing, setups = TRACING_SETUPS),
        file.path(OUT_DIR, "tracing_tier_sweep.rds"))


## =============================================================================
## Plot
## =============================================================================
png_path <- file.path(OUT_DIR, "tracing_tier_sweep.png")
grDevices::png(png_path, width = 1250, height = 620, res = 110)
op <- graphics::par(mfrow = c(1, 2), mar = c(5.5, 5, 4, 1.2), mgp = c(3.2, 0.8, 0))

ylim <- range(res$final_size)
panel <- function(presympt, title, fill) {
  sub <- res[res$presymptomatic == presympt, ]
  vals <- split(sub$final_size, sub$setup)
  graphics::boxplot(vals, log = "y", ylim = ylim, col = fill, border = "grey25",
                    outline = FALSE, las = 1, xaxt = "n",
                    ylab = "final outbreak size (log scale)", xlab = "")
  graphics::axis(1, at = seq_along(vals), labels = names(vals), las = 2, cex.axis = 0.95)
  ## Raw replicates, jittered, so the spread and any cap pile-up are visible.
  for (i in seq_along(vals)) {
    graphics::points(jitter(rep(i, length(vals[[i]])), amount = 0.13), vals[[i]],
                     pch = 16, col = grDevices::adjustcolor("grey15", alpha.f = 0.45), cex = 0.7)
  }
  graphics::abline(h = CAP, lty = 3, col = "firebrick")
  graphics::mtext(sprintf("cap = %d", CAP), side = 4, at = CAP, las = 1,
                  cex = 0.65, col = "firebrick", line = -2.2)
  graphics::title(main = title, adj = 0, cex.main = 1.05)
  graphics::mtext("contacts traced (widening from highest-risk tier)",
                  side = 1, line = 4.2, cex = 0.85)
}
panel(TRUE,  "Presymptomatic transmission ON",  "#9ecae1")
panel(FALSE, "Presymptomatic transmission OFF", "#a1d99b")
graphics::par(op)
grDevices::dev.off()

cat(sprintf("\nwrote %s\n", png_path))
cat(sprintf("wrote %s\n", file.path(OUT_DIR, "tracing_tier_sweep.rds")))
cat("\nNote: R0 is calibrated to", R0_TARGET, "at t = 0 with tracing OFF, so every\n")
cat("setup starts from the same transmissibility and the sweep isolates tracing.\n")
