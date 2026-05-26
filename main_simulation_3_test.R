# Light/fast test of FSEM simulation + fitting with two DEPENDENT functional
# latent variables (eta1, eta2). eta2 is generated as before; eta1 is generated
# CONDITIONAL on eta2 via one of three functional regression scenarios:
#
#   Scenario "fixed"      : eta1     = alpha1 * eta2 + e
#                           (alpha1 is a scalar; supplied by parameter_set1 = 0.7)
#
#   Scenario "concurrent" : eta1(t)  = alpha2(t) * eta2(t) + e(t)
#                           (alpha2(t) supplied by parameter_set1
#                            = 1 + 0.5*cos(pi*t*sqrt(1)*1/m) )
#
#   Scenario "historical" : eta1(t)  = integral_0^t alpha3(s,t) eta2(s) ds + e(t)
#                           (alpha3(s,t) supplied by parameter_set1
#                            = 1 + 0.5*sin(pi*(s+t)*sqrt(1)*1/m) )
#
# This script keeps N, M and n.sim small so the run finishes quickly. Its purpose
# is to verify that the modified simulation/fitting pipeline supports each of the
# three eta1|eta2 dependency structures end to end (no MSE table is produced).

library(fda)
library(Matrix)
library(mvtnorm)
library(tmvtnorm)
library(MASS)

source("Core/FSEM.R")
source("Core/simualtion_utilities.R")

# ---- light settings -----------------------------------------------------------
n.sim   <- 2
N       <- 30   # small sample size for speed
M       <- 10   # few observation points per curve
n.b     <- 4    # small basis
n.em    <- 5    # very few EM iterations (smoke-test)
n.monte <- 20   # small Monte Carlo size in EM

# unbuffered logging to file + stdout
log_file <- "main_simulation_3_test.log"
if (file.exists(log_file)) file.remove(log_file)
msg <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s); flush.console()
  cat(s, file = log_file, append = TRUE)
}

scenarios <- c("fixed", "concurrent", "historical")
test_results <- list()

for (sc in scenarios) {
  msg("\n=========================================\n")
  msg("Scenario: eta1 ~ eta2 with effect =", sc, "\n")
  msg("=========================================\n")

  # ---- TRUE (data-generating) model -----------------------------------------
  # eta2 has indicators z4,z5 (concurrent) and z6 (historical) and intercept-only
  # regression. eta1 has indicators z1,z2 (concurrent) and z3 (historical) and a
  # regression on eta2 with the chosen functional effect type `sc`.
  # NOTE: each indicator has a single effect (concurrent OR historical) to match
  # the well-tested pattern in main_simulation.R.
  model.sim <-
    fsem(eta1 ~~ z1 + z2,      effectType = "concurrent") %+%
    fsem(eta1 ~~ z3,           effectType = "historical") %+%
    fsem(eta1 ~ -1 + eta2,     effectType = sc, scalar.covariate = FALSE,
         latent.covariate = "eta2") %+%
    fsem(eta2 ~~ z4 + z5,      effectType = "concurrent") %+%
    fsem(eta2 ~~ z6,           effectType = "historical") %+%
    fsem(eta2 ~ -1)

  # ---- FIT model (well specified: same eta1~eta2 dependency) ----------------
  model.fit <-
    fsem(eta1 ~~ z1,           effectType = "fixed_concurrent") %+%
    fsem(eta1 ~~ z2,           effectType = "concurrent") %+%
    fsem(eta1 ~~ z3,           effectType = "historical") %+%
    fsem(eta1 ~ -1 + eta2,     effectType = sc, scalar.covariate = FALSE,
         latent.covariate = "eta2") %+%
    fsem(eta2 ~~ z4,           effectType = "fixed_concurrent") %+%
    fsem(eta2 ~~ z5,           effectType = "concurrent") %+%
    fsem(eta2 ~~ z6,           effectType = "historical") %+%
    fsem(eta2 ~ -1)

  scenario_runs <- list()
  for (i in seq_len(n.sim)) {
    msg("  replicate", i, "/", n.sim, "... ")
    t0 <- proc.time()[3]

    sim_ok <- tryCatch({
      simm <- simulation(
        model      = model.sim,
        design     = "regular",
        n.t        = M,
        n.b        = n.b,
        r          = 1,
        rho        = 0.3,
        SNR        = 4,
        n.sample   = N,
        Matern.sem = TRUE,
        Matern.fac = FALSE,
        parameters = "parameter_set1"
      )
      TRUE
    }, error = function(e) { msg("SIM ERROR:", conditionMessage(e), "\n"); FALSE })

    if (!isTRUE(sim_ok)) next

    fit_ok <- tryCatch({
      initial.parameter <- initial.param(model = model.fit, data = simm$data, n.b = n.b)
      params.estimated  <- em.estimation(
        model             = model.fit,
        data              = simm$data,
        initial.parameter = initial.parameter,
        n.b               = n.b,
        n.em              = n.em,
        n.monte           = n.monte,
        s.p               = "min",
        range.min         = NULL,
        range.max         = NULL,
        design            = "regular"
      )
      TRUE
    }, error = function(e) { msg("FIT ERROR:", conditionMessage(e), "\n"); FALSE })

    dt <- proc.time()[3] - t0
    msg(sprintf("done in %.1fs (sim_ok=%s, fit_ok=%s)\n", dt, sim_ok, fit_ok))

    scenario_runs[[i]] <- list(
      sim_ok = sim_ok, fit_ok = fit_ok,
      elapsed = dt,
      estimation = if (isTRUE(fit_ok)) params.estimated else NULL
    )
  }
  test_results[[sc]] <- scenario_runs
}

# ---- summary ------------------------------------------------------------------
msg("\n\n========= TEST SUMMARY =========\n")
for (sc in scenarios) {
  runs <- test_results[[sc]]
  ok   <- sum(vapply(runs, function(r) isTRUE(r$sim_ok) && isTRUE(r$fit_ok), logical(1)))
  tt   <- sum(vapply(runs, function(r) if (is.null(r$elapsed)) 0 else r$elapsed, numeric(1)))
  msg(sprintf("  effect=%-11s  ok=%d/%d  total=%.1fs\n", sc, ok, length(runs), tt))
}
