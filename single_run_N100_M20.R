## single_run_N100_M20.R
## ---------------------------------------------------------------------------
## Single-run FSEM estimation for the sim3 (two dependent latents) setting at
## N = 100, M = 20, and a side-by-side performance report.
##
## What it does:
##   * simulates ONE data set from the sim3 generating model;
##   * fits TWO models:
##       (A) "naive"    = main_simulation_3.R's model.fit.w (all `concurrent`,
##                        NO scale anchor)  -> reproduces the poor performance;
##       (B) "anchored" = same model but with `fixed_concurrent` pinning the
##                        first indicator of each latent (the lambda fix), and
##                        estimated with em.estimation.optimised.revised so the
##                        structural gamma also escapes the gamma=0 trap.
##   * prints per-coefficient metrics (RMS ratio, correlation, MSE) for the
##     standardized loadings (lambda) and standardized structural coef (gamma),
##     plus the residual-covariance eigenvalues.
##
## NOTE: run this yourself, e.g.
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" single_run_N100_M20.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
})

source("Core/FSEM.R")
source("Core/simualtion_utilities.R")
source("Core/FSEM_optimised.R")
source("Core/fsem_kernels_fix.R")         # collision-free Rcpp kernel loader (fixes Error 322)
source("Core/FSEM_optimised_revised.R")   # provides em.estimation.optimised.revised

## ---- settings -------------------------------------------------------------
N      <- 100      # number of subjects
M      <- 20       # number of observation points per curve
n.b    <- 6        # number of basis functions
n.em   <- 100      # EM iterations
n.monte<- 100      # Monte-Carlo draws in the E-step
set.seed(2024)     # reproducible single run

## ---- generating model (sim3) ----------------------------------------------
model.sim <- fsem(eta1~~z1+z2, effectType="concurrent")          %+%
  fsem(eta1~~z3, effectType="historical")                        %+%
  fsem(eta2~~z4+z5+z6, effectType="concurrent")                  %+%
  fsem(eta1~-1+eta2, effectType="concurrent",
       latent.covariate="eta2", scalar.covariate=FALSE)          %+%
  fsem(eta2~-1)

## ---- (A) naive fit: exactly main_simulation_3.R's model.fit.w --------------
fit.naive <- model.sim

## ---- (B) anchored fit: pin one loading per latent (the lambda fix) ---------
fit.anchor <- fsem(eta1~~z1, effectType="fixed_concurrent")      %+%
  fsem(eta1~~z2, effectType="concurrent")                        %+%
  fsem(eta1~~z3, effectType="historical")                        %+%
  fsem(eta2~~z4, effectType="fixed_concurrent")                  %+%
  fsem(eta2~~z5+z6, effectType="concurrent")                     %+%
  fsem(eta1~-1+eta2, effectType="concurrent",
       latent.covariate="eta2", scalar.covariate=FALSE)          %+%
  fsem(eta2~-1)

## ---- simulate one data set -------------------------------------------------
cat(sprintf("Simulating one data set: N=%d, M=%d, n.b=%d ...\n", N, M, n.b))
simm <- simulation(model=model.sim, design="regular", n.t=M, n.b=n.b, r=1,
                   rho=0.3, SNR=4, n.sample=N, Matern.sem=TRUE, Matern.fac=FALSE,
                   parameters="parameter_set1")

## ---------------------------------------------------------------------------
## reporting helpers
## ---------------------------------------------------------------------------
.metrics <- function(tv, ev) {
  n  <- min(length(tv), length(ev))
  tv <- tv[1:n]; ev <- ev[1:n]
  list(rmsratio = sqrt(mean(ev^2) / mean(tv^2)),
       corr     = suppressWarnings(cor(tv, ev)),
       mse      = mean((ev - tv)^2))
}

report.coef <- function(fit.eval, true.eval, field, title) {
  cf.e <- fit.eval[[field]]
  cf.t <- true.eval[[field]]
  cat(sprintf("\n  -- %s (%s) --\n", title, field))
  if (is.null(cf.t) || length(cf.t) == 0) { cat("    (none)\n"); return(invisible()) }
  for (res in names(cf.t)) {
    for (cov in names(cf.t[[res]])) {
      tv <- as.vector(cf.t[[res]][[cov]])
      ev <- as.vector(cf.e[[res]][[cov]])
      if (is.null(ev)) next
      m <- .metrics(tv, ev)
      cat(sprintf("    %-5s <- %-5s  RMSratio=%6.3f  corr=%+.3f  MSE=%.4f\n",
                  res, cov, m$rmsratio, m$corr, m$mse))
    }
  }
}

# report.eig <- function(fit.eval, true.eval, field) {
#   ev <- fit.eval[[field]]; tv <- true.eval[[field]]
#   if (is.null(ev) || is.null(tv)) return(invisible())
#   cat(sprintf("\n  -- %s --\n", field))
#   cat(sprintf("    true: %s\n", paste(sprintf("%.3f", as.vector(tv)), collapse=", ")))
#   cat(sprintf("    est : %s\n", paste(sprintf("%.3f", as.vector(ev)), collapse=", ")))
# }

report.eig <- function(fit.eval, true.eval, field) {
  ev <- fit.eval[[field]]; tv <- true.eval[[field]]
  if (is.null(ev) || is.null(tv)) return(invisible())
  cat(sprintf("\n  -- %s --\n", field))
  # Handle both numeric and list-of-numeric
  print_eig <- function(x, label) {
    if (is.list(x)) {
      for (i in seq_along(x)) {
        cat(sprintf("    %s[%d]: %s\n", label, i, paste(sprintf("%.3f", as.vector(x[[i]])), collapse=", ")))
      }
    } else {
      cat(sprintf("    %s: %s\n", label, paste(sprintf("%.3f", as.vector(x)), collapse=", " )))
    }
  }
  print_eig(tv, "true")
  print_eig(ev, "est ")
}

report.all <- function(fit, label) {
  fe <- fit$result$params.estimated.eval
  te <- simm$params.org.eval
  cat(sprintf("\n================ %s ================\n", label))
  report.coef(fe, te, "coef.fac.std", "standardized loadings  (lambda)")
  report.coef(fe, te, "coef.sem.std", "standardized structural (gamma)")
  report.eig (fe, te, "sigma.fac.eigval.std")
  report.eig (fe, te, "sigma.sem.eigval.std")
}

## ---------------------------------------------------------------------------
## (A) naive fit  (original em.estimation.optimised, no anchor)
## ---------------------------------------------------------------------------
cat("\nFitting (A) NAIVE model (no scale anchor) ...\n")
init.A <- initial.param(model=fit.naive, data=simm$data, n.b=n.b)
est.A  <- em.estimation.optimised(model=fit.naive, data=simm$data,
                                  initial.parameter=init.A, n.b=n.b,
                                  n.em=n.em, n.monte=n.monte, s.p="min",
                                  design="regular", plot.progress=FALSE)

## ---------------------------------------------------------------------------
## (B) anchored fit  (revised estimator: gamma-jitter + scale anchor)
## ---------------------------------------------------------------------------
cat("\nFitting (B) ANCHORED model (fixed_concurrent + gamma jitter) ...\n")
init.B <- initial.param(model=fit.anchor, data=simm$data, n.b=n.b)
est.B  <- em.estimation.optimised.revised(model=fit.anchor, data=simm$data,
                                          initial.parameter=init.B, n.b=n.b,
                                          n.em=n.em, n.monte=n.monte, s.p="min",
                                          design="regular", plot.progress=FALSE,
                                          gamma.init.sd=0.3, seed=1)

## ---------------------------------------------------------------------------
## report
## ---------------------------------------------------------------------------
report.all(est.A, "A) NAIVE  (model.fit.w from main_simulation_3.R)")
report.all(est.B, "B) ANCHORED + revised (lambda fix + gamma fix)")

cat("\nInterpretation: for well-recovered parameters RMSratio should be near 1.0,\n")
cat("corr near +1, and MSE small. The naive fit (A) collapses the loadings\n")
cat("(RMSratio ~ 0); the anchored fit (B) restores them.\n")
cat("\nDONE\n")
