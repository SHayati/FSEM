## Verify em.estimation.optimised.revised fixes the gamma-stuck bug.
suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
  source("Core/FSEM.R")
  source("Core/simualtion_utilities.R")
  source("Core/FSEM_optimised.R")
  source("Core/FSEM_optimised_revised.R")
})

sr <- readRDS("single_run_results.rds")
sim <- sr$sim.data

## Rebuild model (same as main_simulation_3.R well-specified)
model.fit.w <- fsem(eta1~~z1+z2, effectType="concurrent") %+%
  fsem(eta1~~z3, effectType="historical") %+%
  fsem(eta2~~z4+z5+z6, effectType="concurrent") %+%
  fsem(eta1~-1+eta2, effectType="concurrent", latent.covariate="eta2", scalar.covariate=FALSE) %+%
  fsem(eta2~-1)
## Use only first 50 samples to keep runtime modest.
keep.ids <- sort(unique(sim$data$.id))[1:50]
data.sub <- sim$data[sim$data$.id %in% keep.ids, ]
data.sub$.id <- as.integer(factor(data.sub$.id))
n.b <- 6
init <- initial.param(model.fit.w, n.b = n.b, data = data.sub)

NEM <- 20; NMC <- 20

cat("Running ORIGINAL optimised (gamma=0 init)...\n")
t0 <- Sys.time()
fit.orig <- em.estimation.optimised(
  model = model.fit.w, data = data.sub, x.data = NULL,
  initial.parameter = init, n.b = n.b,
  n.em = NEM, n.monte = NMC,
  range.min = 0, range.max = 1, design = "regular",
  plot.progress = FALSE
)
cat("  took", round(as.numeric(Sys.time() - t0), 1), "s\n")

cat("\nRunning REVISED optimised (gamma jittered, seed=1)...\n")
t0 <- Sys.time()
fit.rev <- em.estimation.optimised.revised(
  model = model.fit.w, data = data.sub, x.data = NULL,
  initial.parameter = init, n.b = n.b,
  n.em = NEM, n.monte = NMC,
  range.min = 0, range.max = 1, design = "regular",
  plot.progress = FALSE,
  gamma.init.sd = 0.3, seed = 1
)
cat("  took", round(as.numeric(Sys.time() - t0), 1), "s\n")

true.g <- sim$params.org.eval$coef.sem.std[["eta1"]][["eta2"]]
g.orig <- fit.orig$result$params.estimated.eval$coef.sem.std[["eta1"]][["eta2"]]
g.rev  <- fit.rev$result$params.estimated.eval$coef.sem.std[["eta1"]][["eta2"]]

cat("\n=== Gamma eta1<-eta2 (curve over 200 grid points) ===\n")
cat(sprintf("True : range [%6.3f, %6.3f]   mean|.|= %6.3f\n",
            min(true.g), max(true.g), mean(abs(true.g))))
cat(sprintf("Orig : range [%6.3f, %6.3f]   mean|.|= %6.3f   MSE vs true = %6.3f\n",
            min(g.orig), max(g.orig), mean(abs(g.orig)), mean((g.orig - true.g)^2)))
cat(sprintf("Rev  : range [%6.3f, %6.3f]   mean|.|= %6.3f   MSE vs true = %6.3f\n",
            min(g.rev),  max(g.rev),  mean(abs(g.rev)),  mean((g.rev  - true.g)^2)))

## Also compare lambda MSE (factor loadings should also improve).
true.l <- sim$params.org.eval$coef.fac.std
est.l.o <- fit.orig$result$params.estimated.eval$coef.fac.std
est.l.r <- fit.rev$result$params.estimated.eval$coef.fac.std

mse_lambda <- function(est, tru) {
  vals <- c()
  for (ii in seq_along(tru)) {
    for (jj in seq_along(tru[[ii]])) {
      tv <- tru[[ii]][[jj]]; ev <- est[[ii]][[jj]]
      if (is.numeric(tv) && length(tv) == length(ev) && length(tv) > 0)
        vals <- c(vals, mean((ev - tv)^2))
    }
  }
  mean(vals, na.rm = TRUE)
}
cat(sprintf("\nLambda MSE  Orig: %.4f    Rev: %.4f\n",
            mse_lambda(est.l.o, true.l), mse_lambda(est.l.r, true.l)))
