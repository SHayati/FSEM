## Try more EM iters + larger jitter to confirm convergence.
suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
  source("Core/FSEM.R"); source("Core/simualtion_utilities.R")
  source("Core/FSEM_optimised.R"); source("Core/FSEM_optimised_revised.R")
})
sr <- readRDS("single_run_results.rds"); sim <- sr$sim.data
model.fit.w <- fsem(eta1~~z1+z2, effectType="concurrent") %+%
  fsem(eta1~~z3, effectType="historical") %+%
  fsem(eta2~~z4+z5+z6, effectType="concurrent") %+%
  fsem(eta1~-1+eta2, effectType="concurrent", latent.covariate="eta2", scalar.covariate=FALSE) %+%
  fsem(eta2~-1)
keep.ids <- sort(unique(sim$data$.id))[1:50]
data.sub <- sim$data[sim$data$.id %in% keep.ids, ]; data.sub$.id <- as.integer(factor(data.sub$.id))
n.b <- 6
init <- initial.param(model.fit.w, n.b = n.b, data = data.sub)
true.g <- sim$params.org.eval$coef.sem.std[["eta1"]][["eta2"]]
true.l <- sim$params.org.eval$coef.fac.std
mse_lambda <- function(est, tru) {
  vals <- c()
  for (ii in seq_along(tru)) for (jj in seq_along(tru[[ii]])) {
    tv <- tru[[ii]][[jj]]; ev <- est[[ii]][[jj]]
    if (is.numeric(tv) && length(tv) == length(ev) && length(tv) > 0)
      vals <- c(vals, mean((ev - tv)^2))
  }; mean(vals, na.rm = TRUE)
}

NEM <- 60; NMC <- 30
for (sd in c(0, 0.3, 1.0, 2.0)) {
  cat(sprintf("\n--- gamma.init.sd = %.2f  (n.em=%d, n.monte=%d) ---\n", sd, NEM, NMC))
  t0 <- Sys.time()
  fit <- em.estimation.optimised.revised(
    model = model.fit.w, data = data.sub, x.data = NULL,
    initial.parameter = init, n.b = n.b,
    n.em = NEM, n.monte = NMC,
    range.min = 0, range.max = 1, design = "regular",
    plot.progress = FALSE, gamma.init.sd = sd, seed = 1
  )
  cat(sprintf("  fit took %.1fs\n", as.numeric(Sys.time() - t0)))
  g <- fit$result$params.estimated.eval$coef.sem.std[["eta1"]][["eta2"]]
  l <- fit$result$params.estimated.eval$coef.fac.std
  cat(sprintf("  gamma range [%6.3f, %6.3f]  mean|.|= %6.3f  MSE = %6.3f\n",
              min(g), max(g), mean(abs(g)), mean((g - true.g)^2)))
  cat(sprintf("  lambda MSE = %6.4f\n", mse_lambda(l, true.l)))
}
