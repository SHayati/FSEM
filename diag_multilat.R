## Decisive test: do factor loadings collapse with TWO INDEPENDENT latents
## (no structural coupling)? If yes, the bug is multi-latent E-step, not gamma.
suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
  source("Core/FSEM.R"); source("Core/simualtion_utilities.R"); source("Core/FSEM_optimised.R")
})
set.seed(123)

## (A) single latent (like sim2) -- baseline that is known to work
modelA <- fsem(eta1~~z1+z2+z3, effectType="concurrent") %+% fsem(eta1~-1)

## (B) two INDEPENDENT latents, no coupling
modelB <- fsem(eta1~~z1+z2+z3, effectType="concurrent") %+% fsem(eta1~-1) %+%
          fsem(eta2~~z4+z5+z6, effectType="concurrent") %+% fsem(eta2~-1)

rms_ratio <- function(fit, truth) {
  out <- c()
  cf.e <- fit$result$params.estimated.eval$coef.fac.std
  cf.t <- truth$params.org.eval$coef.fac.std
  for (ind in names(cf.t)) for (fac in names(cf.t[[ind]])) {
    tv <- as.vector(cf.t[[ind]][[fac]]); ev <- as.vector(cf.e[[ind]][[fac]])
    n <- min(length(tv), length(ev))
    out <- c(out, sqrt(mean(ev[1:n]^2)/mean(tv[1:n]^2)))
  }
  out
}

run_one <- function(model, label, N=40, M=10, NEM=30, NMC=20) {
  simm <- simulation(model=model, design="regular", n.t=M, n.b=6, r=1, rho=0.3,
                     SNR=4, n.sample=N, Matern.sem=TRUE, Matern.fac=FALSE,
                     parameters="parameter_set1")
  init <- initial.param(model=model, data=simm$data, n.b=6)
  fit <- em.estimation.optimised(model=model, data=simm$data, initial.parameter=init,
                                 n.b=6, n.em=NEM, n.monte=NMC, s.p="min", design="regular")
  r <- rms_ratio(fit, simm)
  cat(sprintf("\n[%s] RMS(est)/RMS(true) per loading (1.0 = perfect scale):\n  %s\n",
              label, paste(sprintf("%.3f", r), collapse="  ")))
  invisible(r)
}

run_one(modelA, "A: single latent")
run_one(modelB, "B: two INDEPENDENT latents")
cat("\nDONE\n")
