## verify_gpu.R
## ---------------------------------------------------------------------------
## Validate em.estimation.optimised.gpu against em.estimation.optimised.revised
## on a small sim3 (two dependent latents) problem.
##
## Results are NOT expected to be bit-identical (GPU reduction order + RNG),
## but the recovered standardized loadings/structural coefficients should be
## close (high correlation, RMSE ratio ~ 1).
##
## Run:
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" verify_gpu.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
})

source("Core/FSEM.R")
source("Core/simualtion_utilities.R")
source("Core/FSEM_optimised.R")
source("Core/fsem_kernels_fix.R")
source("Core/FSEM_optimised_revised.R")
source("Core/FSEM_gpu.R")

## ---- small settings for a quick comparison --------------------------------
N       <- 50
M       <- 10
n.b     <- 6
n.em    <- 8
n.monte <- 20
set.seed(2024)

model.sim <- fsem(eta1~~z1+z2, effectType="concurrent")          %+%
  fsem(eta1~~z3, effectType="historical")                        %+%
  fsem(eta2~~z4+z5+z6, effectType="concurrent")                  %+%
  fsem(eta1~-1+eta2, effectType="concurrent",
       latent.covariate="eta2", scalar.covariate=FALSE)          %+%
  fsem(eta2~-1)

fit.anchor <- fsem(eta1~~z1, effectType="fixed_concurrent")      %+%
  fsem(eta1~~z2, effectType="concurrent")                        %+%
  fsem(eta1~~z3, effectType="historical")                        %+%
  fsem(eta2~~z4, effectType="fixed_concurrent")                  %+%
  fsem(eta2~~z5+z6, effectType="concurrent")                     %+%
  fsem(eta1~-1+eta2, effectType="concurrent",
       latent.covariate="eta2", scalar.covariate=FALSE)          %+%
  fsem(eta2~-1)

cat(sprintf("Simulating: N=%d, M=%d, n.b=%d (n.em=%d, n.monte=%d)\n",
            N, M, n.b, n.em, n.monte))
simm <- simulation(model=model.sim, design="regular", n.t=M, n.b=n.b, r=1,
                   rho=0.3, SNR=4, n.sample=N, Matern.sem=TRUE, Matern.fac=FALSE,
                   parameters="parameter_set1")

init <- initial.param(model=fit.anchor, data=simm$data, n.b=n.b)

cat("\n=== CPU: em.estimation.optimised.revised ===\n")
t.cpu <- system.time(
  est.cpu <- em.estimation.optimised.revised(
    model=fit.anchor, data=simm$data, initial.parameter=init, n.b=n.b,
    n.em=n.em, n.monte=n.monte, s.p="min", design="regular",
    plot.progress=FALSE, gamma.init.sd=0.3, seed=1)
)
cat(sprintf("CPU time: %.1f s\n", t.cpu[["elapsed"]]))

cat("\n=== GPU: em.estimation.optimised.gpu ===\n")
t.gpu <- system.time(
  est.gpu <- em.estimation.optimised.gpu(
    model=fit.anchor, data=simm$data, initial.parameter=init, n.b=n.b,
    n.em=n.em, n.monte=n.monte, s.p="min", design="regular",
    plot.progress=FALSE, gamma.init.sd=0.3, seed=1)
)
cat(sprintf("GPU time: %.1f s\n", t.gpu[["elapsed"]]))

## ---- compare standardized estimates ---------------------------------------
cmp <- function(field, label) {
  a <- est.cpu$result$params.estimated.eval[[field]]
  b <- est.gpu$result$params.estimated.eval[[field]]
  te <- simm$params.org.eval[[field]]
  cat(sprintf("\n-- %s --\n", label))
  for (res in names(a)) {
    for (cov in names(a[[res]])) {
      av <- as.vector(a[[res]][[cov]]); bv <- as.vector(b[[res]][[cov]])
      if (is.null(bv) || length(av) == 0) next
      n <- min(length(av), length(bv))
      d <- sqrt(mean((av[1:n] - bv[1:n])^2))
      rr <- sqrt(mean(bv[1:n]^2) / mean(av[1:n]^2))
      cc <- suppressWarnings(cor(av[1:n], bv[1:n]))
      tv <- if (!is.null(te[[res]][[cov]])) as.vector(te[[res]][[cov]]) else NULL
      tcc <- if (!is.null(tv)) suppressWarnings(cor(tv[1:min(n,length(tv))], bv[1:min(n,length(tv))])) else NA
      cat(sprintf("   %-5s <- %-6s  RMSE(cpu,gpu)=%.4f  ratio=%.3f  corr(cpu,gpu)=%+.4f  corr(gpu,true)=%+.3f\n",
                  res, cov, d, rr, cc, tcc))
    }
  }
}

cmp("coef.fac.std", "standardized loadings  (lambda)  CPU vs GPU")
cmp("coef.sem.std", "standardized structural (gamma)  CPU vs GPU")

cat("\nDONE\n")
