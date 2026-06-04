## verify_optimised.R  -- correctness + timing check for em.estimation.optimised
suppressMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
})
source("Core/FSEM.R")
source("Core/simualtion_utilities.R")
source("Core/FSEM_optimised.R")

set.seed(20240520)
model.sim <- fsem(eta1~~z1+z2, effectType="concurrent") %+%
  fsem(eta1~~z3, effectType="historical") %+%
  fsem(eta1~-1)
model.fit <- fsem(eta1~~z1, effectType="fixed_concurrent") %+%
  fsem(eta1~~z2, effectType="concurrent") %+%
  fsem(eta1~~z3, effectType="historical") %+%
  fsem(eta1~-1)

simm <- simulation(model=model.sim, design="regular", n.t=10, n.b=6, r=1, rho=0.3,
                   SNR=4, n.sample=25, Matern.sem=TRUE, Matern.fac=FALSE,
                   parameters="parameter_set1")
initial.parameter <- initial.param(model=model.fit, data=simm$data, n.b=6)

SEED <- 999
## n.em=1 gives a clean correctness check: the eta draws are bit-identical
## (same RNG stream + identical sampling arithmetic), so any remaining
## difference is purely the algebraic reformulation (~1e-10). For n.em>1 the
## two runs legitimately diverge because mvrnorm()'s eigen step amplifies
## ~1e-12 round-off into different (equally valid) random draws.
NEM <- 1; NMC <- 3

set.seed(SEED)
t0 <- Sys.time()
orig <- em.estimation(model=model.fit, data=simm$data, initial.parameter=initial.parameter,
                      n.b=6, n.em=NEM, n.monte=NMC, s.p="min", design="regular")
t_orig <- as.numeric(difftime(Sys.time(), t0, units="secs"))

set.seed(SEED)
t0 <- Sys.time()
opt <- em.estimation.optimised(model=model.fit, data=simm$data, initial.parameter=initial.parameter,
                               n.b=6, n.em=NEM, n.monte=NMC, s.p="min", design="regular")
t_opt <- as.numeric(difftime(Sys.time(), t0, units="secs"))

cat("\n=== Rcpp kernels in use:", isTRUE(.fsem_have_rcpp), "===\n")
cat(sprintf("timing: original=%.2fs  optimised=%.2fs  speedup=%.2fx\n",
            t_orig, t_opt, t_orig / t_opt))

## compare standardized evaluated estimates
ev1 <- orig$result$params.estimated.eval
ev2 <- opt$result$params.estimated.eval

rdiff <- function(a, b) {
  a <- unlist(a); b <- unlist(b)
  if (length(a) != length(b)) return(NA_real_)
  max(abs(a - b))
}
flds <- intersect(names(ev1), names(ev2))
cat("\n=== max abs diff per evaluated field (original vs optimised) ===\n")
for (f in flds) {
  d <- tryCatch(rdiff(ev1[[f]], ev2[[f]]), error=function(e) NA_real_)
  cat(sprintf("  %-28s %.3e\n", f, d))
}

p1 <- orig$result$params.estimated$totparam
p2 <- opt$result$params.estimated$totparam
cat("\n=== max abs diff: lambda / gamma / sigma.error ===\n")
cat(sprintf("  lambda.param  %.3e\n", rdiff(p1$lambda.param, p2$lambda.param)))
cat(sprintf("  gamma.param   %.3e\n", rdiff(p1$gamma.param,  p2$gamma.param)))
cat(sprintf("  sigma.error   %.3e\n", rdiff(p1$sigma.error,  p2$sigma.error)))
cat("\nDONE\n")
