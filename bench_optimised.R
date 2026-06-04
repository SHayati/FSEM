## bench_optimised.R -- timing benchmark at a realistic scale (single EM iter)
suppressMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
})
source("Core/FSEM.R")
source("Core/simualtion_utilities.R")
source("Core/FSEM_optimised.R")

set.seed(7)
model.sim <- fsem(eta1~~z1+z2, effectType="concurrent") %+%
  fsem(eta1~~z3, effectType="historical") %+% fsem(eta1~-1)
model.fit <- fsem(eta1~~z1, effectType="fixed_concurrent") %+%
  fsem(eta1~~z2, effectType="concurrent") %+%
  fsem(eta1~~z3, effectType="historical") %+% fsem(eta1~-1)

N <- 100; M <- 20; NB <- 10; NMC <- 20
simm <- simulation(model=model.sim, design="regular", n.t=M, n.b=6, r=1, rho=0.3,
                   SNR=4, n.sample=N, Matern.sem=TRUE, Matern.fac=FALSE,
                   parameters="parameter_set1")
ip <- initial.param(model=model.fit, data=simm$data, n.b=NB)

cat(sprintf("Benchmark: N=%d M=%d n.b=%d n.monte=%d n.em=1\n", N, M, NB, NMC))
cat("Rcpp kernels:", isTRUE({.fsem_load_kernels(); .fsem_have_rcpp}), "\n\n")

set.seed(11); t0 <- Sys.time()
o <- em.estimation(model=model.fit, data=simm$data, initial.parameter=ip,
                   n.b=NB, n.em=1, n.monte=NMC, s.p="min", design="regular")
to <- as.numeric(difftime(Sys.time(), t0, units="secs"))

set.seed(11); t0 <- Sys.time()
p <- em.estimation.optimised(model=model.fit, data=simm$data, initial.parameter=ip,
                             n.b=NB, n.em=1, n.monte=NMC, s.p="min", design="regular")
tp <- as.numeric(difftime(Sys.time(), t0, units="secs"))

cat(sprintf("\noriginal   = %.2fs\noptimised  = %.2fs\nspeedup    = %.2fx\n", to, tp, to/tp))
d <- max(abs(unlist(o$result$params.estimated$totparam$lambda.param) -
             unlist(p$result$params.estimated$totparam$lambda.param)))
cat(sprintf("lambda max abs diff = %.3e\n", d))
