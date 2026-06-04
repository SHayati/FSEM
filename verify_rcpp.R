## Verify em.estimation.rcpp matches em.estimation at n.em=1, then benchmark
## all three implementations (original / optimised / rcpp).

suppressMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
})
source("Core/FSEM.R")
source("Core/simualtion_utilities.R")
source("Core/FSEM_optimised.R")
source("Core/FSEM_rcpp.R")

## -- Small seeded simulation matching verify_optimised.R for correctness ---
make_models <- function() {
  model.sim <- fsem(eta1 ~~ z1 + z2, effectType = "concurrent") %+%
               fsem(eta1 ~~ z3, effectType = "historical") %+%
               fsem(eta1 ~ -1)
  model.fit <- fsem(eta1 ~~ z1, effectType = "fixed_concurrent") %+%
               fsem(eta1 ~~ z2, effectType = "concurrent") %+%
               fsem(eta1 ~~ z3, effectType = "historical") %+%
               fsem(eta1 ~ -1)
  list(sim = model.sim, fit = model.fit)
}

NB  <- 6
SEED <- 4242

cat("=== Correctness check: n.em = 1, small N ===\n")
mods <- make_models()
set.seed(SEED)
simm  <- simulation(model = mods$sim, design = "regular", n.t = 10, n.b = NB,
                    r = 1, rho = 0.3, SNR = 4, n.sample = 25,
                    Matern.sem = TRUE, Matern.fac = FALSE, parameters = "parameter_set1")
ip <- initial.param(model = mods$fit, data = simm$data, n.b = NB)

set.seed(SEED)
t0 <- Sys.time()
res.orig <- em.estimation(model = mods$fit, data = simm$data,
                          initial.parameter = ip, n.b = NB,
                          n.em = 1, n.monte = 3, s.p = "min",
                          design = "regular")
t.orig.sm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

set.seed(SEED)
t0 <- Sys.time()
res.rcpp <- em.estimation.rcpp(model = mods$fit, data = simm$data,
                               initial.parameter = ip, n.b = NB,
                               n.em = 1, n.monte = 3, s.p = "min",
                               design = "regular")
t.rcpp.sm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("\nRcpp kernels loaded:", isTRUE(.fsem_have_rcpp_full), "\n")
cat(sprintf("small-scale time: original=%.2fs  rcpp=%.2fs\n\n", t.orig.sm, t.rcpp.sm))

## Compare params.estimated.eval fields field-by-field
fields <- names(res.orig$result$params.estimated.eval)
cat("=== max abs diff per evaluated field (original vs rcpp) ===\n")
for (fn in fields) {
  a <- res.orig$result$params.estimated.eval[[fn]]
  b <- res.rcpp$result$params.estimated.eval[[fn]]
  d <- tryCatch(max(abs(unlist(a) - unlist(b)), na.rm = TRUE), error = function(e) NA)
  cat(sprintf("  %-30s %.3e\n", fn, d))
}
cat("\n=== max abs diff: lambda / sigma.error / gamma ===\n")
la <- unlist(res.orig$result$params.estimated$totparam$lambda.param)
lb <- unlist(res.rcpp$result$params.estimated$totparam$lambda.param)
cat(sprintf("  lambda.param  %.3e\n", max(abs(la - lb))))
ga <- unlist(res.orig$result$params.estimated$totparam$gamma.param)
gb <- unlist(res.rcpp$result$params.estimated$totparam$gamma.param)
cat(sprintf("  gamma.param   %.3e\n", if (length(ga) > 0) max(abs(ga - gb)) else 0))
sa <- unlist(res.orig$result$params.estimated$totparam$sigma.error)
sb <- unlist(res.rcpp$result$params.estimated$totparam$sigma.error)
cat(sprintf("  sigma.error   %.3e\n", max(abs(sa - sb))))

## -- Realistic-scale benchmark ---------------------------------------------
cat("\n=== Benchmark: N=100 M=20 n.b=10 n.monte=20 n.em=1 ===\n")
set.seed(SEED)
simm  <- simulation(model = mods$sim, design = "regular", n.t = 20, n.b = 10,
                    r = 1, rho = 0.3, SNR = 4, n.sample = 100,
                    Matern.sem = TRUE, Matern.fac = FALSE, parameters = "parameter_set1")
ip <- initial.param(model = mods$fit, data = simm$data, n.b = 10)

set.seed(SEED); t0 <- Sys.time()
b.orig <- em.estimation(model = mods$fit, data = simm$data,
                        initial.parameter = ip, n.b = 10,
                        n.em = 1, n.monte = 20, s.p = "min",
                        design = "regular")
t.orig <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

set.seed(SEED); t0 <- Sys.time()
b.opt  <- em.estimation.optimised(model = mods$fit, data = simm$data,
                                  initial.parameter = ip, n.b = 10,
                                  n.em = 1, n.monte = 20, s.p = "min",
                                  design = "regular")
t.opt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

set.seed(SEED); t0 <- Sys.time()
b.rcpp <- em.estimation.rcpp(model = mods$fit, data = simm$data,
                             initial.parameter = ip, n.b = 10,
                             n.em = 1, n.monte = 20, s.p = "min",
                             design = "regular")
t.rcpp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("\noriginal   = %7.2fs\n",   t.orig))
cat(sprintf("optimised  = %7.2fs   speedup vs original = %.2fx\n", t.opt,  t.orig / t.opt))
cat(sprintf("rcpp       = %7.2fs   speedup vs original = %.2fx   vs optimised = %.2fx\n",
            t.rcpp, t.orig / t.rcpp, t.opt / t.rcpp))

la <- unlist(b.orig$result$params.estimated$totparam$lambda.param)
lr <- unlist(b.rcpp$result$params.estimated$totparam$lambda.param)
lo <- unlist(b.opt$result$params.estimated$totparam$lambda.param)
cat(sprintf("\nlambda max abs diff (rcpp vs original)    = %.3e\n", max(abs(la - lr))))
cat(sprintf("lambda max abs diff (optimised vs original)= %.3e\n", max(abs(la - lo))))

cat("\nDONE\n")
