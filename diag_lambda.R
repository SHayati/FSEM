## Localize the lambda problem: compare est vs true standardized lambda curves,
## raw lambda, and latent variance (sigma.sem) for the single run.
suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
  source("Core/FSEM.R"); source("Core/simualtion_utilities.R")
})
sr  <- readRDS("single_run_results.rds")
sim <- sr$sim.data
ev.t <- sim$params.org.eval
ev.e <- sr$params.estimated.well.specified$result$params.estimated.eval

cat("=== Standardized factor loadings: est vs true (per indicator) ===\n")
for (ind in names(ev.t$coef.fac.std)) {
  for (fac in names(ev.t$coef.fac.std[[ind]])) {
    tt <- ev.t$coef.fac.std[[ind]][[fac]]
    ee <- ev.e$coef.fac.std[[ind]][[fac]]
    if (is.null(ee)) { cat(sprintf("  %s<-%s : EST MISSING\n", ind, fac)); next }
    tv <- as.vector(tt); evv <- as.vector(ee)
    n <- min(length(tv), length(evv)); tv <- tv[1:n]; evv <- evv[1:n]
    corr <- suppressWarnings(cor(tv, evv))
    scale.ratio <- sqrt(mean(evv^2) / mean(tv^2))
    mse  <- mean((evv - tv)^2)
    mse.signflip <- mean((-evv - tv)^2)
    cat(sprintf("  %-4s<-%-4s  corr=%+.3f  RMS(est)/RMS(true)=%.3f  MSE=%.3f  MSE(sign-flip)=%.3f\n",
                ind, fac, corr, scale.ratio, mse, mse.signflip))
  }
}

cat("\n=== Latent innovation variance diag(sigma.sem eval), est vs true ===\n")
## true uses ker.sem (or sigma.sem.std diag); compare evaluated kernels
for (res in names(ev.t$coef.sem.std)) {
  # nothing; handled below via eigval
}
cat("sigma.fac eigenvalues (est vs true, first 4):\n")
for (i in seq_along(ev.t$sigma.fac.eigval.std)) {
  te <- ev.t$sigma.fac.eigval.std[[i]]
  es <- ev.e$sigma.fac.eigval.std[[i]]
  if (is.null(es)) next
  k <- min(4, length(te), length(es))
  cat(sprintf("  indicator %d  true: %s\n", i, paste(signif(te[1:k],3), collapse=" ")))
  cat(sprintf("               est : %s\n", paste(signif(es[1:k],3), collapse=" ")))
}

cat("\n=== sigma.sem eigenvalues (latent innovation), est vs true ===\n")
for (res in names(ev.t$sigma.sem.eigval.std)) {
  te <- ev.t$sigma.sem.eigval.std[[res]]
  es <- ev.e$sigma.sem.eigval.std[[res]]
  if (is.null(es)) next
  k <- min(4, length(te), length(es))
  cat(sprintf("  %-5s true: %s\n", res, paste(signif(te[1:k],3), collapse=" ")))
  cat(sprintf("        est : %s\n", paste(signif(es[1:k],3), collapse=" ")))
}
