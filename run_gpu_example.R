## run_gpu_example.R
## ---------------------------------------------------------------------------
## Minimal runnable example for em.estimation.optimised.gpu() on an RTX GPU.
##
## Same model/settings as single_run_N100_M20.R's anchored fit (B), but the
## structural+factor M-step runs on the GPU (CUDA) instead of the CPU.
##
## Prerequisite (one-time): install the CUDA libtorch backend. For an RTX
## 50-series (Blackwell, sm_120) GPU you MUST pin CUDA 12.8 and raise the
## download timeout (the backend is ~3.4 GB):
##
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" -e \
##     'options(timeout=7200); Sys.setenv(CUDA="12.8"); \
##      torch::install_torch(type="cuda")'
##
## Run:
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" run_gpu_example.R
## ---------------------------------------------------------------------------

Sys.setenv(CUDA = "12.8")   # required so torch loads the cu128 backend on Blackwell

suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
})

source("Core/FSEM.R")
source("Core/simualtion_utilities.R")
source("Core/FSEM_optimised.R")
source("Core/fsem_kernels_fix.R")
source("Core/FSEM_optimised_revised.R")
source("Core/FSEM_gpu.R")               # provides em.estimation.optimised.gpu

## ---- settings -------------------------------------------------------------
N       <- 100
M       <- 20
n.b     <- 6
n.em    <- 100
n.monte <- 100
set.seed(2024)

## ---- generating + anchored fit model (sim3, two dependent latents) --------
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

cat(sprintf("Simulating one data set: N=%d, M=%d, n.b=%d ...\n", N, M, n.b))
simm <- simulation(model=model.sim, design="regular", n.t=M, n.b=n.b, r=1,
                   rho=0.3, SNR=4, n.sample=N, Matern.sem=TRUE, Matern.fac=FALSE,
                   parameters="parameter_set1")

init <- initial.param(model=fit.anchor, data=simm$data, n.b=n.b)

## ---- GPU fit --------------------------------------------------------------
cat("\nFitting on GPU (em.estimation.optimised.gpu) ...\n")
est.gpu <- em.estimation.optimised.gpu(
  model=fit.anchor, data=simm$data, initial.parameter=init, n.b=n.b,
  n.em=n.em, n.monte=n.monte, s.p="min", design="regular",
  plot.progress=FALSE, gamma.init.sd=0.3, seed=1,
  gpu.device="cuda")            # set gpu.device="cpu" to force CPU torch tensors

## ---- quick recovery report vs the truth -----------------------------------
.metrics <- function(tv, ev) {
  n <- min(length(tv), length(ev)); tv <- tv[1:n]; ev <- ev[1:n]
  c(rmsratio = sqrt(mean(ev^2)/mean(tv^2)),
    corr     = suppressWarnings(cor(tv, ev)))
}
report <- function(field, title) {
  fe <- est.gpu$result$params.estimated.eval[[field]]
  te <- simm$params.org.eval[[field]]
  cat(sprintf("\n  -- %s --\n", title))
  for (res in names(te)) for (cov in names(te[[res]])) {
    ev <- as.vector(fe[[res]][[cov]]); tv <- as.vector(te[[res]][[cov]])
    if (is.null(ev)) next
    m <- .metrics(tv, ev)
    cat(sprintf("    %-5s <- %-6s  RMSratio=%6.3f  corr=%+.3f\n",
                res, cov, m["rmsratio"], m["corr"]))
  }
}
report("coef.fac.std", "standardized loadings  (lambda)")
report("coef.sem.std", "standardized structural (gamma)")

cat("\nDONE\n")
