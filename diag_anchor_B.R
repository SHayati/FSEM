## Confirm: anchored fit (fixed_concurrent) recovers standardized loadings on sim3 data.
suppressPackageStartupMessages({
  library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
  source("Core/FSEM.R"); source("Core/simualtion_utilities.R"); source("Core/FSEM_optimised.R")
})
set.seed(2024)
model.sim <- fsem(eta1~~z1+z2,effectType="concurrent")%+%
  fsem(eta1~~z3,effectType="historical")%+%
  fsem(eta2~~z4+z5+z6,effectType="concurrent")%+%
  fsem(eta1~-1+eta2,effectType="concurrent",latent.covariate="eta2",scalar.covariate=FALSE) %+%
  fsem(eta2~-1)
fit.anchor <- fsem(eta1~~z1,effectType="fixed_concurrent")%+%
  fsem(eta1~~z2,effectType="concurrent")%+%
  fsem(eta1~~z3,effectType="historical")%+%
  fsem(eta2~~z4,effectType="fixed_concurrent")%+%
  fsem(eta2~~z5+z6,effectType="concurrent")%+%
  fsem(eta1~-1+eta2,effectType="concurrent",latent.covariate="eta2",scalar.covariate=FALSE) %+%
  fsem(eta2~-1)
simm <- simulation(model=model.sim,design="regular",n.t=20,n.b=6,r=1,rho=0.3,SNR=4,
                   n.sample=50,Matern.sem=TRUE,Matern.fac=FALSE,parameters="parameter_set1")
report <- function(fit, label) {
  cf.e <- fit$result$params.estimated.eval$coef.fac.std
  cf.t <- simm$params.org.eval$coef.fac.std
  cat(sprintf("\n[%s]\n", label))
  for (ind in names(cf.t)) for (fac in names(cf.t[[ind]])) {
    tv <- as.vector(cf.t[[ind]][[fac]]); ev <- as.vector(cf.e[[ind]][[fac]])
    n <- min(length(tv), length(ev))
    cat(sprintf("  %-4s<-%-4s RMSratio=%.3f corr=%+.3f MSE=%.3f\n",
        ind, fac, sqrt(mean(ev[1:n]^2)/mean(tv[1:n]^2)),
        suppressWarnings(cor(tv[1:n],ev[1:n])), mean((ev[1:n]-tv[1:n])^2)))
  }
}
initB <- initial.param(model=fit.anchor, data=simm$data, n.b=6)
fitB <- em.estimation.optimised(model=fit.anchor, data=simm$data, initial.parameter=initB,
                                n.b=6, n.em=40, n.monte=30, s.p="min", design="regular", plot.progress=FALSE)
report(fitB, "WITH fixed_concurrent anchor")
cat("\nDONE\n")
