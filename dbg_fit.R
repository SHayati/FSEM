options(warn = 1)
library(fda); library(Matrix); library(mvtnorm); library(tmvtnorm); library(MASS)
source("Core/FSEM.R"); source("Core/simualtion_utilities.R")

N <- 30; M <- 10; n.b <- 4

sc <- "historical"
model.sim <-
  fsem(eta1 ~~ z1 + z2,      effectType = "concurrent") %+%
  fsem(eta1 ~~ z3,           effectType = "historical") %+%
  fsem(eta1 ~ -1 + eta2,     effectType = sc, scalar.covariate = FALSE,
       latent.covariate = "eta2") %+%
  fsem(eta2 ~~ z4 + z5,      effectType = "concurrent") %+%
  fsem(eta2 ~~ z6,           effectType = "historical") %+%
  fsem(eta2 ~ -1)

model.fit <-
  fsem(eta1 ~~ z1,           effectType = "fixed_concurrent") %+%
  fsem(eta1 ~~ z2,           effectType = "concurrent") %+%
  fsem(eta1 ~~ z3,           effectType = "historical") %+%
  fsem(eta1 ~ -1 + eta2,     effectType = sc, scalar.covariate = FALSE,
       latent.covariate = "eta2") %+%
  fsem(eta2 ~~ z4,           effectType = "fixed_concurrent") %+%
  fsem(eta2 ~~ z5,           effectType = "concurrent") %+%
  fsem(eta2 ~~ z6,           effectType = "historical") %+%
  fsem(eta2 ~ -1)

cat("model.fit$mod$regression:\n")
for (nm in names(model.fit$mod$regression)) {
  r <- model.fit$mod$regression[[nm]]
  cat(" ", nm, " resp=", r$response, " cov=[", paste(r$covariate, collapse=","),
      "] eff=[", paste(r$effect, collapse=","), "]\n", sep="")
}

set.seed(1)
simm <- simulation(model=model.sim, design="regular", n.t=M, n.b=n.b, r=1, rho=0.3, SNR=4,
                   n.sample=N, Matern.sem=TRUE, Matern.fac=FALSE, parameters="parameter_set1")
cat("sim ok\n")

initial.parameter <- tryCatch(
  initial.param(model = model.fit, data = simm$data, n.b = n.b),
  error = function(e) { cat("initial.param ERROR:", conditionMessage(e), "\n"); traceback(); NULL })

if (is.null(initial.parameter)) quit(status=1)
cat("initial.param ok\n")
str(initial.parameter, max.level=2)

cat("\n--- em.estimation with traceback ---\n")
options(error = function() { traceback(3); quit(status=2) })
res <- em.estimation(model=model.fit, data=simm$data, initial.parameter=initial.parameter,
              n.b=n.b, n.em=1, n.monte=10, s.p="min",
              range.min=NULL, range.max=NULL, design="regular")

cat("done\n")
