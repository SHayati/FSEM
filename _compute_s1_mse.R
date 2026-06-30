## Stream S1 (regular.missing, N=100, M=8) files -> MSE of lambda_1,2,3 (factor loadings)
## and gamma_1^x, gamma_2^x (regression coefs). Mirrors Core/MSE_tab3.R (trim=0.1 mean,
## 1/200 normalization for linear/concurrent, sign-flip not needed for these columns).
setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
files <- list.files("outputs", pattern="^result_parameter_set1_regular\\.missing_[0-9]+_N_100_M_20\\.rds$", full.names=TRUE)
cat("S1 files:", length(files), "\n")

lam <- matrix(NA, length(files), 3)   # z1,z2,z3 loadings on eta
gam <- matrix(NA, length(files), 2)   # x1,x2 regression coefs

for (i in seq_along(files)) {
  res <- readRDS(files[i])
  est <- res$estimation$result$params.estimated.eval
  org <- res$simulation$params.org.eval
  nam.fac <- res$model.fit$var$indicators           # z1 z2 z3
  nam.sem <- res$model.fit$var$latents              # eta
  # factor loadings: each indicator j has a list over factors; take the eta loading
  for (j in 1:3) {
    ce <- est$coef.fac.std[[j]]; co <- org$coef.fac.std[[j]]
    if (length(ce) > 0) {
      # effect: z1 fixed_concurrent, z2/z3 concurrent -> all linear-type, 1/200 norm
      k <- 1
      lam[i, j] <- 1/200 * sum((ce[[k]] - co[[k]])^2)
    }
  }
  # regression coefs eta ~ x1 + x2 (linear) -> coef.sem.std$eta$x1,$x2
  cse <- est$coef.sem.std[[nam.sem[1]]]; cso <- org$coef.sem.std[[nam.sem[1]]]
  for (k in 1:2) {
    gam[i, k] <- 1/200 * sum((cse[[k]] - cso[[k]])^2)
  }
  rm(res, est, org); if (i %% 20 == 0) gc()
}

mse <- c(apply(lam, 2, function(x) mean(x, trim=0.1, na.rm=TRUE)),
         apply(gam, 2, function(x) mean(x, trim=0.1, na.rm=TRUE)))
names(mse) <- c("lambda1","lambda2","lambda3","gamma1x","gamma2x")
cat("\n===== TABLE S1 MSE (trim=0.1) =====\n")
print(round(mse, 4))
saveRDS(list(lam=lam, gam=gam, mse=mse), "_s1_mse_acc.rds")
cat("DONE\n")
