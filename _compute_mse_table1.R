## Streaming MSE computation for Table 1 (main manuscript, Section 3)
## Loads each per-iteration result file ONE AT A TIME (memory efficient).
## Replicates Core/MSE_tab1.R logic (plain mean over Monte Carlo runs,
## standardized/identifiable parameters, eigenvector sign-flip).
setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
suppressMessages({library(Matrix)})
options(width=220)

designs <- c("regular","irregular")
combs <- expand.grid(N=c(50,100), M=c(10,20))
combs <- combs[order(combs$N),]

## accumulator: keyed by design|N|M ; stores running sums of per-run MSEs + count
acc <- list()
keyf <- function(d,N,M) paste(d,N,M,sep="|")

per_file_mse <- function(res){
  nam.fac <- res$model.fit$var$indicators      # z1 z2 z3
  nam.sem <- res$model.fit$var$latents          # eta1
  n.coef.fac <- length(nam.fac)
  est  <- res$estimation$result$params.estimated.eval
  org  <- res$simulation$params.org.eval
  estR <- res$estimation$result$params.estimated
  ## sign-flip factor eigenvectors (cols 1 & 2)
  n.eig.fac <- length(est$sigma.fac.eigvec.std)
  for (k in 1:n.eig.fac){
    d <- dim(est$sigma.fac.eigvec.std[[k]])
    for (cc in 1:2){
      s <- mean(est$sigma.fac.eigvec.std[[k]][,cc]*org$sigma.fac.eigvec.std[[k]][,cc])
      if(s<0) est$sigma.fac.eigvec.std[[k]][,cc] <- -est$sigma.fac.eigvec.std[[k]][,cc]
    }
  }
  out <- vector("list", n.coef.fac)
  for (j in 1:n.coef.fac){
    ## factor-loading MSE (single factor per indicator)
    namm <- names(est$coef.fac.std[[j]])
    # find effect for this indicator/factor
    eff <- NA
    for (kk in seq_along(res$model.fit$mod$factorModel)){
      fm <- res$model.fit$mod$factorModel[[kk]]
      if(fm$indicator==nam.fac[j]) eff <- fm$effect[1]
    }
    if(identical(eff,"historical")){
      lam <- 1/(200^2)*sum((est$coef.fac.std[[j]][[1]]-org$coef.fac.std[[j]][[1]])^2)
    } else {
      lam <- 1/200*sum((est$coef.fac.std[[j]][[1]]-org$coef.fac.std[[j]][[1]])^2)
    }
    beta <- 1/200*sum((est$intercept.std[[j]]-org$intercept.std[[j]])^2)
    phi1 <- 1/200*sum((est$sigma.fac.eigvec.std[[j]][,1]-org$sigma.fac.eigvec.std[[j]][,1])^2)
    nu1  <- (est$sigma.fac.eigval.std[[j]][1]-org$sigma.fac.eigval.std[[j]][1])^2
    sigE <- (estR$totparam$sigma.error[[j]]-org$sigma.error[[j]])^2
    out[[j]] <- c(beta=beta, lambda=lam, phi1=phi1, nu1=nu1, sigma=sigE)
  }
  out
}

getNM <- function(res){
  N <- length(unique(res$simulation$data$.id))
  M <- max(sapply(1:N,function(j){max(sapply(1:length(res$model.fit$mod$factorModel),
        function(jj){length(res$simulation$data$.t[res$simulation$data$.id==j & res$simulation$data$.ind==jj])}))}))
  c(N=N,M=M)
}

files <- list.files("outputs", pattern="^result_parameter_set1_(regular|irregular)_[0-9]+_N_[0-9]+_M_[0-9]+\\.rds$", full.names=TRUE)
cat("total files:", length(files), "\n")
t0 <- Sys.time()
for (fi in seq_along(files)){
  f <- files[fi]
  design <- if(grepl("_irregular_", f)) "irregular" else "regular"
  res <- readRDS(f)
  nm <- getNM(res)
  key <- keyf(design, nm["N"], nm["M"])
  m <- per_file_mse(res)
  if(is.null(acc[[key]])){
    acc[[key]] <- list(sum=lapply(m, function(x) x*0), count=0)
  }
  for (j in seq_along(m)) acc[[key]]$sum[[j]] <- acc[[key]]$sum[[j]] + m[[j]]
  acc[[key]]$count <- acc[[key]]$count + 1L
  rm(res); 
  if(fi %% 50 == 0){ gc(); cat(sprintf("  %d/%d  (%.1f min)\n", fi, length(files), as.numeric(difftime(Sys.time(),t0,units="mins")))) }
}

## emit results
cat("\n================ TABLE 1 MSE RESULTS ================\n")
for (design in designs){
  cat("\n#### design:", design, "####\n")
  for (r in 1:nrow(combs)){
    N <- combs$N[r]; M <- combs$M[r]
    key <- keyf(design,N,M)
    if(is.null(acc[[key]])){ cat(sprintf("N=%d M=%d : MISSING\n",N,M)); next }
    cnt <- acc[[key]]$count
    for (j in 1:3){
      mse <- acc[[key]]$sum[[j]]/cnt
      cat(sprintf("FM(%d) N=%3d M=%2d (n=%3d): beta=%.4f lambda=%.4f phi=%.4f nu=%.4f sigma=%.5f\n",
                  j,N,M,cnt, mse["beta"],mse["lambda"],mse["phi1"],mse["nu1"],mse["sigma"]))
    }
  }
}
saveRDS(acc, "_mse_table1_acc.rds")
cat("\nDONE. total time:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
