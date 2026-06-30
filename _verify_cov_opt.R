setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
options(width=200)
suppressMessages({library(Matrix); library(fda)})

f <- "outputs/result_parameter_set1_regular_1_N_50_M_10.rds"
re <- readRDS(f)
mont <- 100
no.f <- length(re$model.fit$mod$factorModel)
lat2 <- re$model.fit$var$latents
samp <- unique(re$simulation$data$.id)
no.it <- 10

## ---- ORIGINAL inner factor computation (verbatim core) ----
orig_one <- function(re, samples){
  out <- vector("list", no.f)
  for (j in 1:no.f){
    ff<-list(); zz<-list(); ff2<-list()
    for (m in 1:mont){
      ff[[m]]<-0; zz[[m]]<-0; ff2[[m]]<-0
      for (i1 in samples){
        ff[[m]]<-rbind(ff[[m]], re$estimation$hist$hist.fac[[j]]$f[[m]][[paste0(i1)]])
        zz[[m]]<-c(zz[[m]], re$estimation$hist$hist.fac[[j]]$generators2.z[[paste0(i1)]][[j]][[m]])
        ff2[[m]]<-c(ff2[[m]], re$estimation$hist$hist.fac[[j]]$f2[[m]][[paste0(i1)]])
      }
      ff[[m]]<-ff[[m]][-1,]; zz[[m]]<-zz[[m]][-1]; ff2[[m]]<-ff2[[m]][-1]
    }
    sig <- re$estimation$result$params.estimated$sigma.fac[[j]]
    deltaa <- re$estimation$hist$hist.fac[[j]]$deltaa
    pen.fac <- re$estimation$hist$hist.fac[[j]]$pen.fac
    for (it in 1:no.it){
      sol <- solve(kronecker(diag(samples), sig))
      lam.sum<-0; lam1.sum<-0
      for (m in 1:mont){
        lam <- t(ff[[m]])%*%sol%*%ff[[m]] + deltaa*pen.fac
        lam1 <- t(ff[[m]])%*%sol%*%(as.vector(zz[[m]])-ff2[[m]])
        lam.sum<-lam.sum+1/mont*lam; lam1.sum<-lam1.sum+1/mont*lam1
      }
      lambda2.1 <- solve(lam.sum)%*%lam1.sum
      sigma1<-0
      for (i1 in samples){
        ss2<-0
        for (m in 1:mont){
          ss <- as.vector(re$estimation$hist$hist.fac[[j]]$generators2.z[[paste0(i1)]][[j]][[m]]) -
                re$estimation$hist$hist.fac[[j]]$f2[[m]][[paste0(i1)]] -
                re$estimation$hist$hist.fac[[j]]$f[[m]][[paste0(i1)]]%*%lambda2.1
          ss2 <- ss2 + 1/mont*(ss%*%t(ss))
        }
        sigma1 <- sigma1 + 1/length(samples)*ss2
      }
      sig <- sigma1
    }
    out[[j]] <- lambda2.1
  }
  out
}

## ---- OPTIMISED: pre-extract per-id blocks once, pre-allocate, kron identity ----
opt_one <- function(re, samples){
  out <- vector("list", no.f)
  ns <- length(samples)
  ids <- paste0(samples)
  for (j in 1:no.f){
    hf <- re$estimation$hist$hist.fac[[j]]
    # block dims from first id
    f0 <- hf$f[[1]][[ids[1]]]
    kdim <- nrow(f0); pdim <- ncol(f0)
    # build ff[[m]] (ns*kdim x pdim), zz/ff2 vectors, by pre-extracting
    ff<-vector("list",mont); zz<-vector("list",mont); ff2<-vector("list",mont)
    # cache per (m,id) blocks
    for (m in 1:mont){
      rows <- ns*kdim
      ffm <- matrix(0, rows, pdim)
      zzm <- numeric(rows); ff2m <- numeric(rows)
      r0 <- 0
      for (s in 1:ns){
        id <- ids[s]
        blk <- hf$f[[m]][[id]]
        ffm[(r0+1):(r0+kdim),] <- blk
        zzm[(r0+1):(r0+kdim)] <- as.vector(hf$generators2.z[[id]][[j]][[m]])
        ff2m[(r0+1):(r0+kdim)] <- hf$f2[[m]][[id]]
        r0 <- r0 + kdim
      }
      ff[[m]]<-ffm; zz[[m]]<-zzm; ff2[[m]]<-ff2m
    }
    sig <- re$estimation$result$params.estimated$sigma.fac[[j]]
    deltaa <- hf$deltaa; pen.fac <- hf$pen.fac
    inv.samp <- 1/samples
    for (it in 1:no.it){
      siginv <- solve(sig)
      # sol = kron(diag(1/samples), siginv); t(ff)%*%sol%*%ff done block-wise
      lam.sum<-0; lam1.sum<-0
      for (m in 1:mont){
        ffm<-ff[[m]]; resid <- zz[[m]]-ff2[[m]]
        # weighted: each block s scaled by inv.samp[s]; apply siginv within block
        WA <- matrix(0, nrow(ffm), pdim)  # sol %*% ff
        Wb <- numeric(nrow(ffm))          # sol %*% resid
        r0<-0
        for (s in 1:ns){
          idx <- (r0+1):(r0+kdim)
          WA[idx,] <- inv.samp[s]*(siginv %*% ffm[idx,,drop=FALSE])
          Wb[idx]  <- inv.samp[s]*(siginv %*% resid[idx])
          r0<-r0+kdim
        }
        lam <- crossprod(ffm, WA) + deltaa*pen.fac
        lam1 <- crossprod(ffm, Wb)
        lam.sum<-lam.sum+lam/mont; lam1.sum<-lam1.sum+lam1/mont
      }
      lambda2.1 <- solve(lam.sum, lam1.sum)
      # sigma update
      sigma1 <- matrix(0,kdim,kdim)
      for (m in 1:mont){
        ffm<-ff[[m]]; resid <- zz[[m]]-ff2[[m]] - ffm%*%lambda2.1
        r0<-0
        for (s in 1:ns){
          idx <- (r0+1):(r0+kdim)
          v <- resid[idx]
          sigma1 <- sigma1 + tcrossprod(v)/mont
          r0<-r0+kdim
        }
      }
      sig <- sigma1/ns
    }
    out[[j]] <- lambda2.1
  }
  out
}

set.seed(123); samples1 <- sort(sample(samp, replace=TRUE))
t0<-Sys.time(); o1 <- orig_one(re, samples1); t_orig<-as.numeric(difftime(Sys.time(),t0,units="secs"))
t0<-Sys.time(); o2 <- opt_one(re, samples1); t_opt<-as.numeric(difftime(Sys.time(),t0,units="secs"))
for(j in 1:no.f){
  cat(sprintf("factor %d: max|orig-opt|=%.3e  (orig norm=%.4f)\n", j,
      max(abs(o1[[j]]-o2[[j]])), sqrt(sum(o1[[j]]^2))))
}
cat(sprintf("\norig: %.2f s/boot-rep   opt: %.3f s/boot-rep   speedup=%.1fx\n", t_orig, t_opt, t_orig/t_opt))
cat(sprintf("=> full boot=200: orig %.0f s/run, opt %.0f s/run\n", t_orig*200, t_opt*200))
