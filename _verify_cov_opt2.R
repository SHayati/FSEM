setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
options(width=200)
suppressMessages({library(Matrix); library(fda)})
f <- "outputs/result_parameter_set1_regular_1_N_100_M_20.rds"  # slow case
re <- readRDS(f)
mont <- 100; no.f <- length(re$model.fit$mod$factorModel); no.it<-10
samp <- unique(re$simulation$data$.id)

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
    deltaa <- re$estimation$hist$hist.fac[[j]]$deltaa; pen.fac <- re$estimation$hist$hist.fac[[j]]$pen.fac
    for (it in 1:no.it){
      sol <- solve(kronecker(diag(samples), sig)); lam.sum<-0; lam1.sum<-0
      for (m in 1:mont){
        lam <- t(ff[[m]])%*%sol%*%ff[[m]] + deltaa*pen.fac
        lam1 <- t(ff[[m]])%*%sol%*%(as.vector(zz[[m]])-ff2[[m]])
        lam.sum<-lam.sum+1/mont*lam; lam1.sum<-lam1.sum+1/mont*lam1
      }
      lambda2.1 <- solve(lam.sum)%*%lam1.sum; sigma1<-0
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

opt_one <- function(re, samples){
  out <- vector("list", no.f); ns <- length(samples); ids <- paste0(samples); inv.samp <- 1/samples
  for (j in 1:no.f){
    hf <- re$estimation$hist$hist.fac[[j]]
    f0 <- hf$f[[1]][[ids[1]]]; kdim <- nrow(f0); pdim <- ncol(f0)
    ff<-vector("list",mont); zz<-vector("list",mont); ff2<-vector("list",mont)
    for (m in 1:mont){
      ffm <- matrix(0, ns*kdim, pdim); zzm <- numeric(ns*kdim); ff2m <- numeric(ns*kdim); r0<-0
      for (s in 1:ns){
        id<-ids[s]; idx<-(r0+1):(r0+kdim)
        ffm[idx,]<-hf$f[[m]][[id]]; zzm[idx]<-as.vector(hf$generators2.z[[id]][[j]][[m]]); ff2m[idx]<-hf$f2[[m]][[id]]
        r0<-r0+kdim
      }
      ff[[m]]<-ffm; zz[[m]]<-zzm; ff2[[m]]<-ff2m
    }
    sig <- re$estimation$result$params.estimated$sigma.fac[[j]]; deltaa<-hf$deltaa; pen.fac<-hf$pen.fac
    for (it in 1:no.it){
      siginv <- solve(sig); lam.sum<-0; lam1.sum<-0
      for (m in 1:mont){
        ffm<-ff[[m]]; resid<-zz[[m]]-ff2[[m]]
        WA <- matrix(0, ns*kdim, pdim)
        for (c in 1:pdim) WA[,c] <- as.vector(sweep(siginv %*% matrix(ffm[,c],kdim,ns),2,inv.samp,'*'))
        Wb <- as.vector(sweep(siginv %*% matrix(resid,kdim,ns),2,inv.samp,'*'))
        lam.sum <- lam.sum + (crossprod(ffm,WA)+deltaa*pen.fac)/mont
        lam1.sum <- lam1.sum + crossprod(ffm,Wb)/mont
      }
      lambda2.1 <- solve(lam.sum, lam1.sum); sigma1<-matrix(0,kdim,kdim)
      for (m in 1:mont){
        resid <- zz[[m]]-ff2[[m]] - ff[[m]]%*%lambda2.1
        Rm <- matrix(resid, kdim, ns); sigma1 <- sigma1 + tcrossprod(Rm)/mont
      }
      sig <- sigma1/ns
    }
    out[[j]] <- lambda2.1
  }
  out
}

set.seed(7); s1 <- sort(sample(samp, replace=TRUE))
t0<-Sys.time(); o1<-orig_one(re,s1); to<-as.numeric(difftime(Sys.time(),t0,units="secs"))
t0<-Sys.time(); o2<-opt_one(re,s1); tp<-as.numeric(difftime(Sys.time(),t0,units="secs"))
for(j in 1:no.f) cat(sprintf("factor %d max|diff|=%.3e\n", j, max(abs(o1[[j]]-o2[[j]]))))
cat(sprintf("N=100,M=20: orig %.2fs  opt %.3fs  speedup %.1fx  -> opt full boot200 = %.0fs/run\n", to,tp,to/tp,tp*200))
