setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
options(width=200)
cat("fregion installed:", requireNamespace("fregion", quietly=TRUE), "\n")
suppressMessages({library(Matrix); library(fda)})

per_file_cov_fac <- function(re, boot=200, no.it=10){
  mont <- 100
  ind <- seq(1,200,by=2)
  no.f <- length(re$model.fit$mod$factorModel)
  lat2 <- re$model.fit$var$latents
  n.b <- length(re$estimation$result$params.estimated$intercept[[1]])
  basis <- create.bspline.basis(nbasis=n.b, rangeval=c(0,1))
  times <- (seq(0,1,length.out=200))[ind]
  eval1 <- t(eval.basis(times, basis))
  samp <- unique(re$simulation$data$.id)
  lambda1 <- vector("list", boot)
  for (k in 1:boot){
    lambda1[[k]] <- vector("list", no.f)
    samples <- sort(sample(samp, replace=TRUE))
    for (j in 1:no.f){
      count3 <- which(lat2==re$model.fit$mod$factorModel[[j]]$factor)
      ff<-vector("list",mont); zz<-vector("list",mont); ff2<-vector("list",mont)
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
      d.sem2 <- t(eval1)%*%re$estimation$result$params.estimated$sigma.fac[[count3]]%*%eval1  # placeholder; replaced below
      # NB sig.sem for factor count uses sem sigma; main sim factor==latent so use sem
      lambda1[[k]][[j]] <- lambda2.1
    }
  }
  lambda1
}

f <- "outputs/result_parameter_set1_regular_1_N_50_M_10.rds"
re <- readRDS(f)
t0 <- Sys.time()
tmp <- per_file_cov_fac(re, boot=5, no.it=10)
dt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat(sprintf("5 bootstrap reps took %.1f sec -> full boot=200 ~ %.1f sec/run -> 800 runs ~ %.1f hours\n",
    dt, dt/5*200, dt/5*200*800/3600))
