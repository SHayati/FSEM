#################extracting norm values##################
library(Matrix)
library(tidyverse)
if (requireNamespace("writexl", quietly = TRUE)) library(writexl)

## The original one-factor procedural computation below operates on an in-memory
## `results` object and the main-simulation design (N=100, M=8).  Set
## FSEM_MSE_DEFINE_ONLY <- TRUE before sourcing this file to SKIP it and only load
## the generalised functions defined at the bottom (used by the two new
## two-factor simulations, whose `results` have a different structure/design).
if (!exists("FSEM_MSE_DEFINE_ONLY") || !isTRUE(FSEM_MSE_DEFINE_ONLY)) {

dat_list <-results
n<-length(dat_list)
n.coef.fac<-length(dat_list[[1]]$model.fit$var$indicators)
n.coef.sem<-length(dat_list[[1]]$model.fit$var$factors)
nam.fac<-dat_list[[1]]$model.fit$var$indicators
nam.sem<-dat_list[[1]]$model.fit$var$latents
n.sigma.fac.vec<-n.coef.fac
n.sigma.fac.val<-n.coef.fac
n.sigma.sem.vec<-n.coef.sem
n.sigma.sem.val<-n.coef.sem
n.intercept<-n.coef.fac
n.sigma.err<-n.coef.fac
no.f<-length(dat_list[[1]]$model.fit$mod$factorModel)
norm<-list()
for (i in 1:n) {
  res <- dat_list[[i]]
  N<-length(unique(dat_list[[1]]$simulation$data$.id)) 
  M<-max(sapply(1:N,function(j){max(sapply(1:length(dat_list[[1]]$model.fit$mod$factorModel), function(jj){length(dat_list[[1]]$simulation$data$.t[dat_list[[1]]$simulation$data$.id==j&dat_list[[1]]$simulation$data$.ind==jj])}))}))
  n.samp<-N
  n.t<-M
  n.eigen.fac<-length(res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std)
  for (k in 1:n.eigen.fac) {
    d<-dim(res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]])
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][kk,1]*
        res$simulation$params.org.eval$sigma.fac.eigvec.std[[k]][kk,1]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,1]<--1*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,1]
    }
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][kk,2]*
        res$simulation$params.org.eval$sigma.fac.eigvec.std[[k]][kk,2]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,2]<--1*res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[k]][,2]
    }
  }
  
  n.eigen.sem<-length(res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std)
  for (k in 1:n.eigen.sem) {
    d<-dim(res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]])
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][kk,1]*
        res$simulation$params.org.eval$sigma.sem.eigvec.std[[k]][kk,1]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,1]<--1*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,1]
    }
    s<-sum(sapply(1:d[1], function(kk){1/d[1]*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][kk,2]*
        res$simulation$params.org.eval$sigma.sem.eigvec.std[[k]][kk,2]}))
    if(s<0){
      res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,2]<--1*res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[k]][,2]
    }
  }
  
  norm[[paste0("res",i)]]<-list(N=N,M=M)
  for (j in 1:n.coef.fac) {
    namm<-names(res$estimation$result$params.estimated.eval$coef.fac.std[[j]])
    if(length(res$estimation$result$params.estimated.eval$coef.fac.std[[j]])>0){
      m<-length(namm)
      for (k in 1:m) {
        for (kk in 1:length(res$model.fit$mod$factorModel)) {
          if(res$model.fit$mod$factorModel[[kk]]$indicator==nam.fac[j]){
            nn<-kk
          } 
          for (kk2 in 1:length(res$model.fit$mod$factorModel[[k]]$factor)) {
            if(res$model.fit$mod$factorModel[[kk]]$factor[kk2]==namm[k]){
              mm<-kk2
            } 
          }
        }
        if(res$model.fit$mod$factorModel[[nn]]$effect[mm]=="historical"){
          norm[[paste0("res",i)]]$coef.fac[[nam.fac[j]]][[namm[k]]]<-1/(200^2)*sum((res$estimation$result$params.estimated.eval$coef.fac.std[[j]][[k]]-
                                                                                  res$simulation$params.org.eval$coef.fac.std[[j]][[k]])^2)
        }else{
          norm[[paste0("res",i)]]$coef.fac[[nam.fac[j]]][[namm[k]]]<-1/200*sum((res$estimation$result$params.estimated.eval$coef.fac.std[[j]][[k]]-
                                                                                  res$simulation$params.org.eval$coef.fac.std[[j]][[k]])^2)
        }
      }
    }  
  }
  for (j in 1:n.coef.sem) {
    namm<-names(res$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
    if(length(namm)>0){
      m<-length(namm)
      for (k in 1:m) {
        for (kk in 1:length(res$model.fit$mod$regression)) {
          if(res$model.fit$mod$regression[[kk]]$response==nam.sem[j]){
            nn<-kk
          } }
        for (kk2 in 1:length(res$model.fit$mod$regression[[nn]]$covariate)) {
          if(res$model.fit$mod$regression[[kk]]$covariate[kk2]==namm[k]){
            mm<-kk2
          } 
        }
        
        if(res$model.fit$mod$regression[[nn]]$effect[mm]=="historical"){
          norm[[paste0("res",i)]]$coef.sem[[nam.sem[j]]][[namm[k]]]<-1/(200^2)*sum((res$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]][[k]]-
                                                                                  res$simulation$params.org.eval$coef.sem.std[[nam.sem[j]]][[k]])^2)
        }else{
          norm[[paste0("res",i)]]$coef.sem[[nam.sem[j]]][[namm[k]]]<-1/200*sum((res$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]][[k]]-
                                                                                  res$simulation$params.org.eval$coef.sem.std[[nam.sem[j]]][[k]])^2)
        }
      }
    }
  }
  for (j in 1:n.intercept) {
    norm[[paste0("res",i)]]$intercept[[nam.fac[j]]]<-1/200*sum((res$estimation$result$params.estimated.eval$intercept.std[[j]]-
                                                                  res$simulation$params.org.eval$intercept.std[[j]])^2)
  }
  for (j in 1:n.sigma.fac.vec) {
    d<-dim(res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[nam.fac[j]]])[2]
    norm[[paste0("res",i)]]$sigma.fac.eigvec[[nam.fac[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.fac.eigvec[[nam.fac[j]]][k]<-1/200*sum((res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[j]][,k]-
                                                                              res$simulation$params.org.eval$sigma.fac.eigvec.std[[j]][,k])^2)  
    }
  }
  for (j in 1:n.sigma.sem.vec) {
    norm[[paste0("res",i)]]$sigma.sem.eigvec[[nam.sem[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.sem.eigvec[[nam.sem[j]]][k]<-1/200*sum((res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[j]][,k]-
                                                                              res$simulation$params.org.eval$sigma.sem.eigvec.std[[j]][,k])^2)  
    } 
  }
  for (j in 1:n.sigma.fac.val) {
    norm[[paste0("res",i)]]$sigma.fac.eigval[[nam.fac[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.fac.eigval[[nam.fac[j]]][k]<-(res$estimation$result$params.estimated.eval$sigma.fac.eigval.std[[j]][k]-
                                                                    res$simulation$params.org.eval$sigma.fac.eigval.std[[j]][k])^2  
    }
  }
  for (j in 1:n.sigma.sem.val) {
    norm[[paste0("res",i)]]$sigma.sem.eigval[[nam.sem[j]]]<-c()
    for (k in 1:d) {
      norm[[paste0("res",i)]]$sigma.sem.eigval[[nam.sem[j]]][k]<-(res$estimation$result$params.estimated.eval$sigma.sem.eigval.std[[j]][k]-
                                                                    res$simulation$params.org.eval$sigma.sem.eigval.std[[j]][k])^2  
    } 
  }
  for (j in 1:n.sigma.err) {
    norm[[paste0("res",i)]]$sigma.error[[nam.fac[j]]]<-(res$estimation$result$params.estimated$totparam$sigma.error[[j]]-
                                                          res$simulation$params.org.eval$sigma.error[[j]])^2
  }
}



comb<-expand.grid(N=100,M=8)
comb<-comb[order(comb$N),]
mse<-list()
ress<-dat_list[[1]]
for (l in 1:dim(comb)[1]) {
  mse[[l]]<-list(N=comb$N[l],M=comb$M[l])
  for (j in 1:n.coef.fac) {
    namm<-names(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])
    mse[[l]][[nam.fac[j]]]<-data.frame(NA)
    if(length(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])>0){
      m<-length(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])
      for (k in 1:m) {
        count<-0
        cc<-c()
        for (i in 1:n) {
          if(N==comb[l,]$N&M==comb[l,]$M){
            count<-count+1
            cc[count]<-norm[[paste0("res",i)]]$coef.fac[[nam.fac[j]]][[namm[k]]]
          }
        }
        mse[[l]][[nam.fac[j]]][[namm[k]]]<-mean(cc,trim=0.1)
      }
    }
    count<-0
    for (i in 1:n) {
      if(N==comb[l,]$N&M==comb[l,]$M){
        count<-count+1
        cc[count]<-norm[[paste0("res",i)]]$intercept[[nam.fac[j]]]
      }
    }
    mse[[l]][[nam.fac[j]]]$intercept<-mean(cc,trim=0.1)
    for (k in 1:2) {
      count<-0
      cc1<-c()
      cc2<-c()
      for (i in 1:n) {
        if(N==comb[l,]$N&M==comb[l,]$M){
          count<-count+1
          cc1[count]<-norm[[paste0("res",i)]]$sigma.fac.eigvec[[nam.fac[j]]][k]
          cc2[count]<-norm[[paste0("res",i)]]$sigma.fac.eigval[[nam.fac[j]]][k]
        }
      }
      mse[[l]][[nam.fac[j]]][[paste0("eig.fuc.",k)]]<-mean(cc1,trim=0.1)
      mse[[l]][[nam.fac[j]]][[paste0("eig.val.",k)]]<-mean(cc2,trim=0.1)
    }
    count<-0
    cc<-c()
    for (i in 1:n) {
      if(N==comb[l,]$N&M==comb[l,]$M){
        count<-count+1
        cc[count]<-norm[[paste0("res",i)]]$sigma.err[[nam.fac[j]]]
      }
    }
    mse[[l]][[nam.fac[j]]]$sig.err<-mean(cc,trim=0.1)
  } 
  
  for (j in 1:n.coef.sem) {
    namm<-names(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
    mse[[l]][[nam.sem[j]]]<-data.frame(NA)
    if(length(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])>0){
      m<-length(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
      for (k in 1:m) {
        count<-0
        cc<-c()
        for (i in 1:n) {
          if(N==comb[l,]$N&M==comb[l,]$M){
            count<-count+1
            cc[count]<-norm[[paste0("res",i)]]$coef.sem[[nam.sem[j]]][[namm[k]]]
          }
        }
        mse[[l]][[nam.sem[j]]][[namm[k]]]<-mean(cc,trim=0.1)
      }
    }
    for (k in 1:2) {
      count<-0
      cc1<-c()
      cc2<-c()
      for (i in 1:n) {
        if(N==comb[l,]$N&M==comb[l,]$M){
          count<-count+1
          cc1[count]<-norm[[paste0("res",i)]]$sigma.sem.eigvec[[nam.sem[j]]][k]
          cc2[count]<-norm[[paste0("res",i)]]$sigma.sem.eigval[[nam.sem[j]]][k]
        }
      }
      mse[[l]][[nam.sem[j]]][[paste0("eig.fuc.",k)]]<-mean(cc1,trim=0.1)
      mse[[l]][[nam.sem[j]]][[paste0("eig.val.",k)]]<-mean(cc2,trim=0.1)
    }
  } 
  mse[[l]]$numb.runs<-count
}

#################extracting tables###############################
tabl.fac<-data.frame(N=comb$N,M=comb$M)
for (l in 1:dim(comb)[1]) {
  for (j in 1:n.coef.fac) {
    tabl.fac[[paste0("beta",j)]][l]<-mse[[l]][[nam.fac[j]]]$intercept
    namm<-names(ress$estimation$result$params.estimated.eval$coef.fac.std[[j]])
    for (m in 1:length(namm)) {
      tabl.fac[[paste0("lambda",j,m)]][l]<-mse[[l]][[nam.fac[j]]][[namm[m]]]
    }
    tabl.fac[[paste0("phi",j,1)]][l]<-mse[[l]][[nam.fac[j]]]$eig.fuc.1
    tabl.fac[[paste0("phi",j,2)]][l]<-mse[[l]][[nam.fac[j]]]$eig.fuc.2
    tabl.fac[[paste0("nu",j,1)]][l]<-mse[[l]][[nam.fac[j]]]$eig.val.1
    tabl.fac[[paste0("nu",j,2)]][l]<-mse[[l]][[nam.fac[j]]]$eig.val.2
    tabl.fac[[paste0("sigma",j)]][l]<-mse[[l]][[nam.fac[j]]]$sig.err
  }  
}


tabl.sem<-data.frame(N=comb$N,M=comb$M)
for (l in 1:dim(comb)[1]) {
  for (j in 1:n.coef.sem) {
    namm<-names(ress$estimation$result$params.estimated.eval$coef.sem.std[[nam.sem[j]]])
    if(length(namm)>0){
      for (m in 1:length(mse[[l]][[nam.sem[j]]])) {
        tabl.sem[[paste0("gamma",j,m)]][l]<-mse[[l]][[nam.sem[j]]][[namm[m]]]
      } 
    }
  }
  for (j in 1:n.coef.sem) {
    tabl.sem[[paste0("psi",j,1)]][l]<-mse[[l]][[nam.sem[j]]]$eig.fuc.1
    tabl.sem[[paste0("psi",j,2)]][l]<-mse[[l]][[nam.sem[j]]]$eig.fuc.2
  }
  for (j in 1:n.coef.sem) {
    tabl.sem[[paste0("mu",j,1)]][l]<-mse[[l]][[nam.sem[j]]]$eig.val.1
    tabl.sem[[paste0("mu",j,2)]][l]<-mse[[l]][[nam.sem[j]]]$eig.val.2    
  }
} 

}  # end  if(!FSEM_MSE_DEFINE_ONLY)  -- original one-factor procedural block


## ===========================================================================
## EXTENSION  --  MSE for the two NEW two-factor simulation studies
## ---------------------------------------------------------------------------
## The block above reproduces the one-factor table of the main simulation.
## The functions below GENERALISE that computation to models with more than one
## latent factor and with latent-on-latent (structural) regressions, so the same
## squared-error summaries can be produced for:
##     (A) two correlated indicators per factor, UNCORRELATED factors  -> outputs/
##     (B) DEPENDENT factors, well-specified (.w) / mis-specified (.m)  -> outputs3/
##
## NB. No simulation is re-run here.  Each stored result object is 300-500 MB, so
## the folder driver streams the .rds files ONE AT A TIME (read -> summarise ->
## discard -> gc()).  The two main_simulation_*_parallel.R scripts call the
## in-memory entry point on the `results` list they have just produced.
## ===========================================================================

## sign-align estimated eigenvectors to the truth (column-wise)
.mse_align_eig <- function(est, truth) {
  if (length(est) == 0 || length(truth) == 0) return(est)
  for (i in seq_along(est)) {
    em <- est[[i]]; tm <- truth[[i]]
    if (is.null(dim(em)) || is.null(dim(tm))) next
    for (cc in seq_len(min(ncol(em), ncol(tm)))) {
      ip <- mean(em[, cc] * tm[, cc])
      if (!is.na(ip) && ip < 0) est[[i]][, cc] <- -est[[i]][, cc]
    }
  }
  est
}

## look up the effect type ("concurrent"/"historical"/...) of a coefficient
.mse_effect_type <- function(model_list, resp, cov, rf, cf) {
  for (e in model_list) {
    if (!identical(e[[rf]], resp)) next
    mi <- match(cov, e[[cf]])
    if (!is.na(mi)) return(e$effect[[mi]])
  }
  NA_character_
}

## per-replicate squared-error "entry" for one result object
.mse_one_entry <- function(res, fit_key, est_key) {
  mfit <- res[[fit_key]]
  est  <- res[[est_key]]$result$params.estimated.eval
  tru  <- res$simulation$params.org.eval
  nam.fac <- mfit$var$indicators
  nam.sem <- mfit$var$latents

  est$sigma.fac.eigvec.std <- .mse_align_eig(est$sigma.fac.eigvec.std, tru$sigma.fac.eigvec.std)
  est$sigma.sem.eigvec.std <- .mse_align_eig(est$sigma.sem.eigvec.std, tru$sigma.sem.eigvec.std)

  N <- length(unique(res$simulation$data$.id))
  M <- max(sapply(unique(res$simulation$data$.id), function(id)
        max(sapply(seq_along(mfit$mod$factorModel), function(ii)
          length(res$simulation$data$.t[res$simulation$data$.id == id & res$simulation$data$.ind == ii])))))
  e <- list(N = N, M = M, coef.fac = list(), intercept = list(),
            sigma.fac.eigvec = list(), sigma.fac.eigval = list(), sigma.error = list(),
            coef.sem = list(), sigma.sem.eigvec = list(), sigma.sem.eigval = list())

  for (fi in seq_along(nam.fac)) {
    ind <- nam.fac[[fi]]
    fn <- names(est$coef.fac.std[[fi]])
    e$coef.fac[[ind]] <- list()
    if (length(fn) > 0) for (k in seq_along(fn)) {
      eff <- .mse_effect_type(mfit$mod$factorModel, ind, fn[[k]], "indicator", "factor")
      den <- if (identical(eff, "historical")) 200^2 else 200
      e$coef.fac[[ind]][[fn[[k]]]] <- sum((est$coef.fac.std[[fi]][[k]] - tru$coef.fac.std[[fi]][[k]])^2) / den
    }
    e$intercept[[ind]] <- mean((est$intercept.std[[fi]] - tru$intercept.std[[fi]])^2)
    nc <- ncol(est$sigma.fac.eigvec.std[[fi]])
    e$sigma.fac.eigvec[[ind]] <- vapply(seq_len(nc), function(cc)
      mean((est$sigma.fac.eigvec.std[[fi]][, cc] - tru$sigma.fac.eigvec.std[[fi]][, cc])^2), numeric(1))
    e$sigma.fac.eigval[[ind]] <- vapply(seq_along(est$sigma.fac.eigval.std[[fi]]), function(cc)
      (est$sigma.fac.eigval.std[[fi]][[cc]] - tru$sigma.fac.eigval.std[[fi]][[cc]])^2, numeric(1))
    e$sigma.error[[ind]] <- (res[[est_key]]$result$params.estimated$totparam$sigma.error[[fi]] -
                             tru$sigma.error[[fi]])^2
  }
  for (si in seq_along(nam.sem)) {
    lat <- nam.sem[[si]]
    sn <- names(est$coef.sem.std[[lat]])
    e$coef.sem[[lat]] <- list()
    if (length(sn) > 0) for (k in seq_along(sn)) {
      eff <- .mse_effect_type(mfit$mod$regression, lat, sn[[k]], "response", "covariate")
      den <- if (identical(eff, "historical")) 200^2 else 200
      e$coef.sem[[lat]][[sn[[k]]]] <- sum((est$coef.sem.std[[lat]][[k]] - tru$coef.sem.std[[lat]][[k]])^2) / den
    }
    nc <- ncol(est$sigma.sem.eigvec.std[[si]])
    e$sigma.sem.eigvec[[lat]] <- vapply(seq_len(nc), function(cc)
      mean((est$sigma.sem.eigvec.std[[si]][, cc] - tru$sigma.sem.eigvec.std[[si]][, cc])^2), numeric(1))
    e$sigma.sem.eigval[[lat]] <- vapply(seq_along(est$sigma.sem.eigval.std[[si]]), function(cc)
      (est$sigma.sem.eigval.std[[si]][[cc]] - tru$sigma.sem.eigval.std[[si]][[cc]])^2, numeric(1))
  }
  e
}

## average a list of entries into tabl.fac / tabl.sem (mean over replicates per N,M)
.mse_aggregate <- function(entries, template_fit, template_eval) {
  nam.fac <- template_fit$var$indicators
  nam.sem <- template_fit$var$latents
  combos <- unique(do.call(rbind, lapply(entries, function(e) data.frame(N = e$N, M = e$M))))
  combos <- combos[order(combos$N, combos$M), , drop = FALSE]; rownames(combos) <- NULL

  tabl.fac <- data.frame(N = combos$N, M = combos$M)
  tabl.sem <- data.frame(N = combos$N, M = combos$M)
  for (r in seq_len(nrow(combos))) {
    mt <- Filter(function(e) e$N == combos$N[r] && e$M == combos$M[r], entries)
    for (fi in seq_along(nam.fac)) {
      ind <- nam.fac[[fi]]
      tabl.fac[[paste0("beta", fi)]][r] <- mean(vapply(mt, function(e) e$intercept[[ind]], numeric(1)))
      fn <- names(template_eval$coef.fac.std[[fi]])
      for (k in seq_along(fn))
        tabl.fac[[paste0("lambda", fi, k)]][r] <- mean(vapply(mt, function(e) e$coef.fac[[ind]][[fn[[k]]]], numeric(1)))
      nc <- ncol(template_eval$sigma.fac.eigvec.std[[fi]])
      for (cc in seq_len(nc)) {
        tabl.fac[[paste0("phi", fi, cc)]][r] <- mean(vapply(mt, function(e) e$sigma.fac.eigvec[[ind]][[cc]], numeric(1)))
        tabl.fac[[paste0("nu", fi, cc)]][r]  <- mean(vapply(mt, function(e) e$sigma.fac.eigval[[ind]][[cc]], numeric(1)))
      }
      tabl.fac[[paste0("sigma", fi)]][r] <- mean(vapply(mt, function(e) e$sigma.error[[ind]], numeric(1)))
    }
    for (si in seq_along(nam.sem)) {
      lat <- nam.sem[[si]]
      sn <- names(template_eval$coef.sem.std[[lat]])
      for (k in seq_along(sn))
        tabl.sem[[paste0("gamma", si, k)]][r] <- mean(vapply(mt, function(e) e$coef.sem[[lat]][[sn[[k]]]], numeric(1)))
      nc <- ncol(template_eval$sigma.sem.eigvec.std[[si]])
      for (cc in seq_len(nc)) {
        tabl.sem[[paste0("psi", si, cc)]][r] <- mean(vapply(mt, function(e) e$sigma.sem.eigvec[[lat]][[cc]], numeric(1)))
        tabl.sem[[paste0("mu", si, cc)]][r]  <- mean(vapply(mt, function(e) e$sigma.sem.eigval[[lat]][[cc]], numeric(1)))
      }
    }
  }
  list(tabl.fac = tabl.fac, tabl.sem = tabl.sem)
}

## ---- public entry points --------------------------------------------------
## In-memory: summarise a `results` list already held in memory (used by the
## main_simulation_*_parallel.R scripts).  `fit_model` may supply a fitted model
## that was not stored inside each element (as for the dependent simulation).
fsem_mse_from_list <- function(results, fit_key = "model.fit", est_key = "estimation",
                               fit_model = NULL) {
  entries <- vector("list", length(results))
  tf <- NULL; te <- NULL
  for (i in seq_along(results)) {
    res <- results[[i]]
    if (!is.null(fit_model)) res[[fit_key]] <- fit_model
    if (is.null(tf)) { tf <- res[[fit_key]]; te <- res[[est_key]]$result$params.estimated.eval }
    entries[[i]] <- .mse_one_entry(res, fit_key, est_key)
  }
  .mse_aggregate(entries, tf, te)
}

## Streaming: summarise every result .rds in `folder`, one file at a time.
fsem_mse_from_folder <- function(folder, fit_key = "model.fit", est_key = "estimation",
                                 fit_model = NULL, verbose = TRUE) {
  files <- list.files(folder, pattern = "^result_\\d+_N_\\d+_M_\\d+\\.rds$", full.names = TRUE)
  entries <- vector("list", length(files))
  tf <- NULL; te <- NULL
  for (i in seq_along(files)) {
    res <- readRDS(files[i])
    if (!is.null(fit_model)) res[[fit_key]] <- fit_model
    if (is.null(tf)) { tf <- res[[fit_key]]; te <- res[[est_key]]$result$params.estimated.eval }
    entries[[i]] <- .mse_one_entry(res, fit_key, est_key)
    rm(res); gc(verbose = FALSE)
    if (verbose && i %% 25 == 0) cat("  ", i, "/", length(files), "files\n")
  }
  .mse_aggregate(entries, tf, te)
}

## ---- driver: regenerate the new-simulation MSE CSVs from stored results -----
## Set RUN_NEW_SIM_MSE <- TRUE before sourcing this file to (re)build the tables
## by streaming the outputs/ and outputs3/ folders.  Off by default so that
## sourcing this script for the one-factor table above stays cheap.
if (exists("RUN_NEW_SIM_MSE") && isTRUE(RUN_NEW_SIM_MSE)) {
  suppressMessages({ library(fda); library(Matrix) })
  source("Core/FSEM.R"); source("Core/simualtion_utilities.R")

  ## fitted models for the dependent simulation are NOT stored in outputs3/
  model.fit.m <- fsem(eta1~~z1, effectType="fixed_concurrent") %+%
    fsem(eta1~~z2, effectType="concurrent")                   %+%
    fsem(eta1~~z3, effectType="historical")                   %+%
    fsem(eta2~~z4, effectType="fixed_concurrent")             %+%
    fsem(eta2~~z5+z6, effectType="concurrent")                %+%
    fsem(eta1~-1, effectType="concurrent")                    %+%
    fsem(eta2~-1)
  model.fit.w <- fsem(eta1~~z1, effectType="fixed_concurrent") %+%
    fsem(eta1~~z2, effectType="concurrent")                   %+%
    fsem(eta1~~z3, effectType="historical")                   %+%
    fsem(eta2~~z4, effectType="fixed_concurrent")             %+%
    fsem(eta2~~z5+z6, effectType="concurrent")                %+%
    fsem(eta1~-1+eta2, effectType="concurrent", latent.covariate="eta2", scalar.covariate=FALSE) %+%
    fsem(eta2~-1)

  cat("=== (A) Uncorrelated two-factor: streaming MSE ===\n")
  a <- fsem_mse_from_folder("outputs", "model.fit", "estimation")
  write.csv(a$tabl.fac, "outputs/MSE_uncor_fac.csv", row.names = FALSE)
  write.csv(a$tabl.sem, "outputs/MSE_uncor_sem.csv", row.names = FALSE)

  cat("=== (B) Dependent two-factor (well-specified): streaming MSE ===\n")
  w <- fsem_mse_from_folder("outputs3", "model.fit.w", "estimation.w", fit_model = model.fit.w)
  write.csv(w$tabl.fac, "outputs3/MSE_dep_well_fac.csv", row.names = FALSE)
  write.csv(w$tabl.sem, "outputs3/MSE_dep_well_sem.csv", row.names = FALSE)

  cat("=== (B) Dependent two-factor (mis-specified): streaming MSE ===\n")
  m <- fsem_mse_from_folder("outputs3", "model.fit.m", "estimation.m", fit_model = model.fit.m)
  write.csv(m$tabl.fac, "outputs3/MSE_dep_miss_fac.csv", row.names = FALSE)
  write.csv(m$tabl.sem, "outputs3/MSE_dep_miss_sem.csv", row.names = FALSE)
  cat("New-simulation MSE tables written.\n")
}
