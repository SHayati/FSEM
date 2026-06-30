setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
options(width=200)
suppressMessages({library(Matrix); library(fda)})
stopifnot(requireNamespace("fregion", quietly=TRUE))

f <- "outputs/result_parameter_set1_regular_3_N_50_M_10.rds"
re <- readRDS(f)
mont<-100; no.it<-10; no.f<-length(re$model.fit$mod$factorModel)
lat2<-re$model.fit$var$latents
n.b<-length(re$estimation$result$params.estimated$intercept[[1]])
ind<-seq(1,200,by=2)
basis<-create.bspline.basis(nbasis=n.b,rangeval=c(0,1))
times<-(seq(0,1,length.out=200))[ind]
eval1<-t(eval.basis(times,basis))
samp<-unique(re$simulation$data$.id)
lat.count<-which(lat2==re$model.fit$mod$regression[[1]]$response)

## sig.sem (no-cov branch) verbatim
sigsem_orig <- function(re,samples){
  sigma1<-0
  for(i1 in samples){ ss2<-0
    for(m in 1:mont){ ss<-re$estimation$hist$hist.fac[[1]]$eta.l[[m]][[i1]][[lat.count]]; ss2<-ss2+1/mont*ss%*%t(ss) }
    sigma1<-sigma1+1/length(samples)*ss2 }
  sigma1
}
## factor inner verbatim (original)
lam_orig <- function(re,samples,j){
  ff<-list();zz<-list();ff2<-list()
  for(m in 1:mont){ ff[[m]]<-0;zz[[m]]<-0;ff2[[m]]<-0
    for(i1 in samples){
      ff[[m]]<-rbind(ff[[m]],re$estimation$hist$hist.fac[[j]]$f[[m]][[paste0(i1)]])
      zz[[m]]<-c(zz[[m]],re$estimation$hist$hist.fac[[j]]$generators2.z[[paste0(i1)]][[j]][[m]])
      ff2[[m]]<-c(ff2[[m]],re$estimation$hist$hist.fac[[j]]$f2[[m]][[paste0(i1)]]) }
    ff[[m]]<-ff[[m]][-1,];zz[[m]]<-zz[[m]][-1];ff2[[m]]<-ff2[[m]][-1] }
  sig<-re$estimation$result$params.estimated$sigma.fac[[j]]
  deltaa<-re$estimation$hist$hist.fac[[j]]$deltaa; pen.fac<-re$estimation$hist$hist.fac[[j]]$pen.fac
  for(it in 1:no.it){ sol<-solve(kronecker(diag(samples),sig)); lam.sum<-0;lam1.sum<-0
    for(m in 1:mont){ lam<-t(ff[[m]])%*%sol%*%ff[[m]]+deltaa*pen.fac; lam1<-t(ff[[m]])%*%sol%*%(as.vector(zz[[m]])-ff2[[m]])
      lam.sum<-lam.sum+lam/mont; lam1.sum<-lam1.sum+lam1/mont }
    lambda2.1<-solve(lam.sum)%*%lam1.sum; sigma1<-0
    for(i1 in samples){ ss2<-0
      for(m in 1:mont){ ss<-as.vector(re$estimation$hist$hist.fac[[j]]$generators2.z[[paste0(i1)]][[j]][[m]])-
        re$estimation$hist$hist.fac[[j]]$f2[[m]][[paste0(i1)]]-re$estimation$hist$hist.fac[[j]]$f[[m]][[paste0(i1)]]%*%lambda2.1
        ss2<-ss2+1/mont*(ss%*%t(ss)) }
      sigma1<-sigma1+1/length(samples)*ss2 }
    sig<-sigma1 }
  lambda2.1
}
## optimized inner
lam_opt <- function(re,samples,j){
  ns<-length(samples); ids<-paste0(samples); inv.samp<-1/samples
  hf<-re$estimation$hist$hist.fac[[j]]; f0<-hf$f[[1]][[ids[1]]]; kdim<-nrow(f0); pdim<-ncol(f0)
  ff<-vector("list",mont);zz<-vector("list",mont);ff2<-vector("list",mont)
  for(m in 1:mont){ ffm<-matrix(0,ns*kdim,pdim);zzm<-numeric(ns*kdim);ff2m<-numeric(ns*kdim);r0<-0
    for(s in 1:ns){ id<-ids[s];idx<-(r0+1):(r0+kdim)
      ffm[idx,]<-hf$f[[m]][[id]];zzm[idx]<-as.vector(hf$generators2.z[[id]][[j]][[m]]);ff2m[idx]<-hf$f2[[m]][[id]];r0<-r0+kdim }
    ff[[m]]<-ffm;zz[[m]]<-zzm;ff2[[m]]<-ff2m }
  sig<-re$estimation$result$params.estimated$sigma.fac[[j]];deltaa<-hf$deltaa;pen.fac<-hf$pen.fac
  for(it in 1:no.it){ siginv<-solve(sig);lam.sum<-0;lam1.sum<-0
    for(m in 1:mont){ ffm<-ff[[m]];resid<-zz[[m]]-ff2[[m]]
      WA<-matrix(0,ns*kdim,pdim)
      for(c in 1:pdim) WA[,c]<-as.vector(sweep(siginv%*%matrix(ffm[,c],kdim,ns),2,inv.samp,'*'))
      Wb<-as.vector(sweep(siginv%*%matrix(resid,kdim,ns),2,inv.samp,'*'))
      lam.sum<-lam.sum+(crossprod(ffm,WA)+deltaa*pen.fac)/mont; lam1.sum<-lam1.sum+crossprod(ffm,Wb)/mont }
    lambda2.1<-solve(lam.sum,lam1.sum); sigma1<-matrix(0,kdim,kdim)
    for(m in 1:mont){ resid<-zz[[m]]-ff2[[m]]-ff[[m]]%*%lambda2.1; Rm<-matrix(resid,kdim,ns); sigma1<-sigma1+tcrossprod(Rm)/mont }
    sig<-sigma1/ns }
  lambda2.1
}

transform_lambda <- function(re,j,sigsem_k,lambda2.1){
  count3<-which(lat2==re$model.sim$mod$factorModel[[j]]$factor)
  d.sem2<-t(eval1)%*%sigsem_k%*%eval1
  if(re$model.fit$mod$factorModel[[j]]$effect%in%c("fixed_concurrent","fixed_historical")){
    as.vector(diag(d.sem2)^(1/2))
  }else{
    if(re$model.sim$mod$factorModel[[j]]$effect=="concurrent"){
      ev<-diag(d.sem2)^(1/2)*t(eval.basis(times,basis))
      as.vector(t(ev)%*%lambda2.1[(n.b+1):length(lambda2.1)])
    }else{
      evv_ind<-t(eval.basis((seq(0,times[5],length.out=200))[ind],basis))
      ev<-diag(d.sem2)^(1/2)*t(eval.basis(times,basis))
      evv<-(diag(d.sem2)^(1/2))[5]*kronecker(ev[,5],evv_ind)
      as.vector(t(evv)%*%lambda2.1[(n.b+1):length(lambda2.1)])
    }
  }
}

per_fac_compute <- function(re, lambda1, samples_list, boot){
  out<-numeric(no.f)
  for(j in 1:no.f){
    lambda1.1<-as.matrix(t(sapply(1:boot,function(k) lambda1[[k]][[j]])))
    count3<-which(lat2==re$model.fit$mod$factorModel[[j]]$factor)
    res<-re$model.fit$var$indicators[j]; cov<-lat2[count3]
    eff<-re$model.fit$mod$factorModel[[j]]$effect
    if(eff%in%c("fixed_concurrent","fixed_historical") || re$model.sim$mod$factorModel[[j]]$effect=="concurrent"){
      lam.fac<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind]
      lam.org<-re$simulation$params.org.eval$coef.fac.std[[res]][[cov]][ind]
    }else{
      lam.fac<-re$estimation$result$params.estimated$lambda.param.std[[j]][[count3]][ind,ind][,5]
      lam.org<-re$simulation$params.org.eval$coef.fac.std[[res]][[cov]][ind,ind][,5]
    }
    lam.fac<-as.vector(lam.fac); cov.fac<-as.matrix(cov(lambda1.1))
    fb<-fregion::fregion.band(lam.fac,cov.fac,type="BEc",conf.level=0.95)
    out[j]<-mean(fb[,3]<lam.org & lam.org<fb[,2])
  }
  out
}

boot<-40
set.seed(99); samples_list<-lapply(1:boot,function(k) sort(sample(samp,replace=TRUE)))
# ORIG
lam_o<-vector("list",boot)
for(k in 1:boot){ sc<-samples_list[[k]]; ss<-sigsem_orig(re,sc); lam_o[[k]]<-vector("list",no.f)
  for(j in 1:no.f){ l21<-lam_orig(re,sc,j); lam_o[[k]][[j]]<-transform_lambda(re,j,ss,l21) } }
per_o<-per_fac_compute(re,lam_o,samples_list,boot)
# OPT (skip inner for fixed_concurrent j=1)
lam_p<-vector("list",boot)
for(k in 1:boot){ sc<-samples_list[[k]]; ss<-sigsem_orig(re,sc); lam_p[[k]]<-vector("list",no.f)
  for(j in 1:no.f){ if(re$model.fit$mod$factorModel[[j]]$effect%in%c("fixed_concurrent","fixed_historical")){ l21<-NULL }else{ l21<-lam_opt(re,sc,j) }
    lam_p[[k]][[j]]<-transform_lambda(re,j,ss,l21) } }
per_p<-per_fac_compute(re,lam_p,samples_list,boot)
cat("per.fac ORIG:", round(per_o,4), "\n")
cat("per.fac OPT :", round(per_p,4), "\n")
cat("max|diff|:", max(abs(per_o-per_p)), "\n")
