## Shared coverage function for Table 2 (one Monte Carlo file -> coverage per factor).
suppressMessages({library(Matrix); library(fda)})

BOOT <- 200; MONT <- 100; NOIT <- 10
PARTDIR <- "_cov2_parts"

cov_one_file <- function(f, boot=BOOT){
  re <- readRDS(f)
  mont<-MONT; no.it<-NOIT
  no.f<-length(re$model.fit$mod$factorModel)
  lat2<-re$model.fit$var$latents
  n.b<-length(re$estimation$result$params.estimated$intercept[[1]])
  ind<-seq(1,200,by=2)
  basis<-create.bspline.basis(nbasis=n.b,rangeval=c(0,1))
  times<-(seq(0,1,length.out=200))[ind]
  eval1<-t(eval.basis(times,basis))
  Bt<-t(eval.basis(times,basis))
  evv_ind<-t(eval.basis((seq(0,times[5],length.out=200))[ind],basis))
  samp<-unique(re$simulation$data$.id)
  lat.count<-which(lat2==re$model.fit$mod$regression[[1]]$response)
  N<-length(samp)
  M<-max(sapply(samp,function(j){max(sapply(1:no.f,function(jj){length(re$simulation$data$.t[re$simulation$data$.id==j & re$simulation$data$.ind==jj])}))}))

  hf1<-re$estimation$hist$hist.fac[[1]]
  Eall<-vector("list",mont)
  for(m in 1:mont){ Em<-matrix(0,n.b,N); for(s in 1:N) Em[,s]<-hf1$eta.l[[m]][[samp[s]]][[lat.count]]; Eall[[m]]<-Em }

  lam_opt <- function(samples,j){
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
  transform_lambda <- function(j,sigsem_k,lambda2.1){
    d.sem2<-t(eval1)%*%sigsem_k%*%eval1
    if(re$model.fit$mod$factorModel[[j]]$effect%in%c("fixed_concurrent","fixed_historical")){
      as.vector(diag(d.sem2)^(1/2))
    }else if(re$model.sim$mod$factorModel[[j]]$effect=="concurrent"){
      ev<-diag(d.sem2)^(1/2)*Bt; as.vector(t(ev)%*%lambda2.1[(n.b+1):length(lambda2.1)])
    }else{
      ev<-diag(d.sem2)^(1/2)*Bt; evv<-(diag(d.sem2)^(1/2))[5]*kronecker(ev[,5],evv_ind)
      as.vector(t(evv)%*%lambda2.1[(n.b+1):length(lambda2.1)])
    }
  }

  lambda1<-vector("list",boot)
  for(k in 1:boot){
    samples<-sort(sample(samp,replace=TRUE)); sidx<-match(samples,samp)
    sigsem<-matrix(0,n.b,n.b)
    for(m in 1:mont){ Em<-Eall[[m]][,sidx,drop=FALSE]; sigsem<-sigsem+tcrossprod(Em) }
    sigsem<-sigsem/(length(samples)*mont)
    lambda1[[k]]<-vector("list",no.f)
    for(j in 1:no.f){
      l21<-if(re$model.fit$mod$factorModel[[j]]$effect%in%c("fixed_concurrent","fixed_historical")) NULL else lam_opt(samples,j)
      lambda1[[k]][[j]]<-transform_lambda(j,sigsem,l21)
    }
  }
  per<-numeric(no.f)
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
    per[j]<-mean(fb[,3]<lam.org & lam.org<fb[,2])
  }
  design<-if(grepl("_irregular_",f)) "irregular" else "regular"
  data.frame(design=design,N=N,M=M,per1=per[1],per2=per[2],per3=per[3],file=basename(f),stringsAsFactors=FALSE)
}
