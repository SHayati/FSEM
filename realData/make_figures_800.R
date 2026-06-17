## make_figures_800.R
## Regenerates manuscript Figures 2, 3 and 4 (HRS GFP case study) from the
## re-estimated n=800 results (realData/result_800.rds).
## This is the "production" version of main_realdata.R: instead of drawing to the
## R device it writes the PDF images used by interactapasample_revised.tex into
## manuscript/images/.  Estimation is NOT re-run (result_800.rds is loaded).
##
## Confidence bands: reproduces the subject-bootstrap of covRate.R, with two fixes
##   (1) the foreach output is captured (the repo copy discards it and instead
##       reads a machine-specific file that does not exist here), and
##   (2) the GLS weighting uses kronecker(I_n, sigma) -- identity weighting, matching
##       the authoritative final-estimation M-step in Core/FSEM.R (lines 1149/1173).
##       The repo's covRate.R used kronecker(diag(<resampled-id-vector>), sigma),
##       which weights subjects by their id value -- a copy/paste bug.
##   The 4956x4956 solve is avoided via the exact sparse block form
##       sol = kronecker(Diagonal(n_boot), solve(sigma)).
setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
suppressMessages({library(Matrix); library(fda)})
stopifnot(requireNamespace("fregion", quietly = TRUE))
source("Core/FSEM.R"); source("Core/simualtion_utilities.R")

out   <- "manuscript/images"
boot  <- 200
seed  <- 20240617
results <- readRDS("realData/result_800.rds")
d.f.se  <- readRDS("realData/d.f.se_n800.rds")

## ----------------------------------------------------------------------------
## Figure 2 : convergence (tolerances) -- no bootstrap needed
## ----------------------------------------------------------------------------
co <- results$tolerance$tol[-1]       # total tolerance of functional parameters
sg <- results$tolerance$tol_sig[-1]   # tolerance of measurement-error variances

pdf(file.path(out, "tol_coef2.pdf"), width = 575/72, height = 418/72)
par(mgp = c(0, 2, 0))
plot(seq_along(co), co, "l", ylab="", xlab="", cex.axis=3, cex.lab=3, lwd=3, family="serif")
dev.off()

pdf(file.path(out, "tol_sig_2.pdf"), width = 573/72, height = 358/72)
par(mgp = c(0, 2, 0))
plot(seq_along(sg), sg, "l", ylab="", xlab="", cex.axis=3, cex.lab=3, lwd=3, family="serif")
dev.off()
cat("[fig2] convergence PDFs written\n"); flush.console()

## ----------------------------------------------------------------------------
## Model specification (same 1-factor GFP model as realDataSEFF.R)
## ----------------------------------------------------------------------------
model.dat.se<-fsem(eta~~z1,effectType="fixed_concurrent")%+%
  fsem(eta~~z2,effectType="concurrent")%+%fsem(eta~~z3,effectType="concurrent")%+%
  fsem(eta~~z4,effectType="concurrent")%+%fsem(eta~~z5,effectType="concurrent")%+%
  fsem(eta~-1+gender,effectType="linear",scalar.covariate=TRUE)%+%
  fsem(eta~-1+cancer,effectType="linear",scalar.covariate=TRUE)%+%
  fsem(eta~-1+diabet,effectType="linear",scalar.covariate=TRUE)%+%
  fsem(eta~-1+heartcon,effectType="linear",scalar.covariate=TRUE)%+%
  fsem(eta~-1+genCancer,effectType="linear",scalar.covariate=TRUE)%+%
  fsem(eta~-1+genDiabet,effectType="linear",scalar.covariate=TRUE)%+%
  fsem(eta~-1+genHeart,effectType="linear",scalar.covariate=TRUE)

re <- list(model.fit=model.dat.se, model.sim=model.dat.se, estimation=results, simmulation=d.f.se)
mont<-100; no.it<-10; ind<-seq(1,200,by=2)
no.f<-length(re$model.fit$mod$factorModel); no.r<-length(re$model.fit$mod$regression)
lat2<-re$model.fit$var$latents; lat2.cov<-re$model.fit$var$observed
n.b<-length(re$estimation$result$params.estimated$intercept[[1]])
basis<-create.bspline.basis(nbasis=n.b,rangeval=c(0,1))
times<-(seq(0,1,length.out=200))[ind]
eval1<-t(eval.basis(times,basis))
samp<-unique(re$simmulation$data$.id)
HSEM <- re$estimation$hist$hist.sem
HFAC <- re$estimation$hist$hist.fac
PE   <- re$estimation$result$params.estimated
eff.fac <- sapply(1:no.f, function(j) re$model.fit$mod$factorModel[[j]]$effect)

boot_once <- function(){
  samples<-sort(sample(samp,replace=TRUE)); nbk<-length(samples); Ds<-Diagonal(nbk)
  gamma1<-list(); lambda1<-list(); sig.sem<-list()
  for (j in 1:no.r) {
    lat.count<-which(lat2==re$model.fit$mod$regression[[j]]$response)
    gamma1[[lat.count]]<-list()
    if(!is.null(re$model.fit$mod$regression[[j]]$covariate)){
      la<-length(re$model.fit$mod$regression[[j]]$covariate)
      gg <- do.call(rbind, lapply(samples, function(i1) HSEM[[lat.count]]$g.x[[i1]]))
      ETT<- vapply(1:mont, function(m) unlist(lapply(samples, function(i1) HFAC[[1]]$eta.l[[m]][[i1]][[lat.count]])), numeric(nrow(gg)))
      Gm <- lapply(1:mont, function(m) do.call(rbind, lapply(samples, function(i1) HSEM[[lat.count]]$g[[m]][[i1]])))
      sig<-PE$sigma.sem[[lat.count]]; deltaa<-HSEM[[lat.count]]$deltaa; pen.sem<-HSEM[[lat.count]]$pen.sem
      ett_bar<-rowMeans(ETT)
      for (it in 1:no.it){
        sol<-kronecker(Ds, solve(sig)); tgs<-t(gg)%*%sol
        gam.sum <- as.matrix(tgs%*%gg)+deltaa*pen.sem
        gam1.sum<- as.matrix(tgs%*%ett_bar)
        gamma2.1<- solve(gam.sum)%*%gam1.sum
        S<-matrix(0,n.b,n.b)
        for (m in 1:mont){ r<-ETT[,m]-as.numeric(Gm[[m]]%*%gamma2.1); Rm<-matrix(r,nrow=n.b); S<-S+Rm%*%t(Rm) }
        sig<-S/(nbk*mont)
      }
      sig.sem[[lat.count]]<-sig
      d.sem1<-t(eval1)%*%sig.sem[[lat.count]]%*%eval1; count<-0
      for (jj in 1:la){
        count3<-which(lat2.cov==re$model.sim$mod$regression[[j]]$covariate[[jj]])
        ga<-gamma2.1[(count+1):(count+n.b)]; count<-count+n.b
        gamma1[[lat.count]][[count3]]<-as.vector(diag(diag(d.sem1)^(-1/2))%*%t(eval1)%*%ga)
      }
    }
  }
  for (j in 1:no.f){
    count3<-which(lat2==re$model.sim$mod$factorModel[[j]]$factor)
    ffm<-lapply(1:mont, function(m) do.call(rbind, lapply(samples, function(i1) HFAC[[j]]$f[[m]][[paste0(i1)]])))
    L<-nrow(ffm[[1]])
    zzm<-vapply(1:mont, function(m) unlist(lapply(samples, function(i1) HFAC[[j]]$generators2.z[[paste0(i1)]][[j]][[m]])), numeric(L))
    f2m<-lapply(1:mont, function(m) rep_len(unlist(lapply(samples, function(i1) HFAC[[j]]$f2[[m]][[paste0(i1)]])), L))
    sig<-PE$sigma.fac[[j]]; deltaa<-HFAC[[j]]$deltaa; pen.fac<-HFAC[[j]]$pen.fac
    for (it in 1:no.it){
      sol<-kronecker(Ds, solve(sig)); lam.sum<-0; lam1.sum<-0
      for (m in 1:mont){ tfs<-t(ffm[[m]])%*%sol; lam.sum<-lam.sum+as.matrix(tfs%*%ffm[[m]]); lam1.sum<-lam1.sum+as.matrix(tfs%*%(zzm[,m]-f2m[[m]])) }
      lam.sum<-lam.sum/mont+deltaa*pen.fac; lam1.sum<-lam1.sum/mont
      lambda2.1<-solve(lam.sum)%*%lam1.sum
      S<-matrix(0,n.b,n.b)
      for (m in 1:mont){ r<-zzm[,m]-f2m[[m]]-as.numeric(ffm[[m]]%*%lambda2.1); Rm<-matrix(r,nrow=n.b); S<-S+Rm%*%t(Rm) }
      sig<-S/(nbk*mont)
    }
    d.sem2<-t(eval1)%*%sig.sem[[count3]]%*%eval1
    if(eff.fac[j]%in%c("fixed_concurrent","fixed_historical")){
      lambda1[[j]]<-as.vector(diag(d.sem2)^(1/2))
    }else{
      ev<-diag(d.sem2)^(1/2)*t(eval.basis(times,basis))
      lambda1[[j]]<-as.vector(t(ev)%*%lambda2.1[(n.b+1):length(lambda2.1)])
    }
  }
  list(gamma1=gamma1, lambda1=lambda1)
}

## ----------------------------------------------------------------------------
## Bootstrap (reuse cached draws if present)
## ----------------------------------------------------------------------------
ress.file <- "realData/res_parallel_con_800.rds"
if(file.exists(ress.file)){
  ress<-readRDS(ress.file); cat("[boot] loaded cached", length(ress), "draws\n")
}else{
  set.seed(seed)
  ress<-vector("list", boot); t0<-Sys.time()
  for(k in 1:boot){
    ress[[k]]<-boot_once()
    if(k%%10==0) cat(sprintf("[boot] %d/%d  (%.1f min elapsed)\n", k, boot,
                             as.numeric(difftime(Sys.time(),t0,units="mins")))) ; flush.console()
  }
  saveRDS(ress, ress.file)
  cat(sprintf("[boot] done %d draws in %.1f min\n", boot, as.numeric(difftime(Sys.time(),t0,units="mins"))))
}

## ----------------------------------------------------------------------------
## Confidence bands (BEc, 95%) -- centres from the point estimates
## ----------------------------------------------------------------------------
fregion.bands.fac <- vector("list", no.f)
for(j in 1:no.f){
  lam.center <- as.vector(PE$lambda.param.std[[j]][[1]][ind])
  lam.draws  <- t(sapply(1:boot, function(k) ress[[k]]$lambda1[[j]]))
  fregion.bands.fac[[j]] <- fregion::fregion.band(lam.center, cov(lam.draws), type=c("BEc"), conf.level=c(0.95))
}
fregion.bands.sem <- vector("list", length(lat2.cov))
for(cc in 1:length(lat2.cov)){
  gam.center <- as.vector(PE$gamma.param.x.std[[1]][[cc]][ind])
  gam.draws  <- t(sapply(1:boot, function(k) ress[[k]]$gamma1[[1]][[cc]]))
  fregion.bands.sem[[cc]] <- fregion::fregion.band(gam.center, cov(gam.draws), type=c("BEc"), conf.level=c(0.95))
}
cat("[bands] computed\n"); flush.console()

## ----------------------------------------------------------------------------
## Figure 3 : factor loadings (lambda1..lambda5)  -- ylim auto-scaled to n=800
## Figure 4 : regression coefficients (gender, cancer)
## All plotted values negated, matching main_realdata.R styling exactly.
## ----------------------------------------------------------------------------
tt<-seq(2008,2022,length.out=100)
ylim_auto<-function(B){ r<-range(-B[,1],-B[,2],-B[,3]); pad<-0.06*diff(r); c(r[1]-pad, r[2]+pad) }
draw_band<-function(file,B,W=573/72,H=358/72){
  pdf(file.path(out,file), width=W, height=H)
  par(mgp=c(0,2,0))
  plot(tt,-B[,1],col="blue","l",ylim=ylim_auto(B),xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
  polygon(c(tt,rev(tt)), c(-B[,3],rev(-B[,2])), col=rgb(0.1,0.1,1,alpha=0.1), border=NA)
  dev.off()
}

load.files<-paste0("load",1:5,"_2.pdf")
for(j in 1:5) draw_band(load.files[j], fregion.bands.fac[[j]])

draw_band("gender2.pdf", fregion.bands.sem[[1]])  # covariate 1 = gender
draw_band("cancer2.pdf", fregion.bands.sem[[2]])  # covariate 2 = cancer
cat("[fig3/4] loading & regression PDFs written\n")
cat("ALL DONE\n")
