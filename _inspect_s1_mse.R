setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
f <- list.files("outputs", pattern="^result_parameter_set1_regular\\.missing_1_N_100_M_20\\.rds$", full.names=TRUE)
res <- readRDS(f)
est <- res$estimation$result$params.estimated.eval
org <- res$simulation$params.org.eval
cat("indicators:", res$model.fit$var$indicators, "\n")
cat("latents:", res$model.fit$var$latents, "\n")
cat("factorModel effects:\n")
for(j in seq_along(res$model.fit$mod$factorModel)) cat("  j",j,res$model.fit$mod$factorModel[[j]]$indicator, res$model.fit$mod$factorModel[[j]]$effect, "\n")
cat("\ncoef.fac.std length:", length(est$coef.fac.std), "\n")
for(j in seq_along(est$coef.fac.std)){
  cat(" j",j,"names:", names(est$coef.fac.std[[j]]), " nsub:", length(est$coef.fac.std[[j]]), "\n")
  if(length(est$coef.fac.std[[j]])>0){
    for(k in seq_along(est$coef.fac.std[[j]])){
      ce<-est$coef.fac.std[[j]][[k]]; co<-org$coef.fac.std[[j]][[k]]
      cat("   k",k,"dim est:", paste(dim(as.matrix(ce)),collapse="x"), " len:",length(ce),
          " mse:", round(1/200*sum((ce-co)^2),4),
          " est[1:3]:", paste(round(head(as.vector(ce),3),3),collapse=","),
          " org[1:3]:", paste(round(head(as.vector(co),3),3),collapse=","), "\n")
    }
  }
}
cat("\ncoef.sem.std names:", names(est$coef.sem.std), "\n")
cse<-est$coef.sem.std[[res$model.fit$var$latents[1]]]
cso<-org$coef.sem.std[[res$model.fit$var$latents[1]]]
cat("sub names:", names(cse), " n:", length(cse), "\n")
for(k in seq_along(cse)){
  cat(" k",k,"len:",length(cse[[k]])," mse:", round(1/200*sum((cse[[k]]-cso[[k]])^2),4),
      " est[1:3]:", paste(round(head(as.vector(cse[[k]]),3),3),collapse=","),
      " org[1:3]:", paste(round(head(as.vector(cso[[k]]),3),3),collapse=","),"\n")
}
cat("\nregression effect:", sapply(res$model.fit$mod$regression, function(r) r$effect), "\n")
