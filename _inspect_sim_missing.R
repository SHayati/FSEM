setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
options(width=200)
fs <- list.files("outputs", pattern="regular.missing_.*_N_100_M_20.rds", full.names=TRUE)
cat("num regular.missing files:", length(fs), "\n")
f1 <- fs[1]
cat("file:", basename(f1), "\n")
r <- readRDS(f1)
cat("indicators:", r$model.fit$var$indicators, "\n")
cat("latents:", r$model.fit$var$latents, "\n")
cat("factors:", r$model.fit$var$factors, "\n")
cat("n factorModel:", length(r$model.fit$mod$factorModel), " n regression:", length(r$model.fit$mod$regression), "\n")
cat("\n regression structure:\n")
for(j in seq_along(r$model.fit$mod$regression)){
  rg<-r$model.fit$mod$regression[[j]]
  cat(" reg",j,"response=",rg$response," covariate=",paste(unlist(rg$covariate),collapse="/")," effect=",paste(unlist(rg$effect),collapse="/"),"\n")
}
cat("\n factorModel effects:\n")
for(j in seq_along(r$model.fit$mod$factorModel)){
  fm<-r$model.fit$mod$factorModel[[j]]
  cat(" fm",j,"indicator=",fm$indicator," factor=",paste(unlist(fm$factor),collapse="/")," effect=",paste(unlist(fm$effect),collapse="/"),"\n")
}
cat("\n coef.sem.std (est) structure:\n"); str(r$estimation$result$params.estimated.eval$coef.sem.std, max.level=3)
cat("\n coef.sem.std (org) structure:\n"); str(r$simulation$params.org.eval$coef.sem.std, max.level=3)
cat("\n coef.fac.std (est):\n"); str(r$estimation$result$params.estimated.eval$coef.fac.std, max.level=2)
cat("\n N ids:", length(unique(r$simulation$data$.id)), "\n")
M<-max(sapply(unique(r$simulation$data$.id),function(j){max(sapply(1:length(r$model.fit$mod$factorModel), function(jj){length(r$simulation$data$.t[r$simulation$data$.id==j&r$simulation$data$.ind==jj])}))}))
cat("M (max time pts):", M, "\n")
cat("gamma.param.x.std present?", !is.null(r$estimation$result$params.estimated$gamma.param.x.std), "\n")
cat("size MB:", round(file.size(f1)/1e6,1), "\n")
