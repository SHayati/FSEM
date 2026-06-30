setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
files <- list.files("outputs", pattern="^result_parameter_set1_regular\\.missing_[0-9]+_N_100_M_20\\.rds$", full.names=TRUE)
# per-run MSE distribution for lambda2, lambda3, gamma1, gamma2
L <- matrix(NA, length(files), 3); G <- matrix(NA, length(files), 2)
for(i in seq_along(files)){
  res<-readRDS(files[i]); est<-res$estimation$result$params.estimated.eval; org<-res$simulation$params.org.eval
  for(j in 1:3){ ce<-est$coef.fac.std[[j]][[1]]; co<-org$coef.fac.std[[j]][[1]]; L[i,j]<-mean((ce-co)^2) }
  cse<-est$coef.sem.std[["eta"]]; cso<-org$coef.sem.std[["eta"]]
  for(k in 1:2){ G[i,k]<-mean((cse[[k]]-cso[[k]])^2) }
  rm(res); if(i%%25==0) gc()
}
cat("lambda1 quantiles:", round(quantile(L[,1]),3),"\n")
cat("lambda2 quantiles:", round(quantile(L[,2]),3),"\n")
cat("lambda3 quantiles:", round(quantile(L[,3]),3),"\n")
cat("gamma1  quantiles:", round(quantile(G[,1]),3),"\n")
cat("gamma2  quantiles:", round(quantile(G[,2]),3),"\n")
cat("\nfraction lambda2 < 0.1:", mean(L[,2]<0.1), " lambda3<0.1:", mean(L[,3]<0.1),"\n")
saveRDS(list(L=L,G=G),"_s1_diag.rds")
