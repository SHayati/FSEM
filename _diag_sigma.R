setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
options(width=200)
# inspect sigma.error per-run for FM3 & FM2, regular N=50 M=10
vals <- matrix(NA, nrow=100, ncol=3)
for (i in 1:100){
  f <- sprintf("outputs/result_parameter_set1_regular_%d_N_50_M_10.rds", i)
  if(!file.exists(f)) next
  r <- readRDS(f)
  for(j in 1:3){
    est <- r$estimation$result$params.estimated$totparam$sigma.error[[j]]
    org <- r$simulation$params.org.eval$sigma.error[[j]]
    vals[i,j] <- (est-org)^2
  }
  rm(r)
}
for(j in 1:3){
  cat(sprintf("FM%d: mean=%.4f  trim10=%.5f  median=%.5f  max=%.4f  #>0.01=%d\n",
      j, mean(vals[,j],na.rm=TRUE), mean(vals[,j],trim=0.1,na.rm=TRUE),
      median(vals[,j],na.rm=TRUE), max(vals[,j],na.rm=TRUE), sum(vals[,j]>0.01,na.rm=TRUE)))
}
cat("\nFM3 sorted top 10:\n"); print(round(sort(vals[,3],decreasing=TRUE)[1:10],4))
cat("\nFM2 sorted top 10:\n"); print(round(sort(vals[,2],decreasing=TRUE)[1:10],4))
