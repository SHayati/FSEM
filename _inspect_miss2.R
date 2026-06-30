setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
r2 <- readRDS("outputs/result_parameter_set1_regular.missing_1_N_100_M_20.rds")
d <- r2$simulation$data
N <- length(unique(d$.id)); no.ind <- length(unique(d$.ind))
cat("N:", N, " indicators:", no.ind, "\n")
tab <- table(d$.id, d$.ind)
cat("obs per (subject,indicator): mean", round(mean(tab),2), " range", paste(range(tab),collapse="-"), "\n")
cat("total rows:", nrow(d), " full would be N*8*ind =", N*8*no.ind, " frac observed:", round(nrow(d)/(N*8*no.ind),3), "\n")
e2<-r2$estimation$result$params.estimated.eval; o2<-r2$simulation$params.org.eval
cat("\nfactor eigval.std sum est:", round(sum(sapply(e2$sigma.fac.eigval.std,sum)),3),
    " org:", round(sum(sapply(o2$sigma.fac.eigval.std,sum)),3),"\n")
cat("intercept est mean:", round(mean(e2$intercept.std[[1]]),3), " org:", round(mean(o2$intercept.std[[1]]),3),"\n")
# Compare: do MAIN sim files also have missingness? check main file obs density
cat("\n=== MAIN regular_1_N_100_M_20 obs density ===\n")
rm(r2,d); gc()
r1 <- readRDS("outputs/result_parameter_set1_regular_1_N_100_M_20.rds")
d1<-r1$simulation$data; tab1<-table(d1$.id,d1$.ind)
cat("main obs per (subj,ind): mean", round(mean(tab1),2)," rows:", nrow(d1)," ids:", length(unique(d1$.id)),"\n")
