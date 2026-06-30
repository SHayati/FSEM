setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
# Compare z2 standardized loading: main sim (eta~-1) vs S1 (eta~x1+x2)
cat("=== MAIN sim regular_1_N_100_M_20 (z2 concurrent) ===\n")
r1 <- readRDS("outputs/result_parameter_set1_regular_1_N_100_M_20.rds")
e1<-r1$estimation$result$params.estimated.eval; o1<-r1$simulation$params.org.eval
for(j in 1:3){ ce<-e1$coef.fac.std[[j]][[1]]; co<-o1$coef.fac.std[[j]][[1]]
  cat("FM",j,"est mean:",round(mean(ce),3)," org mean:",round(mean(co),3)," mse:",round(mean((ce-co)^2),4),"\n") }
cat("eigval.std fac org:", round(sapply(o1$sigma.fac.eigval.std,function(x)x[1]),3),"\n")
cat("eigval.std fac est:", round(sapply(e1$sigma.fac.eigval.std,function(x)x[1]),3),"\n")

cat("\n=== S1 regular.missing_1_N_100_M_20 (z2,z3 concurrent, eta~x1+x2) ===\n")
r2 <- readRDS("outputs/result_parameter_set1_regular.missing_1_N_100_M_20.rds")
e2<-r2$estimation$result$params.estimated.eval; o2<-r2$simulation$params.org.eval
for(j in 1:3){ ce<-e2$coef.fac.std[[j]][[1]]; co<-o2$coef.fac.std[[j]][[1]]
  cat("FM",j,"est mean:",round(mean(ce),3)," org mean:",round(mean(co),3)," mse:",round(mean((ce-co)^2),4),"\n") }
cat("eigval.std fac org:", round(sapply(o2$sigma.fac.eigval.std,function(x)x[1]),3),"\n")
cat("eigval.std fac est:", round(sapply(e2$sigma.fac.eigval.std,function(x)x[1]),3),"\n")
cat("eigval.std sem org:", round(sapply(o2$sigma.sem.eigval.std,function(x)x[1]),3),"\n")
cat("eigval.std sem est:", round(sapply(e2$sigma.sem.eigval.std,function(x)x[1]),3),"\n")
# Is there a NON-std version that matches better?
cat("\n--- non-std coef.fac (S1) ---\n")
cat("names params.estimated:", names(r2$estimation$result$params.estimated)[grepl('lambda|coef|fac',names(r2$estimation$result$params.estimated))],"\n")
le<-r2$estimation$result$params.estimated$lambda.param.std
cat("lambda.param.std length:", length(le), "names j1:", names(le[[1]]),"\n")
for(j in 1:3){ m<-le[[j]][[1]]; cat("FM",j,"lambda.param.std dim:",paste(dim(as.matrix(m)),collapse='x')," mean:",round(mean(m),3),"\n") }
co<-r2$simulation$params.org.eval$coef.fac.std
