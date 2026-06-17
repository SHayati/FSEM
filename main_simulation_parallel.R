
#initial libraries
library(fda)
library(Matrix)
library(mvtnorm)
library(tmvtnorm)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)


source("Core/FSEM.R")
source("Core/simualtion_utilities.R")
source("Core/FSEM_optimised.R")
source("Core/fsem_kernels_fix.R")         # collision-free Rcpp kernel loader (fixes Error 322)
source("Core/FSEM_optimised_revised.R")   # provides em.estimation.optimised.revised

# number of cores
n.cores <- 20#max(1, detectCores() - 1)
cl <- makeCluster(n.cores)
registerDoParallel(cl)

n.sim = n.cores*4
n.b    <- 6
n.em   <- 100      # EM iterations
n.monte<- 100      # Monte-Carlo draws in the E-step


# load packages on workers
clusterEvalQ(cl, {
  library(fda)
  library(Matrix)
  library(mvtnorm)
  library(tmvtnorm)
  library(MASS)
})

# source required files on workers
clusterExport(
  cl,
  varlist = c("fsem", "simulation", "initial.param"),
  envir = environment()
)


clusterEvalQ(cl, {
  setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
  source("Core/FSEM.R")
  source("Core/simualtion_utilities.R")
  source("Core/FSEM_optimised.R")
  source("Core/fsem_kernels_fix.R")         # collision-free Rcpp kernel loader (fixes Error 322)
  source("Core/FSEM_optimised_revised.R")   # provides em.estimation.optimised.revised
})

#Table 1:regular design
##fitting procedure
### set n.sim (number of monte carlo simulations)
n.sim = 100
results<-list()
comb<-expand.grid(N=c(50,100),M=c(10,20))
count<-0
for (j in 1:nrow(comb)) {
  N<-comb[j,]$N
  M<-comb[j,]$M
  # for (i in 1:n.sim) { # use parallel computing if possible
  res_j <- foreach(i = 1:n.sim,
                   .packages = c("fda","Matrix","mvtnorm","tmvtnorm","MASS"),
                   .export = c("N", "M", "n.b", "n.em", "n.monte"),
                   # CRITICAL: do NOT let foreach auto-export the FSEM functions from the
                   # master. em.estimation.optimised(.revised) closes over the master's
                   # .fsem_kernel_env, which holds Rcpp native symbol pointers that are
                   # only valid in the master process. Serializing those closures into a
                   # worker ships dead/NULL pointers that shadow the live kernels the
                   # worker compiled in clusterEvalQ -> the wrapper .fsem_xtSx() then
                   # raises "NULL value passed as symbol address". Listing them in
                   # .noexport forces each worker to use its own clusterEvalQ-sourced
                   # copies (with live kernels).
                   .noexport = c("fsem", "%+%", "simulation", "initial.param",
                                 "em.estimation.optimised", "em.estimation.optimised.revised",
                                 ".fsem_xtSx", ".fsem_xtSz", ".fsem_kernel_env", ".fsem_have_rcpp"),
                   .combine = "list") %dopar% {

    model.sim<-fsem(eta1~~z1+z2,effectType="concurrent")%+%
      fsem(eta1~~z3,effectType="historical")%+%
      fsem(eta1~-1)
    
    model.fit<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
      fsem(eta1~~z2,effectType="concurrent")%+%
      fsem(eta1~~z3,effectType="historical")%+%
      fsem(eta1~-1)
    
    simm<-simulation(model=model.sim,design="regular",n.t=M,n.b=6,r=1,rho = 0.3,SNR=4,
                     n.sample = N,Matern.sem = TRUE,Matern.fac = FALSE,parameters = "parameter_set1") 
    initial.parameter<-initial.param(model = model.fit,data=simm$data,n.b=6)
    params.estimated<-em.estimation.optimised.revised(model=model.fit,data=simm$data,
                                    initial.parameter=initial.parameter,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="regular",gamma.init.sd = 0.3)
    result<-list(simulation=simm,model.fit=model.fit,model.sim=model.sim,estimation=params.estimated)
    

    file_name = paste0("result_parameter_set1_regular_",i,"_N_",N,"_M_",M,".rds")
    if(!dir.exists("outputs"))dir.create("outputs")
    
    saveRDS(result, file.path("outputs",file_name), compress = FALSE)
    result
    
 }  
  
  for (k in seq_along(res_j)) {
    count <- count + 1
    results[[count]] <- res_j[[k]]
  }
}

if(!dir.exists("outputs"))dir.create("outputs")
saveRDS(results, file.path("outputs","results_parameter_set1_regular.rds"), compress = FALSE)

#extracting MSE values
source("Core/MSE_tab1.R")
print(tabl.fac)
print(tabl.sem)

#Table 2:regular design
##extracting coverage rates
source("Core/CovRate_tab2.R")
print(cov.table.fac)
print(cov.table.sem)

summary_parameter_set1 = list(tabl.fac=tabl.fac, tabl.sem=tabl.sem,cov.table.fac=cov.table.fac, cov.table.sem=cov.table.sem)
if(!dir.exists("outputs"))dir.create("outputs")
saveRDS(summary_parameter_set1, file.path("outputs","summary_parameter_set1.rds"), compress = FALSE)

#Table 1:irregular design
##fitting procedure
###set n.sim (number of monte carlo simulations)
results<-list()
comb<-expand.grid(N=c(50,100),M=c(10,20))
count<-0
for (j in 1:nrow(comb)) {
  N<-comb[j,]$N
  M<-comb[j,]$M
# for (i in 1:n.sim) { # use parallel computing if possible
  res_j <- foreach(i = 1:n.sim,
                   .packages = c("fda","Matrix","mvtnorm","tmvtnorm","MASS"),
                   .export = c("N", "M", "n.b", "n.em", "n.monte"),
                   # CRITICAL: do NOT let foreach auto-export the FSEM functions from the
                   # master. em.estimation.optimised(.revised) closes over the master's
                   # .fsem_kernel_env, which holds Rcpp native symbol pointers that are
                   # only valid in the master process. Serializing those closures into a
                   # worker ships dead/NULL pointers that shadow the live kernels the
                   # worker compiled in clusterEvalQ -> the wrapper .fsem_xtSx() then
                   # raises "NULL value passed as symbol address". Listing them in
                   # .noexport forces each worker to use its own clusterEvalQ-sourced
                   # copies (with live kernels).
                   .noexport = c("fsem", "%+%", "simulation", "initial.param",
                                 "em.estimation.optimised", "em.estimation.optimised.revised",
                                 ".fsem_xtSx", ".fsem_xtSz", ".fsem_kernel_env", ".fsem_have_rcpp"),
                   .combine = "list") %dopar% {
  model.sim<-fsem(eta1~~z1+z2,effectType="concurrent")%+%
    fsem(eta1~~z3,effectType="historical")%+%
    fsem(eta1~-1)
  
  model.fit<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
    fsem(eta1~~z2,effectType="concurrent")%+%
    fsem(eta1~~z3,effectType="historical")%+%
    fsem(eta1~-1)
  
  simm<-simulation(model=model.sim,design="irregular",n.t=M,n.b=6,r=1,rho = 0.3,SNR=4,
                   n.sample = N,Matern.sem = TRUE,Matern.fac = FALSE,parameters = "parameter_set1") 
  initial.parameter<-initial.param(model = model.fit,data=simm$data,n.b=6)
  params.estimated<-em.estimation.optimised.revised(model=model.fit,data=simm$data,
                                  initial.parameter=initial.parameter,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="irregular",gamma.init.sd = 0.3)
  result<-list(simulation=simm,model.fit=model.fit,model.sim=model.sim,estimation=params.estimated)
  
  file_name = paste0("result_parameter_set1_irregular_",i,"_N_",N,"_M_",M,".rds")
  if(!dir.exists("outputs"))dir.create("outputs")
  
  saveRDS(result, file.path("outputs",file_name), compress = FALSE)
  result
  
                   }  
  
  for (k in seq_along(res_j)) {
    count <- count + 1
    results[[count]] <- res_j[[k]]
  }
}

if(!dir.exists("outputs"))dir.create("outputs")
saveRDS(results, file.path("outputs","results_parameter_set1_irregular.rds"), compress = FALSE)

#extracting MSE values
source("Core/MSE_tab1.R")
print(tabl.fac)
print(tabl.sem)

#Table 2:irregular design
##extracting coverage rates
source("Core/CovRate_tab2.R")
print(cov.table.fac)
print(cov.table.sem)

summary_parameter_set1 = list(tabl.fac=tabl.fac, tabl.sem=tabl.sem,cov.table.fac=cov.table.fac, cov.table.sem=cov.table.sem)
if(!dir.exists("outputs"))dir.create("outputs")
saveRDS(summary_parameter_set1, file.path("outputs","summary_parameter_set1_irregular.rds"), compress = FALSE)

#Table 3:regular missing at random design for N=100 and M=8
##fitting procedure
###set n.sim (number of monte carlo simulations)
results<-list()
results <- foreach(i = 1:n.sim,
                 .packages = c("fda","Matrix","mvtnorm","tmvtnorm","MASS"),
                 .export = c("N", "M", "n.b", "n.em", "n.monte"),
                 # CRITICAL: do NOT let foreach auto-export the FSEM functions from the
                 # master. em.estimation.optimised(.revised) closes over the master's
                 # .fsem_kernel_env, which holds Rcpp native symbol pointers that are
                 # only valid in the master process. Serializing those closures into a
                 # worker ships dead/NULL pointers that shadow the live kernels the
                 # worker compiled in clusterEvalQ -> the wrapper .fsem_xtSx() then
                 # raises "NULL value passed as symbol address". Listing them in
                 # .noexport forces each worker to use its own clusterEvalQ-sourced
                 # copies (with live kernels).
                 .noexport = c("fsem", "%+%", "simulation", "initial.param",
                               "em.estimation.optimised", "em.estimation.optimised.revised",
                               ".fsem_xtSx", ".fsem_xtSz", ".fsem_kernel_env", ".fsem_have_rcpp"),
                 .combine = "list") %dopar% {
  model.sim<-fsem(eta~~z1,effectType = "concurrent")%+%
    fsem(eta~~z2,effectType="concurrent")%+%
    fsem(eta~~z3,effectType="concurrent")%+%
    fsem(eta~-1+x1,effectType="linear",scalar.covariate = TRUE)%+%
    fsem(eta~-1+x2,effectType="linear",scalar.covariate = TRUE)
  
  model.fit<-fsem(eta~~z1,effectType = "fixed_concurrent")%+%
    fsem(eta~~z2,effectType="concurrent")%+%
    fsem(eta~~z3,effectType="concurrent")%+%
    fsem(eta~-1+x1,effectType="linear",scalar.covariate = TRUE)%+%
    fsem(eta~-1+x2,effectType="linear",scalar.covariate = TRUE)
  
  simm<-simulation(model=model.sim,design="regular.missing",n.t=8,n.b=6,r=1,rho = 0.3,SNR=4,
                   n.sample = 100,Matern.sem = TRUE,Matern.fac = FALSE,parameters = "parameter_set2") 
  initial.parameter<-initial.param(model = model.fit,data=simm$data,n.b=6)
  params.estimated<-em.estimation.optimised.revised(model=model.fit,data=simm$data,x.data=simm$covariate,
                                  initial.parameter=initial.parameter,n.b=6,n.em=2,n.monte=2,s.p="min",range.min=NULL,range.max=NULL,design="regular.missing",gamma.init.sd = 0.3)
  result<-list(simulation=simm,model.fit=model.fit,model.sim=model.sim,estimation=params.estimated)
  
  
  file_name = paste0("result_parameter_set1_regular.missing_",i,"_N_",N,"_M_",M,".rds")
  if(!dir.exists("outputs"))dir.create("outputs")
  
  saveRDS(result, file.path("outputs",file_name), compress = FALSE)
  result
  
}  

if(!dir.exists("outputs"))dir.create("outputs")
saveRDS(results, file.path("outputs","results_parameter_set1_regular.missing.rds"), compress = FALSE)

#extracting MSE values
source("Core/MSE_tab3.R")
print(tabl.fac)
print(tabl.sem)

#Table 3:regular missing at random design for N=100 and M=8
##extracting coverage rates
source("Core/CovRate_tab3.R")
print(cov.table.fac)
print(cov.table.sem)

summary_parameter_set1 = list(tabl.fac=tabl.fac, tabl.sem=tabl.sem,cov.table.fac=cov.table.fac, cov.table.sem=cov.table.sem)
if(!dir.exists("outputs"))dir.create("outputs")
saveRDS(summary_parameter_set1, file.path("outputs","summary_parameter_set1_regular.missing.rds"), compress = FALSE)

stopCluster(cl)

