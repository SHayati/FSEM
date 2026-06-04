
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


#Table 1:regular design
##fitting procedure
### set n.sim (number of monte carlo simulations)
n.sim = 100
results<-list()
comb<-expand.grid(N=c(100),M=c(10,20))
count<-0

# number of cores
n.cores <- 20#max(1, detectCores() - 1)
cl <- makeCluster(n.cores)
registerDoParallel(cl)

n.sim = n.cores*4

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



for (j in 1:nrow(comb)) {
  N<-comb[j,]$N
  M<-comb[j,]$M
  n.b    <- 6
  n.em   <- 100      # EM iterations
  n.monte<- 100      # Monte-Carlo draws in the E-step

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


                     # count<-count+1
                     model.sim <- fsem(eta1~~z1+z2, effectType="concurrent")          %+%
                       fsem(eta1~~z3, effectType="historical")                        %+%
                       fsem(eta2~~z4+z5+z6, effectType="concurrent")                  %+%
                       fsem(eta1~-1+eta2, effectType="concurrent",
                            latent.covariate="eta2", scalar.covariate=FALSE)          %+%
                       fsem(eta2~-1)
                     
                     
                     ## pin one loading per latent (the lambda fix) - misspecified 
                     model.fit.m <- fsem(eta1~~z1, effectType="fixed_concurrent")      %+%
                       fsem(eta1~~z2, effectType="concurrent")                        %+%
                       fsem(eta1~~z3, effectType="historical")                        %+%
                       fsem(eta2~~z4, effectType="fixed_concurrent")                  %+%
                       fsem(eta2~~z5+z6, effectType="concurrent")                     %+%
                       fsem(eta1~-1, effectType="concurrent")          %+%
                       fsem(eta2~-1)
                     
                     ## pin one loading per latent (the lambda fix) - wellspecified                      
                     model.fit.w <- fsem(eta1~~z1, effectType="fixed_concurrent")      %+%
                       fsem(eta1~~z2, effectType="concurrent")                        %+%
                       fsem(eta1~~z3, effectType="historical")                        %+%
                       fsem(eta2~~z4, effectType="fixed_concurrent")                  %+%
                       fsem(eta2~~z5+z6, effectType="concurrent")                     %+%
                       fsem(eta1~-1+eta2, effectType="concurrent", latent.covariate="eta2", scalar.covariate=FALSE)          %+%
                       fsem(eta2~-1)
                     
                     
                     simm <- simulation(model=model.sim, design="regular", n.t=M, n.b=n.b, r=1,
                                        rho=0.3, SNR=4, n.sample=N, Matern.sem=TRUE, Matern.fac=FALSE,
                                        parameters="parameter_set1")
                     
                     ## model fitting  (revised estimator: gamma-jitter + scale anchor)
                     init.m <- initial.param(model=model.fit.m, data=simm$data, n.b=n.b)
                     params.estimated.m  <- em.estimation.optimised.revised(model=model.fit.m, data=simm$data,
                                                                            initial.parameter=init.m, n.b=n.b,
                                                                            n.em=n.em, n.monte=n.monte, s.p="min",
                                                                            design="regular", plot.progress=FALSE,
                                                                            gamma.init.sd=0.3, seed=1)
                     
                     init.w <- initial.param(model=model.fit.w, data=simm$data, n.b=n.b)
                     params.estimated.w  <- em.estimation.optimised.revised(model=model.fit.w, data=simm$data,
                                                                            initial.parameter=init.w, n.b=n.b,
                                                                            n.em=n.em, n.monte=n.monte, s.p="min",
                                                                            design="regular", plot.progress=FALSE,
                                                                            gamma.init.sd=0.3, seed=1)
                     
                     
                     res<-list(simulation=simm,model.sim=model.sim,estimation.m=params.estimated.m,estimation.w=params.estimated.w)
                     
                     
                     file_name = paste0("result_",i,"_N_",N,"_M_",M,".rds")
                     if(!dir.exists("outputs3"))dir.create("outputs3")
                     
                     saveRDS(res, file.path("outputs3",file_name), compress = FALSE)
                     res
                   }
  
  
  for (k in seq_along(res_j)) {
    count <- count + 1
    results[[count]] <- res_j[[k]]
  }
  
}

stopCluster(cl)
if(!dir.exists("outputs3"))dir.create("outputs3")
saveRDS(results, file.path("outputs3","results3.rds"), compress = FALSE)

#extracting MSE values (two-factor model, DEPENDENT factors)
# The fitted models are only created inside the worker loop, so rebuild them here
# for the MSE summary.  The generalised streaming engine lives at the bottom of
# Core/MSE_tab3.R; FSEM_MSE_DEFINE_ONLY loads the helper functions only.
FSEM_MSE_DEFINE_ONLY <- TRUE
source("Core/MSE_tab3.R")

model.fit.m <- fsem(eta1~~z1, effectType="fixed_concurrent")      %+%
  fsem(eta1~~z2, effectType="concurrent")                        %+%
  fsem(eta1~~z3, effectType="historical")                        %+%
  fsem(eta2~~z4, effectType="fixed_concurrent")                  %+%
  fsem(eta2~~z5+z6, effectType="concurrent")                     %+%
  fsem(eta1~-1, effectType="concurrent")                         %+%
  fsem(eta2~-1)
model.fit.w <- fsem(eta1~~z1, effectType="fixed_concurrent")      %+%
  fsem(eta1~~z2, effectType="concurrent")                        %+%
  fsem(eta1~~z3, effectType="historical")                        %+%
  fsem(eta2~~z4, effectType="fixed_concurrent")                  %+%
  fsem(eta2~~z5+z6, effectType="concurrent")                     %+%
  fsem(eta1~-1+eta2, effectType="concurrent", latent.covariate="eta2", scalar.covariate=FALSE) %+%
  fsem(eta2~-1)

mse.well <- fsem_mse_from_list(results, fit_key = "model.fit.w",
                               est_key = "estimation.w", fit_model = model.fit.w)
mse.miss <- fsem_mse_from_list(results, fit_key = "model.fit.m",
                               est_key = "estimation.m", fit_model = model.fit.m)
write.csv(mse.well$tabl.fac, file.path("outputs3","MSE_dep_well_fac.csv"), row.names = FALSE)
write.csv(mse.well$tabl.sem, file.path("outputs3","MSE_dep_well_sem.csv"), row.names = FALSE)
write.csv(mse.miss$tabl.fac, file.path("outputs3","MSE_dep_miss_fac.csv"), row.names = FALSE)
write.csv(mse.miss$tabl.sem, file.path("outputs3","MSE_dep_miss_sem.csv"), row.names = FALSE)
print(mse.well$tabl.fac); print(mse.well$tabl.sem)
print(mse.miss$tabl.fac); print(mse.miss$tabl.sem)


