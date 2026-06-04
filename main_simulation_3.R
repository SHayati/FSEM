
# Simulating and fitting functional structural equation models (FSEM) for regular design with 6 covariates and 2 dependent functional latent variables.

#initial libraries
library(fda)
library(Matrix)
library(mvtnorm)
library(tmvtnorm)
library(MASS)

source("Core/FSEM.R")
source("Core/simualtion_utilities.R")


#Table 1:regular design
##fitting procedure
### set n.sim (number of monte carlo simulations)
n.sim = 100
results<-list()
comb<-expand.grid(N=c(50,100),M=c(20))
count<-0
for (j in 1:nrow(comb)) {
  N<-comb[j,]$N
  M<-comb[j,]$M
  for (i in 1:n.sim) { # use parallel computing if possible
    count<-count+1

    model.sim<-fsem(eta1~~z1+z2,effectType="concurrent")%+%
      fsem(eta1~~z3,effectType="historical")%+%
      fsem(eta2~~z4+z5+z6,effectType="concurrent")%+%
      fsem(eta1~-1+eta2,effectType="concurrent",latent.covariate = "eta2",scalar.covariate = FALSE) %+% 
      fsem(eta2~-1)
    

    #well specified
    model.fit.w <- model.sim
    
    #misspecified
    model.fit.m<-fsem(eta1~~z1+z2,effectType="concurrent")%+%
      fsem(eta1~~z3,effectType="historical")%+%
      fsem(eta2~~z4+z5+z6,effectType="concurrent")%+%
      fsem(eta1~-1,effectType="concurrent",latent.covariate = "eta2",scalar.covariate = FALSE) %+% 
      fsem(eta2~-1)
    
    
    simm<-simulation(model=model.sim,design="regular",n.t=M,n.b=6,r=1,rho = 0.3,SNR=4,
                     n.sample = N,Matern.sem = TRUE,Matern.fac = FALSE,parameters = "parameter_set1") 
    initial.parameter<-initial.param(model = model.sim,data=simm$data,n.b=6)
    params.estimated.w<-em.estimation(model=model.fit.w,data=simm$data,
                                    initial.parameter=initial.parameter,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="regular")
    
    params.estimated.m<-em.estimation(model=model.fit.m,data=simm$data,
                                    initial.parameter=initial.parameter,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="regular")
    
    
    # single_run = list(sim.data=simm,params.estimated.well.specified=params.estimated.w, params.estimated.misspecified = params.estimated.m)
    # saveRDS(single_run, "single_run_results.rds")
    
    results[[count]]<-list(
      simulation = simm,
      model.fit.w = model.fit.w,
      model.fit.m = model.fit.m,
      model.sim = model.sim,
      estimation.w = params.estimated.w,
      estimation.m = params.estimated.m
    )
  }  
}


#extracting MSE values
source("Core/MSE_report.R")
print(tabl.fac.w)
print(tabl.sem.w)
print(tabl.fac.m)
print(tabl.sem.m)


