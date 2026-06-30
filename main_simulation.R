
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
comb<-expand.grid(N=c(50,100),M=c(10,20))
count<-0
for (j in 1:nrow(comb)) {
  N<-comb[j,]$N
  M<-comb[j,]$M
  for (i in 1:n.sim) { # use parallel computing if possible
    count<-count+1
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
    params.estimated<-em.estimation(model=model.fit,data=simm$data,
                                    initial.parameter=initial.parameter,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="regular")
    results[[count]]<-list(simulation=simm,model.fit=model.fit,model.sim=model.sim,estimation=params.estimated)
  }  
}

#extracting MSE values
source("Core/MSE_tab1.R")
print(tabl.fac)
print(tabl.sem)

#Table 2:regular design
##extracting coverage rates
source("Core/CovRate_tab2.R")
print(cov.table.fac)
print(cov.table.sem)

#Table 1:irregular design
##fitting procedure
###set n.sim (number of monte carlo simulations)
results<-list()
comb<-expand.grid(N=c(50,100),M=c(10,20))
count<-0
for (j in 1:nrow(comb)) {
  N<-comb[j,]$N
  M<-comb[j,]$M
for (i in 1:n.sim) { # use parallel computing if possible
  count<-count+1
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
  params.estimated<-em.estimation(model=model.fit,data=simm$data,
                                  initial.parameter=initial.parameter,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="irregular")
  results[[count]]<-list(simulation=simm,model.fit=model.fit,model.sim=model.sim,estimation=params.estimated)
 }
}

#Table 1: extracting MSE values
source("Core/MSE_tab1.R")
print(tabl.fac)
print(tabl.sem)

#Table 2:irregular design
##extracting coverage rates
source("Core/CovRate_tab2.R")
print(cov.table.fac)
print(cov.table.sem)

#Table S1:regular missing at random design for N=100 and M=8
##fitting procedure
###set n.sim (number of monte carlo simulations)
results<-list()
for (i in 1:n.sim) { # use parallel computing if possible
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
  params.estimated<-em.estimation(model=model.fit,data=simm$data,x.data=simm$covariate,
                                  initial.parameter=initial.parameter,n.b=6,n.em=2,n.monte=2,s.p="min",range.min=NULL,range.max=NULL,design="regular.missing")
  results[[i]]<-list(simulation=simm,model.fit=model.fit,model.sim=model.sim,estimation=params.estimated)
}

#extracting MSE values for Table S1:
source("Core/MSE_tab3.R")
print(tabl.fac)
print(tabl.sem)

##extracting coverage rates for Table S1:
source("Core/CovRate_tab3.R")
print(cov.table.fac)
print(cov.table.sem)

