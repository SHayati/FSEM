
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
      fsem(eta1~-1) %+% 
      fsem(eta2~~z4+z5+z6,effectType="concurrent")%+%
      fsem(eta2~-1)
    
    model.fit<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
      fsem(eta1~~z2,effectType="concurrent")%+%
      fsem(eta1~~z3,effectType="historical")%+%
      fsem(eta1~-1)%+%
      fsem(eta2~~z4,effectType="fixed_concurrent")%+%
      fsem(eta2~~z5+z6,effectType="concurrent")%+%
      fsem(eta2~-1)
    
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


