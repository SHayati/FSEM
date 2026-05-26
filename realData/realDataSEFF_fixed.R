library(dplyr)
library(fda)
library(Matrix)
library(mvtnorm)
library(tmvtnorm)
library(MASS)

dir.path = "realData"

source("Core/FSEM.R")
source("Core/simualtion_utilities.R")


n<-400
d.f.se = readRDS(file.path(dir.path, paste0("d.f.se_n",n,".rds")))

model.dat.se<-fsem(eta~~z1+z2+z3+z4+z5,effectType = "fixed_concurrent")%+%
  fsem(eta~-1+gender+cancer+diabet+heartcon+genCancer+genDiabet+genHeart,effectType="fixed_concurrent",scalar.covariate = TRUE)
  


initial.parameter.se<-initial.param(model = model.dat.se,data=d.f.se$data,n.b=6)
params.estimated.se<-em.estimation(model=model.dat.se,data=d.f.se$data,x.data=d.f.se$covariate,
                                   initial.parameter=initial.parameter.se,n.b=6,n.em=100,n.monte=100,s.p="min",range.min=NULL,range.max=NULL,design="regular.missing")
results<-params.estimated.se
