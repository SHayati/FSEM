
#setwd("D:\\FSEM2\\FSEM\\FSEM2")
setwd("D:\\FSEM2\\FSEM\\FSEM2\\copies_final\\copies_final")


debugSource("utilities03022023.3.1 - Copy.R")           #Utilities
debugSource("FSEM24012023 - Copy - Copy - Copy.2 - cross3.R")    #Fitting procedure

parameter<-parameter1
model1.sim<-fsem(eta1~~z1+z2,effectType="concurrent")%+%
  fsem(eta1~~z3,effectType="concurrent")%+%
  fsem(eta1~-1)


model1<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
  fsem(eta1~~z2,effectType="concurrent")%+%
  fsem(eta1~~z3,effectType="concurrent")%+%
  fsem(eta1~-1)

#parameter<-parameter2
model2.sim<-fsem(eta1~~z1+z2+z3+z4,effectType="concurrent")%+%
  fsem(eta2~~z5+z6+z7+z8,effectType="concurrent")%+%
  fsem(eta3~~z9+z10+z11+z12,effectType="concurrent")%+%
  fsem(eta3~-1+eta1,effectType="concurrent",latent.covariate = "eta1",scalar.covariate = FALSE)%+%
  fsem(eta3~-1+eta2,effectType="concurrent",latent.covariate = "eta2",scalar.covariate = FALSE)%+%
  fsem(eta1~-1+eta2,effectType="concurrent",latent.covariate = "eta2",scalar.covariate = FALSE)%+%
  fsem(eta2~-1)

model2<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
  fsem(eta1~~z2+z3+z4,effectType="concurrent")%+%
  fsem(eta2~~z5,effectType="fixed_concurrent")%+%
  fsem(eta2~~z6+z7+z8,effectType="concurrent")%+%
  fsem(eta3~~z9,effectType="fixed_concurrent")%+%
  fsem(eta3~~z10+z11+z12,effectType="concurrent")%+%
  fsem(eta3~-1+eta1,effectType="concurrent",latent.covariate = "eta1",scalar.covariate = FALSE)%+%
  fsem(eta3~-1+eta2,effectType="concurrent",latent.covariate = "eta2",scalar.covariate = FALSE)%+%
  fsem(eta1~-1+eta2,effectType="concurrent",latent.covariate = "eta2",scalar.covariate = FALSE)%+%
  fsem(eta2~-1)

# model3<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
#   fsem(eta1~~z2+z3,effectType="concurrent")%+%
#   fsem(eta2~~z4,effectType="fixed_concurrent")%+%
#   fsem(eta2~~z5+z6,effectType="concurrent")%+%
#   fsem(eta2~-1+eta1,effectType="concurrent",latent.covariate = "eta1",scalar.covariate = FALSE)%+%
#   fsem(eta1~-1)
# 
# model4<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
#   fsem(eta1~~z2+z3,effectType="concurrent")%+%
#   fsem(eta2~~z4,effectType="fixed_concurrent")%+%
#   fsem(eta2~~z5+z6,effectType="concurrent")%+%
#   fsem(eta2~-1)%+%
#   fsem(eta1~-1)
# 
# model5<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
#   fsem(eta1~~z2+z3,effectType="concurrent")%+%
#   fsem(eta2~~z4,effectType="fixed_concurrent")%+%
#   fsem(eta2~~z5+z6,effectType="concurrent")%+%
#   fsem(eta2~-1+x,effectType="smooth",scalar.covariate=TRUE)%+%
#   fsem(eta2~-1+eta1,effectType="concurrent",latent.covariate = "eta1",scalar.covariate = FALSE)%+%
#   fsem(eta1~-1)
# 
# model6<-fsem(eta1~~z1,effectType="fixed_concurrent")%+%
#   fsem(eta1~~z2,effectType="historical")%+%
#   fsem(eta1~~z3,effectType="concurrent")
# fsem(eta1~-1)

simm<-simulation(model=model1.sim,design="regular",n.t=20,n.b.sim=100,r=1.5,rho = 0.1,SNR=10,
                 sample = 50,Matern.sem = TRUE,Matern.fac = FALSE) 
initial.parameter<-initial.param(model = model1,data=simm$data,n.b=12,n.b.sim=100,r=1.5,rho=0.1,SNR=10,Matern.fac=FALSE,Matern.sem=TRUE)

params.estimated<-em.estimation(model=model1,model.sim=model1.sim,param=initial.parameter,
                                range.min=NULL,range.max=NULL,
                                n.em=100,n.monte=100,
                                n.b=12,n.b.sim=100,
                                data=simm$data,s.p=NULL,
                                Matern.sem=TRUE,Matern.fac=FALSE,
                                r=1.5,rho=0.1,SNR=10,
                                x.data =simm$covariate)
params.org.eval<-parameter.org.evaluted(model=model1.sim,r=1.5,rho=0.1,n.b.sim=100,range.min=NULL,range.max=NULL,SNR=10,x.data=simm$covariate,Matern.sem=TRUE,Matern.fac=FALSE)

# tt<-seq(0,1,length.out=length(params.estimated$result$params.estimated.eval$coef.sem$eta2$eta1))
# plot(tt,params.estimated$result$params.estimated.eval$coef.sem$eta2$eta1,"l",col="red",ylim=c(-1,1))
# lines(tt,params.org.eval$coef.sem$eta2$eta1,col="blue")http://127.0.0.1:40637/graphics/bf61f250-3245-4f24-bd18-d8ec6de85001.png
# 
