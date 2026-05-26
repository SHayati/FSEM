
#initial libraries
# source("realData/realDataSEFF.R")
source("realData/realDataSEFF_fixed.R")
source("realData/covRate.R")

t<-seq(2008,2022,length.out=200)
tt<-seq(2008,2022,length.out=100)
par(mfrow=c(1,1))
res<-results
params.estimated.se<-res  

#Figure 3
##factor loadings
par(mgp = c(0, 2, 0))
plot(tt,-fregion.bands.fac[[1]][[1]][,1],col="blue","l", ylim=c(-0.3,0),
     xlab="",
     #ylab=expression(lambda[1](t)),
     ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.fac[[1]][[1]][,3] , rev(-fregion.bands.fac[[1]][[1]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

plot(tt,-fregion.bands.fac[[1]][[2]][,1],col="blue","l", ylim=c(0,0.8),
     xlab="",
     #ylab=expression(lambda[1](t)),
     ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.fac[[1]][[2]][,3] , rev(-fregion.bands.fac[[1]][[2]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

plot(tt,-fregion.bands.fac[[1]][[3]][,1],col="blue","l", ylim=c(-0.1,0.8),
     xlab="",
     #ylab=expression(lambda[1](t)),
     ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.fac[[1]][[3]][,3] , rev(-fregion.bands.fac[[1]][[3]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

plot(tt,-fregion.bands.fac[[1]][[4]][,1],col="blue","l", ylim=c(-0.1,0.8),
     xlab="",
     #ylab=expression(lambda[1](t)),
     ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.fac[[1]][[4]][,3] , rev(-fregion.bands.fac[[1]][[4]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

plot(tt,-fregion.bands.fac[[1]][[5]][,1],col="blue","l", ylim=c(-0.5,0.5),
     xlab="",
     #ylab=expression(lambda[1](t)),
     ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.fac[[1]][[5]][,3] , rev(-fregion.bands.fac[[1]][[5]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

#Figure 4
##tolerances plots
par(mgp = c(0, 2, 0))
plot(1:length(params.estimated.se$tolerance$tol_sig[-1]),params.estimated.se$tolerance$tol_sig[-1],"l",
     #ylab=expression(tol[sigma^2]),
     ylab="",xlab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")

par(mgp = c(0, 2, 0))
plot(1:length(params.estimated.se$tolerance$tol[-1]),params.estimated.se$tolerance$tol[-1],"l",
     #ylab=expression(tol[coef]),
     ylab="",xlab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")

#Figure 5
##regression coefficients
par(mgp = c(0, 2, 0))
plot(tt,-fregion.bands.sem[[1]][[1]][[1]][,1],col="blue","l",ylim=c(-2.5,1.5),
     xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.sem[[1]][[1]][[1]][,3] , rev(-fregion.bands.sem[[1]][[1]][[1]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

par(mgp = c(0, 2, 0))
plot(tt,-fregion.bands.sem[[1]][[1]][[2]][,1],col="blue","l",ylim=c(-2.5,1.5),
     xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.sem[[1]][[1]][[2]][,3] , rev(-fregion.bands.sem[[1]][[1]][[2]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

par(mgp = c(0, 2, 0))
plot(tt,-fregion.bands.sem[[1]][[1]][[3]][,1],col="blue","l",ylim=c(-2.5,1.5),
     xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.sem[[1]][[1]][[3]][,3] , rev(-fregion.bands.sem[[1]][[1]][[3]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

plot(tt,-fregion.bands.sem[[1]][[1]][[4]][,1],col="blue","l",ylim=c(-2.5,1.5),
     xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.sem[[1]][[1]][[4]][,3] , rev(-fregion.bands.sem[[1]][[1]][[4]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

plot(tt,-fregion.bands.sem[[1]][[1]][[5]][,1],col="blue","l",ylim=c(-2.5,2.5),
     xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.sem[[1]][[1]][[5]][,3] , rev(-fregion.bands.sem[[1]][[1]][[5]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

plot(tt,-fregion.bands.sem[[1]][[1]][[6]][,1],col="blue","l",ylim=c(-2.5,1.5),
     xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.sem[[1]][[1]][[6]][,3] , rev(-fregion.bands.sem[[1]][[1]][[6]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)
plot(tt,-fregion.bands.sem[[1]][[1]][[7]][,1],col="blue","l",ylim=c(-2.5,1.5),
     xlab="",ylab="",cex.axis=3,cex.lab=3,lwd=3,family="serif")
polygon(c(tt, rev(tt)), c(-fregion.bands.sem[[1]][[1]][[7]][,3] , rev(-fregion.bands.sem[[1]][[1]][[7]][,2])), col=rgb(0.1 , 0.1 , 1 , alpha=0.1), border=NA)

#Figure 6 and Table 4
##goodness of fit indices
source("realData/goodnessOfFit.R")
