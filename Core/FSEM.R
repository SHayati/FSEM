
#########time points##################################
timeEval<-function(n.b,range.min,range.max,data){
  samples<-length(unique(data$.id))
  no.f<-length(unique(data$.ind))
  m<-list()
  M<-list()
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))
  ## Split the time points by (.id,.ind) ONCE (O(total observations)) instead
  ## of re-scanning the full .t vector with `data$.id==i & data$.ind==j` inside
  ## an O(samples*no.f) loop (which is O(N^2) in the number of subjects and was
  ## the dominant cost when timeEval is called repeatedly).  The within-group
  ## order follows the original row order, so unique() yields identical results.
  grp <- split(data$.t, list(data$.id, data$.ind), drop = FALSE, sep = "\r")
  tp <- function(i, j) { v <- grp[[paste0(i, "\r", j)]]; if (is.null(v)) numeric(0) else unique(v) }
  evall<-list()
  tt<-list()
  for (i in 1:samples) {
    tpts <- lapply(1:no.f, function(j) tp(i, j))
    M[[i]]    <- sapply(tpts, length)
    tt[[i]]   <- tpts
    evall[[i]]<- lapply(tpts, function(v) t(eval.basis(v, basis)))
  }
  m$time.no<-M
  m$eval<-evall
  m$time.point<-tt
  m
}
#########omega matrices##################################
omega.matrix<-function(n.b,range.min,range.max,data){
  eval.t<-timeEval(n.b,range.min,range.max,data)
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))
  tt<-seq(range.min,range.max,length.out=100)
  eval<-t(eval.basis(tt,basis))
  samples=length(unique(data$.id))
  n.t<-max(c(sapply(1:samples,function(i){eval.t$time.no[[i]]})))
  evall2<-t(eval.basis(seq(0,1,length.out=n.t),basis))
  no.fac<-length(unique(data$.ind))
  samples<-length(unique(data$.id))
  m<-list()
  for (i in 1:samples) {
    m$omega1[[i]]<-list() 
    for (j in 1:no.fac) {
      m$omega1[[i]][[j]]<-0
      for (k in 1:eval.t$time.no[[i]][[j]]) {
        omegaa<-eval.t$eval[[i]][[j]][,k]%*%t(eval.t$eval[[i]][[j]][,k])
        m$omega1[[i]][[j]]<-rbind(m$omega1[[i]][[j]],omegaa)
      }
      m$omega1[[i]][[j]]<-m$omega1[[i]][[j]][-1,]
    }
  }
  m$omega1.sem<-0
  for (k in 1:n.t) {
    omegaa<-evall2[,k]%*%t(evall2[,k])
    m$omega1.sem<-rbind(m$omega1.sem,omegaa)
  }
  m$omega1.sem<-m$omega1.sem[-1,]
  
  for (i in 1:samples) {
    m$omega2[[i]]<-list()
    m$omega2.gen[[i]]<-list()
    m$delta[[i]]<-list()
    for (j in 1:no.fac) {
      m$omega2[[i]][[j]]<-0
      m$omega2.gen[[i]][[j]]<-0
      m$delta[[i]][[j]]<-0
      for (k in 1:eval.t$time.no[[i]][[j]]) {
        om1<-0
        om2<-0
        om3<-0
        tt<-seq(range.min,eval.t$time.point[[i]][[j]][k],length.out=100)
        eval<-t(eval.basis(tt,basis))
        for (s in 1:100) {
          om<-eval[,s]%*%kronecker(t(eval.t$eval[[i]][[j]][,k]),t(eval[,s]))
          omm<-kronecker(eval.t$eval[[i]][[j]][,k],eval[,s])%*%t(eval[,s])
          omm2<-t(eval[,s])
          om1<-om1+1/100*om
          om2<-om2+1/100*omm
          om3<-om3+1/100*omm2
        }
        m$omega2[[i]][[j]]<-rbind(m$omega2[[i]][[j]],om1)
        m$omega2.gen[[i]][[j]]<-rbind(m$omega2.gen[[i]][[j]],om2)
        m$delta[[i]][[j]]<-rbind(m$delta[[i]][[j]],om3)
      }
      m$omega2[[i]][[j]]<-m$omega2[[i]][[j]][-1,]
      m$omega2.gen[[i]][[j]]<-m$omega2.gen[[i]][[j]][-1,]
      m$delta[[i]][[j]]<-m$delta[[i]][[j]][-1,]
    }
  }
  m$omega2.sem<-0
  m$omega2.gen.sem<-0
  for (k in 1:n.t) {
    om1<-0
    om2<-0
    se<-seq(0,seq(0,1,length.out=n.t)[k],length.out=100)
    evall2<-t(eval.basis(se,basis))
    for (s in 1:100) {
      om<-eval[,s]%*%kronecker(t(evall2[,k]),t(eval[,s]))
      omm<-kronecker(evall2[,k],eval[,s])%*%t(eval[,s])
      om1<-om1+1/100*om
      om2<-om2+1/100*omm
    }
    m$omega2.sem<-rbind(m$omega2.sem,om1)
    m$omega2.gen.sem<-rbind(m$omega2.gen.sem,om2)
  }
  m$omega2.sem<-m$omega2.sem[-1,]
  m$omega2.gen.sem<-m$omega2.gen.sem[-1,]
  m
}

#########Fitting Procedure##########################################################
weight.matrix<-function(model,n.b,sample){
  weight<-list()
  no.lat<-length(model$var$latents)
  d.eta<-matrix(0,n.b*sample*no.lat,n.b*sample*no.lat)
  for (j in 1:no.lat) {
    d<-matrix(0,n.b*sample,sample*n.b*no.lat)
    for (i in 1:sample) {
      a<-diag(n.b)
      d[((i-1)*n.b+1):(i*n.b),((i-1)*no.lat*n.b+1+n.b*(j-1)):((i-1)*no.lat*n.b+j*n.b)]<-a
    }
    d.eta[((j-1)*n.b*sample+1):(j*n.b*sample),]<-d
  }
  weight[["eta"]]<-d.eta
  
  if(!is.null(model$var$observed)){
    no.obs<-length(model$var$observed)
    n.bb<-c()
    for (j in 1:no.obs) {
      n.bb[j]<-ifelse(model$var$observedScalar[j,2]==TRUE,n.b,n.b)
    }
    d.x<-matrix(0,sum(n.bb)*sample,sum(n.bb)*sample)
    for (j in 1:no.obs) {
      d<-matrix(0,n.bb[j]*sample,sum(n.bb)*sample)
      for (i in 1:sample) {
        a<-diag(n.bb[j])
        if(j>1){
          d[((i-1)*n.bb[j]+1):(i*n.bb[j]),((i-1)*sum(n.bb)+1+n.bb[1]*(j-1)):((i-1)*sum(n.bb)+sum(n.bb[1:j]))]<-a
        }else{
          d[((i-1)*n.bb[j]+1):(i*n.bb[j]),((i-1)*sum(n.bb)+1):((i-1)*sum(n.bb)+n.bb[j])]<-a
        }
      }
      if(j>1){
        d.x[(n.bb[1]*(j-1)*sample+1):(sum(n.bb[1:j])*sample),]<-d  
      }else{
        d.x[(1:(n.bb[j]*sample)),]<-d
      }
    }  
    weight[["x"]]<-d.x 
  }
  weight
}


################test##################################
# n.b<-10
# no.lat<-length(model$var$latents)
# sample<-30
# weightt<-weight.matrix(model,n.b,sample)
# ett<-lapply(1:no.lat,function(i){matrix(runif(sample*n.b),sample,n.b)})
# eta<-NULL
# eta.star<-NULL
# for (i in 1:no.lat) {
#   et<-NULL
#   for (j in 1:sample) {
#     et<-c(et,ett[[i]][j,])
#   }
#   eta<-c(eta,et)
# }
# for (j in 1:sample) {
#   et<-NULL
#   for (i in 1:no.lat) {
#     et<-c(et,ett[[i]][j,])
#   }
#   eta.star<-c(eta.star,et)
# }
# all(eta==weightt$eta%*%eta.star)
# 
# no.obs<-length(model$var$observed)
# n.bb<-c()
# for (j in 1:no.obs) {
#   n.bb[j]<-ifelse(model$var$observedScalar[j,2]==TRUE,1,n.b)
# }
# xx<-lapply(1:no.obs,function(i){matrix(runif(sample*n.bb[i]),sample,n.bb[i])})
# x<-NULL
# x.star<-NULL
# for (i in 1:no.obs) {
#   et<-NULL
#   for (j in 1:sample) {
#     et<-c(et,xx[[i]][j,])
#   }
#   x<-c(x,et)
# }
# for (j in 1:sample) {
#   et<-NULL
#   for (i in 1:no.obs) {
#     et<-c(et,xx[[i]][j,])
#   }
#   x.star<-c(x.star,et)
# }
# all(x==weightt$x%*%x.star)

################################################################################
diag.matrix<-function(model,n.b){
  no.f<-length(model$mod$factorModel)
  no.r<-length(model$mod$regression)
  dd<-list()
  for (j in 1:no.f) {
    d<-0
    a<-model$var$latents%in%model$mod$factorModel[[j]]$factor
    for (k in 1:length(a)) {
      d<-Matrix::bdiag(d,a[k]*diag(n.b))
    }
    d<-d[-1,-1]
    dd[[paste0("f",j)]]<-d
  }
  if(no.r>0)
    for (j in 1:no.r) {
      d<-0
      d1<-0
      a<-model$var$factors%in%model$mod$regression[[j]]$covariate
      for (k in 1:length(a)) {
        d<-Matrix::bdiag(d,a[k]*diag(n.b))
      }
      if(length(model$var$observed)>0){
        a1<-model$var$observed%in%model$mod$regression[[j]]$covariate
        for (k in 1:length(a1)) {
          if(isTRUE(a1[k])&&model$mod$regression[[j]]$effect[which(model$mod$regression[[j]]$covariate==model$var$observed[k])]=="linear"){
            d1<-Matrix::bdiag(d1,a1[k])  
          }else{
            d1<-Matrix::bdiag(d1,a1[k]*diag(n.b))
          }
        }
        d1<-d1[-1,-1]
        dd[[paste0("r",j,".x")]]<-d1
      }
      d<-d[-1,-1]
      dd[[paste0("r",j,".eta")]]<-d
    }
  dd
}


diag.matrix2<-function(model,n.b){
  no.f<-length(model$mod$factorModel)
  no.r<-length(model$mod$regression)
  dd<-list()
  for (j in 1:no.f) {
    d<-diag(n.b)
    a<-model$var$factors%in%model$mod$factorModel[[j]]$factor
    ww<-which(a)
    for (k in 1:length(a)) {
      if(length(ww)>0){
        if(k%in%ww){
          if(model$mod$factorModel[[j]]$effect=="historical"){
            d<-Matrix::bdiag(d,a[k]*diag(n.b*n.b))  
          }
          if(model$mod$factorModel[[j]]$effect=="concurrent"){
            d<-Matrix::bdiag(d,a[k]*diag(n.b)) 
          }
          if(model$mod$factorModel[[j]]$effect=="fixed_concurrent"){
            d<-Matrix::bdiag(d,0*diag(n.b)) 
          }
          if(model$mod$factorModel[[j]]$effect=="fixed_historical"){
            d<-Matrix::bdiag(d,0*diag(n.b)) 
          }
        }else{
          d<-Matrix::bdiag(d,a[k]*diag(n.b))
        }
      }else{
        d<-Matrix::bdiag(d,a[k]*diag(n.b))
      }
    }
    dd[[paste0("f",j)]]<-d 
  }
  regs<-model$var$latents
  no.reg<-which(sapply(1:no.r,function(i2){model$mod$regression[[i2]]$response%in%regs}))
  for (j in no.reg) {
    d<-0
    d1<-0
    a<-model$var$latents%in%model$mod$regression[[j]]$covariate
    ww<-which(a)
    for (k in 1:length(a)) {
      if(length(ww)>0){
        if(k%in%ww){
          if(model$mod$regression[[j]]$effect[model$mod$regression[[j]]$covariate==model$var$latents[k]]=="historical"){
            d<-Matrix::bdiag(d,a[k]*diag(n.b*n.b))  
          }
          if(model$mod$regression[[j]]$effect[model$mod$regression[[j]]$covariate==model$var$latents[k]]=="fixed"){
            d<-Matrix::bdiag(d,0*diag(n.b))
          }
          if(model$mod$regression[[j]]$effect[model$mod$regression[[j]]$covariate==model$var$latents[k]]=="concurrent"){
            d<-Matrix::bdiag(d,a[k]*diag(n.b))
          } 
        }else{
          d<-Matrix::bdiag(d,a[k]*diag(n.b))
        }
      }else{
        d<-Matrix::bdiag(d,a[k]*diag(n.b))
      }
    }
    if(!is.null(model$var$observed)){
      a1<-model$var$observed%in%model$mod$regression[[j]]$covariate
      ww1<-which(a1)
      for (k in 1:length(a1)) {
        if(length(ww)>0){
          if(k%in%ww1){
            if(model$mod$regression[[j]]$effect[model$mod$regression[[j]]$covariate==model$var$observed[k]]%in%c("smooth","historical")){
              d1<-Matrix::bdiag(d1,a1[k]*diag(n.b*n.b))  
            }else{
              d1<-Matrix::bdiag(d1,a1[k]*diag(n.b))
            }  
          }else{
            d1<-Matrix::bdiag(d1,a1[k]*diag(n.b))
          }
        }else{
          d1<-Matrix::bdiag(d1,a1[k]*diag(n.b))
        }
      }
      d1<-d1[-1,-1]
      dd[[paste0("r",j,".x")]]<-d1
    }
    d<-d[-1,-1]
    dd[[paste0("r",j,".eta")]]<-d
  } 
  
  dd
}

#diagMatrix<-diag.matrix2(model,n.b)
##################################################################################
###############generator distributions#############################################

initial.param<-function(model,n.b,data){
  
  range.min = 0
  range.max = 1
  Posdef <- function (n, ev)
  {
    Z <- matrix(ncol=n, runif(n^2,10,11))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
  }
  no.r<-length(model$mod$regression)
  no.obs<-length(model$var$observed)
  no.f<-length(model$mod$factorModel)
  init<-list()
  lat2<-model$var$latents
  obs2<-model$var$observed
  init$gamma.param<-list()
  init$sigma.sem<-list()
  count<-0
  samples<-length(unique(data$.id))
  fr<-length(model$mod$factorModel)
  no.lat<-length(lat2)
  var<-list()
  for (i in 1:fr) {
    lambda.fac<-NULL
    init$lambda.param[[i]]<-list()
    init$sigma.fac[[i]]<-list()
    init$sigma.error[[i]]<-0.5*var(data$.value[data$.id==i])
    if(model$mod$factorModel[[i]]$intercept==TRUE){
      intercept<-rep(0,n.b) #fourier coefficients of intercepts
    }
    for (j in 1:no.lat) {
      fac<-which(model$mod$factorModel[[i]]$factor==model$var$latents[j])
      if(length(fac)==0){
        lambda<-rep(0,n.b)   #fourier coefficients of factor loadings
      }else{
        if(model$mod$factorModel[[i]]$effect[fac]=="concurrent"){
          lambda<-rep(0,n.b)
        }
        if(model$mod$factorModel[[i]]$effect[fac]=="historical"){
          lambda<-rep(0,n.b*n.b)
        }
        if(model$mod$factorModel[[i]]$effect[fac]=="fixed_concurrent"){
          lambda<-rep(0,n.b)
        }
        if(model$mod$factorModel[[i]]$effect[fac]=="fixed_historical"){
          lambda<-rep(0,n.b)
        }
      }
      lambda.fac<-c(lambda.fac,lambda)
    }
    init$lambda.param[[i]]<-c(intercept,lambda.fac)
    
    init$sigma.fac[[i]]<-diag(n.b)   #the matrix sigma in factor model
  }
  for (i in 1:no.r) {
    gamma.eta<-NULL
    gamma.x<-NULL
    count<-which(lat2==model$mod$regression[[i]]$response)
    for (j in 1:length(lat2)){
      fac<-which(model$mod$regression[[i]]$covariate==lat2[j])
      if(length(fac)==0){
        gam.eta<-rep(0,n.b)
      }else{
        if(model$mod$regression[[i]]$effect[fac]=="concurrent"){
          gam.eta<-rep(0,n.b)  #fourier coefficients of coefficients of latents in sem
        }
        if(model$mod$regression[[i]]$effect[fac]=="historical"){
          gam.eta<-rep(0,n.b*n.b)  #fourier coefficients of coefficients of latents in sem
        }
        if(model$mod$regression[[i]]$effect[fac]=="fixed"){
          gam.eta<-rep(0,n.b)
        }
      }
      gamma.eta<-c(gamma.eta,gam.eta)
    }
    if(no.obs>0){
      for (j in 1:length(obs2)) {
        fac<-which(model$mod$regression[[i]]$covariate==obs2[j])
        if(length(fac)==0){
          gam.x<-rep(0,n.b)
        }else{
          if(model$mod$regression[[i]]$effect[fac]=="concurrent"){
            gam.x<-rep(0,n.b)  #fourier coefficients of coefficients of observed vs in sem
          }
          if(model$mod$regression[[i]]$effect[fac]=="historical"){
            gam.x<-rep(0,n.b*n.b)  #fourier coefficients of coefficients of observed vs in sem
          }
          if(model$mod$regression[[i]]$effect[fac]=="smooth"){
            gam.x<-rep(0,n.b*n.b)      #fourier coefficients of coefficients of observed vs in sem
          }
          if(model$mod$regression[[i]]$effect[fac]=="constant"){
            gam.x<-rep(0,n.b)   #fourier coefficients of coefficients of observed vs in sem
          }
          if(model$mod$regression[[i]]$effect[fac]=="linear"){
            gam.x<-rep(0,n.b)   #fourier coefficients of coefficients of observed vs in sem
          }
        }
        gamma.x<-c(gamma.x,gam.x)
      }
    }
    init$gamma.param[[count]]<-if(no.obs>0){c(gamma.eta,gamma.x)}else{gamma.eta}
    init$sigma.sem[[count]]<-diag(n.b)   #the matrix sigma in sem
  }
  init
}


eta.distribution<-function(model,param,n.b,range.min,range.max,data,x.data=NULL,omegaMatrix=NULL){
  eval.t<-timeEval(n.b,range.min,range.max,data)
  diagMatrix<-diag.matrix(model,n.b)
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))
  samples<-length(unique(data$.id))
  n.t<-max(c(sapply(1:samples,function(i){eval.t$time.no[[i]]})))
  evall2<-t(eval.basis(seq(0,1,length.out=n.t),basis))
  soll<-ginv(evall2%*%
               t(evall2))%*%evall2
  weightMatrix<-weight.matrix(model,n.b,sample=samples)
  extra<-list()
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))  
  ome<-inprod(basis,basis)
  ##################eta distribution####################################
  no.r<-length(model$mod$regression)
  no.obs<-length(model$var$observed)
  no.f<-length(model$mod$factorModel)
  no.obs2<-model$var$observed
  lat<-model$var$latents
  lat2<-lat
  obs2<-no.obs2
  b<-0
  upsilon<-0
  b2<-0
  upsilon2<-list()
  for (is in 1:samples) {
    upsilon2[[is]]<-0 
  }
  covv<-0
  covv.2<-0
  extra$mu.et<-list()
  extra$mu.cov<-list()
  for (ii in 1:length(lat2)) {
    r.eta1<-0
    count2<-0
    i<-which(sapply(1:no.r,function(kk){model$mod$regression[[kk]]$response==lat2[ii]}))
    count<-which(lat==model$mod$regression[[i]]$response)
    for (j in 1:length(lat2)){
      fac<-which(model$mod$regression[[i]]$covariate==lat2[j])
      if(length(fac)==0){
        omega<-matrix(0,n.t*n.b,n.b)
        gam.eta<-param$gamma.param[[count]][(count2+1):(count2+n.b)]
        count2<-count2+n.b
      }else{
        if(model$mod$regression[[i]]$effect[fac]=="concurrent"){
          gam.eta<-param$gamma.param[[count]][(count2+1):(count2+n.b)]  #fourier coefficients of coefficients of latents in sem
          omega<-omegaMatrix$omega1.sem
          count2<-count2+n.b
        }
        if(model$mod$regression[[i]]$effect[fac]=="historical"){
          gam.eta<-param$gamma.param[[count]][(count2+1):(count2+n.b*n.b)]  #fourier coefficients of coefficients of latents in sem
          omega<-omegaMatrix$omega2.gen.sem
          count2<-count2+n.b*n.b
        }
      }
      r.eta1<-cbind(r.eta1,t(kronecker(diag(n.t),gam.eta))%*%omega)
    }
    r.eta<-soll%*%r.eta1[,-1]
    extra$mu.et[[ii]]<-r.eta
    r.x1<-0
    if(length(obs2)>0){
      xx<-list()
      for (j in 1:length(obs2)) {
        fac<-which(model$mod$regression[[i]]$covariate==obs2[j])
        if(length(fac)==0){
          omega<-matrix(0,n.t*n.b,n.b)
          gam.x<-param$gamma.param[[count]][(count2+1):(count2+n.b)]
          r.xx<-t(kronecker(diag(n.t),gam.x))%*%omega
          count2<-count2+n.b
        }else{
          if(model$mod$regression[[i]]$effect[fac]=="concurrent"){
            gam.x<-param$gamma.param[[count]][(count2+1):(count2+n.b)]  #fourier coefficients of coefficients of observed vs in sem
            omega<-omegaMatrix$omega1.sem
            r.xx<-t(kronecker(diag(n.t),gam.x))%*%omega
            count2<-count2+n.b
            
            sm<-matrix(NA,n.b,samples)
            for (tt in 1:samples) {
              sm[,tt]<-Data2fd(argvals=seq(range.min,range.max,length.out=length(x.data[[no.obs2[j]]][[tt]])),y=x.data[[no.obs2[j]]][[tt]],basisobj = basis,lambda = 0.5)$coefs
            }
            xx[[j]]<-sm
          }
          if(model$mod$regression[[i]]$effect[fac]=="historical"){
            gam.x<-param$gamma.param[[count]][(count2+1):(count2+n.b*n.b)]  #fourier coefficients of coefficients of observed vs in sem
            omega<-omegaMatrix$omega2.gen.sem
            r.xx<-t(kronecker(diag(n.t),gam.x))%*%omega
            count2<-count2+n.b*n.b
            
            sm<-matrix(NA,n.b,samples)
            for (tt in 1:samples) {
              sm[,tt]<-Data2fd(argvals=seq(range.min,range.max,length.out=length(x.data[[no.obs2[j]]][[tt]])),y=x.data[[no.obs2[j]]][[tt]],basisobj = basis,lambda = 0.5)$coefs
            }
            xx[[j]]<-sm
          }
          if(model$mod$regression[[i]]$effect[fac]=="smooth"){
            Delta<-matrix(param$gamma.param[[count]][(count2+1):(count2+n.b*n.b)],n.b,n.b) #fourier coefficients of coefficients of observed vs in sem
            r.xx<-t(evall2)%*%Delta
            count2<-count2+n.b*n.b
            
            basis3<-create.bspline.basis(rangeval = c(min(x.data[[no.obs2[j]]]),max(x.data[[no.obs2[j]]])),nbasis=n.b)
            smm<-t(eval.basis(x.data[[no.obs2[j]]],basis3))
            xx[[j]]<-smm
          }
          if(model$mod$regression[[i]]$effect[fac]=="linear"){
            delta<-param$gamma.param[[count]][(count2+1):(count2+n.b)] #fourier coefficients of coefficients of observed vs in sem
            r.xx<-t(evall2)%*%delta
            count2<-count2+n.b
            
            xx[[j]]<-x.data[[no.obs2[j]]]
          }
        }
        r.x1<-cbind(r.x1,r.xx)
      }
      x.data12<-list()
      for (is in 1:samples) {
        x.data11<-NULL
        for (j in 1:length(obs2)) {
          if(is.matrix(xx[[j]])){
            x.data11<-c(x.data11,xx[[j]][,is]) 
          }else{
            x.data11<-c(x.data11,xx[[j]][is]) 
          }
        }
        x.data12[[is]]<-x.data11
      }
      x.data13<-NULL
      for (is in 1:samples) {
        x.data13<-c(x.data13,x.data12[[is]]) 
      }
      r.x1<-soll%*%r.x1[,-1]
      ## full-sample kronecker(diag(samples),...) build removed: it fed only the
      ## dead extra$eta$mean/covariance. Per-subject upsilon2 is kept below.
      for (is in 1:samples) {
        upsilon2[[is]]<-rbind(upsilon2[[is]],(r.x1%*%diagMatrix[[paste0("r",i,".x")]]%*%x.data12[[is]]))
      }
      extra$mu.cov[[ii]]<-r.x1
    }else{
      extra$mu.cov[[ii]]<-0
    }
    ## full-sample b/covv kronecker builds removed (fed only dead big mean/cov);
    ## per-subject b2 and covv.2 are retained for the per-subject E-step.
    b2<-rbind(b2,(r.eta%*%diagMatrix[[paste0("r",i,".eta")]]))
    sig.lm<-param$sigma.sem[[count]] #the matrix sigma in sem
    cov.u2<-sig.lm
    covv.2<-Matrix::bdiag(covv.2,cov.u2)
  }
  
  if(length(obs2)>0){
    for (is in 1:samples) {
      upsilon2[[is]]<-upsilon2[[is]][-1,]
    }
    uu.2<-upsilon2
  }
  un<-length(unique(model$var$latents))
  n<-n.b*un*samples
  n2<-n.b*un
  sol2<-solve(diag(n2)-b2[-1,])
  ## NOTE: the full (n.b * q * samples)-dimensional mean and covariance below
  ## are never consumed downstream -- params.z1() uses only the per-subject
  ## extra$mu.eta[[i]] and the per-subject extra$cov.eta (both built from sol2).
  ## Forming the dense (I - B)^{-1} of dimension n.b*q*samples is an O(N^3)
  ## operation that dominates the E-step for large N, so it is skipped here.
  extra$eta$mean<-NULL
  extra$eta$covariance<-NULL
  extra$mu.eta<-list()
  for (i in 1:samples) {
    extra$mu.eta[[i]]<-if(length(obs2)>0){sol2%*%uu.2[[i]]}else{rep(0,n2)}
  }
  extra$cov.eta<-as.matrix(sol2%*%covv.2[-1,-1]%*%t(sol2))
  # extra$cov<-if(length(obs2)>0){t(t(x))}else{rep(0,n.b*samples)}
  extra
}

params.z1<-function(model,param,n.b,range.min,range.max,eta.dist,data,omegaMatrix=NULL,
                    design=c("regular", "irregular", "regular.truncated", "regular.missing")){
  eval.t<-timeEval(n.b,range.min,range.max,data)
  diagMatrix<-diag.matrix(model,n.b)
  samples<-length(unique(data$.id))
  design = match.arg(design)
  if(design=="regular"){
    design.regular=TRUE
  }else{
    design.regular=FALSE
  }
  weightMatrix<-weight.matrix(model,n.b,sample=samples)
  extra<-list()
  ##################params of z####################################
  no.f<-length(model$mod$factorModel)
  fr<-no.f
  no.lat<-length(model$var$latents)
  if(design.regular==TRUE){
    soll<-ginv(eval.t$eval[[1]][[1]]%*%
                 t(eval.t$eval[[1]][[1]]))%*%
      eval.t$eval[[1]][[1]]  
  }
  for (i in 1:samples) {
    sample.no<-i
    lambda.param<-list()
    alpha<-0
    a1<-0
    cov<-0
    for (ii in 1:fr) {
      if(design.regular==FALSE){
        soll<-ginv(eval.t$eval[[sample.no]][[ii]]%*%
                     t(eval.t$eval[[sample.no]][[ii]]))%*%
          eval.t$eval[[sample.no]][[ii]] 
      }
      
      if(model$mod$factorModel[[ii]]$intercept==TRUE){
        intercept<-param$lambda.param[[ii]][1:n.b] #fourier coefficients of intercepts
        alpha<-c(alpha,intercept) 
      }else{
        stop("Intercept must be included")
      }
      theta1<-0
      count.fac<-0
      lambda.param[[ii]]<-param$lambda.param[[ii]][-(1:n.b)]
      for (j in 1:no.lat) {
        fac<-which(model$mod$factorModel[[ii]]$factor==model$var$latents[j])
        if(length(fac)==0){
          lambda<-lambda.param[[ii]][(count.fac+1):(count.fac+n.b)]   #fourier coefficients of factor loadings
          omega<-matrix(0,eval.t$time.no[[sample.no]][[ii]]*n.b,n.b)
          theta<-soll%*%
            t(kronecker(diag(eval.t$time.no[[sample.no]][[ii]]),lambda))%*%omega
          count.fac<-count.fac+n.b
        }else{
          if(model$mod$factorModel[[ii]]$effect[fac]=="concurrent"){
            lambda<-lambda.param[[ii]][(count.fac+1):(count.fac+n.b)]   #fourier coefficients of factor loadings
            omega<-omegaMatrix$omega1[[sample.no]][[ii]]
            theta<-soll%*%
              t(kronecker(diag(eval.t$time.no[[sample.no]][[ii]]),lambda))%*%omega
            count.fac<-count.fac+n.b
          }
          if(model$mod$factorModel[[ii]]$effect[fac]=="historical"){
            lambda<-lambda.param[[ii]][(count.fac+1):(count.fac+n.b*n.b)]   #fourier coefficients of factor loadings
            omega<-omegaMatrix$omega2.gen[[sample.no]][[ii]] 
            theta<-soll%*%
              t(kronecker(diag(eval.t$time.no[[sample.no]][[ii]]),lambda))%*%omega
            count.fac<-count.fac+n.b*n.b
          }
          if(model$mod$factorModel[[ii]]$effect[fac]=="fixed_concurrent"){
            theta<-diag(n.b)
            count.fac<-count.fac+n.b
          }
          if(model$mod$factorModel[[ii]]$effect[fac]=="fixed_historical"){
            theta<-soll%*%
              omegaMatrix$delta[[sample.no]][[ii]] 
            count.fac<-count.fac+n.b
          }
        }
        theta1<-cbind(theta1,theta) 
      }
      theta2<-theta1[,-1]  #the matrix of theta in (22)
      a1<-rbind(a1,theta2%*%diagMatrix[[paste0("f",ii)]])
      sig.j<-param$sigma.fac[[ii]]  #the matrix sigma in (22)
      cov.u<-sig.j   
      cov<-Matrix::bdiag(cov,cov.u) #the covariance of u in z=alpha+AD.eta+u
    }
    a<-a1[-1,]  #the matrix of A in z=alpha+AD.eta+u
    alpha1<-alpha[-1]  #the matrix of alpha in z=alpha+AD.eta+u
    cova<-cov[-1,-1] #the covariance of u in z=alpha+AD.eta+u
    
    extra[[i]]<-list()
    samp<-i
    m.s<-sapply(1:fr,function(ii){sum(sapply(samp, function(ii2){eval.t$time.no[[ii2]][ii]}))})
    err<-0
    ev<-0
    for (ii in 1:fr) {
      err<-Matrix::bdiag(err,param$sigma.error[[ii]]*diag(m.s[ii]))
      ev<-Matrix::bdiag(ev,t(eval.t$eval[[samp]][[ii]]))
    }
    err2<-err[-1,-1]
    ev<-ev[-1,-1]
    weight<-weightMatrix$eta[(no.lat*n.b*(samp-1)+1):(no.lat*n.b*samp),]
    mu.eta<-eta.dist$mu.eta[[samp]]
    sigma.eta<-eta.dist$cov.eta
    mu.z<-alpha1+a%*%mu.eta  #mean of z in z=alpha+AD.eta+u
    ee<-a%*%sigma.eta%*%t(a)
    ee<-(ee+t(ee))/2
    cov.z<-cova+ee #covariance of z in z=alpha+AD.eta+u
    eval2<-ev
    mu.w<-eval2%*%t(t(mu.z))
    cov.w<-eval2%*%cov.z%*%t(eval2)
    cov.w<-(cov.w+t(cov.w))/2
    cov.w2<-cov.w+err2
    extra[[i]]$params$a.matrix<-as.matrix(a)
    extra[[i]]$params$mean.w<-as.vector(mu.w)
    extra[[i]]$params$sigma.w<-as.matrix(cov.w2)
    extra[[i]]$params$mean.z<-as.vector(mu.z)
    extra[[i]]$params$sigma.z<-as.matrix(cov.z)
    extra[[i]]$params$mu.eta<-mu.eta
    extra[[i]]$params$sigma.eta<-sigma.eta
    extra[[i]]$params$weight<-weight
  }
  extra
}


params.z<-function(model,param,n.b,sample.no,w,range.min,range.max,data,eval.t=NULL){
  ## eval.t depends only on (n.b,range.min,range.max,data) and is constant
  ## across the EM iterations and across subjects, so the caller may pass a
  ## precomputed value to avoid an O(N^2)-per-call timeEval inside the
  ## per-subject E-step loop (which made the loop O(N^3)).
  if(is.null(eval.t)) eval.t<-timeEval(n.b,range.min,range.max,data)
  i<-sample.no
  no.f<-length(model$mod$factorModel)
  fr<-no.f
  samples<-length(unique(data$.id))
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))  
  samp<-i
  ev<-0
  for (jj in 1:no.f) {
    ev<-Matrix::bdiag(ev,t(eval.t$eval[[samp]][[jj]]))
  }
  ev<-ev[-1,-1]
  eval2<-ev
  mu.zw<-t(t(param[[i]]$params$mean.z))+param[[i]]$params$sigma.z%*%t(eval2)%*%
    solve(param[[i]]$params$sigma.w)%*%t(t(as.vector(w)-as.vector(param[[i]]$params$mean.w)))
  sigma.zw<-param[[i]]$params$sigma.z-param[[i]]$params$sigma.z%*%t(eval2)%*%
    solve(param[[i]]$params$sigma.w)%*%eval2%*%param[[i]]$params$sigma.z
  sigma.zw<-(sigma.zw+t(sigma.zw))/2
  count<-0
  extra<-list()
  for (j in 1:fr) {
    extra$mu.zw.j[[j]]<-mu.zw[(count+1):(count+n.b)]
    extra$sigma.zw.j[[j]]<-list()
    count1<-0
    for (jj in 1:fr) {
      extra$sigma.zw.j[[j]][[jj]]<-sigma.zw[(count+1):(count+n.b),(count1+1):(count1+n.b)]
      count1<-count1+n.b
    }
    count<-count+n.b
  }
  for (j in 1:fr) {
    if(j==1){
      extra$sigma.zw2[[j]]<-extra$sigma.zw.j[[j]][[j]] 
      extra$sigma.zw2.t[[j]]<-extra$sigma.zw.j[[j]][[j]]  
    }else{
      if(j==2){
        extra$sigma.zw2[[j]]<-extra$sigma.zw.j[[j]][[j-1]] 
        extra$sigma.zw2.t[[j]]<-extra$sigma.zw.j[[j-1]][[j]] 
      }else{
        extra$sigma.zw2[[j]]<-0
        extra$sigma.zw2.t[[j]]<-0
        for (k in 1:(j-1)) {
          s<-extra$sigma.zw.j[[j]][[k]]
          s2<-extra$sigma.zw.j[[k]][[j]]
          extra$sigma.zw2[[j]]<-cbind(extra$sigma.zw2[[j]],s)
          extra$sigma.zw2.t[[j]]<-rbind(extra$sigma.zw2.t[[j]],s2)
        }
        extra$sigma.zw2[[j]]<-extra$sigma.zw2[[j]][,-1]
        extra$sigma.zw2.t[[j]]<-extra$sigma.zw2.t[[j]][-1,]
      } 
    }
  }
  extra$sigma.zw3<-list()
  for (j in 1:(fr-1)) {
    if(j>=2){
      aa<-NULL
      for (jj in 1:j) {
        a<-NULL
        for (jj2 in 1:j) {
          a<-cbind(a,extra$sigma.zw.j[[jj]][[jj2]]) 
        }
        aa<-rbind(aa,a)
      }
      extra$sigma.zw3[[j]]<-aa
    }else{
      extra$sigma.zw3[[j]]<-extra$sigma.zw.j[[j]][[j]]
    }}
  extra
}

generator.z<-function(model,param,indicator.no,zw=NULL,mu.zw=NULL){
  extra<-list()
  j<-indicator.no
  if(j==1){
    mu.zw1<-param$mu.zw.j[[j]]
    sigma.zw1<-param$sigma.zw.j[[j]][[j]]
    zw<-mvrnorm(1,mu.zw1,sigma.zw1)
  }else{
    mu.zw1<-param$mu.zw.j[[j]]+param$sigma.zw2[[j]]%*%solve(param$sigma.zw3[[j-1]])%*%(zw-mu.zw)
    sigma.zw1<-param$sigma.zw.j[[j]][[j]]-param$sigma.zw2[[j]]%*%solve(param$sigma.zw3[[j-1]])%*%param$sigma.zw2.t[[j]]
    sigma.zw1<-(sigma.zw1+t(sigma.zw1))/2
    zw<-as.vector(mvrnorm(1,mu.zw1,sigma.zw1))
  }
  extra$zw.tot<-zw
  extra
}

###############EM algorithm#################################################
em.estimation<-function(model,data,x.data=NULL,
                        initial.parameter,n.b,n.em=100,n.monte=100,s.p=c("min","sd"),range.min=NULL,range.max=NULL,
                        design=c("regular", "irregular", "regular.truncated", "regular.missing"),
                        plot.progress=FALSE,
                        save.every =NULL #list(steps=100,dir_path="path directory to save", file_name_format = "simul_[step].rds",attached_info=NULL)
                        ){
  
  s.p = match.arg(s.p)
  # TODO: Change appropriately
  if(is.null(range.min))
    range.min=0
  if(is.null(range.max))
    range.max=1
  
  if(is.null(s.p))
    s.p="min"
  omegaMatrix<-omega.matrix(n.b,range.min,range.max,data)
  result<-list()
  no.em<-n.em
  no.it<-10 
  parameterr<-initial.parameter
  diagMatrix<-diag.matrix(model,n.b)
  diagMatrix2<-diag.matrix2(model,n.b)
  samples<-length(unique(data$.id))
  weightMatrix<-weight.matrix(model,n.b,sample=samples)
  no.fac<-length(model$mod$factorModel)
  no.lat<-length(model$var$latents)
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))  
  pen<-getbasispenalty(basis)
  ome<-inprod(basis,basis)
  eval.t<-timeEval(n.b,range.min,range.max,data)
  n.t<-max(c(sapply(1:samples,function(i){eval.t$time.no[[i]]})))
  evall2<-t(eval.basis(seq(0,1,length.out=n.t),basis))
  mont<-n.monte
  design = match.arg(design)
  if(design=="regular"){
    design.regular=TRUE
  }else{
    design.regular=FALSE
  }
  
  pb<-txtProgressBar(0,no.em,style = 3)
  pbk<-0
  
  res.list=list()
  cond<-list()
  res.k=0
  
  likeli.fac<-c()
  penal.fac<-c()
  likeli.sem<-c()
  penal.sem<-c()
  toler_mean<-c()
  toler_mean_sig<-c()
  ###############Algorithm#####################################
  for (i in 1:no.em) {
    
    eta.distt<-eta.distribution(model,param=parameterr,n.b,range.min,range.max,data,x.data,omegaMatrix =  omegaMatrix)
    paramz1<-params.z1(model,param=parameterr,n.b,range.min,range.max,eta.dist=eta.distt,data,omegaMatrix= omegaMatrix,design=design)
    time.start<-Sys.time()
    
    ################fac parameter estimations#######################
    generators2.z<-list()
    generators2.w<-list()
    for (i1 in 1:samples) {
      samp<-i1
      gg<-list()
      mu.gg<-list()
      generators2.z[[paste0(samp)]]<-list()
      generators2.w[[paste0(samp)]]<-list()
      zz<-0
      w.tot<-data$.value[data$.id==samp]
      paramz<-params.z(model,param=paramz1,n.b=n.b,
                       range.min=range.min,range.max=range.max,
                       sample.no=i1,w=w.tot,data)
      paramz11<-list()
      for (j1 in 1:no.fac) {
        generators2.z[[paste0(samp)]][[j1]]<-list()
        generators2.w[[paste0(samp)]][[j1]]<-list()
        generators3.w<-data$.value[data$.ind==j1&data$.id==samp]
        for (m in 1:mont) {
          if(j1==1){
            generators.z<-generator.z(model,param=paramz,
                                      indicator.no =j1)
            gg[[j1]]<-generators.z$zw.tot
            mu.gg[[j1]]<-as.vector(paramz$mu.zw.j[[j1]])
          }else{
            generators.z<-generator.z(model,param=paramz,indicator.no =j1,
                                      zw=gg[[j1-1]],mu.zw=mu.gg[[j1-1]])
            gg[[j1]]<-c(as.vector(gg[[j1-1]]),as.vector(generators.z$zw.tot))
            mu.gg[[j1]]<-c(mu.gg[[j1-1]],as.vector(paramz$mu.zw.j[[j1]]))
          }
          generators2.z[[paste0(samp)]][[j1]][[m]]<-generators.z$zw.tot
          generators2.w[[paste0(samp)]][[j1]][[m]]<-generators3.w
        }
      }
    }
    
    beta.s<-list()
    for (m in 1:mont) {
      beta.s[[m]]<-list()
      for (j1 in 1:no.lat) {
        beta.s[[m]][[j1]]<-list()
        for (i1 in 1:samples){
          beta.s[[m]][[j1]][[i1]]<-Inf
        }
      }
    }
    
    paramzz<-paramz1
    xi<-list()
    for (i1 in 1:samples) {
      samp<-i1
      xi[[paste0(samp)]]<-list()
      zz<-c(sapply(1:no.fac, function(j1){as.vector(rowMeans(sapply(1:mont,function(ii){generators2.z[[paste0(samp)]][[j1]][[ii]]})))}))
      for (m in 1:mont) {
        mu.eta.z<-t(t(paramzz[[i1]]$params$mu.eta))+paramzz[[i1]]$params$sigma.eta%*%
          t(paramzz[[i1]]$params$a.matrix)%*%solve(paramzz[[i1]]$params$sigma.z)%*%
          t(t(zz-paramzz[[i1]]$params$mean.z))
        sigma.eta.z<-paramzz[[i1]]$params$sigma.eta-paramzz[[i1]]$params$sigma.eta%*%
          t(paramzz[[i1]]$params$a.matrix)%*%solve(paramzz[[i1]]$params$sigma.z)%*%
          paramzz[[i1]]$params$a.matrix%*%paramzz[[i1]]$params$sigma.eta
        sigma.eta.z <- 0.5*(sigma.eta.z+t(sigma.eta.z))
        loop.count<-0
        repeat{
          eta.zz<-as.vector(mvrnorm(1,mu.eta.z,sigma.eta.z))
          a<-c()
          count2<-0
          for (i2 in 1:no.lat) {
            etaa<-eta.zz[(count2+1):(count2+n.b)]
            count2<-count2+n.b
            a[i2]<-t(etaa)%*%pen%*%etaa<=beta.s[[m]][[i2]][[samp]]
          }
          loop.count<-loop.count+1
          if(all(a)){
            break
          }
          if(loop.count>20){
            break
          }
        }
        xi[[paste0(samp)]][[m]]<-eta.zz
      }
    }
    
    eta.l<-list()
    for (m in 1:mont) {
      eta.l[[m]]<-list()
      for (k in 1:samples) {
        eta.l[[m]][[k]]<-list()
        count2<-0
        for (i2 in 1:no.lat) {
          eta.l[[m]][[k]][[i2]]<-xi[[paste0(k)]][[m]][(count2+1):(count2+n.b)]
          count2<-count2+n.b
        }
      }
    }
    
    for (m in 1:mont) {
      for (i2 in 1:no.lat) {
        mean<-rowMeans(sapply(1:samples, function(k){eta.l[[m]][[k]][[i2]]}))
        for (k in 1:samples) {
          eta.l[[m]][[k]][[i2]]<-eta.l[[m]][[k]][[i2]]-mean
        }
      }
    }
    
    if(design.regular==TRUE){
      soll<-ginv(eval.t$eval[[1]][[1]]%*%
                   t(eval.t$eval[[1]][[1]]))%*%
        eval.t$eval[[1]][[1]]
    }
    delta.fac<-c()
    hist.fac<-list()
    for (j in 1:no.fac) {
      aa<-c()
      for (k in 1:dim(diagMatrix2[[paste0("f",j)]])[1]) {
        if(all(diagMatrix2[[paste0("f",j)]][k,]==0)){
          aa[k]<-0
        }else{
          aa[k]<-1
        }
      }
      pen.fac<-getbasispenalty(basis)
      for (jj in 1:length(model$mod$factorModel[[j]]$effect)) {
        if(model$mod$factorModel[[j]]$effect[jj]%in%c("fixed_concurrent","fixed_historical")){
          pe<-NULL
        }
        if(model$mod$factorModel[[j]]$effect[jj]=="concurrent"){
          pe<-getbasispenalty(basis)
        }
        if(model$mod$factorModel[[j]]$effect[jj]=="historical"){
          P = matrix(0,n.b^2, n.b^2)
          Ps = matrix(0,n.b^2, n.b^2)
          Pt = matrix(0,n.b^2, n.b^2)
          
          times<-seq(0,1,length.out=100)
          E<-t(eval.basis(times,basis))
          E1 <- t(eval.basis(times,basis,Lfdobj=1))
          E2 <- t(eval.basis(times,basis,Lfdobj=2))
          for(l1 in seq_along(times)){
            for(l2 in seq_along(times)){
              kr = matrix(kronecker(E1[,l1],E1[,l2]),ncol = 1)
              P = P + kr%*%t(kr)
              
              kr = matrix(kronecker(E2[,l1],E[,l2]),ncol = 1)
              Ps = Ps + kr%*%t(kr)
              
              kr = matrix(kronecker(E[,l1],E2[,l2]),ncol = 1)
              Pt = Pt + kr%*%t(kr)
            }
          }
          P  = P  / (length(times)^2)
          Ps = Ps / (length(times)^2)
          Pt = Pt / (length(times)^2)
          
          P = P + Ps + Pt
          pe<-P
        }
        if(is.null(pe)){
          next
        }else{
          pen.fac<-Matrix::bdiag(pen.fac,pe) 
        }
      }
      f<-list()
      f2<-list()
      for (m in 1:mont) {
        f[[m]]<-list()
        f2[[m]]<-list()
        for (k in 1:samples) {
          if(design.regular==FALSE){
            soll<-ginv(eval.t$eval[[k]][[j]]%*%
                         t(eval.t$eval[[k]][[j]]))%*%
              eval.t$eval[[k]][[j]] 
          }
          f[[m]][[paste0(k)]]<-diag(n.b)
          for (i2 in 1:no.lat) {
            et<-kronecker(diag(eval.t$time.no[[k]][[j]]),t(eta.l[[m]][[k]][[i2]]))
            fac<-which(model$mod$factorModel[[j]]$factor==model$var$latents[i2])
            if(!length(fac)==0){
              if(model$mod$factorModel[[j]]$effect[fac]=="concurrent"){
                ef<-et%*%omegaMatrix$omega1[[k]][[j]]
                ef<- soll%*%ef
                f2[[m]][[paste0(k)]]<-0
              }
              if(model$mod$factorModel[[j]]$effect[fac]=="historical"){
                ef<-et%*%omegaMatrix$omega2[[k]][[j]] 
                ef<- soll%*%ef
                f2[[m]][[paste0(k)]]<-0
              } 
              if(model$mod$factorModel[[j]]$effect[fac]=="fixed_concurrent"){
                ef<-NULL
                f2[[m]][[paste0(k)]]<-eta.l[[m]][[k]][[i2]]
              }
              if(model$mod$factorModel[[j]]$effect[fac]=="fixed_historical"){
                ef<-NULL
                f2[[m]][[paste0(k)]]<-soll%*%
                  omegaMatrix$delta[[k]][[j]]%*%
                  eta.l[[m]][[k]][[i2]] 
              }
            }else{
              ef<-NULL
            }
            f[[m]][[paste0(k)]]<-cbind(f[[m]][[paste0(k)]],ef)
          }
        }
      }
      
      if(i>0){
        delt_sta<-c(0.001,0.01,0.05,0.1)*log(samples)
        delt <-delt_sta
        samp.train<-sort(sample(1:samples,round(0.7*samples)))
        samp.test<-sort((1:samples)[-samp.train]) 
        sig<-parameterr$sigma.fac[[j]]
        ff<-list()
        zz<-list()
        ff2<-list()
        for (m in 1:mont) {
          ff[[m]]<-0
          zz[[m]]<-0
          ff2[[m]]<-0
          for (i1 in samp.train) {
            ff[[m]]<-rbind(ff[[m]],f[[m]][[paste0(i1)]])
            zz[[m]]<-c(zz[[m]],generators2.z[[paste0(i1)]][[j]][[m]])
            ff2[[m]]<-c(ff2[[m]],f2[[m]][[paste0(i1)]])
          }
          ff[[m]]<-ff[[m]][-1,]
          zz[[m]]<-zz[[m]][-1]
          ff2[[m]]<-ff2[[m]][-1]
        }
        ff.test<-list()
        zz.test<-list()
        ff2.test<-list()
        for (m in 1:mont) {
          ff.test[[m]]<-0
          zz.test[[m]]<-0
          ff2.test[[m]]<-0
          for (i1 in samp.test) {
            ff.test[[m]]<-rbind(ff.test[[m]],f[[m]][[paste0(i1)]])
            zz.test[[m]]<-c(zz.test[[m]],generators2.z[[paste0(i1)]][[j]][[m]])
            ff2.test[[m]]<-c(ff2.test[[m]],f2[[m]][[paste0(i1)]])
          }
          ff.test[[m]]<-ff.test[[m]][-1,]
          zz.test[[m]]<-zz.test[[m]][-1]
          ff2.test[[m]]<-ff2.test[[m]][-1]
        }
        per<-c()
        for (l in 1:length(delt)) {
          for (it in 1:no.it) {
            sol<-solve(kronecker(diag(length(samp.train)),sig))
            lam.sum<-0
            lam1.sum<-0
            for (m in 1:mont) {
              lam<-t(ff[[m]])%*%sol%*%ff[[m]]+delt[l]*pen.fac
              lam1<-t(ff[[m]])%*%sol%*%(as.vector(zz[[m]])-ff2[[m]])
              lam.sum<-lam.sum+1/mont*lam
              lam1.sum<-lam1.sum+1/mont*lam1
            }
            lambda1<-solve(lam.sum)%*%lam1.sum
            lambda2.1<-lambda1
            
            sigma1<-0
            for (i1 in samp.train) {
              ss2<-0
              for (m in 1:mont) {
                ss<-as.vector(generators2.z[[paste0(i1)]][[j]][[m]])-
                  f2[[m]][[paste0(i1)]]-f[[m]][[paste0(i1)]]%*%lambda2.1
                ss2<-ss2+1/mont*(ss%*%t(ss))
              }
              sigma1<-sigma1+1/length(samp.train)*ss2
            }
            sig<-sigma1
          }
          sol<-solve(kronecker(diag(length(samp.test)),sig))
          likeli<-0
          for (m in 1:mont) {
            vall<-(-t(lambda2.1)%*%t(ff.test[[m]])%*%sol%*%(as.vector(zz.test[[m]])-ff2.test[[m]]))+
              t(lambda2.1)%*%((t(ff.test[[m]])%*%sol%*%ff.test[[m]])+delt[l]*pen.fac)%*%lambda2.1
            likeli<-likeli+vall/mont
          }
          per[l]<-as.vector(likeli)
        }
        delt.min<-delt[which(per==min(per))]
        ex1<-c()
        for (m in 1:mont) {
          vall<-(-t(lambda2.1)%*%t(ff.test[[m]])%*%sol%*%(as.vector(zz.test[[m]])-ff2.test[[m]]))+
            t(lambda2.1)%*%((t(ff.test[[m]])%*%sol%*%ff.test[[m]])+delt.min*pen.fac)%*%lambda2.1
          ex1[m]<-as.vector(vall)
        }
        bound<-mean(ex1)+sd(ex1)
        per<-c()
        for (l in 1:length(delt)) {
          likeli<-0
          for (m in 1:mont) {
            vall<-(-t(lambda2.1)%*%t(ff.test[[m]])%*%sol%*%(as.vector(zz.test[[m]])-ff2.test[[m]]))+
              t(lambda2.1)%*%((t(ff.test[[m]])%*%sol%*%ff.test[[m]])+delt[l]*pen.fac)%*%lambda2.1
            likeli<-likeli+1/mont*as.vector(vall)
          }
          if(likeli<=bound){
            per[l]<-TRUE
          }else{
            per[l]<-FALSE
          }
        }
        d.t<-delt[which(per)]
        if(s.p=="min"){
          deltaa<-delt.min
        }else{
          deltaa<-max(d.t)
        }
      }else{
        deltaa<-ifelse(n.bb>2,0.01,0)
      }
      delta.fac[j]<-deltaa
      ff<-list()
      zz<-list()
      ff2<-list()
      for (m in 1:mont) {
        ff[[m]]<-0
        zz[[m]]<-0
        ff2[[m]]<-0
        for (i1 in 1:samples) {
          ff[[m]]<-rbind(ff[[m]],f[[m]][[paste0(i1)]])
          zz[[m]]<-c(zz[[m]],generators2.z[[paste0(i1)]][[j]][[m]])
          ff2[[m]]<-c(ff2[[m]],f2[[m]][[paste0(i1)]])
        }
        ff[[m]]<-ff[[m]][-1,]
        zz[[m]]<-zz[[m]][-1]
        ff2[[m]]<-ff2[[m]][-1]
      }
      sig<-parameterr$sigma.fac[[j]]
      for (it in 1:no.it) {
        sol<-solve(kronecker(diag(samples),sig))
        lam.sum<-0
        lam1.sum<-0
        for (m in 1:mont) {
          lam<-t(ff[[m]])%*%sol%*%ff[[m]]+deltaa*pen.fac
          lam1<-t(ff[[m]])%*%sol%*%(as.vector(zz[[m]])-ff2[[m]])
          lam.sum<-lam.sum+1/mont*lam
          lam1.sum<-lam1.sum+1/mont*lam1
        }
        lambda1<-solve(lam.sum)%*%lam1.sum
        lambda2.1<-lambda1
        
        sigma1<-0
        for (i1 in 1:samples) {
          ss2<-0
          for (m in 1:mont) {
            ss<-as.vector(generators2.z[[paste0(i1)]][[j]][[m]])-
              f2[[m]][[paste0(i1)]]-f[[m]][[paste0(i1)]]%*%lambda2.1
            ss2<-ss2+1/mont*(ss%*%t(ss))
          }
          sigma1<-sigma1+1/samples*ss2
        }
        sig<-sigma1
      }
      
      lambda2<-c()
      count<-1
      for (k in 1:length(aa)) {
        if(aa[k]==1){
          lambda2[k]<-lambda2.1[count]
          count<-count+1
        }else{
          lambda2[k]<-0
        }
      }
      sig.err<-0
      for (i1 in 1:samples) {
        zsum<-0
        for (m in 1:mont) {
          sigg<-generators2.w[[paste0(i1)]][[j]][[1]]-t(eval.t$eval[[i1]][[j]])%*%generators2.z[[paste0(i1)]][[j]][[m]]
          zsum<-zsum+1/mont*(t(sigg)%*%sigg)
        }
        sig.err<-sig.err+zsum 
      }
      sig.err1<-as.vector(sig.err/sum(sapply(1:samples,function(i1)eval.t$time.no[[i1]][[j]])))
      parameterr$lambda.param[[j]]<-as.vector(lambda2)
      parameterr$sigma.fac[[j]]<-sig 
      parameterr$sigma.error[[j]]<-sig.err1
      hist.fac[[j]]<-list(generators2.z=generators2.z,generators2.w=generators2.w,
                          eta.l=eta.l,f=f,f2=f2,deltaa=deltaa,pen.fac=pen.fac)
    }
    ################sem parameter estimations##################  
    no.reg<-length(model$mod$regression)
    lat2<-model$var$latents
    obs<-model$var$observed
    no.lat2<-length(lat2)
    no.obs2<-length(obs)
    soll<-ginv(evall2%*%t(evall2))%*%evall2
    hist.sem<-list()
    for (j in 1:no.reg) {
      co<-model$mod$regression[[j]]$covariate
      no.cov<-length(co)
      lat.count<-which(lat2==model$mod$regression[[j]]$response)
      aa<-c()
      if(length(obs)>0){
        dd<-dim(Matrix::bdiag(diagMatrix2[[paste0("r",j,".eta")]],diagMatrix2[[paste0("r",j,".x")]]))
        for (kk in 1:dd[1]) {
          if(all(Matrix::bdiag(diagMatrix2[[paste0("r",j,".eta")]],diagMatrix2[[paste0("r",j,".x")]])[kk,]==0)){
            aa[kk]<-0
          }else{
            aa[kk]<-1
          }
        }
      }else{
        dd<-dim(diagMatrix2[[paste0("r",j,".eta")]])
        for (kk in 1:dd[1]) {
          if(all(diagMatrix2[[paste0("r",j,".eta")]][kk,]==0)){
            aa[kk]<-0
          }else{
            aa[kk]<-1
          }
        }
      }
      if(length(which(aa==1))>0){
        pen.sem<-0
        for (jj in 1:length(model$mod$regression[[j]]$effect)) {
          if(model$mod$regression[[j]]$effect[jj]=="concurrent"){
            pe<-getbasispenalty(basis)
          }
          if(model$mod$regression[[j]]$effect[jj]=="historical"){
            P = matrix(0,n.b^2, n.b^2)
            Ps = matrix(0,n.b^2, n.b^2)
            Pt = matrix(0,n.b^2, n.b^2)
            
            times<-seq(0,1,length.out=100)
            E<-t(eval.basis(times,basis))
            E1 <- t(eval.basis(times,basis,Lfdobj=1))
            E2 <- t(eval.basis(times,basis,Lfdobj=2))
            for(l1 in seq_along(times)){
              for(l2 in seq_along(times)){
                kr = matrix(kronecker(E1[,l1],E1[,l2]),ncol = 1)
                P = P + kr%*%t(kr)
                
                kr = matrix(kronecker(E2[,l1],E[,l2]),ncol = 1)
                Ps = Ps + kr%*%t(kr)
                
                kr = matrix(kronecker(E[,l1],E2[,l2]),ncol = 1)
                Pt = Pt + kr%*%t(kr)
              }
            }
            P  = P  / (length(times)^2)
            Ps = Ps / (length(times)^2)
            Pt = Pt / (length(times)^2)
            
            P = P + Ps + Pt
            pe<-P
          }
          if(model$mod$regression[[j]]$effect[jj]=="smooth"){
            Ps = matrix(0,n.b^2, n.b^2)
            basis3<-create.bspline.basis(rangeval = range(x.data[[model$mod$regression[[j]]$covariate[jj]]]),nbasis=n.b)
            sm<-t(eval.basis(x.data[[model$mod$regression[[j]]$covariate[jj]]],basis3))
            
            times<-seq(0,1,length.out=100)
            E2 <- t(eval.basis(times,basis,Lfdobj=2))
            for (d in 1:samples) {
              Pss<-matrix(0,n.b^2, n.b^2)
              for(l1 in seq_along(times)){
                kr = matrix(kronecker(E2[,l1],sm[,d]),ncol = 1)
                Pss = Pss + kr%*%t(kr)
              } 
              Ps<-Ps+Pss/ (length(times))
            }
            Ps = Ps 
            pe<-Ps
          }
          if(model$mod$regression[[j]]$effect[jj]=="linear"){
            pe<-getbasispenalty(basis)
          }
          pen.sem<-Matrix::bdiag(pen.sem,pe)
        }
        pen.sem<-pen.sem[-1,-1]
      }
      g.eta<-NULL
      g.x<-NULL
      if(no.cov>0){
        if(any(co%in%model$var$latents)){
          g.eta<-list()
          for (m in 1:mont) {
            g.eta[[m]]<-list()
            for (i1 in 1:samples) {
              g.eta.tmp=NULL
              for (jj in 1:no.lat2) {
                fac<-which(model$mod$regression[[j]]$covariate==lat2[jj])
                if(length(fac)==0){
                  gg.eta<-NULL
                }else{
                  if(model$mod$regression[[j]]$effect[fac]=="concurrent"){
                    gg.eta<-soll%*%
                      kronecker(diag(n.t),t(eta.l[[m]][[i1]][[jj]]))%*%
                      omegaMatrix$omega1.sem
                  }
                  if(model$mod$regression[[j]]$effect[fac]=="historical"){
                    gg.eta<-soll%*%
                      kronecker(diag(n.t),t(eta.l[[m]][[i1]][[jj]]))%*%
                      omegaMatrix$omega2.sem
                  }
                }
                g.eta.tmp<-cbind(g.eta.tmp,gg.eta)
              }
              g.eta[[m]][[i1]]=g.eta.tmp 
            }    
          }
          
        }else{
          g.eta<-lapply(1:mont,function(m){lapply(1:samples, function(i1)(g.eta[[m]][[i1]]=0))})
          et.ex<-"NULL"
        }
        if(length(obs)>0){
          if(any(co%in%model$var$observed)){
            g.x<-list()
            xx<-list()
            for (jj in 1:no.obs2) {
              fac<-which(model$mod$regression[[j]]$covariate==obs[jj])
              if(model$mod$regression[[j]]$effect[fac]%in%c("concurrent","historical")){
                sm<-matrix(NA,n.b,samples)
                for (tt in 1:samples) {
                  sm[,tt]<-Data2fd(argvals=seq(range.min,range.max,length.out=length(x.data[[obs[jj]]][[tt]])),y=x.data[[obs[jj]]][[tt]],basisobj = basis,lambda = 0.5)$coefs
                }
                meann<-rowMeans(sm)
                xx[[jj]]<-sm-meann     
              }
              if(model$mod$regression[[j]]$effect[fac]=="smooth"){
                basis3<-create.bspline.basis(rangeval = range(x.data[[model$var$observed[jj]]]),nbasis=n.b)
                sm<-t(eval.basis(x.data[[model$var$observed[jj]]],basis3))
                meann<-rowMeans(sm)
                xx[[jj]]<-sm-meann
              }
              if(model$mod$regression[[j]]$effect[fac]=="linear"){
                if(length(unique(x.data[[model$var$observed[jj]]]))==2){
                  xx[[jj]]<-x.data[[model$var$observed[jj]]]
                }else{
                  meann<-mean(x.data[[model$var$observed[jj]]])
                  xx[[jj]]<-x.data[[model$var$observed[jj]]]-meann
                }
              }
            }
            for (i1 in 1:samples) {
              g.x.tmp<-NULL
              for (jj in 1:no.obs2) {
                fac<-which(model$mod$regression[[j]]$covariate==obs[jj])
                if(length(fac)==0){
                  gg.x<-NULL
                }else{
                  if(model$mod$regression[[j]]$effect[fac]=="concurrent"){
                    omega<-omegaMatrix$omega1.sem
                    gg.x<-soll%*%
                      kronecker(diag(n.t),t(xx[[jj]][,i1]))%*%omega
                  }
                  if(model$mod$regression[[j]]$effect[fac]=="historical"){
                    omega<-omegaMatrix$omega2.sem
                    gg.x<-soll%*%
                      kronecker(diag(n.t),t(xx[[jj]][,i1]))%*%omega
                  }
                  if(model$mod$regression[[j]]$effect[fac]=="smooth"){
                    gg.x<-soll%*%
                      t(kronecker(evall2,xx[[jj]][,i1]))
                  }
                  if(model$mod$regression[[j]]$effect[fac]=="linear"){
                    gg.x<-as.vector(xx[[jj]][i1])*diag(n.b)
                  }
                }
                g.x.tmp<-cbind(g.x.tmp,gg.x)
                g.x[[i1]]<-g.x.tmp
              }
            }
          }else{
            g.x<-lapply(1:samples, function(i1)(g.x[[i1]]=0))
            x.ex<-"NULL"
          }
          
        }
      } 
      
      sig<-parameterr$sigma.sem[[lat.count]]
      if(!is.null(g.eta)|!is.null(g.x)){
        if(i>0){
          delt_sta<-c(0.001,0.01,0.05,0.1)*log(samples)
          delt <-delt_sta
          samp.train<-sort(sample(1:samples,round(0.7*samples)))
          samp.test<-sort((1:samples)[-samp.train])
          g<-list()  
          gg<-list()
          ett<-list()
          if(!is.null(g.eta)|!is.null(g.x)){
            for (m in 1:mont) {
              g[[m]]<-list()
              gg[[m]]<-0
              ett[[m]]<-0
              for (i1 in samp.train) {
                g[[m]][[i1]]<-cbind(g.eta[[m]][[i1]],g.x[[i1]])
                if(exists("et.ex")){
                  g[[m]][[i1]]<-g[[m]][[i1]][,-1]  
                }
                if(exists("x.ex")){
                  g[[m]][[i1]]<-g[[m]][[i1]][,-ncol(g[[m]][[i1]])]
                }
                gg[[m]]<-rbind(gg[[m]],g[[m]][[i1]])
                ett[[m]]<-c(ett[[m]],eta.l[[m]][[i1]][[lat.count]])
              }
              gg[[m]]<-gg[[m]][-1,]
              ett[[m]]<-ett[[m]][-1]
            }
            g.test<-list()  
            gg.test<-list()
            ett.test<-list()
            if(!is.null(g.eta)|!is.null(g.x))
              for (m in 1:mont) {
                g.test[[m]]<-list()
                gg.test[[m]]<-0
                ett.test[[m]]<-0
                for (i1 in samp.test) {
                  g.test[[m]][[i1]]<-cbind(g.eta[[m]][[i1]],g.x[[i1]])
                  if(exists("et.ex")){
                    g.test[[m]][[i1]]<-g.test[[m]][[i1]][,-1]  
                  }
                  if(exists("x.ex")){
                    g.test[[m]][[i1]]<-g.test[[m]][[i1]][,-ncol(g.test[[m]][[i1]])]
                  }
                  gg.test[[m]]<-rbind(gg.test[[m]],g.test[[m]][[i1]])
                  ett.test[[m]]<-c(ett.test[[m]],eta.l[[m]][[i1]][[lat.count]])
                }
                gg.test[[m]]<-gg.test[[m]][-1,]
                ett.test[[m]]<-ett.test[[m]][-1]
              }   
          }
          
          per<-c()
          for (l in 1:length(delt)) {
            for (it in 1:no.it) {
              sol<-solve(kronecker(diag(length(samp.train)),sig))
              gam.sum<-0
              gam1.sum<-0
              for (m in 1:mont) {
                gam<-t(gg[[m]])%*%sol%*%gg[[m]]+delt[l]*pen.sem
                gam1<-t(gg[[m]])%*%sol%*%ett[[m]]
                gam.sum<-gam.sum+1/mont*gam
                gam1.sum<-gam1.sum+1/mont*gam1
              }
              gamma1<-solve(gam.sum)%*%gam1.sum
              gamma2.1<-gamma1
              sigma1<-0
              for (i1 in samp.train) {
                ss2<-0
                for (m in 1:mont) {
                  if(is.null(g.eta) & is.null(g.x)){
                    ss<-eta.l[[m]][[i1]][[lat.count]]  
                  }else{
                    ss<-eta.l[[m]][[i1]][[lat.count]]- g[[m]][[i1]]%*%gamma2.1
                  }
                  ss2<-ss2+1/mont*ss%*%t(ss)
                }
                sigma1<-sigma1+1/length(samp.train)*ss2
              }
              sig<-sigma1
            }
            sol<-solve(kronecker(diag(length(samp.test)),sig))
            val<-0
            for (m in 1:mont) {
              vall<-(-t(gamma2.1)%*%t(gg.test[[m]])%*%sol%*%as.vector(ett.test[[m]]))+
                t(gamma2.1)%*%((t(gg.test[[m]])%*%sol%*%gg.test[[m]])+delt[l]*pen.sem)%*%gamma2.1 
              val<-val+1/mont*vall
            }
            per[l]<-as.vector(val)
          }
          delt.min<-delt[which(per==min(per))]
          ex1<-c()
          for (m in 1:mont) {
            vall<-(-t(gamma2.1)%*%t(gg.test[[m]])%*%sol%*%as.vector(ett.test[[m]]))+
              t(gamma2.1)%*%((t(gg.test[[m]])%*%sol%*%gg.test[[m]])+delt.min*pen.sem)%*%gamma2.1
            ex1[m]<-as.vector(vall)
          }
          bound<-mean(ex1)+sd(ex1)
          per<-c()
          for (l in 1:length(delt)) {
            likeli<-0
            for (m in 1:mont) {
              vall<-(-t(gamma2.1)%*%t(gg.test[[m]])%*%sol%*%(as.vector(ett.test[[m]])))+
                t(gamma2.1)%*%((t(gg.test[[m]])%*%sol%*%gg.test[[m]])+delt[l]*pen.sem)%*%gamma2.1
              likeli<-likeli+1/mont*as.vector(vall)
            }
            if(likeli<=bound){
              per[l]<-TRUE
            }else{
              per[l]<-FALSE
            }
          }
          d.t<-delt[which(per)]
          if(s.p=="min"){deltaa<-delt.min}else{deltaa<-max(d.t)}
        }else{
          deltaa<-0
        }
        delta.sem<-deltaa
        g<-list()  
        gg<-list()
        ett<-list()
        if(!is.null(g.eta)|!is.null(g.x))
          for (m in 1:mont) {
            g[[m]]<-list()
            gg[[m]]<-0
            ett[[m]]<-0
            for (i1 in 1:samples) {
              g[[m]][[i1]]<-cbind(g.eta[[m]][[i1]],g.x[[i1]])
              if(exists("et.ex")){
                g[[m]][[i1]]<-g[[m]][[i1]][,-1]  
              }
              if(exists("x.ex")){
                g[[m]][[i1]]<-g[[m]][[i1]][,-ncol(g[[m]][[i1]])]
              }
              gg[[m]]<-rbind(gg[[m]],g[[m]][[i1]])
              ett[[m]]<-c(ett[[m]],eta.l[[m]][[i1]][[lat.count]])
            }
            gg[[m]]<-gg[[m]][-1,]
            ett[[m]]<-ett[[m]][-1]
          }
        sig<-parameterr$sigma.sem[[lat.count]]
        sol<-solve(kronecker(diag(samples),sig))
        for (it in 1:no.it) {
          sol<-solve(kronecker(diag(samples),sig))
          gam.sum<-0
          gam1.sum<-0
          for (m in 1:mont) {
            gam<-t(gg[[m]])%*%sol%*%gg[[m]]+deltaa*pen.sem
            gam1<-t(gg[[m]])%*%sol%*%ett[[m]]
            gam.sum<-gam.sum+1/mont*gam
            gam1.sum<-gam1.sum+1/mont*gam1
          }
          gamma1<-solve(gam.sum)%*%gam1.sum
          gamma2.1<-gamma1
          sigma1<-0
          for (i1 in 1:samples) {
            ss2<-0
            for (m in 1:mont) {
              if(is.null(g.eta) & is.null(g.x)){
                ss<-eta.l[[m]][[i1]][[lat.count]]  
              }else{
                ss<-eta.l[[m]][[i1]][[lat.count]]- g[[m]][[i1]]%*%gamma2.1
              }
              ss2<-ss2+1/mont*ss%*%t(ss)
            }
            sigma1<-sigma1+1/samples*ss2
          }
          sig<-sigma1
        }
      }else{
        sigma1<-0
        for (i1 in 1:samples) {
          ss2<-0
          for (m in 1:mont) {
            if(is.null(g.eta) & is.null(g.x)){
              ss<-eta.l[[m]][[i1]][[lat.count]]  
            }else{
              ss<-eta.l[[m]][[i1]][[lat.count]]- g[[m]][[i1]]%*%gamma2.1
            }
            ss2<-ss2+1/mont*ss%*%t(ss)
          }
          sigma1<-sigma1+1/samples*ss2
        }
        sig<-sigma1
      }
      
      gamma2<-c()
      count<-1
      for (k in 1:length(aa)) {
        if(aa[k]==1){
          gamma2[k]<-gamma2.1[count]
          count<-count+1
        }else{
          gamma2[k]<-0
        }
      }
      parameterr$gamma.param[[lat.count]]<-as.vector(gamma2)
      parameterr$sigma.sem[[lat.count]]<-sig  
      if(length(which(aa==1))>0){
        hist.sem[[lat.count]]<-list(g.eta=g.eta,g.x=g.x,g=g,deltaa=deltaa,pen.sem=pen.sem)  
      }else{
        hist.sem[[lat.count]]<-list(g.eta=g.eta,g.x=g.x,deltaa=deltaa) 
      }
    }
    
    pbk<-pbk+1
    setTxtProgressBar(pb,pbk)
    time.end<-Sys.time()
    parameters<-ex.estimations(model,estimate=parameterr,n.b=n.b,range.min=range.min,range.max=range.max,data=data,x.data = x.data)
    parameterr<-parameters$totparam
    params.eval<-parameter.estimated.evaluted(model,param=parameters,n.b=n.b,n.monte=n.monte,range.min=range.min,range.max=range.max,data=data,x.data=x.data)
    
    ###convergence checking##################
    if(i>1){
      count<-0
      toll<-c()
      for (j in 1:no.fac) {
        ind<-model$mod$factorModel[[j]]$indicator
        for (jj in 1:length(model$mod$factorModel[[j]]$factor)) {
          fac<-model$mod$factorModel[[j]]$factor[[jj]]
          par<-params.eval$coef.fac.std[[paste0(ind)]][[paste0(fac)]]
          par.p<-param_prev$coef.fac.std[[paste0(ind)]][[paste0(fac)]]
          rim<-sqrt(mean((par-par.p)^2))
          count<-count+1
          toll[count]<-rim
        }
        par<-params.eval$intercept.std[[paste0(ind)]]
        par.p<-param_prev$intercept.std[[paste0(ind)]]
        rim<-sqrt(mean((par-par.p)^2))
        count<-count+1
        toll[count]<-rim
      }
      for (j in 1:length(model$mod$regression)) {
        res<-model$mod$regression[[j]]$response
        if(length(model$mod$regression[[j]]$covariate)>0){
          for (jj in 1:length(model$mod$regression[[j]]$covariate)) {
            cov<-model$mod$regression[[j]]$covariate[[jj]]
            par<-params.eval$coef.sem.std[[paste0(res)]][[paste0(cov)]]
            par.p<-param_prev$coef.sem.std[[paste0(res)]][[paste0(cov)]]
            rim<-sqrt(mean((par-par.p)^2))
            count<-count+1
            toll[count]<-rim
          }
          
        }
      }
      toler_mean[i]<-max(toll)
      toll_sig<-c()
      for (j in 1:length(model$mod$factorModel)) {
        tol<-abs(parameterr$sigma.error[[j]]-sig_prev[[j]])
        toll_sig[j]<-tol
      }
      toler_mean_sig[i]<-max(toll_sig)
      if(plot.progress){
        plot(1:length(toler_mean_sig[-1]),toler_mean_sig[-1],"l")
        plot(1:length(toler_mean[-1]),toler_mean[-1],"l")  
      }
      
    }
    #####################################
    
    
    if(plot.progress){
      t<-seq(0,1,length.out=200)
      par(mfrow=c(2,2))
      plot(t,params.eval$coef.fac.std$z1$eta,col="red","l")
      plot(t,params.eval$coef.fac.std$z2$eta,col="red","l")
      plot(t,params.eval$coef.fac.std$z3$eta,col="red","l")
      plot(t,params.eval$coef.fac.std$z4$eta,col="red","l")
      plot(t,params.eval$coef.fac.std$z5$eta,col="red","l")
      
      par(mfrow=c(2,2))
      
      plot(t,params.eval$coef.sem.std$eta$gender,col="red","l")
      plot(t,params.eval$coef.sem.std$eta$cancer,col="red","l")
      plot(t,params.eval$coef.sem.std$eta$genCancer,col="red","l")
      
      par(mfrow=c(2,2))
      
      plot(t,params.eval$sigma.sem.eigvec.std$eta[,1],col="red","l")
      plot(t,params.eval$sigma.sem.eigvec.std$eta[,2],col="red","l")
      plot(t,params.eval$sigma.sem.eigvec.std$eta[,3],col="red","l")
    }
    
    ##########################################################################
    param_prev<-params.eval
    sig_prev<-parameterr$sigma.error
    res.k=res.k+1
    res.list[[res.k]]=list(parameters=parameters,parameters.eval=params.eval,timing=c(time.start=time.start,time.end=time.end))
    
    if(plot.progress){
      print(parameterr$sigma.error)
    }
    
    if(!is.null(save.every)){
      #save.every = list(steps=100,dir_path="path directory to save", file_name_format = "simul_[step].rds",attached_info=NULL)
      if(!is.null(save.every$steps) && (i %% save.every$steps)==0){
        #save current step
        file_name = gsub("[step]",i,tolower(save.every$file_name_format),fixed = TRUE)
        file_path = file.path(dir_path,file_name)
        
        result2=result
        result2$params.estimated<-parameters
        result2$params.estimated.eval<-params.eval
        res_to_save<-list(result=result2,
                    tolerance=list(tol=toler_mean,
                                   tol_sig=toler_mean_sig),
                    steps=i,
                    hist=list(hist.fac=hist.fac,hist.sem=hist.sem,eta=eta.l),
                    attached_info = save.every$attached_info
                    )
        saveRDS(res_to_save, file_path)
        
      }
    }
  }
  result$params.estimated<-parameters
  result$params.estimated.eval<-params.eval
  final<-list(result=result,
              tolerance=list(tol=toler_mean,tol_sig=toler_mean_sig),steps=i,hist=list(hist.fac=hist.fac,hist.sem=hist.sem,eta=eta.l))
  final
}

############extracting estimations#################################

ex.estimations<-function(model,estimate,n.b,range.min,range.max,data,x.data=NULL,n.b.sim=NULL){
  estimations<-estimate
  no.obs<-length(model$var$observed)
  no.f<-length(model$mod$factorModel)
  samples<-length(unique(data$.id))
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))
  if(is.null(n.b.sim))
    n.b.sim=100
  
  n.t<-max(100,2*n.b.sim)
  M<-n.t
  eval3<-t(eval.basis(seq(range.min,range.max,length.out=M),basis))
  eval1<-eval3
  eval<-eval1
  ome<-inprod(basis,basis)
  a<-sapply(1:n.b,function(ig){1/100*sum(eval3[ig,1:99])})
  init<-list()
  no.r2<-length(model$mod$regression)
  lat2<-model$var$latents
  obs2<-model$var$observed
  init$gamma.param.eta<-list()
  init$gamma.param.x<-list()
  init$sigma.sem<-list()
  init$sigma.sem.eigvec<-list()
  init$sigma.sem.eigval<-list()
  init$ker.fac<-list()
  init$ker.sem<-list()
  init$totparam$gamma.param<-list()
  init$totparam$sigma.sem<-list()
  
  init$gamma.param.eta.std<-list()
  init$gamma.param.x.std<-list()
  init$sigma.sem.std<-list()
  init$sigma.sem.eigvec.std<-list()
  init$sigma.sem.eigval.std<-list()
  init$ker.fac.std<-list()
  init$ker.sem.std<-list()
  weightMatrix<-weight.matrix(model,n.b,sample=samples)
  soll<-ginv(eval1%*%t(eval1))%*%eval1
  x.data1<-list()
  no.obs2<-model$var$observed
  if(length(no.obs2)>0){
    if(any(c("concurrent","smooth","historical")%in%sapply(1:length(model$mod$regression),function(ii){model$mod$regression[[ii]]$effect}))){
      x.data1<-NULL
      for (j in 1:length(no.obs2)) {
        dimm<-is.list(x.data[[no.obs2[j]]])
        if(dimm){
          sm<-matrix(NA,n.b,samples)
          for (tt in 1:samples) {
            sm[,tt]<-Data2fd(argvals=seq(range.min,range.max,length.out=length(x.data[[no.obs2[j]]][[tt]])),y=x.data[[no.obs2[j]]][[tt]],basisobj = basis,lambda = 0.5)$coefs
          }
          meann<-rowMeans(sm)
          xx<-sm-meann          
        }else{
          basis3<-create.bspline.basis(rangeval = range(x.data[[model$var$observed[j]]]),nbasis=n.b)
          sm<-c(t(eval.basis(x.data[[model$var$observed[j]]],basis3)))
          xx<-sm
        }
        x.data1<-c(x.data1,xx)
      }
      xx.data<-weightMatrix$x%*%x.data1 
    }
  }
  for (i in 1:no.r2) {
    count<-which(model$var$latents==model$mod$regression[[i]]$response)
    d.sem1<-t(eval1)%*%estimations$sigma.sem[[count]]%*%eval1
    count2<-0
    ind<-c()
    count3<-0
    init$gamma.param.eta[[count]]<-list()
    init$totparam$gamma.param[[count]]<-0
    init$gamma.param.eta.std[[count]]<-list()
    for (j in 1:length(lat2)){
      fac<-which(model$mod$regression[[i]]$covariate==lat2[j])
      count3<-count3+1
      d.sem2<-t(eval1)%*%estimations$sigma.sem[[count3]]%*%eval1
      if(length(fac)==0){
        gam.eta<-estimations$gamma.param[[count]][(count2+1):(count2+n.b)]
        count2<-count2+n.b
        gam.eta.std<-gam.eta
      }else{
        if(model$mod$regression[[i]]$effect[fac]=="concurrent"){
          gam.eta<-estimations$gamma.param[[count]][(count2+1):(count2+n.b)]
          count2<-count2+n.b
          gam.eta.std<-diag(diag(d.sem1)^(-1/2))%*%
            diag(diag(d.sem2))^(1/2)%*%t(eval)%*%gam.eta
        }
        if(model$mod$regression[[i]]$effect[fac]=="historical"){
          gam.eta<-estimations$gamma.param[[count]][(count2+1):(count2+n.b*n.b)]
          count2<-count2+n.b*n.b
          gam.eta2<-matrix(gam.eta,n.b,n.b)
          tt<-seq(range.min,range.max,length.out=M)
          evv<-t(eval.basis(tt,basis))
          gam.est<-matrix(NA,M,M)
          evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
          for (ii in 1:M) {
            la<-t(evv_ind[[ii]])%*%matrix(gam.eta2,n.b,n.b)%*%evv[,ii]
            gam.est[,ii]<-la
          }
          gam.eta.std<-t(diag(d.sem1)^(-1/2)*
                           diag(d.sem2)^(1/2)*t(gam.est))
        }
      }
      init$gamma.param.eta[[count]][[count3]]<-gam.eta
      init$totparam$gamma.param[[count]]<-c(init$totparam$gamma.param[[count]],gam.eta)
      init$gamma.param.eta.std[[count]][[count3]]<-gam.eta.std
    }
    init$totparam$gamma.param[[count]]<-init$totparam$gamma.param[[count]][-1]
    
    count3<-0
    if(length(obs2>0)){
      init$gamma.param.x[[count]]<-list()
      init$gamma.param.x.std[[count]]<-list()
      for (j in 1:length(obs2)) {
        fac<-which(model$mod$regression[[i]]$covariate==obs2[j])
        count3<-count3+1
        if(length(fac)==0){
          gam.x<-estimations$gamma.param[[count]][(count2+1):(count2+n.b)]
          count2<-count2+n.b
          gam.x.std<-gam.x
        }else{
          if(model$mod$regression[[i]]$effect[fac]=="concurrent"){
            gam.x<-estimations$gamma.param[[count]][(count2+1):(count2+n.b)]
            count2<-count2+n.b
            gam.x.std<-diag(diag(d.sem1)^(-1/2))%*%t(eval)%*%gam.x
          }
          if(model$mod$regression[[i]]$effect[fac]=="historical"){
            gam.x<-estimations$gamma.param[[count]][(count2+1):(count2+n.b*n.b)]
            count2<-count2+n.b*n.b
            gam.x2<-matrix(gam.x,n.b,n.b)
            tt<-seq(range.min,range.max,length.out=M)
            evv<-t(eval.basis(tt,basis))
            gam.est<-matrix(NA,M,M)
            evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
            for (ii in 1:M) {
              la<-t(evv_ind[[ii]])%*%matrix(gam.x2,n.b,n.b)%*%evv[,ii]
              gam.est[,ii]<-la
            }
            gam.x.std<-t(diag(d.sem1)^(-1/2)*
                           t(gam.est))
          }
          if(model$mod$regression[[i]]$effect[fac]=="smooth"){
            cova<-model$mod$regression[[i]]$covariate[fac]
            gam.x<-estimations$gamma.param[[count]][(count2+1):(count2+n.b*n.b)]
            count2<-count2+n.b*n.b
            gam.x2<-matrix(gam.x,n.b,n.b)
            gam.x.std<-diag(d.sem1)^(-1/2)*t(eval)%*%gam.x2%*%
              matrix(x.data1[[cova]],n.b)
          }
          if(model$mod$regression[[i]]$effect[fac]=="linear"){
            gam.x<-estimations$gamma.param[[count]][(count2+1):(count2+n.b)]
            count2<-count2+n.b
            gam.x.std<-diag(diag(d.sem1)^(-1/2))%*%t(eval)%*%gam.x
          }
          if(model$mod$regression[[i]]$effect[fac]=="constant"){
            gam.x<-estimations$gamma.param[[count]][(count2+1):(count2+n.b)]
            count2<-count2+n.b
            gam.x.std<-t(diag(diag(d.sem1)^(-1/2))%*%t(eval)%*%a%*%gam.x)
          }
        }
        init$gamma.param.x[[count]][[count3]]<-gam.x
        init$totparam$gamma.param[[count]]<-c(init$totparam$gamma.param[[count]],gam.x)
        init$gamma.param.x.std[[count]][[count3]]<-gam.x.std
      }
    }
    init$ker.sem[[count]]<-d.sem1
    init$sigma.sem[[count]]<-as.matrix(estimations$sigma.sem[[count]])
    init$totparam$sigma.sem[[count]]<-init$sigma.sem[[count]]
    init$ker.sem.std[[count]]<-diag(diag(d.sem1)^(-1/2))%*%d.sem1%*%diag(diag(d.sem1)^(-1/2))
    init$sigma.sem.std[[count]]<-soll%*%
      diag(diag(d.sem1)^(-1/2))%*%d.sem1%*%
      diag(diag(d.sem1)^(-1/2))%*%t(soll)
    
  }
  
  fr<-length(model$mod$factorModel)
  no.lat<-length(model$var$latents)
  for (i in 1:fr) {
    d.sem1<-t(eval1)%*%estimations$sigma.fac[[i]]%*%eval1
    count2<-0
    init$lambda.param[[i]]<-list()
    init$sigma.fac[[i]]<-list()
    init$lambda.param.std[[i]]<-list()
    init$sigma.fac.std[[i]]<-list()
    if(model$mod$factorModel[[i]]$intercept==TRUE){
      intercept<-estimations$lambda.param[[i]][(count2+1):(count2+n.b)]
      count2<-count2+n.b
      init$intercept[[i]]<-intercept
      init$intercept.std[[i]]<-intercept
    }
    count3<-0
    lam<-0
    lam.std<-0
    for (j in 1:no.lat) {
      fac<-which(model$mod$factorModel[[i]]$factor==model$var$latents[j])
      count3<-count3+1
      d.sem2<-t(eval1)%*%estimations$sigma.sem[[count3]]%*%eval1
      if(length(fac)==0){
        lambda<-estimations$lambda.param[[i]][(count2+1):(count2+n.b)]
        count2<-count2+n.b
        lambda.std<-lambda
      }else{
        if(model$mod$factorModel[[i]]$effect[fac]=="concurrent"){
          lambda<-estimations$lambda.param[[i]][(count2+1):(count2+n.b)]
          count2<-count2+n.b
          lambda.std<-diag(d.sem2)^(1/2)*(t(eval)%*%lambda)
        }
        if(model$mod$factorModel[[i]]$effect[fac]=="historical"){
          lambda<-estimations$lambda.param[[i]][(count2+1):(count2+n.b*n.b)]
          count2<-count2+n.b*n.b
          lambda2<-matrix(lambda,n.b,n.b)
          tt<-seq(range.min,range.max,length.out=M)
          evv<-t(eval.basis(tt,basis))
          lam.est<-matrix(NA,M,M)
          evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
          for (ii in 1:M) {
            la<-t(evv_ind[[ii]])%*%matrix(lambda2,n.b,n.b)%*%evv[,ii]
            lam.est[,ii]<-la
          }
          lambda.std<-t(diag(d.sem2)^(1/2)*t(lam.est))
        }
        if(model$mod$factorModel[[i]]$effect[fac]=="fixed_concurrent"){
          lambda<-estimations$lambda.param[[i]][(count2+1):(count2+n.b)]
          count2<-count2+n.b
          lambda.std<-diag(d.sem2)^(1/2)
        }
        if(model$mod$factorModel[[i]]$effect[fac]=="fixed_historical"){
          lambda<-estimations$lambda.param[[i]][(count2+1):(count2+n.b)]
          count2<-count2+n.b
          lambda.std<-diag(d.sem2)^(1/2)
        }
      }
      init$lambda.param[[i]][[count3]]<-lambda
      lam<-c(lam,lambda)
      init$lambda.param.std[[i]][[count3]]<-lambda.std
    }
    init$ker.fac[[i]]<-d.sem1
    init$totparam$lambda.param[[i]]<-c(init$intercept[[i]],lam[-1])
    init$sigma.fac[[i]]<-as.matrix(estimations$sigma.fac[[i]])
    init$totparam$sigma.error[[i]]<-estimations$sigma.error[[i]]
    init$totparam$sigma.fac[[i]]<-init$sigma.fac[[i]]
    init$ker.fac.std[[i]]<-d.sem1
    init$sigma.fac.std[[i]]<-init$sigma.fac[[i]]
  }
  init
}

parameter.estimated.evaluted<-function(model,param,n.b,n.b.sim=NULL,n.monte,range.min,range.max,data,x.data=NULL){
  estimates<-param
  no.fac<-length(model$mod$factorModel)
  no.r<-length(model$mod$regression)
  latents<-model$var$latents
  samples<-length(unique(data$.id))
  
  if(is.null(n.b.sim))
    n.b.sim=100
  
  n.t<-max(100,2*n.b.sim)
  M<-n.t
  t<-seq(range.min,range.max,length.out=M)
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))
  eval<-t(eval.basis(seq(range.min,range.max,length.out=M),basis))
  eval3<-t(eval.basis(seq(range.min,range.max,length.out=100),basis))
  ome<-inprod(basis,basis)
  ome.val<-eigen(ome)$values
  ome.vec<-eigen(ome)$vectors
  ome.half<-ome.vec%*%(diag(ome.val)^0.5)%*%t(ome.vec)
  soll.om<-solve(ome)
  a<-sapply(1:n.b,function(ig){1/100*sum(eval3[ig,1:99])})
  weightMatrix<-weight.matrix(model,n.b,sample=samples)
  x.data1<-list()
  no.obs2<-model$var$observed
  if(length(no.obs2)>0){
    if(any(c("concurrent","smooth","historical")%in%sapply(1:length(model$mod$regression),function(ii){model$mod$regression[[ii]]$effect}))){
      x.data1<-NULL
      for (j in 1:length(no.obs2)) {
        dimm<-is.list(x.data[[no.obs2[j]]])
        if(dimm){
          sm<-matrix(NA,n.b,samples)
          for (tt in 1:samples) {
            sm[,tt]<-Data2fd(argvals=seq(range.min,range.max,length.out=length(x.data[[no.obs2[j]]][[tt]])),y=x.data[[no.obs2[j]]][[tt]],basisobj = basis,lambda = 0.5)$coefs
          }
          meann<-rowMeans(sm)
          xx<-sm-meann          
        }else{
          basis3<-create.bspline.basis(rangeval = range(x.data[[model$var$observed[j]]]),nbasis=n.b)
          sm<-c(t(eval.basis(x.data[[model$var$observed[j]]],basis3)))
          xx<-sm
        }
        x.data1<-c(x.data1,xx)
      }
      xx.data<-weightMatrix$x%*%x.data1 
    }
  }
  init<-list()
  lat<-model$var$latents
  no.r2<-no.r
  lat2<-latents
  obs2<-no.obs2
  for (j in 1:no.r2) {
    count<-which(model$var$latents==model$mod$regression[[j]]$response)
    res<-model$mod$regression[[j]]$response
    init$coef.sem[[res]]<-list()
    init$coef.sem.std[[res]]<-list()
    ind<-c()
    count1<-0
    for (jj in 1:length(lat2)) {
      count1<-count1+1
      cov<-which(model$mod$regression[[j]]$covariate==lat2[jj])
      if(!length(cov)==0){
        eff<-model$mod$regression[[j]]$effect[cov]
        cova<-model$mod$regression[[j]]$covariate[cov]
        if(eff=="concurrent"){
          init$coef.sem[[res]][[cova]]<-t(eval)%*%estimates$gamma.param.eta[[count]][[count1]]
          init$coef.sem.std[[res]][[cova]]<-estimates$gamma.param.eta.std[[count]][[count1]]
        }
        if(eff=="historical"){
          tt<-seq(range.min,range.max,length.out=M)
          evv<-t(eval.basis(tt,basis))
          gam.est<-matrix(NA,M,M)
          evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
          for (ii in 1:M) {
            la<-t(evv_ind[[ii]])%*%matrix(estimates$gamma.param.eta[[count]][[count1]],n.b,n.b)%*%evv[,ii]
            gam.est[,ii]<-la
          }
          init$coef.sem[[res]][[cova]]<-gam.est
          init$coef.sem.std[[res]][[cova]]<-estimates$gamma.param.eta.std[[count]][[count1]]
        }
      }
    }
    ind<-c()
    count1<-0
    if(length(obs2)>0)
      for (jj in 1:length(obs2)) {
        count1<-count1+1
        cov<-which(model$mod$regression[[j]]$covariate==obs2[jj])
        if(!length(cov)==0){
          eff<-model$mod$regression[[j]]$effect[cov]
          cova<-model$mod$regression[[j]]$covariate[cov]
          if(eff=="concurrent"){
            init$coef.sem[[res]][[cova]]<-t(eval)%*%estimates$gamma.param.x[[count]][[count1]]
            init$coef.sem.std[[res]][[cova]]<-estimates$gamma.param.x.std[[count]][[count1]]
          }
          if(eff=="historical"){
            tt<-seq(range.min,range.max,length.out=M)
            evv<-t(eval.basis(tt,basis))
            gam.est<-matrix(NA,M,M)
            evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
            for (ii in 1:M) {
              la<-t(evv_ind[[ii]])%*%matrix(estimates$gamma.param.x[[count]][[count1]],n.b,n.b)%*%evv[,ii]
              gam.est[,ii]<-la
            }
            init$coef.sem[[res]][[cova]]<-gam.est
            init$coef.sem.std[[res]][[cova]]<-estimates$gamma.param.x.std[[count]][[count1]]
          }
          if(eff=="smooth"){
            init$coef.sem[[res]][[cova]]<-t(eval)%*%
              matrix(estimates$gamma.param.x[[count]][[count1]],n.b,n.b)%*%
              matrix(x.data1[[cova]],n.b)
            init$coef.sem.std[[res]][[cova]]<-estimates$gamma.param.x.std[[count]][[count1]]
          }
          if(eff=="linear"){
            init$coef.sem[[res]][[cova]]<-t(eval)%*%estimates$gamma.param.x[[count]][[count1]]
            init$coef.sem.std[[res]][[cova]]<-estimates$gamma.param.x.std[[count]][[count1]]
          }
          if(eff=="constant"){
            init$coef.sem[[res]][[cova]]<-t(eval)%*%a%*%
              t(estimates$gamma.param.x[[count]][[count1]])%*%
              matrix(x.data1[[cova]],n.b)
            init$coef.sem.std[[res]][[cova]]<-estimates$gamma.param.x.std[[count]][[count1]]
          }
        }
      }
    eig<-eigen(ome.half%*%
                 estimates$sigma.sem[[count]]%*%
                 ome.half)
    init$ker.sem[[res]]<-estimates$ker.sem[[count]]
    init$sigma.sem.eigvec[[res]]<-t(eval)%*%soll.om%*%ome.half%*%eig$vectors
    init$sigma.sem.eigval[[res]]<-eig$values
    eig.std<-eigen(ome.half%*%
                     estimates$sigma.sem.std[[count]]%*%
                     ome.half)
    init$ker.sem.std[[res]]<-estimates$ker.sem.std[[count]]
    init$sigma.sem.eigvec.std[[res]]<-t(eval)%*%soll.om%*%ome.half%*%eig.std$vectors
    init$sigma.sem.eigval.std[[res]]<-eig.std$values
  }
  for (j in 1:no.fac) {
    res<-model$mod$factorModel[[j]]$indicator
    init$coef.fac[[res]]<-list()
    init$intercept[[res]]<-t(eval)%*%estimates$intercept[[j]]
    init$coef.fac.std[[res]]<-list()
    init$intercept.std[[res]]<-t(eval)%*%estimates$intercept.std[[j]]
    count1<-0
    for (jj in 1:length(latents)) {
      count1<-count1+1
      cov<-which(model$mod$factorModel[[j]]$factor==latents[jj])
      if(!length(cov)==0){
        cova<-model$mod$factorModel[[j]]$factor[cov]
        eff<-model$mod$factorModel[[j]]$effect[cov]
        if(eff=="concurrent"){
          init$coef.fac[[res]][[cova]]<-t(eval)%*%estimates$lambda.param[[j]][[count1]]
          init$coef.fac.std[[res]][[cova]]<-estimates$lambda.param.std[[j]][[count1]]
        }
        if(eff=="historical"){
          tt<-seq(range.min,range.max,length.out=M)
          evv<-t(eval.basis(tt,basis))
          lam.est<-matrix(NA,M,M)
          evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
          for (ii in 1:M) {
            la<-t(evv_ind[[ii]])%*%matrix(estimates$lambda.param[[j]][[count1]],n.b,n.b)%*%evv[,ii]
            lam.est[,ii]<-la
          }
          init$coef.fac[[res]][[cova]]<-lam.est
          init$coef.fac.std[[res]][[cova]]<-estimates$lambda.param.std[[j]][[count1]]
        }
        if(eff=="fixed_concurrent"){
          init$coef.fac[[res]][[cova]]<-c()
          init$coef.fac.std[[res]][[cova]]<-estimates$lambda.param.std[[j]][[count1]]
        }
        if(eff=="fixed_historical"){
          init$coef.fac[[res]][[cova]]<-c()
          init$coef.fac.std[[res]][[cova]]<-estimates$lambda.param.std[[j]][[count1]]
        }
        break
      }
    }
    eigg<-eigen(ome.half%*%
                  estimates$sigma.fac[[j]]%*%
                  ome.half)
    init$ker.fac[[res]]<-estimates$ker.fac[[j]]
    init$sigma.fac.eigvec[[res]]<-t(eval)%*%soll.om%*%ome.half%*%eigg$vectors
    init$sigma.fac.eigval[[res]]<-eigg$values
    eigg.std<-eigen(ome.half%*%
                      estimates$sigma.fac.std[[j]]%*%
                      ome.half)
    init$ker.fac.std[[res]]<-estimates$ker.fac.std[[j]]
    init$sigma.fac.eigvec.std[[res]]<-t(eval)%*%soll.om%*%ome.half%*%eigg.std$vectors
    init$sigma.fac.eigval.std[[res]]<-eigg.std$values
  }
  init
}

