library(mvtnorm)

###model########################################
unravel.formula<-function(formula){
  a<-list()
  name<-all.names(formula)
  if(name[3]=="~"){
    a[[1]]<-"factorModel"
  }else{
    a[[1]]<-"regression"
  }
  a[[2]]<-name
  plus<-length(which(name=="+"))
  if(a[[1]]=="regression"){
    a[[3]]<-ifelse(name[plus+3]=="-",plus,plus+1)
  }else{
    a[[3]]<-ifelse(name[plus+4]=="-",plus,plus+1)  
  }
  if(a[[1]]=="regression"){
    if(name[plus+3]=="-"){
      a[[4]]<-TRUE 
    }else{
      a[[4]]<-FALSE 
    } 
  }else{
    if(name[plus+4]=="-"){
      a[[4]]<-TRUE 
    }else{
      a[[4]]<-FALSE 
    } 
  }
  a
}

fsem<-function(formula,effectType,latent.covariate=NULL,scalar.covariate=NULL){
  ur<-unravel.formula(formula)
  b<-list()
  l<-length(ur[[2]])
  if(ur[[1]]=="factorModel"){
    b[["mod"]][["factorModel"]]<-list()
    b[["var"]][["indicators"]]<-ur[[2]][(l-ur[[3]]+1):l]
    b[["var"]][["factors"]]<-ur[[2]][[2]]  
    ind<-ur[[3]]
    for (i in 1:ind) {
      b[["mod"]][["factorModel"]][[paste0("f",i)]][["type"]]<-ur[[1]]
      b[["mod"]][["factorModel"]][[paste0("f",i)]][["indicator"]]<-ur[[2]][(l-ur[[3]]+1):l][i] 
      b[["mod"]][["factorModel"]][[paste0("f",i)]][["factor"]]<-ur[[2]][[2]]
      b[["mod"]][["factorModel"]][[paste0("f",i)]][["effect"]]<-effectType
      b[["mod"]][["factorModel"]][[paste0("f",i)]][["intercept"]]<-ifelse(ur[[4]]==TRUE,FALSE,TRUE)
    }
    b[["var"]][["latents"]]<-ur[[2]][[2]]
  }else{
    b[["mod"]][["regression"]]<-list()
    b[["mod"]][["regression"]][["r1"]][["type"]]<-ur[[1]]
    b[["mod"]][["regression"]][["r1"]][["response"]]<-ur[[2]][[2]]
    b[["mod"]][["regression"]][["r1"]][["covariate"]]<-ur[[2]][(l-ur[[3]]+1):l]
    if(ur[[4]]==TRUE&ur[[2]][l]=="-"){
      b[["mod"]][["regression"]][["r1"]][["effect"]]<-NULL
      b[["mod"]][["regression"]][["r1"]][["covariate"]]<-NULL
      b[["mod"]][["regression"]][["r1"]][["scalar.covariate"]]<-NULL
    }else{
      if(is.null(scalar.covariate))
        scalar.covariate<-TRUE
      
      b[["mod"]][["regression"]][["r1"]][["effect"]]<-rep(effectType,length(b[["mod"]][["regression"]][["r1"]][["covariate"]]))
    }
    b[["mod"]][["regression"]][["r1"]][["intercept"]]<-ifelse(ur[[4]]==TRUE,FALSE,TRUE)
    wl<-which(latent.covariate%in%b[["mod"]][["regression"]][["r1"]][["covariate"]])
    b[["mod"]][["regression"]][["r1"]][["scalar.covariate"]]<-scalar.covariate
    if(!(length(wl)==0)){
      b[["var"]][["latents"]]<-latent.covariate[wl]
      if(!(length(b[["mod"]][["regression"]][["r1"]][["covariate"]][-wl])==0))
        b[["var"]][["observed"]]<-b[["mod"]][["regression"]][["r1"]][["covariate"]][-wl]
    }else{
      if(!(length(b[["mod"]][["regression"]][["r1"]][["covariate"]])==0))
        b[["var"]][["observed"]]<-b[["mod"]][["regression"]][["r1"]][["covariate"]]   
    }
  }
  b
}

GenModel<-function(fsem1,fsem2){
  l_fsem1<-length(fsem1[["mod"]])
  count1<-rep(0,l_fsem1)
  count2<-"FALSE"
  for(i in 1:l_fsem1){
    if(!(fsem1[["mod"]][[i]][[1]][[1]]==fsem2[["mod"]][[1]][[1]][[1]])){
      count1[i]<-1
    }else{
      count<-matrix(0,length(fsem1[["mod"]][[i]]),length(fsem2[["mod"]][[1]]))
      if(fsem1[["mod"]][[i]][[1]][[1]]=="factorModel"&fsem2[["mod"]][[1]][[1]][[1]]=="factorModel"){
        for (j in 1:length(fsem1[["mod"]][[i]])) {
          for (k in 1:length(fsem2[["mod"]][[1]])) {
            if(fsem1[["mod"]][[i]][[j]][[2]]==fsem2[["mod"]][[1]][[k]][[2]]){
              fsem1[["mod"]][[i]][[j]][["type"]]<-fsem1[["mod"]][[i]][[j]][[1]]
              fsem1[["mod"]][[i]][[j]][["indicator"]]<-fsem1[["mod"]][[i]][[j]][[2]]
              fsem1[["mod"]][[i]][[j]][["factor"]]<-c(fsem1[["mod"]][[i]][[j]][[3]],fsem2[["mod"]][[1]][[k]][[3]])
              fsem1[["mod"]][[i]][[j]][["effect"]]<-c(fsem1[["mod"]][[i]][[j]][[4]],fsem2[["mod"]][[1]][[k]][[4]])
              fsem1[["mod"]][[i]][[j]][["intercept"]]<-ifelse(fsem1[["mod"]][[i]][[j]][["intercept"]]==TRUE,TRUE,FALSE)
              count[j,k]<-1
              break
            }else{
              next
            } 
          }}
        for (l in 1:dim(count)[2]) {
          if(colSums(count)[l]==0){
            fsem1[["mod"]][[i]][[paste0("f",(length(fsem1[["mod"]][[i]])+1))]]<-fsem2[["mod"]][[1]][[l]]  
          }
        }
      }else{
        for (j in 1:length(fsem1[["mod"]][[i]])) {
          if(fsem1[["mod"]][[i]][[j]][[2]]==fsem2[["mod"]][[1]][[1]][[2]]){
            count2<-"TRUE"
            fsem1[["mod"]][[i]][[j]][["type"]]<-fsem1[["mod"]][[i]][[j]][[1]]
            fsem1[["mod"]][[i]][[j]][["response"]]<-fsem1[["mod"]][[i]][[j]][[2]]
            fsem1[["mod"]][[i]][[j]][["covariate"]]<-c(fsem1[["mod"]][[i]][[j]][[3]],fsem2[["mod"]][[1]][[1]][[3]])
            fsem1[["mod"]][[i]][[j]][["effect"]]<-c(fsem1[["mod"]][[i]][[j]][[4]],fsem2[["mod"]][[1]][[1]][[4]])
            fsem1[["mod"]][[i]][[j]][["intercept"]]<-ifelse(fsem1[["mod"]][[i]][[j]][["intercept"]]==TRUE,TRUE,FALSE)
            fsem1[["mod"]][[i]][[j]][["scalar.covariate"]]<-c(fsem1[["mod"]][[i]][[j]][["scalar.covariate"]],fsem2[["mod"]][[1]][[1]][["scalar.covariate"]])
            break
          }else{
            next
          }}
        if(count2=="FALSE"){
          cr<-length(fsem1[["mod"]][["regression"]])
          fsem1[["mod"]][["regression"]][[paste0("r",(cr+1))]]<-fsem2[["mod"]][[1]][[1]] 
        }}
    }
  }
  if(all(count1==1)){
    fsem1[["mod"]]=c(fsem1[["mod"]],fsem2[["mod"]])
  }
  fsem1[["var"]][["indicators"]]<-c(fsem1[["var"]][["indicators"]],fsem2[["var"]][["indicators"]])
  fsem1[["var"]][["factors"]]<-c(fsem1[["var"]][["factors"]],fsem2[["var"]][["factors"]])
  fsem1[["var"]][["latents"]]<-c(fsem1[["var"]][["factors"]],fsem1[["var"]][["latents"]],fsem2[["var"]][["latents"]])
  fsem1[["var"]][["observed"]]<-c(fsem1[["var"]][["observed"]],fsem2[["var"]][["observed"]])
  fsem1[["var"]][["indicators"]]<-unique(fsem1[["var"]][["indicators"]])
  fsem1[["var"]][["factors"]]<-unique(fsem1[["var"]][["factors"]])
  fsem1[["var"]][["latents"]]<-unique(fsem1[["var"]][["latents"]])
  fsem1[["var"]][["observed"]]<-unique(fsem1[["var"]][["observed"]])
  fsem1.sort<-sort.fsem(fsem1)
  aa<-extraInf(fsem1.sort)
  a<-aa
  a
}

sort.fsem<-function(model){
  lr=length(model$mod$regression)
  if(lr<2)
    return(model)
  
  dep.mat=matrix(0,lr,lr) #Matrix of dependencies
  for(i in 1:lr){
    for(j in 1:lr){
      if(i==j)
        next
      if(model$mod$regression[[i]]$response%in%model$mod$regression[[j]]$covariate)
        dep.mat[i,j]=1 
    }
  }
  rownames(dep.mat)=colnames(dep.mat)=paste0("r",1:lr)
  
  which.are.zero<-function(mat){
    which(colSums(mat)==0)
  }
  
  last.gen.lev=Inf
  
  inds=c()
  for(i in 1:(lr-1)){
    
    wz=which.are.zero(dep.mat)
    
    if(length(wz)==0)
      stop("Model is not consistent! There is a loop within the model")
    
    nms=colnames(dep.mat)[wz]
    ind.condidate=wz[1]
    inds[i]=nms[1]
    
    if(i<(lr-1)){
      dep.mat=dep.mat[-ind.condidate,-ind.condidate]
    }else{
      inds[i+1]=colnames(dep.mat)[-ind.condidate]
    }
    
  }
  
  regs=model$mod$regression
  model$mod$regression=list()
  for(i in 1:lr)
    model$mod$regression[[paste0("r",i)]]=regs[[inds[i]]]
  
  model
}


extraInf<-function(model){
  lr=length(model$mod$regression)
  no.obs<-length(model$var$observed)
  fr<-lr
  if(no.obs>0){
    model$var$observedScalar<-data.frame(observed=model$var$observed,scalar=rep(NA,length(model$var$observed)))
    for (i in 1:fr) {
      a<-which(model$var$observedScalar$observed%in%model$mod$regression[[i]]$covariate==TRUE)
      if(length(a)>0)
        model$var$observedScalar$scalar[a]<-model$mod$regression[[i]]$scalar[which(model$mod$regression[[i]]$covariate%in%model$var$observedScalar$observed==TRUE)]
    }
  }
  model
}

"%+%"<-GenModel


###parameters############################################
isolate.fun<-function(func,lst.vars){
  b=new.env()
  nms=names(lst.vars)
  for(n in nms)
    b[[n]]=lst.vars[[n]]
  
  environment(func)=b
  func
}

parameter1<-function(model,r,rho,n.b.sim,SNR,Matern.fac=FALSE,Matern.sem=FALSE,range.min,range.max){
  if(is.null(range.min))
    range.min=0
  if(is.null(range.max))
    range.max=1
  
  if(is.null(n.b.sim))
    n.b.sim=100
  range<-seq(range.min,range.max,length.out=100)
  J<-n.b.sim
  k<-1:J
  p<-list()
  for (m in 1:length(model[["mod"]])) {
    l<-length(model[["mod"]][[m]])
    if(model[["mod"]][[m]][[1]][["type"]]=="factorModel"){
      for (i in 1:l) {
        lf<-length(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect)
        if(Matern.fac==FALSE){
          p[["factorModel"]][[paste0("f",i)]][["eigen"]][["value"]]<-r*(rho^(k-1))/m
          p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]]<-list()
          ker<-0
          for (n in 1:J) {
            if(k[n]%%2==0){
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*sin(pi*n*i*m*t)},list(i=i,m=m,n=n)) 
            }else{
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*cos(pi*n*i*m*t)},list(i=i,m=m,n=n)) 
            } 
            ker<-ker+p[["factorModel"]][[paste0("f",i)]][["eigen"]][["value"]][n]*
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]](range)^2
          }   
        }else{
          p[["factorModel"]][[paste0("f",i)]]$Matern<-isolate.fun(function(s,t){exp(-1*abs(t-s)*i*m)},list(i=i,m=m)) 
          ker<-sapply(range, function(ii){p[["factorModel"]][[paste0("f",i)]]$Matern(ii,range)}) 
        }
        signal<-sum(sapply(1:length(ker),function(ii){1/length(ker)*(sqrt(ker[ii]))}))
        snr<-SNR
        p[["factorModel"]][[paste0("f",i)]]$error.std<-signal/snr
        p[["factorModel"]][[paste0("f",i)]][["coefficient"]]<-list()
        if(model[["mod"]][["factorModel"]][[paste0("f",i)]][["intercept"]]==FALSE){
          p[["factorModel"]][[paste0("f",i)]][["intercept"]]<-0
        }else{
          p[["factorModel"]][[paste0("f",i)]][["intercept"]]<-isolate.fun(function(t){t^2*i*m},list(i=i,m=m)) 
        }
        if(lf>0)
          for (j in 1:lf){
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="concurrent"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1+1/2*sin(pi*t*sqrt(i)*j/2*m)},list(i=i,j=j,m=m))
            }
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="historical"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*cos(pi*(s+t)*sqrt(i)*j/m)},list(i=i,j=j,m=m))
            }
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="fixed"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]=1
            }
          }
        
      }}else{
        for (i in 1:l) {
          res<-model$mod$regression[[i]]$response
          lat.clount<-which(model$var$latents==res)
          lf<-length(model[["mod"]][["regression"]][[paste0("r",i)]]$covariate)
          if(Matern.sem==FALSE){
            p[["regression"]][[paste0("r",i)]][["eigen"]][["value"]]<-r*(rho^(k-1))*lat.clount/m
            p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]]<-list()
            for (n in 1:J) {
              if(k[n]%%2==0){
                p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*sin(pi*n*lat.clount*m*t)},list(lat.clount=lat.clount,m=m,n=n)) 
              }else{
                p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*cos(pi*n*lat.clount*m*t)},list(lat.clount=lat.clount,m=m,n=n))  
              }  
            } 
          }else{
            p[["regression"]][[paste0("r",i)]]$Matern<-isolate.fun(function(s,t){exp(-1*(t-s)^2*lat.clount*m)},list(lat.clount=lat.clount,m=m)) 
          }
          p[["regression"]][[paste0("r",i)]][["coefficient"]]<-list()
          if(model[["mod"]][["regression"]][[paste0("r",i)]][["intercept"]]==FALSE){
            p[["regression"]][[paste0("r",i)]][["intercept"]]<-0 
          }else{
            p[["regression"]][[paste0("r",i)]][["intercept"]]<-isolate.fun(function(t){t*lat.clount*m},list(lat.clount=lat.clount,m=m)) 
          }
          if(lf>0)
            for (j in 1:lf){
              
              if(model[["mod"]][["regression"]][[paste0("r",i)]]$covariate[j]%in%model$var$latents){
                
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="concurrent"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1+1/2*cos(pi*t*sqrt(lat.clount)*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                }
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="historical"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t)*sqrt(lat.clount)*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                }
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="fixed"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-NULL
                }
              }else{
                if(!(is.null(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]))){
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="constant"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(x){x^{2}*j*lat.clount},list(lat.clount=lat.clount,j=j))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="linear"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){t^{2}*j*lat.clount},list(lat.clount=lat.clount,j=j))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="concurrent"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1-1/2*cos(t*lat.clount*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="historical"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t^2)*lat.clount*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                  }  
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="smooth"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(x,t){x*t^2*lat.clount*j},list(lat.clount=lat.clount,j=j))
                  }
                }}
            } 
        }
      } 
  }
  p
}
parameter2<-function(model,r,rho,n.b.sim,SNR,Matern.fac=FALSE,Matern.sem=FALSE,range.min,range.max){
  if(is.null(range.min))
    range.min=0
  if(is.null(range.max))
    range.max=1
  
  if(is.null(n.b.sim))
    n.b.sim=100
  range<-seq(range.min,range.max,length.out=100)
  J<-n.b.sim
  k<-1:J
  p<-list()
  for (m in 1:length(model[["mod"]])) {
    l<-length(model[["mod"]][[m]])
    if(model[["mod"]][[m]][[1]][["type"]]=="factorModel"){
      for (i in 1:l) {
        lf<-length(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect)
        if(Matern.fac==FALSE){
          p[["factorModel"]][[paste0("f",i)]][["eigen"]][["value"]]<-r*(rho^(k-1))/m
          p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]]<-list()
          ker<-0
          for (n in 1:J) {
            if(k[n]%%2==0){
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*sin(pi*n*i*m*t)},list(i=i,m=m,n=n)) 
            }else{
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*cos(pi*n*i*m*t)},list(i=i,m=m,n=n)) 
            } 
            ker<-ker+p[["factorModel"]][[paste0("f",i)]][["eigen"]][["value"]][n]*
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]](range)^2
          }   
        }else{
          p[["factorModel"]][[paste0("f",i)]]$Matern<-isolate.fun(function(s,t){exp(-1*abs(t-s)*i*m)},list(i=i,m=m)) 
          ker<-sapply(range, function(ii){p[["factorModel"]][[paste0("f",i)]]$Matern(ii,range)}) 
        }
        signal<-sum(sapply(1:length(ker),function(ii){1/length(ker)*(sqrt(ker[ii]))}))
        snr<-SNR
        p[["factorModel"]][[paste0("f",i)]]$error.std<-signal/snr
        p[["factorModel"]][[paste0("f",i)]][["coefficient"]]<-list()
        if(model[["mod"]][["factorModel"]][[paste0("f",i)]][["intercept"]]==FALSE){
          p[["factorModel"]][[paste0("f",i)]][["intercept"]]<-0
        }else{
          p[["factorModel"]][[paste0("f",i)]][["intercept"]]<-isolate.fun(function(t){t*i*m},list(i=i,m=m)) 
        }
        if(lf>0)
          for (j in 1:lf){
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="concurrent"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1+1/2*cos(pi*t*sqrt(i)*j/m)},list(i=i,j=j,m=m))
            }
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="historical"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t)*sqrt(i)*j/m)},list(i=i,j=j,m=m))
            }
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="fixed"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]=1
            }
          }
        
      }}else{
        for (i in 1:l) {
          res<-model$mod$regression[[i]]$response
          lat.clount<-which(model$var$latents==res)
          lf<-length(model[["mod"]][["regression"]][[paste0("r",i)]]$covariate)
          if(Matern.sem==FALSE){
            p[["regression"]][[paste0("r",i)]][["eigen"]][["value"]]<-r*(rho^(k-1))*lat.clount/m
            p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]]<-list()
            for (n in 1:J) {
              if(k[n]%%2==0){
                p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*sin(pi*n*lat.clount*m*t)},list(lat.clount=lat.clount,m=m,n=n)) 
              }else{
                p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*cos(pi*n*lat.clount*m*t)},list(lat.clount=lat.clount,m=m,n=n))  
              }  
            } 
          }else{
            p[["regression"]][[paste0("r",i)]]$Matern<-isolate.fun(function(s,t){exp(-1*abs(t-s)*lat.clount*m)},list(lat.clount=lat.clount,m=m)) 
          }
          p[["regression"]][[paste0("r",i)]][["coefficient"]]<-list()
          if(model[["mod"]][["regression"]][[paste0("r",i)]][["intercept"]]==FALSE){
            p[["regression"]][[paste0("r",i)]][["intercept"]]<-0 
          }else{
            p[["regression"]][[paste0("r",i)]][["intercept"]]<-isolate.fun(function(t){t*lat.clount*m},list(lat.clount=lat.clount,m=m)) 
          }
          if(lf>0)
            for (j in 1:lf){
              
              if(model[["mod"]][["regression"]][[paste0("r",i)]]$covariate[j]%in%model$var$latents){
                
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="concurrent"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1+1/2*cos(pi*t*sqrt(lat.clount)*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                }
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="historical"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t)*sqrt(lat.clount)*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                }
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="fixed"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-NULL
                }
              }else{
                if(!(is.null(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]))){
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="constant"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(x){x^{2}*j*lat.clount},list(lat.clount=lat.clount,j=j))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="linear"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){t^{2}*j*lat.clount},list(lat.clount=lat.clount,j=j))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="concurrent"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1-1/2*cos(t*lat.clount*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="historical"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t^2)*lat.clount*j/m)},list(lat.clount=lat.clount,j=j,m=m))
                  }  
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="smooth"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(x,t){x*t^2*lat.clount*j},list(lat.clount=lat.clount,j=j))
                  }
                }}
            } 
        }
      } 
  }
  p
}
parameter3<-function(model,r,rho,n.b.sim,SNR,Matern.fac=FALSE,Matern.sem=FALSE,range.min,range.max){
  if(is.null(range.min))
    range.min=0
  if(is.null(range.max))
    range.max=1
  
  if(is.null(n.b.sim))
    n.b.sim=100
  range<-seq(range.min,range.max,length.out=100)
  J<-n.b.sim
  k<-1:J
  p<-list()
  for (m in 1:length(model[["mod"]])) {
    l<-length(model[["mod"]][[m]])
    if(model[["mod"]][[m]][[1]][["type"]]=="factorModel"){
      for (i in 1:l) {
        lf<-length(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect)
        if(Matern.fac==FALSE){
          p[["factorModel"]][[paste0("f",i)]][["eigen"]][["value"]]<-r*(rho^(k-1))/m
          p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]]<-list()
          ker<-0
          for (n in 1:J) {
            if(k[n]%%2==0){
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*sin(pi*n*i*m*t)},list(i=i,m=m,n=n))
            }else{
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*cos(pi*n*i*m*t)},list(i=i,m=m,n=n))
            }
            ker<-ker+p[["factorModel"]][[paste0("f",i)]][["eigen"]][["value"]][n]*
              p[["factorModel"]][[paste0("f",i)]][["eigen"]][["function"]][[n]](range)^2
          }
        }else{
          p[["factorModel"]][[paste0("f",i)]]$Matern<-isolate.fun(function(s,t){exp(-1*abs(t-s)*i*m)},list(i=i,m=m))
          ker<-sapply(range, function(ii){p[["factorModel"]][[paste0("f",i)]]$Matern(ii,range)})
        }
        signal<-sum(sapply(1:length(ker),function(ii){1/length(ker)*(sqrt(ker[ii]))}))
        snr<-SNR
        p[["factorModel"]][[paste0("f",i)]]$error.std<-signal/snr
        p[["factorModel"]][[paste0("f",i)]][["coefficient"]]<-list()
        if(model[["mod"]][["factorModel"]][[paste0("f",i)]][["intercept"]]==FALSE){
          p[["factorModel"]][[paste0("f",i)]][["intercept"]]<-0
        }else{
          p[["factorModel"]][[paste0("f",i)]][["intercept"]]<-isolate.fun(function(t){t*i*m},list(i=i,m=m))
        }
        if(lf>0)
          for (j in 1:lf){
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="concurrent"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1+1/2*cos(pi*t*sqrt(i)*j/m)},list(i=i,j=j,m=m))
            }
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="historical"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t)*sqrt(i)*j/m)},list(i=i,j=j,m=m))
            }
            if(model[["mod"]][["factorModel"]][[paste0("f",i)]]$effect[j]=="fixed"){
              p[["factorModel"]][[paste0("f",i)]][["coefficient"]][[j]]=1
            }
          }
        
      }}else{
        for (i in 1:l) {
          lf<-length(model[["mod"]][["regression"]][[paste0("r",i)]]$covariate)
          if(Matern.sem==FALSE){
            p[["regression"]][[paste0("r",i)]][["eigen"]][["value"]]<-r*(rho^(k-1))*i/m
            p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]]<-list()
            for (n in 1:J) {
              if(k[n]%%2==0){
                p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*sin(pi*n*i*m*t)},list(i=i,m=m,n=n))
              }else{
                p[["regression"]][[paste0("r",i)]][["eigen"]][["function"]][[n]]<-isolate.fun(function(t){sqrt(2)*cos(pi*n*i*m*t)},list(i=i,m=m,n=n))
              }
            }
          }else{
            p[["regression"]][[paste0("r",i)]]$Matern<-isolate.fun(function(s,t){exp(-1*abs(t-s)*i*m)},list(i=i,m=m))
          }
          p[["regression"]][[paste0("r",i)]][["coefficient"]]<-list()
          if(model[["mod"]][["regression"]][[paste0("r",i)]][["intercept"]]==FALSE){
            p[["regression"]][[paste0("r",i)]][["intercept"]]<-0
          }else{
            p[["regression"]][[paste0("r",i)]][["intercept"]]<-isolate.fun(function(t){t*i*m},list(i=i,m=m))
          }
          if(lf>0)
            for (j in 1:lf){
              
              if(model[["mod"]][["regression"]][[paste0("r",i)]]$covariate[j]%in%model$var$latents){
                
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="concurrent"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1+1/2*cos(pi*t*sqrt(i)*j/m)},list(i=i,j=j,m=m))
                }
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="historical"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t)*sqrt(i)*j/m)},list(i=i,j=j,m=m))
                }
                if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="fixed"){
                  p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-NULL
                }
              }else{
                if(!(is.null(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]))){
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="constant"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(x){x^{2}*j*i},list(i=i,j=j))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="linear"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){t^{2}*j*lat.clount},list(lat.clount=lat.clount,j=j))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="concurrent"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(t){1+1/2*cos(t*i*j/m)},list(i=i,j=j,m=m))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="historical"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(s,t){1+1/2*sin(pi*(s+t^2)*i*j/m)},list(i=i,j=j,m=m))
                  }
                  if(model[["mod"]][["regression"]][[paste0("r",i)]]$effect[j]=="smooth"){
                    p[["regression"]][[paste0("r",i)]][["coefficient"]][[j]]<-isolate.fun(function(x,t){x*t^2*i*j},list(i=i,j=j))
                  }
                }}
            }
        }
      }
  }
  p
}

####simulation############################################

simulation<-function(model,n.sample,n.t,n.b.sim=100,n.b,r=0.5,rho=0.1,SNR=1,
                     design=c("regular", "irregular", "regular.truncated", "regular.missing"),
                     Matern.fac=FALSE,Matern.sem=TRUE,range.min=NULL,range.max=NULL,n.Rimanian=NULL,
                     parameters = c("parameter_set1","parameter_set3","parameter_set2")){
  
  design = match.arg(design)
  paramters = match.arg(parameters)
  
  params.org = list(
    n.b.sim=n.b.sim,r=r,rho=rho,SNR=SNR,design=design,
    Matern.fac=Matern.fac,Matern.sem=Matern.sem,range.min=range.min,range.max=range.max,n.Rimanian=n.Rimanian,
    parameters=parameters
  )
  
  
  #set of functional regression coefficients parameters and other extra parameters of the model
  parameter =  if(parameters=="parameter_set1"){
    #Defined in Core/utilities.R
    parameter1
  }else if(parameters=="parameter_set2"){
    #Defined in Core/utilities.R
    parameter2
  }else{
    #Defined in Core/utilities.R
    parameter3
  }
  
  if(is.null(range.min))
    range.min=0
  if(is.null(range.max))
    range.max=1
  
  parameters<-parameter(model=model,r=r,rho=rho,n.b.sim=n.b.sim,SNR=SNR,Matern.fac=Matern.fac,Matern.sem=Matern.sem,range.min=range.min,range.max=range.max) 
  if(is.null(n.Rimanian))
    n.Rimanian<-200
  
  basis<-create.bspline.basis(nbasis=n.b.sim,rangeval = c(range.min,range.max))
  n.tt<-max(100,2*n.b.sim)
  tt<-seq(range.min,range.max,length.out=n.tt)
  eval<-t(eval.basis(tt,basis))
  Nr<-n.Rimanian
  tM<-list()
  M<-list()
  no.fac<-length(model$mod$factorModel)
  if(design=="regular"){
    for (i in 1:n.sample) {
      tM[[i]]<-list()
      M[[i]]<-rep(n.t,no.fac)
      tM[[i]]<-lapply(1:no.fac,function(j){tM[[i]][[j]]<-seq(range.min,range.max,length.out=M[[i]][[j]])}) 
    }
  }
  if(design=="irregular"){
    kk<-1/(1-(1/2)^n.t) 
    q<-kk*(1/2)^(n.t-1:n.t+1)
    for (i in 1:n.sample) {
      tM[[i]]<-list()
      M[[i]]<-sample(x =1:n.t,no.fac, replace = T, prob = q)
      tM[[i]]<-lapply(1:no.fac,function(j){tM[[i]][[j]]<-sort(runif(M[[i]][j],range.min,range.max))}) 
    }
  }
  if(design=="regular.truncated"){
    kk<-1/(1-(1/2)^n.t) 
    q<-kk*(1/2)^(n.t-1:n.t+1)
    for (i in 1:n.sample) {
      tM[[i]]<-list()
      M[[i]]<-sample(x =1:n.t,no.fac, replace = T, prob = q)
      tM[[i]]<-lapply(1:no.fac,function(j){tM[[i]][[j]]<-seq(range.min,range.max,length.out=M[[i]][j])}) 
    }
  }
  if(design=="regular.missing"){
    kk<-1/(1-(1/2)^n.t) 
    q<-kk*(1/2)^(n.t-1:n.t+1)
    ex<-sum((1:n.t)*q)
    pp<-1-(ex/n.t)
    for (i in 1:n.sample) {
      tM[[i]]<-list()
      tM[[i]]<-lapply(1:no.fac,function(j){tM[[i]][[j]]<-sort(tt[which(sapply(1:n.t, function(ii){rbinom(1,1,(1-pp))})==1)])})
      M[[i]]<-sapply(1:no.fac,function(j){length(tM[[i]][[j]])}) 
    }
  }
  J<-n.b.sim
  no_r<-length(model$mod$regression)
  no_f<-no.fac
  r<-list()
  f<-list()
  ff<-list()
  if(!is.null(model$var$observed)){
    x.generator<-list()
    for (kk in 1:length(model$var$observed)) {
      nx<-n.sample
      x.generator[[model$var$observed[kk]]]<-rnorm(nx,0,1)
    }
    x.data<-list()
    for (i in 1:length(model$var$observed)) {
      nx<-n.sample
      cov<-model$var$observedScalar$observed[i]
      cov.overwrite.later = FALSE
      if(is.na(model$var$observedScalar$scalar[i])&& cov%in%model$var$latents)
        model$var$observedScalar$scalar[i] = FALSE
      
      if(cov%in%model$var$latents) cov.overwrite.later = TRUE
         
      if(model$var$observedScalar$scalar[i]==TRUE){
        x.data[[cov]]<-x.generator[[cov]]  
      }else{
        x.data[[cov]]<-list()
        attr(x.data[[cov]], "cov.overwrite.later") = cov.overwrite.later
        for (ii in 1:nx) {
          x.data[[cov]][[ii]]<-x.generator[[cov]][ii]*tt ##Will be overwritten if it is a latent variable
        }
      }
    }
  }
  
  if(no_r>0)
    for (kk in 1:no_r) {
      no_co<-length(model$mod$regression[[kk]][["covariate"]])
      res<-model$mod$regression[[kk]][["response"]]
      nn<-n.sample
      mm<-list()
      if(Matern.sem==TRUE){
        sigg<-sapply(tt, function(ii){parameters$regression[[kk]]$Matern(ii,tt)})
        for (n in 1:nn) {
          mm[[n]]<-mvrnorm(1,rep(0,n.tt),sigg)
        }
      }else{
        
        ker = 0
        for (j in 1:n.b.sim) {
          ker<-ker+parameters$regression[[kk]]$eigen$value[j]*outer(tt,tt,function(x,y)
            parameters$regression[[kk]]$eigen$`function`[[j]](x)*
              parameters$regression[[kk]]$eigen$`function`[[j]](y))
        }
        for (n in 1:nn) {
          mm[[n]]<-mvrnorm(1,rep(0,n.tt),ker)
        }
      }
      r[[res]]<-mm
      
      if(no_co>0)
        for (j in 1:no_co) {
          nx<-n.sample
          mx<-list()
          if(!(model$mod$regression[[kk]][["covariate"]][j]%in%model$var$latents)){
            x<-x.data[[model$mod$regression[[kk]]$covariate[j]]]
            if(model$mod$regression[[kk]][["effect"]][j]=="constant"){
              m<-parameters$regression[[kk]][["coefficient"]][[j]](x)
              for (i in 1:nx) {
                mx[[i]]<-rep(m[i],n.tt)  
              }
            }
            if(model$mod$regression[[kk]][["effect"]][j]=="smooth"){
              for (i in 1:nx) {
                m<-parameters$regression[[kk]][["coefficient"]][[j]](x[i],tt)
                mx[[i]]<-m}
            }
            if(model$mod$regression[[kk]][["effect"]][j]=="concurrent"){
              for (i in 1:nx) {
                xx<-x[[i]]
                m<-t(t(xx)*parameters$regression[[kk]][["coefficient"]][[j]](tt))
                mx[[i]]<-m
              }
            }
            if(model$mod$regression[[kk]][["effect"]][j]=="linear"){
              for (i in 1:nx) {
                m<-x[i]*parameters$regression[[kk]][["coefficient"]][[j]](tt)
                mx[[i]]<-m
              }
            }
            if(model$mod$regression[[kk]][["effect"]][j]=="historical"){
              for (j1 in 1:nx) {
                xx<-x[[i]]
                x.coef<-solve(eval%*%t(eval))%*%eval%*%xx
                hist<-sapply(tt,function(t){ind1=seq(range.min,t,length.out=n.Rimanian);evv<-t(eval.basis(ind1,basis));sum(sapply(1:length(ind1),function(s){1/Nr*parameters$factorModel[[k]][["coefficient"]][[j]](ind1[s],t)*(t(evv)[s,]%*%x.coef)}))})
                mx[[j1]]<-hist 
              }
            }
          }else{
            ww<-which(names(r)==model$mod$regression[[kk]][["covariate"]][j])
            if(model$mod$regression[[kk]][["effect"]][j]=="fixed")
              mx<-r[[ww]]
            
            if(model$mod$regression[[kk]][["effect"]][j]=="concurrent")
              for (j1 in 1:nx) {
                m<-t(t(r[[ww]][[j1]])*parameters$regression[[kk]][["coefficient"]][[j]](tt))
                mx[[j1]]<-m
              }
            
            if(model$mod$regression[[kk]][["effect"]][j]=="historical"){
              for (j1 in 1:nx) {
                et.coef<-solve(eval%*%t(eval))%*%eval%*%r[[ww]][[j1]]
                hist<-sapply(tt,function(t){ind1=seq(range.min,t,length.out=n.Rimanian);evv<-t(eval.basis(ind1,basis));sum(sapply(1:length(ind1),function(s){1/Nr*parameters$factorModel[[k]][["coefficient"]][[j]](ind1[s],t)*(t(evv)[s,]%*%et.coef)}))})
                mx[[j1]]<-hist 
              }
            }}
          for (i2 in 1:n.sample) {
            r[[res]][[i2]]<-t(t(as.matrix(r[[res]][[i2]]))+as.vector(mx[[i2]]))
            if(model$mod$regression[[kk]][["intercept"]]==TRUE){
              r[[res]][[i2]]<-t(t(r[[res]][[i2]])+parameters$regression[[kk]][["intercept"]](tt)) 
            }
          }
        }
    }
  ev<-list()
  for (k in 1:no_f) {
    ev[[k]]<-list()
    for (n in 1:n.sample) {
      e<-t(eval.basis(tM[[n]][[k]],basis))
      ev[[k]][[n]]<-t(e)%*%solve(eval%*%t(eval))%*%eval
    }
  }
  for(k in 1:no_f){
    no_co<-length(model$mod$factorModel[[k]][["factor"]])
    res<-paste0("w",k)
    res2<-model$mod$factorModel[[k]][["indicator"]]
    nn<-n.sample
    mm<-list()
    err<-list()
    if(Matern.fac==TRUE){
      sM<-tM
      for (n in 1:nn) {
        sigg<-sapply(sM[[n]][[k]], function(ii){parameters$factorModel[[k]]$Matern(ii,tM[[n]][[k]])})
        mm[[n]]<-mvrnorm(1,rep(0,M[[n]][[k]]),sigg)
        err[[n]]<-mvrnorm(1,rep(0,M[[n]][[k]]),parameters$factorModel[[k]]$error.std^2*diag(M[[n]][[k]]))
      }
    }else{
      
      for (n in 1:nn) {
        ker = 0
        for (j in 1:n.b.sim) {
          ker<-ker+parameters$factorModel[[k]]$eigen$value[j]*outer(tM[[n]][[k]],tM[[n]][[k]],function(x,y)
            parameters$factorModel[[k]]$eigen$`function`[[j]](x)*
              parameters$factorModel[[k]]$eigen$`function`[[j]](y))
        }
        mm[[n]]<-mvrnorm(1,rep(0,M[[n]][[k]]),ker)
        err[[n]]<-mvrnorm(1,rep(0,M[[n]][[k]]),parameters$factorModel[[k]]$error.std^2*diag(M[[n]][[k]]))
      }
    }
    ff[[res2]]<-lapply(1:nn,function(ii){mm[[ii]]})
    f[[res]]<-lapply(1:nn,function(ii){mm[[ii]]+err[[ii]]})
    for (j in 1:no_co) {
      ww<-which(names(r)==model$mod$factorModel[[k]][["factor"]][j])
      nx<-n.sample
      mx<-list()
      ef<-model$mod$factorModel[[k]][["effect"]][j]
      fac<-model$mod$factorModel[[k]][["factor"]][j]
      if(ef=="concurrent"){
        for (j1 in 1:nx) {
          mx[[j1]]<-t(t(ev[[k]][[j1]]%*%r[[fac]][[j1]])*parameters$factorModel[[k]][["coefficient"]][[j]](tM[[j1]][[k]]))
        }
      }
      if(ef=="historical"){
        for (j1 in 1:nx) {
          et.coef<-solve(eval%*%t(eval))%*%eval%*%r[[ww]][[j1]]
          hist<-sapply(tM[[j1]][[k]],function(t){ind1=seq(range.min,t,length.out=n.Rimanian);evv<-t(eval.basis(ind1,basis));sum(sapply(1:length(ind1),function(s){1/Nr*parameters$factorModel[[k]][["coefficient"]][[j]](ind1[s],t)*(t(evv)[s,]%*%et.coef)}))})
          mx[[j1]]<-hist 
        }
      }
      if(ef=="fixed_concurrent"){
        for (j1 in 1:nx) {
          mx[[j1]]<-ev[[k]][[j1]]%*%r[[fac]][[j1]]
        } 
      }
      if(ef=="fixed_historical"){
        for (j1 in 1:nx) {
          et.coef<-solve(eval%*%t(eval))%*%eval%*%r[[ww]][[j1]]
          hist<-sapply(tM[[j1]][[k]],function(t){ind1=seq(range.min,t,length.out=n.Rimanian);evv<-t(eval.basis(ind1,basis));sum(sapply(1:length(ind1),function(s){1/Nr*(t(evv)[s,]%*%et.coef)}))})
          mx[[j1]]<-hist 
        } 
      }
      for (i3 in 1:n.sample) {
        f[[res]][[i3]]<-t(t(as.matrix(f[[res]][[i3]]))+as.vector(mx[[i3]]))
        ff[[res2]][[i3]]<-t(t(as.matrix(ff[[res2]][[i3]]))+as.vector(mx[[i3]]))
        if(model$mod$factorModel[[k]][["intercept"]]==TRUE){
          f[[res]][[i3]]<-t(t(f[[res]][[i3]])+parameters$factorModel[[k]][["intercept"]](tM[[i3]][[k]])) 
          ff[[res2]][[i3]]<-t(t(ff[[res2]][[i3]])+parameters$factorModel[[k]][["intercept"]](tM[[i3]][[k]])) 
        }
      }
    }
  }
  data2<-NULL
  for (k in 1:no_f) {
    res2<-model$mod$factorModel[[k]][["indicator"]]
    for (ii in 1:n.sample) {
      data2<-rbind(data2,data.frame(.id=ii,.ind=k,.t=tM[[ii]][[k]],.value=as.vector(f[[paste0("w",k)]][[ii]]),.value.ind=as.vector(ff[[res2]][[ii]])))
    }
  }
  data.et<-NULL
  for (k in 1:no_r) {
    res<-model$mod$regression[[k]][["response"]]
    lat.count<-which(model$var$latents==model$mod$regression[[k]]$response)
    for (ii in 1:n.sample) {
      data.et<- rbind(data.et,data.frame(.id=ii,.et=lat.count,.value=as.vector(ev[[k]][[n]]%*%r[[res]][[ii]]))) 
    }
  }
  data1<-if(!is.null(model$var$observed)){list(data=data2,data.et=data.et,covariate=x.data)}else{list(data=data2,data.et=data.et)}
  
  
  params.org.eval<-parameter.org.evaluted(model=model,
                                          r=params.org$r,rho=params.org$rho,n.b.sim=params.org$n.b.sim,n.b=n.b,
                                          range.min=params.org$range.min,range.max=params.org$range.max,SNR=params.org$SNR,
                                          x.data=data1$covariate,Matern.sem=params.org$Matern.sem,Matern.fac=params.org$Matern.fac,
                                          parameters = params.org$parameters)
  
  data1$params.org.eval = params.org.eval
  data1$params = list(
    model=model,n.sample=n.sample,n.t=n.t,n.b.sim=n.b.sim,r=r,rho=rho,SNR=SNR,
    design=design,
    Matern.fac=Matern.fac,Matern.sem=Matern.sem,range.min=range.min,
    range.max=range.max,n.Rimanian=n.Rimanian,
    parameters = parameters
  )
  
  data1
}  

################################parameter orginal#######################################
parameter.org.evaluted<-function(model,r,rho,n.b.sim,n.b,SNR,Matern.fac=FALSE,Matern.sem=FALSE,range.min=NULL,range.max=NULL,x.data=NULL,
                                 parameters = c("parameter_set1","parameter_set3","parameter_set2")){
  
  paramters = match.arg(parameters)
  
  #set of functional regression coefficients parameters and other extra parameters of the model
  parameter = if(parameters=="parameter_set1"){
    #Defined in Core/utilities.R
    parameter1
  }else if(parameters=="parameter_set2"){
    #Defined in Core/utilities.R
    parameter2
  }else{
    #Defined in Core/utilities.R
    parameter3
  }
  
  if(is.null(range.min))
    range.min=0
  if(is.null(range.max))
    range.max=1
  
  param<-parameter(model,r=r,rho=rho,n.b=n.b,SNR=SNR,Matern.fac=Matern.fac,
                   Matern.sem=Matern.sem, range.min=range.min,range.max=range.max)
  param2<-parameter(model,r=r,rho=rho,n.b=n.b.sim,SNR=SNR,Matern.fac=Matern.fac,
                    Matern.sem=Matern.sem, range.min=range.min,range.max=range.max)
  basis<-create.bspline.basis(nbasis=n.b,rangeval = c(range.min,range.max))
  n.t<-max(100,2*n.b.sim)
  M<-n.t
  t<-seq(range.min,range.max,length.out=M)
  eval<-t(eval.basis(seq(range.min,range.max,length.out=M),basis))
  ome<-inprod(basis,basis)
  ome.val<-eigen(ome)$values
  ome.vec<-eigen(ome)$vectors
  ome.half<-ome.vec%*%(diag(ome.val)^0.5)%*%t(ome.vec)
  no.reg<-length(model$mod$regression)
  init<-list()
  for (i in 1:length(param)) {
    if(model$mod[[i]][[1]]$type=="factorModel"){
      for (i1 in 1:length(param$factorModel)) {
        res<-model$mod$factorModel[[i1]]$indicator
        cov<-model$mod$factorModel[[i1]]$factor
        init$coef.fac[[res]]<-list()
        if(Matern.fac==TRUE){
          s<-t
          ker<-sapply(s,function(ii){param$factorModel[[i1]]$Matern(ii,t)}) 
        }else{
          ker = 0
          for (j in 1:n.b.sim) {
            ker<-ker+param2$factorModel[[i1]]$eigen$value[j]*outer(t,t,function(x,y)
              param2$factorModel[[i1]]$eigen$`function`[[j]](x)*
                param2$factorModel[[i1]]$eigen$`function`[[j]](y))
          }
        }
        init$ker.fac[[res]]<-ker
        init$ker.fac.std[[res]]<-ker
        d.fac1<-diag(diag(ker)^(-1/2))
        lat<-cov
        fac<-model$mod$factorModel[[i1]]$effect
        co<-which(sapply(1:no.reg, function(ii){model$mod$regression[[paste0("r",ii)]]$response==lat}))
        if(Matern.sem==TRUE){
          s<-t
          ker2<-sapply(s,function(ii){param$regression[[co]]$Matern(ii,t)}) 
        }else{
          ker2<-NULL
          for (tt in 1:M) {
            ktt<-0
            for (j in 1:n.b.sim) {
              ktt<-ktt+param2$regression[[co]]$eigen$value[j]*
                param2$regression[[co]]$eigen$`function`[[j]](t[tt])*
                param2$regression[[co]]$eigen$`function`[[j]](t)
            }
            ker2<-rbind(ker2,ktt)
          }
        }
        d.fac2<-diag(diag(ker2))
        if(fac=="concurrent"){
          init$coef.fac[[res]][[cov]]<-param$factorModel[[i1]]$coefficient[[1]](t)
          init$coef.fac.std[[res]][[cov]]<-d.fac2^(1/2)%*%
            param$factorModel[[i1]]$coefficient[[1]](t)
        }
        if(fac=="historical"){
          tt<-t
          evv<-t(eval.basis(tt,basis))
          lam.org<-matrix(NA,M,M)
          evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
          for (ii in 1:M) {
            la.org<-param$factorModel[[i1]]$coefficient[[1]](seq(0,tt[ii],length.out=M),tt[ii])
            lam.org[,ii]<-la.org
          }
          init$coef.fac[[res]][[cov]]<-lam.org
          init$coef.fac.std[[res]][[cov]]<-t(diag(d.fac2)^(1/2)*t(lam.org))
        }
        if(Matern.fac==TRUE){
          sigg<-(ome.half%*%solve(eval%*%t(eval))%*%eval%*%ker%*%t(eval)%*%solve(eval%*%t(eval))%*%ome.half)
          init$sigma.fac.eigval[[res]]<-eigen(sigg)$values
          init$sigma.fac.eigvec[[res]]<-(t(eval)%*%solve(ome)%*%ome.half%*%eigen(sigg)$vectors)
          init$sigma.fac.eigval.std[[res]]<-init$sigma.fac.eigval[[res]]
          init$sigma.fac.eigvec.std[[res]]<-init$sigma.fac.eigvec[[res]]
        }else{
          sigg<-(ome.half%*%solve(eval%*%t(eval))%*%eval%*%ker%*%t(eval)%*%solve(eval%*%t(eval))%*%ome.half)
          init$sigma.fac.eigval[[res]]<-eigen(sigg)$values
          init$sigma.fac.eigvec[[res]]<-(t(eval)%*%solve(ome)%*%ome.half%*%eigen(sigg)$vectors)
          init$sigma.fac.eigval.std[[res]]<-init$sigma.fac.eigval[[res]]
          init$sigma.fac.eigvec.std[[res]]<-init$sigma.fac.eigvec[[res]]
        }
        init$intercept[[res]]<-param$factorModel[[i1]]$intercept(t)
        init$sigma.error[[res]]<-param$factorModel[[i1]]$error.std^2
        init$intercept.std[[res]]<-init$intercept[[res]]
        init$sigma.error.std[[res]]<-init$sigma.error[[res]]
      }
    }else{
      for (i1 in 1:length(model$mod$regression)) {
        res<-model$mod$regression[[i1]]$response
        cov<-model$mod$regression[[i1]]$covariate
        no.cov<-length(cov)
        if(Matern.sem==TRUE){
          s<-t
          ker<-sapply(s,function(ii){param$regression[[i1]]$Matern(ii,t)}) 
        }else{
          ker<-NULL
          for (tt in 1:M) {
            ktt<-0
            for (j in 1:n.b.sim) {
              ktt<-ktt+param2$regression[[i1]]$eigen$value[j]*
                param2$regression[[i1]]$eigen$`function`[[j]](t[tt])*
                param2$regression[[i1]]$eigen$`function`[[j]](t)
            }
            ker<-rbind(ker,ktt)
          } 
        }
        d.fac1<-diag(diag(ker)^-(1/2))
        init$ker.sem[[res]]<-ker
        init$ker.sem.std[[res]]<-d.fac1%*%ker%*%d.fac1
        if(!is.null(cov))
          for (j in 1:no.cov) {
            fac<-model$mod$regression[[i1]]$effect[j]
            if(cov[j]%in%model$var$latents){
              cc<-which(sapply(1:no.reg, function(ii){model$mod$regression[[paste0("r",ii)]]$response==lat}))
              if(Matern.sem==TRUE){
                s<-t
                ker2<-sapply(s,function(ii){param$regression[[cc]]$Matern(ii,t)}) 
              }else{
                ker2<-NULL
                for (tt in 1:M) {
                  ktt<-0
                  for (jj in 1:n.b.sim) {
                    ktt<-ktt+param2$regression[[cc]]$eigen$value[jj]*
                      param2$regression[[cc]]$eigen$`function`[[jj]](t[tt])*
                      param2$regression[[cc]]$eigen$`function`[[jj]](t)
                  }
                  ker2<-rbind(ker2,ktt)
                } 
              }
              d.fac2<-diag(diag(ker2))
              if(fac=="concurrent"){
                init$coef.sem[[res]][[cov[j]]]<-param$regression[[i1]]$coefficient[[j]](t)
                init$coef.sem.std[[res]][[cov[j]]]<-d.fac1%*%d.fac2^(1/2)%*%
                  param$regression[[i1]]$coefficient[[j]](t)
              }
              if(fac=="historical"){
                tt<-t
                evv<-t(eval.basis(tt,basis))
                gam.org<-matrix(NA,M,M)
                evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
                for (ii in 1:M) {
                  la.org<-param$regression[[i1]]$coefficient[[1]](seq(0,tt[ii],length.out=M),tt[ii])
                  gam.org[,ii]<-la.org
                }
                par<-gam.org
                init$coef.sem[[res]][[cov[j]]]<-par
                init$coef.sem.std[[res]][[cov[j]]]<-t(diag(d.fac1)*diag(d.fac2)^(1/2)*t(par))
              }
            }else{
              if(fac=="concurrent"){
                init$coef.sem[[res]][[cov[j]]]<-param$regression[[i1]]$coefficient[[j]](t)
                init$coef.sem.std[[res]][[cov[j]]]<-d.fac1%*%
                  param$regression[[i1]]$coefficient[[j]](t)
              }
              if(fac=="historical"){
                tt<-t
                evv<-t(eval.basis(tt,basis))
                gam.org<-matrix(NA,M,M)
                evv_ind<-lapply(1:length(tt), function(ii){t(eval.basis(seq(0,tt[ii],length.out=M),basis))})
                for (ii in 1:M) {
                  la.org<-param$regression[[i1]]$coefficient[[j]](seq(0,tt[ii],length.out=M),tt[ii])
                  gam.org[,ii]<-la.org
                }
                par<-gam.org
                init$coef.sem[[res]][[cov[j]]]<-par
                init$coef.sem.std[[res]][[cov[j]]]<-t(diag(d.fac1)*t(par))
              }
              if(fac=="smooth"){
                par<-NULL
                for (j1 in 1:length(x.data[[cov[j]]])) {
                  par<-cbind(par,param$regression[[i1]]$coefficient[[j]](x.data[[cov[j]]][j1],t)) 
                }
                init$coef.sem[[res]][[cov[j]]]<-par
                init$coef.sem.std[[res]][[cov[j]]]<-diag(d.fac1)*par
              }
              if(fac=="linear"){
                init$coef.sem[[res]][[cov[j]]]<-param$regression[[i1]]$coefficient[[j]](t)
                init$coef.sem.std[[res]][[cov[j]]]<-d.fac1%*%
                  param$regression[[i1]]$coefficient[[j]](t)
              }
              if(fac=="constant"){
                init$coef.sem[[res]][[cov[j]]]<-param$regression[[i1]]$coefficient[[j]](x.data[[cov[j]]])
                init$coef.sem.std[[res]][[cov[j]]]<-d.fac1%*%rep(1,M)%*%
                  param$regression[[i1]]$coefficient[[j]](x.data[[cov[j]]])
              } 
            }
          }
        sig.sem<-solve(eval%*%t(eval))%*%eval%*%ker%*%t(eval)%*%solve(eval%*%t(eval))
        sig.sem.std<-solve(eval%*%t(eval))%*%
          eval%*%d.fac1%*%ker%*%d.fac1%*%t(eval)%*%solve(eval%*%t(eval))
        
        eig<-eigen(ome.half%*%sig.sem%*%ome.half)
        eig.std<-eigen(ome.half%*%sig.sem.std%*%ome.half)
        if(Matern.sem==TRUE){
          init$sigma.sem.eigvec[[res]]<-(t(eval)%*%solve(ome)%*%ome.half%*%
                                           eig$vectors)
          init$sigma.sem.eigval[[res]]<-eig$values
          init$sigma.sem.eigvec.std[[res]]<-(t(eval)%*%solve(ome)%*%ome.half%*%
                                               eig.std$vectors)
          init$sigma.sem.eigval.std[[res]]<-eig.std$values
        }else{
          init$sigma.sem.eigvec[[res]]<-(t(eval)%*%solve(ome)%*%ome.half%*%
                                           eig$vectors)
          init$sigma.sem.eigval[[res]]<-eig$values
          init$sigma.sem.eigvec.std[[res]]<-(t(eval)%*%solve(ome)%*%ome.half%*%
                                               eig.std$vectors)
          init$sigma.sem.eigval.std[[res]]<-eig.std$values
        }
      }
    }
  }
  init
} 


plot_lines = function(lbl, tt, y1, y2,add=FALSE, ...){
  sgn = ifelse(sum(y1*y2>0),1,-1)
  if(!add){
    plot(tt,y1,"l",col="red", ylab=lbl, ...)
    lines(tt,y2*sgn,col="blue", ...)
  }else{
    lines(tt,y1,"l",col="red", ...)
    lines(tt,y2*sgn,col="blue", ...)
  }
  
}

