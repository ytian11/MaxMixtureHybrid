## BSE_MS
BSE_MS <- function(y,S,knots=NULL,
                   X=NULL,Sp=NULL,Xp=NULL,
                   init.beta=c(0,0,.0001),
                   init.alpha=.4,
                   init.range=1,
                   init.bw=1,
                   vary=c(T,F,F),
                   pri.mn.range=-2,
                   pri.sd.range=1,
                   pri.mn.bw=0,
                   pri.sd.bw=1,
                   pri.var.a=0.1,
                   pri.var.b=0.1,
                   pri.alpha.a=1,
                   pri.alpha.b=1,
                   pri.sd.beta=10,
                   pri.mn.gev=c(0,0,0),
                   pri.sd.gev=c(10,1,0.25),
                   keep.samples=T,
                   iters=iters,burn=burn,
                   update=50,nthin=1){
  #MODEL FITTING
  
  start.time <- proc.time()
  library(fields)
  library(emulator)
  
  if(is.null(knots)){knots<-S}
  if(is.null(X)){X<-cbind(1,S)}
  
  
  #BOOKEEPING
  X<-as.matrix(X)
  n<-nrow(y)
  np<-nrow(Sp)
  nt<-ncol(y)
  nF<-nrow(knots)
  p<-ncol(X)
  
  d<-rdist(S,S)
  diag(d)<-0
  dw2<-rdist(S,knots)^2
  dw2[dw2<0.0001]<-0
  
  #INITIAL VALUES:
  
  beta<-matrix(init.beta,n,3,byrow=T)
  if(vary[1]){for(j in 1:n){
    beta[j,1]<-get.inits.mu(y[j,],beta[j,2],beta[j,3])
  }}
  mnb<-matrix(0,p,3)
  mnb[1,]<-colMeans(beta)
  taub<-rep(1,3)
  lograngeb<-log(init.range)
  logrange<-log(init.bw)
  lambda.ps <- matrix(2, ncol=nt, nrow=nF)
  u.ps<-matrix(0.5, ncol=nt, nrow=nF)
  u<-matrix(0.5,nF,nt)
  alpha<-init.alpha
  q1<-0.8
  logs<-exp(lambda.ps)
  
  
  #COMPUTE QUANTITIES USED IN THE POSTERIOR
  Qb<-solve(exp(-d/exp(lograngeb)))
  fac<-make.kern(dw2,logrange)
  FAC<-stdKern(fac)
  curll<-matrix(0,n,nt)
  theta<-make.theta(FAC,logs,alpha)
  for(j in 1:n){
    curll[j,]<-loglike(y[j,],beta[j,1],beta[j,2],beta[j,3],theta[j,],alpha)
  }
  
  
  #CREATE PLACES TO STORE THE OUTPUT
  keepers.bp<-array(0,c(iters,np,3))     # PREDICTED GEV PARAMETERS
  keepers.thetap<-array(0,c(iters,np,nt))# PREDICTED THETA
  keepers.FACp<-array(0,c(iters,np,nF))  # PREDICTED FAC
  keepers.others<-matrix(0,iters,3)      # ALPHA,Q1 AND LOGLIKELIHOOD
  
  colnames(keepers.others)<-c("alpha","q","log likelihood")
  GEVs<-c("GEV location","GEV log scale","GEV shape")
  SITEs<-paste("Site",1:n)
  beta.mn<-beta.var<-0*beta
  colnames(beta.mn)<-colnames(beta.var)<-GEVs
  rownames(beta.mn)<-rownames(beta.var)<-SITEs
  locs=NULL
  if(keep.samples){
    locs<-array(0,c(iters,n,3))
    dimnames(locs)<-list(paste("sample",1:iters),SITEs,GEVs)
  }
  
  Yp<-Yp1<-Yp2<-locsp<-beta.mn.p<-beta.var.p<-NULL
  if(!is.null(Sp)){
    SITEs<-paste("Pred site",1:np)
    TIMEs<-paste("Replication",1:nt)
    Yp1<-Yp2<-matrix(0,np,nt)
    llp<-matrix(0,np,nt)
    colnames(Yp1)<-colnames(Yp2)<-TIMEs
    rownames(Yp1)<-rownames(Yp2)<-SITEs
    beta.mn.p<-beta.var.p<-matrix(0,np,3)   
    colnames(beta.mn.p)<-colnames(beta.var.p)<-GEVs
    rownames(beta.mn.p)<-rownames(beta.var.p)<-SITEs
    if(keep.samples){
      locsp<-array(0,c(iters,np,3))
      dimnames(locsp)<-list(paste("sample",1:iters),SITEs,GEVs)
      Yp<-array(0,c(iters,np,nt))
      dimnames(Yp)<-list(paste("sample",1:iters),SITEs,TIMEs)
    }
    dw2p<-rdist(Sp,knots)^2
    dw2p[dw2p<0.0001]<-0
    if(sum(vary)>0){
      d12<-rdist(Sp,S)
      d22<-rdist(Sp,Sp)
      d12[d12<0.0001]<-0
      diag(d22)<-0
    }
  }
  
  #SETTINGS FOR THE M-H CANDIDATE DISTRIBUTION
  attb<-accb<-MHb<-c(0.1,0.02,0.02,0.01,0.02,0.02,1)
  MHu<-3
  MHq<-0.25
  attq<-accq<-1
  atts<-accs<-MHs<-c(3,2,2,2,rep(1,10))
  attp<-accp<-MHp<-c(3,2,2,2,rep(1,10))
  cuts<-seq(0,15,2)
  
  
  #START SAMPLING!
  level.ps<-olds.ps<-lambda.ps
  
  for(i in 1:iters){for(rep in 1:nthin){
    olds.ps<-lambda.ps
    
    ##########################################################
    ##############      Random effects S and U    ############
    ##########################################################
    for(t in 1:nt){
      for(l in 1:nF){
        level.ps[l,t]<-get.level(lambda.ps[l,t],cuts)
      }
    }
    ## update ps
    for(t in 1:nt){ 
      v<-beta[,3]*exp(-beta[,2])*(y[,t]-beta[,1])+1
      v<-v^(-1/(alpha*beta[,3]))      
      ccc<-logd(theta[,t],v)
      for(l in 1:nF){
        W<-FAC[,l]^(1/alpha)
        MH1<-MHp[level.ps[l,t]]
        canlambda<-rnorm(1,lambda.ps[l,t],MH1)
        MH2<-MHp[get.level(canlambda,cuts)]
        
        canlogs<-exp(canlambda)
        cantheta<-theta[,t]+W*(canlogs-logs[l,t])
        canccc<-logd(cantheta,v)
        R<-canccc-ccc+
          h1(canlambda,u.ps[l,t],alpha,log=T)-
          h1(lambda.ps[l,t],u.ps[l,t],alpha,log=T)+
          dnorm(lambda.ps[l,t],canlambda,MH2,log=T)-
          dnorm(canlambda,lambda.ps[l,t],MH1,log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          lambda.ps[l,t]<-canlambda  # A_l
          logs[l,t]<-canlogs
          ccc<-canccc
          theta[,t]<-cantheta # theta
        }}
      }
      curll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],theta[,t],alpha)
    }
    
    for(j in 1:length(MHp)){
      accp[j]<-accp[j]+sum(olds.ps[level.ps==j]!=lambda.ps[level.ps==j])
      attp[j]<-attp[j]+sum(level.ps==j)
    }
    for(j in 1:length(attp)){if(i<burn/2 & attp[j]>50){
      if(accp[j]/attp[j]<0.3){MHp[j]<-MHp[j]*0.9}
      if(accp[j]/attp[j]>0.6){MHp[j]<-MHp[j]*1.1}
      accp[j]<-attp[j]<-0
    }}
    
    canu<-rtnorm(u.ps)
    R<-h1(lambda.ps,canu,alpha,log=T)-
      h1(lambda.ps,u.ps,alpha,log=T)+
      dtnorm(u.ps,canu)-
      dtnorm(canu,u.ps)
    acc<-matrix(runif(nt*nF),nF,nt)
    u.ps<-ifelse(acc<exp(R),canu,u.ps)
    
    ##########################################################
    ##############              alpha             ############
    ##########################################################
    
    
    canalpha<-rnorm(1,alpha,MHb[4])    
    if(i>50 & canalpha>0 & canalpha<1){
      attb[4]<-attb[4]+1
      cantheta<-make.theta(FAC,logs,canalpha)
      canll<-curll  
      for(t in 1:nt){
        canll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta[,t],canalpha)
      }
      R<-sum(canll-curll)+
        sum(h1(lambda.ps,u.ps,canalpha))-
        sum(h1(lambda.ps,u.ps,alpha))+
        dbeta(canalpha,pri.alpha.a,pri.alpha.b,log=T)-
        dbeta(alpha,pri.alpha.a,pri.alpha.b,log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        alpha<-canalpha;curll<-canll;theta<-cantheta;
        accb[4]<-accb[4]+1
      }}           
    } 
    
    ##########################################################
    ##############       KERNEL BANDWIDTH         ############
    ##########################################################
    
    attb[5]<-attb[5]+1
    canlogrange<-rnorm(1,logrange,MHb[5])
    canfac<-make.kern(dw2,canlogrange)
    canFAC<-stdKern(canfac)
    canll<-curll
    cantheta<-make.theta(canFAC,logs,alpha)
    for(t in 1:nt){
      canll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta[,t],alpha)
    }
    R<-sum(canll-curll)+
      dnorm(canlogrange,pri.mn.bw,pri.sd.bw,log=T)-
      dnorm(logrange,pri.mn.bw,pri.sd.bw,log=T)
    if(!is.na(exp(R))){if(runif(1)<exp(R)){
      logrange<-canlogrange;fac<-canfac
      FAC<-canFAC;theta<-cantheta;curll<-canll
      accb[5]<-accb[5]+1
    }}
    
    ##########################################################
    ##############          GEV PARAMETERS        ############
    ##########################################################
    
    ## SPATIALLY VARYING PARAMETERS
    for(l in 1:3){if(i>50 & vary[l]){
      Xb<-as.matrix(X)%*%mnb[,l]
      
      for(j in 1:n){
        VVV<-taub[l]*Qb[j,j]
        MMM<-taub[l]*Qb[j,j]*Xb[j]-
          taub[l]*sum(Qb[-j,j]*(beta[-j,l]-Xb[-j]))
        
        attb[l]<-attb[l]+1
        canb<-beta[j,]
        canb[l]<-rnorm(1,beta[j,l],MHb[l])
        canll<-loglike(y[j,],canb[1],canb[2],canb[3],theta[j,],alpha)
        
        R<-sum(canll-curll[j,])+
          dnorm(canb[l],MMM/VVV,1/sqrt(VVV),log=T)-
          dnorm(beta[j,l],MMM/VVV,1/sqrt(VVV),log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          beta[j,]<-canb
          curll[j,]<-canll
          accb[l]<-accb[l]+1
        }}
      }
      
    }}
    
    ## SPATIALLY CONSTANT PARAMETERS
    for(l in 1:3){if(i>50 & !vary[l]){
      attb[l]<-attb[l]+1
      canb<-beta
      canb[,l]<-rnorm(1,beta[1,l],MHb[l])         
      canb[,l]<-beta[1,l]+MHb[l]*rt(1,df=5)
      canll<-curll
      for(t in 1:nt){
        canll[,t]<-loglike(y[,t],canb[,1],canb[,2],canb[,3],theta[,t],alpha)
      }
      R<-sum(canll-curll)+
        dnorm(canb[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)-
        dnorm(beta[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        beta<-canb
        curll<-canll
        accb[l]<-accb[l]+1
      }}
    }}
    
    
    ##########################################################
    ##############    Spatial hyperparameters     ############
    ##########################################################
    
    if(sum(vary)>0){
      tXQ<-t(X)%*%Qb
      tXQX<-tXQ%*%as.matrix(X)
    }
    
    for(l in 1:3){if(vary[l]){
      #MEAN 
      VVV<-solve(taub[l]*tXQX+(1/pri.sd.beta^2)*diag(p))
      MMM<-taub[l]*tXQ%*%beta[,l]
      mnb[,l]<-VVV%*%MMM+t(chol(VVV))%*%rnorm(p)
      
      #VARIANCE
      SS<-quad.form(Qb,beta[,l]-X%*%mnb[,l])
      taub[l]<-rgamma(1,n/2+pri.var.a,SS/2+pri.var.b)
    }}
    
    #SPATIAL RANGE
    if(sum(vary)>0){
      attb[6]<-attb[6]+1
      canlograngeb<-rnorm(1,lograngeb,MHb[6])
      canQb<-solve(exp(-d/exp(canlograngeb)))
      R<-0.5*sum(vary)*logdet(canQb)-
        0.5*sum(vary)*logdet(Qb)+
        dnorm(canlograngeb,pri.mn.range,pri.sd.range,log=T)-
        dnorm(lograngeb,pri.mn.range,pri.sd.range,log=T)
      for(l in 1:3){if(vary[l]){
        R<-R-0.5*taub[l]*quad.form(canQb,beta[,l]-X%*%mnb[,l])
        R<-R+0.5*taub[l]*quad.form(Qb,beta[,l]-X%*%mnb[,l])
      }}
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        lograngeb<-canlograngeb;Qb<-canQb
        accb[6]<-accb[6]+1
      }}
    }
    
    #########  TUNE THE CANDIDATE DISTRIBUTION  #######
    for(j in 1:length(accb)){if(i<burn/2 & attb[j]>50){
      if(accb[j]/attb[j]<0.3){MHb[j]<-MHb[j]*0.9}
      if(accb[j]/attb[j]>0.6){MHb[j]<-MHb[j]*1.1}
      accb[j]<-attb[j]<-0
    }}
  }#END THIN
    
    
    
    
    #MAKE PREDICTIONS AT NEW LOCATIONS
    if(!is.null(Sp)){
      YYY<-matrix(0,np,nt)
      facp<-make.kern(dw2p,logrange)
      FACp<-stdKern(facp)
      thetap<-make.theta(FACp,logs,alpha)^alpha
      
      bp<-matrix(beta[1,],np,3,byrow=T)
      for(j in 1:3){if(vary[j]){
        RRR<-beta[,j]-X%*%mnb[,j]
        bp[,j]<-Xp%*%mnb[,j]+
          proj.beta(RRR,d12,d22,Qb,taub[j],lograngeb)
      }}
      
      for(t in 1:nt){
        res<-rGEV(np,thetap[,t],alpha*thetap[,t],alpha)
        YYY[,t]<-bp[,1]+exp(bp[,2])*(res^(bp[,3])-1)/bp[,3]
      }
      
      if(i>burn){
        Yp1<-Yp1+YYY/(iters-burn)
        Yp2<-Yp2+YYY*YYY/(iters-burn)
        beta.mn.p<-beta.mn.p+bp/(iters-burn)
        beta.var.p<-beta.mn.p+bp*bp/(iters-burn)
      }
      if(keep.samples){
        Yp[i,,]<-YYY
        locsp[i,,]<-bp
      }
    }

    #KEEP TRACK OF STUFF:
    keepers.bp[i,,]<-bp
    keepers.thetap[i,,]<-thetap
    keepers.FACp[i,,]<-FACp
    keepers.others[i,]<-c(alpha,sum(llp),q1)
    if(i>burn){
      nnn<-iters-burn
      beta.mn<-beta.mn+beta/nnn
      beta.var<-beta.var+beta*beta/nnn
    }
    if(keep.samples){locs[i,,]<-beta}
  }

  
  stop.time <- proc.time()
  
  
  #Return output:
  list(time=stop.time-start.time,
       Yp=Yp,keepers.bp=keepers.bp,
       keepers.thetap=keepers.thetap,
       keepers.FACp=keepers.FACp,
       keepers.others=keepers.others)
}

## BSE_DP
BSE_DP <- function(M=50,y,S,knots=NULL,
                   X=NULL,Sp=NULL,Xp=NULL,
                   init.beta=c(0,0,.0001),
                   init.alpha=.4,
                   init.range=1,
                   init.bw=1,
                   vary=c(T,F,F),
                   pri.mn.range=-2,
                   pri.sd.range=1,
                   pri.mn.bw=0,
                   pri.sd.bw=1,
                   pri.var.a=0.1,
                   pri.var.b=0.1,
                   pri.alpha.a=1,
                   pri.alpha.b=1,
                   pri.sd.beta=10,
                   pri.mn.gev=c(0,0,0),
                   pri.sd.gev=c(10,1,0.25),
                   keep.samples=T,
                   iters=iters,burn=burn,
                   update=50,nthin=1){
  
  
  start.time <- proc.time()
  library(fields)
  library(emulator)
  
  if(is.null(knots)){knots<-S}
  if(is.null(X)){X<-cbind(1,S)}
  
  
  #BOOKEEPING
  X<-as.matrix(X)
  n<-nrow(y)
  np<-nrow(Sp)
  nt<-ncol(y)
  nF<-nrow(knots)
  p<-ncol(X)
  # v1<-c(0,0)
  v2<-rep(1,M)
  b1<-1
  
  d<-rdist(S,S)
  diag(d)<-0
  dw2<-rdist(S,knots)^2
  dw2[dw2<0.0001]<-0
  
  #INITIAL VALUES:
  
  beta<-matrix(init.beta,n,3,byrow=T)
  if(vary[1]){for(j in 1:n){
    beta[j,1]<-get.inits.mu(y[j,],beta[j,2],beta[j,3])
  }}
  mnb<-matrix(0,p,3)
  mnb[1,]<-colMeans(beta)
  taub<-rep(1,3)
  lograngeb<-log(init.range)
  logrange<-log(init.bw)
  lambda.dp <- matrix(2, ncol=M, nrow=nF)
  u.dp<-matrix(0.5, ncol=M, nrow=nF)
  u<-matrix(0.5,nF,nt)
  alpha<-init.alpha
  g<-as.vector(matrix(sample(1:M,nt,replace=T),ncol=nt, nrow=1))
  q1<-0
  q2<-rep(0.5,M)
  q2<-q2/sum(q2)
  logs<- exp(lambda.dp[,g])
  K.names <- as.numeric(names(table(g)))
  
  
  #COMPUTE QUANTITIES USED IN THE POSTERIOR
  Qb<-solve(exp(-d/exp(lograngeb)))
  fac<-make.kern(dw2,logrange)
  FAC<-stdKern(fac)
  curll<-matrix(0,n,nt)
  theta<-make.theta(FAC,logs,alpha)
  for(j in 1:n){
    curll[j,]<-loglike(y[j,],beta[j,1],beta[j,2],beta[j,3],theta[j,],alpha)
  }
  
  
  #CREATE PLACES TO STORE THE OUTPUT
  keepers.bp<-array(0,c(iters,np,3))     # PREDICTED GEV PARAMETERS
  keepers.thetap<-array(0,c(iters,np,nt))# PREDICTED THETA
  keepers.FACp<-array(0,c(iters,np,nF))  # PREDICTED FAC
  keepers.lambda<-array(0,c(iters,nF,M)) # FIXED EFFECT 
  keepers.q2<-matrix(0,iters,M)          # MIXTURE COMPONENTS PROPORTIONS 
  keepers.others<-matrix(0,iters,3)      # ALPHA,Q1 AND LOGLIKELIHOOD
  
  colnames(keepers.others)<-c("alpha","q","log likelihood")
  GEVs<-c("GEV location","GEV log scale","GEV shape")
  SITEs<-paste("Site",1:n)
  beta.mn<-beta.var<-0*beta
  colnames(beta.mn)<-colnames(beta.var)<-GEVs
  rownames(beta.mn)<-rownames(beta.var)<-SITEs
  locs=NULL
  if(keep.samples){
    locs<-array(0,c(iters,n,3))
    dimnames(locs)<-list(paste("sample",1:iters),SITEs,GEVs)
  }
  
  Yp<-Yp1<-Yp2<-locsp<-beta.mn.p<-beta.var.p<-NULL
  if(!is.null(Sp)){
    SITEs<-paste("Pred site",1:np)
    TIMEs<-paste("Replication",1:nt)
    Yp1<-Yp2<-matrix(0,np,nt)
    llp<-matrix(0,np,nt)
    colnames(Yp1)<-colnames(Yp2)<-TIMEs
    rownames(Yp1)<-rownames(Yp2)<-SITEs
    beta.mn.p<-beta.var.p<-matrix(0,np,3)   
    colnames(beta.mn.p)<-colnames(beta.var.p)<-GEVs
    rownames(beta.mn.p)<-rownames(beta.var.p)<-SITEs
    if(keep.samples){
      locsp<-array(0,c(iters,np,3))
      dimnames(locsp)<-list(paste("sample",1:iters),SITEs,GEVs)
      Yp<-array(0,c(iters,np,nt))
      dimnames(Yp)<-list(paste("sample",1:iters),SITEs,TIMEs)
    }
    dw2p<-rdist(Sp,knots)^2
    dw2p[dw2p<0.0001]<-0
    if(sum(vary)>0){
      d12<-rdist(Sp,S)
      d22<-rdist(Sp,Sp)
      d12[d12<0.0001]<-0
      diag(d22)<-0
    }
  }
  
  #SETTINGS FOR THE M-H CANDIDATE DISTRIBUTION
  attb<-accb<-MHb<-c(0.1,0.02,0.02,0.01,0.02,0.02,1)
  MHu<-3
  MHq<-0.25
  attq<-accq<-1
  atts<-accs<-MHs<-c(3,2,2,2,rep(1,10))
  attp<-accp<-MHp<-c(3,2,2,2,rep(1,10))
  cuts<-seq(0,15,2)
  
  
  #START SAMPLING!
  # level.ps<-olds.ps<-lambda.ps
  level.dp<-olds.dp<-lambda.dp
  
  for(i in 1:iters){for(rep in 1:nthin){
    olds.dp<-lambda.dp
    
    ##########################################################
    ##############  CLUSTER LABEL g ##########################
    ##########################################################
    for(t in 1:nt){
      if(g[t] %in% g[-t]){
        cang<-sample((1:M)[-K.names],1,replace=T)
        canlambda<-rep(0,nF)
        for(l in 1:nF){
          canlambda[l]<-log(rstable.posit(alpha))
        }
        canlogs<-exp(canlambda)
        cantheta<-FAC^(1/alpha)%*%(canlogs)
        canll<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta,alpha)
        R<-sum(canll,na.rm=T)-sum(curll[,t],na.rm=T)+log(b1)-log(nt-1)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          g[t]<-cang
          curll[,t]<-canll
          lambda.dp[,cang]<-canlambda  # A_l
          theta[,t]<-cantheta # theta
          logs[,t]<-canlogs
        }}
      }else{
        cang<-sample(g[-t],1,replace=T)
        canlambda<-lambda.dp[,cang]
        canlogs<-exp(canlambda)
        cantheta<-FAC^(1/alpha)%*%(canlogs)
        canll<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta,alpha)
        R<-sum(canll,na.rm=T)-sum(curll[,t],na.rm=T)-log(b1)+log(nt-1)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          g[t]<-cang
          curll[,t]<-canll
          lambda.dp[,cang]<-canlambda  # A_l
          theta[,t]<-cantheta # theta
          logs[,t]<-canlogs
        }}
      }
    }
    
    for(t in 1:nt){
      q2<-rep(0,nt)
      if(g[t] %in% g[-t]){
        q2[as.numeric(names(table(g[-t])))]<-table(g[-t])/(nt-1)
        temp = exp(lambda.dp)
        cantheta = FAC^(1/alpha)%*%(temp)
        g.prop = q2*exp(colSums(loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta,alpha)))
        g.prop[is.na(g.prop)]<-0
        if(sum(g.prop)==0){
          g.prop=rep(1, M)
        }else{
          g.prop=g.prop/sum(g.prop)
        }
        g[t] = sample(1:length(g.prop), 1, replace=T, prob=g.prop)
      }
    }
    
    ##########################################################
    ##############  RANDOM EFFECT S(DP) ######################
    ##########################################################
    ### Update dp
    # lambda: MxL
    # v: nxt
    for(j in 1:M){
      for(l in 1:nF){
        level.dp[l,j]<-get.level(lambda.dp[l,j],cuts)
      }
    }
    
    K.names <- as.numeric(names(table(g)))
    
    for(j in 1:M){
      if(j %in% g){
        for(l in 1:nF){
          W<-FAC[,l]^(1/alpha)
          MH1<-MHs[level.dp[l,j]]
          canlambda<-rnorm(1,lambda.dp[l,j],MH1)
          MH2<-MHs[get.level(canlambda,cuts)]
          
          canlogs<-exp(canlambda)
          cantheta<-theta[,g==j]+W*(canlogs-exp(lambda.dp[l,j]))
          canll<-loglike(y[,g==j],beta[,1],beta[,2],beta[,3],
                         cantheta,alpha)
          R<-sum(canll,na.rm=T)-sum(curll[,g==j],na.rm=T)+
            #h1(canlambda,u.dp[l,j],alpha,log=T)-
            #h1(lambda.dp[l,j],u.dp[l,j],alpha,log=T)+
            #dnorm(lambda.dp[l,j],canlambda,MH2,log=T)-
            #dnorm(canlambda,lambda.dp[l,j],MH1,log=T)
           dnorm(lambda.dp[l,j],canlambda,MH2,log=T)-
           dnorm(canlambda,lambda.dp[l,j],MH1,log=T)+
           dnorm(canlambda,pri.mn.bw,pri.sd.bw,log=T)-
           dnorm(lambda.dp[l,j],pri.mn.bw,pri.sd.bw,log=T)
          if(!is.na(exp(R))){if(runif(1)<exp(R)){
            curll[,g==j]<-canll
            lambda.dp[l,j]<-canlambda  # A_l
            theta[,g==j]<-cantheta # theta
            logs[l,g==j]<-canlogs
          }}
        }
      }else{
        for(l in 1:nF){
          lambda.dp[l,j]<-log(rstable.posit(alpha))
        }
      }
    }
    
    
    for(j in 1:length(MHs)){
      accs[j]<-accs[j]+sum(olds.dp[,K.names][level.dp[,K.names]==j]!=lambda.dp[,K.names][level.dp[,K.names]==j])
      atts[j]<-atts[j]+sum(level.dp[,K.names]==j)
    }
    for(j in 1:length(atts)){if(i<burn/2 & atts[j]>50){
      if(accs[j]/atts[j]<0.3){MHs[j]<-MHs[j]*0.9}
      if(accs[j]/atts[j]>0.6){MHs[j]<-MHs[j]*1.1}
      accs[j]<-atts[j]<-0
    }}    
    
    
    #canu<-rtnorm(u.dp[,K.names])
    #R<-h1(lambda.dp[,K.names],canu,alpha,log=T)-
    #  h1(lambda.dp[,K.names],u.dp[,K.names],alpha,log=T)+
    #  dtnorm(u.dp[,K.names],canu)-
    #  dtnorm(canu,u.dp[,K.names])
    #acc<-matrix(runif(length(K.names)*nF),ncol=length(K.names),nrow=nF)
    #u.dp[,K.names]<-ifelse(acc<exp(R),canu,u.dp[,K.names])
    

    logs<-exp(lambda.dp[,g])
    theta<-FAC^(1/alpha)%*%(logs)
    
    for(t in 1:nt){
      curll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],theta[,t],alpha)
    }
    
    ##########################################################
    ##############              alpha             ############
    ##########################################################
    
    
    canalpha<-rnorm(1,alpha,MHb[4])    
    if(i>50 & canalpha>0 & canalpha<1){
      attb[4]<-attb[4]+1
      cantheta<-make.theta(FAC,logs,canalpha)
      canll<-curll  
      for(t in 1:nt){
        canll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta[,t],canalpha)
      }
      R<-sum(canll-curll)+
        dbeta(canalpha,pri.alpha.a,pri.alpha.b,log=T)-
        dbeta(alpha,pri.alpha.a,pri.alpha.b,log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        alpha<-canalpha;curll<-canll;theta<-cantheta;
        accb[4]<-accb[4]+1
      }}           
    } 
    
    ##########################################################
    ##############       KERNEL BANDWIDTH         ############
    ##########################################################
    
    attb[5]<-attb[5]+1
    canlogrange<-rnorm(1,logrange,MHb[5])
    canfac<-make.kern(dw2,canlogrange)
    canFAC<-stdKern(canfac)
    canll<-curll
    cantheta<-make.theta(canFAC,logs,alpha)
    for(t in 1:nt){
      canll[,t]<-loglike(y[,t],beta[,1],beta[,2],beta[,3],cantheta[,t],alpha)
    }
    R<-sum(canll-curll)+
      dnorm(canlogrange,pri.mn.bw,pri.sd.bw,log=T)-
      dnorm(logrange,pri.mn.bw,pri.sd.bw,log=T)
    if(!is.na(exp(R))){if(runif(1)<exp(R)){
      logrange<-canlogrange;fac<-canfac
      FAC<-canFAC;theta<-cantheta;curll<-canll
      accb[5]<-accb[5]+1
    }}
    
    ##########################################################
    ##############          GEV PARAMETERS        ############
    ##########################################################
    
    ## SPATIALLY VARYING PARAMETERS
    for(l in 1:3){if(i>50 & vary[l]){
      Xb<-as.matrix(X)%*%mnb[,l]
      
      for(j in 1:n){
        VVV<-taub[l]*Qb[j,j]
        MMM<-taub[l]*Qb[j,j]*Xb[j]-
          taub[l]*sum(Qb[-j,j]*(beta[-j,l]-Xb[-j]))
        
        attb[l]<-attb[l]+1
        canb<-beta[j,]
        canb[l]<-rnorm(1,beta[j,l],MHb[l])
        canll<-loglike(y[j,],canb[1],canb[2],canb[3],theta[j,],alpha)
        
        R<-sum(canll-curll[j,])+
          dnorm(canb[l],MMM/VVV,1/sqrt(VVV),log=T)-
          dnorm(beta[j,l],MMM/VVV,1/sqrt(VVV),log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          beta[j,]<-canb
          curll[j,]<-canll
          accb[l]<-accb[l]+1
        }}
      }
      
    }}
    
    ## SPATIALLY CONSTANT PARAMETERS
    for(l in 1:3){if(i>50 & !vary[l]){
      attb[l]<-attb[l]+1
      canb<-beta
      canb[,l]<-rnorm(1,beta[1,l],MHb[l])         
      canb[,l]<-beta[1,l]+MHb[l]*rt(1,df=5)
      canll<-curll
      for(t in 1:nt){
        canll[,t]<-loglike(y[,t],canb[,1],canb[,2],canb[,3],theta[,t],alpha)
      }
      R<-sum(canll-curll)+
        dnorm(canb[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)-
        dnorm(beta[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        beta<-canb
        curll<-canll
        accb[l]<-accb[l]+1
      }}
    }}
    
    
    ##########################################################
    ##############    SPATIAL HYPERPARAMETERS     ############
    ##########################################################
    
    if(sum(vary)>0){
      tXQ<-t(X)%*%Qb
      tXQX<-tXQ%*%as.matrix(X)
    }
    
    for(l in 1:3){if(vary[l]){
      #MEAN 
      VVV<-solve(taub[l]*tXQX+(1/pri.sd.beta^2)*diag(p))
      MMM<-taub[l]*tXQ%*%beta[,l]
      mnb[,l]<-VVV%*%MMM+t(chol(VVV))%*%rnorm(p)
      
      #VARIANCE
      SS<-quad.form(Qb,beta[,l]-X%*%mnb[,l])
      taub[l]<-rgamma(1,n/2+pri.var.a,SS/2+pri.var.b)
    }}
    
    #SPATIAL RANGE
    if(sum(vary)>0){
      attb[6]<-attb[6]+1
      canlograngeb<-rnorm(1,lograngeb,MHb[6])
      canQb<-solve(exp(-d/exp(canlograngeb)))
      R<-0.5*sum(vary)*logdet(canQb)-
        0.5*sum(vary)*logdet(Qb)+
        dnorm(canlograngeb,pri.mn.range,pri.sd.range,log=T)-
        dnorm(lograngeb,pri.mn.range,pri.sd.range,log=T)
      for(l in 1:3){if(vary[l]){
        R<-R-0.5*taub[l]*quad.form(canQb,beta[,l]-X%*%mnb[,l])
        R<-R+0.5*taub[l]*quad.form(Qb,beta[,l]-X%*%mnb[,l])
      }}
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        lograngeb<-canlograngeb;Qb<-canQb
        accb[6]<-accb[6]+1
      }}
    }
    
    #########  TUNE THE CANDIDATE DISTRIBUTION  #######
    for(j in 1:length(accb)){if(i<burn/2 & attb[j]>50){
      if(accb[j]/attb[j]<0.3){MHb[j]<-MHb[j]*0.9}
      if(accb[j]/attb[j]>0.6){MHb[j]<-MHb[j]*1.1}
      accb[j]<-attb[j]<-0
    }}}#END THIN
    
    
    
    #MAKE PREDICTIONS AT NEW LOCATIONS
    if(!is.null(Sp)){
      YYY<-matrix(0,np,nt)
      facp<-make.kern(dw2p,logrange)
      FACp<-stdKern(facp)
      thetap<-make.theta(FACp,logs,alpha)^alpha
      
      bp<-matrix(beta[1,],np,3,byrow=T)
      for(j in 1:3){if(vary[j]){
        RRR<-beta[,j]-X%*%mnb[,j]
        bp[,j]<-Xp%*%mnb[,j]+
          proj.beta(RRR,d12,d22,Qb,taub[j],lograngeb)
      }}
      
      for(t in 1:nt){
        res<-rGEV(np,thetap[,t],alpha*thetap[,t],alpha)
        YYY[,t]<-bp[,1]+exp(bp[,2])*(res^(bp[,3])-1)/bp[,3]
      }

      if(i>burn){
        Yp1<-Yp1+YYY/(iters-burn)
        Yp2<-Yp2+YYY*YYY/(iters-burn)
        beta.mn.p<-beta.mn.p+bp/(iters-burn)
        beta.var.p<-beta.mn.p+bp*bp/(iters-burn)
      }
      if(keep.samples){
        Yp[i,,]<-YYY
        locsp[i,,]<-bp
      }
    }

    #KEEP TRACK OF STUFF:
    keepers.bp[i,,]<-bp
    keepers.thetap[i,,]<-thetap
    keepers.FACp[i,,]<-FACp
    keepers.lambda[i,,]<-lambda.dp
    keepers.others[i,]<-c(alpha,NA,NA)
    if(i>burn){
      nnn<-iters-burn
      beta.mn<-beta.mn+beta/nnn
      beta.var<-beta.var+beta*beta/nnn
    }
    if(keep.samples){locs[i,,]<-beta}
    
    print(i)
  }
  
  stop.time <- proc.time()
  
  #Return output:
  list(time=stop.time-start.time,
       Yp=Yp,keepers.bp=keepers.bp,
       keepers.thetap=keepers.thetap,
       keepers.FACp=keepers.FACp,
       keepers.lambda=keepers.lambda,
       keepers.others=keepers.others)
}



## BSE_MM
BSE_MM <- function(M=50,y,S,knots=NULL,
                   X=NULL,Sp=NULL,Xp=NULL,
                   init.beta=c(0,0,.0001),
                   init.alpha=.4,
                   init.range=1,
                   init.bw=1,
                   vary=c(T,F,F),
                   pri.mn.range=-2,
                   pri.sd.range=1,
                   pri.mn.bw=0,
                   pri.sd.bw=1,
                   pri.var.a=0.1,
                   pri.var.b=0.1,
                   pri.alpha.a=1,
                   pri.alpha.b=1,
                   pri.sd.beta=10,
                   pri.mn.gev=c(0,0,0),
                   pri.sd.gev=c(10,1,0.25),
                   keep.samples=T,
                   iters=iters,burn=burn,
                   update=50,nthin=1){
  
  start.time <- proc.time()
  library(fields)
  library(emulator)
  
  if(is.null(knots)){knots<-S}
  if(is.null(X)){X<-cbind(1,S)}
  
  
  #BOOKEEPING
  X<-as.matrix(X)
  n<-nrow(y)
  np<-nrow(Sp)
  nt<-ncol(y)
  nF<-nrow(knots)
  p<-ncol(X)
  v2<-rep(1,M)
  
  d<-rdist(S,S)
  diag(d)<-0
  dw2<-rdist(S,knots)^2
  dw2[dw2<0.0001]<-0
  
  #INITIAL VALUES:
  
  beta<-matrix(init.beta,n,3,byrow=T)
  if(vary[1]){for(j in 1:n){
    beta[j,1]<-get.inits.mu(y[j,],beta[j,2],beta[j,3])
  }}
  mnb<-matrix(0,p,3)
  mnb[1,]<-colMeans(beta)
  taub<-rep(1,3)
  lograngeb<-log(init.range)
  logrange<-log(init.bw)
  lambda.ps <- matrix(2, ncol=nt, nrow=nF)
  lambda.dp <- matrix(2, ncol=M, nrow=nF)
  u.ps<-matrix(0.5, ncol=nt, nrow=nF)
  u.dp<-matrix(0.5, ncol=M, nrow=nF)
  u<-matrix(0.5,nF,nt)
  alpha<-init.alpha
  g<-as.vector(matrix(sample(1:M,nt,replace=T),ncol=nt, nrow=1))
  q1<-0.5
  q2<-rep(0.5,M)
  q2<-q2/sum(q2)
  logs.ps<-exp(lambda.ps)
  logs.dp<-exp(lambda.dp[,g])
  
  
  #COMPUTE QUANTITIES USED IN THE POSTERIOR
  Qb<-solve(exp(-d/exp(lograngeb)))
  fac<-make.kern(dw2,logrange)
  FAC<-stdKern(fac)
  curll<-matrix(0,n,nt)
  theta.ps<-(FAC^(1/alpha))%*%(logs.ps)
  theta.dp<-(FAC^(1/alpha))%*%(logs.dp)
  for(j in 1:n){
    curll[j,]<-loglikeMM(y[j,],beta[j,1],beta[j,2],beta[j,3],
                       theta.ps[j,],theta.dp[j,],alpha,q1)
  }
  
  
  #CREATE PLACES TO STORE THE OUTPUT
  keepers.bp<-array(0,c(iters,np,3))     # PREDICTED GEV PARAMETERS
  keepers.thetapsp<-array(0,c(iters,np,nt))# PREDICTED THETA
  keepers.thetadpp<-array(0,c(iters,np,nt))# PREDICTED THETA
  keepers.FACp<-array(0,c(iters,np,nF))  # PREDICTED FAC
  keepers.lambda<-array(0,c(iters,nF,M)) # FIXED EFFECT 
  keepers.q2<-matrix(0,iters,M)          # MIXTURE COMPONENTS PROPORTIONS 
  keepers.others<-matrix(0,iters,3)      # ALPHA,Q1 AND LOGLIKELIHOOD
  
  colnames(keepers.others)<-c("alpha","q","log likelihood")
  GEVs<-c("GEV location","GEV log scale","GEV shape")
  SITEs<-paste("Site",1:n)
  beta.mn<-beta.var<-0*beta
  colnames(beta.mn)<-colnames(beta.var)<-GEVs
  rownames(beta.mn)<-rownames(beta.var)<-SITEs
  locs=NULL
  if(keep.samples){
    locs<-array(0,c(iters,n,3))
    dimnames(locs)<-list(paste("sample",1:iters),SITEs,GEVs)
  }
  
  
  Yp<-Yp1<-Yp2<-locsp<-beta.mn.p<-beta.var.p<-NULL
  if(!is.null(Sp)){
    SITEs<-paste("Pred site",1:np)
    TIMEs<-paste("Replication",1:nt)
    Yp1<-Yp2<-matrix(0,np,nt)
    llp<-matrix(0,np,nt)
    colnames(Yp1)<-colnames(Yp2)<-TIMEs
    rownames(Yp1)<-rownames(Yp2)<-SITEs
    beta.mn.p<-beta.var.p<-matrix(0,np,3)   
    colnames(beta.mn.p)<-colnames(beta.var.p)<-GEVs
    rownames(beta.mn.p)<-rownames(beta.var.p)<-SITEs
    if(keep.samples){
      locsp<-array(0,c(iters,np,3))
      dimnames(locsp)<-list(paste("sample",1:iters),SITEs,GEVs)
      Yp<-array(0,c(iters,np,nt))
      dimnames(Yp)<-list(paste("sample",1:iters),SITEs,TIMEs)
    }
    dw2p<-rdist(Sp,knots)^2
    dw2p[dw2p<0.0001]<-0
    if(sum(vary)>0){
      d12<-rdist(Sp,S)
      d22<-rdist(Sp,Sp)
      d12[d12<0.0001]<-0
      diag(d22)<-0
    }
  }
  
  
  #SETTINGS FOR THE M-H CANDIDATE DISTRIBUTION
  attb<-accb<-MHb<-c(0.1,0.02,0.02,0.01,0.02,0.02,1)
  MHu<-3
  MHq<-0.25
  attq<-accq<-1
  atts<-accs<-MHs<-c(3,2,2,2,rep(1,10))
  attp<-accp<-MHp<-c(3,2,2,2,rep(1,10))
  cuts<-seq(0,15,2)
  
  
  #START SAMPLING!
  level.ps<-olds.ps<-lambda.ps
  level.dp<-olds.dp<-lambda.dp

  for(i in 1:iters){for(rep in 1:nthin){
    olds.ps<-lambda.ps
    olds.dp<-lambda.dp
    
    ##########################################################
    ##############  Cluster label g ##########################
    ##########################################################
    
    # update g
    for(t in 1:nt){
      canthetadp = FAC^(1/alpha)%*%exp(lambda.dp)
      g.prop = q2*exp(colSums(loglikeMM(y[,t],beta[,1],beta[,2],beta[,3],
                                      theta.ps[,t],canthetadp,alpha,q1),na.rm=T))
      g.prop[is.na(g.prop)]<-0
      if(sum(g.prop)==0){
        g.prop=rep(1, M)
      }else{
        g.prop=g.prop/sum(g.prop)
      }
      g[t] = sample(1:M, 1, replace=T, prob=g.prop)
    }
    
    # # Update p
    for(j in 1:M){
      v2[j]<-rbeta(1,1+sum(g==j),1+sum(g>j))
    }
    prev_probs<-q2
    q2<-make.probs(v2,M)
    q2[q2<0] <- 0
    if(is.na(sum(q2))){q2<-prev_probs}
    
    logs.dp<-exp(lambda.dp[,g])
    theta.dp<-FAC^(1/alpha)%*%(logs.dp)
    
    for(t in 1:nt){
      curll[,t]<-loglikeMM(y[,t],beta[,1],beta[,2],beta[,3],
                         theta.ps[,t],theta.dp[,t],alpha,q1)
    }
    
    # update q1 (MH)
    canq1<-rtnorm(q1, sd=MHq)
    canll<-curll
    if(i>50){
      attq<-attq+1
      for(t in 1:nt){
        canll[,t]<-loglikeMM(y[,t],beta[,1],beta[,2],beta[,3],
                           theta.ps[,t],theta.dp[,t],alpha,canq1)
      }
      
      R<-sum(canll,na.rm=T)-sum(curll,na.rm=T)+
        dtnorm(canq1,q1,MHq)-
        dtnorm(q1,canq1,MHq)+
        dbeta(canq1,1,2,log=T)-
        dbeta(q1,1,2,log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        q1<-canq1
        curll<-canll
        accq <- accq + 1
      }}}
    
    
    #########  TUNE THE CANDIDATE DISTRIBUTION  #######
    if(i<burn/2 & attq>50){
      if(accq/attq<0.3){MHq<-MHq*0.9}
      if(accq/attq>0.6){MHq<-MHq*1.1}
      accq<-attq<-0
    }
    

    ##########################################################
    ##############      Random effects S and U    ############
    ##########################################################
    for(t in 1:nt){
      for(l in 1:nF){
        level.ps[l,t]<-get.level(lambda.ps[l,t],cuts)
      }
    }
    ## update ps
    for(t in 1:nt){ 
      for(l in 1:nF){
        W<-FAC[,l]^(1/alpha)
        MH1<-MHp[level.ps[l,t]]
        canlambda<-rnorm(1,lambda.ps[l,t],MH1)
        MH2<-MHp[get.level(canlambda,cuts)]
        
        canlogs<-exp(canlambda)
        cantheta<-theta.ps[,t]+W*(canlogs-exp(lambda.ps[l,t]))
        canll<-loglikeMM(y[,t],beta[,1],beta[,2],beta[,3],
                       cantheta,theta.dp[,t],alpha,q1)
        R<-sum(canll)-sum(curll[,t])+
          h1(canlambda,u.ps[l,t],alpha,log=T)-
          h1(lambda.ps[l,t],u.ps[l,t],alpha,log=T)+
          dnorm(lambda.ps[l,t],canlambda,MH2,log=T)-
          dnorm(canlambda,lambda.ps[l,t],MH1,log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          curll[,t]<-canll
          lambda.ps[l,t]<-canlambda  # A_l
          logs.ps[l,t]<-canlogs
          theta.ps[,t]<-cantheta # theta
        }}
      }
    }

    for(j in 1:length(MHp)){
      accp[j]<-accp[j]+sum(olds.ps[level.ps==j]!=lambda.ps[level.ps==j])
      attp[j]<-attp[j]+sum(level.ps==j)
    }
    for(j in 1:length(attp)){if(i<burn/2 & attp[j]>50){
      if(accp[j]/attp[j]<0.3){MHp[j]<-MHp[j]*0.9}
      if(accp[j]/attp[j]>0.6){MHp[j]<-MHp[j]*1.1}
      accp[j]<-attp[j]<-0
    }}
    
    canu<-rtnorm(u.ps)
    R<-h1(lambda.ps,canu,alpha,log=T)-
      h1(lambda.ps,u.ps,alpha,log=T)+
      dtnorm(u.ps,canu)-
      dtnorm(canu,u.ps)
    acc<-matrix(runif(nt*nF),nF,nt)
    u.ps<-ifelse(acc<exp(R),canu,u.ps)
    
    ### Update dp
  # lambda: MxL
  # v: nxt
  for(j in 1:M){
    for(l in 1:nF){
      level.dp[l,j]<-get.level(lambda.dp[l,j],cuts)
    }
  }
  
  # w<-beta[,3]*exp(-beta[,2])*(y-beta[,1])+1
  # w<-w^(1/beta[,3])
  # ccc<-logd(theta,w)
  for(j in 1:M){
    # ccc<-logd(theta.ps[,g==j],theta.dp[,g==j],w[,g==j],q1,alpha)
    for(l in 1:nF){
      W<-FAC[,l]^(1/alpha)
      MH1<-MHs[level.dp[l,j]]
      canlambda<-rnorm(1,lambda.dp[l,j],MH1)
      MH2<-MHs[get.level(canlambda,cuts)]
      
      canlogs<-exp(canlambda)
      cantheta<-theta.dp[,g==j]+W*(canlogs-exp(lambda.dp[l,j]))
      # canccc<-logd(cantheta,w[,g==j])
      # canccc<-logd(theta.ps[,g==j],cantheta,w[,g==j],q1,alpha)
      canll<-loglikeMM(y[,g==j],beta[,1],beta[,2],beta[,3],
                     theta.ps[,g==j],cantheta,alpha,q1)
      R<-sum(canll,na.rm=T)-sum(curll[,g==j],na.rm=T)+
        dnorm(lambda.dp[l,j],canlambda,MH2,log=T)-
        dnorm(canlambda,lambda.dp[l,j],MH1,log=T)+
        dnorm(canlambda,pri.mn.bw,pri.sd.bw,log=T)-
        dnorm(lambda.dp[l,j],pri.mn.bw,pri.sd.bw,log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        curll[,g==j]<-canll
        lambda.dp[l,j]<-canlambda  # A_l
        theta.dp[,g==j]<-cantheta # theta
        logs.dp[l,g==j]<-canlogs
      }}
    }
  }

  
  for(j in 1:length(MHs)){
    accs[j]<-accs[j]+sum(olds.dp[level.dp==j]!=lambda.dp[level.dp==j])
    atts[j]<-atts[j]+sum(level.dp==j)
  }
  for(j in 1:length(atts)){if(i<burn/2 & atts[j]>50){
    if(accs[j]/atts[j]<0.3){MHs[j]<-MHs[j]*0.9}
    if(accs[j]/atts[j]>0.6){MHs[j]<-MHs[j]*1.1}
    accs[j]<-atts[j]<-0
  }}    
    
    
    ##########################################################
    ##############              alpha             ############
    ##########################################################
    
    
    canalpha<-rnorm(1,alpha,MHb[4])    
    if(i>50 & canalpha>0 & canalpha<1){
      attb[4]<-attb[4]+1
      canthetaps<-(FAC^(1/canalpha))%*%(logs.ps)
      canthetadp<-(FAC^(1/canalpha))%*%(logs.dp)
      canll<-curll  
      for(t in 1:nt){
        canll[,t]<-loglikeMM(y[,t],beta[,1],beta[,2],beta[,3],
                           canthetaps[,t],canthetadp[,t],canalpha,q1)
      }
      R<-sum(canll,na.rm=T)-sum(curll,na.rm=T)+
        sum(h1(lambda.dp,u.dp,canalpha))-
        sum(h1(lambda.dp,u.dp,alpha))+
        dbeta(canalpha,pri.alpha.a,pri.alpha.b,log=T)-
        dbeta(alpha,pri.alpha.a,pri.alpha.b,log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        alpha<-canalpha;curll<-canll;
        theta.ps<-canthetaps
        theta.dp<-canthetadp;
        accb[4]<-accb[4]+1
      }}           
    } 
    
    ##########################################################
    ##############       KERNEL BANDWIDTH         ############
    ##########################################################
    
    attb[5]<-attb[5]+1
    canlogrange<-rnorm(1,logrange,MHb[5])
    canfac<-make.kern(dw2,canlogrange)
    canFAC<-stdKern(canfac)
    canll<-curll
    canthetaps<-(canFAC^(1/alpha))%*%(logs.ps)
    canthetadp<-(canFAC^(1/alpha))%*%(logs.dp)
    for(t in 1:nt){
      canll[,t]<-loglikeMM(y[,t],beta[,1],beta[,2],beta[,3],
                         canthetaps[,t],canthetadp[,t],alpha,q1)
    }
    R<-sum(canll,na.rm=T)-sum(curll,na.rm=T)+
      dnorm(canlogrange,pri.mn.bw,pri.sd.bw,log=T)-
      dnorm(logrange,pri.mn.bw,pri.sd.bw,log=T)
    if(!is.na(exp(R))){if(runif(1)<exp(R)){
      logrange<-canlogrange;fac<-canfac
      FAC<-canFAC;
      theta.ps<-canthetaps;
      theta.dp<-canthetadp;
      curll<-canll
      accb[5]<-accb[5]+1
    }}
    
    ##########################################################
    ##############          GEV PARAMETERS        ############
    ##########################################################
    
    ## SPATIALLY VARYING PARAMETERS
    for(l in 1:3){if(i>50 & vary[l]){
      Xb<-as.matrix(X)%*%mnb[,l]
      
      for(j in 1:n){
        VVV<-taub[l]*Qb[j,j]
        MMM<-taub[l]*Qb[j,j]*Xb[j]-
          taub[l]*sum(Qb[-j,j]*(beta[-j,l]-Xb[-j]))
        
        attb[l]<-attb[l]+1
        canb<-beta[j,]
        canb[l]<-rnorm(1,beta[j,l],MHb[l])
        canll<-loglikeMM(y[j,],canb[1],canb[2],canb[3],
                       theta.ps[j,],theta.dp[j,],alpha,q1)
        
        R<-sum(canll,na.rm=T)-sum(curll[j,],na.rm=T)+
          dnorm(canb[l],MMM/VVV,1/sqrt(VVV),log=T)-
          dnorm(beta[j,l],MMM/VVV,1/sqrt(VVV),log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          beta[j,]<-canb
          curll[j,]<-canll
          accb[l]<-accb[l]+1
        }}
      }
      
    }}
    
    ## SPATIALLY CONSTANT PARAMETERS
    for(l in 1:3){if(i>50 & !vary[l]){
      attb[l]<-attb[l]+1
      canb<-beta
      # canb[,l]<-rnorm(1,beta[1,l],MHb[l])         
      canb[,l]<-beta[1,l]+MHb[l]*rt(1,df=5)
      canll<-curll
      for(t in 1:nt){
        canll[,t]<-loglikeMM(y[,t],canb[,1],canb[,2],canb[,3],
                           theta.ps[,t],theta.dp[,t],alpha,q1)
      }
      R<-sum(canll,na.rm=T)-sum(curll,na.rm=T)+
        dnorm(canb[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)-
        dnorm(beta[1,l],pri.mn.gev[l],pri.sd.gev[l],log=T)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        beta<-canb
        curll<-canll
        accb[l]<-accb[l]+1
      }}
    }}
    
    
    ##########################################################
    ##############    Spatial hyperparameters     ############
    ##########################################################
    
    if(sum(vary)>0){
      tXQ<-t(X)%*%Qb
      tXQX<-tXQ%*%as.matrix(X)
    }
    
    for(l in 1:3){if(vary[l]){
      #MEAN 
      VVV<-solve(taub[l]*tXQX+(1/pri.sd.beta^2)*diag(p))
      MMM<-taub[l]*tXQ%*%beta[,l]
      mnb[,l]<-VVV%*%MMM+t(chol(VVV))%*%rnorm(p)
      
      #VARIANCE
      SS<-quad.form(Qb,beta[,l]-X%*%mnb[,l])
      taub[l]<-rgamma(1,n/2+pri.var.a,SS/2+pri.var.b)
    }}
    
    #SPATIAL RANGE
    if(sum(vary)>0){
      attb[6]<-attb[6]+1
      canlograngeb<-rnorm(1,lograngeb,MHb[6])
      canQb<-solve(exp(-d/exp(canlograngeb)))
      R<-0.5*sum(vary)*logdet(canQb)-
        0.5*sum(vary)*logdet(Qb)+
        dnorm(canlograngeb,pri.mn.range,pri.sd.range,log=T)-
        dnorm(lograngeb,pri.mn.range,pri.sd.range,log=T)
      for(l in 1:3){if(vary[l]){
        R<-R-0.5*taub[l]*quad.form(canQb,beta[,l]-X%*%mnb[,l])
        R<-R+0.5*taub[l]*quad.form(Qb,beta[,l]-X%*%mnb[,l])
      }}
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        lograngeb<-canlograngeb;Qb<-canQb
        accb[6]<-accb[6]+1
      }}
    }
    
    #########  TUNE THE CANDIDATE DISTRIBUTION  #######
    for(j in 1:length(accb)){if(i<burn/2 & attb[j]>50){
      if(accb[j]/attb[j]<0.3){MHb[j]<-MHb[j]*0.9}
      if(accb[j]/attb[j]>0.6){MHb[j]<-MHb[j]*1.1}
      accb[j]<-attb[j]<-0
    }}
    
  }#end nthin
    
    #MAKE PREDICTIONS AT NEW LOCATIONS
    if(!is.null(Sp)){
      YYY<-matrix(0,np,nt)
      facp<-make.kern(dw2p,logrange)
      FACp<-stdKern(facp)
      thetapsp<-(FACp^(1/alpha))%*%(logs.ps)
      thetadpp<-(FACp^(1/alpha))%*%(logs.dp)
      
      bp<-matrix(beta[1,],np,3,byrow=T)
      for(j in 1:3){if(vary[j]){
        RRR<-beta[,j]-X%*%mnb[,j]
        bp[,j]<-Xp%*%mnb[,j]+
          proj.beta(RRR,d12,d22,Qb,taub[j],lograngeb)
      }}
      
      for(t in 1:nt){
        YYY[,t]<-rmmGEV(np,bp[,1],exp(bp[1,2]),bp[1,3],
                        thetapsp[,t],thetadpp[,t],alpha,q1)
      }
      
      if(i>burn){
        Yp1<-Yp1+YYY/(iters-burn)
        Yp2<-Yp2+YYY*YYY/(iters-burn)
        beta.mn.p<-beta.mn.p+bp/(iters-burn)
        beta.var.p<-beta.mn.p+bp*bp/(iters-burn)
      }
      if(keep.samples){
        Yp[i,,]<-YYY
        locsp[i,,]<-bp
      }
    }
    
    #KEEP TRACK OF STUFF:
    keepers.bp[i,,]<-bp
    keepers.thetapsp[i,,]<-thetapsp
    keepers.thetadpp[i,,]<-thetadpp
    keepers.FACp[i,,]<-FACp
    keepers.lambda[i,,]<-lambda.dp
    keepers.q2[i,]<-q2
    keepers.others[i,]<-c(alpha,sum(llp,na.rm=T),q1)
    if(i>burn){
      nnn<-iters-burn
      beta.mn<-beta.mn+beta/nnn
      beta.var<-beta.var+beta*beta/nnn
    }
    if(keep.samples){locs[i,,]<-beta}
  }
  
  stop.time <- proc.time()
  
  #Return output:
  list(time=stop.time-start.time,
       Yp=Yp,keepers.bp=keepers.bp,
       keepers.thetapsp=keepers.thetapsp,
       keepers.thetadpp=keepers.thetadpp,
       keepers.FACp=keepers.FACp,
       keepers.lambda=keepers.lambda,
       keepers.q2=keepers.q2,
       keepers.others=keepers.others)
}