library(fields)
library(emulator)
library(BMAmevt)


################################################################
# From Morris(2017)
# Arguments:
#   preds(iters, yp, nt): mcmc predictions at validation
#                         locations
#   probs(nprobs): sample quantiles for scoring
#   validate(np, nt): validation data
#   trans(bool): are the mcmc predictions transposed
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
QuantScore <- function(preds, probs, validate, trans = FALSE) {
  
  nt <- ncol(validate)  # number of prediction days
  np <- nrow(validate)  # number of prediction sites
  nprobs <- length(probs)  # number of quantiles to find quantile score
  
  # we get the predicted quantile for each site nprobs x np
  if (trans) {
    pred.quants <- apply(preds, 3, quantile, probs = probs, na.rm = T)
  } else {
    pred.quants <- apply(preds, 2, quantile, probs=probs, na.rm=T)
  }
  
  pred.quants <- t(pred.quants)  # need np x nprobs for proper matrix subtraction
  scores.sites <- array(NA, dim=c(nprobs, np, nt))
  
  # we need to figure out how many times the site did or didn't exceed the prediction
  
  for (q in 1:nprobs) {
    diff <- pred.quants[, q] - validate
    i <- diff >= 0  # diff >= 0 means qhat is larger
    scores.sites[q, , ] <- 2 * (i - probs[q]) * diff
  }
  
  scores <- apply(scores.sites, 1, mean, na.rm=T)
  
  return(scores)
}

#############################################################
########   FUNCTIONS TO COMPUTE INITIAL VALUES    ###########
#############################################################

get.inits.mu<-function(y,logsig=0,xi=.1){
  m<-median(y,na.rm=T)
  mu<-m-exp(logsig)*(log(2)^(-xi)-1)/xi
  return(mu)}


######################################################
########   THE POSITIVE STABLE DENSITY    ###########
######################################################


h1<-function(logs,u,alpha,log=T){
  s<-exp(logs)
  psi<-pi*u
  c<-(sin(alpha*psi)/sin(psi))^(1/(1-alpha))
  c<-c*sin((1-alpha)*psi)/sin(alpha*psi)
  
  logd<-log(alpha)-log(1-alpha)-(1/(1-alpha))*logs+
    log(c)-c*(1/s^(alpha/(1-alpha)))+
    logs
  if(!log){logd<-exp(logd)}
  logd}


######################################################
##########    FUNCTIONS USED FOR PREDICTION  #########
######################################################

rGEV<-function(n,mu,sig,xi){
  tau<-runif(n)
  x<--1/log(tau)
  x<-x^(xi)-1
  x<-mu+sig*x/xi
  x}

rmmGEV<-function(n,mu=1,sig=1,xi=1,theta1,theta2,alpha=.5,q){
  library(evd)
  library(BMAmevt)
  
  mu1.star<-mu+sig*(theta1^(xi*q)*q^xi-1)/xi
  sig1.star<-alpha*q*sig*theta1^(xi*q)*q^xi
  xi1.star<-alpha*xi*q
  y1<-rGEV(n,mu1.star,sig1.star,xi1.star)
  
  mu2.star<-mu+sig*(theta2^(xi*(1-q))*(1-q)^xi-1)/xi
  sig2.star<-alpha*(1-q)*sig*theta2^(xi*(1-q))*(1-q)^xi
  xi2.star<-alpha*xi*(1-q)
  y2<-rGEV(n,mu2.star,sig2.star,xi2.star)
  
  y<-pmax(y1,y2,na.rm = T)
  
  y
}

proj.beta<-function(B,d12,d22,S11inv,tau,logrho){
  #B<-n-vector of observed beta (minus the mean)
  ns<-nrow(d22)
  rho<-exp(logrho)
  
  S22<-exp(-d22/rho)/tau
  S12<-exp(-d12/rho)/tau
  S11inv<-S11inv*tau
  
  P2<-S12%*%S11inv
  P1<-S22-S12%*%S11inv%*%t(S12)
  P1<-t(chol(P1))    
  
  Bnew<-P2%*%B+P1%*%rnorm(ns)
  return(Bnew)}



######################################################
####  FUNCTION TO COMPUTE THE RANDOM EFFECTS  ########
######################################################

make.theta<-function(FAC,logs,alpha){
  #FAC is nxnF
  #logs is nxnFxnt
  #alpha in (0,1)
  n <- dim(FAC)[1]
  nt <- dim(logs)[2]
  xxx <- matrix(0,n,nt)
  if(length(logs)==1){
    xxx<-(FAC^(1/alpha))*(logs)}
  if(length(logs)>1){
    xxx<-(FAC^(1/alpha))%*%(logs)}
  xxx}  


stdKern<-function(w,single=F){
  if(single){K<-w/sum(w)}   
  if(!single){K<-sweep(w,1,rowSums(w),"/")}
  K}  

make.kern<-function(d2,logrho){
  rho2<-exp(logrho)^2
  w<-exp(-0.5*d2/rho2)
  w}

#Use normal instead of t
#Compute the stick-breaking probabilities
#v: current probability #m: number of clusters
make.probs<-function(v,m){
  if(m==1){
    p=1
  }
  if(m>1){
    p<-v
    p[2:m]<-p[2:m]*cumprod(1-v[2:m-1])
    p[m]<-1-sum(p[2:m-1])
  }
  p
}

get.logs<-function(lambda.ps,lambda.dp,g,nt,nF){
  logs<-matrix(2,nF,nt)
  for(t in 1:nt){
    if(g[t]==0){
      logs[,t] <- lambda.ps[,t]
    }else{
      logs[,t] <- lambda.dp[g[t],]
    }
  }
  logs
}

######################################################
########            OTHER FUNCTIONS        ###########
######################################################

get.level<-function(logs,cuts){
  sum(logs>cuts)+1
}

logdet<-function(X){
  determinant(X)$modulus
}

trunc<-function(x,eps=.1){
  x<-ifelse(x<eps,eps,x)
  x<-ifelse(x>1-eps,1-eps,x)
  x}

rtnorm<-function(mn,sd=.25,fudge=0){
  upper<-pnorm(1-fudge,mn,sd)
  lower<-pnorm(fudge,mn,sd)
  if(is.matrix(mn)){
    U<-matrix(runif(prod(dim(mn)),lower,upper),dim(mn)[1],dim(mn)[2])
  }
  if(!is.matrix(mn)){
    U<-runif(length(mn),lower,upper)
  }
  return(qnorm(U,mn,sd))}

dtnorm<-function(y,mn,sd=.25,fudge=0){
  upper<-pnorm(1-fudge,mn,sd)
  lower<-pnorm(fudge,mn,sd)
  l<-dnorm(y,mn,sd,log=T)-log(upper-lower)
  return(l)}

loglike<-function(y,mu,logsig,xi,theta,alpha){
  missing<-is.na(y)
  theta<-theta^alpha
  mu.star<-mu+exp(logsig)*(theta^xi-1)/xi
  sig.star<-alpha*exp(logsig)*(theta^xi)
  xi.star<-alpha*xi
  ttt<-(1+xi.star*(y-mu.star)/sig.star)^(-1/xi.star)
  lll<- -log(sig.star)+(xi.star+1)*log(ttt)-ttt
  lll[missing]<-0
  lll}

loglikeMM<-function(y,mu,logsig,xi,theta1,theta2,alpha,q){
  missing<-is.na(y)
  theta1<-theta1^alpha
  theta2<-theta2^alpha
  mu1.star<-mu+exp(logsig)*(theta1^(xi*q)*q^xi-1)/xi
  sig1.star<-alpha*q*exp(logsig)*theta1^(xi*q)*q^xi
  xi1.star<-alpha*xi*q
  mu2.star<-mu+exp(logsig)*(theta2^(xi*(1-q))*(1-q)^xi-1)/xi
  sig2.star<-alpha*(1-q)*exp(logsig)*theta2^(xi*(1-q))*(1-q)^xi
  xi2.star<-alpha*xi*(1-q)
  
  ttt1<-(1+xi1.star*(y-mu1.star)/sig1.star)^(-1/xi1.star)
  ttt2<-(1+xi2.star*(y-mu2.star)/sig2.star)^(-1/xi2.star)
  ttt1[is.na(ttt1)]<-0
  ttt2[is.na(ttt2)]<-0
  lll1<- -log(sig1.star)+(xi1.star+1)*log(ttt1)-ttt1-ttt2
  lll2<- -log(sig2.star)+(xi2.star+1)*log(ttt2)-ttt2-ttt1
  lll1[missing]<-0
  lll2[missing]<-0
  lll1[is.na(lll1)]<-0
  lll2[is.na(lll2)]<-0
  lll<-log(exp(lll1)+exp(lll2))
  lll[which(lll==-Inf)]<- -10^4
  #lll[which(lll==-Inf)]<- NA
  lll}


logd<-function(theta,v){
  sum(log(theta)-theta*v) 
}


concdf<-function(y0,y,q2,lambda.dp,FAC,alpha,l){
  if(y0<0){lll<-0}
  else{
    ttt<-y0
    xxx<-(FAC[l,]/ttt)^(1/alpha)%*%exp(lambda.dp)
    pow<- -(xxx)
    lll<-sum(q2*exp(pow))
    if(y0<min(y)){lll[is.na(lll)]<-0}
    if(y0>max(y)){lll[is.na(lll)]<-1}
  }
  lll}

concdfinv<-function(q,y,q2,lambda.dp,FAC,alpha,l){
  uniroot(function(v){concdf(v,y,q2,lambda.dp,FAC,alpha,l)-q}, c(0,2000))$root
}

Bconcdf<-function(y01,y02,y,q2,lambda.dp,FAC,alpha,l1,l2){
  ttt1<-y01
  ttt2<-y02
  xxx<-((FAC[l1,]/ttt1)^(1/alpha)+(FAC[l2,]/ttt2)^(1/alpha))%*%exp(lambda.dp)
  pow<- -(xxx)
  lll<-sum(q2*exp(pow))
  if(y01<min(y) & (y02<min(y))){lll[is.na(lll)]<-0}
  if(y01>max(y) & (y02>max(y))){lll[is.na(lll)]<-1}
  lll}

rmmspatial<-function(nreps,S,knots,mu=1,sig=1,xi=1,alpha=.5,bw=1,M=3,lambda=lambda,q){
  library(evd)
  library(BMAmevt)
  
  n      <- nrow(S)
  nknots <- nrow(knots)
  
  d           <- rdist(S,knots)
  d[d<0.0001] <- 0
  w           <- make.kern(d^2,log(bw))
  K           <- stdKern(w)^(1/alpha)
  
  y1 <- matrix(0,n,nreps)
  for(t in 1:nreps){
    A     <- rep(0,nknots) 
    for(j in 1:nknots){A[j] <- rstable.posit(alpha)}
    theta <- (K%*%A)^alpha
    
    mu1.star<-mu+sig*(theta^(xi*q)-1)/xi
    sig1.star<-alpha*q*sig*(theta^(xi*q))*q^xi
    xi1.star<-alpha*xi*q
    
    y1[,t]    <- rGEV(n,mu1.star,sig1.star,xi1.star)
  }  
  
  y2 <- matrix(0,n,nreps)
  g<-matrix(apply(rmultinom(nreps, 1, c(0.5,0.3,0.2)),2,function(x){which(x==1)}),nreps)
  
  A<-matrix(0,ncol = nreps, nrow=n)
  for(i in 1:nreps){
    A[,i] <- lambda[g[i],]
  }
  
  for(t in 1:nreps){
    theta <- (K%*%A[,t])^alpha
    
    mu2.star<-mu+sig*(theta^(xi*(1-q))-1)/xi
    sig2.star<-alpha*(1-q)*sig*(theta^(xi*(1-q)))*(1-q)^xi
    xi2.star<-alpha*xi*(1-q)
    
    y2[,t]    <- rGEV(n,mu2.star,sig2.star,xi2.star)
  }
  
  y=pmax(y1,y2,na.rm = T)
  
  y
}

rgevspatial<-function(nreps,S,knots,mu=1,sig=1,xi=1,alpha=.5,bw=1){
  library(evd)
  library(BMAmevt)
  
  n      <- nrow(S)
  nknots <- nrow(knots)
  
  d           <- rdist(S,knots)
  d[d<0.0001] <- 0
  w           <- make.kern(d^2,log(bw))
  K           <- stdKern(w)^(1/alpha)
  
  y <- matrix(0,n,nreps)
  for(t in 1:nreps){
    A     <- rep(0,nknots) 
    for(j in 1:nknots){A[j] <- rstable.posit(alpha)}
    theta <- (K%*%A)^alpha
    
    xi_star  <- alpha*xi
    mu_star  <- mu+sig*(theta^xi-1)/xi
    sig_star <- alpha*sig*theta^xi
    
    y[,t]    <- rgev(n,mu_star,sig_star,xi_star)
  }  
  
  return(y)
}


rgaussian<-function(nreps,S,mu=1,bw=1){
  library(evd)
  library(BMAmevt)
  
  n           <- nrow(S)
  d           <- rdist(S,S)
  d[d<0.0001] <- 0
  y           <- matrix(0,n,nreps)
  C           <- exp(-d/bw)
  
  for(t in 1:nreps){
    y[,t]    <- mu+t(chol(C))%*%rnorm(n)
  }  
  
  return(list(y=y,C=C))
}

rgevspatial_dp<-function(nreps,S,knots,mu=1,sig=1,xi=1,alpha=.5,bw=1,M=5,lambda=lambda){
  library(evd)
  library(BMAmevt)
  
  n      <- nrow(S)
  nknots <- nrow(knots)
  
  d           <- rdist(S,knots)
  d[d<0.0001] <- 0
  w           <- make.kern(d^2,log(bw))
  K           <- stdKern(w)^(1/alpha)
  
  y <- matrix(0,n,nreps)
  
  g<-matrix(apply(rmultinom(nreps, 1, c(0.5,0.3,0.2)),2,function(x){which(x==1)}),nreps)
  
  A<-matrix(0,ncol = nreps, nrow=n)
  for(i in 1:nreps){
    A[,i] <- lambda[g[i],]
  }
  
  for(t in 1:nreps){
    theta <- (K%*%A[,t])^alpha
    
    xi_star  <- alpha*xi
    mu_star  <- mu+sig*(theta^xi-1)/xi
    sig_star <- alpha*sig*theta^xi
    
    y[,t]    <- rgev(n,mu_star,sig_star,xi_star)
  }  
  
  y}


rskewt<-function(nreps,S,lambda=3,a=10,b=10,mu=1,bw=1){
  library(evd)
  library(BMAmevt)
  n           <- nrow(S)
  d           <- rdist(S,S)
  d[d<0.0001] <- 0
  y           <- matrix(0,n,nreps)
  C           <- exp(-d/bw)
  for(t in 1:nreps){
    sigma    <- 1/rgamma(1,a/2,b/2)
    z        <- rnorm(1,sd=sigma)
    y[,t]    <- mu+lambda*abs(z)+t(chol(C))%*%rnorm(n)*sigma
  }
  return(list(y=y,C=C))
}


## data generation
data_generate <- function(sim="MS",lambda=NULL){
  #SPECIFY HOW TO GENERATE THE DATA
  S<-expand.grid(1:7,1:7)           #DATA LOCATIONS
  knots<-S                          #KNOT LOCATIONS
  nt<-50                           #NUMBER OF YEARS
  alpha<-.3                         #NUGGET EFFECT
  mu<-rep(S[1,1]/10,49)                      #GEV LOCATION
  sig<-1                            #GEV SCALE
  xi<-.1                            #GEV SHAPE
  bw<-1                             #KERNEL BANDWIDTH
  
  #ASSIGN COVARIATES IN THE PRIOR MEAN OF THE 
  #SPATIALLY-VARYING GEV PARAMETERS
  n<-nrow(S)
  X<-as.matrix(cbind(1,S),n,3)
  
  #WITHHOLD A TEST SET FOR PREDICTION
  test<-runif(n)<.30
  
  #GENERATE DATA FROM THE MODEL
  if(sim=="MS"){
    y<-rgevspatial(nreps=nt,S=S,knots=knots,mu=mu,sig=sig,xi=xi,alpha=alpha,bw=bw)
    Sp<-S[test,]
    res<-list(test=test,X=X,S=S,y=y,knots=knots,alpha=alpha)
  }else if(sim=="DP"){
    y<-rgevspatial_dp(nreps=nt,S=S,knots=knots,mu=mu,sig=sig,xi=xi,alpha=alpha,bw=bw,M=3,lambda=lambda)    
    res<-list(test=test,X=X,S=S,y=y,lambda=lambda,knots=knots,alpha=alpha)
  }else if(sim=="Gaussian"){
    tempy<-rgaussian(nreps=nt,S=S,mu=mu,bw=bw)
    y<-tempy$y
    C<-tempy$C
    res<-list(test=test,X=X,S=S,y=y,C=C,mu=mu)
  }else if(sim=="skewt"){
    tempy<-rskewt(nreps=nt,S=S,lambda=3,a=8,b=2,mu=mu,bw=bw)
    y<-tempy$y
    C<-tempy$C
    res<-list(test=test,X=X,S=S,y=y,C=C,mu=mu)
  }
  
  return(res)
}






# chi-true
chi_true <- function(y,S,Sp,knots,mu=NULL,sig=1,xi=0.1,C=NULL,lambda=NULL,alpha,sim="Sim1"){
  library(fields)
  library(sn)
  library(mvtnorm)
  np<-nrow(Sp)
  if(sim=="Sim1"){
    dw2p<-rdist(Sp,knots)^2
    dw2p[dw2p<0.0001]<-0
    d12<-rdist(Sp,S)
    d22<-rdist(Sp,Sp)
    d12[d12<0.0001]<-0
    diag(d22)<-0
    
    facp<-make.kern(dw2p,1)
    FACp<-stdKern(facp)
    
    chitrue <- matrix(0,nrow=5,ncol=np*(np-1)/2)
    Qtrue <- matrix(0,nrow=5,ncol=np)
    level <- c(0.5,0.9,0.95,0.99,0.995)
    for(i in 1:5){
      temp <- NULL
      for(j in 1:(np-1)){
        for(k in (j+1):np){
          temp <- c(temp, sum(colSums(FACp[c(j,k),]^(1/alpha))^(alpha)))
        }
      }
      chitrue[i,] <- (1-2*level[i]+level[i]^temp)/(1-level[i])
    }
    for(i in 1:5){
      temp <- NULL
      for(j in 1:np){
          temp <- c(temp, qgev(level[i],mu[j],sig,xi))
      }
      Qtrue[i,] <- temp
    }
  }else if(sim=="Sim2"){
    dw2p<-rdist(Sp,knots)^2
    dw2p[dw2p<0.0001]<-0
    d12<-rdist(Sp,S)
    d22<-rdist(Sp,Sp)
    d12[d12<0.0001]<-0
    diag(d22)<-0
    
    facp<-make.kern(dw2p,1)
    FACp<-stdKern(facp)
    
    chitrue <- matrix(0,nrow=5,ncol=np*(np-1)/2)
    Qtrue <- matrix(0,nrow=5,ncol=np)
    level <- c(0.5,0.9,0.95,0.99,0.995)
    for(i in 1:5){
      temp <- NULL
      for(j in 1:(np-1)){
        for(k in (j+1):np){
          qj = concdfinv(level[i],y,c(0.5,0.3,0.2),
                         t(log(lambda)),FACp,alpha,j)
          qk = concdfinv(level[i],y,c(0.5,0.3,0.2),
                         t(log(lambda)),FACp,alpha,k)
          temp = c(temp, Bconcdf(qj,qk,y,c(0.5,0.3,0.2),
                                 t(log(lambda)),FACp,alpha,j,k))
        }
      }
      chitrue[i,] <- (1-2*level[i]+temp)/(1-level[i])
    }
    for(i in 1:5){
      temp <- NULL
      for(j in 1:np){
        temp <- c(temp, concdfinv(level[i],y,c(0.5,0.3,0.2),
                                  t(log(lambda)),FACp,alpha,j))
      }
      Qtrue[i,] <- temp
    }
  }else if(sim=="Sim3"){
    chitrue <- matrix(0,nrow=5,ncol=np*(np-1)/2)
    Qtrue <- matrix(0,nrow=5,ncol=np)
    level <- c(0.5,0.9,0.95,0.99,0.995)
    for(i in 1:5){
      temp <- NULL
      for(j in 1:(np-1)){
        for(k in (j+1):np){
          qj = qnorm(level[i],mu[j],1)
          qk = qnorm(level[i],mu[k],1)
          temp = c(temp, 
                   pmvnorm(upper=c(qj,qk),mean=mu[c(j,k)],
                           sigma=matrix(c(C[j,j],C[j,k],C[k,j],C[k,k]),2,2)))
        }
      }
      chitrue[i,] <- (1-2*level[i]+temp)/(1-level[i])
    }
    for(i in 1:5){
      temp <- NULL
      for(j in 1:np){
        temp <- c(temp, qnorm(level[i],mu[j],1))
      }
      Qtrue[i,] <- temp
    }
    
  }else if(sim=="Sim4"){
    chitrue <- matrix(0,nrow=5,ncol=np*(np-1)/2)
    Qtrue <- matrix(0,nrow=5,ncol=np)
    level <- c(0.5,0.9,0.95,0.99,0.995)
    for(i in 1:5){
      temp <- NULL
      for(j in 1:(np-1)){
        for(k in (j+1):np){
          qj = qst(level[i],xi=mu[j],omega=sqrt(2/8*10),nu=8,alpha=3)
          qk = qst(level[i],xi=mu[k],omega=sqrt(2/8*10),nu=8,alpha=3)
          sigma=matrix(c(C[j,j],C[j,k],C[k,j],C[k,k]),2,2)
		  omega=(sigma+9)/(1+9*sum(solve(sigma)))^2*4
          #alpha=3*(1+9)^(0.5)*(1+9*sum(solve(sigma)))^(-0.5)*rowSums(solve(sigma))
          alpha=3*colSums(solve(sigma))
		  temp = c(temp, pmst(c(qj,qk),xi=mu[c(j,k)],
                              Omega=(sigma+9)*2/8,
                              alpha=alpha,nu=8))
        }
      }
      chitrue[i,] <- (1-2*level[i]+temp)/(1-level[i])
    }
    for(i in 1:5){
      temp <- NULL
      for(j in 1:np){
        temp <- c(temp, qst(level[i],xi=mu[j],omega=sqrt(2/8*10),nu=8,alpha=3))
      }
      Qtrue[i,] <- temp
    }
    
  }
  
  return(list(chitrue=chitrue,Qtrue=Qtrue))  
}

unif<-function(x){
  rank(x)/(length(x) + 1)
}

