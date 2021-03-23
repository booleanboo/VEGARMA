#system.time(source("scriptname.R"))
library(MTS)
library(FKF)
library(orthopolynom)
library(MASS)


Gest<-NULL


for(k in 1:1000)
{
  
  estd1<-vector()
  estu1<-vector()
  estd2<-vector()
  estu2<-vector()
  estrho<-vector()
  aic<-vector()
  
 
  ##############################################################################
  unnormalized.p.list1 <- gegenbauer.polynomials(10, 0.30, normalized=FALSE)
  theta1<-polynomial.values(unnormalized.p.list1,0.6)
  unnormalized.p.list2 <- gegenbauer.polynomials(10, 0.40, normalized=FALSE)
  theta2<-polynomial.values(unnormalized.p.list2,0.8)
  
  
  N=100 #length of simulated series
  Sigma <- matrix(c(1,0.8,0.8,1),ncol=2)
  at<-mvrnorm(n = 5000, rep(0, 2), Sigma, empirical = TRUE)
  xt<-filter(at[,1],theta1,method="convolution",sides=1)
  yt<-filter(at[,2],theta2,method="convolution",sides=1)
  sim<-matrix(c(xt[1000:(1000+N-1)],yt[1000:(1000+N-1)]),ncol=2)
  
  
  ################################################################################
  i=1
  
  VMA1approx <- function(d1,u1,d2,u2,rho) {
    si11<-2*d1*u1
    si12<-2*d2*u2
    Tt <- matrix(c(rep(0,8), 1, 0, 0, 0, 0, 1, 0, 0), ncol = 4)
    Zt <- matrix(c(1, 0, 0, 1,rep(0,4)), ncol = 4)
    ct <- matrix(0,nrow=2)
    dt <- matrix(0, nrow = 4)
    GGt <- matrix(0,nrow=2,ncol=2)
    H <- matrix(c(1,0,si11,0,0,1,0,si12), nrow = 4)%*% matrix(c(1,rho,0,sqrt(1-rho^2)),ncol=2)
    HHt <- H %*% t(H)
    a0 <- c(0, 0, 0, 0)
    P0<- P0 <- diag(1e+07,4,4)
    #P0<-toeplitz(c(si11^2+1,0,si11,0))
    return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                HHt = HHt))
  }
  
  
  ## The objective function passed to 'optim'
  objective <- function(theta, yt) {
    sp <- VMA1approx(theta[1], theta[2],theta[3],theta[4],theta[5])
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
    return(-ans$logLik)
  }
  
  
  theta <- c(0.1,0.2,0.25,0.55,0.1)
  objective(theta,yt=t(sim))
  fit <- optim(theta, objective, method="L-BFGS-B", lower=c(0,0,0,0,0), upper=c(0.5,1,0.5,1,0.99),yt = t(sim), hessian = TRUE)
  estd1[i]=fit$par[1]
  estu1[i]=fit$par[2]
  estd2[i]=fit$par[3]
  estu2[i]=fit$par[4]
  estrho[i]=fit$par[5]
  aic[i]=10+2*fit$value
  ###############################################################################
  
  i=2
  
  VMA2approx <- function(d1,u1,d2,u2,rho) {
    si11<-2*d1*u1
    si12<-2*d2*u2
    si21<-2*d1^2*u1^2+2*d1*u1^2-d1
    si22<-2*d2^2*u2^2+2*d2*u2^2-d2
    Tt <- matrix(c(rep(0,12), 1, rep(0,6), 1, rep(0,6),1,rep(0,6),1,0,0), ncol = 6)
    Zt <- matrix(c(1, 0, 0, 1,rep(0,8)), ncol = 6)
    ct <- matrix(0,nrow=2)
    dt <- matrix(0, nrow = 6)
    GGt <- matrix(0,nrow=2,ncol=2)
    H <- matrix(c(1,0,si11,0,si21,0,0,1,0,si12,0,si22), nrow = 6)%*% matrix(c(1,rho,0,sqrt(1-rho^2)),ncol=2)
    HHt <- H %*% t(H)
    a0 <- c(0, 0, 0, 0, 0, 0)
    P0 <- diag(1e+07,6,6)
    #P0<-toeplitz(c(1+si11^2+si21^2,0,si11+si11*si21,0,si21,0))
    return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                HHt = HHt))
  }
  
  
  ## The objective function passed to 'optim'
  objective <- function(theta, yt) {
    sp <- VMA2approx(theta[1], theta[2],theta[3],theta[4],theta[5])
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
    return(-ans$logLik)
  }
  
  
  theta <- c(0.1,0.2,0.25,0.55,0.1)
  #theta<-c(0.1,0.1,0.1,0.1,0.1)
  objective(theta,yt=t(sim))
  fit <- optim(theta, objective, method="L-BFGS-B", lower=c(0, 0, 0, 0, 0), upper=c(0.5,1,0.5,1,0.99),yt = t(sim), hessian = TRUE)
  estd1[i]=fit$par[1]
  estu1[i]=fit$par[2]
  estd2[i]=fit$par[3]
  estu2[i]=fit$par[4]
  estrho[i]=fit$par[5]
  aic[i]=10+2*fit$value
  
  ##############################################################################
  
  for(i in 3:40)
  { ## |u|<1 0<d<0.5
    ## Create a state space representation out of two parameters
    
    VMAapprox <- function(d1,u1,d2,u2,rho) {
      k=i
      si<-matrix(NA,nrow=k,ncol=2)
      si[1,1]<-2*d1*u1
      si[1,2]<-2*d2*u2
      si[2,1]<-2*d1^2*u1^2+2*d1*u1^2-d1
      si[2,2]<-2*d2^2*u2^2+2*d2*u2^2-d2
      for(j in 3:k)
      {
        si[j,1]=2*u1*((d1-1)/j+1)*si[j-1,1]-(2*(d1-1)/j+1)*si[j-2,1]
        si[j,2]=2*u2*((d2-1)/j+1)*si[j-1,2]-(2*(d2-1)/j+1)*si[j-2,2]
        
      } 
      
      Tt1 <-cbind(matrix(rep(0,4*k),ncol=2),diag(2*k))
      Tt2 <- matrix(rep(0,4*(k+1)),ncol=2*(k+1))
      Tt <-rbind(Tt1,Tt2)
      Zt <- matrix(c(1, 0, 0, 1,rep(0,4*k)), ncol = 2*(k+1))
      ct <- matrix(0,nrow=2)
      dt <- matrix(0, nrow = 2*(k+1))
      GGt <- matrix(0,nrow=2,ncol=2)
      combosi<-cbind(si[,1],0,0,si[,2])
      H <- matrix(c(1,0,as.vector(t(combosi[,1:2])),0,1,as.vector(t(combosi[,3:4]))), ncol=2)%*% matrix(c(1,rho,0,sqrt(1-rho^2)),ncol=2)
      HHt <- H %*% t(H)
      a0 <- rep(0,2*(k+1))
      P0 <- diag(1e+07,2*(k+1),2*(k+1))
      return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                  HHt = HHt))
    }
    
    
    ## The objective function passed to 'optim'
    objective <- function(theta, yt) {
      sp <- VMAapprox(theta[1], theta[2],theta[3],theta[4],theta[5])
      ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
                 Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
      #print(theta)
      return(-ans$logLik)
    }
    
    
    theta <- c(0.1,0.2,0.25,0.55,0.1)
    objective(theta,yt=t(sim))
    fit <- optim(theta, objective, method="L-BFGS-B", lower=c(0, 0, 0, 0, 0), upper=c(0.5,1,0.5,1,0.99),yt = t(sim), hessian = TRUE)
    estd1[i]=fit$par[1]
    estu1[i]=fit$par[2]
    estd2[i]=fit$par[3]
    estu2[i]=fit$par[4]
    estrho[i]=fit$par[5]
    aic[i]=10+2*fit$value
    
  }
  
  temp<-cbind(estd1,estu1,estd2,estu2,estrho)
  
  Gest<-rbind(Gest,temp)
  
}

write.csv(Gest, file = "VMAest.csv",row.names=TRUE)




