library(MTS)
library(FKF)
library(orthopolynom)
library(MASS)


Gest<-NULL

for(z in 1:1000)
  
{
  estd1<-vector()
  estu1<-vector()
  estd2<-vector()
  estu2<-vector()
  estrho<-vector()
  aic<-vector()
  
  ##############################################################################
  unnormalized.p.list1 <- gegenbauer.polynomials(50, 0.30, normalized=FALSE)
  theta1<-polynomial.values(unnormalized.p.list1,0.6)
  unnormalized.p.list2 <- gegenbauer.polynomials(50, 0.40, normalized=FALSE)
  theta2<-polynomial.values(unnormalized.p.list2,0.8)
  

  N=100 #length of simulated series
  Sigma <- matrix(c(1,0.8,0.8,1),ncol=2)
  at<-mvrnorm(n = 5000, rep(0, 2), Sigma, empirical = TRUE)
  xt<-filter(at[,1],theta1,method="convolution",sides=1)
  yt<-filter(at[,2],theta2,method="convolution",sides=1)
  sim<-matrix(c(xt[1000:(1000+N-1)],yt[1000:(1000+N-1)]),ncol=2)
  
  ###############VAR(1)#########################################################
  
  i=1
  
  VAR1approx <- function(d1, u1, d2, u2,rho) {
    phi1 <- 2*d1*u1
    phi2 <- 2*d2*u2
    Tt <- matrix(c(phi1, 0, 0, phi2), ncol = 2)
    Zt <- diag(2)
    ct <- matrix(0, nrow=2)
    dt <- matrix(0, nrow = 2)
    GGt <- matrix(0,nrow=2,ncol=2)
    HHt <- matrix(c(1,rho,rho,1),2,2)
    a0<-c(0,0)
    P0 <- matrix(c(1e+07,0e+00, 0e+00, 1e+07),nrow = 2, ncol = 2)
    return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                HHt = HHt))
  }
  
  ## The objective function passed to 'optim'
  objective <- function(theta, yt) {
    sp <- VAR1approx(theta[1], theta[2], theta[3], theta[4], theta[5])
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt )
    return(-ans$logLik)
  }
  
  theta <- c(0.1,0.2,0.25,0.55,0.1)
  objective(theta,t(sim))
  fit <- optim(theta, objective, yt = t(sim), method="L-BFGS-B", lower=rep(0.01,5), upper=c(0.5, 1, 0.5, 1, 0.99), hessian = TRUE)
  estd1[i]=fit$par[1]
  estu1[i]=fit$par[2]
  estd2[i]=fit$par[3]
  estu2[i]=fit$par[4]
  estrho[i]=fit$par[5]
  aic[i]=10+2*fit$value
  
  #############VAR(2)############################################################
  i=2
  
  VAR2approx <- function(d1, u1, d2, u2,rho) {
    pi11<- 2*d1*u1
    pi12<- 2*d2*u2
    pi21<- -d1-2*u1^2*d1*(d1-1)
    pi22<- -d2-2*u2^2*d2*(d2-1)
    Tt <- matrix(c(pi11,0,pi21,0,0,pi12,0,pi21,1,0,0,0,0,1,0,0), ncol = 4)
    Zt <- matrix(c(1,0,0,1,0,0,0,0),ncol=4)
    ct <- matrix(0, nrow=2)
    dt <- matrix(0, nrow = 4)
    GGt <- matrix(0,nrow=2,ncol=2)
    H <- matrix(c(1,0,0,0,0,1,0,0), nrow = 4)%*% matrix(c(1,rho,0,sqrt(1-rho^2)),ncol=2)
    HHt <- H %*% t(H)
    a0 <- c(0, 0, 0, 0)
    P0 <- diag(1e+07,4,4)
    return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                HHt = HHt))
  }
  
  ## The objective function passed to 'optim'
  objective <- function(theta, yt) {
    sp <- VAR2approx(theta[1], theta[2], theta[3], theta[4],theta[5])
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt )
    
    return(-ans$logLik)
  }
  
  theta <- c(0.1,0.2,0.25,0.55,0.1)
  objective(theta,t(sim))
  fit <- optim(theta, objective, yt = t(sim), method="L-BFGS-B", lower=rep(0.01,5), upper=c(0.5, 1, 0.5, 1, 0.99), hessian = TRUE)
  estd1[i]=fit$par[1]
  estu1[i]=fit$par[2]
  estd2[i]=fit$par[3]
  estu2[i]=fit$par[4]
  estrho[i]=fit$par[5]
  aic[i]=10+2*fit$value
  
  
  ##############################################################################
  
  
  for(i in 3:40)
  { 
    VARapprox<- function(d1,u1,d2,u2,rho) {
      k=i
      pi<-matrix(NA,nrow=k,ncol=2)
      pi[1,1]<--2*d1*u1
      pi[1,2]<--2*d2*u2
      pi[2,1]<-2*d1^2*u1^2-2*d1*u1^2+d1
      pi[2,2]<-2*d2^2*u2^2-2*d2*u2^2+d2
      for(j in 3:k)
      {
        pi[j,1]=2*u1*((-d1-1)/j+1)*pi[j-1,1]-(2*(-d1-1)/j+1)*pi[j-2,1]
        pi[j,2]=2*u2*((-d2-1)/j+1)*pi[j-1,2]-(2*(-d2-1)/j+1)*pi[j-2,2]
        
      } 
      pi=-pi
      combopi<-cbind(pi[,1],0,0,pi[,2])
      
      Tt1 <- t(matrix(as.vector(t(combopi)), ncol = 2*k))
      Tt2 <- rbind(diag(2*(k-1)),matrix(rep(0,4*(k-1)),nrow = 2))
      Tt  <- cbind(Tt1,Tt2)
      Zt <- matrix(c(1,0,0,1,rep(0,4*(k-1))),ncol=2*k)
      ct <- matrix(0, nrow=2)
      dt <- matrix(0, nrow = 2*k)
      GGt <- matrix(0,nrow=2,ncol=2)
      H <- rbind(diag(2),matrix(c(rep(0,4*(k-1))),ncol=2))%*% matrix(c(1,rho,0,sqrt(1-rho^2)),ncol=2)
      HHt <- H %*% t(H)
      a0 <- rep(0,2*k)
      P0 <- diag(1e+07,2*k,2*k)
      return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                  HHt = HHt))
    }
    
    ## The objective function passed to 'optim'
    objective <- function(theta, yt) {
      sp <- VARapprox(theta[1], theta[2],theta[3],theta[4],theta[5])
      ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
                 Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt )
      return(-ans$logLik)
    }
    
    theta <- c(0.1,0.2,0.25,0.55,0.1)
    objective(theta,yt=t(sim))
    fit <- optim(theta, objective, method="L-BFGS-B", lower=rep(0.01,5), upper=c(0.5,1,0.5,1,0.99),yt = t(sim), hessian = TRUE)
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


write.csv(Gest,"VARest.csv")



