library(LongMemoryTS)
library(MTS)
library(FKF)
library(orthopolynom)
library(MASS)
library(longmemo)
library(beyondWhittle)
library(forecast)
library(nlme)
library(msos)


unnormalized.p.list1 <- gegenbauer.polynomials(50, 0.30, normalized=FALSE)
theta1<-polynomial.values(unnormalized.p.list1,0.6)
unnormalized.p.list2 <- gegenbauer.polynomials(50, 0.40, normalized=FALSE)
theta2<-polynomial.values(unnormalized.p.list2,0.8)
  
  
Sigma <- matrix(c(1,0.8,0.8,1),ncol=2)
  at<-mvrnorm(n = 3000, rep(0, 2), Sigma, empirical = TRUE)
  xt<-filter(at[,1],theta1,method="convolution",sides=1)
  yt<-filter(at[,2],theta2,method="convolution",sides=1)
  
  
  
  ###########################################################
  
  n <- 500 
  m<-floor(n/2)
  
  lambdaj <- 2 * pi * (1:m)/n
  #lambdaj <- 2 * pi * (0:(m-1))/n
  G_d<-matrix(rep(0,4),2,2)
  Q_d<-0
  
  varfima.whittle<-function(d1,u1,d2,u2)
  {
    
    weight.mat <- matrix(NA, n, m) #500*m
    for (k in 1:m) {
      weight.mat[, k] <- exp((0+1i) * (1:n) * lambdaj[k]) 
    }
    wx <- 1/sqrt(2 * pi * n) * xt[1000:1499] %*% weight.mat #1*56
    wy <- 1/sqrt(2 * pi * n) * yt[1000:1499] %*% weight.mat  #1*56
    
    
    
    for(j in 1:m)
    {
      
      c11<-abs(2*cos(lambdaj[j])-2*u1)^(-d1)*exp(-1i*(lambdaj[j]-pi)*d1)
      c22<-abs(2*cos(lambdaj[j])-2*u2)^(-d2)*exp(-1i*(lambdaj[j]-pi)*d2)
      
      
      Lamdaj<-matrix(c(c11,0,0,c22),2,2)
      
      
      ###periodogram matrix
      
      wz<-matrix(c(wx[,j],wy[,j]),2,1)
      
      I.lambda <- wz%*%Conj(t(wz))
      
      G_d<-Re(solve(Lamdaj)%*%I.lambda%*%solve(Conj(t(Lamdaj))))+G_d
      
      Q_d<-2*(d1*log(abs(2*cos(lambdaj[j])-2*u1))+d2*log(abs(2*cos(lambdaj[j])-2*u2))) +Q_d
      
      
    }
    
    R_d<-log(det(G_d/m))-Q_d/m
    
    R_d
  }
  
  
  objective <- function(theta) {
    sp <-varfima.whittle(theta[1],theta[2],theta[3],theta[4])
    return(sp) 
  }
  
  theta<-c(0.3,0.6,0.3,0.6)
  objective(theta)
  fit2 <- optim(theta, objective, method="L-BFGS-B", lower=rep(0, 4), upper=c(0.5,0.99,0.5,0.99), hessian = TRUE)
  fit2$par
  
  ######G_0################
  d1<-fit2$par[1]
  u1<-fit2$par[2]
  d2<-fit2$par[3]
  u2<-fit2$par[4]
  
  
  
  G_d<-matrix(rep(0,4),2,2)
  Q_d<-0
  
  weight.mat <- matrix(NA, n, m) #500*m
  for (k in 1:m) {
    weight.mat[, k] <- exp((0+1i) * (1:n) * lambdaj[k]) 
  }
  wx <- 1/sqrt(2 * pi * n) * xt[1000:1499] %*% weight.mat #1*56
  wy <- 1/sqrt(2 * pi * n) * yt[1000:1499] %*% weight.mat  #1*56
  
  
  
  for(j in 1:m)
  {
    c11<-abs(2*cos(lambdaj[j])-2*u1)^(-d1)*exp(-1i*(lambdaj[j]-pi)*d1)
    c22<-abs(2*cos(lambdaj[j])-2*u2)^(-d2)*exp(-1i*(lambdaj[j]-pi)*d2)
    
    
    Lamdaj<-matrix(c(c11,0,0,c22),2,2)
    
    
    ###periodogram matrix
    
    wz<-matrix(c(wx[,j],wy[,j]),2,1)
    
    I.lambda <- wz%*%Conj(t(wz))
    
    G_d<-Re(solve(Lamdaj)%*%I.lambda%*%solve(Conj(t(Lamdaj))))+G_d
    Q_d<-2*(d1*log(abs(2*cos(lambdaj[j])-2*u1))+d2*log(abs(2*cos(lambdaj[j])-2*u2))) +Q_d
   
    
  }
  
  
  G_0<-G_d/m*2*pi
  G_0
  
  
 
