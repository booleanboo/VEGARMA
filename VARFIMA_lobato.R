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


  theta1<-vector()
  theta1[1]<-1
  theta1[2]<-0.2
  
  for(j in 3:1000)
  { 
    theta1[j]=theta1[j-1]*(j-2+0.2)/(j-1)
    
  }
  theta2<-vector()
  theta2[1]<-1
  theta2[2]<-0.4
  
  for(k in 3:1000)
  { 
    theta2[k]=theta2[k-1]*(k-2+0.4)/(k-1)
    
  }
  
  Sigma <- matrix(c(1,0.8,0.8,1),ncol=2)
  at<-mvrnorm(n = 3000, rep(0, 2), Sigma, empirical = TRUE)
  xt<-filter(at[,1],theta1,method="convolution",sides=1)
  yt<-filter(at[,2],theta2,method="convolution",sides=1)
  
  
  
  n <- 500 
  m<-n^0.65
  lambdaj <- 2 * pi * (1:m)/n
  loglambdaj<-sum(log(lambdaj))/m
  G_d<-matrix(rep(0,4),2,2)
  
  varfima.whittle<-function(d1,d2)
  {
    weight.mat <- matrix(NA, n, m)
    for (k in 1:m) {
      weight.mat[, k] <- exp((0+1i) * (1:n) * lambdaj[k])
    }
    wx <- 1/sqrt(2 * pi * n) * xt[1000:1499] %*% weight.mat
    wy <- 1/sqrt(2 * pi * n) * yt[1000:1499] %*% weight.mat
    
    
    for(j in 1:m)
    {
      
      
      Lamdaj<-matrix(c(lambdaj[j]^d1,0,0,lambdaj[j]^d2),2,2)
      
      
      ###periodogram matrix
      
      
      wz<-matrix(c(wx[,j],wy[,j]),2,1)
      I.lambda <- wz%*%Conj(t(wz))
      
      G_d<-Re(Lamdaj%*%I.lambda%*%Lamdaj) +G_d
      
    }
    R_d<-log(det(G_d/m))-2*(d1*loglambdaj+d2*loglambdaj) 
    
    
    R_d
  }
  
  
  objective <- function(theta) {
    sp <-varfima.whittle(theta[1], theta[2])
    return(sp)
  }
  
  theta<-c(0.1,0.1)
  objective(theta)
  fit2 <- optim(theta, objective, method="L-BFGS-B", lower=c(0, 0), upper=c(0.5, 0.5), hessian = TRUE)
  fit2$par
  
  
  ######G_0################
  G_d<-matrix(rep(0,4),2,2)
  d1<-fit2$par[1]
  d2<-fit2$par[2]
  
  weight.mat <- matrix(NA, n, m) #500*m
  for (k in 1:m) {
    weight.mat[, k] <- exp((0+1i) * (1:n) * lambdaj[k]) 
  }
  wx <- 1/sqrt(2 * pi * n) * xt[1000:1499] %*% weight.mat #1*56
  wy <- 1/sqrt(2 * pi * n) * yt[1000:1499] %*% weight.mat  #1*56
  
  
  
  for(j in 1:m)
  {
    
    Lamdaj<-matrix(c(lambdaj[j]^d1,0,0,lambdaj[j]^d2),2,2)
    
    ###periodogram matrix
    
    wz<-matrix(c(wx[,j],wy[,j]),2,1)
    
    I.lambda <- wz%*%Conj(t(wz))
    
    G_d<-Re(Lamdaj%*%I.lambda%*%Lamdaj) +G_d
    
  }
  
  
  G_0<-G_d/m*2*pi
  G_0
  
  