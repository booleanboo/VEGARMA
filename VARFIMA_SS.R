library(multiwave)
library(GeneCycle)
library(arfima)
library(FKF)


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
  sim<-matrix(c(xt[2000:2500],yt[2000:2500]),ncol=2)
  i=20
  arfimass1 <- function(d1, d2,rho) {
    k=i
    pi<-matrix(NA,nrow=k,ncol=2)
    pi[1,1]<- d1
    pi[1,2]<- d2
    for(j in 2:k)
    {
      pi[j,1]<- (j-1-d1)/j*pi[j-1,1]
      pi[j,2]<- (j-1-d2)/j*pi[j-1,2]
    }
    pi<-pi
    combopi<-cbind(pi[,1],0,0,pi[,2])
    
    Tt1 <- cbind(as.vector(t(combopi[,1:2])),as.vector(t(combopi[,3:4])))
    Tt2<-rbind(diag(2*(k-1)),matrix(0,nrow=2,ncol=2*(k-1)))
    Tt<-cbind(Tt1,Tt2)
    Zt <-cbind(diag(2),matrix(0,nrow=2,ncol=2*(k-1)))
    ct <- matrix(0, nrow=2)
    dt <- matrix(0, nrow = 2*k)
    GGt <- matrix(0,nrow=2,ncol=2)
    H <- rbind(diag(2),matrix(0,nrow=2*(k-1),ncol=2))%*% matrix(c(1,rho,0,sqrt(1-rho^2)),ncol=2)
    HHt <- H %*% t(H)
    a0 <- rep(0,2*k)
    P0 <- diag(1e+07,2*k,2*k)
    return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                HHt = HHt))
  }
  
  objective <- function(theta, yt) {
    sp <- arfimass1(theta[1], theta[2],theta[3])
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt )
    return(-ans$logLik)
  }
  
  
  theta <- c(0.1,0.1,0.1)
  objective(theta,yt=t(sim))
  fit <- optim(theta, objective, yt = t(sim), method="L-BFGS-B", lower=c(0, 0, 0), upper=c(0.5,0.5,0.99), hessian = TRUE)
  fit$par