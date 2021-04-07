library(LongMemoryTS)
library(MTS)
library(FKF)
library(orthopolynom)
library(MASS)
library(longmemo)
library(GeneCycle)

lag_transform <- function(x, k= 1){
  
  lagged =  c(rep(NA, k), x[1:(length(x)-k)])
  DF = as.data.frame(cbind(lagged, x))
  colnames(DF) <- c( paste0('x-', k), 'x')
  DF[is.na(DF)] <- 0
  return(DF)
}


unnormalized.p.list1 <- gegenbauer.polynomials(50, -0.4, normalized=FALSE)
theta1<-polynomial.values(unnormalized.p.list1,0.8)
c1<-as.numeric(theta1[2:50])
c1<--c1
  
unnormalized.p.list1 <- gegenbauer.polynomials(50, -0.2, normalized=FALSE)
theta2<-polynomial.values(unnormalized.p.list1,0.2)
c2<-as.numeric(theta2[2:50])
c2<--c2
  
  
at<-rnorm(5000,sd=1)
xt<-filter(at,c1,method="recursive")
yt<-filter(xt,c2,method="recursive")
sim<-yt[1000:1999]
  
 ##############AR approximation#########################
  
  i=51 

  kgarmafit <- function(d1,u1,k1,k2,sigma) {
      d2=k1*d1
      u2=k2*u1
      k=i-1
      j=(i+1)/2
      pi1=vector()
      pi1[1]<-1
      pi1[2]<--2*d2*u2
      pi1[3]<-2*d2^2*u2^2-2*d2*u2^2+d2
      for(g in 4:j)
      {
        pi1[g]=2*u2*((-d2-1)/(g-1)+1)*pi1[g-1]-(2*(-d2-1)/(g-1)+1)*pi1[g-2]
        
      } 
      
      pi2=vector()
      pi2[1]<-1
      pi2[2]<--2*d1*u1
      pi2[3]<-2*d1^2*u1^2-2*d1*u1^2+d1
      
      for(g in 4:j)
      {
        pi2[g]=2*u1*((-d1-1)/(g-1)+1)*pi2[g-1]-(2*(-d1-1)/(g-1)+1)*pi2[g-2]
        
      } 
      
      pi2[(j+1):i]<-0
      
      pi_temp<-matrix(NA, nrow=(i-1)/2, ncol=i)
      
      o<-(i-1)/2
      for(z in 1:o)
      { 
        pi_temp[z,]= lag_transform(pi1[z+1]*pi2, z)[,1]
      }
      
      pi_final<-pi2+colSums(pi_temp)
      pi_final=-pi_final[2:i]
      
      #AR(k)
      Tt <- cbind(pi_final,rbind(diag(k-1),rep(0,k-1)))
      Zt <- matrix(c(1, rep(0,k-1)), ncol = k)
      ct <- matrix(0)
      dt <- matrix(0, nrow = k)
      GGt <- matrix(0)
      H <- matrix(c(1,rep(0,k-1)), nrow = k)*sigma
      HHt <- H %*% t(H)
      a0 <- rep(0,k)
      P0 <- diag(1e+07,nrow = k)
      
      return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                  HHt = HHt))
    }
    
    ## The objective function passed to 'optim'
    objective <- function(theta, yt) {
      sp <- kgarmafit(theta[1], theta[2], theta[3], theta[4], theta[5])
      ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
                 Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
      return(-ans$logLik)
    }
    
    
    theta <- c(0.3,0.6,0.3,0.6,0.1)
    objective(theta,yt=rbind(sim))
    fit <- optim(theta, objective, method="L-BFGS-B", lower=c(0, 0, 0, 0, 0.01), upper=c(0.5,1,1,1,2),yt = rbind(sim), hessian = TRUE)
    