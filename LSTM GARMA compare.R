library(keras)
library(tensorflow)

library(LongMemoryTS)
library(MTS)
library(FKF)
library(orthopolynom)
library(MASS)
library(longmemo)
library(ggplot2)
library(ggfortify)
library(forecast)

start_time <- Sys.time()

MAE1<-vector()
MAE2<-vector()
MSE1<-vector()
MSE2<-vector()
for(j in 1:500)
{
unnormalized.p.list <- gegenbauer.polynomials(50, 0.4, normalized=FALSE)
theta<-polynomial.values(unnormalized.p.list,0.8)

at<-rnorm(2500,sd=1)
#at<-rt(5000,df=3)
xt<-stats::filter(at,theta,method="convolution",sides=1)
data<-xt[1000:1499]



lag_transform <- function(x, k= 1){
  
  lagged =  c(rep(NA, k), x[1:(length(x)-k)])
  DF = as.data.frame(cbind(lagged, x))
  colnames(DF) <- c( paste0('x-', k), 'x')
  DF[is.na(DF)] <- 0
  return(DF)
}

supervised = lag_transform(data, 1)
#head(supervised)

N = nrow(supervised)
n = round(N *0.7, digits = 0)
train = supervised[1:n, ]
test  = supervised[(n+1):N,  ]

#nrow(train)
#nrow(test)



## scale data
scale_data = function(train, test, feature_range = c(0, 1)) {
  x = train
  fr_min = feature_range[1]
  fr_max = feature_range[2]
  std_train = ((x - min(x) ) / (max(x) - min(x)  ))
  std_test  = ((test - min(x) ) / (max(x) - min(x)  ))
  
  scaled_train = std_train *(fr_max -fr_min) + fr_min
  scaled_test = std_test *(fr_max -fr_min) + fr_min
  
  return( list(scaled_train = as.vector(scaled_train), scaled_test = as.vector(scaled_test) ,scaler= c(min =min(x), max = max(x))) )
  
}


Scaled = scale_data(train, test, c(-1, 1))

y_train = Scaled$scaled_train[, 2]
x_train = Scaled$scaled_train[, 1]

y_test = Scaled$scaled_test[, 2]
x_test = Scaled$scaled_test[, 1]

## inverse-transform
invert_scaling = function(scaled, scaler, feature_range = c(0, 1)){
  min = scaler[1]
  max = scaler[2]
  t = length(scaled)
  mins = feature_range[1]
  maxs = feature_range[2]
  inverted_dfs = numeric(t)
  
  for( i in 1:t){
    X = (scaled[i]- mins)/(maxs - mins)
    rawValues = X *(max - min) + min
    inverted_dfs[i] <- rawValues
  }
  return(inverted_dfs)
}


# Reshape the input to 3-dim
dim(x_train) <- c(length(x_train), 1, 1)

# specify required arguments
X_shape2 = dim(x_train)[2]
X_shape3 = dim(x_train)[3]
batch_size = 1                # must be a common factor of both the train and test samples
units = 1                     # can adjust this, in model tuninig phase

#=========================================================================================
#batch_input_shape (batch_size, num_steps, features)


model <- keras_model_sequential() 
model%>%
  layer_lstm(units, batch_input_shape = c(batch_size, X_shape2, X_shape3), stateful= TRUE)%>%
  layer_dense(units = 1)


model %>% compile(
  loss = 'mean_squared_error',
  optimizer = optimizer_adam( lr= 0.02, decay = 1e-6 ),  
  metrics = c('accuracy')
)
summary(model)

#Epochs = 50  
Epochs = 5
for(i in 1:Epochs ){
  model %>% fit(x_train, y_train, epochs=1, batch_size=batch_size, verbose=1, shuffle=FALSE)
  model %>% reset_states()
}

L = length(x_test)
scaler = Scaled$scaler
predictions = numeric(L)

for(i in 1:L){
  X = x_test[i]
  dim(X) = c(1,1,1)
  yhat = model %>% predict(X, batch_size=batch_size)
  # invert scaling
  yhat = invert_scaling(yhat, scaler,  c(-1, 1))
  # invert differencing
  #yhat  = yhat + Series[(n+i)]
  # store
  predictions[i] <- yhat
}


plot(data[351:500],type="l")
lines(predictions, col="red")



#############whittle##########
T<-length(data[1:350])
peri <- per(data[1:350])[-1]
#peri<-spec.pgram(data)$spec[-1]
m=floor((T-1)/2)
#m=floor(1+T^0.65)

Garma.whittle<-function (d,u) 
{ 
  lambda <- 2 * pi/T
  K <- sum(peri[1:m] * (2*abs(cos(lambda * (1:m))-u))^(2*d)) 
  K
}

objective <- function(theta) {
  sp <-Garma.whittle(theta[1], theta[2])
  return(sp)
}

theta<-c(0.1,0.1)
objective(theta)
fit <- optim(theta, objective, method="L-BFGS-B", lower=c(0, 0), upper=c(0.5, 1), hessian = TRUE)
fit$par

#Q.T<-Garma.whittle(fit$par[1],fit$par[2])
Q.T<-objective(fit$par)
sigma.T<-sqrt(4*pi*Q.T/T)




##########Forecast###############################


unnormalized.p.list1 <- gegenbauer.polynomials(51, -fit$par[1], normalized=FALSE)
theta1<-polynomial.values(unnormalized.p.list1,fit$par[2])
c1<-as.numeric(theta1[2:51])
c1<--c1


fc1<-vector()
for(i in 1:150)
{
  fc1[i]=t(rev(c1))%*%data[(300+i):(349+i)]
}
plot(data[351:500],type="l")
lines(fc1, col="red")

MAE1[j]<-sum(abs(fc1-data[351:500]))/150
MSE1[j]<-sum((fc1-data[351:500])^2)/150

MAE2[j]<-sum(abs(predictions-data[351:500]))/150
MSE2[j]<-sum((predictions-data[351:500])^2)/150


}

write.csv(cbind(MAE1,MSE1,MAE2,MSE2),"GARMA_whittle_LSTM_0408.csv")

end_time <- Sys.time()

end_time - start_time