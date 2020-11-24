rm(list = ls())
setwd("C:/Users/Alvin/Desktop")
##' load required packages
require(tseries)
require(forecast)
require(dplyr)
##' read in data for the demo
PHMM.data = read.csv("PHMM.csv")
PHMM.data = PHMM.data[,c(2,3,4)]

colnames(PHMM.data) <- c("M4","M5","M3") 
#M4: magnitude 4
#M5: magnitude 5
#M3: magnitude 3

##' have a quick look at the data
plot(PHMM.data$M4)

##' define how much data we want to use to train the model and how far ahead we want to forecast
train.length = 572
forecast.length = 5

res1 <- matrix(0, nrow=11, ncol=2)     
colnames(res1) <- c("mdl1_mu","mdl1_sd")

for (k in 0:10){
  

##" split into train and test data
train.data = PHMM.data[1: train.length, ]
test.data = PHMM.data[(train.length+1) : (train.length + forecast.length), ]

############################################################

##################  LINEAR REGRESSION   ####################

############################################################

##' find potential lags for covariates to cases
a = ccf(train.data$M4, train.data$M5)
M5.lag = max(forecast.length, a$lag[which(a$acf == max(a$acf))])

a = ccf(train.data$M4, train.data$M3)
M3.lag = max(forecast.length, a$lag[which(a$acf == max(a$acf))])

##' generate columns for these lagged data
lagged.data = PHMM.data
lagged.data$M5.lag = lag(lagged.data$M5, M5.lag)
lagged.data$M3.lag = lag(lagged.data$M3, M3.lag)


##' generate lagged training and testing data set
train.lag.data = lagged.data[1:train.length, ]
test.lag.data = lagged.data[(train.length + 1):(train.length + forecast.length), ]



##' fit lagged lm model
lm.model.lag = lm(M4 ~ M5.lag + M3.lag, 
                  data = train.lag.data)

# if(k > 50)
#    browser()

##' make predictions
lm.pred.lag = predict(lm.model.lag, test.lag.data, se.fit = T)


##' make lower and upper prediction intervals
lower.lm.lag = lm.pred.lag$fit - qnorm(0.975)*sqrt(lm.pred.lag$se.fit^2 + lm.pred.lag$residual.scale^2)
upper.lm.lag = lm.pred.lag$fit + qnorm(0.975)*sqrt(lm.pred.lag$se.fit^2 + lm.pred.lag$residual.scale^2)

##' how to convert cases to logged cases
##' other transformations are available


#   log.cases.lag.data = lagged.data
#   log.cases.lag.data$Cases = log(lagged.data$Cases + 1)

############################################################

####################     ARIMA MODEL    ####################

############################################################

##' check stationarity of case data
adf.test(train.data$M4)  # stationarity if p-value < 0.05
kpss.test(train.data$M4) # stationarity if p-value > 0.05


##' possibly not stationary, try differencing
adf.test(diff(train.data$M4))
kpss.test(diff(train.data$M4))



##' looks to be stationary now
##' can we guess the order of our ARIMA function
acf(diff(train.data$M4))  
pacf(diff(train.data$M4)) 



##' use auto.arima to fit a model and predict forwards
##' this will automatically fit the differencing (d), along with p and q
arima.model = auto.arima(train.data$M4)



##' check which model is fit
arima.model



##' produce predictions 
arima.pred = predict(arima.model, n.ahead = forecast.length)



##' make lower and upper intervals
upper.arima = arima.pred$pred + qnorm(0.975)*arima.pred$se
lower.arima = arima.pred$pred - qnorm(0.975)*arima.pred$se

res1[1+k,1] = arima.pred$pred[forecast.length]
res1[1+k,2] = arima.pred$se[forecast.length]

}

############################################################

####################      PLOT     ####################

############################################################



##' plot a comparison of the data and the forecasts
##' plot the data first 
plot(c(tail(train.lag.data$M4,25), test.lag.data$M4), type = 'o', ylim = c(0, 10), pch = 16, 
     ylab = 'M4', bty = 'n')
abline(v = 20.5)

##' plot point predictions from the arima model
##' and make a polygon for lower and upper intervals
points(21:(20+forecast.length),as.numeric(arima.pred$pred), col = 'red', pch = 16, type = 'o')
polygon(x = c(21:(20+forecast.length),rev(21:(20+forecast.length))), y = c(upper.arima, rev(lower.arima)),
        col=adjustcolor('darkred', 0.25), border=F)

############################################################

####################      ERROR     ####################

############################################################

mse1=0

for (k in 0:10){
  mse1=((res1[1+k,2])^2)+mse1
}

tmse1=sqrt(mse1/11) #MSE for Mdl1

tmse1

write.csv(res1,file="5timeahead_alv.csv")
