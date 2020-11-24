setwd("C:/Users/Alvin/Desktop/PHMM/PHMM")
PHMM.data = read.csv("PHMM.csv")
PHMM.data = PHMM.data[,c(2,3,4)]
colnames(PHMM.data) <- c("M4","M5","M7") 
train.length = 719
forecast.length = 1
train.data = PHMM.data[1: train.length, ]
arima.model = auto.arima(train.data$M4)
arima.pred = predict(arima.model, n.ahead = forecast.length)
arima.pred
