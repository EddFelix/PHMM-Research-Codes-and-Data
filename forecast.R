setwd("~/Documents/PHMM")

library(nltsa)
library(ggplot2)
library(depmixS4)
library(readxl)

data = read_excel("Frequencies 1960-2020.xlsx", sheet="30-day interval")
earth = data.frame(data$to, data$`m>4`)

X=data$`m>4`[1:719]
actual = data$`m>4`[719:737]
n=length(actual)
date = data$from[719:737]

lam1 = 1.645341000
forecast_1 = rep(0,n)
forecast_1[1] = actual[1]
forecast_1[2:n]=lam1

forecast_2 = rep(0,n)

# 2-states
set.seed(90210)
model <- depmix(X ~1, nstates=2, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
if (u[7] <= u[8]) { para.mle = c(u[3:6], exp(u[7]), exp(u[8]))
} else { para.mle = c(u[6:3], exp(u[8]), exp(u[7]))
}
A2 = matrix(para.mle[1:4], byrow=TRUE, nrow=2)
lams2 = para.mle[5:6]
print(lams2)
print(A2)

residual_2 = matrix(0,nrow=2, ncol=n)
for (i in 1:n){
  residual_2[1,i]=abs(actual[i]-lams2[1])
  residual_2[2,i]=abs(actual[i]-lams2[2])
}

forecast_2[2:n]=lams2[1]*A2[1,1]+lams2[2]*A2[1,2]

#three states
set.seed(90210)
model <- depmix(X ~1, nstates=3, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[4:12], exp(u[13]), exp(u[14]), exp(u[15]))

A3 = matrix(para.mle[1:9], byrow=TRUE, nrow=3)
lams3 = rev(para.mle[10:12])
print(lams3)
print(A3)

residual_3 = matrix(0,nrow=3, ncol=n)
rank3 = matrix(0,nrow=3, ncol=n)
forecast_3 = rep(0,n)
for (i in 1:n){
  residual_3[1,i]=abs(actual[i]-lams3[1])
  residual_3[2,i]=abs(actual[i]-lams3[2])
  residual_3[3,i]=abs(actual[i]-lams3[3])
  rank3[,i] = rank(residual_3[,i])
}

for (i in 2:n){
  if (rank3[1,(i-1)]==1){
    forecast_3[i] = lams3[1]*A3[1,1]+lams3[2]*A3[1,2]+lams3[3]*A3[1,3]
  }
  else if (rank3[2,(i-1)]==1){
    forecast_3[i] = lams3[1]*A3[2,1]+lams3[2]*A3[2,2]+lams3[3]*A3[2,3]
  }
  else {
    forecast_3[i] = lams3[1]*A3[3,1]+lams3[2]*A3[3,2]+lams3[3]*A3[3,3]
  }
}

# 4-states
set.seed(90210)
model <- depmix(X ~1, nstates=4, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[5:20], exp(u[21]), exp(u[22]), exp(u[23]), exp(u[24]))

A4 = matrix(para.mle[1:16], byrow=TRUE, nrow=4)
lams4 = para.mle[17:20]
print(lams4)
print(A4)

residual_4 = matrix(0,nrow=4, ncol=n)
rank4 = matrix(0,nrow=4, ncol=n)
forecast_4 = rep(0,n)
for (i in 1:n){
  residual_4[1,i]=abs(actual[i]-lams4[1])
  residual_4[2,i]=abs(actual[i]-lams4[2])
  residual_4[3,i]=abs(actual[i]-lams4[3])
  residual_4[4,i]=abs(actual[i]-lams4[4])
  rank4[,i] = rank(residual_4[,i])
}

for (i in 2:n){
  if (rank4[1,(i-1)]==1){
    forecast_4[i] = lams4[1]*A4[1,1]+lams4[2]*A4[1,2]+lams4[3]*A4[1,3]+lams4[4]*A4[1,4]
  }
  else if (rank4[2,(i-1)]==1){
    forecast_4[i] = lams4[1]*A4[2,1]+lams4[2]*A4[2,2]+lams4[3]*A4[2,3]+lams4[4]*A4[2,4]
  }
  else if (rank4[3,(i-1)]==1){
    forecast_4[i] = lams4[1]*A4[3,1]+lams4[2]*A4[3,2]+lams4[3]*A4[3,3]+lams4[4]*A4[3,4]
  }
  else {
    forecast_4[i] = lams4[1]*A4[4,1]+lams4[2]*A4[4,2]+lams4[3]*A4[4,3]+lams4[4]*A4[4,4]
  }
}

# 5-states
set.seed(90210)
model <- depmix(X ~1, nstates=5, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[6:30], exp(u[31]), exp(u[32]), exp(u[33]), exp(u[34]), exp(u[35]))

A5 = matrix(para.mle[1:25], byrow=TRUE, nrow=5)
lams5 = para.mle[26:30]
print(lams5)
print(A5)

residual_5 = matrix(0,nrow=5, ncol=n)
rank5 = matrix(0,nrow=5, ncol=n)
forecast_5 = rep(0,n)
for (i in 1:n){
  residual_5[1,i]=abs(actual[i]-lams5[1])
  residual_5[2,i]=abs(actual[i]-lams5[2])
  residual_5[3,i]=abs(actual[i]-lams5[3])
  residual_5[4,i]=abs(actual[i]-lams5[4])
  residual_5[5,i]=abs(actual[i]-lams5[5])
  rank5[,i] = rank(residual_5[,i])
}

for (i in 2:n){
  if (rank5[1,(i-1)]==1){
    forecast_5[i] = lams5[1]*A5[1,1]+lams5[2]*A5[1,2]+lams5[3]*A5[1,3]+lams5[4]*A5[1,4]+lams5[5]*A5[1,5]
  }
  else if (rank5[2,(i-1)]==1){
    forecast_5[i] = lams5[1]*A5[2,1]+lams5[2]*A5[2,2]+lams5[3]*A5[2,3]+lams5[4]*A5[2,4]+lams5[5]*A5[2,5]
  }
  else if (rank5[3,(i-1)]==1){
    forecast_5[i] = lams5[1]*A5[3,1]+lams5[2]*A5[3,2]+lams5[3]*A5[3,3]+lams5[4]*A5[3,4]+lams5[5]*A5[3,5]
  }
  else if (rank5[4,(i-1)]==1){
    forecast_5[i] = lams5[1]*A5[4,1]+lams5[2]*A5[4,2]+lams5[3]*A5[4,3]+lams5[4]*A5[4,4]+lams5[5]*A5[4,5]
  }
  else {
    forecast_5[i] = lams5[1]*A5[5,1]+lams5[2]*A5[5,2]+lams5[3]*A5[5,3]+lams5[4]*A5[5,4]+lams5[5]*A5[5,5]
  }
}

# 6-states
set.seed(90210)
model <- depmix(X ~1, nstates=6, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[7:42], exp(u[43]), exp(u[44]), exp(u[45]), exp(u[46]), exp(u[47]), exp(u[48]))

A6 = matrix(para.mle[1:36], byrow=TRUE, nrow=6)
lams6 = para.mle[37:42]
print(lams6)
print(A6)

residual_6 = matrix(0,nrow=6, ncol=n)
rank6 = matrix(0,nrow=6, ncol=n)
forecast_6 = rep(0,n)
for (i in 1:n){
  residual_6[1,i]=abs(actual[i]-lams6[1])
  residual_6[2,i]=abs(actual[i]-lams6[2])
  residual_6[3,i]=abs(actual[i]-lams6[3])
  residual_6[4,i]=abs(actual[i]-lams6[4])
  residual_6[5,i]=abs(actual[i]-lams6[5])
  residual_6[6,i]=abs(actual[i]-lams6[6])
  rank6[,i] = rank(residual_6[,i])
}

for (i in 2:n){
  if (rank6[1,(i-1)]==1){
    forecast_6[i] = lams6[1]*A6[1,1]+lams6[2]*A6[1,2]+lams6[3]*A6[1,3]+lams6[4]*A6[1,4]+lams6[5]*A6[1,5]+lams6[6]*A6[1,6]
  }
  else if (rank6[2,(i-1)]==1){
    forecast_6[i] = lams6[1]*A6[2,1]+lams6[2]*A6[2,2]+lams6[3]*A6[2,3]+lams6[4]*A6[2,4]+lams6[5]*A6[2,5]+lams6[6]*A6[2,6]
  }
  else if (rank6[3,(i-1)]==1){
    forecast_6[i] = lams6[1]*A6[3,1]+lams6[2]*A6[3,2]+lams6[3]*A6[3,3]+lams6[4]*A6[3,4]+lams6[5]*A6[3,5]+lams6[6]*A6[3,6]
  }
  else if (rank6[4,(i-1)]==1){
    forecast_6[i] = lams6[1]*A6[4,1]+lams6[2]*A6[4,2]+lams6[3]*A6[4,3]+lams6[4]*A6[4,4]+lams6[5]*A6[4,5]+lams6[6]*A6[4,6]
  }
  else if (rank5[5,(i-1)]==1){
    forecast_6[i] = lams6[1]*A6[5,1]+lams6[2]*A6[5,2]+lams6[3]*A6[5,3]+lams6[4]*A6[5,4]+lams6[5]*A6[5,5]+lams6[6]*A6[5,6]
  }
  else {
    forecast_6[i] = lams6[1]*A6[6,1]+lams6[2]*A6[6,2]+lams6[3]*A6[6,3]+lams6[4]*A6[6,4]+lams6[5]*A6[6,5]+lams6[6]*A6[6,6]
  }
}

#arima
arima=rep(0,n)
arima[2] = 1.41789
arima[3] = 1.8812
arima[4] = 2.614079
arima[5] = 2.452054
arima[6] = 2.493527
arima[7] = 1.661978
arima[8] = 1.391435
arima[9] = 1.84189
arima[10] = 4.263278
arima[11] = 1.700743
arima[12] = 3.064235
arima[13] = 1.610838
arima[14] = 4.969254
arima[15] = 3.812922
arima[16] = 1.666197
arima[17] = 2.928414
arima[18] = 5.233352
arima[19] = 1.876134

plot(date, actual, type="l", ylim=c(0,30), lwd=2)
lines(date, forecast_1, type="l", col="red")
lines(date, forecast_2, type="l", col="green")
lines(date, forecast_3, type="l", col="blue")
lines(date, forecast_4, type="l", col="orange")
lines(date, forecast_5, type="l", col="violet")
lines(date, forecast_6, type="l", col="brown", lty=2)
lines(date, arima, type="l", col="red", lwd=2, lty=2)
legend("topright", legend=c("actual", "1-state", "2-state", "3-state", "4-state", "5-state", "6-state"), col=c("black", "red", "green", "blue", "orange", "violet", "brown"), lty=1, cex=0.5)


#RMSE
rmse1=sqrt(mean(abs(actual-forecast_1)^2))
print(rmse1)
rmse2=sqrt(mean(abs(actual-forecast_2)^2))
print(rmse2)
rmse3=sqrt(mean(abs(actual-forecast_3)^2))
print(rmse3)
rmse4=sqrt(mean(abs(actual-forecast_4)^2))
print(rmse4)
rmse5=sqrt(mean(abs(actual-forecast_5)^2))
print(rmse5)
rmse6=sqrt(mean(abs(actual-forecast_6)^2))
print(rmse6)
rmse_arima=sqrt(mean(abs(actual-arima)^2))
print(rmse_arima)

#MAE
mae1 = mean(abs(actual-forecast_1))
print(mae1)
mae2 = mean(abs(actual-forecast_2))
print(mae2)
mae3 = mean(abs(actual-forecast_3))
print(mae3)
mae4 = mean(abs(actual-forecast_4))
print(mae4)
mae5 = mean(abs(actual-forecast_5))
print(mae5)
mae6 = mean(abs(actual-forecast_6))
print(mae6)
mae_arima = mean(abs(actual-arima))
print(mae_arima)

#RAE
rae1 = (sum(abs(actual-forecast_1)))/(sum(abs(actual-mean(forecast_1))))
print(rae1)
rae2 = (sum(abs(actual-forecast_2)))/(sum(abs(actual-mean(forecast_2))))
print(rae2)
rae3 = (sum(abs(actual-forecast_3)))/(sum(abs(actual-mean(forecast_3))))
print(rae3)
rae4 = (sum(abs(actual-forecast_4)))/(sum(abs(actual-mean(forecast_4))))
print(rae4)
rae5 = (sum(abs(actual-forecast_5)))/(sum(abs(actual-mean(forecast_5))))
print(rae5)
rae6 = (sum(abs(actual-forecast_6)))/(sum(abs(actual-mean(forecast_6))))
print(rae6)
rae_arima = (sum(abs(actual-arima)))/(sum(abs(actual-mean(arima))))
print(rae_arima)
