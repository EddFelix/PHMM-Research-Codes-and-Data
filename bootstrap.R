setwd("~/Documents/PHMM")

library(nltsa)
library(ggplot2)
library(depmixS4)
library(readxl)

data = read_excel("Frequencies 1960-2020.xlsx", sheet="30-day interval")
earth = data.frame(data$to, data$`m>4`)
#earth = data.frame(data$X__2, data$freq1)
#par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
#layout(matrix(c(1,2, 1,3), nc=2))
plot(earth, type="h", xlab="date", ylab="frequency")
acf(earth)
pacf(earth)


X=data$`m>4`[1:719]
n=length(X)
# 5-states
set.seed(90210)
model <- depmix(X ~1, nstates=5, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[6:30], exp(u[31]), exp(u[32]), exp(u[33]), exp(u[34]), exp(u[35]))

mtrans = matrix(para.mle[1:25], byrow=TRUE, nrow=5)
lams = para.mle[26:30]
print(lams)
print(mtrans)

lams=sort(lams)

#declaring the variables
simul = 10000
k=30
run=floor(n/k)
XX = matrix(0, nrow=simul, ncol=n)
lam = matrix(0, nrow=simul, ncol=5)

for (i in 1:simul){
  for (j in 1:run){
    XX[i,((k*j-29):(k*j))] = sample(X[(k*j-29):(k*j)], replace=TRUE)
  }
  XX[i, ((k*j+1):n)] = sample(X[(k*j+1):n], replace=TRUE)

  model <- depmix(XX[i,] ~1, nstates=5, data=data.frame(XX[i,]), family=poisson())
  fm <- fit(model)
  u <-as.vector(getpars(fm))
  para.mle = c(u[6:30], exp(u[31]), exp(u[32]), exp(u[33]), exp(u[34]), exp(u[35]))
  
  mtrans = matrix(para.mle[1:25], byrow=TRUE, nrow=5)
  lam[i,] = para.mle[26:30]  
  lam[i,] = sort(lam[i,])
}

lam = read.csv("lambda.csv")
