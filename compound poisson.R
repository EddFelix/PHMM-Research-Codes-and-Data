setwd("~/Documents/PHMM")

library(anytime)
library(Hmisc)
library(pastecs)
library(readxl)

#data <- read_excel("miso_minn_24hr_avg.xlsx", sheet=2)
#NN=data$`pois 3sd`
#NN=na.omit(NN)


NN=rpois(1000, 20)
N=sum(NN)

#Y=na.omit(data$`JUMP mean+3sd`)

Y=rexp(1000, 1/6)

sumY = sum(Y)
m = length(Y)
tol=1e-6

h <-  function(y, eta){
  sum = eta
  l=1
  k = 2
  tol=1e-6
  while (l>tol) {
    l = eta^k*y^(k-1) / (factorial(k)*factorial(k-1))
    sum = sum + l
    k = k + 1
  }
  return(sum)
}

dh <-  function(y, eta){
  sum = 1
  l=1
  k = 2
  tol=1e-6
  while (l>tol) {
    l = eta^(k-1)*y^(k-1) / (factorial(k-1)*factorial(k-1))
    sum = sum + l
    k = k + 1
  }
  return(sum)
}

ddh <- function(y, eta){
  ddh = y/eta*h(y, eta)
  return(ddh)
}

hess <- function(y, eta){
  hess = ddh(y, eta)/h(y, eta)- (dh(y, eta)/h(y, eta))^2
  return(hess)
}

G <- function(eta, Y){
  logh=0
  for (i in 1:m){
    logh = logh + log(h(Y[i], eta))
  }
  G = -2*sqrt(eta*N*sumY) + logh
  return(G)
}

dG <- function(eta, Y){
  dlogh = 0
  for (i in 1:m){
    dlogh = dlogh + dh(Y[i], eta) / h(Y[i], eta)
  }
  dG = -eta^(-1/2)*sqrt(N*sumY) + dlogh
  return(dG)
}

ddG <- function(eta, Y){
  sumhess=0
  for (i in 1:m){
    sumhess = sumhess + hess(Y[i], eta)
  }
  ddG = 1/2*eta^(-3/2)*sqrt(N*sumY)+sumhess
  return(ddG)
}

eta0=.01
toleta=1e-6
maxiter=1000
eta1 = eta0-dG(eta0,Y)/ddG(eta0,Y)
diff=abs(eta1-eta0)
iter=1

while ((iter<maxiter)&(diff>toleta)){
  eta0=eta1
  eta1 = eta0-dG(eta0,Y)/ddG(eta0,Y)
  diff=abs(eta1-eta0)
  iter=iter+1
}

eta=eta1
lambda=sqrt(eta*sumY/N)
theta=eta/lambda
theta=1/theta
lambda=1/lambda

print(eta)
print(lambda)
print(theta)