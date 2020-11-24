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
print(mean(X))
print(sd(X))
print(min(X))
print(max(X))
print(mode(X))
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
mode <- getmode(X)
print(mode)



set.seed(90210)
# 1-state model
model <- depmix(X ~1, nstates=1, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[2], exp(u[3]))
lams = para.mle[2]
print(lams)



# 2-states
set.seed(90210)
model <- depmix(X ~1, nstates=2, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
if (u[7] <= u[8]) { para.mle = c(u[3:6], exp(u[7]), exp(u[8]))
} else { para.mle = c(u[6:3], exp(u[8]), exp(u[7]))
}
mtrans = matrix(para.mle[1:4], byrow=TRUE, nrow=2)
lams = para.mle[5:6]
print(lams)
print(mtrans)

pi1 = mtrans[2,1]/(2 - mtrans[1,1] - mtrans[2,2])
pi2 = 1 - pi1
plot(ts(posterior(fm)[,3], start=1900), ylab=expression(hat(pi)[~t]*' (2 | dat
a)'))
abline(h=.5, lty=2)

hist(data$`m>4`, breaks=30, prob=TRUE, main="")
xvals = seq(1,45)
u1 = pi1*dpois(xvals, lams[1])
u2 = pi2*dpois(xvals, lams[2])
lines(xvals, u1, col=3)
lines(xvals, u2, col=4)

#data and states
plot(earth, main="", xlab="date", ylab='frequency', type='h', col='#808080')
text(earth, col=posterior(fm)[,1], labels=posterior(fm)[,1], cex=.9)

# residuals
resid <- earth - lams[posterior(fm)[,1]]
plot(resid, type="h")
acf(resid)
pacf(resid)

#three states
set.seed(90210)
model <- depmix(X ~1, nstates=3, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[4:12], exp(u[13]), exp(u[14]), exp(u[15]))

mtrans = matrix(para.mle[1:9], byrow=TRUE, nrow=3)
lams = rev(para.mle[10:12])
print(lams)
print(mtrans)

pi1 = mtrans[2,1]/(2 - mtrans[1,1] - mtrans[2,2])
pi2 = 1 - pi1
plot(ts(posterior(fm)[,3], start=1900), ylab=expression(hat(pi)[~t]*' (2 | data)'))
abline(h=.5, lty=2)

hist(data$`m>4`, breaks=30, prob=TRUE, main="")
xvals = seq(1,45)
u1 = pi1*dpois(xvals, lams[1])
u2 = pi2*dpois(xvals, lams[2])
lines(xvals, u1, col=3)
lines(xvals, u2, col=4)

S=posterior(fm)
for (i in 1:572){
if (S[i,1]==1){ S[i,1]=3}
else if (S[i,1]==3){S[i,1]=2}
else {S[i,1]=1}
}

#data and states
plot(earth, main="", xlab="date", ylab='frequency', type='h', col='#808080')
text(earth, col=S[,1], labels=S[,1], cex=.9)

# residuals
resid <- earth - lams[posterior(fm)[,1]]
plot(resid, type="h")
acf(resid)
pacf(resid)

# 4-states
set.seed(90210)
model <- depmix(X ~1, nstates=4, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[5:20], exp(u[21]), exp(u[22]), exp(u[23]), exp(u[24]))

mtrans = matrix(para.mle[1:16], byrow=TRUE, nrow=4)
lams = para.mle[17:20]
print(lams)
print(mtrans)

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

# 6-states
set.seed(90210)
model <- depmix(X ~1, nstates=6, data=data.frame(X), family=poisson())
fm <- fit(model)
u <-as.vector(getpars(fm))
para.mle = c(u[7:42], exp(u[43]), exp(u[44]), exp(u[45]), exp(u[46]), exp(u[47]), exp(u[48]))

mtrans = matrix(para.mle[1:36], byrow=TRUE, nrow=6)
lams = para.mle[37:42]
print(lams)
print(mtrans)

# 7-states
set.seed(90210)
model <- depmix(X ~1, nstates=7, data=data.frame(X), family=poisson())
fm <- fit(model)

