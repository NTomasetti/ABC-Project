library(reshape)
library(plyr)
library(moments)
library(ggplot2)
library(tidyr)
library(dplyr)
library(randomForest)
library(caret)
library(matrixStats)
library(abcrf)

T <- 500
n <- 5000
#Simulate g-and-k random variables: Four parameters: a,b,g,k
param <- 4
a <- 3 #mean
b <- 1 #scale
g <- 2 #skewness
k <- .5 # kurtosis

# Step-1: Draw standard normal random variables
x <- rnorm(T,mean=0, sd=1)
# Step-2: Plug these into the following formula

y_var <- a+b*(1+.8*((1-exp(-g*x))/(1+exp(-g*x))))*((1+x^2)^k)*x

# Different choices of Summaries

# Simple moments: Should be very poor for many choices of g and k
s1=function(x, qr) 1/T*sum(x)
s2=function(x, qr) (1/T*sum((x-s1(x))^2))^0.5
s3=function(x, qr) 1/T*sum((x-s1(x))^3)/((1/T*sum((x-s1(x))^2))^1.5)
s4=function(x, qr) (1/T*sum((x-s1(x))^4)/((1/T*sum((x-s1(x))^2))^2))^0.5

# Robust moments: Should be very poor for many choices of g and k
p <- c(.25,.5,.75)
qr = as.numeric(quantile(x,p))
qr_stupid = as.numeric(quantile(x,c(.08,.36,.6,.875)))

s1r = function(y, x) y[4]
s2r = function(y, x) y[6]-y[2]
s3r = function(y, x) (y[6]+y[2]-2*y[4])/(y[6]-y[2])
s4r = function(y, x) (y[7]-y[5]+y[3]-y[1])/(y[6]-y[2])

# Priors for simulating the psuedo-data. 
theta <- matrix(0, nrow=n, ncol=param)
theta[,1] <- runif(n, 0, 10)
theta1 = theta[,1]
theta[,2] <- runif(n, 0, 10)
theta2 = theta[,2]
theta[,3] <- runif(n, 0, 10)
theta3 = theta[,3]
theta[,4] <- runif(n, 0, 10)
theta4= theta[,4]

# Generating simulated data now

u <- matrix( rnorm(T*n,mean=0,sd=1), T, n) 
# Step-2: Plug these into the following formula
z_var <- theta[,1]+theta[,1]*(1+.8*((1-exp(-theta[,3]*u))/(1+exp(-theta[,3]*u))))*((1+u^2)^theta[,4])*u


##
vars <- matrix(0,nrow=n, ncol = 4)
vars[,1] <- apply(z_var,2,s1)
vars[,2] <- apply(z_var,2,s2)
vars[,3] <- apply(z_var,2,s3)
vars[,4] <- apply(z_var,2,s4)
vars[,1] <- colQuantiles(z_var,.5)
vars[,2] <- colQuantiles(z_var,.75)-colQuantiles(z_var,.25)
vars[,3] <- (colQuantiles(z_var,.75)+colQuantiles(z_var,.25)-2*colQuantiles(z_var,.5))/vars[,6]
vars[,4] <- (colQuantiles(z_var,.875)-colQuantiles(z_var,.6)+colQuantiles(z_var,.375)-colQuantiles(z_var,.125))/vars[,6]


df1 = data.frame(theta1,vars)
df2 = data.frame(theta2,vars)
df3 = data.frame(theta3,vars)
df4 = data.frame(theta4,vars)



require(randomForest)
require(varSelRF)
fit=randomForest(theta1~., data=df1,ntree=1000,keep.forest=TRUE,importance=TRUE)
importance(fit,type=1)
varImpPlot(fit)
partialPlot(fit,df1,x.var=X3)
