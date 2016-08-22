library(reshape)
library(plyr)
library(moments)
library(ggplot2)
library(tidyr)
library(dplyr)
library(MASS)

T <- 2000
n <- 5000
percent <- .5
R <- 3 #starting length
M <- 6 #ending length

#Simulate g-and-k random variables: Four parameters: a,b,g,k
a <- 3 #mean
b <- 1 #scale
g <- 2 #skewness
k <- 0.5 # kurtosis
#dimension of theta
p <- 4

# Step-1: Draw standard normal random variables
xy <- rnorm(T,mean=0, sd=1)
qr = as.numeric(quantile(xy,c(.08,.25,.36,.5,.6,.75,.875)))
# Step-2: Plug these into the following formula
y <- a+b*(1+.8*((1-exp(-g*xy))/(1+exp(-g*xy))))*((1+xy^2)^k)*xy


# Different choices of Summaries

# Simple moments: Should be very poor for many choices of g and k
s1=function(x, qr) 1/T*sum(x)
s2=function(x, qr) (1/T*sum((x-s1(x))^2))^0.5
s3=function(x, qr) 1/T*sum((x-s1(x))^3)/((1/T*sum((x-s1(x))^2))^1.5)
s4=function(x, qr) (1/T*sum((x-s1(x))^4)/((1/T*sum((x-s1(x))^2))^2))^0.5
s1r = function(x, qr) qr[4]
s2r = function(x, qr) qr[6]-qr[2]
s3r = function(x, qr) (qr[6]+qr[2]-2*qr[4])/(qr[6]-qr[2])
s4r = function(x, qr) (qr[7]-qr[5]+qr[3]-qr[1])/(qr[6]-qr[2])


#Add RF functions to ABCvec manually, change names if you want
ABCvecmean <- c(s1=s1, s1r=s1r)
ABCvecother <- c(s2=s2, s3=s3, s4=s4, s2r=s2r, s3r=s3r, s4r=s4r)

##Setting up functions for ABC##

##euclidean distance between first and second half of vector, 
##had to rbind summ. stat matrices and use apply to calculate each column separately in ABC function
##replace last line with other choices of d(,) if required
Wdist <- function(x){
  le <- length(x)
  y <- as.numeric(x[1:(le/3)])
  z <- as.numeric(x[(le/3+1):(2*le/3)])
  var <- as.numeric(x[((2*le/3)+1):le])
  if(length(var)==1)
    W <- 1/sqrt(var)
  else
    W <- diag(1/sqrt( var ))
  sqrt( t((y-z)) %*% W %*%(y-z))
}

##Return N smallest values
index.bottom.N = function(x, N){
  o = order(x, na.last=FALSE, decreasing=TRUE)
  o.length = length(o)
  o[((o.length-N+1):o.length)]
}

acceptother <- list()
theta <- matrix(0, ncol=n, nrow=p)
accepting <- matrix(0, nrow=n*percent/100, ncol=3)
theta[1,] <- runif(n, 0, 10)
theta[2,] <- runif(n, 0, 10)
theta[3,] <- runif(n, 0, 10)
theta[4,] <- runif(n, 0, 10)
z <- matrix(0, ncol=n, nrow=T)
x <- matrix(rnorm(T*n), ncol=n)
#loop through simulation of data
for(k in 1:n){
  z[,k] <- theta[1,k]+theta[2,k]*(1+.8*((1-exp(-theta[3,k]*x[,k]))/(1+exp(-theta[3,k]*x[,k]))))*((1+x[,k]^2)^theta[4,k])*x[,k]
}

#use s1, s1r on first parameter
qr = as.numeric(quantile(y,c(.08,.25,.36,.5,.6,.75,.875)))
ssfn <- each(ABCvecmean) #able to apply all statistics in one function 
yss <- ssfn(y, qr) #calculate sample statistics on y
distance <- matrix(0, ncol=n, nrow=3)
zss <- matrix(0, ncol=n, nrow=2)
for(k in 1:n){
  qr = as.numeric(quantile(z[,k],c(.08,.25,.36,.5,.6,.75,.875))) #calculate sample statistics on z
  zss[,k] <- ssfn(z[,k], qr)
}
varz <- apply(zss, 1, var) #variance of a sample statistic across draws of z
for(k in 1:n){
  distance[1,k] <- Wdist(c(yss[1], zss[1,k], varz[1])) #calculate distances
  distance[2,k] <- Wdist(c(yss[2], zss[2,k], varz[2]))
  distance[3,k] <- Wdist(c(yss, zss[,k], varz))
}
lowdist <- apply(distance, 1, index.bottom.N, n*percent/100) #choose lowest distances

for(l in 1:3){
  accepting[,l] <- (theta[1,lowdist[,l]]) #accept theta that correspond to lowest distances
}
colnames(accepting) <- c("s1", "s1r", "s1,s1r") #column names of the statistics used
mean.long <- melt(accepting)[,-1] #convert to long form for ease of plotting
colnames(mean.long) <- c("stat", "value")

#use other s.stats on other parameters
for(j in R:M){
  fn <- combn(ABCvecother, j)
  names <- combn(names(ABCvecother), j)
  keep <- apply(names, 2, function(x) any(grepl("s2", x)) & any(grepl("s3", x)) & any(grepl("s4", x)))
  b <- fn[,keep]
  names <- names[,keep]
  qr = as.numeric(quantile(y,c(.08,.25,.36,.5,.6,.75,.875)))
  if(j ==length(ABCvecother)) {
    ssfn <- each(b)
    yss <- as.matrix(ssfn(y, qr))
    colnames(yss) <- "s2,s3,s4,s2r,s3r,s4r"
    distance <- rep(0, n)
    zss <- matrix(0, ncol=n, nrow=6)
    for(k in 1:n){
      qr = as.numeric(quantile(z[,k],c(.08,.25,.36,.5,.6,.75,.875)))
      zss[,k] <- ssfn(z[,k], qr)
    }
    varz <- apply(zss, 1, var)
    for(k in 1:n){
      jointss <- c(yss, zss[,k], varz)
      distance[k] <- Wdist(jointss)
    }
    lowdist <- index.bottom.N(distance, n*percent/100)
    accepting <- array(0, dim=c(n*percent/100, 2, 3)) #to store accepted values
    accepting[,1,] <- t(theta[2:4,lowdist]) #accept theta that correspond to lowest distances
    colnames(accepting) <- c("s2,s3,s4,s2r,s3r,s4r", NA) #column names of the statistics used
    acceptother[[(j-R+1)]] <- accepting #store different dimension arrays on list
  } else {
    yss <- apply(b, 1:2, function(f) f[[1]](y, qr))
    if(j == 1){
      colnames(yss) <- names
    } else {
      colnames(yss) <- apply(names[1:j,], 2, paste, collapse=",")
    }
    distance <- matrix(0, ncol=n, nrow=ncol(yss)) #store distances to z
    zss <- array(0, dim=c(j, ncol(yss), n))
    for(k in 1:n){ #for each simulation k, calculate sample statistics
      qr = as.numeric(quantile(z[,k],c(.08,.25,.36,.5,.6,.75,.875)))
      zss[,,k] <- apply(b, 1:2, function(f) f[[1]](z[,k], qr))
    }
    varz <- apply(zss, 1:2, var)
    for(k in 1:n){
      jointss <- rbind(yss, zss[,,k], varz) #merge together for apply to work properly
      distance[,k] <- apply(jointss, 2, Wdist) #euclidean distance to y
    }
    lowdist <- apply(distance, 1, index.bottom.N, n*percent/100) #select % lowest
    accepting <- array(0, dim=c(n*percent/100, ncol(yss), 3)) #to store accepted values
    for(l in 1:ncol(yss)){
      accepting[,l,] <- t(theta[2:4,lowdist[,l]]) #accept theta that correspond to lowest distances
    }
    colnames(accepting) <- colnames(yss) #column names of the statistics used
    acceptother[[(j-R+1)]] <- accepting #store different dimension arrays on list
  }
}
other.long <- melt.list(acceptother)[,-1] #convert list to dataframe for ease of plotting
colnames(other.long) <- c("stat", "theta", "value", "length")
other.long <- subset(other.long, !is.na(stat))
other.long <- mutate(other.long, theta = paste0("Theta ", theta+1),
                      length = length + R -1)
#output is four columns - summary statistics used, which theta, value of theta, length of eta

ggplot(data=mean.long) + geom_density(aes(x=value, colour=stat))
ggplot(data=subset(other.long, length==3)) + geom_density(aes(x=value, colour=stat)) + facet_wrap(~theta) +ylim(0,2)
