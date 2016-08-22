library(reshape)
library(plyr)
library(moments)
library(ggplot2)
library(tidyr)
library(dplyr)
library(MASS)

T <- 1000
n <- 15000
percent <- .5
R <- 1 #starting length
M <- 1 #ending length

#Simulate g-and-k random variables: Four parameters: a,b,g,k
a <- 3 #mean
b <- 1 #scale
g <- 2 #skewness
k <- .5 # kurtosis
#dimension of theta
p <- 4

# Step-1: Draw standard normal random variables
xy <- rnorm(T,mean=0, sd=1)
qr = (quantile(y,c(.08,.25,.36,.5,.6,.75,.875)))
# Step-2: Plug these into the following formula
y <- a+b*(1+.8*((1-exp(-g*xy))/(1+exp(-g*xy))))*((1+xy^2)^k)*xy


# Different choices of Summaries

# Simple moments: Should be very poor for many choices of g and k
s1=function(x, qr) 1/T*sum(x)
s2=function(x, qr) (1/T*sum((x-s1(x))^2))
s3=function(x, qr) 1/T*sum((x-s1(x))^3)/((1/T*sum((x-s1(x))^2))^1.5)
s4=function(x, qr) (1/T*sum((x-s1(x))^4)/((1/T*sum((x-s1(x))^2))^2))^0.5

#pr <- 

#qr_stupid = as.numeric(quantile(x,c(.08,.36,.6,.875)))

s1r = function(x, qr) qr[4]
s2r = function(x, qr) qr[6]-qr[2]
s3r = function(x, qr) (qr[6]+qr[2]-2*qr[4])/(qr[6]-qr[2])
s4r = function(x, qr) (qr[7]-qr[5]+qr[3]-qr[1])/(qr[6]-qr[2])


#Add RF functions to ABCvec manually, change names if you want
ABCvec <- c(s1=s1r)#, s2=s2r, s3=s3r, s4=s4r)#, s1r=s1r, s2r=s2r, s3r=s3r, s4r=s4r) c(s1r=s1r, s2r=s2r, s3r=s3r, s4r=s4r)

##Setting up functions for ABC##

##euclidean distance between first and second half of vector, 
##had to rbind summ. stat matrices and use apply to calculate each column separately in ABC function
##replace last line with other choices of d(,) if required
Wdist <- function(x){
  le <- length(x)
  y <- as.numeric(x[1:(le/3)])
  z <- as.numeric(x[(le/3+1):(2*le/3)])
  var <- as.numeric(x[((2*le/3)+1):le])
  W <- diag(1/sqrt( var ))
  ( t((y-z)) %*%(y-z))
}

##Return N smallest values
index.bottom.N = function(x, N){
  o = order(x, na.last=FALSE, decreasing=TRUE)
  o.length = length(o)
  o[((o.length-N+1):o.length)]
}

#storing results, different dimensions as j changes so requires a list
accept <- list()
theta <- matrix(0, ncol=n, nrow=p)
theta[1,] <- runif(n, 0, 10)
#theta[2,] <- runif(n, 0, 10)
theta[2,] <- rep(1,n)
#theta[3,] <- runif(n, 0, 10)
#theta[3,] <- 2*rep(1,n)
theta[3,] <- 2*rep(1,n)
theta[4,] <- .5*rep(1,n)
#theta[4,] <- runif(n, 0, 10)
z <- matrix(0, ncol=n, nrow=T)
x <- matrix(rnorm(T*n), ncol=n)
  #loop through simulation of data
for(k in 1:n){
  z[,k] <- theta[1,k]+theta[2,k]*(1+.8*((1-exp(-theta[3,k]*x[,k]))/(1+exp(-theta[3,k]*x[,k]))))*((1+x[,k]^2)^theta[4,k])*x[,k]
}
for(j in R:M){ #loop through different lengths of eta
  #calculate y sample statistics
  fn <- combn(ABCvec, j)
  names <- combn(names(ABCvec), j)
  keep <- apply(names, 2, function(x) any(grepl("s1", x)) & any(grepl("s2", x)) & any(grepl("s3", x)) & any(grepl("s4", x)))
  b <- fn[,keep]
  names <- names[,keep]
  qr = as.numeric(quantile(y,c(.08,.25,.36,.5,.6,.75,.875)))
  if(j ==length(ABCvec)) {
    ssfn <- each(b)
    yss <- as.matrix(ssfn(y, qr))
    colnames(yss) <- "s1,s2,s3,s4,s1r,s2r,s3r,s4r"
    distance <- rep(0, n)
    zss <- matrix(0, ncol=n, nrow=8)
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
    accepting <- array(0, dim=c(n*percent/100, 2, p)) #to store accepted values
    accepting[,1,] <- t(theta[,lowdist]) #accept theta that correspond to lowest distances
    colnames(accepting) <- c("s1,s2,s3,s4,s1r,s2r,s3r,s4r", NA) #column names of the statistics used
    accept[[(j-R+1)]] <- accepting #store different dimension arrays on list
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
    accepting <- array(0, dim=c(n*percent/100, ncol(yss), p)) #to store accepted values
    for(l in 1:ncol(yss)){
      accepting[,l,] <- t(theta[,lowdist[,l]]) #accept theta that correspond to lowest distances
    }
    colnames(accepting) <- colnames(yss) #column names of the statistics used
    accept[[(j-R+1)]] <- accepting #store different dimension arrays on list
  }
}
accept.long <- melt.list(accept)[,-1] #convert list to dataframe for ease of plotting
colnames(accept.long) <- c("stat", "theta", "value", "length")
accept.long <- subset(accept.long, !is.na(stat))
accept.long <- mutate(accept.long, theta = paste0("Theta ", theta),
                        length = length + R -1)
#output is four columns - summary statistics used, which theta, value of theta, length of eta


##Example plots. First: All summary statistics, Second: Only shows eta of length 2, 
##Third: Only shows eta containing s2, Fourth: Only shows eta containing s3 and s2

ggplot(data=accept.long) + geom_density(aes(x=value, colour=stat)) + ylim(0,1.5) + facet_wrap(~theta, ncol=2)
#ggplot(data=subset(results, length==4)) + geom_density(aes(x=value, colour=stat)) + facet_wrap(~theta, ncol=2)
#ggplot(data=subset(results, grepl("s2", results$stat))) + geom_density(aes(x=value, colour=stat)) + facet_wrap(~theta, ncol=2)
#ggplot(data=subset(accept.long, grepl("s1", accept.long$stat) & grepl("s2", accept.long$stat))) + geom_density(aes(x=value, colour=stat)) + facet_wrap(~theta, ncol=2)
