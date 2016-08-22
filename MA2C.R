##MA(2) Model
library(Rcpp)
library(ggplot2)
library(plyr)
library(dplyr)
library(RcppArmadillo)
library(microbenchmark)
library(reshape)
sourceCpp("MA2C.cpp")

#ABC details
n <- 50000
percent <- 1

##Data
T <- 250
##True thetas are 0.6 and 0.2
y <- genMA2(T, c(0.6, 0.2))  
##Dimension of theta
p <- 2



##Return indeces of N smallest values in a vector
index.bottom.N = function(x, N){
  o = order(x, na.last=FALSE, decreasing=TRUE)
  o.length = length(o)
  o[((o.length-N+1):o.length)]
}

##ABC inputs: number of reps, data, percent of reps kept, starting/finishing value for j, 
ABC <- function(n, y, percent, R, M){ #j is the number of sample statistics to use at once
  accept <- list() #storing results, different dimensions as j changes so requires a list
  theta <- gentheta(n) #generates from U(0,1) priors, accepts only stationary draws
  z <- apply(theta, 1, genMA2, n = T) #One MA2 model simulated per theta pair
  ABCvec <- c(s0=s0, s1=s1, s2=s2, s3=s3) #summary statistics used, name=function pairs
  abcss <- each(ABCvec) #combined into one function
  ABCnames <- names(ABCvec) #for ease of labels later
  for(j in R:M){ #loop through different lengths of eta
     yss <- combn(abcss(y), j)   #calculate y sample statistics
     names <- combn(ABCnames, j) #names of each statistic in yss
     if(j == 1){    ##assign column names to yss that match the sample statistics used
       colnames(yss) <- names
     } else if(j == length(ABCvec)){
       colnames(yss) <- apply(names, 2, paste, collapse=",")
     } else {
       colnames(yss) <- apply(names[1:j,], 2, paste, collapse=",")
     }
     distance <- matrix(0, ncol=n, nrow=ncol(yss)) #store distances between yss, zss
     for(k in 1:n){
       zss <- combn(abcss(z[,k]), j) #for each simulation z(k), calculate sample statistics
       jointss <- rbind(yss, zss) #merge together for apply to work properly
       distance[,k] <- apply(jointss, 2, eucdist) #euclidean distance between zss/yss
     }
     lowdist <- apply(distance, 1, index.bottom.N, n*percent/100) #select % lowest
     accepting <- array(0, dim=c(n*percent/100, ncol(yss), p)) #to store accepted values
     for(l in 1:ncol(yss)){
       accepting[,l,] <- theta[lowdist[,l],] #accept theta that correspond to lowest distances
     }
     colnames(accepting) <- colnames(yss) #column names of the statistics used
     accept[[(j-R+1)]] <- accepting #store different dimension arrays on list
  }
  accept.long <- melt.list(accept)[,-1] #convert list to dataframe for ease of plotting
  colnames(accept.long) <- c("stat", "theta", "value", "length")
  accept.long <- mutate(accept.long, theta = paste0("Theta ", theta),
                        length = length + R -1)
  return(accept.long) #output is four columns - summary statistics used, which theta, value of theta, length of eta
}
##Running ABC and plotting##

results <- ABC(n, y, percent, p, 4)
##ABC with 50000 simulations, matching against data vector y, keeping lowest 1% distance
##Starting with eta of length equal theta, ending with eta with length 4 (contains all s.stats)

##Example plots. First: All summary statistics, Second: Only shows eta containing s1 and s2
ggplot(data=results) + geom_density(aes(x=value, colour=stat)) + facet_grid(.~theta)
ggplot(data=subset(results, grepl("s1", results$stat) & grepl("s2", results$stat))) + geom_density(aes(x=value, colour=stat)) + facet_wrap(~theta)


