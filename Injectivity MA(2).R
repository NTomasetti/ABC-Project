source("bindma2.R")

##Step One - Initialise theta grid##

T <- 15000

d_t <- 2
d_e <- 4

theta1 <- seq(-0.95, .95, .5)
theta2 <- seq(-0.95, .95, .5)
theta <- expand.grid(theta1, theta2)

##Step Two - Generate errors##
eps <- rnorm(T)

##Step Three - Calculate numerical partial derivatives##
h = 0.0000001

score <- array(0, dim=c(d_e, d_t, nrow(theta)))
for(i in 1:nrow(theta)){
  score1 <- rep(0, 4)
  score2 <- rep(0, 4)
  grid1 <- c(theta[i,1]+0.5*h, theta[i,1]-0.5*h)
  score1 <- (bind.ma2(T, grid1[1], theta[i,2], eps) - bind.ma2(T, grid1[2], theta[i,2], eps))/h
  grid2 <- c(theta[i,2]+0.5*h, theta[i,2]-0.5*h)
  score2 <- (bind.ma2(T, theta[i,1], grid2[1], eps) - bind.ma2(T, theta[i,1], grid2[2], eps))/h
  score[,,i] <- cbind(score1, score2)
}

##Step Four - Take every combination of 2x2 matrices##
b <- combn(d_e, d_t)
mats2x2 <- array(0, dim=c(d_t, d_t, nrow(theta), ncol(b)))
for(i in 1:ncol(b)){
  mats2x2[,,,i] <- score[b[,i],,]   ##First two dims - 2x2 matrix, third dim - theta grid, fourth dim - b combination
}

##Step Five - Check eigenvalues##
posdef <- apply(mats2x2, 3:4, function(x) all(Re(eigen(x)$value)>0 & all(Im(eigen(x)$value)==0)))
b[,apply(posdef, 2, function(x) all(x == TRUE))]
