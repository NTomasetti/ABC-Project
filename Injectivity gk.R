source("bindgk.R")


##Step One - Initialise theta grid##

T <- 25000

d_t <- 4
d_e <- 4

theta1 <- seq(-1, 4, 1)
theta2 <- seq(0.5, 4, 1)
theta3 <- seq(-1, 1, 0.5)
theta4 <- seq(0, 1, 0.5)
theta <- expand.grid(theta1, theta2, theta3, theta4)

##Step Two - Generate errors##
eps <- rnorm(T)

##Step Three - Calculate numerical partial derivatives##
h = 0.0000001

score <- array(0, dim=c(d_e, d_t, nrow(theta)))

par.deriv <- function(theta, T, eps){
  score1 <- rep(0, 4)
  score2 <- rep(0, 4)
  score3 <- rep(0, 4)
  score4 <- rep(0, 4)
  grid1 <- c(theta[1]+0.5*h, theta[1]-0.5*h)
  score1 <- (bind.gk(T, grid1[1], theta[2], theta[3], theta[4], eps) - bind.gk(T, grid1[2], theta[2], theta[3], theta[4],eps))/h
  grid2 <- c(theta[2]+0.5*h, theta[2]-0.5*h)
  score2 <- (bind.gk(T, theta[1], grid2[1],theta[3],theta[4],eps) - bind.gk(T, theta[1], grid2[2],theta[3],theta[4],eps))/h
  grid3 <- c(theta[3]+0.5*h, theta[3]-0.5*h)
  score3 <- (bind.gk(T, theta[1], theta[2], grid3[1], theta[4],eps) - bind.gk(T, theta[1], theta[2], grid3[2], theta[4],eps))/h
  grid4 <- c(theta[4]+0.5*h, theta[4]-0.5*h)
  score4 <- (bind.gk(T, theta[1], theta[2],theta[3],grid4[1],eps) - bind.gk(T, theta[1], theta[2],theta[3],grid4[2],eps))/h
  out <- matrix(c(score1, score2, score3, score4), ncol=4, byrow=TRUE)
}

score <- array(apply(theta, 1, par.deriv, T, eps), dim=c(4,4,nrow(theta)))

mat1x1 <- score[1,1,]
mat3x3 <- score[2:4,2:4,]

##Step Four - Check eigenvalues##
posdef1 <- sapply(mat1x1, function(x) all(Re(eigen(x)$value)>0 & all(Im(eigen(x)$value)==0)))
posdef3 <- apply(mat3x3, 3, function(x) all(Re(eigen(x)$value)>0 & all(Im(eigen(x)$value)==0)))
posdef <- apply(score, 3, function(x) all(Re(eigen(x)$value)>0 & all(Im(eigen(x)$value)==0)))
all(posdef1 == TRUE)
all(posdef3 == TRUE)
all(posdef == posdef3)

#variable combinations where it is pos.def
#theta[posdef,]
