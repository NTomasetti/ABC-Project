library(R2Cuba)

true1 <- 0.6   #setup
true2 <- 0.2
T <- 250

y <- arima.sim(list(ma=c(true1, true2)), T)

theta1 <- seq(0.19, .99, 0.02)
theta2 <- seq(-.21, .59, 0.02)
theta <- expand.grid(theta1, theta2)

#kalman filter
kalman <- function(theta){
  eps <- matrix(0, ncol=T, nrow=3) #eps(t|t)
  epsl <- matrix(0, ncol=T, nrow=3) #eps(t|t-1)
  P <- array(0, dim=c(3,3,T)) #P(t|t)
  Pl <- array(0, dim=c(3,3,T)) #P(t|t-1)
  R <- matrix(0, ncol=3, nrow=3) #formerly F
  R[2,1] <- 1
  R[3,2] <- 1
  H <- c(1, theta[1], theta[2])
  Q <- matrix(c(1, rep(0,8)), 3,3)
  v <- rep(0,T) #y - mean(y)
  F <- rep(0,T) #var(y)
  M <- matrix(0, ncol=T, nrow=3)
  
  
  #time 1 - #epsl[,1] starts at 0, Pl[,,1] starts at Identity(3)
  Pl[,,1] <- diag(rep(1,3))
  v[1] <- y[1] - t(H)%*%epsl[,1]
  F[1] <- t(H)%*%Pl[,,1]%*%H
  M[,1] <- Pl[,,1]%*%H
  eps[,1] <- epsl[,1] + M[,1]%*%solve(F[1])%*%v[1]
  P[,,1] <- Pl[,,1] - M[,1]%*%solve(F[1])%*%t(M[,1])

  #update loop
  for(t in 2:T){
    epsl[,t] <-  R%*%eps[,(t-1)]
    Pl[,,t] <- R%*%P[,,(t-1)]%*%t(R)+Q
    v[t] <- y[t] - t(H)%*%epsl[,t]
    F[t] <- t(H)%*%Pl[,,t]%*%H
    M[,t] <- Pl[,,t]%*%H
    eps[,t] <- epsl[,t] + M[,t]%*%solve(F[t])%*%v[t]
    P[,,t] <- Pl[,,t] - M[,t]%*%solve(F[t])%*%t(M[,t])
  }
  prob <- rep(0,T)
  prob[1] <- dnorm(y[1])
  for(t in 2:T) {
    prob[t] <- dnorm(y[t], mean=t(H)%*%epsl[,t], sd=sqrt(F[t]))
  }
  lik <- prod(prob)
  lik
}

Denom <- cuhre(2, 1, kalman, lower=c(-0.99,-0.99), upper=c(0.99,0.99))[[5]] #Integrate over region Theta

lik <- matrix(apply(theta, 1, kalman), ncol=length(theta2))/(2*Denom)

pos.theta1 <- rep(0, length(theta1)) #Trapezoid integration
for(i in 1:length(theta1)){
  pos.theta1[i] <- sum(lik[i,2:length(theta1)], lag(lik)[i,2:length(theta1)])*((max(theta1)-min(theta1))/(2*(length(theta1)-1)))
}
#verify marginal density integrates to 1
sum(pos.theta1[2:length(pos.theta1)], lag(pos.theta1)[2:length(pos.theta1)])*((max(theta1)-min(theta1))/(2*(length(pos.theta1)-1)))
pos.theta2 <- rep(0, length(theta2))
for(i in 1:length(theta2)){
  pos.theta2[i] <- sum(lik[2:length(theta2),i], lag(lik)[2:length(theta2),i])*((max(theta2)-min(theta2))/(2*(length(theta2)-1)))
}

posterior <- data.frame(value = c(pos.theta1, pos.theta2), support = c(theta1, theta2), theta=c(rep("Theta 1", length(theta1)), rep("Theta 2", length(theta2))))
ggplot() + geom_line(data=posterior, aes(x=support, y=value)) + facet_wrap(~theta) #+ geom_density(data=results, aes(x=value, colour=stat))
#where results is the output from the ABC procedure