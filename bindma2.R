##Binding Funciton##
s0=function(x) 1/T*sum(x^2)
s1=function(x) 1/T*sum(x[1:(T-1)]*x[2:T])
s2=function(x) 1/T*sum(x[1:(T-2)]*x[3:T])
s3=function(x) 1/T*sum(x[1:(T-3)]*x[4:T])

bind.ma2 <- function(T, theta1, theta2, eps){
  z <- rep(0, T)
  z[1] <- eps[1]
  z[2] <- eps[2]+theta1*eps[1]
  for(t in 3:T){
    z[t] <- eps[t] + theta1*eps[(t-1)]+theta2*eps[(t-2)]
  }
  e0 <- s0(z)
  e1 <- s1(z)
  e2 <- s2(z)
  e3 <- s3(z)
  return(c(e0, e1, e2, e3))
}