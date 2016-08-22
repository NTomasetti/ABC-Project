##Binding Funciton##
s1=function(x) 1/T*sum(x)
s2=function(x) (1/T*sum((x-s1(x))^2))^0.5
s3=function(x) 1/T*sum((x-s1(x))^3)/((1/T*sum((x-s1(x))^2))^1.5)
s4=function(x) (1/T*sum((x-s1(x))^4)/((1/T*sum((x-s1(x))^2))^2))^0.5

bind.gk <- function(T, theta1, theta2, theta3, theta4, eps){
  z <- theta1+theta2*(1+.8*((1-exp(-theta3*eps))/(1+exp(-theta3*eps))))*((1+eps^2)^theta4)*eps
  e1 <- s1(z)
  e2 <- s2(z)
  e3 <- s3(z)
  e4 <- s4(z)
  return(c(e1, e2, e3, e4))
}