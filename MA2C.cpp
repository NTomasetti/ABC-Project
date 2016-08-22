// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <functional>   
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
vec genMA2(int n, vec theta) { //Generates the errors of a length n MA2 model
  vec eps(n);        //then loops through the formula to create it
  eps = randn<vec>(n);         //Theta parameter should be length 2
  vec out(n);
  out[0] = eps[0];
  out[1] = theta[0]*eps[0] + eps[1];
  for(int i = 2; i < n; ++i){
    out[i] = theta[0]*eps[i-1] + theta[1]*eps[i-2] + eps[i];
  }
  return out;
}

// [[Rcpp::export]]
mat gentheta(int n) {  //n is the number of theta pairs
  mat out(n,2);
  for(int i = 0; i < n; ++i) {
    bool check = FALSE;  //Stationary flag
    while(!check) {
      out(i,0) = randu<vec>(1)[1]; //Generates from uniform priors
      out(i,1) = randu<vec>(1)[1];
      if(out(i,0) + out(i,1) < 1 & out(i,0)-out(i,1) < 1) {
        check = TRUE;  //Checking the stationary conditions of an MA2 model
      }                //We only accept stationary values for theta
    }
  }
  return(out);
}

// [[Rcpp::export]]
double s0(vec x) {  //We want 1/T sum(x^2)
  int T = x.size();   
  double out =0;
  for(int i = 0; i < T; ++i) {
    out += x[i]*x[i];
  }
  return out/T;
}

// [[Rcpp::export]]
double s1(vec x) { //We want 1/T sum(x[1:T-1]*x[2:T])
  int T = x.size();  
  double out;
  for(int i = 0; i < T-1; ++i) {
    out += x[i+1]*x[i];
  }
  return out/T;
}

// [[Rcpp::export]]
double s2(vec x) { //We want 1/T sum(x[1:T-2]*x[3:T])
  int T = x.size(); 
  double out;
  for(int i = 0; i < T-2; ++i) {
    out += x[i+2]*x[i];
  }
  return out/T;
}

// [[Rcpp::export]]
double s3(vec x) { //We want 1/T sum(x[1:T-3]*x[4:T])
  int T = x.size();
  double out;
  for(int i = 0; i < T-3; ++i) {
    out += x[i+3]*x[i];
  }
  return out/T;
}

// [[Rcpp::export]]
double eucdist(vec x) { //Euclidean distance
  int n = x.size();                //Because of the way R's apply works, we need zss and yss on the same input vector. So we use the
  double out = 0;                  //first n/2 elements as yss, second n/2 as zss then take the euclidean distance sqrt(sum((x-y)^2))
  for(int i=0; i < n/2; ++i) {     //all handled by the loop.
    out += (x[i]-x[n/2+i])*(x[i]-x[n/2+i]); 
  }
  return sqrt(out);
}




