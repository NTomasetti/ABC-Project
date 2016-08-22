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
vec gengk(int T, vec theta) {
  double a = theta[0];
  double b = theta[1];
  double g = theta[2];
  double k = theta[3];
  vec out(T);
  vec eps(T);
  eps = randn<vec>(T);
  for(int i = 0; i < T; ++i) {
    double expo = exp(-g*eps[i]);
    out(i) = a + b*(1+0.8*((1-expo)/(1+expo)))*pow(1+eps[i]*eps[i],k)*eps[i];
  }
  return out;
}

// [[Rcpp::export]]
double s1(vec x, vec qr) {
  int T = x.size();
  double out = 0;
  for(int i = 0; i < T; ++i) {
    out += x[i]*x[i];
  }
  return out/T;
}

// [[Rcpp::export]]
double s2(vec x, vec qr) {
  int T = x.size();
  double out = 0;
  double mean = s1(x, qr);
  for(int i = 0; i < T; ++i) {
    out += pow((x[i]-mean), 2.0);
  }
  return sqrt(out/T);
}

// [[Rcpp::export]]
double s3(vec x, vec qr) {
  int T = x.size();
  double num = 0;
  double mean = s1(x, qr);
  for(int i = 0; i < T; ++i) {
    num += pow((x[i]-mean), 3.0);
  }
  double denom = pow(s2(x, qr), 3.0);
  return num/(T*denom);
}

// [[Rcpp::export]]
double s4(vec x, vec qr) {
  int T = x.size();
  double num = 0;
  double mean = s1(x, qr);
  for(int i = 0; i < T; ++i) {
    num += pow((x[i]-mean), 4.0);
  }
  double denom = pow(s2(x, qr), 4.0);
  return sqrt(num/(T*denom));
}

// [[Rcpp::export]]
double s1r(vec x, vec qr) {
  return qr[3];
}

// [[Rcpp::export]]
double s2r(vec x, vec qr) {
  return qr[5]-qr[1];
}

// [[Rcpp::export]]
double s3r(vec x, vec qr) {
  return (qr[5]+qr[1]-2*qr[3])/(qr[5]-qr[1]);
}

// [[Rcpp::export]]
double s4r(vec x, vec qr) {
  return (qr[6]-qr[4]+qr[2]-qr[0])/(qr[5]-qr[1]);
}

// [[Rcpp::export]]
double Wdist(vec x){
  int n = x.size();
  vec y(n/3);
  vec z(n/3);
  vec var(n/3);
  mat W(n/3, n/3);
  for(int i = 0; i < n/3; ++i) {
    y[i] = x[i];
    z[i] = x[n/3+i];
    var[i] = x[2*n/3+i];
  }
  for(int i = 0; i < n/3; ++i) {
      W(i,i) = 1/sqrt(var[i]);
  }
  double out = 0;
  vec diff = y - z;
  for(int i = 0; i < n/3; ++i) {
    for(int j = 0; j < n/3; ++j)  {
      out += diff[i]*diff[j]*W(i,j);
    }
  }
  return sqrt(out);
}

// [[Rcpp::export]]
mat gentheta(int n, int lower, int upper) {
  mat out(n, 4);
  for(int i = 0; i < n; ++i) {
    out.row(i) = trans(lower + (upper-lower)*randu<vec>(4));
  }
  return out;
}
