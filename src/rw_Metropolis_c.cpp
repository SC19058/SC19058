#include <Rcpp.h>
using namespace Rcpp;
//' @title Random walk Metropolis using Rcpp
//' @description Implement the random walk version of the Metropolis using Rcpp
//' @param sigma the variance
//' @param x0 the intial location
//' @param N the total steps
//' @return a list conclude x and the number of accpetance times \code{accept}
//' @export
// [[Rcpp::export]]
List rw_Metropolis_c(double x0, double sigma, int N) {
  NumericVector x(N);
  as<DoubleVector>(x)[0] = x0;
  NumericVector u(N);
  u = as<DoubleVector>(runif(N));
  List out(2);
  int accept = 1;
  for(int i=1;i<N;i++){
    double y = as<double>(rnorm(1,x[i-1],sigma));
    if(u[i] <= exp(abs(x[i-1])-abs(y))){
        x[i] = y;accept = accept + 1;
    }
    else{
        x[i] = x[i-1];
    }
  }  
  out[0] = x;
  out[1] = accept;
  return(out);
}
