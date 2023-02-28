#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// ==============================================================================
//' @title Generate test statistic from bootstrap samples
//' 
//' @description This function simulates the DGP under the null hypothesis and computes the test statistic for each sample
//' 
//' @return vector of simulated test statistics
//' 
//' @export
// [[Rcpp::export]]
arma::vec bootNullDist_ex(double param_null, double std, int n, int B){
  arma::vec tau_B_tmp(B, arma::fill::zeros);
  arma::vec y_t_sim;
  double gamma_hat_sim;
  double std_sim;
  double tau_sim;
  for (int xb = 0; xb < B; xb++){
    y_t_sim = param_null + arma::randn<arma::vec>(n)*std;
    gamma_hat_sim = mean(y_t_sim);
    std_sim = sqrt(sum(pow(y_t_sim-gamma_hat_sim,2))/(n-1));
    tau_sim = (gamma_hat_sim - param_null)/(std_sim/sqrt(n));
    tau_B_tmp(xb) = tau_sim;
  }
  return(tau_B_tmp);
}


