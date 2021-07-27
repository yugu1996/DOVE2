#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <ctime>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec BS(double t, arma::vec knots) {
  // define a row arma::vector {B_1(t), ..., B_K(t)}, based on the knots 
  // the slope of the last piece is 0 if constantVE = True
  int npc = knots.n_elem+1;
  arma::rowvec B;
  
  B.set_size(npc);
  B(0) = t;
  for (int i=1; i<npc; ++i) {
    B(i) = (t>knots(i-1) ? t-knots(i-1) : 0);
  }
  B = B/30.5;
  return B;
}
