#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace arma;

// Sparse A
// [[Rcpp::export]]
arma::vec expAv_cpp(const arma::sp_mat& A, const arma::vec& v, double tol, int renorm_freq, bool forward) {
  const double rho=-A.min();
  if (rho==0) {
    return v;
  }
  const int nr=A.n_rows, nc=A.n_cols;
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  unsigned int n;
  const arma::sp_mat P = A+rho*arma::speye(nr, nc);
    
  arma::vec term = v, ans = v;
  double log_scale = 0.0;
  double s = 0.0;
  // Rcout << "rho: " << rho << std::endl;
  // Rcout << "N: " << niters << std::endl;

  n=1;
  while (n<=niters) {
    term/=n;
    if(!forward) { term = P*term;
    } else { term = term*P; }
    
    ans = ans + term;
    
    if(n % renorm_freq == 0){
      s = sum(abs(ans));
      term/=s;
      ans/=s;
      log_scale += log(s);
    }
    n++;
  }
  ans = exp(-rho + log_scale) * ans;
  return ans;
}

// arma::vec expAv_gr(const arma::sp_mat& A, const arma::vec& v, double tol, int renorm_freq, bool forward) {
  // const double rho=-A.min();
  // if (rho==0) {
    // return v;
  // }
  // const int nr=A.n_rows, nc=A.n_cols;
  // const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  // unsigned int n;
  // const arma::sp_mat P = A+rho*arma::speye(nr, nc);
    
  // arma::vec term = v, ans = v;
  // double log_scale = 0.0;
  // double s = 0.0;
  //// Rcout << "rho: " << rho << std::endl;
  //// Rcout << "N: " << niters << std::endl;

  // n=1;
  // while (n<=niters) {
    // term/=n;
    // if(!forward) { term = P*term;
    // } else { term = term*P; }
    
    // ans = ans + term;
    
    // if(n % renorm_freq == 0){
      // s = sum(abs(ans));
      // term/=s;
      // ans/=s;
      // log_scale += log(s);
    // }
    // n++;
  // }
  // ans = exp(-rho + log_scale) * ans;
  // return ans;
// }
