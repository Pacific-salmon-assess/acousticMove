#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH)]]
#include <Rcpp.h>
using namespace Rcpp;

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

// Sparse A
// [[Rcpp::export]]
arma::vec expAv_gr(const arma::sp_mat& A, const arma::vec& v, double tol, int renorm_freq, bool forward) {
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


expAv_gr <- function(Q, v, dQ, t = 1, rescale_freq = 50){
  N <- nrow(Q)
  g <- max(-diag(Q))
  gt <- g*t
  S <- g*t*diag(N) + Q*t
  D1 <- Q
  D2 <- Q
  D1@x <- t*dQ[,1]
  D2@x <- t*dQ[,2]
  w <- v
  expQv <- v
  delta <- Matrix(0, nrow = N, ncol = 2)
  dexpQv <- Matrix(0, nrow = N, ncol = 2)
  Nmax <- qpois(1e-8, gt, lower.tail = FALSE)
  log_scale <- 0

  for( n in 1:Nmax ){
    delta[,1] <- (D1 %*% w + S %*% delta[,1])/n
    delta[,2] <- (D2 %*% w + S %*% delta[,2])/n
    dexpQv[,1] <- dexpQv[,1] + delta[,1]
    dexpQv[,2] <- dexpQv[,2] + delta[,2]
    w <- S %*% w/n
    expQv <- expQv + w    
    if(!(n %% rescale_freq)){
      s <- sum(abs(w))
      w <- w/s
      delta <- delta/s
      expQv <- expQv/s
      dexpQv <- dexpQv/s
      log_scale <- log_scale + log(s)
    }    
  }
  return(list(expQv*exp(-gt + log_scale), dexpQv*exp(-gt + log_scale)))
}
