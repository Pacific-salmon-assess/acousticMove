#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace arma;

//' Compute expAv with the iters and uniform matrix pre-computed.
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v Vector (column vector) of to do exp(A)*v
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @details Implements Sherlock 2021 uniformization method.
// [[Rcpp::export]]
arma::vec expAv_light(const arma::sp_mat& A, const arma::vec& v, const double rho, const unsigned int niters, int renorm_freq){

  double log_scale = 0.0;
  double s = 0.0;  
  arma::vec term = v, ans = v;
  unsigned int n=1;
  while (n<=niters) {
    term/=n;
    term = A*term;
    
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

//' Compute Expected nij in CTMC
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v1 Vector (column vector) of to do v*exp(A)
//' @param v2 Vector (column vector) of to do exp(A)*v
//' @param itoj Integer matrix that indicates indices of non zero elements of A.
//' @param tau Time period between steps.
//' @param delta Time step for midpoint rule.
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param tau_only Logical to compute expected time only vs expected time and expected number of jumps.
//' @details Implements Sherlock 2021 uniformization method.
// [[Rcpp::export]]
arma::vec expected_n(const arma::sp_mat& A, const arma::vec& v1, const arma::vec& v2, const arma::umat& itoj, const double& tau, double delta, double tol, int renorm_freq, bool tau_only) {
  double rho=-A.min();
  if (rho==0) {
    return v1;
  }

  //' Setup midopint rule: v1*exp(As)exp(A(t-s))v2ds
  int ndiv = floor(tau / delta);
  delta = tau/ndiv;

  const int nr=A.n_rows, nc=A.n_cols, nij = itoj.n_rows;

  double rho_1 = rho*delta/2;
  double rho_2 = rho*delta;  
  const unsigned int niters_1 = ::Rf_qpois( tol, rho_1, false, false);
  const unsigned int niters_2 = ::Rf_qpois( tol, rho_2, false, false);

  // Key values to return:
  arma::vec expected_nij;
  if(!tau_only){ 
    expected_nij.set_size(nij);
    expected_nij.zeros();
  }
  arma::vec expected_tau(nr);

  arma::sp_mat P = A + rho*arma::speye(nr, nc);
  arma::sp_mat Pt = P.t();
  
  arma::dmat pforward(nr, ndiv);  
  pforward.col(0) = expAv_light(Pt*delta/2, v1, rho_1, niters_1, renorm_freq);
  for( int i=1; i<ndiv; ++i){
    pforward.col(i) =  expAv_light(Pt*delta, pforward.col(i-1), rho_2, niters_2, renorm_freq);
  }
  // The final step to complete expAv over the small final step:
  arma::vec pall = expAv_light(Pt*delta/2, pforward.col(ndiv-1), rho_1, niters_1, renorm_freq); 
  
  arma::vec preverse = expAv_light(P*delta/2, v2, rho_1, niters_1, renorm_freq);
  // .elem is ordered by stacking columns.
  if(!tau_only){ expected_nij = delta * pforward.elem((ndiv-1)*nr + itoj.col(0) - 1) % preverse(itoj.col(1) - 1); } // First step of integral using midpoint
  expected_tau = delta * pforward.col(ndiv-1) % preverse; // First step of integral using midpoint
  
  for( int i=1; i < ndiv; ++i) {
    preverse = expAv_light(P*delta, preverse, rho_2, niters_2, renorm_freq);   
    if(!tau_only){ expected_nij += delta * pforward.elem((ndiv-i-1)*nr + itoj.col(0) - 1) % preverse(itoj.col(1) - 1); }
    expected_tau += delta * pforward.col(ndiv-i-1) % preverse;
  }
  double pobs = sum(pall % v2);

  expected_tau /= pobs;
  if(tau_only){
    return expected_tau;    
  }
  expected_nij /= pobs;
  for(int i = 0; i<nij; ++i){
    expected_nij(i) *= A(itoj(i,0)-1, itoj(i,1)-1);
  }
  arma::vec out = join_cols(expected_tau, expected_nij);
  return out;
}

//' Compute int v1 expAs expA(t-s)  v2 ds
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v1 Vector (column vector)
//' @param v2 Vector (column vector)
//' @param tol Tolerance default = 1e-8.
//' @details Implements Sherlock 2021 uniformization method.
// [[Rcpp::export]]
arma::vec v1expAv2_old(const arma::sp_mat& A, const double& tau, const arma::vec& v1, const arma::vec& v2, const arma::umat& itoj, double tol) { 
  double rho=-A.min();
  if (rho==0) {
    return v1;
  }
  const int nr=A.n_rows, nc=A.n_cols, nij = itoj.n_rows;
  const int niters = ::Rf_qpois( tol, rho*tau, false, false);
  arma::sp_mat P = (A/rho+arma::speye(nr, nc));

  // Key values to return:
  arma::vec expected_nij(nij, fill::zeros);
  arma::vec expected_tau(nr, fill::zeros);
    
  arma::mat forward(nr, niters+1, arma::fill::none);
  arma::mat reverse(nr, niters+1, arma::fill::none);
  forward.col(0) = v1;
  reverse.col(0) = v2;
  arma::sp_mat Pt = P.t();

  // Run through products:
  for ( int n=1; n <= niters; ++n ) {
    forward.col(n) = Pt * forward.col(n-1);
    reverse.col(n) = P * reverse.col(n-1);
  }

  double dp = exp(-rho*tau);
  expected_tau = v1 % v2 * dp * tau;
  expected_nij = (v1.elem(itoj.col(0) - 1) % v2.elem(itoj.col(1) - 1)) * dp * tau;
  
  for( int n = 1; n <= niters; ++n ){
    dp *= rho*tau/n;
    if(dp > 1e-12){
      arma::vec inner_sum_tau = forward.col(0) % reverse.col(n);
      arma::vec inner_sum_nij = forward.elem(itoj.col(0) - 1) % reverse.elem(n*nr + itoj.col(1) - 1);
      for( int j = 1; j <= n; ++j ){
        inner_sum_tau += forward.col(j) % reverse.col(n-j);
        inner_sum_nij += forward.elem(j*nr + itoj.col(0) - 1) % reverse.elem((n-j)*nr + itoj.col(1) - 1);
      }
      expected_tau += inner_sum_tau * dp * tau/(n+1);
      expected_nij += inner_sum_nij * dp * tau/(n+1);
    }
  }
  // Figure out pobs:
  double pobs = sum(expected_tau)/tau;
  // Rcpp::Rcout << "pobs: " << pobs << std::endl;
  expected_tau /= pobs;
  expected_nij /= pobs;
  for(int i = 0; i<nij; ++i){
    expected_nij(i) *= A(itoj(i,0)-1, itoj(i,1)-1);
  }
  
  arma::vec out = join_cols(expected_tau, expected_nij);
  return out;
}

//' Compute int v1 expAs expA(t-s)  v2 ds
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v1 Vector (column vector)
//' @param v2 Vector (column vector)
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @details Implements Sherlock 2021 uniformization method.
// [[Rcpp::export]]
arma::vec EN(const arma::sp_mat& A, const double& tau, const arma::vec& v1, const arma::vec& v2, const arma::umat& itoj, double tol, bool tau_only) { 
  double rho=-A.min();
  if (rho==0) {
    return v1;
  }
  const int nr=A.n_rows, nc=A.n_cols, nij = itoj.n_rows;
  const int niters = ::Rf_qpois( tol, rho*tau, false, false);
  arma::sp_mat P = (A/rho+arma::speye(nr, nc));

  // Key values to return:
  arma::vec expected_nij;
  if(!tau_only){ 
    expected_nij.set_size(nij);
    expected_nij.zeros();
  }
  arma::vec expected_tau(nr, fill::zeros);
  
  arma::vec forward = v1;
  arma::mat reverse(nr, niters+1, arma::fill::none);
  reverse.col(0) = v2;
  arma::sp_mat Pt = P.t();

  arma::vec scale(niters+1);
  double lambda = rho*tau;
  scale(0) = exp(-lambda)*tau;
  for ( int n=1; n <= niters; ++n ) {
    reverse.col(n) = P * reverse.col(n-1);
    scale(n) = scale(n-1)*lambda/(n+1);
  }

  for( int n = 0; n <= niters; ++n ){
    if(n > 0){ 
      forward = Pt*forward;
    }
    for( int j = 0; j <= (niters-n); ++j){
      if(scale(n+j) > 1e-12){
        expected_tau += (forward % reverse.col(j)) * scale(n+j);
        if(!tau_only){ expected_nij += (forward.elem(itoj.col(0) - 1) % reverse.elem(j*nr + itoj.col(1) - 1)) * scale(n+j); }
      }
    }
  }
  // Figure out pobs:
  double pobs = sum(expected_tau)/tau;
  // Rcpp::Rcout << "pobs: " << pobs << std::endl;
  expected_tau /= pobs;
  if(tau_only){
    return expected_tau;
  }
  expected_nij /= pobs;
  for(int i = 0; i<nij; ++i){
    expected_nij(i) *= A(itoj(i,0)-1, itoj(i,1)-1);
  }
  arma::vec out = join_cols(expected_tau, expected_nij);
  return out;
}
