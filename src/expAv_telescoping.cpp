#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace arma;

//' Overwrite a sparse matrix with its transpose
//' @param X A sparse matrix from R
//' @details Basic function that forces a transpose on a `arma::sp_mat` matrix without generating two full in memory.
// [[Rcpp::export]]
void inplace_transpose_sparse(arma::sp_mat& X) {
    // Step 1: Create a temporary sparse matrix holding the transpose
    arma::sp_mat tmp = X.t();
    
    // Step 2: Overwrite the original matrix memory
    X = std::move(tmp); 
}

//' Compute expAv
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v Vector (column vector) of to do exp(A)*v
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param trans Logical to confirm if we should do transpose of A.
//' @details Implements Sherlock 2021 uniformization method.
// [[Rcpp::export]]
arma::vec expAv_cpp(const arma::sp_mat& A, const arma::vec& v, double tol, int renorm_freq, bool trans = false) {
  const double rho=-A.min();
  if (rho==0) {
    return v;
  }
  const int nr=A.n_rows, nc=A.n_cols;
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  unsigned int n;
  arma::sp_mat P = A+rho*arma::speye(nr, nc);
    
  arma::vec term = v, ans = v;

  if(trans) inplace_transpose_sparse(P);

  double log_scale = 0.0;
  double s = 0.0;
  // Rcout << "rho: " << rho << std::endl;
  // Rcout << "N: " << niters << std::endl;

  n=1;
  while (n<=niters) {
    term/=n;
    term = P*term;
    
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


// [[Rcpp::export]]
arma::sp_mat expAv_sparsev_cpp(const arma::sp_mat& A, const arma::sp_mat& v, double tol, int renorm_freq, bool trans = false) {
  const double rho=-A.min();
  if (rho==0) {
    return v;
  }
  const int nr=A.n_rows, nc=A.n_cols;
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  unsigned int n;
  arma::sp_mat P = A+rho*arma::speye(nr, nc);
    
  arma::sp_mat term = v, ans = v;

  if(trans) inplace_transpose_sparse(P);

  double log_scale = 0.0;
  double s = 0.0;
  // Rcout << "rho: " << rho << std::endl;
  // Rcout << "N: " << niters << std::endl;

  n=1;
  while (n<=niters) {
    term/=n;
    term = P*term;
    
    ans = ans + term;
    
    if(n % renorm_freq == 0){
      s = arma::as_scalar(arma::sum(abs(ans.col(0))));
      term/=s;
      ans/=s;
      log_scale += log(s);
    }
    n++;
  }
  ans = exp(-rho + log_scale) * ans;
  return ans;
}


//' Compute expAv
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v Vector (column vector) of to do exp(A)*v
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param trans Logical to confirm if we should do transpose of A.
//' @details Implements Sherlock 2021 uniformization method.
// [[Rcpp::export]]
arma::vec expAv_telescope_cpp(const arma::sp_mat& A, const arma::vec& v, double tol, int renorm_freq, bool trans = false, const double& vmin = 0.0) {
  const double rho=-A.min();
  const int nr=A.n_rows, nc=A.n_cols;

  if (rho==0) {
    return v;
  }
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  unsigned int n;
  arma::sp_mat P = A+rho*arma::speye(nr, nc);
    
  arma::vec term = v, ans = v;

  if(trans) inplace_transpose_sparse(P);

  const double* values = P.values;
  const arma::uword* row_indx = P.row_indices;
  const arma::uword* col_ptr = P.col_ptrs;

  double log_scale = 0.0;
  double s = 0.0;
  // Rcout << "rho: " << rho << std::endl;
  // Rcout << "N: " << niters << std::endl;
  arma::vec all_zero(nr, fill::zeros);
  arma::vec new_term = all_zero;
  
  n=1;
  while (n<=niters) {
    new_term = all_zero;
    for(int i=0; i<nr; ++i){
      if(term(i) > vmin){
        term(i)/=n;
        for(int k = col_ptr[i]; k < col_ptr[i+1]; ++k){
          new_term[row_indx[k]] += values[k] * term[i];
        }        
      }
    }
    term = new_term;
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