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

//' Fast Multiply of multiple sparse matrix vector.
//' @param dA Dense matrix with each column the non-zero values of a sparse matrix.
//' @param term Vector to multiply against matrix.
//' @param row_indx Row index from column-sparse matrix "dgCMatrix" \code{A@i}.
//' @param col_ptr Column pointer from column-sparse matrix "dgCMatrix" \code{A@p}.
// [[Rcpp::export]]
arma::mat fast_multiply(const arma::mat& dA, const arma::vec& term, const std::vector<int>& row_indx, const std::vector<int>& col_ptr){
  int npars = dA.n_cols;
  int ncols = term.n_elem; // Square matrix assumed with term length cols and rows.
  arma::mat dQw(ncols, npars, fill::zeros);

  for(int j = 0; j < ncols; ++j){
    for(int k = col_ptr[j]; k < col_ptr[j+1]; ++k){
      dQw.row(row_indx[k]) += dA.row(k) * term[j];
    }
  }
  return dQw;
}

//' Compute expAv and gradient of expAv jointly
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param dA A matrix that has row of non-zero gradient of A and columns for each parameter.
//' @param v Vector (column vector) of to do exp(A)*v
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param trans Logical to confirm if we should do transpose of A.
//' @param row_indx Row index from column-sparse matrix "dgCMatrix" \code{A@i}.
//' @param col_ptr Column pointer from column-sparse matrix "dgCMatrix" \code{A@p}.
//' @details Computes joint expAv and the gradient as implemented in 'Differentiated uniformization: a new method for inferring Markov chains on combinatorial state spaces including stochastic epidemic models'.
// [[Rcpp::export]]
Rcpp::List expAv_gr_cpp(const arma::sp_mat& A, const arma::mat& dA, const arma::vec& v, double tol, int renorm_freq, const std::vector<int>& row_indx, const std::vector<int>& col_ptr) {
  const double rho=-A.min();

  const int nr=A.n_rows, nc=A.n_cols;
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  unsigned int n;
  const arma::sp_mat P = A+rho*arma::speye(nr, nc);

  arma::vec term = v, ans = v;
  double log_scale = 0.0;
  double s = 0.0;

  const int npars = dA.n_cols; 

  arma::mat delta(nr, npars, fill::zeros);
  arma::mat gr_ans(nr, npars, fill::zeros);
  arma::mat d_term;

  n=1;
  while (n<=niters) {
    term/=n;
    d_term = fast_multiply(dA, term, row_indx, col_ptr);
    for(int i=0; i<npars; ++i){
      delta.col(i) = d_term.col(i) + P*delta.col(i)/n;
      gr_ans.col(i) += delta.col(i);
    }
    term = P*term;    
    ans = ans + term;
    
    if(n % renorm_freq == 0){
      s = sum(abs(ans));
      term/=s;
      ans/=s;
      delta/=s;
      gr_ans/=s;
      log_scale += log(s);
    }
    n++;
  }
  double correction = exp(-rho + log_scale);
  ans *= correction;
  gr_ans *= correction;
  return Rcpp::List::create(Rcpp::Named("ans") = ans, Rcpp::Named("gr_ans") = gr_ans);
}

//' Compute Linear approximation of expAv and gradient of expAv jointly
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v Vector (column vector) of to do exp(A)*v.
//' @param lambda Vector (column vector) of detection rate.
//' @param tol Tolerance default = 1e-8.
//' @param trans Logical to confirm if we should do transpose of A.
//' @details Computes joint expAv and the gradient as implemented in 'Differentiated uniformization: a new method for inferring Markov chains on combinatorial state spaces including stochastic epidemic models'.
// [[Rcpp::export]]
arma::vec expAv_approx_cpp(const arma::sp_mat& A, const arma::vec& v, 
                               const arma::vec& lambda, double tol) {
  const double rho=-A.min();

  const int nr=A.n_rows, nc=A.n_cols;
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  // const unsigned int niters = ceil(rho*2);
  unsigned int n;
  const arma::sp_mat R = A/niters+arma::speye(nr, nc);
  arma::vec term = v;
  const arma::vec pdet = 1-lambda/niters;//exp(-lambda / niters);
  
  n=1;
  while (n<=niters) {
    term %= pdet;
    term = R*term;    
    n++;
  }
  term %= pdet;
  return term;
}


// //' Compute Linear approximation of expAv and gradient of expAv jointly
// //' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
// //' @param dA A matrix that has row of non-zero gradient of A and columns for each parameter.
// //' @param v Vector (column vector) of to do exp(A)*v
// //' @param pdet Vector (column vector) of exp(-Lambda)
// //' @param dpdet Vector (column vector) of gradient of exp(-Lambda)
// //' @param tol Tolerance default = 1e-8.
// //' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
// //' @param trans Logical to confirm if we should do transpose of A.
// //' @param row_indx Row index from column-sparse matrix "dgCMatrix" \code{A@i}.
// //' @param col_ptr Column pointer from column-sparse matrix "dgCMatrix" \code{A@p}.
// //' @details Computes joint expAv and the gradient as implemented in 'Differentiated uniformization: a new method for inferring Markov chains on combinatorial state spaces including stochastic epidemic models'.
//// [[Rcpp::export]]
// Rcpp::List expAv_gr_approx_cpp(const arma::sp_mat& A, arma::mat dA, const arma::vec& v, 
                               // const arma::vec& pdet, const arma::vec& dpdet, double tol,
                               // int renorm_freq, const std::vector<int>& row_indx, const std::vector<int>& col_ptr) {
  // const double rho=-A.min();

  // const int nr=A.n_rows, nc=A.n_cols;
  // const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  // unsigned int n;
  // const arma::sp_mat R = A/niters+arma::speye(nr, nc);
  // dA /= niters;
  // arma::vec term = v;

  // const int npars = dA.n_cols; 

  // arma::mat delta(nr, npars, fill::zeros);
  // arma::mat gr_ans(nr, npars, fill::zeros);
  // arma::mat d_term;

  // n=1;
  // while (n<=niters) {
    // term*=pdet;
    // d_term = fast_multiply(dA, term, row_indx, col_ptr); // Need to loop through all theta except logitq.
    // for(int i=0; i<npars; ++i){
      // delta.col(i) = d_term.col(i) + R*delta.col(i);  // Need to add dpdet to this. Will be something like: dR*p*term + R*dp*v + Rp dv, where dv is accumulated.
    // }
    // term = R*term;    
    // n++;
  // }
  // return Rcpp::List::create(Rcpp::Named("ans") = term, Rcpp::Named("gr_ans") = delta);
// }
