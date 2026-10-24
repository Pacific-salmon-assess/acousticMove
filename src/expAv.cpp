#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::export]]
void transpose_in_place(arma::sp_mat& X) {
  arma::inplace_trans(X); 
}

// [[Rcpp::export]]
void kahan_sum(arma::vec& full_sum, const arma::mat& A) {
  const int npar = A.n_cols;
  full_sum.set_size(npar);
  for(int i = 0; i < npar; ++i){
    double ans = 0.0;
    volatile double c = 0.0; // Prevent compiler optimization
    for (double val : A.col(i)) {
      volatile double y = val - c;
      volatile double t = ans + y;
      c = (t - ans) - y;
      ans = t;
    }
    full_sum(i) = ans;
  }
}

// [[Rcpp::export]]
void expAv_internal(arma::vec& ans, const arma::sp_mat& A, const arma::vec& v, const double tol, const int renorm_freq, const bool trans = false, double vmin = 0.0, int state = 0) {
  const double rho=-A.min();
  const int nr=A.n_rows, nc=A.n_cols;

  if (rho==0) {
    return;
  }
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  unsigned int n;
  arma::sp_mat P = A+rho*arma::speye(nr, nc);
    
  arma::vec term = v;
  ans = v;

  if(trans) transpose_in_place(P);

  const double* values = P.values;
  const arma::uword* row_indx = P.row_indices;
  const arma::uword* col_ptr = P.col_ptrs;

  double log_scale = 0.0;
  double s = 0.0;
  // Rcout << "rho: " << rho << std::endl;
  // Rcout << "N: " << niters << std::endl;
  arma::vec next_term(nr, fill::zeros);
  
  n=1;
  while (n<=niters) {
    next_term.zeros();
    for(int i=0; i<nr; ++i){
      if(term(i) > vmin){
        const double x = term(i) / static_cast<double>(n);
        for(arma::uword k = col_ptr[i]; k < col_ptr[i+1]; ++k){
          next_term[row_indx[k]] += values[k] * x;
        }
      }
    }
    term.swap(next_term);
    ans = ans + term;
    
    if(n % renorm_freq == 0){
      s = sum(abs(ans));
      term/=s;
      ans/=s;
      log_scale += log(s);
    }
    n++;
  }
  if(state > 0){
    ans = exp(-rho + log_scale) * ans(state-1);   
  }else{
    ans = exp(-rho + log_scale) * ans;
  }
}

//' Compute expAv via telescope
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v Vector (column vector) of to do exp(A)*v
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param trans Logical to confirm if we should do transpose of A.
//' @param vmin Additional approximation by excluding some of the animal's path. Suggest keeping 0.
//' @details Implements Sherlock 2021 uniformization method, differs from expAv regular by taking advantage of sparsity of v.
// [[Rcpp::export]]
arma::vec expAv_cpp(const arma::sp_mat& A, const arma::vec& v, const double tol, const int renorm_freq, const bool trans = false, double vmin = 0.0) {
  arma::vec ans;
  expAv_internal(ans, A, v, tol, renorm_freq, trans, vmin, 0);
  return ans;
}

// [[Rcpp::export]]
void expAv_gr_internal(arma::vec& ans, arma::mat& gr_ans, const arma::sp_mat& A, const arma::mat& dA, const arma::vec& v, const double tol, const int renorm_freq, double vmin = 0.0, int state = 0) {
  const double rho=-A.min();

  const int nr=A.n_rows, nc=A.n_cols;
  const unsigned int niters = ::Rf_qpois( tol, rho, false, false);
  unsigned int n;
  const arma::sp_mat P = A+rho*arma::speye(nr, nc);

  arma::vec term = v;
  ans = v;
  
  double log_scale = 0.0;
  arma::vec next_term(nr, arma::fill::zeros);

  const int npars = dA.n_cols; 

  arma::mat delta(nr, npars, fill::zeros);
  arma::mat next_delta(nr, npars, arma::fill::zeros);  
  gr_ans.set_size(nr, npars);
  gr_ans.zeros();

  // Telescoping components via sparseness:
  const double* values = P.values;
  const arma::uword* row_indx = P.row_indices;
  const arma::uword* col_ptr = P.col_ptrs;
  
  const arma::uword* Arow_indx = A.row_indices;
  const arma::uword* Acol_ptr = A.col_ptrs;    
    
  n=1;
  while (n<=niters) {
    // const double inv_n = 1.0 / static_cast<double>(n);    
    
    next_term.zeros();
    next_delta.zeros();
    
    for (int i = 0; i < nc; ++i) {
      const arma::uword begin = col_ptr[i];
      const arma::uword end   = col_ptr[i + 1];

      // Step 1: P*delta.col(j) for each npars
      for( int j = 0; j < npars; ++j){
        if(delta(i,j) != 0.0){
          const double x = delta(i,j) / static_cast<double>(n);
          for (arma::uword k = begin; k < end; ++k) {
              const arma::uword r = row_indx[k];
              next_delta(r,j) += values[k] * x;
          }
        }
      }

      // // if ti is negative things have gone very poorly.
      if (term[i] <= vmin) {
          continue;
      }
      const double ti = term[i] / static_cast<double>(n);
      // Now do P*term:
      for (arma::uword k = begin; k < end; ++k) {
          const arma::uword r = row_indx[k];
          const double p = values[k];

          next_term[r] += p * ti;
      }
      // Sparsity changes bewteen A and P. So have to
      // to adapt for dA * term.
      const arma::uword abegin = Acol_ptr[i];
      const arma::uword aend = Acol_ptr[i + 1];
      
      for (arma::uword k = abegin; k < aend; ++k) {
        const arma::uword Ar = Arow_indx[k];
         next_delta.row(Ar) += dA.row(k) * ti;
      }
    }
    term.swap(next_term);
    delta.swap(next_delta);
    
    gr_ans += delta;
    ans += term;
    
    if(n % renorm_freq == 0){
      const double s = arma::accu(arma::abs(ans));
      term/=s;
      ans/=s;
      delta/=s;
      gr_ans/=s;
      log_scale += std::log(s);
    }
    n++;
  }
  const double correction = exp(-rho + log_scale);
  
  if(state > 0){
    ans = correction * ans(state-1);  
    gr_ans = correction * gr_ans.row(state-1);
  }else{
    ans *= correction;
    gr_ans *= correction;
  }
}

//' Compute expAv and gradient of expAv jointly via telescope
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param dA A matrix that has row of non-zero gradient of A and columns for each parameter.
//' @param v Vector (column vector) of to do exp(A)*v
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param trans Logical to confirm if we should do transpose of A.
//' @param vmin Additional approximation by excluding some of the animal's path.
//' @details Computes joint expAv and the gradient as implemented in 'Differentiated uniformization: a new method for inferring Markov chains on combinatorial 
// 'state spaces including stochastic epidemic models'. Differs from expAv regular by taking advantage of sparsity of v.
// [[Rcpp::export]]
Rcpp::List expAv_gr_cpp(const arma::sp_mat& A, const arma::mat& dA, const arma::vec& v, const double tol, const int renorm_freq, double vmin = 0.0) {
  arma::vec ans;
  arma::mat gr_ans;
  
  expAv_gr_internal(ans, gr_ans, A, dA, v, tol, renorm_freq, vmin, 0);

  return Rcpp::List::create(Rcpp::Named("ans") = ans, Rcpp::Named("gr_ans") = gr_ans);
}

//' Compute Compute survival function for ctmc with known states, exp(A)_ij
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param v_idx What state is the animal in, relating to j.
//' @param v_other What state is the animal coming from, or going to, (e.g. i).
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param vmin Additional approximation by excluding some of the animal's path.
//' @param log Return the log or the real scale (default log = FALSE).
//' @details Implements Sherlock 2021 uniformization method, differs from expAv regular by taking advantage of sparsity of v. Reduces memory as a wrapper to fetch ij. When
//' v_other = 0, we assume right censoring and sum values.
// [[Rcpp::export]]
double survival_ctmc(const arma::sp_mat& A, int v_idx = 0, int v_other = 0, double tol = 1e-8, int renorm_freq = 10, double vmin = 0.0, bool log = false) {
  const int nr = A.n_rows;
  arma::vec v(nr, arma::fill::zeros);
  v[v_idx - 1] = 1;
  arma::vec ans;

  expAv_internal(ans, A, v, tol, renorm_freq, false, vmin, v_other);
  double surv;
    
  if(v_other > 0){
    surv = ans(0);
    if(log) surv = std::log(surv);
    return surv;
  }
  // Strictly positive but can be some very big values But there can be some very small values.
  arma::vec log_ans = arma::log(ans);
  double max_val = log_ans.max();
  double log_surv = max_val + std::log(arma::sum(arma::exp(log_ans - max_val)));
  if(log) return log_surv;
  surv = exp(log_surv);
  return surv;
}

//' Compute Compute survival function for ctmc with known states, exp(A)_ij and the gradient.
//' @param A A sparse matrix from R of \code{class(A) == "dgCMatrix"}.
//' @param dA A matrix that has row of non-zero gradient of A and columns for each parameter.
//' @param v_idx What state is the animal in, relating to j.
//' @param v_other What state is the animal coming from, or going to, (e.g. i).
//' @param tol Tolerance default = 1e-8.
//' @param renorm_freq Renoramlization frequency to avoid computational overload from rho^n/n!.
//' @param vmin Additional approximation by excluding some of the animal's path.
//' @param log Return the log or the real scale (default log = FALSE).
//' @details Computes joint expAv and the gradient as implemented in 'Differentiated uniformization: a new method for inferring Markov chains on combinatorial 
// 'state spaces including stochastic epidemic models'. Differs from expAv regular by taking advantage of sparsity of v. When
//' v_other = 0, we assume right censoring and sum values. Otherwise returns the exp(A)ij and the gradient.
// [[Rcpp::export]]
Rcpp::List survival_gr_ctmc(const arma::sp_mat& A, const arma::mat& dA, int v_idx = 0, int v_other = 0, double tol = 1e-8, int renorm_freq = 10, double vmin = 0.0, bool log = false) {
  const int nr = A.n_rows;
  arma::vec v(nr, arma::fill::zeros);
  v[v_idx - 1] = 1;
  
  arma::vec ans;
  arma::mat gr_ans;

  expAv_gr_internal(ans, gr_ans, A, dA, v, tol, renorm_freq, vmin, v_other);
  
  if(v_other > 0){
    if(log) {
      gr_ans /= ans(0);
      ans = arma::log(ans);
    }
    return Rcpp::List::create(Rcpp::Named("ans") = ans, Rcpp::Named("gr_ans") = gr_ans);  
  }
  // log sum exponential: for small positive values
  arma::vec log_ans = arma::log(ans);
  double max_val = log_ans.max();
  double log_surv = max_val + std::log(arma::sum(arma::exp(log_ans - max_val)));
  double surv = exp(log_surv);

  // Column sum of the matrix. Use Kahan Sum for small double values. Matters... but maybe not totally neccesary.
  arma::vec surv_gr;
  kahan_sum(surv_gr, gr_ans);
  // arma::rowvec surv_gr = arma::sum(gr_ans,0);
  if(log){
    surv_gr /= surv;
    surv = log_surv;
  }

  return Rcpp::List::create(Rcpp::Named("ans") = surv, Rcpp::Named("gr_ans") = surv_gr);  
}

