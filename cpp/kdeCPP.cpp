// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec kde_multivariate(const arma::mat& data, const arma::vec& point, const arma::mat& H) {
  int n = data.n_rows;
  int d = data.n_cols;
  arma::vec kde_values(n);
  arma::mat H_inv = arma::inv(H);
  double det_H = arma::det(H);
  double norm_const = pow(2 * M_PI, -0.5 * d) * pow(det_H, -0.5);
  
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    arma::vec diff = data.row(i).t() - point;
    kde_values(i) = exp(-0.5 * arma::as_scalar(diff.t() * H_inv * diff));
  }
  
  return norm_const * kde_values;
}

// [[Rcpp::export]]
double kde_multivariate_sum( const arma::mat& data, 
                             const arma::vec& point, 
                             const arma::mat& H ) {
  int n = data.n_rows;
  int d = data.n_cols;
  arma::vec kde_values(n);
  // double out_sum = 0;
  arma::mat H_inv = arma::inv(H);
  double det_H = arma::det(H);
  double norm_const = pow(2 * M_PI, -0.5 * d) * pow(det_H, -0.5);
  
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    arma::vec diff = data.row(i).t() - point;
    // out_sum += exp(-0.5 * arma::as_scalar(diff.t() * H_inv * diff));
    kde_values(i) = exp(-0.5 * arma::as_scalar(diff.t() * H_inv * diff));
  }
  return norm_const * arma::accu(kde_values) / n;
}

// [[Rcpp::export]]
arma::vec estimate_kde_for_pts( const arma::mat& pts,
                                const arma::mat& data,
                                const arma::mat& H )
{
  int n_pt  = pts.n_rows;
  int n_dat = data.n_rows;
  int dim   = pts.n_cols; // should be == data.n_cols
  arma::vec estimated_vals(n_pt);
  // arma::vec kde_values(n_dat);
  
  arma::mat H_inv = arma::inv(H);
  double det_H = arma::det(H);
  double norm_const = pow(2 * M_PI, -0.5 * dim) * pow(det_H, -0.5);

#pragma omp parallel for  
  for( int pid = 0; pid < n_pt; ++pid )
  {
    arma::vec point = pts.row(pid).t();
    arma::vec kde_values(n_dat);

    for (int did = 0; did < n_dat; ++did) {
      arma::vec diff = data.row(did).t() - point;
      // out_sum += exp(-0.5 * arma::as_scalar(diff.t() * H_inv * diff));
      kde_values(did) = exp(-0.5 * arma::as_scalar(diff.t() * H_inv * diff));
    }
    estimated_vals(pid) = norm_const * arma::accu(kde_values) / n_dat;    
  }
  return estimated_vals;
}
// with weights
// [[Rcpp::export]]
arma::vec estimate_kde_for_pts_wgt( const arma::mat& pts,
                                    const arma::mat& data,
                                    const arma::mat& H,
                                    const arma::vec& weights ) // Add weights parameter
{
  int n_pt  = pts.n_rows;
  int n_dat = data.n_rows;
  int dim   = pts.n_cols; // should be == data.n_cols
  arma::vec estimated_vals(n_pt);
  
  arma::mat H_inv = arma::inv(H);
  double det_H = arma::det(H);
  double norm_const = pow(2 * M_PI, -0.5 * dim) * pow(det_H, -0.5);
  
#pragma omp parallel for  
  for(int pid = 0; pid < n_pt; ++pid)
  {
    arma::vec point = pts.row(pid).t();
    arma::vec kde_values(n_dat);
    
    for(int did = 0; did < n_dat; ++did) {
      arma::vec diff = data.row(did).t() - point;
      kde_values(did) = weights(did) * exp(-0.5 * arma::as_scalar(diff.t() * H_inv * diff)); // Apply weights
    }
    estimated_vals(pid) = norm_const * arma::accu(kde_values) / arma::accu(weights); // Normalize by sum of weights
  }
  return estimated_vals;
}