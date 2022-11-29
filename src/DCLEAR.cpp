#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat dist_w_ham(const arma::mat X, const arma::vec InfoW) {

  const int n_cells = X.n_rows;
  const int n_loc = X.n_cols;
  const int n_InfoW = InfoW.size();

  arma::mat D = arma::zeros(n_cells, n_cells);
  
  for(int i = 0; i < n_cells ; i++ ) {
    for(int j = i + 1; j < n_cells ; j++) {
      for(int k = 0; k < n_loc ; k++) {
	if( X(i,k) != X(j,k)) {
	  D(i,j) += InfoW[ X(i,k)-1 ] * InfoW[ X(j,k)-1 ];
	}
      }
    }
  }
  
  return D;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat dist_w_ham2(const arma::mat X, const arma::vec InfoW) {

  const int n_cells = X.n_rows;
  const int n_loc = X.n_cols;
  const int n_InfoW = InfoW.size();

  arma::mat D = arma::zeros(n_cells, n_cells);

  
  for(int i = 0; i < n_cells ; i++ ) {
    for(int j = i + 1; j < n_cells ; j++) {

      int flag_i = 0;
      int flag_j = 0;
      int m_i = 0;
      int m_j = 0;
      int d_c = 0;
      int n_nonmissing = 0;
      for(int k = 0; k < n_loc ; k++) {

	if((flag_i == 0) && (X(i,k) == 2) && (X(j,k) != 2) ) {
	  flag_i = 1;
	  D(i,j) += InfoW[1];
	}
	if((flag_j == 0) && (X(j,k) == 2) && (X(i,k) != 2) ) {
	  flag_j = 1;
	  D(i,j) += InfoW[1];
	}
	if((flag_i == 1) && (X(i,k) != 2)) {
	  flag_i = 0;
	  m_i += 1;
	}
	if((flag_j == 1) && (X(j,k) != 2)) {
	  flag_j = 0;
	  m_j += 1;
	}

	if((X(i,k) != 2) && (X(j,k) != 2)) {
	  n_nonmissing += 1;
	}
	
	if( (X(i,k) != X(j,k)) && (X(i,k) != 2) && (X(j,k) !=2) ) {
	  D(i,j) += InfoW[ X(i,k)-1 ] * InfoW[ X(j,k)-1 ];
	}

	if(n_nonmissing + m_i + m_j == 0) {
	  n_nonmissing = n_loc/2;
	}
		
      }
    }
  }
  
  return D;
}

