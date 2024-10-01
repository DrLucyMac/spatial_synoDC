#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

// [[Rcpp::export]]
std::vector<arma::umat> get_idx_cpp(
    arma::mat& X
) {
    size_t ncells = X.max(); 
    std::vector<arma::umat> res(ncells); 
    for (int i = 1; i <= ncells; i++) {
        res[i-1] = ind2sub(size(X), find(X == i)) + 1; 
    }
    return res;
}
