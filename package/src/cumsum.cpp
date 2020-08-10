#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix apply_cumsum_col(NumericMatrix m) {
    for (int j = 0; j < m.ncol(); ++j) {
        for (int i = 1; i < m.nrow(); ++i) {
            m(i, j) += m(i - 1, j);
        }
    }
    return m;
}
