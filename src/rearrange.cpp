#include <Rcpp.h>
#include <TreeTools.h>
using namespace std;
using namespace Rcpp;

//  [[Rcpp::export]]
List asan_error (const IntegerMatrix x) {
  const int
    n_edge = x.nrow(),
    n_internal = n_edge / 2,
    n_tip = n_internal + 1,
    n_vert = n_internal + n_tip,
    root_node = n_tip + 1
  ;
  
  Rf_error("Oh dear.");
  return List::create();
}
