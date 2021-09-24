#include <Rcpp.h>
#include <TreeTools.h>
using namespace std;
using namespace Rcpp;

//  [[Rcpp::export]]
List asan_error (const IntegerMatrix x) {
  Rf_error("Oh dear.");
  return List::create();
}
