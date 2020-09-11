#include <Rcpp.h>
/// [ [Rcpp::depends(TreeTools)]]
#include "TreeTools.h"
using namespace Rcpp;
/// [ [Rcpp::interfaces(r, cpp)]]

typedef int_fast16_t int16;
const int16 UNDEFINED = -1;

int16 sample_one(int16 len) {
  // TODO use more sophisticated RNG
  return 0;
}

// Assumptions: 
//  * Tree is bifurcating and rooted on a tip.
//  [[Rcpp::export]]
IntegerMatrix nni(IntegerMatrix edge,
                  IntegerVector randomEdge,
                  IntegerVector whichSwitch) {
  const int16 n_edge = edge.nrow(),
    chosen_edge = randomEdge[0],
    chosen_switch = whichSwitch[0] % 2,
    n_tip = (n_edge / 2) + 1;
    
  int16 n_samplable = 0;
  int16 *samplable = new int16[n_edge];
  for (int16 i = n_edge; --i; ) {
    if (edge(i, 1) > n_tip) {
      samplable[n_samplable++] = i;
    }
  }
  if (!n_samplable) {
    throw std::length_error("Not enough edges to allow NNI rearrangement");
  }
   
  const int16
    edge_to_break = samplable[chosen_edge % n_samplable],
    end1 = edge(edge_to_break, 0),
    end2 = edge(edge_to_break, 1);
  delete[] samplable;
  int16
    ind1 = UNDEFINED,
    ind2 = UNDEFINED;
  
  for (int16 i = n_edge; i--; ) {
    if (i != edge_to_break && edge(i, 0) == end1) {
      ind1 = i;
      break;
    }
  }
  
  for (int16 i = n_edge; i--; ) {
    if (edge(i, 0) == end2) {
      if (ind2 != UNDEFINED || chosen_switch) {
        ind2 = i;
        break;
      } else {
        ind2 = i;
      }
    }
  }
  
  IntegerMatrix ret = clone(edge);
  ret(ind1, 1) = edge(ind2, 1);
  ret(ind2, 1) = edge(ind1, 1);
  
  return preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
}
