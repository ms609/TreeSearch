#include <Rcpp.h>

extern "C" {
#include "mpl.h"
#include "RMorphy.h"
}

using namespace Rcpp;

// [[Rcpp::export]]
int preorder_morphy(IntegerMatrix edge, SEXP MorphyHandl) {
  Morphy handl = R_ExternalPtrAddr(MorphyHandl);
  const int
    n_tip = mpl_get_numtaxa(handl),
    n_internal = mpl_get_num_internal_nodes(handl),
    n_vertex = n_tip + n_internal,
    root_node = n_tip
  ;
  
  IntegerVector parent_of(n_vertex);
  IntegerVector left_child(n_internal);
  IntegerVector right_child(n_internal);
  
  for (int i = edge.nrow(); i--; ) {
    const int
      parent = edge(i, 0) - 1,
      child = edge(i, 1) - 1
    ;
    parent_of[child] = parent;
    if (right_child[parent - n_tip]) {
      left_child[parent - n_tip] = child;
    } else {
      right_child[parent - n_tip] = child;
    }
  }
  parent_of[0] = root_node;
  
  const int 
    /* INTEGER gives pointer to first element of an R vector */
    *ancestor = parent_of.begin(),
    *left = left_child.begin(),
    *right = right_child.begin()
  ; 
  
  
  /* Initialize return variables */
  int score = 0;
  morphy_length(ancestor, left, right, handl, &score); /* Updates score */
  return score;
}

// [[Rcpp::export]]
IntegerVector preorder_morphy_by_char(IntegerMatrix edge, List MorphyHandls) {
  Morphy handl = R_ExternalPtrAddr(MorphyHandls[0]);
  const int
    n_tip = mpl_get_numtaxa(handl),
    n_internal = mpl_get_num_internal_nodes(handl),
    n_vertex = n_tip + n_internal,
    root_node = n_tip
  ;
  
  IntegerVector parent_of(n_vertex);
  IntegerVector left_child(n_internal);
  IntegerVector right_child(n_internal);
  
  for (int i = edge.nrow(); i--; ) {
    const int
      parent = edge(i, 0) - 1,
      child = edge(i, 1) - 1
    ;
    parent_of[child] = parent;
    if (right_child[parent - n_tip]) {
      left_child[parent - n_tip] = child;
    } else {
      right_child[parent - n_tip] = child;
    }
  }
  parent_of[0] = root_node;
  
  const int 
    /* INTEGER gives pointer to first element of an R vector */
    *ancestor = parent_of.begin(),
    *left = left_child.begin(),
    *right = right_child.begin()
  ; 
  
  
  /* Initialize return variables */
  IntegerVector ret (MorphyHandls.length());
  for (int i = MorphyHandls.length(); i--; ) {
    int score = 0;
    Morphy handl = R_ExternalPtrAddr(MorphyHandls[i]);
    morphy_length(ancestor, left, right, handl, &score); /* Updates score */
    ret[i] = score;
  }
  return ret;
}
