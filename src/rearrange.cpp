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
IntegerMatrix nni(const IntegerMatrix edge,
                  const IntegerVector randomEdge,
                  const IntegerVector whichSwitch) {
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

// #TODO move to TreeTools and import
// edge must be in preorder
//  [[Rcpp::export]]
IntegerMatrix root_on_node(const IntegerMatrix edge, int outgroup) {
  if (edge(0, 1) == outgroup || edge(1, 1) == outgroup) return edge;
  
  const int16 n_edge = edge.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    root_node = n_tip + 1,
    max_node = edge(n_edge - 1, 0);
  if (outgroup > max_node) throw std::range_error("outgroup exceeds number of nodes");
  if (outgroup == root_node) return edge;
  int16* edge_above = new int16[max_node + 1];
  int16 root_edges[2] = {0, 0};
  for (int16 i = n_edge; i--; ) {
    edge_above[edge(i, 1)] = i;
    if (edge(i, 0) == root_node) {
      root_edges[root_edges[1] ? 0 : 1] = i;
    }
  }
  
  IntegerMatrix ret = clone(edge);
  int16 invert_next = edge_above[outgroup];
  
  // We'll later add an edge from the now-unallocated root node to the outgroup.
  ret(invert_next, 0) = root_node;
  ret(invert_next, 1) = edge(invert_next, 0);
  do {
    Rcout << "\n - Inverted edge " << (invert_next + 1) << " with parent " <<
      edge(invert_next, 0) << " which is child of " << edge_above[edge(invert_next, 0)];
    invert_next = edge_above[edge(invert_next, 0)];
    ret(invert_next, 0) = edge(invert_next, 1);
    ret(invert_next, 1) = edge(invert_next, 0);
  } while (edge(invert_next, 0) != root_node);
  Rcout << "\n * Invert next left at: " << invert_next << "\n";
  
  delete[] edge_above;
  
  
  // second root i.e. 16 -- 24 must be replaced with root -> outgroup.
  int16 spare_edge = (ret(root_edges[0], 0) == root_node ? 0 : 1);
  Rcout << "Spare edge: " << (root_edges[spare_edge] + 1) << "\n";
  ret(invert_next, 1) = edge(root_edges[spare_edge], 1);
  ret(root_edges[spare_edge], 1) = outgroup;
  
  for (int i = 0; i < n_edge; i++) {
    Rcout << "\n["<< (i + 1)<< ",] " << ret(i, 0) << " -- " << ret(i, 1);
  }
  return preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
  return ret;
}


// Assumptions: 
//  * Tree is bifurcating, in preorder; first two edges have root as parent.
//  [[Rcpp::export]]
IntegerMatrix spr (const IntegerMatrix edge,
                   const IntegerMatrix randomEdge,
                   const IntegerMatrix mergeEdge) {
  const int16
    n_edge = edge.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    root_node = n_tip + 1,
    broken_option = randomEdge[0] % (n_edge - 1), // Don't break root sister
    broken_edge = broken_option == 0 ? 0 : broken_option + 1,
    chosen_regraft = mergeEdge[0]
  ;
  if (n_edge < 5) return edge;
  if (edge(0, 0) != root_node) throw std::invalid_argument("First edge of tree must connect root to leaf. Try Preorder(root(tree)).");
  if (edge(1, 0) != root_node) throw std::invalid_argument("Second edge of tree must connect root to leaf. Try Preorder(root(tree)).");
  
  const int16 broken_edge_parent = edge(broken_edge, 0);
  
  IntegerMatrix ret = clone(edge);
  
  if (broken_edge) { // We are breaking a non-root edge
    int16 edge_above_broken = 0, edge_beside_broken = 0;
    
    int16* merge_options = new int16[n_edge];
    int16 n_merge_options = 0, i = 0;
    bool adrift = false;
    for (i = 2; i != n_edge; i++) {
      if (edge(i, 1) == broken_edge_parent) {
        edge_above_broken = i;
        continue;
      }
      if (i == broken_edge) {
        adrift = true;
        continue;
      }
      if (adrift) {
        if (edge(i, 0) == broken_edge_parent) {
          break; // Now we know that all remaining edges will be potential merge sites
        }
      } else {
        if (edge(i, 0) == broken_edge_parent) {
          edge_beside_broken = i;
        } else {
          merge_options[n_merge_options++] = i;
        }
      }
    }
    if (!edge_beside_broken) {
      edge_beside_broken = i;
    }
    while (i != n_edge) {
      merge_options[n_merge_options++] = i++;
    }
    const int16 merge_edge = merge_options[chosen_regraft % n_merge_options];
    delete[] merge_options;
    
    ret(edge_beside_broken, 0) = edge(edge_above_broken, 0);
    ret(edge_above_broken, 0) = edge(merge_edge, 0);
    ret(merge_edge, 0) = broken_edge_parent;
    return preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
  } else { // We are breaking the root edge
    
    
    int16 root_daughter_2 = 0, merge_options_considered = 0, merge_edge = -1;
    const int16
      second_root_child = root_node + 1,
      chosen_merge_option = chosen_regraft % (n_edge - 4)
    ;
    
    for (int16 i = 3; 
         merge_options_considered < chosen_merge_option && root_daughter_2;
         i++) {
      if (edge(i, 0) == second_root_child) {
        root_daughter_2 = i;
      } else {
        if (chosen_merge_option == merge_options_considered++) {
          merge_edge = i;
        }
      }
    }
    if (merge_edge == -1) throw std::invalid_argument("Merge location not found. Sack the programmer.");
    
    ret(2, 0) = broken_edge_parent;
    ret(root_daughter_2, 0) = broken_edge_parent;
    
    //child [brokenEdgeSister] <- child[mergeEdge]
    ret(1, 1) = edge(merge_edge, 1);
    //parent[brokenEdge | brokenEdgeSister] <- spareNode
    ret(0, 0) = 2;
    ret(1, 0) = 2;
    // child[mergeEdge] <- spareNode
    ret(merge_edge, 1) = 2;
    ret = preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
    return root_on_node(ret, 1);
  }
}
