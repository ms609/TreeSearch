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
  
  if (edge(0, 1) == outgroup) return edge;
  
  const int16 n_edge = edge.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    root_node = n_tip + 1,
    max_node = n_node + n_tip;
  
  if (outgroup < 1) throw std::range_error("`outgroup` must be a positive integer");
  if (outgroup > max_node) throw std::range_error("`outgroup` exceeds number of nodes");
  if (outgroup == root_node) return edge;
  
  int16* edge_above = new int16[max_node + 1];
  int16 root_edges[2] = {0, 0};
  
  for (int16 i = n_edge; i--; ) {
    
    edge_above[edge(i, 1)] = i;
    
    if (edge(i, 0) == root_node) {
      if (edge(i, 1) == outgroup) {
        delete[] edge_above;
        return edge;
      }
      root_edges[root_edges[1] ? 0 : 1] = i;
    }
    
  }
  
  IntegerMatrix ret = clone(edge);
  int16 invert_next = edge_above[outgroup];
  
  // We'll later add an edge from the now-unallocated root node to the outgroup.
  ret(invert_next, 0) = root_node;
  ret(invert_next, 1) = edge(invert_next, 0);
  
  do {
    invert_next = edge_above[edge(invert_next, 0)];
    ret(invert_next, 0) = edge(invert_next, 1);
    ret(invert_next, 1) = edge(invert_next, 0);
  } while (edge(invert_next, 0) != root_node);
  
  delete[] edge_above;
  
  // second root i.e. 16 -- 24 must be replaced with root -> outgroup.
  int16 spare_edge = (ret(root_edges[0], 0) == root_node ? 0 : 1);
  ret(invert_next, 1) = edge(root_edges[spare_edge], 1);
  ret(root_edges[spare_edge], 1) = outgroup;
  
  return preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
}

// edge must be in preorder
//  [[Rcpp::export]]
IntegerMatrix spr_moves(const IntegerMatrix edge) {
  const int16
    n_edge = edge.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    root_node = n_tip + 1
  ;
  if (n_edge < 5) return IntegerMatrix(0, 0);
  if (edge(0, 0) != root_node) throw std::invalid_argument("edge[1,] must connect root to leaf. Try Preorder(root(tree)).");
  if (edge(1, 0) != root_node) throw std::invalid_argument("edge[2,] must connect root to leaf. Try Preorder(root(tree)).");
  
  int16* prune = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* graft = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* above = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* bside = new int16[(n_edge - 1) * (n_edge - 3)];
  int16 n_moves = 0, root_daughter_2 = 0;
  const int16
    second_root_child = root_node + 1
  ;
  // Root edge first
  for (int16 i = 3; i != n_edge; i++) {
    if (edge(i, 0) == second_root_child) {
      //Rcout << "Root daughter edges are 3 and " << (1+i) << "\n";
      root_daughter_2 = i;
    } else {
      //Rcout << "\n _ Logging graft option, 1 -> "  << (i + 1) << "\n";
      prune[n_moves] = 0;
      graft[n_moves] = i;
      ++n_moves;
    }
  }
  
  for (int16 i = 0; i != n_moves; i++) {
    above[i] = -1;
    bside[i] = root_daughter_2;
  }
  
  for (int16 prune_candidate = 2; prune_candidate != n_edge; prune_candidate++) {
    const int16
      prune_parent = edge(prune_candidate, 0),
      first_prune_move = n_moves;
    ;
    int16 edge_above = 0, edge_beside = 0, i = 0;
    bool adrift = false;
    
    if (edge(1, 1) == prune_parent) edge_above = 1;
    for (i = 2; i != n_edge; i++) {
      if (edge(i, 1) == prune_parent) {
        //Rcout << "\n - Edge above broken is " << (1 + i);
        edge_above = i;
        continue;
      }
      if (i == prune_candidate) {
        if (edge(i, 1) <= n_tip) {
          ++i;
          break;
        }
        //Rcout << "\n - We're adrift! " << (1 + i);
        adrift = true;
        continue;
      }
      if (adrift) {
        if (edge(i, 0) == prune_parent) {
          //Rcout << "\n ...   Back to shore! " << (1 + i);
          break; // Now we know that all remaining edges will be potential merge sites
        }
        //Rcout << "\n ...   Still adrift! " << (1 + i);
      } else {
        if (edge(i, 0) == prune_parent) {
          //Rcout << "\n - Edge beside broken = " << (1 + i);
          edge_beside = i;
        } else {
          prune[n_moves] = prune_candidate;
          graft[n_moves] = i;
          ++n_moves;
        }
      }
    }
    if (!edge_beside) {
      //Rcout << "\n - Oo Err. Edge beside broken = " << (1 + i);
      edge_beside = i++;
    }
    if (i != n_edge + 1) while (i != n_edge) {
      prune[n_moves] = prune_candidate;
      graft[n_moves] = i;
      ++n_moves;
      ++i;
    }
    for (int16 j = first_prune_move; j != n_moves; j++) {
      above[j] = edge_above;
      bside[j] = edge_beside;
    }
  }
  
  
  IntegerMatrix ret(n_moves, 4);
  for (int16 i = n_moves; i--; ) {
    ret(i, 0) = prune[i];
    ret(i, 1) = graft[i];
    ret(i, 2) = above[i];
    ret(i, 3) = bside[i];
  }
  delete[] prune;
  delete[] graft;
  delete[] above;
  delete[] bside;
  return (ret);
}


// Assumptions: 
//  * Tree is bifurcating, in preorder; first two edges have root as parent.
//  [[Rcpp::export]]
IntegerMatrix spr (const IntegerMatrix edge,
                   const IntegerVector move) {
  const IntegerMatrix move_list = spr_moves(edge);
  const int16
    n_edge = edge.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    root_node = n_tip + 1,
    move_id = move[0] % move_list.nrow(),
    prune_edge = move_list(move_id, 0),
    graft_edge = move_list(move_id, 1),
    broken_edge_parent = edge(prune_edge, 0)
  ;
  
  if (n_edge < 5) return edge;
  if (edge(0, 0) != root_node) throw std::invalid_argument("edge[1,] must connect root to leaf. Try Preorder(root(tree)).");
  if (edge(1, 0) != root_node) throw std::invalid_argument("edge[2,] must connect root to leaf. Try Preorder(root(tree)).");
  
  IntegerMatrix ret = clone(edge);
  
  if (prune_edge) { // We are breaking a non-root edge
    const int16
      edge_above = move_list(move_id, 2),
      edge_beside = move_list(move_id, 3)
    ;
    
    ret(edge_beside, 0) = edge(edge_above, 0);
    ret(edge_above, 0) = edge(graft_edge, 0);
    ret(graft_edge, 0) = broken_edge_parent;
  } else { // We are breaking the root edge
    ret(2, 0) = broken_edge_parent;
    ret(move_list(move_id, 3), 0) = broken_edge_parent;
    
    //child [brokenEdgeSister] <- child[mergeEdge]
    ret(1, 1) = edge(graft_edge, 1);
    //parent[brokenEdge | brokenEdgeSister] <- spareNode
    const int spare_node = edge(1, 1);
    ret(0, 0) = spare_node;
    ret(1, 0) = spare_node;
    // child[mergeEdge] <- spareNode
    ret(graft_edge, 1) = spare_node;
  }
  ret = preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
  return root_on_node(ret, 1);
}
