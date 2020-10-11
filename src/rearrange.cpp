#include <Rcpp.h>
// [ [Rcpp::depends(TreeTools)]]
#include <TreeTools.h>
#include <memory> /* for unique_ptr */
using namespace std;
using namespace Rcpp;

typedef int_fast16_t int16;
const int16 UNDEFINED = -1;

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
  
  return TreeTools::preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
}

// edge must be in preorder
//  [[Rcpp::export]]
IntegerMatrix spr_moves(const IntegerMatrix edge) {
  const int16
    n_edge = edge.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    root_node = n_tip + 1,
    second_root_child = root_node + 1
  ;
  if (n_edge < 5) return IntegerMatrix(0, 0);
  if (edge(0, 0) != root_node) throw std::invalid_argument("edge[1,] must connect root to leaf. Try Preorder(root(tree)).");
  if (edge(1, 0) != root_node) throw std::invalid_argument("edge[2,] must connect root to leaf. Try Preorder(root(tree)).");
  
  int16* prune = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* graft = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* above = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* bside = new int16[(n_edge - 1) * (n_edge - 3)];
  int16 n_moves = 0, root_daughter_2 = 0;
  
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
  
  if (n_edge < 5) throw std::invalid_argument("No SPR rearrangements possible on a tree with < 5 edges");
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
  ret = TreeTools::preorder_edges_and_nodes(ret(_, 0), ret(_, 1));
  return TreeTools::root_binary(ret, 1);
}

// Assumptions: 
//  * Tree is bifurcating, in preorder; first two edges have root as parent.
//  [[Rcpp::export]]
IntegerMatrix tbr_moves(const IntegerMatrix edge) {
  const int16
    n_edge = edge.nrow(),
    n_node = n_edge / 2,
    n_tip = n_node + 1,
    root_node = n_tip + 1,
    second_root_child = root_node + 1
  ;
  if (n_edge < 5) throw std::invalid_argument("No TBR rearrangements possible on a tree with < 5 edges");
  if (edge(0, 0) != root_node) throw std::invalid_argument("edge[1,] must connect root to leaf. Try Preorder(root(tree)).");
  if (edge(1, 0) != root_node) throw std::invalid_argument("edge[2,] must connect root to leaf. Try Preorder(root(tree)).");
  
  std::unique_ptr<int16[]> n_edges_above = std::make_unique<int16[]>(n_edge);
  std::unique_ptr<int16[]> probibited_parent = std::make_unique<int16[]>(n_edge);
  std::unique_ptr<int16[]> probibited_sibling = std::make_unique<int16[]>(n_edge);
  
  int16* prune = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* graft = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* above = new int16[(n_edge - 1) * (n_edge - 3)];
  int16* bside = new int16[(n_edge - 1) * (n_edge - 3)];
  int16 n_moves = 0, root_daughter_2 = 0;
  
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
  
  for (int16 bisect = 1; bisect != n_edge; bisect++) {
    
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

//  [[Rcpp::export]]
IntegerMatrix tbr (const IntegerMatrix edge,
                   const IntegerVector move) {
  const IntegerMatrix move_list = tbr_moves(edge);
  // Actually do TBR move
  return IntegerMatrix(0, 0);
}

inline void set_child(unique_ptr<int16[]> &side, const int16 parent, 
               const int16 value, const int16 n_tip) {
  side[parent - 1 - n_tip] = value;
}

inline int16 get_child(unique_ptr<int16[]> &side, const int16 parent, const int16 n_tip) {
  return side[parent - 1 - n_tip];
}

inline int16 count_children(unique_ptr<int16[]> &n_children, const int16 vert) {
  return n_children[vert - 1];
}

inline void add_children(unique_ptr<int16[]> &n_children,
                         const int16 parent, const int16 child) {
  n_children[parent - 1] += n_children[child - 1];
}

inline int16 edge_above(unique_ptr<int16[]> &parent_edge, const int16 vert) {
  return parent_edge[vert - 1];
}

// Assumptions: 
//  * Tree is bifurcating, in preorder; first two edges have root as parent.
//  [[Rcpp::export]]
List all_tbr (const IntegerMatrix edge,
              const IntegerVector break_order) {
  const int16
    n_edge = edge.nrow(),
    n_internal = n_edge / 2,
    n_tip = n_internal + 1,
    n_vert = n_internal + n_tip,
    root_node = n_tip + 1,
    second_root_child = root_node + 1
  ;
  if (n_edge < 5) throw std::invalid_argument("No TBR rearrangements possible on a tree with < 5 edges");
  if (edge(0, 0) != root_node) throw std::invalid_argument("edge[1,] must connect root to leaf. Try Preorder(root(tree)).");
  if (edge(1, 0) != root_node) throw std::invalid_argument("edge[2,] must connect root to leaf. Try Preorder(root(tree)).");
  
  IntegerVector break_seq;
  if (break_order.length()) {
    break_seq = clone(break_order);
  } else {
    IntegerVector tmp (n_edge - 2);
    break_seq = tmp;
    for (int16 i = n_edge - 2; i--; ) {
      break_seq[i] = i + 3;
    }
  }

  unique_ptr<int16[]> n_children = make_unique<int16[]>(n_vert);
  unique_ptr<int16[]> left_node = make_unique<int16[]>(n_internal);
  unique_ptr<int16[]> right_node = make_unique<int16[]>(n_internal);
  unique_ptr<int16[]> left_edge = make_unique<int16[]>(n_internal);
  unique_ptr<int16[]> right_edge = make_unique<int16[]>(n_internal);
  unique_ptr<int16[]> parent_edge = make_unique<int16[]>(n_vert);
  for (int16 i = n_tip; i--; ) {
    n_children[i] = 1;
  }
  
  for (int16 i = n_edge; i--; ) {
    const int parent = edge(i, 0);
    const int child = edge(i, 1);
    parent_edge[child - 1] = i;
    if (get_child(left_node, parent, n_tip)) {
      set_child(right_node, parent, child, n_tip);
      set_child(right_edge, parent, i, n_tip);
    } else {
      set_child(left_node, parent, child, n_tip);
      set_child(left_edge, parent, i, n_tip);
    }
  }
  
  
  List ret = List::create();
  
  // Let's go.
  for (int16 i = break_seq.length(); i--; ) {
    IntegerMatrix two_bits = clone(edge);
    const int16
      break_edge = break_seq[i] - 1,
      break_parent = edge(break_edge, 0),
      break_child = edge(break_edge, 1),
      spare_node = break_parent,
      fragment_leaves = count_children(n_children, break_child),
      fragment_edges = fragment_leaves + fragment_leaves - 1,
      fragment_min_edge = break_edge,
      fragment_max_edge = break_edge + fragment_edges - 1,
      base_leaves = n_tip - fragment_leaves,
      base_edges = n_edge - fragment_edges - 1 // -1 for the broken edge
    ;
    const bool broken_on_left = get_child(left_edge, break_parent, n_tip) == break_edge;
    const int16
      spare_edge = broken_on_left ?
        get_child(right_edge, break_parent, n_tip) :
        get_child(left_edge, break_parent, n_tip)
    ;
    
    two_bits(edge_above(parent_edge, break_parent), 1) = broken_on_left ?
      get_child(right_node, break_parent, n_tip) :
      get_child(left_node, break_parent, n_tip);
    if (fragment_leaves == 1) {
      for (int16 graft_edge = n_edge; graft_edge--; ) {
        if (graft_edge == fragment_max_edge) {
          // Remember: graft_location will be decremented after continue
          graft_edge = fragment_min_edge;
          continue; 
        } else if (broken_on_left && graft_edge == get_child(right_node, break_parent, n_tip)) {
          // Remember: graft_location will be decremented after continue
          graft_edge = parent_edge[break_parent];
          continue;
        } else if (graft_edge == parent_edge[break_parent]) {
          continue;
        }
        IntegerMatrix new_tree = clone(two_bits);
        
        new_tree(spare_edge, 1) = two_bits(graft_edge, 1);
        new_tree(graft_edge, 1) = spare_node;
        new_tree(break_edge, 0) = spare_node;
        new_tree = TreeTools::postorder_edges(new_tree);
        ret.push_back(new_tree);
      }
    } else {
      for (int16 loose_root = fragment_edges; loose_root--; ) {
        for (int16 graft_location = base_edges; graft_location--; ) {
          
        }
      }
    }
    
    
  }
  return ret;
}