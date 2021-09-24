#include <Rcpp.h>
#include <TreeTools.h>
using namespace std;
using namespace Rcpp;

inline void add_children(unique_ptr<int16[]> &n_children,
                         const int16 parent, const int16 child) {
  n_children[parent - 1] += n_children[child - 1];
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

inline int16 edge_above(const int16 vert, unique_ptr<int16[]> &parent_edge) {
  return parent_edge[vert - 1];
}

inline IntegerMatrix fuse(const IntegerMatrix& tree_bits,
                          const int16* graft_edge, const int16* break_edge, 
                          const int16* spare_edge, const int16* spare_node) {
  IntegerMatrix new_tree = clone(tree_bits);
  new_tree(*spare_edge, 1) = tree_bits(*graft_edge, 1);
  new_tree(*graft_edge, 1) = *spare_node;
  new_tree(*break_edge, 0) = *spare_node;
  return TreeTools::preorder_edges_and_nodes(new_tree(_, 0), new_tree(_, 1));
}

//  [[Rcpp::export]]
List asan_error (const IntegerMatrix edge, const IntegerVector break_order) {
  const int
    n_edge = edge.nrow(),
    n_internal = n_edge / 2,
    n_tip = n_internal + 1,
    n_vert = n_internal + n_tip,
    root_node = n_tip + 1
  ;
  
  if (n_edge < 5000) {
    Rf_error("Oh dear.");
  }
  
  IntegerVector break_seq;
  if (break_order.length()) {
    break_seq = clone(break_order);
  } else {
    IntegerVector tmp (n_edge - 1);
    break_seq = tmp;
    for (int16 i = n_edge - 1; i--; ) {
      break_seq[i] = i + 2;
    }
  }
  
  unique_ptr<int16[]> n_children = make_unique<int16[]>(n_vert);
  // if both internal, left_node < right_node
  unique_ptr<int16[]> left_node = make_unique<int16[]>(n_internal);
  unique_ptr<int16[]> right_node = make_unique<int16[]>(n_internal);
  // left_edge is plotted on top with ape::plot.phylo
  unique_ptr<int16[]> left_edge = make_unique<int16[]>(n_internal);
  unique_ptr<int16[]> right_edge = make_unique<int16[]>(n_internal);
  unique_ptr<int16[]> parent_edge = make_unique<int16[]>(n_vert);
  for (int16 i = n_tip; i--; ) {
    n_children[i] = 1;
  }
  
  for (int16 i = n_edge; i--; ) {
    const int parent = edge(i, 0);
    const int child = edge(i, 1);
    add_children(n_children, parent, child);
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
    if (break_seq[i] > n_edge) {
      Rf_warning("Ignoring SPR break locations that exceed number of edges in tree.\n");
      continue;
    }
    if (break_seq[i] < 2) {
      Rf_warning("Ignoring break locations < 2");
      continue;
    }
    const int16
      break_edge = break_seq[i] - 1,
        break_parent = edge(break_edge, 0),
        break_child = edge(break_edge, 1),
        spare_node = break_parent,
        
        fragment_root = break_child,
        fragment_leaves = count_children(n_children, break_child),
        fragment_edges = fragment_leaves + fragment_leaves - 1,
        fragment_min_edge = break_edge,
        fragment_max_edge = break_edge + fragment_edges - 1
    ;
    const bool broken_on_left = get_child(left_edge, break_parent, n_tip) == break_edge;
    const int16
      spare_edge = broken_on_left ?
    get_child(right_edge, break_parent, n_tip) :
      get_child(left_edge, break_parent, n_tip)
      ;
    
    two_bits(edge_above(break_parent, parent_edge), 1) = broken_on_left ?
    get_child(right_node, break_parent, n_tip) :
      get_child(left_node, break_parent, n_tip);
    if (break_edge == 1) {
      const int16
      fragment_base_right = 2,
        fragment_base_left = get_child(left_edge, fragment_root, n_tip);
      ;
      
      for (int16 insertion_point = fragment_min_edge + 2;
           insertion_point != fragment_max_edge + 1; insertion_point++) {
        if (insertion_point == fragment_base_left) {
          continue;
        }
        
        int16 invert_next = insertion_point;
        IntegerMatrix rerooted = clone(two_bits);
        
        rerooted(invert_next, 0) = break_child; // Borrow fragment-root node id
        rerooted(invert_next, 1) = two_bits(invert_next, 0);
        
        do {
          invert_next = edge_above(two_bits(invert_next, 0), parent_edge);
          rerooted(invert_next, 0) = two_bits(invert_next, 1);
          rerooted(invert_next, 1) = two_bits(invert_next, 0);
        } while (two_bits(invert_next, 0) != fragment_root);
        
        const bool new_root_on_right = invert_next == fragment_base_right;
        const int16 repurposed_edge = new_root_on_right ?
        fragment_base_left : 
          fragment_base_right;
        rerooted(invert_next, 1) = two_bits(repurposed_edge, 1);
        rerooted(repurposed_edge, 1) = two_bits(insertion_point, 1);
        rerooted = TreeTools::preorder_edges_and_nodes(rerooted(_, 0), rerooted(_, 1));
        ret.push_back(rerooted);
      }
    } else {
      for (int16 graft_edge = n_edge - 1; graft_edge; graft_edge--) {
        if (graft_edge == fragment_max_edge) {
          graft_edge = fragment_min_edge;
          continue; 
        } else if (broken_on_left && graft_edge == get_child(right_edge, break_parent, n_tip)) {
          graft_edge = edge_above(break_parent, parent_edge);
          continue;
        } else if (graft_edge == spare_edge) {
          continue;
        } else if (graft_edge == edge_above(break_parent, parent_edge)) {
          continue;
        }
        ret.push_back(fuse(two_bits, &graft_edge, &break_edge, &spare_edge,
                           &spare_node));
        if (graft_edge < 0) break; // TODO REMOVE
      }
    }
  }
  return ret;
}