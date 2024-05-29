#include <Rcpp.h>
#include <TreeTools/renumber_tree.h>
using namespace Rcpp;

inline bool multiple_bits(int n) {
  bool found_one = false;
  while (n) {
    if (n & 1) {
      if (found_one) return true;
      found_one = true;
    }
    n >>= 1;
  }
  return false;
}

inline int first_bit(int n) {
  int ret = 1;
  while(n) {
    if (n & 1) return ret;
    n >>= 1;
    ret <<= 1;
  }
  return 0;
}

// rows in states correspond to sequence of tree labels, columns to characters
// each entry in `states` is a binary map of tokens present.
// [[Rcpp::export]]
List character_regions(const List tree, const IntegerMatrix states,
                       const LogicalVector acctrans) {
  IntegerMatrix edge = tree["edge"];
  IntegerVector parent = edge(_, 0);
  IntegerVector child = edge(_, 1);
  const IntegerVector postorder = TreeTools::postorder_order(edge);

  const int n_node = tree["Nnode"];
  const int n_tip = CharacterVector(tree["tip.label"]).length();
  const int n_vert = n_node + n_tip;
  
  const bool acctran = acctrans[0];
  /*
  edge <- tree[["edge"]][postorder, ]
  parent <- edge[, 1]
  child <- edge[, 2]*/
  
  IntegerVector left(n_vert + 1);
  IntegerVector right(n_vert + 1);
  IntegerVector parent_of(n_vert + 1);
  IntegerVector postorder_nodes(n_node);
  
  // Read structure of tree
  for (int i = postorder.length(); i--; ) {
    int e = postorder[i];
    int pa = parent[e];
    int ch = child[e];
    if (i % 2) {
      postorder_nodes[i / 2] = pa;
    }
    parent_of[ch] = pa;
    if (right[pa]) {
      left[pa] = ch;
    } else {
      right[pa] = ch;
    }
  }
  const int root_node = postorder_nodes[n_node - 1];
  parent_of[root_node] = root_node;
  const int n_patterns = states.ncol();
  
  IntegerMatrix state(n_vert, n_patterns);
  for (int i = n_vert; i--; ) {
    for (int j = n_patterns; j--; ) {
      state(i, j) = states(i, j);
    }
  }
  
  // Fitch optimization
  for (int i = 0; i != n_node; i++) {
    const int node = postorder_nodes[i];
    const IntegerVector left_state = state(left[node], _);
    const IntegerVector right_state = state(right[node], _);
    const IntegerVector common = left_state & right_state;
    state(node, _) = ifelse(common, common, left_state | right_state);
  }
  
  for (int i = n_patterns; i--; ) {
    state[root_node, i] = first_bit(state[root_node, i]);
  }
  
  for (int i = n_node; i--; ) {
    const int node = postorder_nodes[i];
    const int node_state = state(node, _);
    const int anc_state = state(parent_of[node], _);
    const int left_state = state(left[node], _);
    const int right_state = state(right[node], _);
    const int left_pipe_right = left_state | right_state;
    const int inherited = node_state & anc_state;
    state(node, _) = ifelse(
      inherited == anc_state,
      inherited,
      acctran ? node_state : ifelse(
        left_state & right_state,
        ifelse(anc_state & left_pipe_right, anc_state, left_pipe_right),
        anc_state
      )
    );
  }
  
  const int n_pattern = states.ncol();
  List ret(n_pattern);
  IntegerMatrix true_label(n_node, n_pattern);
  for (int pat = n_pattern; pat--; ) {
    int last_label = 1;
    true_label[root_node] = last_label;
    
    for (int i = n_node; i--; ) {
      const int node = postorder_nodes[i];
      const int anc = parent_of[node];
      const int final_state = state(node, pat);
      const int anc_state = state(anc, pat);
      true_label(node, pat) = final_state == anc_state ? 
        true_label(anc, pat) : ++last_label;
    }
    
    IntegerMatrix n_with_label(n_node);
    for (int i = n_tip; i--; ) {
      const int tip = i + 1; // C++ -> R
      tip_state = state(tip, pat);
      // Only consider unambiguous leaf labels
      if (!multiple_bits(tip_state)) {
        const int anc = parent_of[tip];
        const int anc_state = state(nodeParent, pat);
        const int this_label = anc_state == tip_state ? 
          true_label(anc, pat) : ++last_label;
        ++(n_with_label[this_label]);
      }
    }
    ret[pat] = IntegerVector(n_with_label.begin(), n_with_label.begin() + last_label);
  }
  return ret;
}
