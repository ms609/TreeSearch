#include <Rcpp.h>
#include <cassert>
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

inline IntegerVector AND(const IntegerVector a, const IntegerVector b) {
  const int n = a.length();
  assert(b.length() == n);
  IntegerVector ret(n);
  for (int i = n; i--; ) {
    ret[i] = a[i] & b[i];
  }
  return ret;
}

inline IntegerVector OR(const IntegerVector a, const IntegerVector b) {
  const int n = a.length();
  assert(b.length() == n);
  IntegerVector ret(n);
  for (int i = n; i--; ) {
    ret[i] = a[i] | b[i];
  }
  return ret;
}


// rows in states correspond to sequence of tree labels, columns to characters
// each entry in `states` is a binary map of tokens present.
// [[Rcpp::export]]
List character_regions(const List tree, const IntegerMatrix states,
                       const LogicalVector acctrans) {
  // Retain R numbering - needed for postorder_order
  const IntegerMatrix edge = tree["edge"];
  
  // We'll number everything else from zero as soon as it lands in C++
  const IntegerVector parent = edge(_, 0) - 1;
  const IntegerVector child = edge(_, 1) - 1;
  const IntegerVector postorder = TreeTools::postorder_order(edge) - 1;

  const int n_node = tree["Nnode"];
  const int n_tip = CharacterVector(tree["tip.label"]).length();
  const int n_vert = n_node + n_tip;
  const int n_edge = postorder.length();
  
  const bool acctran = acctrans[0];
  
  IntegerVector left(n_vert);
  IntegerVector right(n_vert);
  IntegerVector parent_of(n_vert);
  IntegerVector postorder_nodes(n_node);
  
  int next_node = 0;
  
  // Read structure of tree
  for (int i = 0; i != n_edge; i++) {
    int e = postorder[i];
    // Rcout << "Edge " << e <<": ";
    int pa = parent[e];
    int ch = child[e];
    // Rcout << pa << " - " << ch <<"\n";
    if (ch > n_tip) {
      postorder_nodes[next_node++] = ch;
    }
    parent_of[ch] = pa;
    if (right[pa]) {
      left[pa] = ch;
    } else {
      right[pa] = ch;
    }
  }
  const int last_edge = n_edge - 1;
  const int root_node = parent[postorder[last_edge]];
  // Rcout << "Root node = " << root_node <<  ".\n";
  parent_of[root_node] = root_node;
  postorder_nodes[next_node] = root_node;
  
  const int n_pattern = states.ncol();
  IntegerMatrix state(n_vert, n_pattern);
  for (int i = n_vert; i--; ) {
    for (int j = n_pattern; j--; ) {
      state(i, j) = states(i, j);
    }
  }
  
  // Fitch optimization
  for (int i = 0; i != n_node; i++) {
    const int node = postorder_nodes[i];
    const IntegerVector left_state = state(left[node], _);
    const IntegerVector right_state = state(right[node], _);
    const IntegerVector common = AND(left_state, right_state);
    for (int pat = n_pattern; pat--; ) {
      state(node, pat) = common[pat] ?
        common[pat] : left_state[pat] | right_state[pat];
    }
  }
  
  for (int i = n_pattern; i--; ) {
    state(root_node, i) = first_bit(state(root_node, i));
  }
    
  for (int i = n_node; i--; ) {
    const int node = postorder_nodes[i];
    const IntegerVector node_state = state(node, _);
    const IntegerVector anc_state = state(parent_of[node], _);
    const IntegerVector left_state = state(left[node], _);
    const IntegerVector right_state = state(right[node], _);
    const IntegerVector inherited = AND(node_state, anc_state);
    
    for (int pat = n_pattern; pat--; ) {
      const IntegerVector left_pipe_right = OR(left_state, right_state);
      if (inherited[pat] == anc_state[pat]) {
        state(node, pat) = inherited[pat];
      } else if (acctran) {
        state(node, pat) = node_state[pat];
      } else if (left_state[pat] & right_state[pat]) {
        const int left_pipe_right = left_state[pat] | right_state[pat];
        state(node, pat) = anc_state[pat] & left_pipe_right ?
          anc_state[pat] : left_pipe_right;
      } else {
        state(node, pat) = anc_state[pat];
      }
    }
  }
  
  List ret(n_pattern);
  IntegerMatrix true_label(n_node, n_pattern);
  for (int pat = n_pattern; pat--; ) {
    int last_label = 0;
    true_label(root_node, pat) = last_label;
    
    for (int i = n_node; i--; ) {
      const int node = postorder_nodes[i];
      const int anc = parent_of[node];
      const int final_state = state(node, pat);
      const int anc_state = state(anc, pat);
      // Rcout << "- Node " << node << " (state " << final_state << "), anc = "
      //       << anc << " (" << anc_state << ").\n";
      true_label(node, pat) = final_state == anc_state ? 
        true_label(anc, pat) : ++last_label;
      // Rcout << " - Labelling true label(" << node << ", " << pat <<") = "
      //       << (final_state == anc_state ? "T" : "F") <<  true_label(node, pat)
      //       << "\n";
    }
    
    IntegerVector n_with_label(n_vert);
    for (int i = n_tip; i--; ) {
      const int tip = i;
      const int tip_state = state(tip, pat);
      // Only consider unambiguous leaf labels
      if (!multiple_bits(tip_state)) {
        const int anc = parent_of[tip];
        const int anc_state = state(anc, pat);
        const int this_label = anc_state == tip_state ? 
          true_label(anc, pat) : ++last_label;
        // Rcout << " This label: " << this_label << "; last = " << last_label <<"\n";
        ++n_with_label[this_label];
      }
    }
    // Rcout << "Last label = " << last_label << ".\n";
    ret[pat] = IntegerVector(n_with_label.begin(), n_with_label.begin() + last_label + 1);
  }
  return ret;
}
