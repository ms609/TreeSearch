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
  IntegerMatrix edge = tree["edge"];
  IntegerVector parent = edge(_, 0);
  IntegerVector child = edge(_, 1);
  const IntegerVector postorder = TreeTools::postorder_order(edge);

  const int n_node = tree["Nnode"];
  const int n_tip = CharacterVector(tree["tip.label"]).length();
  const int n_vert = n_node + n_tip;
  const int n_edge = postorder.length();
  
  const bool acctran = acctrans[0];
  /*
  edge <- tree[["edge"]][postorder, ]
  parent <- edge[, 1]
  child <- edge[, 2]*/
  
  // These will hold the properties of each node based on its R label, i.e.
  // their 0th entry will be unused
  IntegerVector left(n_vert + 1);
  IntegerVector right(n_vert + 1);
  IntegerVector parent_of(n_vert + 1);
  IntegerVector postorder_nodes(n_node);
  
  int next_node = 0;
  
  // Read structure of tree
  for (int i = 0; i != n_edge; i++) {
    int e = postorder[i] - 1;
    Rcout << "Edge " << (e + 1) <<": ";
    int pa = parent[e];
    int ch = child[e];
    Rcout << pa << " - " << ch <<"\n";
    if (ch > n_tip) {
      Rcout << "  Postorder_nodes[" << next_node << "] = " << ch << ".\n";
      postorder_nodes[next_node++] = ch;
    }
    parent_of[ch] = pa;
    if (right[pa]) {
      left[pa] = ch;
    } else {
      right[pa] = ch;
    }
  }
  const int root_node = parent[postorder[n_edge - 1] - 1];
  Rcout << "Root node = " << root_node <<  ".\n";
  parent_of[root_node] = root_node;
    Rcout << "  Postorder_nodes[" << next_node << "] = " << root_node << ".\n";
  postorder_nodes[next_node] = root_node;
  
  const int n_patterns = states.ncol();
  // state uses C++ numbering, i.e. subtract 1 to get the right row
  IntegerMatrix state(n_vert, n_patterns);
  for (int i = n_vert; i--; ) {
    for (int j = n_patterns; j--; ) {
      state(i, j) = states(i, j);
    }
  }
  
  Rcout << "Loaded states.";
  
  // Fitch optimization
  for (int i = 0; i != n_node; i++) {
    const int node = postorder_nodes[i];
    const IntegerVector left_state = state(left[node] - 1, _);
    const IntegerVector right_state = state(right[node] - 1, _);
    const IntegerVector common = AND(left_state, right_state);
    const IntegerVector left_pipe_right = OR(left_state, right_state);
    state(node - 1, _) = ifelse(common != 0, common, left_pipe_right);
  }
  
  Rcout << "   * 120 \n";
  
  for (int i = n_patterns; i--; ) {
    state(root_node - 1, i) = first_bit(state(root_node - 1, i));
  }
  
  Rcout << "   * 126 \n";
    
  for (int i = n_node; i--; ) {
    const int node = postorder_nodes[i];
    const IntegerVector node_state = state(node - 1, _);
    const IntegerVector anc_state = state(parent_of[node] - 1, _);
    const IntegerVector left_state = state(left[node] - 1, _);
    const IntegerVector right_state = state(right[node] - 1, _);
    const IntegerVector left_pipe_right = OR(left_state, right_state);
    const IntegerVector inherited = AND(node_state, anc_state);
    state(node - 1, _) = ifelse(
      inherited == anc_state,
      inherited,
      acctran ? node_state : ifelse(
        AND(left_state, right_state) != 0,
        ifelse(AND(anc_state, left_pipe_right) != 0, anc_state, left_pipe_right),
        anc_state
      )
    );
  }
  
  Rcout << "   * 147 \n";
  const int n_pattern = states.ncol();
  List ret(n_pattern);
  // true_label uses C++ node numbering so we need a -1
  IntegerMatrix true_label(n_node, n_pattern);
  for (int pat = n_pattern; pat--; ) {
    int last_label = 1;
    true_label[root_node] = last_label;
    
    for (int i = n_node; i--; ) {
      const int node = postorder_nodes[i];
      const int anc = parent_of[node];
      const int final_state = state(node - 1, pat);
      const int anc_state = state(anc - 1, pat);
      Rcout << "- Node " << node << " (state " << final_state << "), anc = "
            << anc << " (" << anc_state << ").\n";
      true_label(node - 1, pat) = final_state == anc_state ? 
        true_label(anc - 1, pat) : ++last_label;
      Rcout << "- Labelling true label(" << (node - 1) << ", " << pat <<") = "
            << (final_state == anc_state ? "T" : "F") <<  true_label(node - 1, pat)
            << "\n";
    }
    
  Rcout << "   * 165 \n\n";
    IntegerVector n_with_label(n_vert);
    for (int i = n_tip; i--; ) {
      const int tip = i + 1; // C++ -> R
      const int tip_state = state(tip - 1, pat);
      // Only consider unambiguous leaf labels
      if (!multiple_bits(tip_state)) {
        const int anc = parent_of[tip];
        const int anc_state = state(anc - 1, pat);
        const int this_label = anc_state == tip_state ? 
          true_label(anc - 1, pat) : ++last_label;
        Rcout << " This label: " << this_label << "; last = " << last_label <<"\n";
        ++(n_with_label[this_label]);
      }
  Rcout << "   * 178 \n";
    }
  Rcout << "   * 180 \n";
    Rcout << "Last label = " << last_label << ".\n";
    ret[pat] = IntegerVector(n_with_label.begin(), n_with_label.begin() + last_label);
  Rcout << "   * 181 \n";
  }
  return ret;
}
