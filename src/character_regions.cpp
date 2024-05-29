#include <Rcpp.h>
#include <TreeTools/renumber_tree.h>
using namespace Rcpp;


// rows in dataset must correspond to sequence of tree labels
// [[Rcpp::export]]
List character_regions(const List tree, const List dataset,
                       const LogicalVector acctrans) {
  IntegerMatrix edge = tree["edge"];
  IntegerVector parent = edge(_, 0);
  IntegerVector child = edge(_, 1);
  IntegerVector postorder = TreeTools::postorder_order(edge);

  int n_node = tree["Nnode"];
  int n_tip = tree["tip.label"].length();
  /*
  edge <- tree[["edge"]][postorder, ]
  parent <- edge[, 1]
  child <- edge[, 2]*/
  
  int left = n_node + n_tip;
  int right = left;
  IntegerVector parent_of(left + 1);
  IntegerVector postorder_nodes(n_node);
  
  for (int i = postorder.length(); i--; ) {
    int e = postorder[i];
    int pa = parent[e];
    int child = child[e];
    if (i %% 2) {
      postorder_nodes[i / 2] = pa;
    }
    parent_of[ch] = pa;
    if (right[pa]) {
      left[pa] = ch;
    } else {
      right[pa] = ch;
    }
  }
  int root_node = postorder_nodes[n_node - 1];
  parent_of[root_node] = root_node;
}
