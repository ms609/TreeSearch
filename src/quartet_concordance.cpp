#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List quartet_concordance(const LogicalMatrix splits, const IntegerMatrix characters) {
  const int n_splits = splits.ncol();
  const int n_chars = characters.ncol();
  const int n_taxa = splits.nrow();
  
  NumericMatrix concordant(n_splits, n_chars);
  NumericMatrix decisive(n_splits, n_chars);

  // Pre-allocate vectors to avoid repeated allocations in the inner loop
  // We'll use these to store the counts for each state.
  // 32 is usually plenty for states; if you expect more, use a larger constant.
  std::vector<int> n0(32), n1(32);
  std::vector<int> active_states;
  active_states.reserve(32);

  for (int c = 0; c < n_chars; ++c) {
    std::vector<int> char_col(n_taxa);
    for (int t = 0; t < n_taxa; ++t) char_col[t] = characters(t, c);

    for (int s = 0; s < n_splits; ++s) {
      active_states.clear();

      for (int t = 0; t < n_taxa; ++t) {
        int state = char_col[t];
        if (IntegerVector::is_na(state)) continue;
        
        if (state >= (int)n0.size()) {
          n0.resize(state + 1, 0);
          n1.resize(state + 1, 0);
        }
        
        if (n0[state] == 0 && n1[state] == 0) {
          active_states.push_back(state);
        }

        if (splits(t, s)) {
          n1[state]++;
        } else {
          n0[state]++;
        }
      }

      // 2. Combinatorics logic
      double total_choose_n1 = 0;
      double sum_prod = 0;
      double sum_sq_prod = 0;

      for (int state : active_states) {
        double cnt0 = n0[state];
        double cnt1 = n1[state];
        
        total_choose_n1 += (cnt1 * (cnt1 - 1.0) / 2.0);
        double p = cnt0 * cnt1;
        sum_prod += p;
        sum_sq_prod += p * p;
      }

      double conc_val = 0;
      for (int state : active_states) {
        double cnt0 = n0[state];
        double cnt1 = n1[state];
        double cn0 = cnt0 * (cnt0 - 1.0) / 2.0;
        double cn1 = cnt1 * (cnt1 - 1.0) / 2.0;
        conc_val += cn0 * (total_choose_n1 - cn1);
      }

      concordant(s, c) = conc_val;
      decisive(s, c) = conc_val + (0.5 * (sum_prod * sum_prod - sum_sq_prod));
      
      // Cleanup using the tracked active states
      for (int state : active_states) {
        n0[state] = 0;
        n1[state] = 0;
      }
    }
  }

  return List::create(_["concordant"] = concordant, _["decisive"] = decisive);
}
