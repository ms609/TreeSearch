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

  // Pre-allocate buffers
  std::vector<int> n0(32, 0), n1(32, 0);
  std::vector<int> active_states;
  active_states.reserve(32);

  for (int c = 0; c < n_chars; ++c) {
    // Cache character column and find max state in one pass
    std::vector<int> char_col(n_taxa);
    int max_state = 0;
    for (int t = 0; t < n_taxa; ++t) {
      int state = characters(t, c);
      char_col[t] = state;
      if (!IntegerVector::is_na(state) && state > max_state) max_state = state;
    }
    // Hoist resize outside split loop: only reallocate when a new character
    // has states beyond the current buffer capacity.
    if (max_state >= (int)n0.size()) {
      n0.resize(max_state + 1, 0);
      n1.resize(max_state + 1, 0);
    }

    for (int s = 0; s < n_splits; ++s) {
      active_states.clear();

      for (int t = 0; t < n_taxa; ++t) {
        int state = char_col[t];
        if (IntegerVector::is_na(state)) continue;

        if (n0[state] == 0 && n1[state] == 0) {
          active_states.push_back(state);
        }

        if (splits(t, s)) {
          n1[state]++;
        } else {
          n0[state]++;
        }
      }

      if (active_states.size() >= 2) {
        double total_choose_n1 = 0;
        double sum_prod = 0;
        double sum_sq_prod = 0;

        for (int state : active_states) {
          double cnt1 = n1[state];
          total_choose_n1 += (cnt1 * (cnt1 - 1.0) / 2.0);
          
          double p = (double)n0[state] * cnt1;
          sum_prod += p;
          sum_sq_prod += p * p;
        }

        double conc_val = 0;
        for (int state : active_states) {
          double cn0 = (double)n0[state] * (n0[state] - 1.0) / 2.0;
          double cn1 = (double)n1[state] * (n1[state] - 1.0) / 2.0;
          conc_val += cn0 * (total_choose_n1 - cn1);
        }

        concordant(s, c) = conc_val;
        decisive(s, c) = conc_val + (0.5 * (sum_prod * sum_prod - sum_sq_prod));
      }

      // Reset buffers for next split
      for (int state : active_states) {
        n0[state] = 0;
        n1[state] = 0;
      }
    }
  }

  return List::create(_["concordant"] = concordant, _["decisive"] = decisive);
}
