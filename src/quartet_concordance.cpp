#include <Rcpp.h>
#include <map>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List quartet_concordance(const LogicalMatrix splits, const IntegerMatrix characters) {
  int n_splits = splits.ncol();
  int n_chars = characters.ncol();
  int n_taxa = splits.nrow();
  
  // Storage for results
  // Using NumericMatrix to safely hold potentially large counts (double precision)
  NumericMatrix concordant(n_splits, n_chars);
  NumericMatrix decisive(n_splits, n_chars);
  
  // Loop over every Split (s) and every Character (c)
  for (int s = 0; s < n_splits; ++s) {
    for (int c = 0; c < n_chars; ++c) {
      
      // 1. Build Contingency Table manually
      // Map: Character State (int) -> pair<count_in_split_FALSE, count_in_split_TRUE>
      // We use a map to handle arbitrary state integers (0, 1, 2, or 10, 20...).
      std::map<int, std::pair<int, int>> counts; 
      
      for (int t = 0; t < n_taxa; ++t) {
        int char_state = characters(t, c);
        
        // Skip NAs (ambiguous characters)
        if (IntegerVector::is_na(char_state)) continue;
        
        // Check split membership (FALSE=0, TRUE=1)
        if (splits(t, s)) {
          counts[char_state].second++;
        } else {
          counts[char_state].first++;
        }
      }
      
      // If fewer than 2 states are present, no quartets can be decisive
      if (counts.size() < 2) {
        concordant(s, c) = 0;
        decisive(s, c) = 0;
        continue;
      }
      
      // 2. Calculate Combinatorics
      // We need to sum over combinations without O(N^2) loops.
      
      // Variables for Concordance: sum( choose(n0i, 2) * sum(choose(n1j, 2)) ) where j != i
      // This is equivalent to: sum( choose(n0i, 2) * (Total_Choose_N1 - choose(n1i, 2)) )
      double total_choose_n1 = 0;
      std::vector<double> choose_n0_vec;
      std::vector<double> choose_n1_vec;
      choose_n0_vec.reserve(counts.size());
      choose_n1_vec.reserve(counts.size());
      
      // Variables for Discordance: sum( n0i*n1i * n0j*n1j ) for i < j
      // This is equivalent to: 0.5 * ( (sum P)^2 - sum(P^2) ) where P_k = n0k * n1k
      double sum_prod = 0;
      double sum_sq_prod = 0;
      
      for (auto const& [state, cnt] : counts) {
        double n0 = cnt.first;
        double n1 = cnt.second;
        
        // Concordance Pre-calcs
        double cn0 = n0 * (n0 - 1) / 2.0;
        double cn1 = n1 * (n1 - 1) / 2.0;
        
        choose_n0_vec.push_back(cn0);
        choose_n1_vec.push_back(cn1);
        total_choose_n1 += cn1;
        
        // Discordance Pre-calcs
        double p = n0 * n1;
        sum_prod += p;
        sum_sq_prod += p * p;
      }
      
      // Finalize Concordance
      double conc_val = 0;
      for (size_t k = 0; k < choose_n0_vec.size(); ++k) {
        // Pairs of (0,0) from state K and (1,1) from any OTHER state
        conc_val += choose_n0_vec[k] * (total_choose_n1 - choose_n1_vec[k]);
      }
      
      // Finalize Discordance
      double disc_val = 0.5 * (sum_prod * sum_prod - sum_sq_prod);
      
      // Store results
      concordant(s, c) = conc_val;
      // Decisive = Concordant + Discordant
      decisive(s, c) = conc_val + disc_val;
    }
  }
  
  return List::create(
    _["concordant"] = concordant, 
    _["decisive"] = decisive
  );
}
