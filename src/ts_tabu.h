#ifndef TS_TABU_H
#define TS_TABU_H

// Fixed-size circular buffer of topology hashes for tabu search.
//
// During TBR plateau exploration (accept_equal moves), the search can
// cycle between previously visited topologies. The tabu list prevents
// this by remembering recent topology hashes and rejecting moves that
// would revisit them.
//
// Linear scan for contains() — capacity is small (typically 100–1000)
// and data is contiguous, so this is faster than unordered_set.

#include <vector>
#include <cstdint>

namespace ts {

class TabuList {
  std::vector<uint64_t> buf_;
  int capacity_;
  int pos_ = 0;
  int count_ = 0;

public:
  explicit TabuList(int capacity)
      : buf_(capacity > 0 ? capacity : 0, 0), capacity_(capacity > 0 ? capacity : 0) {}

  // Add a hash (overwrites oldest entry when full).
  void insert(uint64_t hash) {
    if (capacity_ <= 0) return;
    buf_[pos_] = hash;
    pos_ = (pos_ + 1) % capacity_;
    if (count_ < capacity_) ++count_;
  }

  // Check whether a hash is in the list.
  bool contains(uint64_t hash) const {
    for (int i = 0; i < count_; ++i) {
      if (buf_[i] == hash) return true;
    }
    return false;
  }

  void clear() { pos_ = 0; count_ = 0; }

  int size() const { return count_; }
  int capacity() const { return capacity_; }
  bool active() const { return capacity_ > 0; }
};

} // namespace ts

#endif // TS_TABU_H
