// Measure std::getenv cost on this machine's environment (ucrt linear scan).
#include <cstdlib>
#include <chrono>
#include <cstdio>
int main() {
  const long N = 5'000'000;
  volatile int sink = 0;
  // warm
  for (long i = 0; i < 100000; ++i) sink += (std::getenv("TS_FREE_HTU_PROBE") != nullptr);
  double best = 1e9;
  for (int rep = 0; rep < 7; ++rep) {
    auto t0 = std::chrono::steady_clock::now();
    for (long i = 0; i < N; ++i) sink += (std::getenv("TS_FREE_HTU_PROBE") != nullptr);
    auto t1 = std::chrono::steady_clock::now();
    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / N;
    if (ns < best) best = ns;
  }
  printf("std::getenv: %.1f ns/call (best of 7, N=%ld); sink=%d\n", best, N, sink);
  printf("3 getenv/pick x 3840 picks = %.4f s\n", 3.0 * 3840 * best / 1e9);
  return 0;
}
