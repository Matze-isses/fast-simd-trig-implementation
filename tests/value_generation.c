#include "clock_utils.h"
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

// Helper: generate a uniform random double in [0, 1]
static inline double uniform01(void) {
  unsigned long long high = (unsigned long long)rand();
  unsigned long long low = (unsigned long long)rand();
  unsigned long long combined = (high << 31) | low;
  return (double)combined / (double)((1ULL << 62) - 1);
}

// Standard normal N(0,1) via Box–Muller
static inline double normal01(void) {
  // avoid log(0)
  double u1 = uniform01();
  if (u1 < 1e-15) u1 = 1e-15;
  double u2 = uniform01();
  return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// Main function: fill vector with uniform random doubles in [lower, upper]
void fill_uniform(double lower, double upper, size_t n, double *vec) {
  if (upper < lower) {
    fprintf(stderr, "Error: upper bound < lower bound.\n");
    return;
  }
  double range = upper - lower;
  for (size_t i = 0; i < n; i++) {
    vec[i] = lower + uniform01() * range;
  }

  vec[0] = lower;
  vec[n - 1] = upper;
}

void fill_dense_pi_over_2(double lower, double upper, size_t n, double *vec, double sigma) {
  if (upper < lower) {
    fprintf(stderr, "Error: upper bound < lower bound.\n");
    return;
  }
  if (!(sigma > 0.0)) {
    fprintf(stderr, "Error: sigma must be > 0.\n");
    return;
  }

  const double step = M_PI / 2.0;

  // k such that k*(pi/2) lies in [lower, upper]
  long long k_min = (long long)ceil(lower / step);
  long long k_max = (long long)floor(upper / step);

  // If the interval is smaller than one step and contains no center,
  // fall back to uniform (or you could instead bias around nearest center).
  if (k_max < k_min) {
    fill_uniform(lower, upper, n, vec);
    return;
  }

  unsigned long long count = (unsigned long long)(k_max - k_min + 1);

  for (size_t i = 0; i < n; i++) {
    // cap attempts to avoid pathological infinite loops for tiny intervals + tiny sigma
    for (int attempt = 0; attempt < 10000; attempt++) {
      // choose a center uniformly among available k's
      unsigned long long idx = (unsigned long long)(uniform01() * (double)count);
      if (idx >= count) idx = count - 1;

      long long k = k_min + (long long)idx;
      double center = (double)k * step;

      // sample around center
      double x = center + sigma * normal01();

      if (x >= lower && x <= upper) {
        vec[i] = x;
        break;
      }

      // fallback if we failed too often
      if (attempt == 9999) {
        vec[i] = lower + uniform01() * (upper - lower);
      }
    }
  }
}
