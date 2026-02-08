// bench.c
#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>

#ifndef N
#define N 1000000000ULL
#endif

static volatile uint64_t g_sink = 0;

static inline uint64_t ns_now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

// Simple fast PRNG (xorshift32). Ensures non-zero outputs.
static inline uint32_t xorshift32(uint32_t *state) {
  uint32_t x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return x ? x : 1u; // force non-zero
}

__attribute__((noinline)) static uint32_t func_a(uint32_t x) {

  return x;
}

__attribute__((noinline)) static uint32_t func_b(uint32_t x) {
  return x;
}
// --------------------

typedef uint32_t (*fn_t)(uint32_t);

static uint64_t run_one(const uint32_t *arr, size_t n, fn_t f) {
  // Accumulate result to force the loop to matter.
  uint64_t acc = 0;
  for (size_t i = 0; i < n; i++) {
    acc += f(arr[i]);
  }
  // Make it observable (prevents dead-code elimination).
  g_sink ^= acc;
  return acc;
}

static uint64_t time_one(const char *name, const uint32_t *arr, size_t n, fn_t f) {
  // Warm-up (optional, helps reduce first-touch effects)
  run_one(arr, n / 100, f);

  uint64_t t0 = ns_now();
  uint64_t acc = run_one(arr, n, f);
  uint64_t t1 = ns_now();

  uint64_t dt = t1 - t0;
  double sec = (double)dt / 1e9;
  double ns_per = (double)dt / (double)n;

  printf("%s: %.3f s, %.2f ns/elem, checksum=%llu\n",
         name, sec, ns_per, (unsigned long long)acc);
  return dt;
}

int main(void) {
  const size_t n = (size_t)N;

  // 4GB for 1e9 uint32_t. Ensure your machine has enough RAM.
  uint32_t *arr = NULL;
  int rc = posix_memalign((void**)&arr, 64, n * sizeof(uint32_t));
  if (rc != 0) {
    fprintf(stderr, "posix_memalign failed: %s\n", strerror(rc));
    return 1;
  }

  // Fill with non-zero values.
  uint32_t s = 1u;
  for (size_t i = 0; i < n; i++) {
    arr[i] = xorshift32(&s);
  }

  // Time both (same input).
  uint64_t dt_a = time_one("func_a", arr, n, func_a);
  uint64_t dt_b = time_one("func_b", arr, n, func_b);

  // Simple ratio
  if (dt_b) {
    printf("Speed ratio (a/b): %.3f\n", (double)dt_a / (double)dt_b);
  }

  // Print volatile sink so it can't be “assumed unused”
  printf("sink=%llu\n", (unsigned long long)g_sink);

  free(arr);
  return 0;
}
