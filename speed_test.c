// ######################## Test with uniform numbers for speed #######################################

// Helper: generate a uniform random double in [0, 1]
static inline double uniform01(void) {
  unsigned long long high = (unsigned long long)rand();
  unsigned long long low = (unsigned long long)rand();
  unsigned long long combined = (high << 31) | low;
  return (double)combined / (double)((1ULL << 62) - 1);
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
}

int random_test() {
  srand((unsigned)time(NULL));

  size_t n = 1000000000;
  double lower = 6.5;
  double upper = 1000000000.0;

  double *vec = malloc(n * sizeof(double));

  if (!vec) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower, upper, n, vec);

  START_CLOCK;

  for (int i = 0; i < n; i++) {
    get_quadrant(vec[i]);
  }

  END_CLOCK("Quadrant Calculation WARMUP");

  START_CLOCK;

  for (int i = 0; i < n; i++) {
    get_quadrant(vec[i]);
  }

  END_CLOCK("Quadrant Calculation       ");

  free(vec);
}
