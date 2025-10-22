#include "tests/clock_utils.h"
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "range_reduction.h"
#include "bit_printing.h"

#define MAX_EXACT_INT (9007199254740992)
#define LONG_RANGE_END (57952155664616982739.0)

const double END_FRIENDLY_RANGE = M_PI * 2.0 * INT_MAX;

const double RANGE_MAX = M_PI * 2.0;
const double TWO_POW_NEG_49 = pow(2, -49);
const double TWO_POW_NEG_99 = pow(2, -99);
const double ONE_OVER_RANGE = 1 / RANGE_MAX;
const double ONE_OVER_PI_2 = 1 / M_PI_2;
const int MAX_SIMD_DOUBLES = (int)(SIMD_LENGTH / 64);
const int TAYLOR_DEGREE = 8;



typedef struct {
  double value;
  const char *name;
  int quadrant;
} UglyCaseTriplet;


int get_quadrant(double x) {
  int quadrant;
  double reduced_range;
  int n;

  n = floor(x * ONE_OVER_RANGE);
  reduced_range = x - n * RANGE_MAX;
  quadrant = floor(reduced_range * ONE_OVER_PI_2);
  quadrant = (quadrant < 0) ? ((quadrant + 4) % 4) : quadrant;

  return quadrant;
}


// When using the -O2 this is faster. For non optimized code this is slower.
void get_simd_quadrant_double(double *src, double *quad, double *range) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(src);
    SDOUBLE two_pi = LOAD_DOUBLE(RANGE_MAX);
    SDOUBLE one_over_2_pi = LOAD_DOUBLE(ONE_OVER_RANGE);
    SDOUBLE one_over_pi_2 = LOAD_DOUBLE(ONE_OVER_PI_2);

    // works but is potentially negative
    SDOUBLE n = MUL_DOUBLE_S(x, one_over_2_pi);
    n = FLOOR_DOUBLE_S(n);


    SDOUBLE range_multiple = MUL_DOUBLE_S(n, two_pi);
    SDOUBLE in_range = SUB_DOUBLE_S(x, range_multiple); // in [0, 2 * pi]

    n = MUL_DOUBLE_S(in_range, one_over_pi_2);
    n = FLOOR_DOUBLE_S(n);

    SIMD_TO_DOUBLE_VEC(quad, n);
    SIMD_TO_DOUBLE_VEC(range, in_range);
}

void get_simd_quadrant_float(float *src, float *quad, float *range) {
  SFLOAT x = LOAD_FLOAT_VEC(src);
  SFLOAT two_pi = LOAD_FLOAT((float) RANGE_MAX);
  SFLOAT one_over_2_pi = LOAD_FLOAT((float) ONE_OVER_RANGE);
  SFLOAT one_over_pi_2 = LOAD_FLOAT((float) ONE_OVER_PI_2);

  SFLOAT n = MUL_FLOAT_S(x, one_over_pi_2);
  n = FLOOR_FLOAT_S(n);

  SFLOAT range_multiple = MUL_FLOAT_S(n, two_pi);
  SFLOAT in_range = SUB_FLOAT_S(x, range_multiple); // in [0, 2 * pi]

  n = MUL_FLOAT_S(in_range, one_over_pi_2);
  n = FLOOR_FLOAT_S(n);

  SIMD_TO_FLOAT_VEC(quad, n);
  SIMD_TO_FLOAT_VEC(range, in_range);
}

/*
 * It is assumed that the src vector has the correct length given the AVX version and that
 * the ranges are already reduced to [0, pi/2]. This function only works for calculating the 
 * values in the first quadrant.
 */
void taylor_sin(double *src, double *res) {
}

void sin_simd_float(float *x, float *res, size_t n, float prec) {
  float *quadrants = malloc(MAX_SIMD_DOUBLES * sizeof(int));
  float *reduced_range = malloc(MAX_SIMD_DOUBLES * sizeof(float));
  float *partial_values = malloc(MAX_SIMD_DOUBLES * sizeof(float));
  
  for (int i = 0; i < n; i+= MAX_SIMD_DOUBLES) {
    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      partial_values[j] = x[i+j];
    }

    get_simd_quadrant_float(partial_values, quadrants, reduced_range);
    
    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      res[i+j] = quadrants[j];
    }
  }

  free(quadrants);
  free(reduced_range);
  free(partial_values);
}


void sin_simd(double *x, double *res, size_t n, float prec) {
  double *quadrants = malloc(MAX_SIMD_DOUBLES * sizeof(double));
  double *reduced_range = malloc(MAX_SIMD_DOUBLES * sizeof(double));
  double *partial_values = malloc(MAX_SIMD_DOUBLES * sizeof(double));
  
  for (int i = 0; i < n; i+= MAX_SIMD_DOUBLES) {
    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      partial_values[j] = x[i+j];
    }

    get_simd_quadrant_double(partial_values, quadrants, reduced_range);
    
    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      res[i+j] = quadrants[j];
    }
  }

  free(quadrants);
  free(reduced_range);
  free(partial_values);
}

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
  double *res = malloc(n * sizeof(double));

  if (!vec) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower, upper, n, vec);

  START_CLOCK;

  for (int i = 0; i < n; i++) {
    res[i] = get_quadrant(vec[i]);
  }

  END_CLOCK("Quadrant Calculation WARMUP");

  START_CLOCK;

  for (int i = 0; i < n; i++) {
    res[i] = get_quadrant(vec[i]);
  }

  END_CLOCK("Quadrant Calculation       ");

  START_CLOCK;

  sin_simd(vec, res, n, 0.1);

  END_CLOCK("Quadrant Calculation SIMD");

  free(vec);
}


// ######################## Test with ugly values #######################################

// Example usage
int main(void) {

  UglyCaseTriplet ugly_tests[] = {
    {0.0, "0", 0},
    {M_PI_2, "pi / 2", 1},
    {M_PI, "pi", 2},
    {M_PI * 1.5, "1.5 * pi", 3},
    {M_PI * 2.0, "2.0 * pi", 0},

    {-M_PI_2, "- pi / 2", 3},
    {-M_PI, "- pi", 2},
    {-M_PI * 1.5, "- 1.5 * pi", 1},
    {-M_PI * 2.0, "- 2.0 * pi", 0},

    {nextafter(0.0, -INFINITY), "0 - eps", 3},
    {nextafter(M_PI_2, -INFINITY), "pi / 2 - eps", 0},
    {nextafter(M_PI, -INFINITY), "pi - eps", 1},
    {nextafter(M_PI * 1.5, -INFINITY), "1.5 * pi - eps", 2},
    {nextafter(M_PI * 2.0, -INFINITY), "2.0 * pi - eps", 3},

    {nextafter(0.0, INFINITY),  "0 + eps", 0},
    {nextafter(M_PI_2, INFINITY), "pi / 2 + eps", 1},
    {nextafter(M_PI, INFINITY), "pi + eps", 2},
    {nextafter(M_PI * 1.5, INFINITY), "1.5 * pi + eps", 3},
    {nextafter(M_PI * 2.0, INFINITY), "2.0 * pi + eps", 0},

    {nextafter(-M_PI_2, INFINITY), "-pi / 2 + eps", 3},
    {nextafter(-M_PI, INFINITY), "-pi + eps", 2},
    {nextafter(-M_PI * 1.5, INFINITY), "-1.5 * pi + eps", 1},
    {nextafter(-M_PI * 2.0, INFINITY), "-2.0 * pi + eps", 0},

    {nextafter(M_PI * 16.5, INFINITY), "16.5 * pi + eps", 1},
    {nextafter(M_PI * 17.0, INFINITY), "17.0 * pi + eps", 2},
    {nextafter(M_PI * 17.5, INFINITY), "17.5 * pi + eps", 3},
    {nextafter(M_PI * 18.0, INFINITY), "18.0 * pi + eps", 0},


    {nextafter(- M_PI * 16.5, INFINITY), "- 16.5 * pi + eps", 3},
    {nextafter(- M_PI * 17.0, INFINITY), "- 17.0 * pi + eps", 2},
    {nextafter(- M_PI * 17.5, INFINITY), "- 17.5 * pi + eps", 1},
    {nextafter(- M_PI * 18.0, INFINITY), "- 18.0 * pi + eps", 0},

    {44.5, "44.5", 0},
    {46.25, "46.25", 1},
    {47.5, "47.5", 2},
    {49.75, "49.75", 3},

    {255.0, "255.0", 2},

    {302.0, "302.0", 0},
    {305.5, "305.5", 2},

    // The following examples have at the start of all 8 columns as long values
    // 0, which should make the calculation of t exact.
    {6316927.0, "6316927", 0},
    {6324095.0, "6324095", 3},
    {1077952576.0, "1077952576", 3},

    // The following have at the start of some 8 columns as long a 1,
    // which should destroy the calculation of t
    {6356612.0, "6356612.0", 0},
    {6356612.0, "6356612.0", 0},
    {6356630.0, "6356630.0", 0},
    {6324223.0, "6324223.0", 1},
    {6356632.0, "6356632.0", 1},
    {6356991.0, "6356991.0", 2},
    {6356608.0, "6356608.0", 2},
    {6356628.0, "6356628.0", 3},
    {6356635.0, "6356635.0", 3},
    {23788.5, "23788.5", 0},
  };

  int n = sizeof(ugly_tests) / sizeof(ugly_tests[0]);
  bool correct_results_own[n];
  int quadrants[n];

  for (int i = 0; i < n; i++) {
    int q = get_quadrant(ugly_tests[i].value);
    quadrants[n] = q;
    printf("  %-25s is in quadrant %d", ugly_tests[i].name, q);

    if (q == ugly_tests[i].quadrant) {
      printf(" Correct:  True (%d)\n", ugly_tests[i].quadrant); // ]]
      correct_results_own[i] = true;
    } else {
      correct_results_own[i] = false;
      printf(" Correct: \033[31mFalse\033[0m (%d)\n", ugly_tests[i].quadrant); // ]]
    }
  }

  
  float* test_values_float = malloc(n * sizeof(float));
  float* test_results_float = malloc(n * sizeof(float));

  double *test_values = malloc(n * sizeof(double));
  double *test_results = malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    test_values[i] = ugly_tests[i].value;
  }

  sin_simd(test_values, test_results, n, 0.1);

  for (int i = 0; i < n; i++) {
    int q = (int)round(test_results[i]);
    printf("(SIMD DOUBLE)  %-25s is in quadrant %d", ugly_tests[i].name, q);

    if (q == ugly_tests[i].quadrant) {
      printf(" Correct:  True (%d)\n", ugly_tests[i].quadrant); // ]]
      correct_results_own[i] = true;
    } else {
      correct_results_own[i] = false;
      printf(" Correct: \033[31mFalse\033[0m (%d)\n", ugly_tests[i].quadrant); // ]]
    }
  }

  free(test_values);
  free(test_results);


  printf("\n\nHere\n\n");

  for (int i = 0; i < n; i++) {
    test_values_float[i] = (float)ugly_tests[i].value;
  }

  sin_simd_float(test_values_float, test_results_float, n, 0.1);

  for (int i = 0; i < n; i++) {
    int q = (int)round(test_results_float[i]);
    printf("(SIMD FLOAT)  %-25s is in quadrant %d", ugly_tests[i].name, q);

    if (q == ugly_tests[i].quadrant) {
      printf(" Correct:  True (%d)\n", ugly_tests[i].quadrant); // ]]
      correct_results_own[i] = true;
    } else {
      correct_results_own[i] = false;
      printf(" Correct: \033[31mFalse\033[0m (%d)\n", ugly_tests[i].quadrant); // ]]
    }
  }

  random_test();

  free(test_values_float);
  free(test_results_float);
  return 0;
}
// gcc range_reduction.c bit_printing.c -o range_reduction -lm -mavx -O2 && ./range_reduction
