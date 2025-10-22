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
void get_simd_quadrant(double *src, double *quad, double *range) {
    SDOUBLE x   = LOAD_DOU_VEC(src);
    SDOUBLE two_pi = LOAD_DOU(RANGE_MAX);
    SDOUBLE one_over_2_pi = LOAD_DOU(ONE_OVER_RANGE);
    SDOUBLE one_over_pi_2 = LOAD_DOU(ONE_OVER_PI_2);

    // works but is potentially negative
    SDOUBLE n = MUL_DOU(x, one_over_2_pi);
    n = FLOOR_S(n);


    SDOUBLE range_multiple = MUL_DOU(n, two_pi);
    SDOUBLE in_range = SUB_DOUBLE_S(x, range_multiple); // in [0, 2 * pi]

    n = MUL_DOU(in_range, one_over_pi_2);
    n = FLOOR_S(n);

    SIMD_TO_DOUBLE_VEC(quad, n);
    SIMD_TO_DOUBLE_VEC(range, in_range);
}

void sin_simd(double *x, double *res, size_t n, float prec) {
  double *quadrants = malloc(MAX_SIMD_DOUBLES * sizeof(int));
  double *reduced_range = malloc(MAX_SIMD_DOUBLES * sizeof(double));

  double *partial_values = malloc(MAX_SIMD_DOUBLES * sizeof(double));
  
  for (int i = 0; i < n; i+= MAX_SIMD_DOUBLES) {
    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      partial_values[j] = x[i+j];
    }
    get_simd_quadrant(partial_values, quadrants, reduced_range);
    
    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      res[i+j] = quadrants[j];
    }
  }
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

    //  {pow(2, 57) * M_PI, "pi * 2^58", 0},

    //  // Those are the edge cases where it stops working
    //  {nextafter(M_PI * 2.0 * INT_MAX, -INFINITY), "2 * pi * intmax - eps", 3},

    //  {M_PI * 2.0 * INT_MAX, "2 * pi * intmax", 0},
    //  {nextafter(M_PI * 2.0 * INT_MAX, INFINITY), "2 * pi * intmax + eps", 0},
    //  {M_PI * 2.0 * INT_MAX + 1, "2 * pi * intmax + 1", 0},
    //  {M_PI * 2.0 * INT_MAX + M_PI, "pi * (2 * intmax + 1)", 2},

    //  // Converted large number to bits then took those and put them into a bit to decimal calculator (https://numeral-systems.com/ieee-754-converter/)
    //  // then put the result into an full precision calculator and checked the quadrant (https://www.mathsisfun.com/calculator-precision.html)
    //  {180899997887870864251481058976852284292172291974620700999680.0, "approx 1.8 * 10^59", 0},
    //  {9999999978010827105118843667605086770705495962419200.0, "approx 9 * 10^51", 1},
    //  {999999997801082606665947196063956106460623011840.0, "approx 9 * 10^47", 1},
    //  {679899997801082676275550042794554161600728662016.0, "approx 6.8 * 10^47", 1},
    //  {579899997887998780515705500125536825728609841840128.0, "approx 5.8 * 10^50", 1},
    //  {3698999978879987753059268883445507785324923154088429879296.0, "approx 3.7 * 10^57", 2},
    //  {1408999978878708665079569693857962203230771643241038251115788338092413288448.0, "approx 1.4 * 10^75", 2},
    //  {140899997887870848106493420455200807689536941778985583455291404091450082399540603781120.0, "approx 1.4 * 10^86", 2},
    //  {61000687969105996493672285664758005575629682910825211579427204492165120.0, "approx 6.1 * 10^70", 2}
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

  
  double *test_values = malloc(n * sizeof(double));
  double *test_results = malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    test_values[i] = ugly_tests[i].value;
  }

  sin_simd(test_values, test_results, n, 0.1);

  for (int i = 0; i < n; i++) {
    int q = (int)round(test_results[i]);
    printf("(SIMD)  %-25s is in quadrant %d", ugly_tests[i].name, q);

    if (q == ugly_tests[i].quadrant) {
      printf(" Correct:  True (%d)\n", ugly_tests[i].quadrant); // ]]
      correct_results_own[i] = true;
    } else {
      correct_results_own[i] = false;
      printf(" Correct: \033[31mFalse\033[0m (%d)\n", ugly_tests[i].quadrant); // ]]
    }
  }

  random_test();
  return 0;
}
// gcc range_reduction.c bit_printing.c -o range_reduction -lm -mavx -O2 && ./range_reduction
