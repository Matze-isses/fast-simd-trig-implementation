#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sin_simd.h"
#include "./util/bit_printing.h"

#define MAX_EXACT_INT (9007199254740992)
#define LONG_RANGE_END (57952155664616982739.0)

const double END_FRIENDLY_RANGE = M_PI * 2.0 * INT_MAX;

const double RANGE_MAX = M_PI * 2.0;
const double SMALL_RANGE = M_PI;
const double RANGE_CENTER = M_PI_2;

const double ONE_OVER_RANGE = 1 / RANGE_MAX;
const double ONE_OVER_SMALL_RANGE = 1 / SMALL_RANGE;

const double ONE_OVER_PI_2 = 1 / M_PI_2;

const int MAX_SIMD_DOUBLES = (int)(SIMD_LENGTH / 64);
const int MAX_SIMD_FLOAT = (int)(SIMD_LENGTH / 32);

const int TAYLOR_DEGREE = 20;
const int TAYLOR_LAST_COEFF = TAYLOR_DEGREE - 1;
const int TAYLOR_LOOP_INTERATIONS = TAYLOR_DEGREE - 2;
const double TAYLOR_COEFF_SIN[] = {
  1,
  6.123233995736766e-17,
  -0.5,
  -1.020538999289461e-17,
  0.041666666666666664,
  5.102694996447305e-19,
  -0.0013888888888888889,
  -1.2149273801065012e-20,
  2.4801587301587302e-05,
  1.6873991390368072e-22,
  -2.7557319223985888e-07,
  -1.5339992173061884e-24,
  2.08767569878681e-09,
  9.8333283160653097e-27,
  -1.1470745597729725e-11,
  -4.6825372933644332e-29,
  4.7794773323873853e-14,
  1.7215210637369241e-31,
  -1.5619206968586225e-16,
  -5.0336873208681989e-34
};


void get_reduced_range(double x, int *quadrant, double *reduced_range) {
  int n;

  n = floor(x * ONE_OVER_RANGE);
  *reduced_range = x - n * RANGE_MAX;
  *quadrant = floor(*reduced_range * ONE_OVER_PI_2);
  *quadrant = (*quadrant < 0) ? ((*quadrant + 4) % 4) : *quadrant;
}

double taylor_eval(double x, double a, const double coeffs[], int n) {
    double result = coeffs[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result = result * (x - a) + coeffs[i];
    }
    return result;
}

void sin_simd(double *input, double *res, size_t n, float prec) {
  const SDOUBLE two_pi = LOAD_DOUBLE(RANGE_MAX);
  const SDOUBLE one_over_2_pi = LOAD_DOUBLE(ONE_OVER_RANGE);
  const SDOUBLE one_over_small_range = LOAD_DOUBLE(ONE_OVER_SMALL_RANGE);

  const SDOUBLE small_range = LOAD_DOUBLE(SMALL_RANGE);
  const SDOUBLE center_point = LOAD_DOUBLE(RANGE_CENTER);

  const SDOUBLE quadrant_multiplier = LOAD_DOUBLE(-2.0);
  const SDOUBLE addition_vector = LOAD_DOUBLE(1.0);
  
  #pragma omp parallel for
  for (int i = 0; i < (int) n; i += MAX_SIMD_DOUBLES) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    // works but is potentially negative
    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_2_pi);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, two_pi);
    const SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple);                       // in [0, 2 * pi]

    const SDOUBLE small_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    const SDOUBLE simd_quadrants = FLOOR_DOUBLE_S(small_ranges_away);                     // used later
    const SDOUBLE small_subtraction_amount = MUL_DOUBLE_S(simd_quadrants, small_range);
    const SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, small_subtraction_amount);      // in smaller range

    const SDOUBLE centered_values = SUB_DOUBLE_S(in_range, center_point);

    SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_SIN[TAYLOR_LAST_COEFF]);

    for (int j = TAYLOR_LOOP_INTERATIONS; j >= 0; --j) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_SIN[j]);
      result = MUL_DOUBLE_S(result, centered_values);
      result = ADD_DOUBLE_S(result, coeff);
    }

    const SDOUBLE multiplied_quadrants = MUL_DOUBLE_S(simd_quadrants, quadrant_multiplier); 
    const SDOUBLE quadrant_evaluation = ADD_DOUBLE_S(multiplied_quadrants, addition_vector);
    const SDOUBLE quadrant_evaluated_result = MUL_DOUBLE_S(result, quadrant_evaluation);

    PRINT_M256D(quadrant_evaluated_result);
    SIMD_TO_DOUBLE_VEC(&res[i], quadrant_evaluated_result); 
  }

 
  int num_left_over = (n % 4);

  #pragma omp parallel for
  for (int i = n - num_left_over; i < (int)n; i++) {
    double reduced_range;
    int quadrant;
    get_reduced_range(input[i], &quadrant, &reduced_range);
    res[i] = taylor_eval(reduced_range, 0.5 * M_PI_2, TAYLOR_COEFF_SIN, TAYLOR_DEGREE);

    if (quadrant == 1) {
      res[i] = sqrt(1 - (res[i] * res[i]));
    } else if (quadrant == 2) {
      res[i] = - res[i];
    } else if (quadrant == 3) {
      res[i] = - sqrt(1 - (res[i] * res[i]));
    }
  }
}

// gcc ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/sin_arb.c ./sin_simd.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 12 0 100 0
