#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "trig_simd.h"
#include "./util/bit_printing.h"

#define MAX_EXACT_INT (9007199254740992)
#define LONG_RANGE_END (57952155664616982739.0)

const double RANGE_MAX = M_PI;
const double SMALL_RANGE = M_PI_2;
const double RANGE_CENTER = 0;

const double ONE_OVER_RANGE = 1 / RANGE_MAX;
const double ONE_OVER_SMALL_RANGE = 1 / SMALL_RANGE;

const double ONE_OVER_PI_2 = 1 / M_PI_2;

const int MAX_SIMD_DOUBLES = (int)(SIMD_LENGTH / 64);
const int MAX_SIMD_FLOAT = (int)(SIMD_LENGTH / 32);

const int TAYLOR_DEGREE = 27;
const int TAYLOR_LAST_COEFF = TAYLOR_DEGREE - 1;
const int TAYLOR_LOOP_INTERATIONS = TAYLOR_DEGREE - 2;

const double TAYLOR_COEFF_TAN[] = {
  0,
  1,
  0,
  0.33333333333333331,
  0,
  0.13333333333333333,
  0,
  0.053968253968253971,
  0,
  0.021869488536155203,
  0,
  0.0088632355299021973,
  0,
  0.0035921280365724811,
  0,
  0.0014558343870513183,
  0,
  0.00059002744094558595,
  0,
  0.00023912911424355248,
  0,
  9.6915379569294509e-05,
  0,
  3.9278323883316833e-05,
  0,
  1.5918905069328964e-05,
  0,
  6.4516892156554306e-06
};


void get_reduced_range(double x, int *quadrant, double *reduced_range) {
  int n;

  n = floor(x * ONE_OVER_RANGE);
  *reduced_range = x - n * RANGE_MAX;
  *quadrant = floor(*reduced_range * ONE_OVER_SMALL_RANGE);
  *quadrant = (*quadrant < 0) ? ((*quadrant + 4) % 4) : *quadrant;
  *reduced_range -= *quadrant * SMALL_RANGE;
}

double taylor_eval(double x, double a, const double coeffs[], int n) {
    double result = coeffs[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result = result * (x - a) + coeffs[i];
    }
    return result;
}

void sin_simd(double *input, double *res, size_t n, float prec){
  // other script just written to prevent header issues
}

void tan_simd(double *input, double *res, size_t n, float prec) {
  const SDOUBLE range_max = LOAD_DOUBLE(RANGE_MAX);
  const SDOUBLE one_over_range_max = LOAD_DOUBLE(ONE_OVER_RANGE);

  const SDOUBLE small_range = LOAD_DOUBLE(SMALL_RANGE);
  const SDOUBLE one_over_small_range = LOAD_DOUBLE(ONE_OVER_SMALL_RANGE);

  const SDOUBLE center_point = LOAD_DOUBLE(RANGE_CENTER);

  const SDOUBLE quadrant_multiplier = LOAD_DOUBLE(-2.0);
  const SDOUBLE addition_vector = LOAD_DOUBLE(1.0);
  
  #pragma omp parallel for
  for (int i = 0; i < (int) n; i += MAX_SIMD_DOUBLES) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    // works but is potentially negative
    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_range_max);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, range_max);
    const SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple);                       // in [0, 2 * pi]

    const SDOUBLE small_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    const SDOUBLE small_range_sub_part = FLOOR_DOUBLE_S(small_ranges_away);                     // used later
    const SDOUBLE move_second_half_vec = MUL_DOUBLE_S(small_range_sub_part, range_max);

    // this is checked and it is the correct range
    const SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, move_second_half_vec);

    SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_TAN[TAYLOR_LAST_COEFF]);

    for (int j = TAYLOR_LOOP_INTERATIONS; j >= 0; --j) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result = MUL_DOUBLE_S(result, in_range);
      result = ADD_DOUBLE_S(result, coeff);
    }


    SIMD_TO_DOUBLE_VEC(&res[i], result); 
  }

 
  int num_left_over = (n % 4);

  #pragma omp parallel for
  for (int i = n - num_left_over; i < (int)n; i++) {
    double reduced_range;
    int quadrant;
    get_reduced_range(input[i], &quadrant, &reduced_range);
    res[i] = taylor_eval(reduced_range, 0, TAYLOR_COEFF_TAN, TAYLOR_DEGREE);

    if (quadrant == 1) {
      res[i] = -res[i];
    }
  }
}

// gcc ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./tan_simd.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 100000000 0 1000 100000000
