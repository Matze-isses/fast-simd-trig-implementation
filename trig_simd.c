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

const double END_FRIENDLY_RANGE = M_PI * 2.0 * INT_MAX;

// ---------- SIN -----------
const double RANGE_MAX_SIN = M_PI * 2.0;
const double MED_RANGE_SIN = M_PI;
const double SMALL_RANGE_SIN = M_PI_2;

const double ONE_OVER_RANGE_SIN = 1 / RANGE_MAX_SIN;
const double ONE_OVER_MED_RANGE_SIN = 1 / MED_RANGE_SIN;
const double ONE_OVER_SMALL_RANGE_SIN = 1 / SMALL_RANGE_SIN;

const double RANGE_CENTER_SIN = 0;


// ---------- TAN -----------
const double RANGE_MAX_TAN = M_PI;
const double SMALL_RANGE_TAN = M_PI_2;
const double RANGE_CENTER_TAN = 0;
const double ONE_OVER_RANGE_TAN = 1 / RANGE_MAX_TAN;
const double ONE_OVER_SMALL_RANGE_TAN = 1 / SMALL_RANGE_TAN;

const double ONE_OVER_PI_2 = 1 / M_PI_2;

const int MAX_SIMD_DOUBLES = (int)(SIMD_LENGTH / 64);
const int MAX_SIMD_FLOAT = (int)(SIMD_LENGTH / 32);

const int TAYLOR_DEGREE = 20;
const int TAYLOR_LAST_COEFF = TAYLOR_DEGREE - 1;
const int TAYLOR_LOOP_INTERATIONS = TAYLOR_DEGREE - 2;
const double TAYLOR_COEFF_SIN[] = {
  0,
  1,
  -0,
  -0.16666666666666666,
  0,
  0.0083333333333333332,
  -0,
  -0.00019841269841269841,
  0,
  2.7557319223985893e-06,
  -0,
  -2.505210838544172e-08,
  0,
  1.6059043836821613e-10,
  -0,
  -7.6471637318198164e-13,
  0,
  2.8114572543455206e-15,
  -0,
  -8.2206352466243295e-18,
  0,
  1.9572941063391263e-20,
  -0,
  -3.8681701706306835e-23,
  0,
  6.4469502843844736e-26,
  -0,
  -9.183689863795546e-29,
  0,
  1.1309962886447718e-31,
  -0,
  -1.2161250415535181e-34,
  0,
  1.1516335620771951e-37,
  -0,
  -9.6775929586318907e-41,
  0,
  7.2654601791530714e-44,
  -0,
  -4.9024697565135435e-47,
  0,
  2.9893108271424046e-50,
  -0,
  -1.6552108677421951e-53,
  0,
  8.3596508471828045e-57,
  -0,
  -3.866628513960594e-60,
  0,
  1.6439747083165791e-63,
};

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


void get_reduced_range(double x, int *quadrant, double *reduced_range, double one_over_range_max, double range_max) {
  int n;

  n = floor(x * one_over_range_max);
  *reduced_range = x - n * range_max;
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

void sin_simd(double *input, double *res, size_t n, int prec) {
// TIME: 1329.5143421658272

  const int taylor_degree = (int)prec;
  const int taylor_last_coeff = taylor_degree - 1;
  const int taylor_loop_iteration = taylor_degree - 2;

  const SDOUBLE spi = LOAD_DOUBLE(M_PI);

  const SDOUBLE two_pi = LOAD_DOUBLE(RANGE_MAX_SIN);
  const SDOUBLE one_over_2_pi = LOAD_DOUBLE(ONE_OVER_RANGE_SIN);

  const SDOUBLE med_range = LOAD_DOUBLE(MED_RANGE_SIN);
  const SDOUBLE one_over_med_range = LOAD_DOUBLE(ONE_OVER_MED_RANGE_SIN);

  const SDOUBLE one_over_small_range = LOAD_DOUBLE(ONE_OVER_SMALL_RANGE_SIN);

  const SDOUBLE small_range = LOAD_DOUBLE(SMALL_RANGE_SIN);
  const SDOUBLE center_point = LOAD_DOUBLE(RANGE_CENTER_SIN);

  const SDOUBLE quadrant_multiplier = LOAD_DOUBLE(-2.0);
  const SDOUBLE ones = LOAD_DOUBLE(1.0);

  for (int i = 0; i < (int) n; i += MAX_SIMD_DOUBLES) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    // works but is potentially negative
    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_2_pi);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, two_pi);
    SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple);

    // Gives Sign of the Result
    const SDOUBLE medium_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_med_range);
    const SDOUBLE sign = FLOOR_DOUBLE_S(medium_ranges_away);

    // Gives Quadrant of the result
    const SDOUBLE small_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    const SDOUBLE q = FLOOR_DOUBLE_S(small_ranges_away);

    // [0, 1, 2, 3] -> [1, 0, 1, 2]
    const SDOUBLE q1 = ABS_PD(SUB_DOUBLE_S(q, ones));

    // [1, 0, 1, 2] -> [1, 0, 1, 4]
    const SDOUBLE q2 = MUL_DOUBLE_S(q1, q1);

    // [1, 0, 1, 4] -> [0, 1, 0, 3]
    const SDOUBLE q3 = ABS_PD(SUB_DOUBLE_S(q2, ones));

    // q3 * pi gives the mirroring points
    const SDOUBLE mirroring = MUL_DOUBLE_S(spi, q3);

    // all values mirroringare eighter in the first or in the thierd quadrant
    in_outer_range = SUB_DOUBLE_S(mirroring, in_outer_range);
    in_outer_range = ABS_PD(in_outer_range);

    const SDOUBLE small_ranges_away1 = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    const SDOUBLE q11 = FLOOR_DOUBLE_S(small_ranges_away1);
    const SDOUBLE initial_move = MUL_DOUBLE_S(small_range, q11);
    const SDOUBLE small_subtraction_amount = MUL_DOUBLE_S(q11, small_range);
    const SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, small_subtraction_amount);

    const SDOUBLE centered_values = SUB_DOUBLE_S(in_range, center_point);
    SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_SIN[taylor_last_coeff]);

    for (int j = taylor_loop_iteration; j >= 0; --j) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_SIN[j]);
      result = MUL_DOUBLE_S(result, centered_values);
      result = ADD_DOUBLE_S(result, coeff);
    }

    const SDOUBLE multiplied_quadrants = MUL_DOUBLE_S(sign, quadrant_multiplier);
    const SDOUBLE quadrant_evaluation = ADD_DOUBLE_S(multiplied_quadrants, ones);
    const SDOUBLE quadrant_evaluated_result = MUL_DOUBLE_S(result, quadrant_evaluation);

    SIMD_TO_DOUBLE_VEC(&res[i], quadrant_evaluated_result);
  }


  int num_left_over = (n % 4);

#pragma omp parallel for
  for (int i = n - num_left_over; i < (int)n; i++) {
    double reduced_range;
    int quadrant;
    get_reduced_range(input[i], &quadrant, &reduced_range, ONE_OVER_SMALL_RANGE_SIN, RANGE_MAX_SIN);
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


void tan_simd(double *input, double *res, size_t n, int prec) {
  const double cutoff = 0.07;
  const int taylor_degree = (int)prec;
  const int taylor_last_coeff = taylor_degree - 1;
  const int taylor_loop_iteration = taylor_degree - 2;

  double *approx_check = malloc(MAX_SIMD_DOUBLES * sizeof(double));

  const SDOUBLE pi_2 = LOAD_DOUBLE(M_PI_2);

  const SDOUBLE range_max = LOAD_DOUBLE(RANGE_MAX_TAN);
  const SDOUBLE one_over_range_max = LOAD_DOUBLE(ONE_OVER_RANGE_TAN);

  const SDOUBLE small_range = LOAD_DOUBLE(SMALL_RANGE_TAN);
  const SDOUBLE one_over_small_range = LOAD_DOUBLE(ONE_OVER_SMALL_RANGE_TAN);

  const SDOUBLE center_point = LOAD_DOUBLE(RANGE_CENTER_TAN);

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

    const SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, move_second_half_vec);
    SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_TAN[taylor_last_coeff]);

    for (int j = taylor_loop_iteration; j >= 0; --j) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result = MUL_DOUBLE_S(result, in_range);
      result = ADD_DOUBLE_S(result, coeff);
    }
    // PRINT_M256D(result);

    SIMD_TO_DOUBLE_VEC(approx_check, in_range);
    SIMD_TO_DOUBLE_VEC(&res[i], result);

    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      if (M_PI_2 - fabs(approx_check[j]) < cutoff) {
        res[i + j] = 1 / approx_check[j];
      }
    }
  }


  int num_left_over = (n % 4);

  #pragma omp parallel for
  for (int i = n - num_left_over; i < (int)n; i++) {
    double reduced_range;
    int quadrant;
    get_reduced_range(input[i], &quadrant, &reduced_range, ONE_OVER_SMALL_RANGE_TAN, RANGE_MAX_TAN);
    res[i] = taylor_eval(reduced_range, 0, TAYLOR_COEFF_TAN, taylor_degree);

    if (quadrant == 1) {
      res[i] = -res[i];
    }

    if (M_PI_2 - fabs(reduced_range) < cutoff) {
      res[i] = 1 / reduced_range;
    }
  }
}

// gcc ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/sin_arb.c ./sin_simd.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 12 0 100 0
