#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include "trig_simd.h"

#include <string.h>   // just used for the error calculation
                      
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

const int SIZE_TAYLOR_COEFF = 20;
const double TAYLOR_COEFF_SIN[] = {
  1,
  -0.16666666666666666,
  0.0083333333333333332,
  -0.00019841269841269841,
  2.7557319223985893e-06,
  -2.505210838544172e-08,
  1.6059043836821613e-10,
  -7.6471637318198164e-13,
  2.8114572543455206e-15,
  -8.2206352466243295e-18,
  1.9572941063391263e-20,
  -3.8681701706306835e-23,
  6.4469502843844736e-26,
  -9.183689863795546e-29,
  1.1309962886447718e-31,
  -1.2161250415535181e-34,
  1.1516335620771951e-37,
  -9.6775929586318907e-41,
  7.2654601791530714e-44,
  -4.9024697565135435e-47,
  2.9893108271424046e-50,
};

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

static inline int taylor_degree_from_prec(double max_input, double req_prec) {
    uint64_t bits;
    memcpy(&bits, &max_input, sizeof bits);

    // extract 11-bit exponent field
    int exp_bits = (int)((bits >> 52) & 0x7FF);

    // subtract IEEE-754 bias (1023) + 1 due to the lost bit
    int exponent =  exp_bits - 1023;
    int bit_loss = 52 - exponent;
    double item_loss = pow(2, -bit_loss);

    // get taylor loss if used 0 degree
    double taylor_loss = 0;
    for (int i = 0; i < SIZE_TAYLOR_COEFF; i++) { taylor_loss += fabs(TAYLOR_COEFF_SIN[i]); }

    int used_coeffs = 0;
    double total_loss = item_loss + taylor_loss;

    while (total_loss > req_prec) {
      total_loss -= fabs(TAYLOR_COEFF_SIN[used_coeffs]);
      used_coeffs += 1;

      if (used_coeffs >= SIZE_TAYLOR_COEFF) {
        printf("[WARNING] The required precision of sin_simd cannot be reached! The best possible precision is %.17g, which is used!", total_loss);
        break;
      }
    }

    used_coeffs += 1; // safety ensurance

    // printf("Precision: %.17g -> Taylor Degree: %d\n", total_loss, used_coeffs);
    return used_coeffs;
}

void sin_simd(double *input, double *res, size_t n, double prec) {
  double max_element = input[0];
  // for (size_t i = 1; i < n; ++i) { if (input[i] > max_element) max_element = input[i]; }
  // int taylor_degree = taylor_degree_from_prec(max_element, prec);
  int taylor_degree = 19;

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

    // q3 * pi gives the mirroring points, where 0 and 2 do not need to be mirrored
    const SDOUBLE mirroring = MUL_DOUBLE_S(spi, q3);

    // all values mirroringare eighter in the first or in the thierd quadrant
    in_outer_range = SUB_DOUBLE_S(mirroring, in_outer_range);
    in_outer_range = ABS_PD(in_outer_range);

    const SDOUBLE small_ranges_away1 = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    const SDOUBLE q11 = FLOOR_DOUBLE_S(small_ranges_away1);
    const SDOUBLE initial_move = MUL_DOUBLE_S(small_range, q11);
    const SDOUBLE small_subtraction_amount = MUL_DOUBLE_S(q11, small_range);
    const SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, small_subtraction_amount);

    SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_SIN[taylor_last_coeff]);
    const SDOUBLE x_square = MUL_DOUBLE_S(in_range, in_range);

    for (int j = taylor_loop_iteration; j >= 0; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_SIN[j]);
      result = FMADD_PD(result, x_square, coeff);
    }

    // to uneven the degrees
    result = MUL_DOUBLE_S(result, in_range);

    const SDOUBLE multiplied_quadrants = MUL_DOUBLE_S(sign, quadrant_multiplier);
    const SDOUBLE quadrant_evaluation = ADD_DOUBLE_S(multiplied_quadrants, ones);
    const SDOUBLE quadrant_evaluated_result = MUL_DOUBLE_S(result, quadrant_evaluation);

    SIMD_TO_DOUBLE_VEC(&res[i], quadrant_evaluated_result);
  }


  int num_left_over = (n % 4);

  for (int i = n - num_left_over; i < (int)n; i++) {
    res[i] = sin(input[i]);
  }
}

// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/trig_arb_comparison.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 1000000000 -8 8 1000000
//
// Compilation and run on laptop
//
//
// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/trig_arb_comparison.c -o test -lm -mavx -mavx2 -mfma -O2 -Wextra $(pkg-config --cflags --libs flint)  && ./test 1000000 -8 8 100000


//
