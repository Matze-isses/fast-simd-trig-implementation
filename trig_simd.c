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

const double RANGE_CENTER_SIN = M_PI_2;


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
  -5.0336873208681989e-34,
  4.1103176233121648e-19,
  1.198496981159095e-36,
  -8.8967913924505741e-22,
  -2.3685711090100691e-39,
  1.6117375710961184e-24,
  3.9476185150167817e-42,
  -2.4795962632247972e-27,
  -5.6233881980296044e-45,
  3.2798892370698378e-30,
  6.9253549236817788e-48,
  -3.7699876288159054e-33,
  -7.4466181975072882e-51,
  3.8003907548547441e-36,
  7.0517217779425086e-54,
  -3.3871575355211618e-39,
  -5.9258166201197543e-57,
  2.688220266286636e-42,
  4.4488112763661821e-60,
  -1.911963205040282e-45,
  -3.0018969476155074e-63,
  1.225617439128386e-48,
  1.8304249680582362e-66,
  -7.117406731291439e-52,
  -1.0135243455471962e-69,
  3.7618428812322616e-55,
  5.1188098259959402e-73,
  -1.817315401561479e-58,
  -2.367627116556864e-76,
  8.0554760707512382e-62,
  1.0066441822095511e-79
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

    const SDOUBLE medium_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_med_range);
    // Sign of the sin
    const SDOUBLE sign = FLOOR_DOUBLE_S(medium_ranges_away);

    const SDOUBLE small_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    // quadrant
    const SDOUBLE q = FLOOR_DOUBLE_S(small_ranges_away);

    PRINT_M256D(q);

    SDOUBLE q1 = SUB_DOUBLE_S(q, ones);
    q1 = ABS_PD(q1); 
    // PRINT_M256D(q1);

    SDOUBLE q2 = MUL_DOUBLE_S(q1, q1);
    // PRINT_M256D(q2);

    SDOUBLE q3 = SUB_DOUBLE_S(q2, ones);
    // Those makes the quadrant to a vector for the mirroring points
    // 0 -> 0
    // 1 -> 1
    // 2 -> 0
    // 3 -> 3
    q3 = ABS_PD(q3);

    // PRINT_M256D(q3);


    const SDOUBLE mirroring = MUL_DOUBLE_S(spi, q3);

    // all values are eighter in the first or in the thierd quadrant
    in_outer_range = SUB_DOUBLE_S(mirroring, in_outer_range);
    PRINT_M256D(in_outer_range);

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
