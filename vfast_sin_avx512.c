#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include "trig_simd.h"

#include <string.h>   // just used for the error calculation
// ---------- SIN -----------

// const double RANG_REDUCTION_CORRECTION = 3.8981718325193755e-17;
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


void vfast_sin(double *input, double *res, size_t n) {

//  int taylor_degree = 19;

//  const int taylor_last_coeff = taylor_degree - 1;
//  const int taylor_loop_iteration = taylor_degree - 2;

//  SET1_PD(range_reduction_correction, RANG_REDUCTION_CORRECTION);

//  SET1_PD(spi, M_PI);

//  SET1_PD(two_pi, RANGE_MAX_SIN);
//  SET1_PD(one_over_2_pi, ONE_OVER_RANGE_SIN);

//  SET1_PD(med_range, MED_RANGE_SIN);
//  SET1_PD(one_over_med_range, ONE_OVER_MED_RANGE_SIN);

//  SET1_PD(one_over_small_range, ONE_OVER_SMALL_RANGE_SIN);

//  SET1_PD(small_range, SMALL_RANGE_SIN);
//  SET1_PD(center_point, RANGE_CENTER_SIN);

//  SET1_PD(quadrant_multiplier, -2.0);
//  SET1_PD(ones, 1.0);

//  SET1_PD(taylor_coeff0, sin_tp0);
//  SET1_PD(taylor_coeff1, sin_tp1);
//  SET1_PD(taylor_coeff2, sin_tp2);
//  SET1_PD(taylor_coeff3, sin_tp3);
//  SET1_PD(taylor_coeff4, sin_tp4);
//  SET1_PD(taylor_coeff5, sin_tp5);
//  SET1_PD(taylor_coeff6, sin_tp6);
//  SET1_PD(taylor_coeff7, sin_tp7);
//  SET1_PD(taylor_coeff8, sin_tp8);
//  SET1_PD(taylor_coeff9, sin_tp9);
//  SET1_PD(taylor_coeff10, sin_tp10);
//  SET1_PD(taylor_coeff11, sin_tp11);
//  SET1_PD(taylor_coeff12, sin_tp12);
//  SET1_PD(taylor_coeff13, sin_tp13);
//  SET1_PD(taylor_coeff14, sin_tp14);
//  SET1_PD(taylor_coeff15, sin_tp15);
//  SET1_PD(taylor_coeff16, sin_tp16);
//  SET1_PD(taylor_coeff17, sin_tp17);
//  SET1_PD(taylor_coeff18, sin_tp18);
//  SET1_PD(taylor_coeff19, sin_tp19);

//  for (int i = 0; i < (int) n; i += SIMD_DOUBLES) {
//      SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

//      // works but is potentially negative
//      const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_2_pi);
//      const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
//      const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, two_pi);

//      SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple);
//      SDOUBLE correction_term = MUL_DOUBLE_S(x, range_reduction_correction);
//      in_outer_range = SUB_DOUBLE_S(in_outer_range, correction_term);

//      // Gives Sign of the Result
//      const SDOUBLE medium_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_med_range);
//      const SDOUBLE sign = FLOOR_DOUBLE_S(medium_ranges_away);

//      // Gives Quadrant of the result
//      const SDOUBLE small_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
//      const SDOUBLE q = FLOOR_DOUBLE_S(small_ranges_away);
//      const SDOUBLE q1 = ABS_PD(SUB_DOUBLE_S(q, ones));
//      const SDOUBLE q2 = MUL_DOUBLE_S(q1, q1);
//      const SDOUBLE q3 = ABS_PD(SUB_DOUBLE_S(q2, ones));

//      // q3 * pi gives the mirroring points, where 0 and 2 do not need to be mirrored
//      const SDOUBLE mirroring = MUL_DOUBLE_S(spi, q3);

//      // all values mirroringare eighter in the first or in the thierd quadrant
//      in_outer_range = SUB_DOUBLE_S(mirroring, in_outer_range);
//      in_outer_range = ABS_PD(in_outer_range);

//      const SDOUBLE small_ranges_away1 = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
//      const SDOUBLE q11 = FLOOR_DOUBLE_S(small_ranges_away1);
//      const SDOUBLE initial_move = MUL_DOUBLE_S(small_range, q11);
//      const SDOUBLE small_subtraction_amount = MUL_DOUBLE_S(q11, small_range);
//      const SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, small_subtraction_amount);

//      SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_SIN[taylor_last_coeff]);
//      const SDOUBLE x_square = MUL_DOUBLE_S(in_range, in_range);

//      /* ---- Taylor Loop ---- */
//      const SDOUBLE result_q0_t17 = FMADD_PD(taylor_coeff18, x_square, taylor_coeff17);
//      const SDOUBLE result_q0_t16 = FMADD_PD(result_q0_t17,  x_square, taylor_coeff16);
//      const SDOUBLE result_q0_t15 = FMADD_PD(result_q0_t16,  x_square, taylor_coeff15);
//      const SDOUBLE result_q0_t14 = FMADD_PD(result_q0_t15,  x_square, taylor_coeff14);
//      const SDOUBLE result_q0_t13 = FMADD_PD(result_q0_t14,  x_square, taylor_coeff13);
//      const SDOUBLE result_q0_t12 = FMADD_PD(result_q0_t13,  x_square, taylor_coeff12);
//      const SDOUBLE result_q0_t11 = FMADD_PD(result_q0_t12,  x_square, taylor_coeff11);
//      const SDOUBLE result_q0_t10 = FMADD_PD(result_q0_t11,  x_square, taylor_coeff10);
//      const SDOUBLE result_q0_t9  = FMADD_PD(result_q0_t10,  x_square, taylor_coeff9);
//      const SDOUBLE result_q0_t8  = FMADD_PD(result_q0_t9,   x_square, taylor_coeff8);
//      const SDOUBLE result_q0_t7  = FMADD_PD(result_q0_t8,   x_square, taylor_coeff7);
//      const SDOUBLE result_q0_t6  = FMADD_PD(result_q0_t7,   x_square, taylor_coeff6);
//      const SDOUBLE result_q0_t5  = FMADD_PD(result_q0_t6,   x_square, taylor_coeff5);
//      const SDOUBLE result_q0_t4  = FMADD_PD(result_q0_t5,   x_square, taylor_coeff4);
//      const SDOUBLE result_q0_t3  = FMADD_PD(result_q0_t4,   x_square, taylor_coeff3);
//      const SDOUBLE result_q0_t2  = FMADD_PD(result_q0_t3,   x_square, taylor_coeff2);
//      const SDOUBLE result_q0_t1  = FMADD_PD(result_q0_t2,   x_square, taylor_coeff1);
//      const SDOUBLE result_q0_t0  = FMADD_PD(result_q0_t1,   x_square, taylor_coeff0);

//      // to uneven the degrees
//      result = MUL_DOUBLE_S(result_q0_t0, in_range);

//      const SDOUBLE multiplied_quadrants = MUL_DOUBLE_S(sign, quadrant_multiplier);
//      const SDOUBLE quadrant_evaluation = ADD_DOUBLE_S(multiplied_quadrants, ones);
//      const SDOUBLE quadrant_evaluated_result = MUL_DOUBLE_S(result, quadrant_evaluation);

//      SIMD_TO_DOUBLE_VEC(&res[i], quadrant_evaluated_result);
//  }

//  int num_left_over = (n % 4);

//  for (int i = n - num_left_over; i < (int)n; i++) {
//      res[i] = sin(input[i]);
//  }
}

// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/trig_arb_comparison.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 1000000000 -8 8 1000000
//
// Compilation and run on laptop
//
//
// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/trig_arb_comparison.c -o test -lm -mavx -mavx2 -mfma -O2 -Wextra $(pkg-config --cflags --libs flint)  && ./test 1000000 -8 8 100000


//
