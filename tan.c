#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <stdint.h>

#include <stdlib.h>
#include "trig_simd.h"

#include <string.h>   // just used for the error calculation
                      
#include "./util/bit_printing.h"



double TAYLOR_COEFF_TAN[] = {
  0.3333333333333333,
  0.13333333333333333,
  0.05396825396825397,
  0.021869488536155203,
  0.008863235529902197,
  0.003592128036572481,
  0.0014558343870513183,
  0.000590027440945586,
  0.00023912911424355248,
  9.691537956929451e-05,
  3.927832388331683e-05,
  1.5918905069328964e-05,
  6.451689215655431e-06,  /* 6.451689215655431e-06 | 27-th degree*/
  2.6147711512907546e-06 /* 2.6147711512907546e-06 */ // <-- Current last
};


void tan_simd(double *input, double *res, size_t n) {
  int simd_doubles = SIMD_LENGTH / 64;

  const SDOUBLE pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  const SDOUBLE one_over_pi_2 = LOAD_DOUBLE(1/M_PI_2);

  const SDOUBLE q2_bitshift = LOAD_DOUBLE(-pow(2, -52));
  const SDOUBLE b_correction = LOAD_DOUBLE(MIN_POSITIVE_COS_VALUE);

  const SDOUBLE neg_one = LOAD_DOUBLE(-1.0);
  const SDOUBLE neg_half = LOAD_DOUBLE(-0.5);
  const SDOUBLE zero = SET_ZERO();
  const SDOUBLE half = LOAD_DOUBLE(0.5);
  const SDOUBLE one = LOAD_DOUBLE(1.0);
  const SDOUBLE two = LOAD_DOUBLE(2.0);
  const SDOUBLE three = LOAD_DOUBLE(3.0);
  
  for (int i = 0; i < (int) n; i += SIMD_DOUBLES) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    const SDOUBLE abs_x = ABS_PD(x);
    const SDOUBLE x_negative = DIV_DOUBLE_S(x, abs_x);

    const SDOUBLE x_negative01_0 = MUL_DOUBLE_S(x_negative, neg_one);
    const SDOUBLE x_negative01_1 = ADD_DOUBLE_S(x_negative01_0, one);
    const SDOUBLE x_negative01 = MUL_DOUBLE_S(x_negative01_1, half);

    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_pi_2);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, pi_2);

    x = SUB_DOUBLE_S(x, range_multiple);

    // here problems occure because negative singularities need a ceil not a floor
    const SDOUBLE singularities_away0 = ADD_DOUBLE_S(num_ranges_away, x_negative01);
    const SDOUBLE singularities_away1 = ABS_PD(singularities_away0);
    const SDOUBLE singularities_away = MUL_DOUBLE_S(singularities_away1, half);
    const SDOUBLE num_singularities_away = FLOOR_DOUBLE_S(singularities_away);

    const SDOUBLE constants_away0 = MUL_DOUBLE_S(num_singularities_away, two);
    const SDOUBLE constants_away = ADD_DOUBLE_S(constants_away0, one);

    // Check if even
    //  Default Range Reduction
    const SDOUBLE sign_adjust_0 = MUL_DOUBLE_S(num_ranges_away, half);
    const SDOUBLE sign_adjust_1 = FLOOR_DOUBLE_S(sign_adjust_0);
    const SDOUBLE sign_adjust_2 = MUL_DOUBLE_S(sign_adjust_1, two);

    const SDOUBLE in_odd_range  = SUB_DOUBLE_S(num_ranges_away, sign_adjust_2);
    const SDOUBLE in_even_range = SUB_DOUBLE_S(one, in_odd_range);

    const SDOUBLE sign_adjust_3 = MUL_DOUBLE_S(in_odd_range, two);
    const SDOUBLE sign_adjust = SUB_DOUBLE_S(one, sign_adjust_3);

    SDOUBLE from_behind = SUB_DOUBLE_S(pi_2, x);

    SDOUBLE in_odd_range_reduction = from_behind;
    in_odd_range_reduction = MUL_DOUBLE_S(in_odd_range_reduction, in_odd_range);
    in_odd_range_reduction = SUB_DOUBLE_S(in_odd_range_reduction, x);
    x = FMADD_PD(in_odd_range_reduction, in_odd_range, x);

    from_behind = SUB_DOUBLE_S(pi_2, x);

    const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
    const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);

    /* obtaining bool vectors for the each quadrant */
    // 1 if quadrant == 0 else 0
    SDOUBLE in_q0 = SUB_DOUBLE_S(quadrant, two);
    in_q0 = ABS_PD(in_q0);
    in_q0 = MUL_DOUBLE_S(in_q0, half);
    in_q0 = FLOOR_DOUBLE_S(in_q0);

    // 1 if quadrant == 1 else 0
    SDOUBLE in_q1 = SUB_DOUBLE_S(quadrant, one);
    in_q1 = ABS_PD(in_q1);
    in_q1 = SUB_DOUBLE_S(in_q1, two);
    in_q1 = ABS_PD(in_q1);
    in_q1 = MUL_DOUBLE_S(in_q1, half);
    in_q1 = FLOOR_DOUBLE_S(in_q1);

    // 1 if quadrant == 2 else 0
    SDOUBLE in_q2 = SUB_DOUBLE_S(quadrant, two);
    in_q2 = ABS_PD(in_q2);
    in_q2 = SUB_DOUBLE_S(in_q2, two);
    in_q2 = ABS_PD(in_q2);
    in_q2 = MUL_DOUBLE_S(in_q2, half);
    in_q2 = FLOOR_DOUBLE_S(in_q2);

    // 1 if quadrant == 3 else 0
    SDOUBLE in_q3 = SUB_DOUBLE_S(quadrant, one);
    in_q3 = ABS_PD(in_q3);
    in_q3 = MUL_DOUBLE_S(in_q3, half);
    in_q3 = FLOOR_DOUBLE_S(in_q3);

    // Mirror it to move it to range 1
    SDOUBLE q2_reduction = from_behind;
    q2_reduction = MUL_DOUBLE_S(q2_reduction, in_q2);
    q2_reduction = SUB_DOUBLE_S(q2_reduction, x);
    x = FMADD_PD(q2_reduction, in_q2, x);
    // x = q2_reduction if q2_reduction != 0 else x

    // reduce to move it to range 0
    SDOUBLE q1_reduction = MUL_DOUBLE_S(x, neg_half);
    x = FMADD_PD(q1_reduction, in_q1, x);
    x = FMADD_PD(q1_reduction, in_q2, x);

    // move q3 in q0
    SDOUBLE q3_reduction = SUB_DOUBLE_S(from_behind, x);
    x = FMADD_PD(q3_reduction, in_q3, x);
    
    /* ---- Taylor Loop ---- */
    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);
    const SDOUBLE result_q0_t0 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[12]);

    const SDOUBLE coeff12 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[11]);
    const SDOUBLE result_q0_t1 = FMADD_PD(result_q0_t0, x_square, coeff12);

    const SDOUBLE coeff11 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[10]);
    const SDOUBLE result_q0_t2 = FMADD_PD(result_q0_t1, x_square, coeff11);

    const SDOUBLE coeff10 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[9]);
    const SDOUBLE result_q0_t3 = FMADD_PD(result_q0_t2, x_square, coeff10);

    const SDOUBLE coeff9 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[8]);
    const SDOUBLE result_q0_t4 = FMADD_PD(result_q0_t3, x_square, coeff9);

    const SDOUBLE coeff8 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[7]);
    const SDOUBLE result_q0_t5 = FMADD_PD(result_q0_t4, x_square, coeff8);

    const SDOUBLE coeff7 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[6]);
    const SDOUBLE result_q0_t6 = FMADD_PD(result_q0_t5, x_square, coeff7);

    const SDOUBLE coeff6 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[5]);
    const SDOUBLE result_q0_t7 = FMADD_PD(result_q0_t6, x_square, coeff6);

    const SDOUBLE coeff5 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[4]);
    const SDOUBLE result_q0_t8 = FMADD_PD(result_q0_t7, x_square, coeff5);

    const SDOUBLE coeff4 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[3]);
    const SDOUBLE result_q0_t9 = FMADD_PD(result_q0_t8, x_square, coeff4);

    const SDOUBLE coeff3 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[2]);
    const SDOUBLE result_q0_t10 = FMADD_PD(result_q0_t9, x_square, coeff3);

    const SDOUBLE coeff2 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[1]);
    const SDOUBLE result_q0_t11 = FMADD_PD(result_q0_t10, x_square, coeff2);

    const SDOUBLE coeff1 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[0]);
    const SDOUBLE result_q0_t12 = FMADD_PD(result_q0_t11, x_square, coeff1);

    const SDOUBLE correction_sign_1 = SUB_DOUBLE_S(in_even_range, in_odd_range);
    const SDOUBLE correction_sign  = MUL_DOUBLE_S(correction_sign_1, x_negative);

    const SDOUBLE correction_0 = MUL_DOUBLE_S(b_correction, constants_away);
    const SDOUBLE correction = MUL_DOUBLE_S(correction_0, correction_sign);

    const SDOUBLE one_over_from_behind = DIV_DOUBLE_S(one, from_behind);

    const SDOUBLE correction_term = MUL_DOUBLE_S(correction, one_over_from_behind);

    const SDOUBLE coeff0 = FMADD_PD(correction_term, in_q3, one);
    const SDOUBLE result_q0_1 = FMADD_PD(result_q0_t12, x_square, coeff0);

    const SDOUBLE result_q0 = MUL_DOUBLE_S(result_q0_1, x);
    /* ---- End first Calculation ---- */

    /* ---- Readjusting for the second range ---- */
    const SDOUBLE nominator = MUL_DOUBLE_S(two, result_q0);
    const SDOUBLE result_q0_square = MUL_DOUBLE_S(result_q0, result_q0);
    const SDOUBLE denominator = SUB_DOUBLE_S(one, result_q0_square);

    /* Obtaining the interval results */
    const SDOUBLE result_q1 = DIV_DOUBLE_S(nominator, denominator);
    const SDOUBLE result_q2 = DIV_DOUBLE_S(one, result_q1);
    const SDOUBLE result_q3 = DIV_DOUBLE_S(one, result_q0);

    // Add quadrant results together
    const SDOUBLE partial_result_0 = MUL_DOUBLE_S(result_q0, in_q0);
    const SDOUBLE partial_result_1 = FMADD_PD(result_q1, in_q1, partial_result_0);
    const SDOUBLE partial_result_2 = FMADD_PD(result_q2, in_q2, partial_result_1);

    const SDOUBLE partial_result_31 = FMADD_PD(q2_bitshift, in_q2, partial_result_2);
    const SDOUBLE partial_result_3  = FMADD_PD(result_q3, in_q3, partial_result_31);

    const SDOUBLE result = MUL_DOUBLE_S(sign_adjust, partial_result_3);
    SIMD_TO_DOUBLE_VEC(&res[i], result);
  }

//for (size_t i = 0; i < n; i++) {
//    if (isnan(res[i])) {
//        printf("THERE IS NAN AT: %25.17e", input[i]);
//    }
//}

  int num_left_over = (n % simd_doubles);

  for (size_t i = n - num_left_over; i < n; i++) {
    res[i] = tan(input[i]);
  }
}


// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
// Compilation and run on laptop
//
//
// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 $(pkg-config --cflags --libs flint) -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
