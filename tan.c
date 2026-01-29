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
  1.000000000000000,
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
  //printf("%.17g %.17g", nextafter(M_PI/2.0, INFINITY), nextafter(3.0 * M_PI/2.0, -INFINITY));
  int simd_doubles = SIMD_LENGTH / 64;

  const SDOUBLE pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  const SDOUBLE one_over_pi_2 = LOAD_DOUBLE(1/M_PI_2);

  //const SDOUBLE correction = LOAD_DOUBLE(TAN_CORRECTION);
  //const SDOUBLE range_reduction_correction = LOAD_DOUBLE(RANG_REDUCTION_CORRECTION);

  // -1 => -0.0000000000000000612323399570
  //  0 =>  0.0000000000000000612323399570
  //  1 =>  0.0000000000000001836970199000
  //  2 =>  0.0000000000000003061616998040
  //
  // 
  // 3/2 pi => 1 => -0.0000000000000000995799000000000
  // 5/2 pi => 2 =>  0.0000000000000000114423776000000
  // 
  //
  // 0.000000000000000122464679943 * num_ranges_away + 0.0000000000000000612323399570
  //
  double q2_bitshift_double = -pow(2, -52);
  const SDOUBLE q2_bitshift = LOAD_DOUBLE(q2_bitshift_double);


  double tan_correction_a = 0.0000000000000001224646799430000;
  // double tan_correction_a =-0.000000000000000244931000000000;
  double tan_correction_b = 0.00000000000000006123233995736766;

  //tan_correction_b = 0.0;
  //tan_correction_a = 0.0;

  const SDOUBLE a_correction = LOAD_DOUBLE(tan_correction_a);
  const SDOUBLE b_correction = LOAD_DOUBLE(tan_correction_b);
  const SDOUBLE range_reduction_correction = LOAD_DOUBLE(0.0);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = last_taylor_coeff - 1;

  const SDOUBLE zero = SET_ZERO();

  const SDOUBLE neg_one = LOAD_DOUBLE(-1.0);
  const SDOUBLE neg_half = LOAD_DOUBLE(-0.5);
  const SDOUBLE half = LOAD_DOUBLE(0.5);
  const SDOUBLE one = LOAD_DOUBLE(1.0);
  const SDOUBLE two = LOAD_DOUBLE(2.0);
  const SDOUBLE three = LOAD_DOUBLE(3.0);
  
  for (int i = 0; i < (int) n; i += simd_doubles) {
    SDOUBLE result = LOAD_DOUBLE(0.0);
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    SDOUBLE x_negative = DIV_DOUBLE_S(x, ABS_PD(x));
    SDOUBLE x_negative01 = MUL_DOUBLE_S(x_negative, neg_one);
    x_negative01 = ADD_DOUBLE_S(x_negative01, one);
    x_negative01 = MUL_DOUBLE_S(x_negative01, half);

    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_pi_2);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, pi_2);
    SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple);

    // here problems occure because negative singularities need a ceil not a floor
    const SDOUBLE singularities_away0 = ADD_DOUBLE_S(num_ranges_away, x_negative01);
    const SDOUBLE singularities_away1 = ABS_PD(singularities_away0);
    const SDOUBLE singularities_away = MUL_DOUBLE_S(singularities_away1, half);
    const SDOUBLE num_singularities_away = FLOOR_DOUBLE_S(singularities_away);
    SDOUBLE constants_away = MUL_DOUBLE_S(num_singularities_away, two);
    constants_away = ADD_DOUBLE_S(constants_away, one);

    SDOUBLE range_reduction_correction_term = MUL_DOUBLE_S(x, range_reduction_correction);
    x = SUB_DOUBLE_S(in_outer_range, range_reduction_correction_term);

    // Check if even
    //  Default Range Reduction
    const SDOUBLE sign_adjust0 = MUL_DOUBLE_S(num_ranges_away, half);
    const SDOUBLE sign_adjust1 = FLOOR_DOUBLE_S(sign_adjust0);
    const SDOUBLE sign_adjust2 = MUL_DOUBLE_S(sign_adjust1, two);
    const SDOUBLE in_odd_range = SUB_DOUBLE_S(num_ranges_away, sign_adjust2);
    const SDOUBLE in_even_range = SUB_DOUBLE_S(one, in_odd_range);
    const SDOUBLE sign_adjust4 = MUL_DOUBLE_S(in_odd_range, two);
    const SDOUBLE sign_adjust = SUB_DOUBLE_S(one, sign_adjust4);
    // PRINT_M256D(sign_adjust);

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
    
    /* ---- Calculation for first range ---- */
    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);
    SDOUBLE result_q0 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    for (int j = taylor_loop_iteration; j >= 1; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q0 = FMADD_PD(result_q0, x_square, coeff);
    }

    SDOUBLE c_sign0 = SUB_DOUBLE_S(in_even_range, in_odd_range);
    SDOUBLE c_sign = MUL_DOUBLE_S(c_sign0, x_negative);
    
    SDOUBLE correction = MUL_DOUBLE_S(b_correction, constants_away);
    correction = MUL_DOUBLE_S(correction, c_sign);

//  PRINT_M256D(constants_away);
//  PRINT_M256D(c_sign);
//  PRINT_FULL_M256D(correction);
    // PRINT_M256D(num_ranges_away);
    // PRINT_FULL_M256D(correction);

    SDOUBLE one_over_from_behind = DIV_DOUBLE_S(one, from_behind);
    SDOUBLE correction_term = MUL_DOUBLE_S(correction, one_over_from_behind);
    SDOUBLE first_coeff = FMADD_PD(correction_term, in_q3, one);

    result_q0 = FMADD_PD(result_q0, x_square, first_coeff);
    result_q0 = MUL_DOUBLE_S(result_q0, x);
    /* ---- End first Calculation ---- */

    /* ---- Readjusting for the second range ---- */
    SDOUBLE nominator = MUL_DOUBLE_S(two, result_q0);
    SDOUBLE result_q0_square = MUL_DOUBLE_S(result_q0, result_q0);
    SDOUBLE denominator = SUB_DOUBLE_S(one, result_q0_square);
    SDOUBLE result_q1 = DIV_DOUBLE_S(nominator, denominator);

    /* ---- Calculation for thierd range ---- */
    SDOUBLE result_q2 = DIV_DOUBLE_S(one, result_q1);
    
    /* ---- Calculation for fourth range ---- */
    SDOUBLE result_q3 = DIV_DOUBLE_S(one, result_q0);

    // Add up of the results
    result = FMADD_PD(result_q0, in_q0, result);
    result = FMADD_PD(result_q1, in_q1, result);
    result = FMADD_PD(result_q2, in_q2, result);

    // ULP of q2 is between +5 and -3 shifting one up  gives +4 -4 => smaller abs ULP
    result = FMADD_PD(q2_bitshift, in_q2, result);

    result = FMADD_PD(result_q3, in_q3, result);
    
    result = MUL_DOUBLE_S(sign_adjust, result);
    SIMD_TO_DOUBLE_VEC(&res[i], result);
  }

//for (size_t i = 0; i < n; i++) {
//  if (M_PI_4 < input[i] && input[i] < 3.0 * M_PI_8) {
//    res[i] = nextafter(res[i], -INFINITY);
//  }
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
