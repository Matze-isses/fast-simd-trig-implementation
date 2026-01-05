#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
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
  int simd_doubles = 4;

  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE correction = LOAD_DOUBLE(TAN_CORRECTION);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = last_taylor_coeff - 1;

  const SDOUBLE neg_half = LOAD_DOUBLE(-0.5);
  const SDOUBLE half = LOAD_DOUBLE(0.5);
  const SDOUBLE one = LOAD_DOUBLE(1.0);
  const SDOUBLE two = LOAD_DOUBLE(2.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE result = LOAD_DOUBLE(0.0);
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    //TODO: RANGE REDUCTION MISSING

    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);

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

    SDOUBLE not_in_q3 = SUB_DOUBLE_S(one, in_q3);

    SDOUBLE one_over_from_behind = DIV_DOUBLE_S(one, from_behind);
    SDOUBLE correction_term = MUL_DOUBLE_S(correction, one_over_from_behind);
    SDOUBLE adjusted_first = ADD_DOUBLE_S(one, correction_term);
    SDOUBLE first_coeff = FMADD_PD(adjusted_first, in_q3, not_in_q3);
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

    result = FMADD_PD(result_q0, in_q0, result);
    result = FMADD_PD(result_q1, in_q1, result);
    result = FMADD_PD(result_q2, in_q2, result);
    result = FMADD_PD(result_q3, in_q3, result);


    SIMD_TO_DOUBLE_VEC(&res[i], result);

  }

  int num_left_over = (n % 4);

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
