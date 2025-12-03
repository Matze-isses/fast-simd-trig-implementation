#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include "trig_simd.h"

#include <string.h>   // just used for the error calculation
                      
#include "./util/bit_printing.h"


double CORRECTION = 0.00000000000000006123233995736765;
double M_PI_8 = M_PI / 8;

// Working until 0.4
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
  6.451689215655431e-06,  /* 6.451689215655431e-06 */
  2.6147711512907546e-06, /* 2.6147711512907546e-06 */ // <-- Current last
  1.0597268320104656e-06,
  4.2949110782738063e-07,
  1.7406618963571645e-07,
  7.054636946400967e-08,
  2.859136662305253e-08,
  1.1587644432798851e-08,
  4.696295398230901e-09,
  1.9033368339312746e-09,
  7.713933635359061e-10,
  3.1263395458920874e-10
};


#define N_FIRST_MID 17

static const double X_FIRST_MID[17] = {
    0.39269908169872414,
    0.3974311529539063,
    0.4054030240311012,
    0.4286863854186324,
    0.4576011079513348,
    0.487426901464288,
    0.5407991673806207,
    0.5723175108611952,
    0.6115367603429416,
    0.6435891799153872,
    0.6768825409954158,
    0.7002130267439945,
    0.7267833123145859,
    0.7490809891502078,
    0.7665916163094255,
    0.7796607724399203,
    0.7853981633974483
};

static const double Y_FIRST_MID[17] = {
    0.41421356237309503,
    0.4197684582474178,
    0.42917669979674156,
    0.4570320684430114,
    0.4924645633683938,
    0.5300875167665825,
    0.6005164649453173,
    0.6442430603066724,
    0.7012087722833138,
    0.7501376202189168,
    0.8035182551172894,
    0.8426526043282521,
    0.8891412771435833,
    0.9298813332611141,
    0.9630769477832821,
    0.9885905533499666,
    1.0000000000000000
};

static double LAGRANGE_DEN_FIRST_MID[N_FIRST_MID][N_FIRST_MID];

void init_first_mid_lagrange_table(void)
{
    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        for (size_t j = 0; j < N_FIRST_MID; ++j) {
            if (i == j) {
                LAGRANGE_DEN_FIRST_MID[i][j] = 0.0;
            } else {
                LAGRANGE_DEN_FIRST_MID[i][j] = 1.0 / (X_FIRST_MID[i] - X_FIRST_MID[j]);
            }
        }
    }
}

void start_of_range(double input, double *res) {
    double taylor = input;
    double x_square = taylor * taylor;
    double result = TAYLOR_COEFF_TAN[13];

    for (int j = 12; j >= 0; j-=1) {
      double coeff = TAYLOR_COEFF_TAN[j];
      result = result * x_square + coeff;
    }

    *res = result * taylor;

}

void end_of_range(double input, double *res) {
    double from_behind = M_PI_2 - input;
    double x_square = from_behind * from_behind;

    TAYLOR_COEFF_TAN[0] += CORRECTION * 1/from_behind;
    double result = TAYLOR_COEFF_TAN[13];

    for (int j = 12; j >= 0; j-=1) {
      double coeff = TAYLOR_COEFF_TAN[j];
      result = result * x_square + coeff;
    }

    result = result * from_behind;
    *res = 1 / result;
    TAYLOR_COEFF_TAN[0] = 1.0;
}

void first_mid_range(double input, double *res) {
    double result = 0.0;

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        double Li = 1.0;

        for (size_t j = 0; j < N_FIRST_MID; ++j) {
            if (j == i) continue;
            Li *= (input - X_FIRST_MID[j]) * LAGRANGE_DEN_FIRST_MID[i][j];
        }

        result += Y_FIRST_MID[i] * Li;
    }

    *res = result;
}

void sec_mid_range(double input, double *res) {
    double result = 0.0;
    double from_end = M_PI_2 - input;

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        double Li = 1.0;

        for (size_t j = 0; j < N_FIRST_MID; ++j) {
            if (j == i) continue;
            Li *= (from_end - X_FIRST_MID[j]) * LAGRANGE_DEN_FIRST_MID[i][j];
        }

        result += Y_FIRST_MID[i] * Li;
    }

    *res = 1/result;
}


void tan_simd(double *input, double *res, size_t n) {
  int simd_doubles = 4;

  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE correction = LOAD_DOUBLE(CORRECTION);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = last_taylor_coeff - 1;

  const SDOUBLE neg_half = LOAD_DOUBLE(-0.5);
  const SDOUBLE half = LOAD_DOUBLE(0.5);
  const SDOUBLE one = LOAD_DOUBLE(1.0);
  const SDOUBLE two = LOAD_DOUBLE(2.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE result = LOAD_DOUBLE(0.0);
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);
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
   if (input[i] < M_PI_8) {
      start_of_range(input[i], &res[i]);

    } else if (input[i] < 2 * M_PI_8) {
      first_mid_range(input[i], &res[i]);

    } else if (input[i] < 3 * M_PI_8) {
      sec_mid_range(input[i], &res[i]);

    } else if (input[i] > 3 * M_PI_8){
      end_of_range(input[i], &res[i]);
    }
  }

}


// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
// Compilation and run on laptop
//
//
// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 $(pkg-config --cflags --libs flint) -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
