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
  1.0000000000000000,
  0.33333333333333331,
  0.13333333333333333,
  0.053968253968253971,
  0.021869488536155203,
  0.0088632355299021973,
  0.0035921280365724811,
  0.0014558343870513183,
  0.00059002744094558595,
  0.00023912911424355248,
  9.6915379569294509e-05,
  3.9278323883316833e-05,
  1.5918905069328964e-05,
  6.4516892156554306e-06
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


void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;
  init_first_mid_lagrange_table();

  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE correction = LOAD_DOUBLE(CORRECTION);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = 12;
  double test_vec[4] = {0.0, 1.0, 2.0, 3.0};

  const SDOUBLE half = LOAD_DOUBLE(0.5);
  const SDOUBLE one = LOAD_DOUBLE(1.0);
  const SDOUBLE two = LOAD_DOUBLE(2.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);

    const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
    const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);
    // const SDOUBLE quadrant = LOAD_DOUBLE_VEC(test_vec);

    SDOUBLE in_q0 = SUB_DOUBLE_S(quadrant, two);
    in_q0 = ABS_PD(in_q0);
    in_q0 = MUL_DOUBLE_S(in_q0, half);
    in_q0 = FLOOR_DOUBLE_S(in_q0);

    SDOUBLE in_q1 = SUB_DOUBLE_S(quadrant, one);
    in_q1 = ABS_PD(in_q1);
    in_q1 = SUB_DOUBLE_S(in_q1, two);
    in_q1 = ABS_PD(in_q1);
    in_q1 = MUL_DOUBLE_S(in_q1, half);
    in_q1 = FLOOR_DOUBLE_S(in_q1);

    SDOUBLE in_q2 = SUB_DOUBLE_S(quadrant, two);
    in_q2 = ABS_PD(in_q2);
    in_q2 = SUB_DOUBLE_S(in_q2, two);
    in_q2 = ABS_PD(in_q2);
    in_q2 = MUL_DOUBLE_S(in_q2, half);
    in_q2 = FLOOR_DOUBLE_S(in_q2);

    SDOUBLE in_q3 = SUB_DOUBLE_S(quadrant, one);
    in_q3 = ABS_PD(in_q3);
    in_q3 = MUL_DOUBLE_S(in_q3, half);
    in_q3 = FLOOR_DOUBLE_S(in_q3);

    /*
    PRINT_M256D(in_q0);
    PRINT_M256D(in_q1);
    PRINT_M256D(in_q2);
    PRINT_M256D(in_q3);
    */


    SDOUBLE result = LOAD_DOUBLE(0.0);
    
    /* ---- Calculation for first range ---- */
    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);
    SDOUBLE result_q0 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    for (int j = taylor_loop_iteration; j >= 0; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q0 = FMADD_PD(result_q0, x_square, coeff);
    }

    result_q0 = MUL_DOUBLE_S(result_q0, x);


    /* ---- Calculation for second range ---- */
    SDOUBLE result_q1 = LOAD_DOUBLE(0.0);

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q1 = ADD_DOUBLE_S(result_q1, eval_inner);
    }

    /* ---- Calculation for thierd range ---- */
    SDOUBLE result_q2 = LOAD_DOUBLE(0.0);

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(from_behind, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(from_behind, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q2 = ADD_DOUBLE_S(result_q2, eval_inner);
    }

    result_q2 = DIV_DOUBLE_S(one, result_q2);
    
    /* ---- Calculation for fourth range ---- */
    SDOUBLE from_behind_square = MUL_DOUBLE_S(from_behind, from_behind);

    // Calculation of the correction term
    SDOUBLE one_over_from_behind = DIV_DOUBLE_S(one, from_behind);
    SDOUBLE correction_term = MUL_DOUBLE_S(correction, one_over_from_behind);
    SDOUBLE first_coeff = ADD_DOUBLE_S(one, correction_term);

    SDOUBLE result_q3 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    // Important that first term is treated with correction
    for (int j = taylor_loop_iteration; j >= 1; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q3 = FMADD_PD(result_q3, from_behind_square, coeff);
    }

    // add the corrected first coefficiant
    result_q3 = FMADD_PD(result_q3, from_behind_square, first_coeff);
    result_q3 = MUL_DOUBLE_S(result_q3, from_behind);
    result_q3 = DIV_DOUBLE_S(one, result_q3);


    /* ---- Final Addup ---- */

    result_q1 = REMOVE_INF(result_q1);
    result_q2 = REMOVE_INF(result_q2);

    /*
    PRINT_M256D(result_q0);
    PRINT_M256D(result_q1);
    PRINT_M256D(result_q2);
    PRINT_M256D(result_q3);
    */

    result = FMADD_PD(result_q0, in_q0, result);
    result = FMADD_PD(result_q1, in_q1, result);
    result = FMADD_PD(result_q2, in_q2, result);
    result = FMADD_PD(result_q3, in_q3, result);


    SIMD_TO_DOUBLE_VEC(&res[i], result);

  }

  int num_left_over = (n % 4);

  for (size_t i = n - num_left_over; i < (int)n; i++) {
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

  for (int i = 0; i < n; i++) {
      if (isnan(res[i])) {
        printf("NaN at index %d, input = %.17g\n", i, input[i]);
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
