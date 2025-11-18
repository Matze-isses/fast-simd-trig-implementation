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
    double taylor = M_PI_2 - input;
    double x_square = taylor * taylor;

    TAYLOR_COEFF_TAN[0] += CORRECTION * 1/taylor;
    double result = TAYLOR_COEFF_TAN[13];

    for (int j = 12; j >= 0; j-=1) {
      double coeff = TAYLOR_COEFF_TAN[j];
      result = result * x_square + coeff;
    }

    result = result * taylor;
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


void tan_non_simd(double *input, double *res, size_t n, int prec) {
  printf("M_PI_2: %.17g\n", 4 * M_PI_8);
  init_first_mid_lagrange_table();

  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  
  for (int i = 0; i < (int) n; i++) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
    const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);


    /*
    if (input[i] < M_PI_8) {
      start_of_range(input[i], &res[i]);

    } else if (input[i] < 2 * M_PI_8) {
      first_mid_range(input[i], &res[i]);

    } else if (input[i] < 3 * M_PI_8) {
      sec_mid_range(input[i], &res[i]);

    } else if (input[i] > 3 * M_PI_8){
      end_of_range(input[i], &res[i]);
    }
    */
  }
}


