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
double M_PI_8 = M_PI /8;

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
  
}


void tan_simd(double *input, double *res, size_t n, int prec) {
  printf("M_PI_2: %.17g\n", M_PI_2/4);
  
  for (int i = 0; i < (int) n; i++) {
    if (input[i] < M_PI_8) {
      start_of_range(input[i], &res[i]);
    }
    else if (input[i] > 3 * M_PI_8){
      end_of_range(input[i], &res[i]);
    }
  }

}

// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/trig_arb_comparison.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 1000000000 -8 8 1000000
//
// Compilation and run on laptop
//
//
//gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/trig_arb_comparison.c -o test -lm -mavx -mavx2 -mfma -O2 -Wextra $(pkg-config --cflags --libs flint)  && ./test 1000000 -8 8 100000


//
