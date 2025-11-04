#include "clock_utils.h"
#include <omp.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "test_interface.h"
#include "../trig_simd.h"
#include "../cmeasure/cmeasure.h"
#include "../cmeasure/CrystalClockInC.h"

int main(int argc, char *argv[]) {

  if (argc < 3) {
    fprintf(stderr, "Usage: %s n lower upper\n", argv[0]);
    return 1;
  }

  int n = atoi(argv[1]);
  double lower = atof(argv[2]);
  double upper = atof(argv[3]);

  bool eval_glibc = false;

  int test_size = (argc > 4) ? atoi(argv[4]) : n;

  printf("\n Test Setup \n----------------------------------------------------------------------------------------------------\n");
  printf("Input Interval:                       [%f, %f]\n", lower, upper);
  printf("Number of inputs (speed test):        n=%d\n", n);
  printf("Number of inputs (error calculation): n=%d\n", test_size);
  printf("----------------------------------------------------------------------------------------------------\n\n");

  srand((unsigned)time(NULL));

  double *test_values = malloc(n * sizeof(double));

  double *correct_results = malloc(n * sizeof(double));
  double *own_results = malloc(n * sizeof(double));
  double *glibc_results = malloc(n * sizeof(double));


  if (!test_values) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower, upper, n, test_values);
  
  // user information for the current state of the script
  printf("Test values are generated! Starting calculation of correct results.\n");

  // Cristal Clock Setup
  uint64_t own_execution_cycles_warmup = 0, own_execution_cycles = 0;
  double own_execution_ms_warmup = 0, own_execution_ms = 0;
  struct CrystalClock clk;

  printf("\n -------- Own Script WARMUP Execution ---------- \n\n");
  clk._begin = current();

  sin_simd(test_values, own_results, n, 20);

  clk._end   = current();
  printf("\n ------- End Own Script WARMUP Execution ------- \n");

  own_execution_cycles_warmup = own_execution_cycles_warmup + cycles1(clk);
  own_execution_ms_warmup = own_execution_ms_warmup + duration_ms1(clk);
  printf("WARMUP Execution Time Cristal Clock MS: %.17g\n", own_execution_ms_warmup);


  START_TCLOCK;
  sin_simd(test_values, own_results, n, 20);
  double own_time = GET_TCLOCK;


  clk._begin = current();
  sin_simd(test_values, own_results, n, 20);
  clk._end   = current();

  own_execution_ms = duration_ms1(clk);
  printf("TIME OC: %.17g\n", own_time);
  printf("TIME CC: %.17g\n", own_execution_ms);

  if (eval_glibc) {
    START_CLOCK;
    for (int i = 0; i < n; i++) { glibc_results[i] = sin(test_values[i]); }
    END_CLOCK("Time needed by glibc (second)               ");
  }


  /*
  double total_time = 0.0;
  for (int i = 0; i < 10; i++) {
    START_TCLOCK;
    sin_simd(test_values, own_results, n, 0.1);
    total_time += GET_TCLOCK;    
  }
  printf("\n10 Evaluations for the own implementation took on average %.17g ms\n", total_time/10);
  */

  // precision test
  if (test_size > 0) {
    // user information for the current state of the script
    printf("\nThe results are obtained! Starting error calculation.\n");

    double *correct_results_partial = malloc(test_size * sizeof(double));
    double *glibc_results_partial = malloc(test_size * sizeof(double));
    double *own_results_partial = malloc(test_size * sizeof(double));

    double *test_values_partial = malloc(test_size * sizeof(double));

    for (int i = 0; i < test_size; i++) {
      correct_results_partial[i] = correct_results[i];
      glibc_results_partial[i] = glibc_results[i];
      own_results_partial[i] = own_results[i];
      test_values_partial[i] = test_values[i];
    }

    double abs_error = compare_results_sin(test_values_partial, own_results_partial, test_size);
    double abs_error_glibc = compare_results_sin(test_values_partial, glibc_results_partial, test_size);

    printf("Accumulated Absolut Error Own   Results: %.17g\n", abs_error);
    printf("Accumulated Absolut Error glibc Results: %.17g\n", abs_error_glibc);

    printf("\nAbsolut Error Own   Results: %.17g\n", abs_error/test_size);
    printf("Absolut Error glibc Results: %.17g\n", abs_error_glibc/test_size);

    free(correct_results_partial);
    free(glibc_results_partial);
    free(own_results_partial);
    free(test_values_partial);
  }

  free(test_values);
  free(correct_results);
  free(glibc_results);
  free(own_results);


  return 0;
}

// gcc ./tests/test_interface_sin.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./trig_simd.c ./tests/test_object.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 12 0 100 0
