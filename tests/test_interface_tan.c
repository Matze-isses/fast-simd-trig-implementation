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

  const bool eval_glibc = true;

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
  test_values[0] = upper;
  
  // user information for the current state of the script
  printf("Test values are generated! Starting calculation of correct results.\n");

  uint64_t own_execution_cycles_warmup = 0, own_execution_cycles = 0;
  double own_execution_ms_warmup = 0, own_execution_ms = 0;
  double own_time_old_clock = -1.0;
  struct CrystalClock clk;
  clk._freq = frequency();


  if (test_size == 0) {
    printf("\n -------- Own Script WARMUP Execution ---------- \n\n");
    clk._begin = current();

    tan_simd(test_values, own_results, n, 20);

    clk._end   = current();
    printf("\n ------- End Own Script WARMUP Execution ------- \n\n");

    own_execution_ms = duration_ms1(clk);
    printf("WARMUP time CC: %.17g\n", own_execution_ms);


    printf("\n -------- Own Script OC Execution ---------- \n\n");
    START_TCLOCK;

    tan_simd(test_values, own_results, n, 20);

    own_time_old_clock = GET_TCLOCK;
    printf("\n ------- End Own Script OC Execution ------- \n\n");

  }

  printf("\n -------- Own Script Execution ---------- \n\n");
  clk._begin = current();

  tan_simd(test_values, own_results, n, 20);

  clk._end   = current();
  printf("\n ------- End Own Script Execution ------- \n\n");

  own_execution_ms = duration_ms1(clk);
  printf("OWN TIME OC: %.17g\n\n", own_time_old_clock);
  printf("OWN TIME CC: %.17g\n\n", own_execution_ms);

  // Evaluate glibc if local constant is set
  if (eval_glibc) {
    clk._begin = current();

    for (int i = 0; i < n; i++) { glibc_results[i] = tan(test_values[i]); }

    clk._end   = current();
    own_execution_ms = duration_ms1(clk);
    printf("GLIBC TIME CC: %.17g\n", own_execution_ms);
  }

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
    double abs_error, abs_error_glibc, max_error, max_error_glibc, value_max_error, value_max_error_glibc;

    compare_results_tan(test_values_partial, own_results_partial, &abs_error, &max_error, &value_max_error, test_size);
    compare_results_tan(test_values_partial, glibc_results_partial, &abs_error_glibc, &max_error_glibc, &value_max_error_glibc, test_size);

    printf("Max Error Own:   %.17g;   At Value %.17g\n", max_error, value_max_error);
    printf("Max Error glibc: %.17g;   At Value %.17g\n\n", max_error_glibc, value_max_error_glibc);

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
