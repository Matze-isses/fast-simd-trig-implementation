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


void speed_test(double *test_values, double *own_results, size_t test_size) {
  uint64_t own_execution_cycles_warmup = 0, own_execution_cycles = 0;
  double own_execution_ms_warmup = 0, own_execution_ms = 0, glibc_execution_ms = 0;
  double own_time_old_clock = -1.0, glibc_time_old_clock = -1.0;
  struct CrystalClock clk;
  clk._freq = frequency();

  struct CrystalClock clk1;
  clk1._freq = frequency();


  /* ---- WARMUP ---- */ 
  START_TCLOCK;
  sin_simd(test_values, own_results, test_size);
  own_execution_ms_warmup = GET_TCLOCK;

  printf("WARMUP Time OC: %.17g\n", own_execution_ms_warmup);


  clk._begin = current();
  sin_simd(test_values, own_results, test_size);
  clk._end   = current();

  own_execution_ms_warmup = duration_ms1(clk);
  printf("WARMUP Time CC: %.17g\n\n", own_execution_ms_warmup);


  /* ---- Own Results ---- */
  START_TCLOCK;
  sin_simd(test_values, own_results, test_size);
  own_time_old_clock = GET_TCLOCK;

  printf("Own   Time OC: %.17g\n", own_time_old_clock);


  clk._begin = current();
  sin_simd(test_values, own_results, test_size);
  clk._end   = current();

  own_execution_ms = duration_ms1(clk);
  printf("Own   Time CC: %.17g\n\n", own_execution_ms);



  /* ---- GLIBC Results ---- */
  // Written last to obtain here the best case for glibc, where mine could still be in warmup
  double *glibc_results = malloc(test_size * sizeof(double));

  START_TCLOCK;
  for (size_t i = 0; i < test_size; i++) { glibc_results[i] = sin(test_values[i]); }
  glibc_time_old_clock = GET_TCLOCK;

  printf("GLIBC Time OC: %.17g\n", glibc_time_old_clock);

  clk1._begin = current();
  for (size_t i = 0; i < test_size; i++) { glibc_results[i] = sin(test_values[i]); }
  clk1._end   = current();

  glibc_execution_ms = duration_ms1(clk1);

  printf("GLIBC Time CC: %.17g\n", glibc_execution_ms);

  free(glibc_results);
}


void quadrant_error_test(size_t n) {

  /* bounds: [0, pi/8], [pi/8, pi/4], [pi/4, 3*pi/8], [3*pi/8, pi/2] */
  double bounds[5] = {
    0.0,
    M_PI / 8.0,
    M_PI / 4.0,
    3.0 * M_PI / 8.0,
    M_PI / 2.0
  };

  const char *interval_names[4] = {
    "[0, pi/8]",
    "[pi/8, pi/4]",
    "[pi/4, 3pi/8]",
    "[3pi/8, pi/2]"
  };

  printf("==============================================================================================================\n");
  printf("Quadrant Error Analysis for n = %d\n", (int)n);
  printf("==============================================================================================================\n\n");

  printf("+----------------+----------------+---------------------------+---------------------------+\n");
  printf("| Interval       | Implementation | Max Error                 | Avg Abs Error             |\n");
  printf("+----------------+----------------+---------------------------+---------------------------+\n");

  for (int q = 0; q < 4; ++q) {
    double lower = bounds[q];
    double upper = bounds[q + 1];

    double *test_values    = malloc((size_t)n * sizeof(double));
    double *own_results    = malloc((size_t)n * sizeof(double));
    double *glibc_results  = malloc((size_t)n * sizeof(double));

    if (!test_values || !own_results || !glibc_results) {
      perror("malloc");
      free(test_values);
      free(own_results);
      free(glibc_results);
      return;
    }

    /* generate random inputs in the given sub-interval */
    fill_uniform(lower, upper, n, test_values);

    /* compute results */
    sin_simd(test_values, own_results, n);
    for (size_t i = 0; i < n; ++i) {
      glibc_results[i] = sin(test_values[i]);
    }

    /* error metrics */
    double abs_error_own, max_error_own, value_max_error_own;
    double abs_error_glibc, max_error_glibc, value_max_error_glibc;

    compare_results_sin(test_values, own_results,
                        &abs_error_own, &max_error_own, &value_max_error_own, n);
    compare_results_sin(test_values, glibc_results,
                        &abs_error_glibc, &max_error_glibc, &value_max_error_glibc, n);

    double avg_error_own   = abs_error_own   / (double)n;
    double avg_error_glibc = abs_error_glibc / (double)n;

    /* print table rows */
    printf("| %-14s | %-14s | %25.17e | %25.17e |\n",
           interval_names[q], "own",   max_error_own,   avg_error_own);
    printf("| %-14s | %-14s | %25.17e | %25.17e |\n",
           "",               "glibc", max_error_glibc, avg_error_glibc);

    printf("+----------------+----------------+---------------------------+---------------------------+\n");

    free(test_values);
    free(own_results);
    free(glibc_results);
  }

  printf("\n");
}


int main(int argc, char *argv[]) {

  if (argc < 3) {
    fprintf(stderr, "Usage: %s n lower upper\n", argv[0]);
    return 1;
  }

  int n = atoi(argv[1]);
  if (n % 4 != 0) { n = n + (4 - (n % 4)); }

  double lower = atof(argv[2]);
  double upper = atof(argv[3]);


  printf("\n Test Setup \n----------------------------------------------------------------------------------------------------\n");
  printf("Input Interval:                       [%f, %f]\n", lower, upper);
  printf("Number of inputs (speed test):        n=%d\n", n);
  printf("----------------------------------------------------------------------------------------------------\n\n");

  size_t test_size = (argc > 4) ? atoi(argv[4]) : (size_t)n;

  /* NEW: special mode – only quadrant error test when lower == upper */
  if (lower == upper) {
    quadrant_error_test(test_size);
  }


  srand((unsigned)time(NULL));

  double *test_values = malloc((size_t)n * sizeof(double));

  double *correct_results = malloc((size_t)n * sizeof(double));
  double *own_results = malloc((size_t)n * sizeof(double));
  double *glibc_results = malloc((size_t)n * sizeof(double));


  if (!test_values) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower, upper, test_size, test_values);
  
  if (n != 0) {
    speed_test(test_values, own_results, n);
  }

  // precision test
  if (test_size > 0 && lower != upper) {
    test_values[0] = upper;

    sin_simd(test_values, own_results, test_size);
    for (size_t i = 0; i < test_size; i++) { glibc_results[i] = sin(test_values[i]); }

    // user information for the current state of the script
    printf("\nThe results are obtained! Starting error calculation.\n");

    double *correct_results_partial = malloc(test_size * sizeof(double));
    double *glibc_results_partial = malloc(test_size * sizeof(double));
    double *own_results_partial = malloc(test_size * sizeof(double));

    double *test_values_partial = malloc(test_size * sizeof(double));

    for (size_t i = 0; i < test_size; i++) {
      correct_results_partial[i] = correct_results[i];
      glibc_results_partial[i] = glibc_results[i];
      own_results_partial[i] = own_results[i];
      test_values_partial[i] = test_values[i];
    }
    double abs_error, abs_error_glibc, max_error, max_error_glibc, value_max_error, value_max_error_glibc;

    compare_results_sin(test_values_partial, own_results_partial, &abs_error, &max_error, &value_max_error, test_size);
    compare_results_sin(test_values_partial, glibc_results_partial, &abs_error_glibc, &max_error_glibc, &value_max_error_glibc, test_size);

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
