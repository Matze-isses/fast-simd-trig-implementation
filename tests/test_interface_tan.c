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


  /* ---- WARMUP ---- */ 
  START_TCLOCK;
  tan_simd(test_values, own_results, test_size);
  own_execution_ms_warmup = GET_TCLOCK;

  printf("WARMUP Time OC: %.17g\n", own_execution_ms_warmup);


  clk._begin = current();
  tan_simd(test_values, own_results, test_size);
  clk._end   = current();

  own_execution_ms_warmup = duration_ms1(clk);
  printf("WARMUP Time CC: %.17g\n\n", own_execution_ms_warmup);


  /* ---- Own Results ---- */
  START_TCLOCK;
  tan_simd(test_values, own_results, test_size);
  own_time_old_clock = GET_TCLOCK;

  printf("Own   Time OC: %.17g\n", own_time_old_clock);


  clk._begin = current();
  tan_simd(test_values, own_results, test_size);
  clk._end   = current();

  own_execution_ms = duration_ms1(clk);
  printf("Own   Time CC: %.17g\n\n", own_execution_ms);



  /* ---- GLIBC Results ---- */
  // Written last to obtain here the best case for glibc, where mine could still be in warmup
  double *glibc_results = malloc(test_size * sizeof(double));

  START_TCLOCK;
  for (size_t i = 0; i < test_size; i++) { glibc_results[i] = tan(test_values[i]); }
  glibc_time_old_clock = GET_TCLOCK;

  printf("GLIBC Time OC: %.17g\n", glibc_time_old_clock);

  clk._begin = current();
  for (size_t i = 0; i < test_size; i++) { glibc_results[i] = tan(test_values[i]); }
  clk._end   = current();

  glibc_execution_ms = duration_ms1(clk);
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
    tan_simd(test_values, own_results, n);
    for (size_t i = 0; i < n; ++i) {
      glibc_results[i] = tan(test_values[i]);
    }

    /* error metrics */
    double abs_error_own, max_error_own, value_max_error_own;
    double abs_error_glibc, max_error_glibc, value_max_error_glibc;

    compare_results_tan(test_values, own_results,
                        &abs_error_own, &max_error_own, &value_max_error_own, n);
    compare_results_tan(test_values, glibc_results,
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

  printf("\033[1A\033[2K\033[1A\033[2K\033[1A\033[2K");

  if (argc < 3) {
    fprintf(stderr, "Usage: %s n lower upper\n", argv[0]);
    return 1;
  }

  int speed_test_size = atoi(argv[1]);
  if (speed_test_size % 4 != 0) { speed_test_size = speed_test_size + (4 - (speed_test_size % 4)); }

  double lower = atof(argv[2]);
  double upper = atof(argv[3]);
  size_t accuracy_test_size = (argc > 4) ? atoi(argv[4]) : speed_test_size;


  printf("\n Test Setup \n----------------------------------------------------------------------------------------------------\n");
  printf("Input Interval:                       [%f, %f]\n", lower, upper);
  printf("Number of inputs (speed test):        n=%d\n", speed_test_size);
  printf("Number of inputs (error calculation): n=%d\n", (int)accuracy_test_size);
  printf("----------------------------------------------------------------------------------------------------\n\n");

  if (accuracy_test_size > 0) {
    quadrant_error_test(accuracy_test_size);
  }

  srand((unsigned)time(NULL));



  
  if (speed_test_size != 0) {
    double *speed_test_values = malloc(speed_test_size * sizeof(double));
    double *own_speed_results = malloc(speed_test_size * sizeof(double));

    fill_uniform(lower, upper, speed_test_size, speed_test_values);
    speed_test(speed_test_values, own_speed_results, speed_test_size);

    free(speed_test_values);
    free(own_speed_results);
  }


  // precision test
  if (accuracy_test_size > 0 && speed_test_size == 0) {
    double *test_values = malloc(accuracy_test_size * sizeof(double));

    double *correct_results = malloc(accuracy_test_size * sizeof(double));
    double *own_results = malloc(accuracy_test_size * sizeof(double));
    double *glibc_results = malloc(accuracy_test_size * sizeof(double));

    fill_uniform(lower, upper, accuracy_test_size, test_values);
    test_values[0] = upper;

    tan_simd(test_values, own_results, accuracy_test_size);
    for (size_t i = 0; i < accuracy_test_size; i++) { glibc_results[i] = tan(test_values[i]); }

    // user information for the current state of the script
    printf("\nThe results are obtained! Starting error calculation.\n");
    double abs_error, abs_error_glibc, max_error, max_error_glibc, value_max_error, value_max_error_glibc;

    compare_results_tan(test_values, own_results, &abs_error, &max_error, &value_max_error, accuracy_test_size);
    compare_results_tan(test_values, glibc_results, &abs_error_glibc, &max_error_glibc, &value_max_error_glibc, accuracy_test_size);

    printf("Max Error Own:   %.17g;   At Value %.17g\n", max_error, value_max_error);
    printf("Max Error glibc: %.17g;   At Value %.17g\n\n", max_error_glibc, value_max_error_glibc);

    printf("Accumulated Absolut Error Own   Results: %.17g\n", abs_error);
    printf("Accumulated Absolut Error glibc Results: %.17g\n", abs_error_glibc);

    printf("\nAbsolut Error Own   Results: %.17g\n", abs_error/accuracy_test_size);
    printf("Absolut Error glibc Results: %.17g\n", abs_error_glibc/accuracy_test_size);


    free(test_values);
    free(correct_results);
    free(glibc_results);
    free(own_results);
  }

  return 0;
}

// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
// Compilation and run on laptop
//
//
// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 $(pkg-config --cflags --libs flint) -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
