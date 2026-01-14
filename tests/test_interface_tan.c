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


static void run_accuracy_test(size_t accuracy_test_size) {
  if (accuracy_test_size > 0) {
    quadrant_error_test(accuracy_test_size);
  }
}

static void run_speed_test(double lower, double upper, int speed_test_size) {
  if (speed_test_size == 0) return;

  double *speed_test_values = (double *)malloc((size_t)speed_test_size * sizeof(double));
  double *own_speed_results = (double *)malloc((size_t)speed_test_size * sizeof(double));

  if (!speed_test_values || !own_speed_results) {
    fprintf(stderr, "Allocation failed in run_speed_test.\n");
    free(speed_test_values);
    free(own_speed_results);
    return;
  }

  fill_uniform(lower, upper, (size_t)speed_test_size, speed_test_values);
  speed_test(speed_test_values, own_speed_results, speed_test_size);

  free(speed_test_values);
  free(own_speed_results);
}

static void run_precision_test(double lower, double upper,
                               size_t accuracy_test_size,
                               int speed_test_size) {
  double *test_values    = (double *)malloc(accuracy_test_size * sizeof(double));
  double *correct_results = (double *)malloc(accuracy_test_size * sizeof(double)); /* kept to match your original */
  double *own_results    = (double *)malloc(accuracy_test_size * sizeof(double));
  double *glibc_results  = (double *)malloc(accuracy_test_size * sizeof(double));

  if (!test_values || !correct_results || !own_results || !glibc_results) {
    fprintf(stderr, "Allocation failed in run_precision_test.\n");
    free(test_values);
    free(correct_results);
    free(own_results);
    free(glibc_results);
    return;
  }

  fill_uniform(lower, upper, accuracy_test_size, test_values);
  test_values[0] = upper;

  tan_simd(test_values, own_results, accuracy_test_size);
  for (size_t i = 0; i < accuracy_test_size; i++) {
    glibc_results[i] = tan(test_values[i]);
  }

  printf("\nThe results are obtained! Starting error calculation.\n");

  double abs_error = 0.0, abs_error_glibc = 0.0;
  double max_error = 0.0, max_error_glibc = 0.0;
  double value_max_error = 0.0, value_max_error_glibc = 0.0;

  compare_results_tan(test_values, own_results,
                      &abs_error, &max_error, &value_max_error,
                      accuracy_test_size);

  compare_results_tan(test_values, glibc_results,
                      &abs_error_glibc, &max_error_glibc, &value_max_error_glibc,
                      accuracy_test_size);

  printf("Max Error Own:   %.17g;   At Value %.17g\n", max_error, value_max_error);
  printf("Max Error glibc: %.17g;   At Value %.17g\n\n", max_error_glibc, value_max_error_glibc);

  printf("Accumulated Absolut Error Own   Results: %.17g\n", abs_error);
  printf("Accumulated Absolut Error glibc Results: %.17g\n", abs_error_glibc);

  printf("\nAbsolut Error Own   Results: %.17g\n", abs_error / (double)accuracy_test_size);
  printf("Absolut Error glibc Results: %.17g\n", abs_error_glibc / (double)accuracy_test_size);

  free(test_values);
  free(correct_results);
  free(own_results);
  free(glibc_results);
}


static void plot_error_behavior(double lower, double upper, size_t accuracy_test_size) {
  if (accuracy_test_size == 0) return;

  double *test_values = (double *)malloc(accuracy_test_size * sizeof(double));
  double *own_results = (double *)malloc(accuracy_test_size * sizeof(double));
  double *err         = (double *)malloc(accuracy_test_size * sizeof(double));

  if (!test_values || !own_results || !err) {
    fprintf(stderr, "Allocation failed in plot_error_behavior.\n");
    free(test_values);
    free(own_results);
    free(err);
    return;
  }

  /* linspace(lower, upper, accuracy_test_size) with inclusive endpoints */
  if (accuracy_test_size == 1) {
    test_values[0] = lower;
  } else {
    double step = (upper - lower) / (double)(accuracy_test_size - 1);
    for (size_t i = 0; i < accuracy_test_size; i++) {
      test_values[i] = lower + step * (double)i;
    }
    /* ensure exact endpoints (numerical hygiene) */
    test_values[0] = lower;
    test_values[accuracy_test_size - 1] = upper;
  }

  /* run your implementation */
  tan_simd(test_values, own_results, accuracy_test_size);

  /* compute absolute errors: err[i] = |tan(test_values[i]) - own_results[i]| */
  compare_results_tan_err(test_values, own_results, err, accuracy_test_size);

  /* write data for python (simple whitespace-separated columns: x err) */
  const char *data_path = "tan_error_behavior.tsv";
  FILE *f = fopen(data_path, "w");
  if (!f) {
    fprintf(stderr, "Failed to open %s for writing.\n", data_path);
    free(test_values);
    free(own_results);
    free(err);
    return;
  }

  /* header (optional, python-friendly) */
  fprintf(f, "x\terr\n");

  for (size_t i = 0; i < accuracy_test_size; i++) {
    fprintf(f, "%.17g\t%.17g\n", test_values[i], err[i]);
  }
  fclose(f);

  printf("Wrote data: %s\n", data_path);

  free(test_values);
  free(own_results);
  free(err);
}


int main(int argc, char *argv[]) {

  if (argc < 4) {
    fprintf(stderr, "Usage: %s n lower upper [accuracy_test_size]\n", argv[0]);
    return 1;
  }

  int speed_test_size = atoi(argv[1]);
  if (speed_test_size % 4 != 0) {
    speed_test_size = speed_test_size + (4 - (speed_test_size % 4));
  }

  double lower = atof(argv[2]);
  double upper = atof(argv[3]);
  size_t accuracy_test_size = (argc > 4) ? (size_t)atoi(argv[4]) : (size_t)speed_test_size;

  printf("\n----------------------------------------------------------------------------------------------------\n");
  printf(" Test Setup \n");
  printf("----------------------------------------------------------------------------------------------------\n");
  printf("Input Interval:                       [%f, %f]\n", lower, upper);
  printf("Number of inputs (speed test):        n=%d\n", speed_test_size);
  printf("Number of inputs (error calculation): n=%d\n", (int)accuracy_test_size);
  printf("----------------------------------------------------------------------------------------------------\n\n");

  //run_accuracy_test(accuracy_test_size);

  srand((unsigned)time(NULL));

  // run_speed_test(lower, upper, speed_test_size);
  // run_precision_test(lower, upper, accuracy_test_size, speed_test_size);
  
  plot_error_behavior(lower, upper, accuracy_test_size);

  return 0;
}


// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
// Compilation and run on laptop
//
//
// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 $(pkg-config --cflags --libs flint) -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
