/* test_interface.c */
#include "clock_utils.h"
#include <omp.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

#include "test_interface.h"
#include "../trig_simd.h"

#include "../cmeasure/cmeasure.h"
#include "../cmeasure/cbind_to_hw_thread.h"
#include "../cmeasure/CrystalClockInC.h"

typedef enum {
  F_TAN = 0,
  F_SIN = 1
} func_kind_t;

typedef void (*simd_fn_t)(double*, double*, size_t);
typedef double (*libm_fn_t)(double);

static inline const char* func_name(func_kind_t k) { return (k == F_SIN) ? "sin" : "tan"; }
static inline simd_fn_t   pick_simd(func_kind_t k) { return (k == F_SIN) ? sin_simd : tan_simd; }
static inline libm_fn_t   pick_libm(func_kind_t k) { return (k == F_SIN) ? sin : tan; }

static inline void compare_results_generic(func_kind_t fk,
                                           double *x, double *y,
                                           double *cum_error, double *max_error, double *value_max_error,
                                           double *cum_ulp_error, double *max_ulp_error, double *value_max_ulp_error,
                                           size_t n)
{
  if (fk == F_SIN) {
    compare_results_sin(x, y, cum_error, max_error, value_max_error,
                        cum_ulp_error, max_ulp_error, value_max_ulp_error, n);
  } else {
    compare_results_tan(x, y, cum_error, max_error, value_max_error,
                        cum_ulp_error, max_ulp_error, value_max_ulp_error, n);
  }
}

static inline void compare_results_err_generic(func_kind_t fk,
                                               double *x, double *y, double *err, size_t n)
{
  if (fk == F_SIN) compare_results_sin_err(x, y, err, n);
  else            compare_results_tan_err(x, y, err, n);
}

static inline void compare_results_ulp_signed_generic(func_kind_t fk,
                                                      double *x, double *y, int64_t *ulp_err, size_t n)
{
  if (fk == F_SIN) compare_results_sin_ulp_err_signed(x, y, ulp_err, n);
  else            compare_results_tan_ulp_err_signed(x, y, ulp_err, n);
}

static void speed_test(func_kind_t fk, double *test_values, double *own_results, size_t test_size)
{
  simd_fn_t own_fn = pick_simd(fk);
  libm_fn_t libm   = pick_libm(fk);

  struct CrystalClock clk;
  clk._freq = frequency();

  double own_execution_ms_warmup = 0, own_execution_ms = 0;
  double glibc_execution_ms_warmup = 0, glibc_execution_ms = 0;
  double own_time_old_clock = -1.0, glibc_time_old_clock = -1.0;

  /* ---- WARMUP (own) ---- */
  START_TCLOCK;
  own_fn(test_values, own_results, test_size);
  own_execution_ms_warmup = GET_TCLOCK;
  printf("WARMUP %s Time OC: %.17g\n", func_name(fk), own_execution_ms_warmup);

  clk._begin = current();
  own_fn(test_values, own_results, test_size);
  clk._end   = current();
  own_execution_ms_warmup = duration_ms1(clk);
  printf("WARMUP %s Time CC: %.17g\n\n", func_name(fk), own_execution_ms_warmup);

  /* ---- Own Results ---- */
  START_TCLOCK;
  own_fn(test_values, own_results, test_size);
  own_time_old_clock = GET_TCLOCK;
  printf("Own   %s Time OC: %.17g\n", func_name(fk), own_time_old_clock);

  clk._begin = current();
  own_fn(test_values, own_results, test_size);
  clk._end   = current();
  own_execution_ms = duration_ms1(clk);
  printf("Own   %s Time CC: %.17g\n\n", func_name(fk), own_execution_ms);

  /* ---- GLIBC Warmup Results ---- */
  double *glibc_results = (double*)malloc(test_size * sizeof(double));
  if (!glibc_results) {
    fprintf(stderr, "malloc failed in speed_test (glibc_results)\n");
    return;
  }

  START_TCLOCK;
  for (size_t i = 0; i < test_size; i++) glibc_results[i] = libm(test_values[i]);
  glibc_time_old_clock = GET_TCLOCK;
  printf("Warmup GLIBC %s Time OC: %.17g\n", func_name(fk), glibc_time_old_clock);

  clk._begin = current();
  for (size_t i = 0; i < test_size; i++) glibc_results[i] = libm(test_values[i]);
  clk._end   = current();
  glibc_execution_ms_warmup = duration_ms1(clk);
  printf("Warmup GLIBC %s Time CC: %.17g\n\n", func_name(fk), glibc_execution_ms_warmup);

  /* ---- GLIBC Results ---- */
  START_TCLOCK;
  for (size_t i = 0; i < test_size; i++) glibc_results[i] = libm(test_values[i]);
  glibc_time_old_clock = GET_TCLOCK;
  printf("GLIBC %s Time OC: %.17g\n", func_name(fk), glibc_time_old_clock);

  clk._begin = current();
  for (size_t i = 0; i < test_size; i++) glibc_results[i] = libm(test_values[i]);
  clk._end   = current();
  glibc_execution_ms = duration_ms1(clk);

  printf("GLIBC %s Time CC: %.17g\n\n", func_name(fk), glibc_execution_ms);

  printf("Speedup: (GLIBC/OWN): %.17g\n\n", glibc_execution_ms / own_execution_ms);

  free(glibc_results);
}

static void quadrant_error_test(func_kind_t fk, size_t n)
{
  simd_fn_t own_fn = pick_simd(fk);
  libm_fn_t libm   = pick_libm(fk);

  const char *interval_names_sin[4] = {
    "Q1 (0, pi/2)",
    "Q2 (pi/2, pi)",
    "Q3 (pi, 3pi/2)",
    "Q4 (3pi/2, 2pi)"
  };

  const char *interval_names_tan[4] = {
    "[0, pi/8]",
    "[pi/8, pi/4]",
    "[pi/4, 3pi/8]",
    "[3pi/8, pi/2]"
  };

  double bounds_sin[5] = {
    nextafter(0.0, 1.0),
    nextafter(M_PI / 2.0, 0.0),
    nextafter(M_PI, 0.0),
    nextafter(3.0 * M_PI / 2.0, M_PI),
    nextafter(2.0 * M_PI, 0.0)
  };

  /* For tan: stay in (0, pi/2) but split into 4 sub-intervals */
  double bounds_tan[5] = {
    1e-11,                /* avoid exactly 0 if you want (optional) */
    M_PI / 8.0,
    M_PI / 4.0,
    3.0 * M_PI / 8.0,
    nextafter(M_PI / 2.0, 0.0) /* avoid pole */
  };

  const char **interval_names = (fk == F_TAN) ? interval_names_tan : interval_names_sin;
  double *bounds              = (fk == F_TAN) ? bounds_tan         : bounds_sin;

  printf("==============================================================================================================\n");
  printf("Quadrant Error Analysis for %s, n = %d\n", func_name(fk), (int)n);
  printf("==============================================================================================================\n\n");

  printf("+----------------+----------------+---------------------------+---------------------------+---------------------------+---------------------------+\n");
  printf("| Interval       | Implementation | Max Error                 | Avg Abs Error             | Max ULP Error             | Avg ULP Error             |\n");
  printf("+----------------+----------------+---------------------------+---------------------------+---------------------------+---------------------------+\n");

  for (int q = 0; q < 4; ++q) {
    double lower = bounds[q];
    double upper = bounds[q + 1];

    double *test_values   = (double*)malloc(n * sizeof(double));
    double *own_results   = (double*)malloc(n * sizeof(double));
    double *glibc_results = (double*)malloc(n * sizeof(double));

    if (!test_values || !own_results || !glibc_results) {
      perror("malloc");
      free(test_values);
      free(own_results);
      free(glibc_results);
      return;
    }

    fill_uniform(lower, upper, n, test_values);

    own_fn(test_values, own_results, n);
    for (size_t i = 0; i < n; ++i) glibc_results[i] = libm(test_values[i]);

    /* error metrics (own) */
    double cum_error_own=0.0, max_error_own=0.0, value_max_error_own=0.0;
    double cum_ulp_error_own=0.0, max_ulp_error_own=0.0, value_max_ulp_error_own=0.0;

    compare_results_generic(fk, test_values, own_results,
                            &cum_error_own, &max_error_own, &value_max_error_own,
                            &cum_ulp_error_own, &max_ulp_error_own, &value_max_ulp_error_own,
                            n);

    double avg_error_own     = cum_error_own / (double)n;
    double avg_ulp_error_own = cum_ulp_error_own / (double)n;

    /* error metrics (glibc) */
    double cum_error_glibc=0.0, max_error_glibc=0.0, value_max_error_glibc=0.0;
    double cum_ulp_error_glibc=0.0, max_ulp_error_glibc=0.0, value_max_ulp_error_glibc=0.0;

    compare_results_generic(fk, test_values, glibc_results,
                            &cum_error_glibc, &max_error_glibc, &value_max_error_glibc,
                            &cum_ulp_error_glibc, &max_ulp_error_glibc, &value_max_ulp_error_glibc,
                            n);

    double avg_error_glibc     = cum_error_glibc / (double)n;
    double avg_ulp_error_glibc = cum_ulp_error_glibc / (double)n;

    printf("| %-14s | %-14s | %25.17e | %25.17e | %25.17e | %25.17e |\n",
           interval_names[q], "own",
           max_error_own, avg_error_own, max_ulp_error_own, avg_ulp_error_own);

    printf("| %-14s | %-14s | %25.17e | %25.17e | %25.17e | %25.17e |\n",
           "", "glibc",
           max_error_glibc, avg_error_glibc, max_ulp_error_glibc, avg_ulp_error_glibc);

    printf("+----------------+----------------+---------------------------+---------------------------+---------------------------+---------------------------+\n");

    free(test_values);
    free(own_results);
    free(glibc_results);
  }

  printf("\n");
}

static void run_accuracy_test(func_kind_t fk, size_t accuracy_test_size)
{
  if (accuracy_test_size > 0) quadrant_error_test(fk, accuracy_test_size);
}

static void run_speed_test(func_kind_t fk, double lower, double upper, int speed_test_size)
{
  if (speed_test_size == 0) return;

  double *speed_test_values = (double*)malloc((size_t)speed_test_size * sizeof(double));
  double *own_speed_results = (double*)malloc((size_t)speed_test_size * sizeof(double));

  if (!speed_test_values || !own_speed_results) {
    fprintf(stderr, "Allocation failed in run_speed_test.\n");
    free(speed_test_values);
    free(own_speed_results);
    return;
  }

  fill_uniform(lower, upper, (size_t)speed_test_size, speed_test_values);
  speed_test(fk, speed_test_values, own_speed_results, (size_t)speed_test_size);

  free(speed_test_values);
  free(own_speed_results);
}

static void run_precision_test(func_kind_t fk, double lower, double upper, size_t accuracy_test_size)
{
  simd_fn_t own_fn = pick_simd(fk);
  libm_fn_t libm   = pick_libm(fk);

  double *test_values     = (double*)malloc(accuracy_test_size * sizeof(double));
  double *correct_results = (double*)malloc(accuracy_test_size * sizeof(double)); /* kept for compatibility */
  double *own_results     = (double*)malloc(accuracy_test_size * sizeof(double));
  double *glibc_results   = (double*)malloc(accuracy_test_size * sizeof(double));

  if (!test_values || !correct_results || !own_results || !glibc_results) {
    fprintf(stderr, "Allocation failed in run_precision_test.\n");
    free(test_values); free(correct_results); free(own_results); free(glibc_results);
    return;
  }

  fill_uniform(lower, upper, accuracy_test_size, test_values);

  own_fn(test_values, own_results, accuracy_test_size);
  for (size_t i = 0; i < accuracy_test_size; i++) glibc_results[i] = libm(test_values[i]);

  printf("\nThe results are obtained! Starting error calculation for %s.\n", func_name(fk));

  double cum_error_own=0.0, max_error_own=0.0, value_max_error_own=0.0;
  double cum_ulp_own=0.0, max_ulp_own=0.0, value_max_ulp_own=0.0;

  double cum_error_glibc=0.0, max_error_glibc=0.0, value_max_error_glibc=0.0;
  double cum_ulp_glibc=0.0, max_ulp_glibc=0.0, value_max_ulp_glibc=0.0;

  compare_results_generic(fk, test_values, own_results,
                          &cum_error_own, &max_error_own, &value_max_error_own,
                          &cum_ulp_own, &max_ulp_own, &value_max_ulp_own,
                          accuracy_test_size);

  compare_results_generic(fk, test_values, glibc_results,
                          &cum_error_glibc, &max_error_glibc, &value_max_error_glibc,
                          &cum_ulp_glibc, &max_ulp_glibc, &value_max_ulp_glibc,
                          accuracy_test_size);

  printf("Max Error Own:   %.17g;   At Value %.17g\n", max_error_own, value_max_error_own);
  printf("Max Error glibc: %.17g;   At Value %.17g\n\n", max_error_glibc, value_max_error_glibc);

  printf("Max ULP Error Own:   %.17g;   At Value %.17g\n", max_ulp_own, value_max_ulp_own);
  printf("Max ULP Error glibc: %.17g;   At Value %.17g\n\n", max_ulp_glibc, value_max_ulp_glibc);

  printf("Accumulated Absolut Error Own   Results: %.17g\n", cum_error_own);
  printf("Accumulated Absolut Error glibc Results: %.17g\n", cum_error_glibc);

  printf("\nAbsolut Error Own   Results: %.17g\n", cum_error_own / (double)accuracy_test_size);
  printf("Absolut Error glibc Results: %.17g\n", cum_error_glibc / (double)accuracy_test_size);

  free(test_values);
  free(correct_results);
  free(own_results);
  free(glibc_results);
}

static void plot_error_behavior(func_kind_t fk, double lower, double upper, size_t accuracy_test_size, int dtype)
{
  if (accuracy_test_size == 0) return;

  simd_fn_t own_fn = pick_simd(fk);

  double *test_values = (double*)malloc(accuracy_test_size * sizeof(double));
  double *own_results = (double*)malloc(accuracy_test_size * sizeof(double));
  double *err         = (double*)malloc(accuracy_test_size * sizeof(double));

  if (!test_values || !own_results || !err) {
    fprintf(stderr, "Allocation failed in plot_error_behavior.\n");
    free(test_values); free(own_results); free(err);
    return;
  }

  if (dtype == 0) {
    if (accuracy_test_size == 1) {
      test_values[0] = lower;
    } else {
      double step = (upper - lower) / (double)(accuracy_test_size - 1);
      for (size_t i = 0; i < accuracy_test_size; i++) test_values[i] = lower + step * (double)i;
      test_values[0] = lower;
      test_values[accuracy_test_size - 1] = upper;
    }
  } else if (dtype == 1) {
    fill_uniform(lower, upper, accuracy_test_size, test_values);
  } else if (dtype == 2) {
    fill_dense_pi_over_2(lower, upper, accuracy_test_size, test_values, 0.01);
  }

  own_fn(test_values, own_results, accuracy_test_size);
  compare_results_err_generic(fk, test_values, own_results, err, accuracy_test_size);

  char data_path[64];
  snprintf(data_path, sizeof(data_path), "%s_error_behavior.tsv", func_name(fk));

  FILE *f = fopen(data_path, "w");
  if (!f) {
    fprintf(stderr, "Failed to open %s for writing.\n", data_path);
    free(test_values); free(own_results); free(err);
    return;
  }

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

static void plot_data_ulp(func_kind_t fk, double lower, double upper, size_t accuracy_test_size, int dtype)
{
  if (accuracy_test_size == 0) return;

  simd_fn_t own_fn = pick_simd(fk);

  double  *test_values = (double*)malloc(accuracy_test_size * sizeof(double));
  double  *own_results = (double*)malloc(accuracy_test_size * sizeof(double));
  int64_t *ulp_err     = (int64_t*)malloc(accuracy_test_size * sizeof(int64_t));

  if (!test_values || !own_results || !ulp_err) {
    fprintf(stderr, "Allocation failed in plot_data_ulp.\n");
    free(test_values); free(own_results); free(ulp_err);
    return;
  }

  if (dtype == 0) {
    if (accuracy_test_size == 1) {
      test_values[0] = lower;
    } else {
      double step = (upper - lower) / (double)(accuracy_test_size - 1);
      for (size_t i = 0; i < accuracy_test_size; i++) test_values[i] = lower + step * (double)i;
      test_values[0] = lower;
      test_values[accuracy_test_size - 1] = upper;
    }
  } else if (dtype == 1) {
    fill_uniform(lower, upper, accuracy_test_size, test_values);
  } else if (dtype == 2) {
    fill_dense_pi_over_2(lower, upper, accuracy_test_size, test_values, 0.01);
  }

  own_fn(test_values, own_results, accuracy_test_size);
  compare_results_ulp_signed_generic(fk, test_values, own_results, ulp_err, accuracy_test_size);

  char data_path[64];
  snprintf(data_path, sizeof(data_path), "%s_ulp_error_behavior.tsv", func_name(fk));

  FILE *f = fopen(data_path, "w");
  if (!f) {
    fprintf(stderr, "Failed to open %s for writing.\n", data_path);
    free(test_values); free(own_results); free(ulp_err);
    return;
  }

  fprintf(f, "x\tulp_err\n");
  for (size_t i = 0; i < accuracy_test_size; i++) {
    fprintf(f, "%.17g\t%" PRId64 "\n", test_values[i], ulp_err[i]);
  }
  fclose(f);

  printf("Wrote data: %s\n", data_path);

  free(test_values);
  free(own_results);
  free(ulp_err);
}

static void interval_accuracy_table(func_kind_t fk, double lower, double upper, size_t n)
{
  simd_fn_t own_fn = pick_simd(fk);
  libm_fn_t libm   = pick_libm(fk);

  double *x            = (double*)malloc(n * sizeof(double));
  double *own_results  = (double*)malloc(n * sizeof(double));
  double *glibc_results= (double*)malloc(n * sizeof(double));

  if (!x || !own_results || !glibc_results) {
    fprintf(stderr, "Allocation failed in interval_accuracy_table.\n");
    free(x); free(own_results); free(glibc_results);
    return;
  }

  fill_uniform(lower, upper, n, x);

  own_fn(x, own_results, n);
  for (size_t i = 0; i < n; ++i) glibc_results[i] = libm(x[i]);

  double cum_error_own=0.0, max_error_own=0.0, value_max_error_own=0.0;
  double cum_ulp_own=0.0,  max_ulp_own=0.0,  value_max_ulp_own=0.0;

  double cum_error_glibc=0.0, max_error_glibc=0.0, value_max_error_glibc=0.0;
  double cum_ulp_glibc=0.0,  max_ulp_glibc=0.0,  value_max_ulp_glibc=0.0;

  compare_results_generic(fk, x, own_results,
                          &cum_error_own, &max_error_own, &value_max_error_own,
                          &cum_ulp_own, &max_ulp_own, &value_max_ulp_own,
                          n);

  compare_results_generic(fk, x, glibc_results,
                          &cum_error_glibc, &max_error_glibc, &value_max_error_glibc,
                          &cum_ulp_glibc, &max_ulp_glibc, &value_max_ulp_glibc,
                          n);

  double avg_error_own     = cum_error_own / (double)n;
  double avg_ulp_error_own = cum_ulp_own   / (double)n;

  double avg_error_glibc     = cum_error_glibc / (double)n;
  double avg_ulp_error_glibc = cum_ulp_glibc   / (double)n;

  printf("==============================================================================================================\n");
  printf("Interval Error Analysis for %s on [%.17g, %.17g], n = %d\n", func_name(fk), lower, upper, (int)n);
  printf("==============================================================================================================\n\n");

  printf("+----------------------+----------------+---------------------------+---------------------------+---------------------------+---------------------------+\n");
  printf("| Interval             | Implementation | Max Error                 | Avg Abs Error             | Max ULP Error             | Avg ULP Error             |\n");
  printf("+----------------------+----------------+---------------------------+---------------------------+---------------------------+---------------------------+\n");

  char interval_label[64];
  snprintf(interval_label, sizeof(interval_label), "[%.6g, %.6g]", lower, upper);

  printf("| %-20s | %-14s | %25.17e | %25.17e | %25.17e | %25.17e |\n",
         interval_label, "own",
         max_error_own, avg_error_own, max_ulp_own, avg_ulp_error_own);

  printf("| %-20s | %-14s | %25.17e | %25.17e | %25.17e | %25.17e |\n",
         "", "glibc",
         max_error_glibc, avg_error_glibc, max_ulp_glibc, avg_ulp_error_glibc);

  printf("+----------------------+----------------+---------------------------+---------------------------+---------------------------+---------------------------+\n\n");

  /* optional: show where worst cases happen */
  printf("Own   worst abs err at x=%.17g;  worst ulp err at x=%.17g\n", value_max_error_own, value_max_ulp_own);
  printf("glibc worst abs err at x=%.17g;  worst ulp err at x=%.17g\n\n", value_max_error_glibc, value_max_ulp_glibc);

  free(x);
  free(own_results);
  free(glibc_results);
}

int main(int argc, char *argv[]) {
  cbind_to_hw_thread(2, 1);

  if (argc < 5) {
    fprintf(stderr, "Usage: %s n lower upper tan|sin [accuracy_test_size]\n", argv[0]);
    return 1;
  }

  int speed_test_size = atoi(argv[1]);
  if (speed_test_size % 4 != 0) speed_test_size += (4 - (speed_test_size % 4));

  double lower = atof(argv[2]);
  double upper = atof(argv[3]);

  func_kind_t fk = F_TAN;
  if (strcmp(argv[4], "sin") == 0) fk = F_SIN;
  else if (strcmp(argv[4], "tan") == 0) fk = F_TAN;
  else {
    fprintf(stderr, "4th arg must be 'tan' or 'sin'\n");
    return 1;
  }

  size_t accuracy_test_size = (argc > 5) ? (size_t)atoi(argv[5]) : (size_t)speed_test_size;

  printf("\n----------------------------------------------------------------------------------------------------\n");
  printf(" Test Setup \n");
  printf("----------------------------------------------------------------------------------------------------\n");
  printf("Function:                             %s\n", func_name(fk));
  printf("Input Interval:                       [%f, %f]\n", lower, upper);
  printf("Number of inputs (speed test):        n=%d\n", speed_test_size);
  printf("Number of inputs (error calculation): n=%d\n", (int)accuracy_test_size);
  printf("----------------------------------------------------------------------------------------------------\n\n");

  srand((unsigned)time(NULL));

  /* Enable/disable as you like */
  /* run_accuracy_test(fk, accuracy_test_size); */

  run_speed_test(fk, lower, upper, speed_test_size);
  interval_accuracy_table(fk, lower, upper, accuracy_test_size);
  // quadrant_error_test(fk, accuracy_test_size);

  /* run_precision_test(fk, lower, upper, accuracy_test_size); */

  /* dtype: 0 linspace, 1 uniform random, 2 dense near pi/2 (mostly for tan) */
  // plot_error_behavior(fk, lower, upper, accuracy_test_size, 1);
  // plot_data_ulp(fk, lower, upper, accuracy_test_size, 1);

  return 0;
}
