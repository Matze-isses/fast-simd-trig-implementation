#include "clock_utils.h"
#include <omp.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "test_interface.h"
#include "../trig_simd.h"

double test_glibc_sin_time(int taylor_degree, double lower_bound, double upper_bound, int test_size) {
  srand((unsigned)time(NULL));

  double *test_values = malloc(test_size * sizeof(double));
  double *glibc_results = malloc(test_size * sizeof(double));


  if (!test_values) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower_bound, upper_bound, test_size, test_values);

  START_TCLOCK;
  for (int i = 0; i < test_size; i++) { glibc_results[i] = sin(test_values[i]); }
  double total_time = GET_TCLOCK;    

  free(test_values);
  free(glibc_results);

  return total_time;
}


double test_sin_time(int taylor_degree, double lower_bound, double upper_bound, int test_size) {
  srand((unsigned)time(NULL));

  double *test_values = malloc(test_size * sizeof(double));
  double *own_results = malloc(test_size * sizeof(double));


  if (!test_values) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower_bound, upper_bound, test_size, test_values);

  START_TCLOCK;
  sin_simd(test_values, own_results, test_size, taylor_degree);
  double total_time = GET_TCLOCK;    

  free(test_values);
  free(own_results);

  return total_time;
}

double test_tan_time(int taylor_degree, double lower_bound, double upper_bound, int test_size) {
  srand((unsigned)time(NULL));

  double *test_values = malloc(test_size * sizeof(double));
  double *own_results = malloc(test_size * sizeof(double));


  if (!test_values) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower_bound, upper_bound, test_size, test_values);

  START_TCLOCK;
  tan_simd(test_values, own_results, test_size, taylor_degree);
  double total_time = GET_TCLOCK;    

  free(test_values);
  free(own_results);

  return total_time;
}

double test_sin_accuracy(int taylor_degree, double lower_bound, double upper_bound, int test_size) {
  srand((unsigned)time(NULL));

  double *test_values = malloc(test_size * sizeof(double));
  double *own_results = malloc(test_size * sizeof(double));

  fill_uniform(lower_bound, upper_bound, test_size, test_values);

  sin_simd(test_values, own_results, test_size, taylor_degree);
  double abs_error = compare_results_sin(test_values, own_results, test_size);
  double mean_abs_error = abs_error / test_size;

  free(test_values);
  free(own_results);
  return mean_abs_error;
}

double test_tan_accuracy(int taylor_degree, double lower_bound, double upper_bound, int test_size) {
  srand((unsigned)time(NULL));

  double *test_values = malloc(test_size * sizeof(double));
  double *own_results = malloc(test_size * sizeof(double));

  fill_uniform(lower_bound, upper_bound, test_size, test_values);

  tan_simd(test_values, own_results, test_size, taylor_degree);
  double abs_error = compare_results_tan(test_values, own_results, test_size);
  double mean_abs_error = abs_error / test_size;

  free(test_values);
  free(own_results);
  return mean_abs_error;
}


// gcc -O3 -fPIC -mavx -mavx2 -mfma -fopenmp -c tests/test_object.c trig_simd.c tests/value_generation.c tests/trig_arb_comparison.c && gcc -shared -o libtest_object.so test_object.o trig_simd.o value_generation.o trig_arb_comparison.o -lm -lflint -fopenmp
