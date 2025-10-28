#include "clock_utils.h"
#include <omp.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "test_interface.h"
#include "../trig_simd.h"

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
