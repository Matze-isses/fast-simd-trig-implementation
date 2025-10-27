#include "clock_utils.h"
#include <omp.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "test_interface.h"
#include "../range_reduction.h"

int main(int argc, char *argv[]) {

  if (argc < 3) {
    fprintf(stderr, "Usage: %s n lower upper\n", argv[0]);
    return 1;
  }

  int n = atoi(argv[1]);
  double lower = atof(argv[2]);
  double upper = atof(argv[3]);

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


  printf("\n -------- Own Script Execution ---------- \n\n");
  START_CLOCK;
  sin_simd(test_values, own_results, n, 0.1);
  END_CLOCK("\n\n ------- End Own Script Execution ------- \n\nTime needed by own implementiation ");

  START_CLOCK;
  for (int i = 0; i < n; i++) { glibc_results[i] = sin(test_values[i]); }
  END_CLOCK("Time needed by glibc               ");

  printf("\n -------- Second Own Script Execution ---------- \n\n");
  START_CLOCK;
  sin_simd(test_values, own_results, n, 0.1);
  END_CLOCK("\n\n ------- Second End Own Script Execution ------- \n\nTime needed by own implementiation (second) ");

  START_CLOCK;
  for (int i = 0; i < n; i++) { glibc_results[i] = sin(test_values[i]); }
  END_CLOCK("Time needed by glibc (second)               ");
  // precision test


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

  double abs_error = compare_results(test_values_partial, own_results_partial, test_size);
  double abs_error_glibc = compare_results(test_values_partial, glibc_results_partial, test_size);

  printf("Accumulated Absolut Error Own   Results: %.17g\n", abs_error);
  printf("Accumulated Absolut Error glibc Results: %.17g\n", abs_error_glibc);

  printf("\nAbsolut Error Own   Results: %.17g\n", abs_error/n);
  printf("Absolut Error glibc Results: %.17g\n", abs_error_glibc/n);

  free(test_values);
  free(correct_results);
  free(glibc_results);
  free(own_results);

  free(correct_results_partial);
  free(glibc_results_partial);
  free(own_results_partial);
  free(test_values_partial);

  return 0;
}

// gcc test_interface.c value_generation.c sin_arb.c -o test -lm -mavx -O2 -lflint -Wall -Wextra && ./test 1000000 0 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
