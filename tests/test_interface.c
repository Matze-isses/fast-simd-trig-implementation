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
  if (argc != 4) {
    fprintf(stderr, "Usage: %s n lower upper\n", argv[0]);
    return 1;
  }

  int n = atoi(argv[1]);
  double lower = atof(argv[2]);
  double upper = atof(argv[3]);

  printf("\nThe test is executed for n=%d test-values within the interval [%f, %f].\n", n, lower, upper);

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
  printf("Test values are generated! Starting calculation.\n");

  // precision test
  #pragma omp parallel for
  for (int i = 0; i < n; i++){
    correct_results[i] = correct_result_sin(test_values[i]);
  }

  printf(" ------- Own Script Execution ------- \n\n");
  START_CLOCK;
  sin_simd(test_values, own_results, n, 0.1);
  END_CLOCK("Time needed by own implementiation ");

  START_CLOCK;
  for (int i = 0; i < n; i++) { glibc_results[i] = sin(test_values[i]); }
  END_CLOCK("Time needed by glibc               ");


  // user information for the current state of the script
  printf("The results are obtained! Starting error calculation.\n");

  double abs_error = 0;
  double abs_error_glibc = 0;
  for (int i = 0; i < n; i++){
    abs_error += fabs(correct_results[i] - own_results[i]);
    abs_error_glibc += fabs(correct_results[i] - glibc_results[i]);
  }

  printf("Accumulated Absolut Error Own   Results: %.17g\n", abs_error);
  printf("Accumulated Absolut Error glibc Results: %.17g\n", abs_error_glibc);

  printf("\nAbsolut Error Own   Results: %.17g\n", abs_error/n);
  printf("Absolut Error glibc Results: %.17g\n", abs_error_glibc/n);
  free(test_values);
  free(correct_results);
  free(own_results);
  free(glibc_results);

  return 0;
}

// gcc test_interface.c value_generation.c sin_arb.c -o test -lm -mavx -O2 -lflint -Wall -Wextra && ./test 1000000 0 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
