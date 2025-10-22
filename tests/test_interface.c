#include "clock_utils.h"
#include <omp.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "test_interface.h"

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
    glibc_results[i] = sin(test_values[i]);
  }

  // user information for the current state of the script
  printf("The results are obtained! Starting error calculation.\n");

  double abs_error = 0;
  for (int i = 0; i < n; i++){
    abs_error += abs(correct_results[i] - glibc_results[i]);
  }

  printf("Absolut Error glibc: %.17g\n", abs_error);

  return 0;
}
