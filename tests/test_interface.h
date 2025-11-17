#ifndef TEST_INTERFACE_H 

#define TEST_INTERFACE_H 
#include <stdint.h> 
#include "flint/arb.h"


void fill_uniform(double lower, double upper, size_t n, double *vec);
void compare_results_sin(double *x, double *y, double *cum_error, double *max_error, size_t n);
void compare_results_tan(double *x, double *y, double *cum_error, double *max_error, size_t n);
double test_sin_time(int taylor_degree, double lower_bound, double upper_bound, int test_size);
double test_tan_accuracy(int taylor_degree, double lower_bound, double upper_bound, int test_size);
double test_tan_time(int taylor_degree, double lower_bound, double upper_bound, int test_size);
double test_sin_accuracy(int taylor_degree, double lower_bound, double upper_bound, int test_size);

double test_glibc_sin_time(int taylor_degree, double lower_bound, double upper_bound, int test_size);


#endif

