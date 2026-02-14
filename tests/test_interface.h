#ifndef TEST_INTERFACE_H
#define TEST_INTERFACE_H

#include <stdint.h>
#include <stddef.h>   /* size_t */
#include "flint/arb.h"

/* input generators */
void fill_uniform(double lower, double upper, size_t n, double *vec);
void fill_dense_pi_over_2(double lower, double upper, size_t n, double *vec, double sigma);

/* Arb-based reference comparisons (absolute + ULP where applicable) */
void compare_results_sin(double *x, double *y,
                         double *cum_error, double *max_error, double *value_max_error,
                         double *cum_ulp_error, double *max_ulp_error, double *value_max_ulp_error,
                         size_t n);

void compare_results_sin_err(double *x, double *y, double *err, size_t n);
void compare_results_sin_ulp_err_signed(double *x, double *y, int64_t *ulp_err, size_t n);

void compare_results_tan(double *x, double *y,
                         double *cum_error, double *max_error, double *value_max_error,
                         double *cum_ulp_error, double *max_ulp_error, double *value_max_ulp_error,
                         size_t n);

void compare_results_tan_err(double *x, double *y, double *err, size_t n);
void compare_results_tan_ulp_err_signed(double *x, double *y, int64_t *ulp_err, size_t n);

/* legacy / other tests (kept as-is) */
double test_sin_time(int taylor_degree, double lower_bound, double upper_bound, int test_size);
double test_tan_accuracy(int taylor_degree, double lower_bound, double upper_bound, int test_size);
double test_tan_time(int taylor_degree, double lower_bound, double upper_bound, int test_size);
double test_sin_accuracy(int taylor_degree, double lower_bound, double upper_bound, int test_size);

double test_glibc_sin_time(int taylor_degree, double lower_bound, double upper_bound, int test_size);

#endif /* TEST_INTERFACE_H */
