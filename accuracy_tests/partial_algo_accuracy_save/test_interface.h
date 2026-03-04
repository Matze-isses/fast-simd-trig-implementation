#ifndef TEST_INTERFACE_H
#define TEST_INTERFACE_H

#include <stdint.h>
#include <stddef.h>   /* size_t */
#include "flint/arb.h"

void compare_results_tan(double *x, double *y,
                         double *cum_error, double *max_error, double *value_max_error,
                         double *cum_ulp_error, double *max_ulp_error, double *value_max_ulp_error,
                         size_t n);

void compare_results_tan_err(double *x, double *y, double *err, size_t n);
void compare_results_tan_ulp_err_signed(double *x, double *y, int64_t *ulp_err, size_t n);


#endif /* TEST_INTERFACE_H */
