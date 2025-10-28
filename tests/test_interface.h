#ifndef TEST_INTERFACE_H 

#define TEST_INTERFACE_H 
#include <stdint.h> 
#include "flint/arb.h"


void fill_uniform(double lower, double upper, size_t n, double *vec);
double compare_results_sin(double *x, double *y, size_t n);
double compare_results_tan(double *x, double *y, size_t n);

#endif

