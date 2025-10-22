#ifndef TEST_INTERFACE_H 

#define TEST_INTERFACE_H 
#include <stdint.h> 
#include "flint/arb.h"


void fill_uniform(double lower, double upper, size_t n, double *vec);
double correct_result_sin(double x);

#endif

