#include "flint/flint.h"
#include "flint/arb.h"
#include "test_interface.h"

const slong PRECISION = 512;


void compare_results_sin(double *x, double *y, double *cum_error, double *max_error, size_t n) {
  // Error should be continously added to and is therefore set before
  // i do not want to go into vector operations with flint...
  arb_t error;
  arb_t arb_max_error;
  
  arb_init(error);
  arb_init(arb_max_error);
  arb_set_d(error, 0.0);
  
  for (int i = 0; i < (int)n; i++) {
    // Initialize all Variables
    arb_t arb_x;
    arb_t arb_y;
    arb_t true_result;
    arb_t difference;

    arb_init(arb_x);
    arb_init(arb_y);
    arb_init(true_result);
    arb_init(difference);
    
    arb_set_d(arb_x, x[i]);
    arb_set_d(arb_y, y[i]);

    // calculate sin
    arb_sin(true_result, arb_x, PRECISION); 
    
    // get the error of the calculation
    arb_sub(difference, true_result, arb_y, PRECISION);
    arb_abs(difference, difference);

    arb_add(error, error, difference, PRECISION);
    arb_max(arb_max_error, arb_max_error, difference, PRECISION);

    // cleanup
    arb_clear(arb_x);
    arb_clear(arb_y);
    arb_clear(true_result);
    arb_clear(difference);
  }

  *cum_error = arf_get_d(arb_midref(error), ARF_RND_NEAR);
  *max_error = arf_get_d(arb_midref(arb_max_error), ARF_RND_NEAR);

  arb_clear(arb_max_error);
  arb_clear(error);
  flint_cleanup();
}


void compare_results_tan(double *x, double *y, double *cum_error, double *max_error, size_t n) {
  // Error should be continously added to and is therefore set before
  // i do not want to go into vector operations with flint...
  arb_t error;
  arb_t arb_max_error;
  
  arb_init(error);
  arb_init(arb_max_error);

  arb_set_d(error, 0.0);

  arb_t neg_one;
  arb_init(neg_one);
  arb_set_d(neg_one, -1.0);
  
  for (int i = 0; i < (int)n; i++) {
    // Initialize all Variables
    arb_t arb_x;
    arb_t arb_y;
    arb_t true_result;
    arb_t difference;

    arb_init(arb_x);
    arb_init(arb_y);
    arb_init(true_result);
    arb_init(difference);
    
    arb_set_d(arb_x, x[i]);
    arb_set_d(arb_y, y[i]);

    // calculate tan
    arb_tan(true_result, arb_x, PRECISION); 
    
    // get the error of the calculation
    arb_sub(difference, true_result, arb_y, PRECISION);

    arb_abs(difference, difference);

    arb_add(error, error, difference, PRECISION);
    arb_max(arb_max_error, arb_max_error, difference, PRECISION);

    // cleanup
    arb_clear(arb_x);
    arb_clear(arb_y);
    arb_clear(true_result);
    arb_clear(difference);
  }

  *cum_error = arf_get_d(arb_midref(error), ARF_RND_NEAR);
  *max_error = arf_get_d(arb_midref(arb_max_error), ARF_RND_NEAR);

  arb_clear(error);
  arb_clear(arb_max_error);

  flint_cleanup();
}
