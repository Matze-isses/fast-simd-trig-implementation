#include "flint/flint.h"
#include "flint/arb.h"
#include "test_interface.h"

const slong PRECISION = 512;


double compare_results(double *x, double *y, size_t n) {
  // Error should be continously added to and is therefore set before
  // i do not want to go into vector operations with flint...
  arb_t error;
  arb_init(error);
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

    // cleanup
    arb_clear(arb_x);
    arb_clear(arb_y);
    arb_clear(true_result);
    arb_clear(difference);
  }

  double res = arf_get_d(arb_midref(error), ARF_RND_NEAR);
  arb_clear(error);
  flint_cleanup();

  return res;
}
