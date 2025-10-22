#include "flint/flint.h"
#include "flint/arb.h"
#include "test_interface.h"

const slong PRECISION = 1024;


double correct_result_sin(double x) {
  arb_t arb_x;
  arb_t result;

  arb_init(arb_x);
  arb_init(result);
  
  arb_set_d(arb_x, x);
  arb_sin(result, arb_x, PRECISION); 

  double res = arf_get_d(arb_midref(result), ARF_RND_NEAR);

  arb_clear(arb_x);
  arb_clear(result);

  return res;
}
