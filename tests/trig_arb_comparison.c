#include "flint/flint.h"
#include "flint/arb.h"
#include "test_interface.h"

#include <stdint.h>
#include <stddef.h>
#include <string.h>   // memcpy
#include <math.h>     // isnan, isinf

// Map a double to an ordered 64-bit integer so that integer difference equals ULP distance.
static inline uint64_t double_to_ordered_u64(double a)
{
    uint64_t u;
    memcpy(&u, &a, sizeof(u));

    // Make ordering monotone with respect to numeric ordering:
    // Negative doubles map to [0 .. 0x7FFF..] in reverse, positives map to [0x8000.. .. 0xFFFF..]
    if (u & 0x8000000000000000ULL)
        u = 0x8000000000000000ULL - u;
    else
        u = u + 0x8000000000000000ULL;

    return u;
}

static inline uint64_t ulp_distance_double(double a, double b)
{
    // Convention: if either is NaN, return max. You can choose a different policy if desired.
    if (isnan(a) || isnan(b)) return UINT64_MAX;

    // If both are infinities of same sign => distance 0; otherwise max.
    if (isinf(a) || isinf(b)) {
        return (a == b) ? 0ULL : UINT64_MAX;
    }

    uint64_t ua = double_to_ordered_u64(a);
    uint64_t ub = double_to_ordered_u64(b);
    return (ua > ub) ? (ua - ub) : (ub - ua);
}

const slong PRECISION = 512;


void compare_results_sin(double *x, double *y, double *cum_error, double *max_error, double *value_max_error, size_t n) {
  // Error should be continously added to and is therefore set before
  // i do not want to go into vector operations with flint...
  arb_t error;
  arb_t arb_max_error;
  
  arb_init(error);
  arb_init(arb_max_error);
  arb_set_d(error, 0.0);

  arb_t arb_value_max_error;
  arb_init(arb_value_max_error);
  arb_set_d(arb_value_max_error, -41);


  arb_t arb_x;
  arb_t arb_y;
  arb_t true_result;
  arb_t difference;
  arb_t diff_to_max;

  arb_init(arb_x);
  arb_init(arb_y);
  arb_init(true_result);
  arb_init(difference);
  arb_init(diff_to_max);
  
  for (int i = 0; i < (int)n; i++) {
    // Initialize all Variables
    arb_set_d(arb_x, x[i]);
    arb_set_d(arb_y, y[i]);

    // calculate sin
    arb_sin(true_result, arb_x, PRECISION); 
    
    // get the error of the calculation
    arb_sub(difference, true_result, arb_y, PRECISION);
    arb_abs(difference, difference);

    arb_add(error, error, difference, PRECISION);
    arb_max(arb_max_error, arb_max_error, difference, PRECISION);

    // get the error of the calculation
    arb_sub(difference, true_result, arb_y, PRECISION);
    arb_abs(difference, difference);
    arb_add(error, error, difference, PRECISION);

    arb_sub(diff_to_max, difference, arb_max_error, PRECISION);
    if (arb_is_positive(diff_to_max)) {
      arb_max(arb_max_error, arb_max_error, difference, PRECISION);
      *value_max_error = x[i];
    }

  }
  // cleanup
  arb_clear(arb_x);
  arb_clear(arb_y);
  arb_clear(true_result);
  arb_clear(difference);
  arb_clear(diff_to_max);

  *cum_error = arf_get_d(arb_midref(error), ARF_RND_NEAR);
  *max_error = arf_get_d(arb_midref(arb_max_error), ARF_RND_NEAR);

  arb_clear(arb_max_error);
  arb_clear(error);
  flint_cleanup();
}


void compare_results_tan(double *x, double *y, double *cum_error, double *max_error, double *value_max_error, double *cum_ulp_error, double *max_ulp_error, double *value_max_ulp_error, size_t n) {
    arb_t error;
    arb_t arb_max_error;

    arb_init(error);
    arb_init(arb_max_error);

    arb_set_d(error, 0.0);
    arb_set_d(arb_max_error, 0.0);

    // init outputs
    if (cum_ulp_error) *cum_ulp_error = 0.0;
    if (max_ulp_error) *max_ulp_error = 0.0;
    if (value_max_ulp_error) *value_max_ulp_error = 0.0;

    for (size_t i = 0; i < n; i++) {
        arb_t arb_x, arb_y, true_result, difference, diff_to_max;
        arb_init(arb_x);
        arb_init(arb_y);
        arb_init(true_result);
        arb_init(difference);
        arb_init(diff_to_max);

        arb_set_d(arb_x, x[i]);
        arb_set_d(arb_y, y[i]);

        // Reference tan in Arb
        arb_tan(true_result, arb_x, PRECISION);

        // Absolute error (your existing logic)
        arb_sub(difference, true_result, arb_y, PRECISION);
        arb_abs(difference, difference);
        arb_add(error, error, difference, PRECISION);

        arb_sub(diff_to_max, difference, arb_max_error, PRECISION);
        if (arb_is_positive(diff_to_max)) {
            arb_set(arb_max_error, difference);
            *value_max_error = x[i];
        }

        // --- ULP error ---
        if (cum_ulp_error || max_ulp_error) {
            // Convert Arb reference midpoint to a double (nearest)
            double ref = arf_get_d(arb_midref(true_result), ARF_RND_NEAR);
            double got = y[i];

            uint64_t ulp = ulp_distance_double(got, ref);

            if (cum_ulp_error) {
                // If ulp is UINT64_MAX (NaN/inf mismatch), this will overflow meaningfully; adjust policy if needed.
                *cum_ulp_error += (double)ulp;
            }
            if (max_ulp_error && (double)ulp > *max_ulp_error) {
                *max_ulp_error = (double)ulp;
                if (value_max_ulp_error) *value_max_ulp_error = x[i];
            }
        }

        arb_clear(arb_x);
        arb_clear(arb_y);
        arb_clear(true_result);
        arb_clear(difference);
        arb_clear(diff_to_max);
    }

    *cum_error = arf_get_d(arb_midref(error), ARF_RND_NEAR);
    *max_error = arf_get_d(arb_midref(arb_max_error), ARF_RND_NEAR);

    arb_clear(error);
    arb_clear(arb_max_error);

    flint_cleanup();
}


void compare_results_tan_err(double *x, double *y, double *err, size_t n) {
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

    // get the difference: tan(x) - y
    arb_sub(difference, true_result, arb_y, PRECISION);

    // store as double
    err[i] = arf_get_d(arb_midref(difference), ARF_RND_NEAR);

    // cleanup
    arb_clear(arb_x);
    arb_clear(arb_y);
    arb_clear(true_result);
    arb_clear(difference);
  }

  flint_cleanup();
}
