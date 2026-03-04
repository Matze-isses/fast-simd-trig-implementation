#include "flint/flint.h"
#include "flint/arb.h"
#include "test_interface.h"

#include <stdint.h>
#include <stddef.h>
#include <string.h>   // memcpy
#include <math.h>     // isnan, isinf
                    
const double M_PI_8 = 0x1.921fb54442d18p-2;

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

static inline int64_t ulp_error_signed_double(double f_hat, double f_ref)
{
    // Policy: if NaN involved, return 0 (or choose a sentinel if you prefer).
    if (isnan(f_hat) || isnan(f_ref)) return 0;

    // If either is inf: return 0 only if equal, else saturate to INT64_{MIN,MAX}.
    if (isinf(f_hat) || isinf(f_ref)) {
        if (f_hat == f_ref) return 0;
        return (f_ref > f_hat) ? INT64_MAX : INT64_MIN;
    }

    uint64_t oh = double_to_ordered_u64(f_hat);
    uint64_t orr = double_to_ordered_u64(f_ref);

    // ord(ref) - ord(hat)
    if (orr >= oh) {
        uint64_t d = orr - oh;
        return (d > (uint64_t)INT64_MAX) ? INT64_MAX : (int64_t)d;
    } else {
        uint64_t d = oh - orr;
        return (d > (uint64_t)INT64_MAX) ? INT64_MIN : -(int64_t)d;
    }
}

const slong PRECISION = 1024;


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

double reference_calculation(double x)
{
    arb_t arb_x, before_floor, floored, reduced_range, pi_2, tmp, q2_tmp, true_result, true_result_tan;
    arb_init(arb_x);
    arb_init(pi_2);
    arb_init(tmp);
    arb_init(before_floor);
    arb_init(floored);

    arb_init(reduced_range);
    arb_init(q2_tmp);
    arb_init(true_result);

    arb_init(true_result_tan);


    arb_set_d(arb_x, x);

    arb_const_pi(tmp, PRECISION);          // tmp = pi
    arb_mul_2exp_si(pi_2, tmp, -1);        // pi_2 = pi / 2
    
    //Flint stuff
    arb_div(before_floor, arb_x, pi_2, PRECISION);
    arb_floor(floored, before_floor, PRECISION);
    arb_sub(reduced_range, arb_x, floored, PRECISION);

    if (x < M_PI_8) {
        arb_mul_2exp_si(true_result, reduced_range, 0);

    } else if (M_PI_8 <= x && x < M_PI_4) {
        arb_mul_2exp_si(true_result, reduced_range, -1);

    } else if (M_PI_4 <= x && x < 3.0 * M_PI_8) {
        arb_sub(q2_tmp, pi_2, arb_x, PRECISION);
        arb_mul_2exp_si(true_result, q2_tmp, -1);

    } else if (3.0 * M_PI_8 <= x && x < M_PI_2) {
        arb_sub(true_result, pi_2, arb_x, PRECISION);
    }

    // double res = arf_get_d(arb_midref(true_result), ARF_RND_NEAR);

    arb_tan(true_result_tan, arb_x, PRECISION);
    double res = arf_get_d(arb_midref(true_result_tan), ARF_RND_NEAR);


    arb_clear(arb_x);
    arb_clear(pi_2);
    arb_clear(tmp);
    arb_clear(before_floor);
    arb_clear(floored);

    arb_clear(reduced_range);
    arb_clear(q2_tmp);
    arb_clear(true_result);

    arb_clear(true_result_tan);

    // If you want to return both, see note below.
    // For now, returning res as requested:
    return res;
}

void compare_results_tan_ulp_err_signed(double *x, double *y, int64_t *ulp_err, size_t n) {

    for (size_t i = 0; i < n; i++) {
        ulp_err[i] = ulp_error_signed_double(y[i], reference_calculation(x[i]));
    }

    flint_cleanup();
}
