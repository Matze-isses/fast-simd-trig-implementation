#include "flint/flint.h"
#include "flint/arb.h"
#include "test_interface.h"

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <limits.h>

const slong PRECISION = 512;

/* Map a double to an ordered 64-bit integer so that integer difference equals ULP distance. */
static inline uint64_t double_to_ordered_u64(double a)
{
    uint64_t u;
    memcpy(&u, &a, sizeof(u));

    /* Make ordering monotone with respect to numeric ordering:
       Negative doubles map to [0 .. 0x7FFF..] in reverse, positives map to [0x8000.. .. 0xFFFF..] */
    if (u & 0x8000000000000000ULL)
        u = 0x8000000000000000ULL - u;
    else
        u = u + 0x8000000000000000ULL;

    return u;
}

static inline uint64_t ulp_distance_double(double a, double b)
{
    if (isnan(a) || isnan(b)) return UINT64_MAX;

    if (isinf(a) || isinf(b)) {
        return (a == b) ? 0ULL : UINT64_MAX;
    }

    uint64_t ua = double_to_ordered_u64(a);
    uint64_t ub = double_to_ordered_u64(b);
    return (ua > ub) ? (ua - ub) : (ub - ua);
}

static inline int64_t ulp_error_signed_double(double f_hat, double f_ref)
{
    if (isnan(f_hat) || isnan(f_ref)) return 0;

    if (isinf(f_hat) || isinf(f_ref)) {
        if (f_hat == f_ref) return 0;
        return (f_ref > f_hat) ? INT64_MAX : INT64_MIN;
    }

    uint64_t oh  = double_to_ordered_u64(f_hat);
    uint64_t orr = double_to_ordered_u64(f_ref);

    /* ord(ref) - ord(hat) */
    if (orr >= oh) {
        uint64_t d = orr - oh;
        return (d > (uint64_t)INT64_MAX) ? INT64_MAX : (int64_t)d;
    } else {
        uint64_t d = oh - orr;
        return (d > (uint64_t)INT64_MAX) ? INT64_MIN : -(int64_t)d;
    }
}

/* ========================= TAN ========================= */

void compare_results_tan(double *x, double *y,
                         double *cum_error, double *max_error, double *value_max_error,
                         double *cum_ulp_error, double *max_ulp_error, double *value_max_ulp_error,
                         size_t n)
{
    arb_t error, arb_max_error;
    arb_init(error);
    arb_init(arb_max_error);

    arb_set_d(error, 0.0);
    arb_set_d(arb_max_error, 0.0);

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

        arb_tan(true_result, arb_x, PRECISION);

        arb_sub(difference, true_result, arb_y, PRECISION);
        arb_abs(difference, difference);
        arb_add(error, error, difference, PRECISION);

        arb_sub(diff_to_max, difference, arb_max_error, PRECISION);
        if (arb_is_positive(diff_to_max)) {
            arb_set(arb_max_error, difference);
            if (value_max_error) *value_max_error = x[i];
        }

        if (cum_ulp_error || max_ulp_error) {
            double ref = arf_get_d(arb_midref(true_result), ARF_RND_NEAR);
            double got = y[i];
            uint64_t ulp = ulp_distance_double(got, ref);

            if (cum_ulp_error) *cum_ulp_error += (double)ulp;
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

    if (cum_error) *cum_error = arf_get_d(arb_midref(error), ARF_RND_NEAR);
    if (max_error) *max_error = arf_get_d(arb_midref(arb_max_error), ARF_RND_NEAR);

    arb_clear(error);
    arb_clear(arb_max_error);

    flint_cleanup();
}

void compare_results_tan_err(double *x, double *y, double *err, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        arb_t arb_x, arb_y, true_result, difference;
        arb_init(arb_x);
        arb_init(arb_y);
        arb_init(true_result);
        arb_init(difference);

        arb_set_d(arb_x, x[i]);
        arb_set_d(arb_y, y[i]);

        arb_tan(true_result, arb_x, PRECISION);
        arb_sub(difference, true_result, arb_y, PRECISION);

        err[i] = arf_get_d(arb_midref(difference), ARF_RND_NEAR);

        arb_clear(arb_x);
        arb_clear(arb_y);
        arb_clear(true_result);
        arb_clear(difference);
    }

    flint_cleanup();
}

void compare_results_tan_ulp_err_signed(double *x, double *y, int64_t *ulp_err, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        arb_t arb_x, true_result;
        arb_init(arb_x);
        arb_init(true_result);

        arb_set_d(arb_x, x[i]);
        arb_tan(true_result, arb_x, PRECISION);

        double ref = arf_get_d(arb_midref(true_result), ARF_RND_NEAR);
        ulp_err[i] = ulp_error_signed_double(y[i], ref);

        arb_clear(arb_x);
        arb_clear(true_result);
    }

    flint_cleanup();
}

/* ========================= SIN ========================= */

void compare_results_sin(double *x, double *y,
                         double *cum_error, double *max_error, double *value_max_error,
                         double *cum_ulp_error, double *max_ulp_error, double *value_max_ulp_error,
                         size_t n)
{
    arb_t error, arb_max_error;
    arb_init(error);
    arb_init(arb_max_error);

    arb_set_d(error, 0.0);
    arb_set_d(arb_max_error, 0.0);

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

        arb_sin(true_result, arb_x, PRECISION);

        arb_sub(difference, true_result, arb_y, PRECISION);
        arb_abs(difference, difference);
        arb_add(error, error, difference, PRECISION);

        arb_sub(diff_to_max, difference, arb_max_error, PRECISION);
        if (arb_is_positive(diff_to_max)) {
            arb_set(arb_max_error, difference);
            if (value_max_error) *value_max_error = x[i];
        }

        if (cum_ulp_error || max_ulp_error) {
            double ref = arf_get_d(arb_midref(true_result), ARF_RND_NEAR);
            double got = y[i];
            uint64_t ulp = ulp_distance_double(got, ref);

            if (cum_ulp_error) *cum_ulp_error += (double)ulp;
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

    if (cum_error) *cum_error = arf_get_d(arb_midref(error), ARF_RND_NEAR);
    if (max_error) *max_error = arf_get_d(arb_midref(arb_max_error), ARF_RND_NEAR);

    arb_clear(error);
    arb_clear(arb_max_error);

    flint_cleanup();
}

void compare_results_sin_err(double *x, double *y, double *err, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        arb_t arb_x, arb_y, true_result, difference;
        arb_init(arb_x);
        arb_init(arb_y);
        arb_init(true_result);
        arb_init(difference);

        arb_set_d(arb_x, x[i]);
        arb_set_d(arb_y, y[i]);

        arb_sin(true_result, arb_x, PRECISION);
        arb_sub(difference, true_result, arb_y, PRECISION);

        err[i] = arf_get_d(arb_midref(difference), ARF_RND_NEAR);

        arb_clear(arb_x);
        arb_clear(arb_y);
        arb_clear(true_result);
        arb_clear(difference);
    }

    flint_cleanup();
}

void compare_results_sin_ulp_err_signed(double *x, double *y, int64_t *ulp_err, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        arb_t arb_x, true_result;
        arb_init(arb_x);
        arb_init(true_result);

        arb_set_d(arb_x, x[i]);
        arb_sin(true_result, arb_x, PRECISION);

        double ref = arf_get_d(arb_midref(true_result), ARF_RND_NEAR);
        ulp_err[i] = ulp_error_signed_double(y[i], ref);

        arb_clear(arb_x);
        arb_clear(true_result);
    }

    flint_cleanup();
}
