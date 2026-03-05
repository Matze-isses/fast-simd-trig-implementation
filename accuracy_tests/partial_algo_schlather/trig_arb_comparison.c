#include "test_interface.h"

#include <stdint.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>   // memcpy
#include <math.h>
#include <stdlib.h>

const double M_PI_8 = 0x1.921fb54442d18p-2;  // pi/8 in double

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

/* No longer used (Arb precision); kept only to minimize diffs */
typedef long slong;
const slong PRECISION = 128;

/* long double pi */
static inline long double pi_ld(void)
{
    return acosl(-1.0L);
}

double reference_calculation(double x)
{
    // Mimic your Arb code computation *as written*, but in long double.
    const long double pi   = pi_ld();
    const long double pi_2 = ldexpl(pi, -1);  // pi/2

    long double arb_x = (long double)x;

    long double true_result = 0.0L;

    if (x < M_PI_8) {
        true_result = arb_x;

    } else if (M_PI_8 <= x && x < M_PI_4) {
        true_result = ldexpl(arb_x, -1);

    } else if (M_PI_4 <= x && x < 3.0 * M_PI_8) {
        long double q2_tmp = pi_2 - arb_x;
        true_result = ldexpl(q2_tmp, -1);

    } else if (3.0 * M_PI_8 <= x && x < M_PI_2) {
        true_result = pi_2 - arb_x;
    }

    long double true_result_tan = tanl(true_result);

    return (double)true_result_tan;
}




double reference_calculation2(double x)
{
    // Mimic your Arb code computation *as written*, but in long double.
    const long double pi   = pi_ld();
    const long double pi_2 = ldexpl(pi, -1);  // pi/2

    long double arb_x = (long double)x;
#define c 1e19
    
    arb_x = truncl(arb_x * c) / c;

    long double true_result = 0.0L;
    

    if (x < M_PI_8) {
        true_result = arb_x;

    } else if (M_PI_8 <= x && x < M_PI_4) {
        true_result = ldexpl(arb_x, -1);

    } else if (M_PI_4 <= x && x < 3.0 * M_PI_8) {
        long double q2_tmp = pi_2 - arb_x;
        true_result = ldexpl(q2_tmp, -1);

    } else if (3.0 * M_PI_8 <= x && x < M_PI_2) {
        true_result = pi_2 - arb_x;
    }

    long double true_result_tan = tanl(true_result);

    return (double)true_result_tan;
}



void compare_results_tan_ulp_err_signed(double *x, double *y, int64_t *ulp_err, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        ulp_err[i] = ulp_error_signed_double(y[i], reference_calculation(x[i]));
        if (abs(ulp_err[i]) >= 2) printf("%ld ", i);
    }
    printf("\n");
}


void compare_results_tan_ulp_err_signed2(double *x, double *y, int64_t *ulp_err, size_t n)
{
    for (size_t i = 0; i < n; i++) {
      ulp_err[i] = ulp_error_signed_double(reference_calculation2(x[i]),
					   reference_calculation(x[i]));
    }
}
