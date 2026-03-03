#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>   // add this include

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>

#include "test_interface.h"
#include "trig_simd.h"

static inline uint64_t dbl_to_u64(double x) {
    uint64_t u;
    memcpy(&u, &x, sizeof u);
    return u;
}

static inline double u64_to_dbl(uint64_t u) {
    double x;
    memcpy(&x, &u, sizeof x);
    return x;
}

/* Change the raw IEEE-754 bit pattern by +/- delta (LSB steps). */
static inline double tweak_lsb(double x, int delta) {
    uint64_t u = dbl_to_u64(x);

    // Optional safety: skip NaN/Inf (all exponent bits = 1)
    if ((u & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL) {
        return x;
    }

    u = (uint64_t)((int64_t)u + (int64_t)delta);
    return u64_to_dbl(u);
}


int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n <= 0) {
        fprintf(stderr, "n must be a positive integer\n");
        return 1;
    }

    double *x = (double*)malloc((size_t)n * sizeof(double));
    double *res = (double*)malloc((size_t)n * sizeof(double));
    int64_t *ulp_error = (int64_t*)malloc((size_t)n * sizeof(int64_t));


    const double a = 0;
    const double b = M_PI_2;

    // const double a = 3 * M_PI / 8.0;
    // const double b = nextafter(M_PI_2, 0);

    srand((unsigned)time(NULL));

    for (int i = 0; i < n; ++i) {
        double u = (double)(i)/(double)(n-1);
        x[i] = a + (b - a) * u;
    }

    for (int i = 0; i < n; i++) {
        res[i] = 1 - x[i] * x[i];
    }

    compare_results_one_minus_x2_ulp_err_signed(x, res, ulp_error, n);
    int cum_ulp_error = 0;

    for (int i = 0; i < n; i++) {
        cum_ulp_error += ulp_error[i];
        printf("[%5.3f, %d],\n", x[i], ulp_error[i]);
    }

    printf("\nAVG ULP ERROR: %.17g\n", (double)cum_ulp_error/(double)n);

    return 0;
}
