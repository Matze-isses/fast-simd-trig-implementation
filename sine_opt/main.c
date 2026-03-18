#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "trig_simd.h"

static uint64_t double_to_ordered_u64(double x)
{
    uint64_t u;
    memcpy(&u, &x, sizeof(u));

    if (u & 0x8000000000000000ULL)
        return ~u;

    return u | 0x8000000000000000ULL;
}

static long long signed_ulp_diff(double approx, double ref)
{
    uint64_t ua = double_to_ordered_u64(approx);
    uint64_t ur = double_to_ordered_u64(ref);

    if (ua >= ur)
        return (long long)(ua - ur);

    return -(long long)(ur - ua);
}

static double reference_sin_as_double(double x)
{
    long double xl = (long double)x;
    long double yl = sinl(xl);
    return (double)yl;
}

static int save_results_tsv(const char *filename,
                            const double *input,
                            const double *res,
                            size_t n)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
        return -1;

    fprintf(fp, "x\tvfast_sin(x)\tref_double\tulp_diff\n");

    for (size_t i = 0; i < n; ++i) {
        double ref = reference_sin_as_double(input[i]);
        long long ulp = signed_ulp_diff(res[i], ref);

        fprintf(fp, "%.17g\t%.17g\t%.17g\t%lld\n",
                input[i], res[i], ref, ulp);
    }

    fclose(fp);
    return 0;
}

int main(void)
{
    const size_t n = 100000;
    const double a = 0;
    const double b = M_PI_2;

    double input[n];
    double res[n];

    for (size_t i = 0; i < n; ++i) {
        input[i] = a + ((double)i / (double)(n - 1)) * (b - a);
    }

    vfast_sin(input, res, n);

    printf("%24s %24s %24s %12s\n", "x", "vfast_sin(x)", "ref_double", "ulp_diff");

    for (size_t i = 0; i < n; ++i) {
        double ref = reference_sin_as_double(input[i]);
        long long ulp = signed_ulp_diff(res[i], ref);

        // printf("%24.17g %24.17g %24.17g %12lld\n", input[i], res[i], ref, ulp);
    }

    if (save_results_tsv("sin_results.tsv", input, res, n) != 0) {
        fprintf(stderr, "Could not write sin_results.tsv\n");
        return 1;
    }

    return 0;
}
