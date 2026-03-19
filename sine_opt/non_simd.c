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

double cos_taylor30(double x) {
    double x2 = x * x;

    double r =
    0x1.0000000000000p+0 +
    x2 * (
    -0x1.0000000000000p-1 +
    x2 * (
     0x1.5555555555555p-5 +
    x2 * (
    -0x1.6c16c16c16c17p-10 +
    x2 * (
     0x1.a01a01a01a01ap-16 +
    x2 * (
    -0x1.27e4fb7789f5cp-22 +
    x2 * (
     0x1.1eed8eff8d898p-29 +
    x2 * (
    -0x1.93974a81d1e3cp-37 +
    x2 * (
     0x1.4f8b588e368f1p-45 +
    x2 * (
    -0x1.1f5d7d9a7b0aap-54 +
    x2 * (
     0x1.0a6c0f5c0b56dp-63 +
    x2 * (
    -0x1.3b3b3b3b3b3b4p-73 +
    x2 * (
     0x1.3b3b3b3b3b3b4p-84 +
    x2 * (
    -0x1.1111111111111p-95 +
    x2 * (
     0x1.745d1745d1746p-108
    ))))))))))))));

    return r;
}

double sin_taylor29(double x)
{
    double x2 = x * x;

    double r =
    0x1.0000000000000p+0 +
    x2 * (
    -0x1.5555555555555p-3 +
    x2 * (
     0x1.1111111111111p-7 +
    x2 * (
    -0x1.a01a01a01a01ap-13 +
    x2 * (
     0x1.71de3a556c734p-19 +
    x2 * (
    -0x1.ae64567f544e4p-26 +
    x2 * (
     0x1.6124613a86d09p-33 +
    x2 * (
    -0x1.ae7f3e733b81fp-41 +
    x2 * (
     0x1.3c6a5b63d3f1fp-49 +
    x2 * (
    -0x1.6c16c16c16c17p-58 +
    x2 * (
     0x1.1eed8eff8d898p-67 +
    x2 * (
    -0x1.27e4fb7789f5cp-77 +
    x2 * (
     0x1.93974a81d1e3cp-88 +
    x2 * (
    -0x1.4f8b588e368f1p-100
    )))))))))))));

    return x * r;
}


int main(void) {
    const size_t n = 100000;
    const double a = 0;
    const double b = M_PI_2;

    double input[n];
    double res[n];

    for (size_t i = 0; i < n; ++i) {
        input[i] = a + ((double)i / (double)(n - 1)) * (b - a);
    }

    for (size_t i = 0; i < n; i++) {
        res[i] = sin_taylor29(input[i]);
    }

    printf("%24s %24s %24s %12s\n", "x", "vfast_sin(x)", "ref_double", "ulp_diff");
    int count_one = 0;
    int count_two = 0;
    int count_severe = 0;

    for (size_t i = 0; i < n; ++i) {
        double ref = reference_sin_as_double(input[i]);
        long long ulp = signed_ulp_diff(res[i], ref);

        if (ulp > 2 || ulp < -2) count_severe++;
        if (ulp == 2 || ulp == -2) count_two++;
        if (ulp == 1 || ulp == -1) count_one++;

        //printf("%24.17g %24.17g %24.17g %12lld\n", input[i], res[i], ref, ulp);
    }
    printf("Total ULP 1 Values: %d\n\n", count_one);
    printf("Total ULP 2 Values: %d\n\n", count_two);
    printf("Total Unacceptable: %d\n\n", count_severe);


}
