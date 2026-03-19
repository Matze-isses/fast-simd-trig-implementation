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

static double reference_cos_as_double(double x)
{
    long double xl = (long double)x;
    long double yl = cosl(xl);
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

double cos_taylor(double x)
{
    double r = 0;
    double x_square = x * x;
    double x_cube = x * x_square;

    r = cos_tp5;

//  r = r * x_square + cos_tp13;
//  r = r * x_square + cos_tp12;
//  r = r * x_square + cos_tp11;
//  r = r * x_square + cos_tp10;
//  r = r * x_square + cos_tp9;
//  r = r * x_square + cos_tp8;
    //r = r * x_square + cos_tp7;
    //r = r * x_square + cos_tp6;
    //r = r * x_square + cos_tp5;
    r = r * x_square + cos_tp4;
    r = r * x_square + cos_tp3;
    r = r * x_square + cos_tp2;
    r = r * x_square + cos_tp1;
    r = r * x_square + cos_tp0;

    return r;
}

double sin_taylor(double x)
{
    double r = 0;
    double x_square = x * x;
    double x_cube = x * x_square;

    r = sin_tp11 * x_square + sin_tp10;

    r = r * x_square + sin_tp9;
    r = r * x_square + sin_tp8;
    r = r * x_square + sin_tp7;
    r = r * x_square + sin_tp6;
    r = r * x_square + sin_tp5;
    r = r * x_square + sin_tp4;
    r = r * x_square + sin_tp3;
    r = r * x_square + sin_tp2;
    r = r * x_square + sin_tp1;

    r = r * x_cube + x;

    return r;
}


int main(void) {
    const size_t n = 100000;
    const double a = 0;
    const double b = 0.15658;

    double input[n];
    double res[n];

    for (size_t i = 0; i < n; ++i) {
        input[i] = a + ((double)i / (double)(n - 1)) * (b - a);
    }

    for (size_t i = 0; i < n; i++) {
        res[i] = cos_taylor(input[i]);
    }

    printf("%24s %24s %24s %12s\n", "x", "vfast_sin(x)", "ref_double", "ulp_diff");
    int count_one = 0;
    int count_two = 0;
    int count_severe = 0;

    for (size_t i = 0; i < n; ++i) {
        double ref = reference_cos_as_double(input[i]);
        long long ulp = signed_ulp_diff(res[i], ref);

        if (ulp > 2 || ulp < -2) count_severe++;
        if (ulp == 2 || ulp == -2) count_two++;
        if (ulp == 1 || ulp == -1) count_one++;

        //printf("%24.17g %24.17g %24.17g %12lld\n", input[i], res[i], ref, ulp);
    }

    printf("Total ULP 0 Values: %d\n", n - count_one - count_two - count_severe);
    printf("Total ULP 1 Values: %d\n", count_one);
    printf("Total ULP 2 Values: %d\n", count_two);
    printf("Total Unacceptable: %d\n\n", count_severe);


}
