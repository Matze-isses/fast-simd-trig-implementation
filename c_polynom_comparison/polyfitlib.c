// polyfitlib.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <string.h>
#include <omp.h>

void generate_random_angles(double *values, size_t N) {
    for (size_t i = 0; i < N; i++) {
        double u = (double) rand() / (double) RAND_MAX;
        values[i] = u * 0.5 * M_PI;
    }
}

void apply_sin_highprec(const double *values, long double *out_ld, size_t N) {

    #pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
        out_ld[i] = sinl((long double)values[i]);
    }
}

long double average_abs_deviation_ld_d(const long double *ref_ld,
                                       const double *other_d,
                                       size_t N) {
    if (N == 0) return 0.0L;
    long double acc = 0.0L;

    for (size_t i = 0; i < N; i++) {
        long double diff = ref_ld[i] - (long double)other_d[i];
        acc += fabsl(diff);
    }
    return acc / (long double)N;
}

void calculate_polynomial_values(double *target, double *x, size_t n, double *coeff, int num_coefficiants) {
    const double pi_over_4 = M_PI / 4.0;

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        target[i] = 0.0;
        for (int j = 0; j < num_coefficiants; j++) {
            target[i] += coeff[j] * pow(x[i] - pi_over_4, j);
        }   
    }
}

/* Exposed function to call from Python */
double compute_mae(int num_polynomial_degree, const double *polynomial_coefficiants) {
    srand((unsigned int) time(NULL));

    size_t N = 500000;

    double *values = malloc(N * sizeof(double));
    double *result = malloc(N * sizeof(double));
    long double *correct_values = malloc(N * sizeof(long double));

    if (!values || !result || !correct_values) {
        if (values) free(values);
        if (result) free(result);
        if (correct_values) free(correct_values);
        return NAN;
    }

    generate_random_angles(values, N);
    apply_sin_highprec(values, correct_values, N);

    /* cast away const only to match existing signature of calculate_polynomial_values */
    calculate_polynomial_values(result, values, N, (double*)polynomial_coefficiants, num_polynomial_degree);

    long double mae_ld = average_abs_deviation_ld_d(correct_values, result, N);

    free(result);
    free(values);
    free(correct_values);

    return (double)mae_ld;
} // gcc -O3 -fopenmp -shared -fPIC polyfitlib.c -o libpolyfit.so -lm
