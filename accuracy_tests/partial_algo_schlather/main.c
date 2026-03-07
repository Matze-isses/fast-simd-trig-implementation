// Compilation and Running: 
//
// gcc ./vfast_tan_avx512.c ./trig_arb_comparison.c ./main.c -O2 -Wextra -lm -mavx512f -mavx512bw -mavx512dq -mavx512vl -mfma -fopenmp -pedantic -o test && ./test 1000000
//
//
//
// Plotting (requires some python packages!):
//
// python ./plotting_stuff.py
//
//
// Data is also saved in data.json

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
//#include "trig_simd.h"



void write_data_json(double *x, int64_t *ulp_error, int n)
{
    FILE *f = fopen("data.json", "w");
    if (!f) {
        perror("fopen(data.json)");
        return;
    }

    double cum_ulp_error = 0.0;

    // JSON: array of [x, ulp_error] pairs
    fputs("[\n", f);
    for (int i = 0; i < n; i++) {
        cum_ulp_error += fabs((double)ulp_error[i]);

        // print comma only between elements (no trailing comma)
        fprintf(f, "  [%0.17g, %ld]%s\n",
                x[i], ulp_error[i], (i + 1 < n) ? "," : "");
    }
    fputs("]\n", f);

    fclose(f);

    // optional: print summary to stdout
    printf("Wrote %d points to data.json\n", n);
    printf("cum_ulp_error = %.17g\n", cum_ulp_error);
 }


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
    int64_t *ulp_error2 = (int64_t*)malloc((size_t)n * sizeof(int64_t));

    int *lsb = (int*)malloc((size_t)n * sizeof(int));


    const double a = M_PI_4 * 1.5 - 0.01;  
    const double b = M_PI_2 - 0.01;
    //const double b = ;

    // const double a = 3 * M_PI / 8.0;
    // const double b = nextafter(M_PI_2, 0);

    //    srand((unsigned)time(NULL)); // martin temporaer

    for (int i = 0; i < n; ++i) {
        double u = (double)(i)/(double)(n);
        x[i] = a + (b - a) * u;
    }


    vfast_tan(x, res, lsb, n);
    compare_results_tan_ulp_err_signed(x, res, ulp_error, n);
    compare_results_tan_ulp_err_signed2(x, res, ulp_error2, n);


    int cum_ulp_error = 0;
    int cum_ulp_error2 = 0;

    for (int i = 0; i < n; i++) {
        cum_ulp_error += abs((int)ulp_error[i]);
        cum_ulp_error2 += abs((int)ulp_error2[i]);
    } 

    write_data_json(x, ulp_error, n);
    // write_data_json(x, ulp_error2, n);

    printf("\nAVG ULP ERROR: %.17g\n", (double)cum_ulp_error/(double)n);
    //    printf("AVG ULP DOUBL: %.17g\n", (double)cum_ulp_error2/(double)n);

    return 0;
}
