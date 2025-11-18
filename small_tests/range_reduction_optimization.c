#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>


#include <stdio.h>
#include <stdlib.h>
#include "../trig_simd.h"
#include "flint/flint.h"
#include "flint/arb.h"
#include <flint/fmpz.h>

#include <string.h>   // just used for the error calculation


void compare_results(double *x, double *y, double *left_over, size_t n) {
    /* Choose a working precision (in bits) for the Arb computations. */
    const slong prec = 256;

    /* Arb variables */
    arb_t pi, two_pi;
    arb_t ax;    /* lifted x[i] */
    arb_t t;     /* x[i] / (2*pi) */
    arb_t tmp;   /* 2*pi * k */
    arb_t r;     /* remainder = x[i] - k * 2*pi */

    /* Integer for the quotient k */
    fmpz_t k;

    /* Initialize Arb and fmpz objects */
    arb_init(pi);
    arb_init(two_pi);
    arb_init(ax);
    arb_init(t);
    arb_init(tmp);
    arb_init(r);
    fmpz_init(k);

    /* two_pi = 2 * pi (computed once) */
    arb_const_pi(pi, prec);          /* pi at 'prec' bits precision */
    arb_mul_ui(two_pi, pi, 2, prec); /* two_pi = 2 * pi */

    for (size_t i = 0; i < n; i++) {
        /* ax = x[i] as an Arb number (exactly that double) */
        arb_set_d(ax, x[i]);

        /* t = ax / (2 * pi) */
        arb_div(t, ax, two_pi, prec);

        /*
         * k = nearest integer to t, using the midpoint of t:
         *    k ≈ round( x[i] / (2*pi) )
         */
        arf_get_fmpz(k, arb_midref(t), ARF_RND_NEAR);

        /* tmp = k * (2 * pi) */
        arb_mul_fmpz(tmp, two_pi, k, prec);

        /* r = ax - tmp = high-precision (x[i] mod 2*pi) in the [-pi,pi]-style sense */
        arb_sub(r, ax, tmp, prec);

        /* Get the midpoint of r as a double (correctly rounded) */
        double true_mod = arf_get_d(arb_midref(r), ARF_RND_NEAR);

        /* Store the error: your result minus the high-precision result */
        left_over[i] = y[i] - true_mod;
    }

    /* Clean up Arb / FLINT objects */
    arb_clear(pi);
    arb_clear(two_pi);
    arb_clear(ax);
    arb_clear(t);
    arb_clear(tmp);
    arb_clear(r);
    fmpz_clear(k);
}


void do_range_reduction(double *input, double *res, size_t n) {
  double CORRECTION = 3.8981718325193755e-17;

  const SDOUBLE one_over_2_pi = LOAD_DOUBLE(1 / (2.0 * M_PI));
  const SDOUBLE two_pi = LOAD_DOUBLE(2.0 * M_PI);
  const SDOUBLE correction = LOAD_DOUBLE(CORRECTION);

  for (int i = 0; i < (int)n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_2_pi);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, two_pi);

    SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple);
    SDOUBLE correction_term = MUL_DOUBLE_S(x, correction);
    in_outer_range = SUB_DOUBLE_S(in_outer_range, correction_term);


    SIMD_TO_DOUBLE_VEC(&res[i], in_outer_range);
  }
}


int main(int argc, char *argv[]) {
  printf("Hello World\n");

  int n = 12;

  double *test_vals = malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    test_vals[i] = 1000000000000000000 * i * M_PI;
  }

  double *res = malloc(n * sizeof(double));
  double *left_over = malloc(n * sizeof(double));
  
  do_range_reduction(test_vals, res, n);
  compare_results(test_vals, res, left_over, n);
  for (int i = 0; i < n; i++) {
    printf("ORIGIN: %.17g; RESULT: %.17g; LEFT-OVER: %.17g\n", test_vals[i], res[i], left_over[i]);
  }
}
