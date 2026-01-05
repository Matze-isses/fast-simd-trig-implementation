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


/*
 * Compare y[i] against a high-precision reduction of x[i] modulo (pi/2),
 * using nearest-integer quotient (so the remainder is centered around 0).
 *
 * That is:
 *    k = round( x / (pi/2) )
 *    r = x - k*(pi/2)
 */
void compare_results(double *x, double *y, double *left_over, size_t n) {
    const slong prec = 512;

    arb_t pi, pi_over_2;
    arb_t ax, t, tmp, r;
    arb_t ay, err;          /* y[i] lifted, and high-precision error */
    fmpz_t k;

    arb_init(pi);
    arb_init(pi_over_2);
    arb_init(ax);
    arb_init(t);
    arb_init(tmp);
    arb_init(r);
    arb_init(ay);
    arb_init(err);
    fmpz_init(k);

    arb_const_pi(pi, prec);
    arb_mul_2exp_si(pi_over_2, pi, -1);   /* pi/2 */

    for (size_t i = 0; i < n; i++) {
        /* ax = x[i] exactly as that double */
        arb_set_d(ax, x[i]);

        /* t = x / (pi/2) */
        arb_div(t, ax, pi_over_2, prec);

        /* k = floor(t) (match your SIMD FLOOR) */
        arf_get_fmpz(k, arb_midref(t), ARF_RND_FLOOR);

        /* r = x - k*(pi/2) */
        arb_mul_fmpz(tmp, pi_over_2, k, prec);
        arb_sub(r, ax, tmp, prec);

        /* ay = y[i] exactly as that double */
        arb_set_d(ay, y[i]);

        /* err = ay - r  (done in high precision) */
        arb_sub(err, ay, r, prec);

        /* Convert the (midpoint of) high-precision error to double */
        left_over[i] = arf_get_d(arb_midref(err), ARF_RND_NEAR);
    }

    arb_clear(pi);
    arb_clear(pi_over_2);
    arb_clear(ax);
    arb_clear(t);
    arb_clear(tmp);
    arb_clear(r);
    arb_clear(ay);
    arb_clear(err);
    fmpz_clear(k);
}


void do_range_reduction(double *input, double *res, size_t n) {
  /*
   * Your original CORRECTION constant was tuned for 2*pi reduction.
   * A reasonable first-order adaptation is to scale it by the same factor as the range:
   *   (pi/2) = (2*pi)/4  => correction_pi_over_2 ≈ CORRECTION / 4
   *
   * If you have a derivation for CORRECTION, you should recompute it for pi/2 directly.
   */
  // const double CORRECTION_A = 3.8981718325193755e-17;
  const double CORRECTION_A = 0;
  const double CORRECTION_B = 3.8981718325193756e-18;

  const SDOUBLE pi_2 = LOAD_DOUBLE(M_PI_2);                   /* 1 / (pi/2) */
  const SDOUBLE one_over_pi_2 = LOAD_DOUBLE(2.0 / M_PI);    /* pi/2 */

  const SDOUBLE correction_a  = LOAD_DOUBLE(CORRECTION_A);
  const SDOUBLE correction_b  = LOAD_DOUBLE(CORRECTION_B);


  for (int i = 0; i < (int)n; i += 4) {
    SDOUBLE x = LOAD_DOUBLE_VEC(&input[i]);

    const SDOUBLE ranges_away     = MUL_DOUBLE_S(x, one_over_pi_2);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple  = MUL_DOUBLE_S(num_ranges_away, pi_2);

    SDOUBLE in_outer_range   = SUB_DOUBLE_S(x, range_multiple);
    SDOUBLE correction = MUL_DOUBLE_S(x, correction_a);
    in_outer_range = SUB_DOUBLE_S(in_outer_range, correction);

    SIMD_TO_DOUBLE_VEC(&res[i], in_outer_range);
  }
}


int main(int argc, char *argv[]) {
  (void)argc; (void)argv;

  printf("Hello World\n");

  int n = 20;

  double *test_vals = (double*)malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    test_vals[i] = 1.0 * (double)i * M_PI_2 + 0.1;
  }

  double *res = (double*)malloc((size_t)n * sizeof(double));
  double *left_over = (double*)malloc((size_t)n * sizeof(double));

  do_range_reduction(test_vals, res, (size_t)n);
  compare_results(test_vals, res, left_over, (size_t)n);

  for (int i = 0; i < n; i++) {
    printf("ORIGIN: %.17g; RESULT: %.17g; LEFT-OVER: %.17g\n",
           test_vals[i], res[i], left_over[i]);
  }

  free(test_vals);
  free(res);
  free(left_over);

  return 0;
}
