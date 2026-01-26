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

#define PRINT_FULL_M256D(reg) do {                           \
    double vals[4];                                          \
    _mm256_storeu_pd(vals, (reg));                           \
    printf(#reg " = [%.19g, %.19g, %.19g, %.19g]\n",         \
           vals[0], vals[1], vals[2], vals[3]);              \
} while (0)

void compare_results(double *x, double *y, double *left_over, double reduction_value, size_t n) {
    const slong prec = 512;

    arb_t range_reduction;
    arb_t arb_x;
    arb_t t;
    arb_t tmp;
    arb_t r;

    fmpz_t k;

    arb_init(range_reduction);
    arb_init(arb_x);
    arb_init(t);
    arb_init(tmp);
    arb_init(r);

    fmpz_init(k);

    arb_const_pi(range_reduction, prec);

    /* divide by 2 → π/2 */
    arb_div_ui(range_reduction, range_reduction, 2, prec);
    // arb_set_d(range_reduction, reduction_value);



    for (size_t i = 0; i < n; i++) {
        arb_set_d(arb_x, x[i]);
        arb_div(t, arb_x, range_reduction, prec);

        arf_get_fmpz(k, arb_midref(t), ARF_RND_FLOOR);

        /* tmp = k * (2 * pi) */
        arb_mul_fmpz(tmp, range_reduction, k, prec);
        arb_sub(r, arb_x, tmp, prec);

        double true_mod = arf_get_d(arb_midref(r), ARF_RND_NEAR);
        left_over[i] = y[i] - true_mod;
    }

    /* Clean up Arb / FLINT objects */
    arb_clear(range_reduction);

    arb_clear(arb_x);
    arb_clear(t);
    arb_clear(tmp);
    arb_clear(r);
    fmpz_clear(k);
}


void do_range_reduction(double *input, double *res, double reduction_value, size_t n) {
  double CORRECTION = 0.0;

  const SDOUBLE one_over_reduction_value = LOAD_DOUBLE(1.0 / reduction_value);
  const SDOUBLE reduced_range = LOAD_DOUBLE(reduction_value);
  const SDOUBLE correction = LOAD_DOUBLE(CORRECTION);

  for (int i = 0; i < (int)n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_reduction_value);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, reduced_range);

    SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple);
    SDOUBLE correction_term = MUL_DOUBLE_S(x, correction);
    in_outer_range = SUB_DOUBLE_S(in_outer_range, correction_term);


    SIMD_TO_DOUBLE_VEC(&res[i], in_outer_range);
  }
}


int main(int argc, char *argv[]) {
  int n = 12;
  double reduction_value = M_PI/2;

  double *test_vals = malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    test_vals[i] = 1 * (i+1);
  }

  double *res = malloc(n * sizeof(double));
  double *left_over = malloc(n * sizeof(double));
  
  do_range_reduction(test_vals, res, reduction_value, n);
  compare_results(test_vals, res, left_over, reduction_value, n);


  double cum_leftover = 0;
  double rel_cum_leftover = 0;

  printf("\n");
  printf("+---------------------------+---------------------------+---------------------------+\n");
  printf("| Test Value                | Left Over                 | Ranges Away               |\n");
  printf("+---------------------------+---------------------------+---------------------------+\n");

  for (int i = 0; i < n; i++) {
    double num_away = test_vals[i] / reduction_value;

    printf("| %25.17e | %25.17e | %25.17e | \n", test_vals[i], left_over[i], num_away);
    cum_leftover += left_over[i];
    rel_cum_leftover += fabs(left_over[i] / test_vals[i]);
  }
  printf("+---------------------------+---------------------------+---------------------------+\n\n");

  double avg_leftover = cum_leftover/n;
  double rel_avg_leftover = rel_cum_leftover/n;

  printf("AVG Leftover: %.17g REL AVG Leftover: %.17g\n", avg_leftover, rel_avg_leftover);
}




