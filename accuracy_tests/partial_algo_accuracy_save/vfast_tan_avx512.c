#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include "trig_simd.h"



void apply_taylor_ld(const double *input, double *res, size_t n)
{
    // Coeffs for tan(x) ≈ x + c1*x^3 + c2*x^5 + ... + c14*x^29
    // Stored as long double hex literals.
    const long double c1  = 0x8p-3L;
    const long double c2  = 0xa.aaaaaaaaaaaaaabp-5L;
    const long double c3  = 0x8.888888888888889p-6L;
    const long double c4  = 0xd.d0dd0dd0dd0dd0ep-8L;
    const long double c5  = 0xb.327a4416087cf99p-9L;
    const long double c6  = 0x9.1371aaf3611e47bp-10L;
    const long double c7  = 0xe.b69e870abeefdbp-12L;
    const long double c8  = 0xb.ed1b2295baf15b5p-13L;
    const long double c9  = 0x9.aac12401b3a2291p-14L;
    const long double c10 = 0xf.abebb9a68b3210dp-16L;
    const long double c11 = 0xc.b3f0c57e57d6451p-17L;
    const long double c12 = 0xa.4bec7751292c99fp-18L;
    const long double c13 = 0x8.589969cd3028276p-19L;
    const long double c14 = 0xd.87b969f71274916p-21L;

    for (size_t i = 0; i < n; ++i) {
        long double x  = (long double)input[i];
        long double x2 = x * x;

        // Evaluate: tan(x) ≈ x * (1 + c1*x2 + c2*x2^2 + ... + c14*x2^14)
        long double p =
            ((((((((((((((c14 * x2 + c13) * x2 + c12) * x2 + c11) * x2 + c10)
            * x2 + c9) * x2 + c8) * x2 + c7) * x2 + c6) * x2 + c5)
            * x2 + c4) * x2 + c3) * x2 + c2) * x2 + c1) * x2 + 1.0L);

        long double y = x * p;

        // One final rounding to nearest double (per current rounding mode, normally FE_TONEAREST)
        res[i] = (double)y;
    }
}



void vfast_tan(double *input, double *res, size_t n) {
    int simd_doubles = SIMD_LENGTH / 64;

    const SDOUBLE pi_2_hi = LOAD_DOUBLE(M_PI_2);
    const SDOUBLE pi_2_lo = LOAD_DOUBLE(0x1.1a62633145c07p-54);

    const SDOUBLE pi_4 = LOAD_DOUBLE(M_PI_4);
    const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);

    const SDOUBLE q2_bitshift_rr1 = LOAD_DOUBLE(pow(2, -54));
    const SDOUBLE q2_bitshift_rr2 = LOAD_DOUBLE(-pow(2, -55));



    const SDOUBLE neg_one = LOAD_DOUBLE(-1.0);
    const SDOUBLE neg_half = LOAD_DOUBLE(-0.5);

    const SDOUBLE zero = _mm512_setzero_pd();
    const SDOUBLE zero_two_five = LOAD_DOUBLE(0.25);

    const SDOUBLE half = LOAD_DOUBLE(0.5);
    const SDOUBLE one = LOAD_DOUBLE(1.0);
    const SDOUBLE two = LOAD_DOUBLE(2.0);
    const SDOUBLE three = LOAD_DOUBLE(3.0);

    const SDOUBLE taylor_coeff1  = LOAD_DOUBLE(tan_tp1);
    const SDOUBLE taylor_coeff2  = LOAD_DOUBLE(tan_tp2);
    const SDOUBLE taylor_coeff3  = LOAD_DOUBLE(tan_tp3);
    const SDOUBLE taylor_coeff4  = LOAD_DOUBLE(tan_tp4);
    const SDOUBLE taylor_coeff5  = LOAD_DOUBLE(tan_tp5);
    const SDOUBLE taylor_coeff6  = LOAD_DOUBLE(tan_tp6);
    const SDOUBLE taylor_coeff7  = LOAD_DOUBLE(tan_tp7);
    const SDOUBLE taylor_coeff8  = LOAD_DOUBLE(tan_tp8);
    const SDOUBLE taylor_coeff9  = LOAD_DOUBLE(tan_tp9);
    const SDOUBLE taylor_coeff10 = LOAD_DOUBLE(tan_tp10);
    const SDOUBLE taylor_coeff11 = LOAD_DOUBLE(tan_tp11);
    const SDOUBLE taylor_coeff12 = LOAD_DOUBLE(tan_tp12);
    const SDOUBLE taylor_coeff13 = LOAD_DOUBLE(tan_tp13);
    
    for (int i = 0; i < (int) n; i += SIMD_DOUBLES) {
        SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

        const SDOUBLE from_behind = SUB_DOUBLE_S(pi_2_hi, x);
        const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
        const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);

        /* obtaining bool vectors for the each quadrant */
        MASK8 m0 = CMP_MASK(quadrant, zero,  _CMP_EQ_OQ);
        MASK8 m1 = CMP_MASK(quadrant, one,   _CMP_EQ_OQ);
        MASK8 m2 = CMP_MASK(quadrant, two,   _CMP_EQ_OQ);
        MASK8 m3 = CMP_MASK(quadrant, three, _CMP_EQ_OQ);

        x = MASK_MUL_PD(x, m1, x, half); // Without loss

        x = MASK_SUB_PD(x, m2, pi_2_hi, x);
        x = MASK_MUL_PD(x, m2, x, half);

        // mit der aktion haut auch q2 hin
        x = MASK_ADD_PD(x, m2, x, q2_bitshift_rr1);

        __mmask8 mx_lt_quarter = _mm512_cmp_pd_mask(x, zero_two_five, _CMP_LT_OQ);
        __mmask8 m2_and_x_lt_quarter = _kand_mask8(m2, mx_lt_quarter);

        x = MASK_ADD_PD(x, m2_and_x_lt_quarter, x, q2_bitshift_rr2);

        // mit dem ding ist auch q3 exakt
        x = MASK_SUB_PD(x, m3, pi_2_hi, x);
        x = MASK_ADD_PD(x, m3, x, pi_2_lo);

        /* ---- Taylor Loop ---- */
        const SDOUBLE x_square = MUL_DOUBLE_S(x, x);

        const SDOUBLE result_q0_t1  = FMADD_PD(taylor_coeff13, x_square, taylor_coeff12);
        const SDOUBLE result_q0_t2  = FMADD_PD(result_q0_t1, x_square, taylor_coeff11);
        const SDOUBLE result_q0_t3  = FMADD_PD(result_q0_t2, x_square, taylor_coeff10);
        const SDOUBLE result_q0_t4  = FMADD_PD(result_q0_t3, x_square, taylor_coeff9);
        const SDOUBLE result_q0_t5  = FMADD_PD(result_q0_t4, x_square, taylor_coeff8);
        const SDOUBLE result_q0_t6  = FMADD_PD(result_q0_t5, x_square, taylor_coeff7);
        const SDOUBLE result_q0_t7  = FMADD_PD(result_q0_t6, x_square, taylor_coeff6);
        const SDOUBLE result_q0_t8  = FMADD_PD(result_q0_t7, x_square, taylor_coeff5);
        const SDOUBLE result_q0_t9  = FMADD_PD(result_q0_t8, x_square, taylor_coeff4);
        const SDOUBLE result_q0_t10 = FMADD_PD(result_q0_t9, x_square, taylor_coeff3);
        const SDOUBLE result_q0_t11 = FMADD_PD(result_q0_t10, x_square, taylor_coeff2);
        const SDOUBLE result_q0_t12 = FMADD_PD(result_q0_t11, x_square, taylor_coeff1);
        const SDOUBLE result_q0_1 = FMADD_PD(result_q0_t12, x_square, one);

        const SDOUBLE result_q0 = MUL_DOUBLE_S(result_q0_1, x);

        SIMD_TO_DOUBLE_VEC(&res[i], result_q0);
    }

    /* Treatment of the left overs with glibc */
    int num_left_over = (n % simd_doubles);
    if (num_left_over != 0) {printf("[WARNING] Test got wrong input number");}
}


void safe_vfast_tan(double *input, double *res, size_t n, double error_threshold) {
  // based on the error of the interval [pi/4, 3 pi/8] 
  // because there it is the largest and is not at singularity
  double maximum_error_coeff = 2.62e-16;

  for (size_t i = 0; i < n; i++) {
    double error_of_element = maximum_error_coeff * fabs(input[i]);
    if (error_of_element > error_threshold) {
      printf("[Warning] Tan calculation exceeds error threshold for Element at index %d!\n");
    }
  }
  vfast_tan(input, res, n);
}
