#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include "trig_simd.h"


void vfast_tan(double *input, double *res, size_t n) {
    int simd_doubles = SIMD_LENGTH / 64;

    const SDOUBLE pi_2          = LOAD_DOUBLE(M_PI_2);
    const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
    const SDOUBLE one_over_pi_2 = LOAD_DOUBLE(1/M_PI_2);

    const SDOUBLE pi_2_hi = LOAD_DOUBLE(M_PI_2);
    const SDOUBLE cor_coeff = LOAD_DOUBLE(COR_COEFF);
    const SDOUBLE one = LOAD_DOUBLE(1.0);

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

    /* Range Reduction */
    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_pi_2);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, pi_2);
    x = SUB_DOUBLE_S(x, range_multiple);

    /* treatment of cotan ranges */
    GEN_MASK_IF_ODD(odd_mask, num_ranges_away);
    x = MASK_SUB_PD(x, odd_mask, pi_2, x);

    /* Get Quadrant */
    const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
    const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);

    /* Generate Mask */
    MASK8 m2 = CMP_MASK(quadrant, one, _CMP_GT_OQ);

    // Kann kein fehler haben da nur exponent geaendert wird
    x = MASK_SUB_PD(x, m2, pi_2_hi, x);
    HALF_PD_FAST(x, x);

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
    const SDOUBLE result_q0_t13 = MUL_DOUBLE_S(result_q0_t12, x_square); // Note that here a one is Missing
                                                                         
    /* Calculate Correction */
    SDOUBLE result_q0_partial   = MUL_DOUBLE_S(result_q0_t13, x);
    SDOUBLE correction          = FMADD_PD(x_square, cor_coeff, cor_coeff);
    SDOUBLE result_q0_corrected = MASK_ADD_PD(result_q0_partial, m2, result_q0_partial, correction);

    /* Obtain result */
    SDOUBLE result_q0           = ADD_DOUBLE_S(result_q0_corrected, x); // Here the one is added

    /* Getting Values for Double Angle */
    DOUBLE_PD_FAST(nominator, result_q0);
    const SDOUBLE result_q0_square = MUL_DOUBLE_S(result_q0, result_q0);
    const SDOUBLE denominator      = SUB_DOUBLE_S(one, result_q0_square);

    /* Obtaining the interval results */
    const SDOUBLE result_q01 = DIV_DOUBLE_S(nominator, denominator);
    const SDOUBLE result_q23 = DIV_DOUBLE_S(denominator, nominator);

    /* Putting together results */
    MASK_MOV_PD(partial_result, m2, result_q01, result_q23);

    /* adjust cotan ranges */
    FLIP_SIGN_IF_MASK_PD(result, odd_mask, partial_result);

    /* save result */
    SIMD_TO_DOUBLE_VEC(&res[i], result);
  }

  /* Treatment of the left overs with glibc */
  int num_left_over = (n % simd_doubles);

  for (size_t i = n - num_left_over; i < n; i++) {
    res[i] = tan(input[i]);
  }
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
