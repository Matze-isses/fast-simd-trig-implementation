#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <stdint.h>

#include <stdlib.h>
#include "trig_simd.h"

#include <string.h>


void vfast_tan(double *input, double *res, size_t n) {
    int simd_doubles = SIMD_LENGTH / 64;

    SET1_PD(pi_2, M_PI_2);
    SET1_PD(one_over_pi_8, 1/M_PI_8);
    SET1_PD(one_over_pi_2, 1/M_PI_2);

    SET1_PD(pi_2_hi, M_PI_2);
    SET1_PD(cor_coeff, COR_COEFF);
    SET1_PD(one, 1.0);

    SET1_PD(taylor_coeff1, tan_tp1);
    SET1_PD(taylor_coeff2, tan_tp2);
    SET1_PD(taylor_coeff3, tan_tp3);
    SET1_PD(taylor_coeff4, tan_tp4);
    SET1_PD(taylor_coeff5, tan_tp5);
    SET1_PD(taylor_coeff6, tan_tp6);
    SET1_PD(taylor_coeff7, tan_tp7);
    SET1_PD(taylor_coeff8, tan_tp8);
    SET1_PD(taylor_coeff9, tan_tp9);
    SET1_PD(taylor_coeff10, tan_tp10);
    SET1_PD(taylor_coeff11, tan_tp11);
    SET1_PD(taylor_coeff12, tan_tp12);
    SET1_PD(taylor_coeff13, tan_tp13);

  for (int i = 0; i < (int) n; i += SIMD_DOUBLES) {
    LOAD_DOUBLE_VEC(x_in, &input[i]);

    /* Range Reduction */
    MUL_DOUBLE_S(ranges_away, x_in, one_over_pi_2);
    FLOOR_DOUBLE_S(num_ranges_away, ranges_away);
    MUL_DOUBLE_S(range_multiple, num_ranges_away, pi_2);
    SUB_DOUBLE_S(x_reduced_range, x_in, range_multiple);

    /* treatment of cotan ranges */
    GEN_MASK_IF_ODD(odd_mask, num_ranges_away);
    MASK_SUB_PD(x_cotan_adjust, x_reduced_range, odd_mask, pi_2, x_reduced_range);

    /* Get Quadrant */
    MUL_DOUBLE_S(not_floored, x_cotan_adjust, one_over_pi_8);
    FLOOR_DOUBLE_S(quadrant, not_floored);

    /* Generate Mask */
    CMP_MASK(m2, quadrant, one, _CMP_GT_OQ);

    // Kann kein fehler haben da nur exponent geaendert wird
    MASK_SUB_PD(x_second_adjust, x_cotan_adjust, m2, pi_2_hi, x_cotan_adjust);
    HALF_PD_FAST(x_half, x_second_adjust);

    /* ---- Taylor Loop ---- */
    MUL_DOUBLE_S(x_square, x_half, x_half);

    FMADD_PD(result_q0_t1, taylor_coeff13, x_square, taylor_coeff12);
    FMADD_PD(result_q0_t2, result_q0_t1, x_square, taylor_coeff11);
    FMADD_PD(result_q0_t3, result_q0_t2, x_square, taylor_coeff10);
    FMADD_PD(result_q0_t4, result_q0_t3, x_square, taylor_coeff9);
    FMADD_PD(result_q0_t5, result_q0_t4, x_square, taylor_coeff8);
    FMADD_PD(result_q0_t6, result_q0_t5, x_square, taylor_coeff7);
    FMADD_PD(result_q0_t7, result_q0_t6, x_square, taylor_coeff6);
    FMADD_PD(result_q0_t8, result_q0_t7, x_square, taylor_coeff5);
    FMADD_PD(result_q0_t9, result_q0_t8, x_square, taylor_coeff4);
    FMADD_PD(result_q0_t10, result_q0_t9, x_square, taylor_coeff3);
    FMADD_PD(result_q0_t11, result_q0_t10, x_square, taylor_coeff2);
    FMADD_PD(result_q0_t12, result_q0_t11, x_square, taylor_coeff1);

    MUL_DOUBLE_S(result_q0_t13, result_q0_t12, x_square); // Note that here a one is Missing
                                                                         
    /* Calculate Correction */
    MUL_DOUBLE_S(result_q0_partial, result_q0_t13, x_half);
    FMADD_PD(correction, x_square, cor_coeff, cor_coeff);
    MASK_ADD_PD(result_q0_corrected, result_q0_partial, m2, result_q0_partial, correction);

    /* Obtain result */
    ADD_DOUBLE_S(result_q0, result_q0_corrected, x_half); // Here the one is added

    /* Getting Values for Double Angle */
    DOUBLE_PD_FAST(nominator, result_q0);
    MUL_DOUBLE_S(result_q0_square, result_q0, result_q0);
    SUB_DOUBLE_S(denominator, one, result_q0_square);

    /* Obtaining the interval results */
    DIV_DOUBLE_S(result_q01, nominator, denominator);
    DIV_DOUBLE_S(result_q23, denominator, nominator);

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
