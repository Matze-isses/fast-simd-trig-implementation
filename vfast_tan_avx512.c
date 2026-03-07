#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include "trig_simd.h"


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

  
  for (int i = 0; i < (int) n; i += 4 * SIMD_DOUBLES) {
    LOAD_DOUBLE_VEC(x_in_l0, &input[i]);
    LOAD_DOUBLE_VEC(x_in_l1, &input[i+SIMD_DOUBLES]);
    LOAD_DOUBLE_VEC(x_in_l2, &input[i+2*SIMD_DOUBLES]);
    LOAD_DOUBLE_VEC(x_in_l3, &input[i+3*SIMD_DOUBLES]);

    /* Range Reduction */
    MUL_DOUBLE_S(ranges_away_l0, x_in_l0, one_over_pi_2);
    FLOOR_DOUBLE_S(num_ranges_away_l0, ranges_away_l0);
    MUL_DOUBLE_S(range_multiple_l0, num_ranges_away_l0, pi_2);

    MUL_DOUBLE_S(ranges_away_l1, x_in_l1, one_over_pi_2);
    FLOOR_DOUBLE_S(num_ranges_away_l1, ranges_away_l1);
    MUL_DOUBLE_S(range_multiple_l1, num_ranges_away_l1, pi_2);

    MUL_DOUBLE_S(ranges_away_l2, x_in_l2, one_over_pi_2);
    FLOOR_DOUBLE_S(num_ranges_away_l2, ranges_away_l2);
    MUL_DOUBLE_S(range_multiple_l2, num_ranges_away_l2, pi_2);

    MUL_DOUBLE_S(ranges_away_l3, x_in_l3, one_over_pi_2);
    FLOOR_DOUBLE_S(num_ranges_away_l3, ranges_away_l3);
    MUL_DOUBLE_S(range_multiple_l3, num_ranges_away_l3, pi_2);

    SUB_DOUBLE_S(x_reduced_range_l0, x_in_l0, range_multiple_l0);
    SUB_DOUBLE_S(x_reduced_range_l1, x_in_l1, range_multiple_l1);
    SUB_DOUBLE_S(x_reduced_range_l2, x_in_l2, range_multiple_l2);
    SUB_DOUBLE_S(x_reduced_range_l3, x_in_l3, range_multiple_l3);

    /* treatment of cotan ranges */
    GEN_MASK_IF_ODD(odd_mask_l0, num_ranges_away_l0);
    GEN_MASK_IF_ODD(odd_mask_l1, num_ranges_away_l1);
    GEN_MASK_IF_ODD(odd_mask_l2, num_ranges_away_l2);
    GEN_MASK_IF_ODD(odd_mask_l3, num_ranges_away_l3);

    MASK_SUB_PD(x_cotan_adjust_l0, x_reduced_range_l0, odd_mask_l0, pi_2, x_reduced_range_l0);
    MASK_SUB_PD(x_cotan_adjust_l1, x_reduced_range_l1, odd_mask_l1, pi_2, x_reduced_range_l1);
    MASK_SUB_PD(x_cotan_adjust_l2, x_reduced_range_l2, odd_mask_l2, pi_2, x_reduced_range_l2);
    MASK_SUB_PD(x_cotan_adjust_l3, x_reduced_range_l3, odd_mask_l3, pi_2, x_reduced_range_l3);

    /* Get Quadrant */
    MUL_DOUBLE_S(not_floored_l0, x_cotan_adjust_l0, one_over_pi_8);
    MUL_DOUBLE_S(not_floored_l1, x_cotan_adjust_l1, one_over_pi_8);
    MUL_DOUBLE_S(not_floored_l2, x_cotan_adjust_l2, one_over_pi_8);
    MUL_DOUBLE_S(not_floored_l3, x_cotan_adjust_l3, one_over_pi_8);

    FLOOR_DOUBLE_S(quadrant_l0, not_floored_l0);
    FLOOR_DOUBLE_S(quadrant_l1, not_floored_l1);
    FLOOR_DOUBLE_S(quadrant_l2, not_floored_l2);
    FLOOR_DOUBLE_S(quadrant_l3, not_floored_l3);

    /* Generate Mask */
    CMP_MASK(m2_l0, quadrant_l0, one, _CMP_GT_OQ);
    CMP_MASK(m2_l1, quadrant_l1, one, _CMP_GT_OQ);
    CMP_MASK(m2_l2, quadrant_l2, one, _CMP_GT_OQ);
    CMP_MASK(m2_l3, quadrant_l3, one, _CMP_GT_OQ);

    // Kann kein fehler haben da nur exponent geaendert wird
    MASK_SUB_PD(x_second_adjust_l0, x_cotan_adjust_l0, m2_l0, pi_2_hi, x_cotan_adjust_l0);
    MASK_SUB_PD(x_second_adjust_l1, x_cotan_adjust_l1, m2_l1, pi_2_hi, x_cotan_adjust_l1);
    MASK_SUB_PD(x_second_adjust_l2, x_cotan_adjust_l2, m2_l2, pi_2_hi, x_cotan_adjust_l2);
    MASK_SUB_PD(x_second_adjust_l3, x_cotan_adjust_l3, m2_l3, pi_2_hi, x_cotan_adjust_l3);

    HALF_PD_FAST(x_half_l0, x_second_adjust_l0);
    HALF_PD_FAST(x_half_l1, x_second_adjust_l1);
    HALF_PD_FAST(x_half_l2, x_second_adjust_l2);
    HALF_PD_FAST(x_half_l3, x_second_adjust_l3);

    /* ---- Taylor Loop ---- */
    MUL_DOUBLE_S(x_square_l0, x_half_l0, x_half_l0);
    MUL_DOUBLE_S(x_square_l1, x_half_l1, x_half_l1);
    MUL_DOUBLE_S(x_square_l2, x_half_l2, x_half_l2);
    MUL_DOUBLE_S(x_square_l3, x_half_l3, x_half_l3);
    
    // Sortierung überarbeiten!!!
    FMADD_PD(result_q0_t1_l0, taylor_coeff13, x_square_l0, taylor_coeff12);
    FMADD_PD(result_q0_t2_l0, result_q0_t1_l0, x_square_l0, taylor_coeff11);
    FMADD_PD(result_q0_t1_l1, taylor_coeff13, x_square_l1, taylor_coeff12);

    FMADD_PD(result_q0_t2_l1, result_q0_t1_l1, x_square_l1, taylor_coeff11);
    FMADD_PD(result_q0_t3_l0, result_q0_t2_l0, x_square_l0, taylor_coeff10);
    FMADD_PD(result_q0_t1_l2, taylor_coeff13, x_square_l2, taylor_coeff12);

    FMADD_PD(result_q0_t2_l2, result_q0_t1_l2, x_square_l2, taylor_coeff11);
    FMADD_PD(result_q0_t3_l1, result_q0_t2_l1, x_square_l1, taylor_coeff10);
    FMADD_PD(result_q0_t1_l3, taylor_coeff13, x_square_l3, taylor_coeff12);

    FMADD_PD(result_q0_t4_l0, result_q0_t3_l0, x_square_l0, taylor_coeff9);
    FMADD_PD(result_q0_t3_l2, result_q0_t2_l2, x_square_l2, taylor_coeff10);
    FMADD_PD(result_q0_t2_l3, result_q0_t1_l3, x_square_l3, taylor_coeff11);

    FMADD_PD(result_q0_t5_l0, result_q0_t4_l0, x_square_l0, taylor_coeff8);
    FMADD_PD(result_q0_t4_l1, result_q0_t3_l1, x_square_l1, taylor_coeff9);
    FMADD_PD(result_q0_t3_l3, result_q0_t2_l3, x_square_l3, taylor_coeff10);

    FMADD_PD(result_q0_t6_l0, result_q0_t5_l0, x_square_l0, taylor_coeff7);
    FMADD_PD(result_q0_t5_l1, result_q0_t4_l1, x_square_l1, taylor_coeff8);
    FMADD_PD(result_q0_t4_l2, result_q0_t3_l2, x_square_l2, taylor_coeff9);

    FMADD_PD(result_q0_t7_l0, result_q0_t6_l0, x_square_l0, taylor_coeff6);
    FMADD_PD(result_q0_t5_l2, result_q0_t4_l2, x_square_l2, taylor_coeff8);
    FMADD_PD(result_q0_t6_l1, result_q0_t5_l1, x_square_l1, taylor_coeff7);

    FMADD_PD(result_q0_t4_l3, result_q0_t3_l3, x_square_l3, taylor_coeff9);
    FMADD_PD(result_q0_t6_l2, result_q0_t5_l2, x_square_l2, taylor_coeff7);
    FMADD_PD(result_q0_t7_l1, result_q0_t6_l1, x_square_l1, taylor_coeff6);

    FMADD_PD(result_q0_t8_l0, result_q0_t7_l0, x_square_l0, taylor_coeff5);
    FMADD_PD(result_q0_t5_l3, result_q0_t4_l3, x_square_l3, taylor_coeff8);
    FMADD_PD(result_q0_t7_l2, result_q0_t6_l2, x_square_l2, taylor_coeff6);

    FMADD_PD(result_q0_t6_l3, result_q0_t5_l3, x_square_l3, taylor_coeff7);
    FMADD_PD(result_q0_t8_l1, result_q0_t7_l1, x_square_l1, taylor_coeff5);
    FMADD_PD(result_q0_t9_l0, result_q0_t8_l0, x_square_l0, taylor_coeff4);

    FMADD_PD(result_q0_t7_l3, result_q0_t6_l3, x_square_l3, taylor_coeff6);
    FMADD_PD(result_q0_t8_l2, result_q0_t7_l2, x_square_l2, taylor_coeff5);
    FMADD_PD(result_q0_t9_l1, result_q0_t8_l1, x_square_l1, taylor_coeff4);

    FMADD_PD(result_q0_t8_l3, result_q0_t7_l3, x_square_l3, taylor_coeff5);
    FMADD_PD(result_q0_t9_l2, result_q0_t8_l2, x_square_l2, taylor_coeff4);
    FMADD_PD(result_q0_t10_l0, result_q0_t9_l0, x_square_l0, taylor_coeff3);

    FMADD_PD(result_q0_t9_l3, result_q0_t8_l3, x_square_l3, taylor_coeff4);
    FMADD_PD(result_q0_t10_l1, result_q0_t9_l1, x_square_l1, taylor_coeff3);
    FMADD_PD(result_q0_t11_l0, result_q0_t10_l0, x_square_l0, taylor_coeff2);

    FMADD_PD(result_q0_t10_l2, result_q0_t9_l2, x_square_l2, taylor_coeff3);
    FMADD_PD(result_q0_t11_l1, result_q0_t10_l1, x_square_l1, taylor_coeff2);
    FMADD_PD(result_q0_t10_l3, result_q0_t9_l3, x_square_l3, taylor_coeff3);

    FMADD_PD(result_q0_t11_l2, result_q0_t10_l2, x_square_l2, taylor_coeff2);
    FMADD_PD(result_q0_t11_l3, result_q0_t10_l3, x_square_l3, taylor_coeff2);
    FMADD_PD(result_q0_t12_l0, result_q0_t11_l0, x_square_l0, taylor_coeff1);

    FMADD_PD(result_q0_t12_l1, result_q0_t11_l1, x_square_l1, taylor_coeff1);
    FMADD_PD(result_q0_t12_l2, result_q0_t11_l2, x_square_l2, taylor_coeff1);
    FMADD_PD(result_q0_t12_l3, result_q0_t11_l3, x_square_l3, taylor_coeff1);


    MUL_DOUBLE_S(result_q0_t13_l0, result_q0_t12_l0, x_square_l0); // Note that here a one is Missing
    MUL_DOUBLE_S(result_q0_t13_l1, result_q0_t12_l1, x_square_l1);
    MUL_DOUBLE_S(result_q0_t13_l2, result_q0_t12_l2, x_square_l2);
    MUL_DOUBLE_S(result_q0_t13_l3, result_q0_t12_l3, x_square_l3);
                                                                         
    /* Calculate Correction */
    MUL_DOUBLE_S(result_q0_partial_l0, result_q0_t13_l0, x_half_l0);
    MUL_DOUBLE_S(result_q0_partial_l1, result_q0_t13_l1, x_half_l1);
    MUL_DOUBLE_S(result_q0_partial_l2, result_q0_t13_l2, x_half_l2);
    MUL_DOUBLE_S(result_q0_partial_l3, result_q0_t13_l3, x_half_l3);

    FMADD_PD(correction_l0, x_square_l0, cor_coeff, cor_coeff);
    FMADD_PD(correction_l1, x_square_l1, cor_coeff, cor_coeff);
    FMADD_PD(correction_l2, x_square_l2, cor_coeff, cor_coeff);
    FMADD_PD(correction_l3, x_square_l3, cor_coeff, cor_coeff);

    MASK_ADD_PD(result_q0_corrected_l0, result_q0_partial_l0, m2_l0, result_q0_partial_l0, correction_l0);
    MASK_ADD_PD(result_q0_corrected_l1, result_q0_partial_l1, m2_l1, result_q0_partial_l1, correction_l1);
    MASK_ADD_PD(result_q0_corrected_l2, result_q0_partial_l2, m2_l2, result_q0_partial_l2, correction_l2);
    MASK_ADD_PD(result_q0_corrected_l3, result_q0_partial_l3, m2_l3, result_q0_partial_l3, correction_l3);

    /* Obtain result */
    ADD_DOUBLE_S(result_q0_l0, result_q0_corrected_l0, x_half_l0); // Here the one is added
    ADD_DOUBLE_S(result_q0_l1, result_q0_corrected_l1, x_half_l1);
    ADD_DOUBLE_S(result_q0_l2, result_q0_corrected_l2, x_half_l2);
    ADD_DOUBLE_S(result_q0_l3, result_q0_corrected_l3, x_half_l3);

    /* Getting Values for Double Angle */
    DOUBLE_PD_FAST(nominator_l0, result_q0_l0);
    DOUBLE_PD_FAST(nominator_l1, result_q0_l1);
    DOUBLE_PD_FAST(nominator_l2, result_q0_l2);
    DOUBLE_PD_FAST(nominator_l3, result_q0_l3);

    MUL_DOUBLE_S(result_q0_square_l0, result_q0_l0, result_q0_l0);
    MUL_DOUBLE_S(result_q0_square_l1, result_q0_l1, result_q0_l1);
    MUL_DOUBLE_S(result_q0_square_l2, result_q0_l2, result_q0_l2);
    MUL_DOUBLE_S(result_q0_square_l3, result_q0_l3, result_q0_l3);

    SUB_DOUBLE_S(denominator_l0, one, result_q0_square_l0);
    SUB_DOUBLE_S(denominator_l1, one, result_q0_square_l1);
    SUB_DOUBLE_S(denominator_l2, one, result_q0_square_l2);
    SUB_DOUBLE_S(denominator_l3, one, result_q0_square_l3);

    /* Obtaining the interval results */
    DIV_DOUBLE_S(result_q01_l0, nominator_l0, denominator_l0);
    DIV_DOUBLE_S(result_q01_l1, nominator_l1, denominator_l1);
    DIV_DOUBLE_S(result_q01_l2, nominator_l2, denominator_l2);
    DIV_DOUBLE_S(result_q01_l3, nominator_l3, denominator_l3);

    DIV_DOUBLE_S(result_q23_l0, denominator_l0, nominator_l0);
    DIV_DOUBLE_S(result_q23_l1, denominator_l1, nominator_l1);
    DIV_DOUBLE_S(result_q23_l2, denominator_l2, nominator_l2);
    DIV_DOUBLE_S(result_q23_l3, denominator_l3, nominator_l3);

    /* Putting together results */
    MASK_MOV_PD(partial_result_l0, m2_l0, result_q01_l0, result_q23_l0);
    MASK_MOV_PD(partial_result_l1, m2_l1, result_q01_l1, result_q23_l1);
    MASK_MOV_PD(partial_result_l2, m2_l2, result_q01_l2, result_q23_l2);
    MASK_MOV_PD(partial_result_l3, m2_l3, result_q01_l3, result_q23_l3);

    /* adjust cotan ranges */
    FLIP_SIGN_IF_MASK_PD(result_l0, odd_mask_l0, partial_result_l0);
    FLIP_SIGN_IF_MASK_PD(result_l1, odd_mask_l1, partial_result_l1);
    FLIP_SIGN_IF_MASK_PD(result_l2, odd_mask_l2, partial_result_l2);
    FLIP_SIGN_IF_MASK_PD(result_l3, odd_mask_l3, partial_result_l3);

    /* save result */
    SIMD_TO_DOUBLE_VEC(&res[i], result_l0);
    SIMD_TO_DOUBLE_VEC(&res[i+SIMD_DOUBLES], result_l1);
    SIMD_TO_DOUBLE_VEC(&res[i+2*SIMD_DOUBLES], result_l2);
    SIMD_TO_DOUBLE_VEC(&res[i+3*SIMD_DOUBLES], result_l3);
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
