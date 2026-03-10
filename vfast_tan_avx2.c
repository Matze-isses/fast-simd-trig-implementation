#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <stdint.h>

#include <stdlib.h>
#include "trig_simd.h"

#include <string.h>
#include "./util/bit_printing.h"


void vfast_tan(double *input, double *res, size_t n) {
    int simd_doubles = SIMD_LENGTH / 64;

    SET1_PD(pi_2, M_PI_2);
    SET1_PD(one_over_pi_8, 1/M_PI_8);
    SET1_PD(one_over_pi_2, 1/M_PI_2);

    SET1_PD(pi_2_hi, M_PI_2);
    SET1_PD(cor_coeff, COR_COEFF);
    SET1_PD(one, 1.0);

    SET1_PD(half, 0.5);
    SET1_PD(two, 2.0);

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

  for (int i = 0; i < (int) n; i += 5 * SIMD_DOUBLES) {
    LOAD_DOUBLE_VEC(x_in_0, &input[i]);
    LOAD_DOUBLE_VEC(x_in_1, &input[i+1*SIMD_DOUBLES]);
    LOAD_DOUBLE_VEC(x_in_2, &input[i+2*SIMD_DOUBLES]);
    LOAD_DOUBLE_VEC(x_in_3, &input[i+3*SIMD_DOUBLES]);
    LOAD_DOUBLE_VEC(x_in_4, &input[i+4*SIMD_DOUBLES]);

    /* Range Reduction */
    MUL_DOUBLE_S(ranges_away_0, x_in_0, one_over_pi_2);
    MUL_DOUBLE_S(ranges_away_1, x_in_1, one_over_pi_2);
    MUL_DOUBLE_S(ranges_away_2, x_in_2, one_over_pi_2);
    MUL_DOUBLE_S(ranges_away_3, x_in_3, one_over_pi_2);
    MUL_DOUBLE_S(ranges_away_4, x_in_4, one_over_pi_2);

    FLOOR_DOUBLE_S(num_ranges_away_0, ranges_away_0);
    FLOOR_DOUBLE_S(num_ranges_away_1, ranges_away_1);
    FLOOR_DOUBLE_S(num_ranges_away_2, ranges_away_2);
    FLOOR_DOUBLE_S(num_ranges_away_3, ranges_away_3);
    FLOOR_DOUBLE_S(num_ranges_away_4, ranges_away_4);

    MUL_DOUBLE_S(range_multiple_0, num_ranges_away_0, pi_2);
    MUL_DOUBLE_S(range_multiple_1, num_ranges_away_1, pi_2);
    MUL_DOUBLE_S(range_multiple_2, num_ranges_away_2, pi_2);
    MUL_DOUBLE_S(range_multiple_3, num_ranges_away_3, pi_2);
    MUL_DOUBLE_S(range_multiple_4, num_ranges_away_4, pi_2);

    SUB_DOUBLE_S(x_reduced_range_0, x_in_0, range_multiple_0);
    SUB_DOUBLE_S(x_reduced_range_1, x_in_1, range_multiple_1);
    SUB_DOUBLE_S(x_reduced_range_2, x_in_2, range_multiple_2);
    SUB_DOUBLE_S(x_reduced_range_3, x_in_3, range_multiple_3);
    SUB_DOUBLE_S(x_reduced_range_4, x_in_4, range_multiple_4);

    /* treatment of cotan ranges */
    MUL_DOUBLE_S(odd_mask_partial_0, num_ranges_away_0, half);
    MUL_DOUBLE_S(odd_mask_partial_1, num_ranges_away_1, half);
    MUL_DOUBLE_S(odd_mask_partial_2, num_ranges_away_2, half);
    MUL_DOUBLE_S(odd_mask_partial_3, num_ranges_away_3, half);
    MUL_DOUBLE_S(odd_mask_partial_4, num_ranges_away_4, half);

    FLOOR_DOUBLE_S(odd_mask_partial_floor_0, odd_mask_partial_0);
    FLOOR_DOUBLE_S(odd_mask_partial_floor_1, odd_mask_partial_1);
    FLOOR_DOUBLE_S(odd_mask_partial_floor_2, odd_mask_partial_2);
    FLOOR_DOUBLE_S(odd_mask_partial_floor_3, odd_mask_partial_3);
    FLOOR_DOUBLE_S(odd_mask_partial_floor_4, odd_mask_partial_4);

    SUB_DOUBLE_S(odd_mask_partial1_0, odd_mask_partial_0, odd_mask_partial_floor_0);
    SUB_DOUBLE_S(odd_mask_partial1_1, odd_mask_partial_1, odd_mask_partial_floor_1);
    SUB_DOUBLE_S(odd_mask_partial1_2, odd_mask_partial_2, odd_mask_partial_floor_2);
    SUB_DOUBLE_S(odd_mask_partial1_3, odd_mask_partial_3, odd_mask_partial_floor_3);
    SUB_DOUBLE_S(odd_mask_partial1_4, odd_mask_partial_4, odd_mask_partial_floor_4);

    MUL_DOUBLE_S(odd_mask_0, odd_mask_partial1_0, two);
    MUL_DOUBLE_S(odd_mask_1, odd_mask_partial1_1, two);
    MUL_DOUBLE_S(odd_mask_2, odd_mask_partial1_2, two);
    MUL_DOUBLE_S(odd_mask_3, odd_mask_partial1_3, two);
    MUL_DOUBLE_S(odd_mask_4, odd_mask_partial1_4, two);

    DOUBLE_PD_FAST(two_x1_0, x_reduced_range_0);
    DOUBLE_PD_FAST(two_x1_1, x_reduced_range_1);
    DOUBLE_PD_FAST(two_x1_2, x_reduced_range_2);
    DOUBLE_PD_FAST(two_x1_3, x_reduced_range_3);
    DOUBLE_PD_FAST(two_x1_4, x_reduced_range_4);

    SUB_DOUBLE_S(pi_2_sub_2x1_0, pi_2_hi, two_x1_0);
    SUB_DOUBLE_S(pi_2_sub_2x1_1, pi_2_hi, two_x1_1);
    SUB_DOUBLE_S(pi_2_sub_2x1_2, pi_2_hi, two_x1_2);
    SUB_DOUBLE_S(pi_2_sub_2x1_3, pi_2_hi, two_x1_3);
    SUB_DOUBLE_S(pi_2_sub_2x1_4, pi_2_hi, two_x1_4);

    MUL_DOUBLE_S(second_adjust1_0, pi_2_sub_2x1_0, odd_mask_0);
    MUL_DOUBLE_S(second_adjust1_1, pi_2_sub_2x1_1, odd_mask_1);
    MUL_DOUBLE_S(second_adjust1_2, pi_2_sub_2x1_2, odd_mask_2);
    MUL_DOUBLE_S(second_adjust1_3, pi_2_sub_2x1_3, odd_mask_3);
    MUL_DOUBLE_S(second_adjust1_4, pi_2_sub_2x1_4, odd_mask_4);

    ADD_DOUBLE_S(x_cotan_adjust_0, x_reduced_range_0, second_adjust1_0);
    ADD_DOUBLE_S(x_cotan_adjust_1, x_reduced_range_1, second_adjust1_1);
    ADD_DOUBLE_S(x_cotan_adjust_2, x_reduced_range_2, second_adjust1_2);
    ADD_DOUBLE_S(x_cotan_adjust_3, x_reduced_range_3, second_adjust1_3);
    ADD_DOUBLE_S(x_cotan_adjust_4, x_reduced_range_4, second_adjust1_4);

    /* Get Quadrant */
    MUL_DOUBLE_S(not_floored_0, x_cotan_adjust_0, one_over_pi_8);
    MUL_DOUBLE_S(not_floored_1, x_cotan_adjust_1, one_over_pi_8);
    MUL_DOUBLE_S(not_floored_2, x_cotan_adjust_2, one_over_pi_8);
    MUL_DOUBLE_S(not_floored_3, x_cotan_adjust_3, one_over_pi_8);
    MUL_DOUBLE_S(not_floored_4, x_cotan_adjust_4, one_over_pi_8);

    FLOOR_DOUBLE_S(quadrant_0, not_floored_0);
    FLOOR_DOUBLE_S(quadrant_1, not_floored_1);
    FLOOR_DOUBLE_S(quadrant_2, not_floored_2);
    FLOOR_DOUBLE_S(quadrant_3, not_floored_3);
    FLOOR_DOUBLE_S(quadrant_4, not_floored_4);

    /* Generate Mask */
    MUL_DOUBLE_S(m2_1_0, quadrant_0, half);
    MUL_DOUBLE_S(m2_1_1, quadrant_1, half);
    MUL_DOUBLE_S(m2_1_2, quadrant_2, half);
    MUL_DOUBLE_S(m2_1_3, quadrant_3, half);
    MUL_DOUBLE_S(m2_1_4, quadrant_4, half);

    FLOOR_DOUBLE_S(neg_m2_0, m2_1_0);
    FLOOR_DOUBLE_S(neg_m2_1, m2_1_1);
    FLOOR_DOUBLE_S(neg_m2_2, m2_1_2);
    FLOOR_DOUBLE_S(neg_m2_3, m2_1_3);
    FLOOR_DOUBLE_S(neg_m2_4, m2_1_4);

    SUB_DOUBLE_S(m2_0, one, neg_m2_0);
    SUB_DOUBLE_S(m2_1, one, neg_m2_1);
    SUB_DOUBLE_S(m2_2, one, neg_m2_2);
    SUB_DOUBLE_S(m2_3, one, neg_m2_3);
    SUB_DOUBLE_S(m2_4, one, neg_m2_4);

    /* Apply Mask */
    DOUBLE_PD_FAST(two_x_0, x_cotan_adjust_0);
    DOUBLE_PD_FAST(two_x_1, x_cotan_adjust_1);
    DOUBLE_PD_FAST(two_x_2, x_cotan_adjust_2);
    DOUBLE_PD_FAST(two_x_3, x_cotan_adjust_3);
    DOUBLE_PD_FAST(two_x_4, x_cotan_adjust_4);

    SUB_DOUBLE_S(pi_2_sub_2x_0, pi_2_hi, two_x_0);
    SUB_DOUBLE_S(pi_2_sub_2x_1, pi_2_hi, two_x_1);
    SUB_DOUBLE_S(pi_2_sub_2x_2, pi_2_hi, two_x_2);
    SUB_DOUBLE_S(pi_2_sub_2x_3, pi_2_hi, two_x_3);
    SUB_DOUBLE_S(pi_2_sub_2x_4, pi_2_hi, two_x_4);

    MUL_DOUBLE_S(second_adjust_0, pi_2_sub_2x_0, neg_m2_0);
    MUL_DOUBLE_S(second_adjust_1, pi_2_sub_2x_1, neg_m2_1);
    MUL_DOUBLE_S(second_adjust_2, pi_2_sub_2x_2, neg_m2_2);
    MUL_DOUBLE_S(second_adjust_3, pi_2_sub_2x_3, neg_m2_3);
    MUL_DOUBLE_S(second_adjust_4, pi_2_sub_2x_4, neg_m2_4);

    ADD_DOUBLE_S(x_second_adjust_0, x_cotan_adjust_0, second_adjust_0);
    ADD_DOUBLE_S(x_second_adjust_1, x_cotan_adjust_1, second_adjust_1);
    ADD_DOUBLE_S(x_second_adjust_2, x_cotan_adjust_2, second_adjust_2);
    ADD_DOUBLE_S(x_second_adjust_3, x_cotan_adjust_3, second_adjust_3);
    ADD_DOUBLE_S(x_second_adjust_4, x_cotan_adjust_4, second_adjust_4);

    // PRINT_FULL_M256D(x_second_adjust);
    MUL_DOUBLE_S(x_half_0, x_second_adjust_0, half);
    MUL_DOUBLE_S(x_half_1, x_second_adjust_1, half);
    MUL_DOUBLE_S(x_half_2, x_second_adjust_2, half);
    MUL_DOUBLE_S(x_half_3, x_second_adjust_3, half);
    MUL_DOUBLE_S(x_half_4, x_second_adjust_4, half);

    /* ---- Taylor Loop ---- */
    MUL_DOUBLE_S(x_square_0, x_half_0, x_half_0);
    MUL_DOUBLE_S(x_square_1, x_half_1, x_half_1);
    MUL_DOUBLE_S(x_square_2, x_half_2, x_half_2);
    MUL_DOUBLE_S(x_square_3, x_half_3, x_half_3);
    MUL_DOUBLE_S(x_square_4, x_half_4, x_half_4);

    FMADD_PD(result_q0_t1_0, taylor_coeff13, x_square_0, taylor_coeff12);
    FMADD_PD(result_q0_t2_0, result_q0_t1_0, x_square_0, taylor_coeff11);
    FMADD_PD(result_q0_t1_1, taylor_coeff13, x_square_1, taylor_coeff12);

    FMADD_PD(result_q0_t3_0, result_q0_t2_0, x_square_0, taylor_coeff10);
    FMADD_PD(result_q0_t2_1, result_q0_t1_1, x_square_1, taylor_coeff11);
    FMADD_PD(result_q0_t1_2, taylor_coeff13, x_square_2, taylor_coeff12);

    FMADD_PD(result_q0_t4_0, result_q0_t3_0, x_square_0, taylor_coeff9);
    FMADD_PD(result_q0_t3_1, result_q0_t2_1, x_square_1, taylor_coeff10);
    FMADD_PD(result_q0_t2_2, result_q0_t1_2, x_square_2, taylor_coeff11);
    FMADD_PD(result_q0_t1_3, taylor_coeff13, x_square_3, taylor_coeff12);

    FMADD_PD(result_q0_t5_0, result_q0_t4_0, x_square_0, taylor_coeff8);
    FMADD_PD(result_q0_t4_1, result_q0_t3_1, x_square_1, taylor_coeff9);
    FMADD_PD(result_q0_t3_2, result_q0_t2_2, x_square_2, taylor_coeff10);
    FMADD_PD(result_q0_t2_3, result_q0_t1_3, x_square_3, taylor_coeff11);
    FMADD_PD(result_q0_t1_4, taylor_coeff13, x_square_4, taylor_coeff12);

    FMADD_PD(result_q0_t6_0, result_q0_t5_0, x_square_0, taylor_coeff7);
    FMADD_PD(result_q0_t5_1, result_q0_t4_1, x_square_1, taylor_coeff8);
    FMADD_PD(result_q0_t4_2, result_q0_t3_2, x_square_2, taylor_coeff9);
    FMADD_PD(result_q0_t3_3, result_q0_t2_3, x_square_3, taylor_coeff10);
    FMADD_PD(result_q0_t2_4, result_q0_t1_4, x_square_4, taylor_coeff11);

    FMADD_PD(result_q0_t7_0, result_q0_t6_0, x_square_0, taylor_coeff6);
    FMADD_PD(result_q0_t6_1, result_q0_t5_1, x_square_1, taylor_coeff7);
    FMADD_PD(result_q0_t5_2, result_q0_t4_2, x_square_2, taylor_coeff8);
    FMADD_PD(result_q0_t4_3, result_q0_t3_3, x_square_3, taylor_coeff9);
    FMADD_PD(result_q0_t3_4, result_q0_t2_4, x_square_4, taylor_coeff10);

    FMADD_PD(result_q0_t8_0, result_q0_t7_0, x_square_0, taylor_coeff5);
    FMADD_PD(result_q0_t7_1, result_q0_t6_1, x_square_1, taylor_coeff6);
    FMADD_PD(result_q0_t6_2, result_q0_t5_2, x_square_2, taylor_coeff7);
    FMADD_PD(result_q0_t5_3, result_q0_t4_3, x_square_3, taylor_coeff8);
    FMADD_PD(result_q0_t4_4, result_q0_t3_4, x_square_4, taylor_coeff9);

    FMADD_PD(result_q0_t9_0, result_q0_t8_0, x_square_0, taylor_coeff4);
    FMADD_PD(result_q0_t8_1, result_q0_t7_1, x_square_1, taylor_coeff5);
    FMADD_PD(result_q0_t7_2, result_q0_t6_2, x_square_2, taylor_coeff6);
    FMADD_PD(result_q0_t6_3, result_q0_t5_3, x_square_3, taylor_coeff7);
    FMADD_PD(result_q0_t5_4, result_q0_t4_4, x_square_4, taylor_coeff8);

    FMADD_PD(result_q0_t10_0, result_q0_t9_0, x_square_0, taylor_coeff3);
    FMADD_PD(result_q0_t9_1, result_q0_t8_1, x_square_1, taylor_coeff4);
    FMADD_PD(result_q0_t8_2, result_q0_t7_2, x_square_2, taylor_coeff5);
    FMADD_PD(result_q0_t7_3, result_q0_t6_3, x_square_3, taylor_coeff6);
    FMADD_PD(result_q0_t6_4, result_q0_t5_4, x_square_4, taylor_coeff7);

    FMADD_PD(result_q0_t11_0, result_q0_t10_0, x_square_0, taylor_coeff2);
    FMADD_PD(result_q0_t10_1, result_q0_t9_1, x_square_1, taylor_coeff3);
    FMADD_PD(result_q0_t9_2, result_q0_t8_2, x_square_2, taylor_coeff4);
    FMADD_PD(result_q0_t8_3, result_q0_t7_3, x_square_3, taylor_coeff5);
    FMADD_PD(result_q0_t7_4, result_q0_t6_4, x_square_4, taylor_coeff6);

    FMADD_PD(result_q0_t12_0, result_q0_t11_0, x_square_0, taylor_coeff1);
    FMADD_PD(result_q0_t11_1, result_q0_t10_1, x_square_1, taylor_coeff2);
    FMADD_PD(result_q0_t10_2, result_q0_t9_2, x_square_2, taylor_coeff3);
    FMADD_PD(result_q0_t9_3, result_q0_t8_3, x_square_3, taylor_coeff4);
    FMADD_PD(result_q0_t8_4, result_q0_t7_4, x_square_4, taylor_coeff5);

    FMADD_PD(result_q0_t12_1, result_q0_t11_1, x_square_1, taylor_coeff1);
    FMADD_PD(result_q0_t11_2, result_q0_t10_2, x_square_2, taylor_coeff2);
    FMADD_PD(result_q0_t10_3, result_q0_t9_3, x_square_3, taylor_coeff3);
    FMADD_PD(result_q0_t9_4, result_q0_t8_4, x_square_4, taylor_coeff4);

    FMADD_PD(result_q0_t12_2, result_q0_t11_2, x_square_2, taylor_coeff1);
    FMADD_PD(result_q0_t11_3, result_q0_t10_3, x_square_3, taylor_coeff2);
    FMADD_PD(result_q0_t10_4, result_q0_t9_4, x_square_4, taylor_coeff3);

    FMADD_PD(result_q0_t12_3, result_q0_t11_3, x_square_3, taylor_coeff1);
    FMADD_PD(result_q0_t11_4, result_q0_t10_4, x_square_4, taylor_coeff2);

    FMADD_PD(result_q0_t12_4, result_q0_t11_4, x_square_4, taylor_coeff1);


    MUL_DOUBLE_S(result_q0_t13_0, result_q0_t12_0, x_square_0); // Note that here a one is Missing
    MUL_DOUBLE_S(result_q0_t13_1, result_q0_t12_1, x_square_1);
    MUL_DOUBLE_S(result_q0_t13_2, result_q0_t12_2, x_square_2);
    MUL_DOUBLE_S(result_q0_t13_3, result_q0_t12_3, x_square_3);
    MUL_DOUBLE_S(result_q0_t13_4, result_q0_t12_4, x_square_4);
                                                                         
    /* Calculate Correction */
    MUL_DOUBLE_S(result_q0_partial_0, result_q0_t13_0, x_half_0);
    MUL_DOUBLE_S(result_q0_partial_1, result_q0_t13_1, x_half_1);
    MUL_DOUBLE_S(result_q0_partial_2, result_q0_t13_2, x_half_2);
    MUL_DOUBLE_S(result_q0_partial_3, result_q0_t13_3, x_half_3);
    MUL_DOUBLE_S(result_q0_partial_4, result_q0_t13_4, x_half_4);

    FMADD_PD(correction_0, x_square_0, cor_coeff, cor_coeff);
    FMADD_PD(correction_1, x_square_1, cor_coeff, cor_coeff);
    FMADD_PD(correction_2, x_square_2, cor_coeff, cor_coeff);
    FMADD_PD(correction_3, x_square_3, cor_coeff, cor_coeff);
    FMADD_PD(correction_4, x_square_4, cor_coeff, cor_coeff);

    MUL_DOUBLE_S(eval_correction_0, neg_m2_0, correction_0);
    MUL_DOUBLE_S(eval_correction_1, neg_m2_1, correction_1);
    MUL_DOUBLE_S(eval_correction_2, neg_m2_2, correction_2);
    MUL_DOUBLE_S(eval_correction_3, neg_m2_3, correction_3);
    MUL_DOUBLE_S(eval_correction_4, neg_m2_4, correction_4);

    ADD_DOUBLE_S(result_q0_corrected_0, result_q0_partial_0, eval_correction_0);
    ADD_DOUBLE_S(result_q0_corrected_1, result_q0_partial_1, eval_correction_1);
    ADD_DOUBLE_S(result_q0_corrected_2, result_q0_partial_2, eval_correction_2);
    ADD_DOUBLE_S(result_q0_corrected_3, result_q0_partial_3, eval_correction_3);
    ADD_DOUBLE_S(result_q0_corrected_4, result_q0_partial_4, eval_correction_4);

    /* Obtain result */
    ADD_DOUBLE_S(result_q0_0, result_q0_corrected_0, x_half_0); // Here the one is added
    ADD_DOUBLE_S(result_q0_1, result_q0_corrected_1, x_half_1); // Here the one is added
    ADD_DOUBLE_S(result_q0_2, result_q0_corrected_2, x_half_2); // Here the one is added
    ADD_DOUBLE_S(result_q0_3, result_q0_corrected_3, x_half_3); // Here the one is added
    ADD_DOUBLE_S(result_q0_4, result_q0_corrected_4, x_half_4); // Here the one is added

    /* Getting Values for Double Angle */
    DOUBLE_PD_FAST(nominator_0, result_q0_0);
    DOUBLE_PD_FAST(nominator_1, result_q0_1);
    DOUBLE_PD_FAST(nominator_2, result_q0_2);
    DOUBLE_PD_FAST(nominator_3, result_q0_3);
    DOUBLE_PD_FAST(nominator_4, result_q0_4);

    MUL_DOUBLE_S(result_q0_square_0, result_q0_0, result_q0_0);
    MUL_DOUBLE_S(result_q0_square_1, result_q0_1, result_q0_1);
    MUL_DOUBLE_S(result_q0_square_2, result_q0_2, result_q0_2);
    MUL_DOUBLE_S(result_q0_square_3, result_q0_3, result_q0_3);
    MUL_DOUBLE_S(result_q0_square_4, result_q0_4, result_q0_4);

    SUB_DOUBLE_S(denominator_0, one, result_q0_square_0);
    SUB_DOUBLE_S(denominator_1, one, result_q0_square_1);
    SUB_DOUBLE_S(denominator_2, one, result_q0_square_2);
    SUB_DOUBLE_S(denominator_3, one, result_q0_square_3);
    SUB_DOUBLE_S(denominator_4, one, result_q0_square_4);

    /* Obtaining the interval results */
    DIV_DOUBLE_S(result_q01_0, nominator_0, denominator_0);
    DIV_DOUBLE_S(result_q01_1, nominator_1, denominator_1);
    DIV_DOUBLE_S(result_q01_2, nominator_2, denominator_2);
    DIV_DOUBLE_S(result_q01_3, nominator_3, denominator_3);
    DIV_DOUBLE_S(result_q01_4, nominator_4, denominator_4);

    DIV_DOUBLE_S(result_q23_0, denominator_0, nominator_0);
    DIV_DOUBLE_S(result_q23_1, denominator_1, nominator_1);
    DIV_DOUBLE_S(result_q23_2, denominator_2, nominator_2);
    DIV_DOUBLE_S(result_q23_3, denominator_3, nominator_3);
    DIV_DOUBLE_S(result_q23_4, denominator_4, nominator_4);

    /* Putting together results */
    MUL_DOUBLE_S(result_first_0, m2_0, result_q01_0);
    MUL_DOUBLE_S(result_first_1, m2_1, result_q01_1);
    MUL_DOUBLE_S(result_first_2, m2_2, result_q01_2);
    MUL_DOUBLE_S(result_first_3, m2_3, result_q01_3);
    MUL_DOUBLE_S(result_first_4, m2_4, result_q01_4);

    MUL_DOUBLE_S(result_sec_0, neg_m2_0, result_q23_0);
    MUL_DOUBLE_S(result_sec_1, neg_m2_1, result_q23_1);
    MUL_DOUBLE_S(result_sec_2, neg_m2_2, result_q23_2);
    MUL_DOUBLE_S(result_sec_3, neg_m2_3, result_q23_3);
    MUL_DOUBLE_S(result_sec_4, neg_m2_4, result_q23_4);

    ADD_DOUBLE_S(partial_result_0, result_first_0, result_sec_0);
    ADD_DOUBLE_S(partial_result_1, result_first_1, result_sec_1);
    ADD_DOUBLE_S(partial_result_2, result_first_2, result_sec_2);
    ADD_DOUBLE_S(partial_result_3, result_first_3, result_sec_3);
    ADD_DOUBLE_S(partial_result_4, result_first_4, result_sec_4);

    /* adjust cotan ranges */
    MUL_DOUBLE_S(odd_mask1_0, odd_mask_0, two);
    MUL_DOUBLE_S(odd_mask1_1, odd_mask_1, two);
    MUL_DOUBLE_S(odd_mask1_2, odd_mask_2, two);
    MUL_DOUBLE_S(odd_mask1_3, odd_mask_3, two);
    MUL_DOUBLE_S(odd_mask1_4, odd_mask_4, two);

    SUB_DOUBLE_S(odd_mask01_0, one, odd_mask1_0);
    SUB_DOUBLE_S(odd_mask01_1, one, odd_mask1_1);
    SUB_DOUBLE_S(odd_mask01_2, one, odd_mask1_2);
    SUB_DOUBLE_S(odd_mask01_3, one, odd_mask1_3);
    SUB_DOUBLE_S(odd_mask01_4, one, odd_mask1_4);

    MUL_DOUBLE_S(result_0, odd_mask01_0, partial_result_0);
    MUL_DOUBLE_S(result_1, odd_mask01_1, partial_result_1);
    MUL_DOUBLE_S(result_2, odd_mask01_2, partial_result_2);
    MUL_DOUBLE_S(result_3, odd_mask01_3, partial_result_3);
    MUL_DOUBLE_S(result_4, odd_mask01_4, partial_result_4);

    /* save result */
    SIMD_TO_DOUBLE_VEC(&res[i],                result_0);
    SIMD_TO_DOUBLE_VEC(&res[i+1*SIMD_DOUBLES], result_1);
    SIMD_TO_DOUBLE_VEC(&res[i+2*SIMD_DOUBLES], result_2);
    SIMD_TO_DOUBLE_VEC(&res[i+3*SIMD_DOUBLES], result_3);
    SIMD_TO_DOUBLE_VEC(&res[i+4*SIMD_DOUBLES], result_4);
  }

  /* Treatment of the left overs with glibc */
  int num_left_over = (n % (5*simd_doubles));

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
