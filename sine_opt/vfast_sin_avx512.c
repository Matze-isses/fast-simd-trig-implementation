#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include "trig_simd.h"
#include "bit_printing.h"


void vfast_sin(double *input, double *res, size_t n) {
    SET1_PD(range_reduction_correction, RANG_REDUCTION_CORRECTION);

    SET1_PD(two_pi, TWO_PI);
    SET1_PD(pi, M_PI);
    SET1_PD(three_pi, 3.0 * M_PI);
    SET1_PD(pi_2, M_PI_2);
    SET1_PD(pi_4, M_PI_4);

    SET1_PD(one_over_two_pi, ONE_OVER_TWO_PI);
    SET1_PD(one_over_pi, ONE_OVER_PI);
    SET1_PD(one_over_pi_2, ONE_OVER_PI_2);
    SET1_PD(cor_coeff, COR_COEFF);


    SET1_PD(pi_lo, 0x1.1a62633145c07p-53);

    SET1_PD(martin_const, 0.0000000000000000074010000000001);

    SET1_PD(one, 1.0);
    SET1_PD(two, 2.0);
    SET1_PD(three, 3.0);
    SET1_PD(neg_one, -1.0);

    SET1_PD(half, 0.5);

    SET1_PD(tc_sin0, sin_tp0);
    SET1_PD(tc_sin1, sin_tp1);
    SET1_PD(tc_sin2, sin_tp2);
    SET1_PD(tc_sin3, sin_tp3);
    SET1_PD(tc_sin4, sin_tp4);
    SET1_PD(tc_sin5, sin_tp5);
    SET1_PD(tc_sin6, sin_tp6);
    SET1_PD(tc_sin7, sin_tp7);
    SET1_PD(tc_sin8, sin_tp8);
    SET1_PD(tc_sin9, sin_tp9);
    SET1_PD(tc_sin10, sin_tp10);
    SET1_PD(tc_sin11, sin_tp11);

    SET1_PD(tc_cos0, cos_tp0);
    SET1_PD(tc_cos1, cos_tp1);
    SET1_PD(tc_cos2, cos_tp2);
    SET1_PD(tc_cos3, cos_tp3);
    SET1_PD(tc_cos4, cos_tp4);
    SET1_PD(tc_cos5, cos_tp5);
    SET1_PD(tc_cos6, cos_tp6);
    SET1_PD(tc_cos7, cos_tp7);
    SET1_PD(tc_cos8, cos_tp8);
    SET1_PD(tc_cos9, cos_tp9);
    SET1_PD(tc_cos10, cos_tp10);
    SET1_PD(tc_cos11, cos_tp11);

    for (int i = 0; i < (int) n; i += SIMD_DOUBLES) {
        LOAD_DOUBLE_VEC(x, &input[i]);

        /* range reduction */
        MUL_DOUBLE_S(ranges_away, x, one_over_two_pi);
        FLOOR_DOUBLE_S(num_ranges_away, ranges_away);
        MUL_DOUBLE_S(range_multiple, num_ranges_away, two_pi);

        SUB_DOUBLE_S(in_outer_range_uncorrected, x, range_multiple);
        MUL_DOUBLE_S(correction_term, x, range_reduction_correction);
        SUB_DOUBLE_S(corrected_range, in_outer_range_uncorrected, correction_term);

        /* get masks */
        MUL_DOUBLE_S(small_ranges_away, corrected_range, one_over_pi_2);
        FLOOR_DOUBLE_S(q, small_ranges_away);

        CMP_MASK(sign_mask, q, one, _CMP_GT_OQ);

        CMP_MASK(q1, q, one, _CMP_EQ_OQ);
        CMP_MASK(q2, q, two, _CMP_EQ_OQ);
        CMP_MASK(q3, q, three, _CMP_EQ_OQ);

        MASK_SUB_PD(in_range_partial, corrected_range, q2 | q3, corrected_range, pi);
        MASK_SUB_PD(in_range_hi, in_range_partial, q1 | q3, pi, in_range_partial);

        MASK_ADD_PD(in_range_partial1, in_range_hi, q1 | q3, in_range_hi, pi_lo);

        CMP_MASK(second_half, in_range_partial1, pi_4, _CMP_GT_OQ) ;
        MASK_SUB_PD(in_range, in_range_partial1, second_half, pi_2, in_range_partial1); 

        MUL_DOUBLE_S(x_square, in_range, in_range);
        MUL_DOUBLE_S(x_cube, x_square, in_range);

        FMADD_PD(sin_res_t10, tc_sin11, x_square, tc_sin10);
        FMADD_PD(sin_res_t9, sin_res_t10, x_square, tc_sin9);
        FMADD_PD(sin_res_t8, sin_res_t9, x_square, tc_sin8);
        FMADD_PD(sin_res_t7, sin_res_t8, x_square, tc_sin7);
        FMADD_PD(sin_res_t6, sin_res_t7, x_square, tc_sin6);
        FMADD_PD(sin_res_t5, sin_res_t6, x_square, tc_sin5);
        FMADD_PD(sin_res_t4, sin_res_t5, x_square, tc_sin4);
        FMADD_PD(sin_res_t3, sin_res_t4, x_square, tc_sin3);
        FMADD_PD(sin_res_t2, sin_res_t3, x_square, tc_sin2);
        FMADD_PD(sin_res_t1, sin_res_t2, x_square, tc_sin1);

        FMADD_PD(result_sin, sin_res_t1, x_cube, in_range);

        FMADD_PD(cos_res_t10, tc_cos11, x_square, tc_cos10);
        FMADD_PD(cos_res_t9, cos_res_t10, x_square, tc_cos9);
        FMADD_PD(cos_res_t8, cos_res_t9, x_square, tc_cos8);
        FMADD_PD(cos_res_t7, cos_res_t8, x_square, tc_cos7);
        FMADD_PD(cos_res_t6, cos_res_t7, x_square, tc_cos6);
        FMADD_PD(cos_res_t5, cos_res_t6, x_square, tc_cos5);
        FMADD_PD(cos_res_t4, cos_res_t5, x_square, tc_cos4);
        FMADD_PD(cos_res_t3, cos_res_t4, x_square, tc_cos3);
        FMADD_PD(cos_res_t2, cos_res_t3, x_square, tc_cos2);
        FMADD_PD(cos_res_t1, cos_res_t2, x_square, tc_cos1);
        FMADD_PD(result_cos, cos_res_t1, x_square, tc_cos0);

        __m512d result = _mm512_mask_blend_pd(second_half, result_sin, result_cos);


        MASK_MUL_PD(sign_adjusted_result, result, sign_mask, result, neg_one);

        SIMD_TO_DOUBLE_VEC(&res[i], sign_adjusted_result);
    }


}
