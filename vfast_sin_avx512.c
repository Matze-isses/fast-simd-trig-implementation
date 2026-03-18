#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include "trig_simd.h"
#include "util/bit_printing.h"



void vfast_sin(double *input, double *res, size_t n) {
    SET1_PD(range_reduction_correction, RANG_REDUCTION_CORRECTION);

    SET1_PD(two_pi, TWO_PI);
    SET1_PD(pi, M_PI);
    SET1_PD(three_pi, 3.0 * M_PI);
    SET1_PD(pi_2, M_PI_2);

    SET1_PD(one_over_two_pi, ONE_OVER_TWO_PI);
    SET1_PD(one_over_pi, ONE_OVER_PI);
    SET1_PD(one_over_pi_2, ONE_OVER_PI_2);
    SET1_PD(cor_coeff, COR_COEFF);

    SET1_PD(pi_lo, 0x1.1a8p-53);

    SET1_PD(one, 1.0);
    SET1_PD(two, 2.0);
    SET1_PD(three, 3.0);
    SET1_PD(neg_one, -1.0);

    SET1_PD(taylor_coeff0, sin_tp0);
    SET1_PD(taylor_coeff1, sin_tp1);
    SET1_PD(taylor_coeff2, sin_tp2);
    SET1_PD(taylor_coeff3, sin_tp3);
    SET1_PD(taylor_coeff4, sin_tp4);
    SET1_PD(taylor_coeff5, sin_tp5);
    SET1_PD(taylor_coeff6, sin_tp6);
    SET1_PD(taylor_coeff7, sin_tp7);
    SET1_PD(taylor_coeff8, sin_tp8);
    SET1_PD(taylor_coeff9, sin_tp9);
    SET1_PD(taylor_coeff10, sin_tp10);
    SET1_PD(taylor_coeff11, sin_tp11);
    SET1_PD(taylor_coeff12, sin_tp12);
    SET1_PD(taylor_coeff13, sin_tp13);
    SET1_PD(taylor_coeff14, sin_tp14);
    SET1_PD(taylor_coeff15, sin_tp15);
    SET1_PD(taylor_coeff16, sin_tp16);
    SET1_PD(taylor_coeff17, sin_tp17);
    SET1_PD(taylor_coeff18, sin_tp18);

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
        MASK_ADD_PD(in_range_t, in_range_hi, q1 | q3, pi_lo, in_range_hi);
        MASK_ADD_PD(in_range, in_range_t, q3, in_range_t, pi_lo);

        MUL_DOUBLE_S(x_square, in_range, in_range);

        FMADD_PD(result_q0_t1, taylor_coeff18, x_square, taylor_coeff17);
        FMADD_PD(result_q0_t2, result_q0_t1, x_square, taylor_coeff16);
        FMADD_PD(result_q0_t3, result_q0_t2, x_square, taylor_coeff15);
        FMADD_PD(result_q0_t4, result_q0_t3, x_square, taylor_coeff14);
        FMADD_PD(result_q0_t5, result_q0_t4, x_square, taylor_coeff13);
        FMADD_PD(result_q0_t6, result_q0_t5, x_square, taylor_coeff12);
        FMADD_PD(result_q0_t7, result_q0_t6, x_square, taylor_coeff11);
        FMADD_PD(result_q0_t8, result_q0_t7, x_square, taylor_coeff10);
        FMADD_PD(result_q0_t9, result_q0_t8, x_square, taylor_coeff9);
        FMADD_PD(result_q0_t10, result_q0_t9, x_square, taylor_coeff8);
        FMADD_PD(result_q0_t11, result_q0_t10, x_square, taylor_coeff7);
        FMADD_PD(result_q0_t12, result_q0_t11, x_square, taylor_coeff6);
        FMADD_PD(result_q0_t13, result_q0_t12, x_square, taylor_coeff5);
        FMADD_PD(result_q0_t14, result_q0_t13, x_square, taylor_coeff4);
        FMADD_PD(result_q0_t15, result_q0_t14, x_square, taylor_coeff3);
        FMADD_PD(result_q0_t16, result_q0_t15, x_square, taylor_coeff2);
        FMADD_PD(result_q0_t17, result_q0_t16, x_square, taylor_coeff1);
        FMADD_PD(result_q0_t18, result_q0_t17, x_square, taylor_coeff0);
        MUL_DOUBLE_S(result, result_q0_t18, in_range);

        MASK_MUL_PD(sign_adjusted_result, result, sign_mask, result, neg_one);

        SIMD_TO_DOUBLE_VEC(&res[i], sign_adjusted_result);
    }

    int num_left_over = (n % 4);

    for (int i = n - num_left_over; i < (int)n; i++) {
        res[i] = sin(input[i]);
    }
}
