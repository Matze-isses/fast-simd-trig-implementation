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

  const SDOUBLE pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE pi_4 = LOAD_DOUBLE(M_PI_4);
  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  const SDOUBLE one_over_pi_2 = LOAD_DOUBLE(1/M_PI_2);

  const SDOUBLE q2_bitshift = LOAD_DOUBLE(-pow(2, -52));
  const SDOUBLE b_correction = LOAD_DOUBLE(MIN_POSITIVE_COS_VALUE);

  const SDOUBLE neg_one = LOAD_DOUBLE(-1.0);
  const SDOUBLE neg_half = LOAD_DOUBLE(-0.5);

  const SDOUBLE zero = _mm512_setzero_pd();

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

    MASK8 m_pos = CMP_MASK(x, zero, _CMP_GT_OQ);
    MASK8 m_neg = CMP_MASK(x, zero, _CMP_LT_OQ);

    SDOUBLE x_negative = _mm512_mask_mov_pd(neg_one, m_pos, one);
    SDOUBLE x_negative01 = _mm512_maskz_mov_pd(m_neg, one);

    const SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_pi_2);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, pi_2);

    x = SUB_DOUBLE_S(x, range_multiple);

    // here problems occure because negative singularities need a ceil not a floor
    const SDOUBLE singularities_away0 = ADD_DOUBLE_S(num_ranges_away, x_negative01);
    const SDOUBLE singularities_away1 = ABS_PD(singularities_away0);
    const SDOUBLE singularities_away = MUL_DOUBLE_S(singularities_away1, half);
    const SDOUBLE num_singularities_away = FLOOR_DOUBLE_S(singularities_away);

    const SDOUBLE constants_away0 = MUL_DOUBLE_S(num_singularities_away, two);
    const SDOUBLE constants_away = ADD_DOUBLE_S(constants_away0, one);

    const SDOUBLE sign_adjust_0 = MUL_DOUBLE_S(num_ranges_away, half);
    const SDOUBLE sign_adjust_1 = FLOOR_DOUBLE_S(sign_adjust_0);
    const SDOUBLE sign_adjust_2 = MUL_DOUBLE_S(sign_adjust_1, two);

    const SDOUBLE in_odd_range  = SUB_DOUBLE_S(num_ranges_away, sign_adjust_2); // 1 if in odd range, else 0
    const SDOUBLE in_even_range = SUB_DOUBLE_S(one, in_odd_range);              // 1 if in even range, else 0

    const SDOUBLE sign_adjust_3 = MUL_DOUBLE_S(in_odd_range, two);
    const SDOUBLE sign_adjust   = SUB_DOUBLE_S(one, sign_adjust_3);

    const SDOUBLE uneval_from_behind       = SUB_DOUBLE_S(pi_2, x);

    const SDOUBLE in_odd_range_reduction_1 = MUL_DOUBLE_S(uneval_from_behind, in_odd_range);
    const SDOUBLE in_odd_range_reduction   = SUB_DOUBLE_S(in_odd_range_reduction_1, x);

    x = FMADD_PD(in_odd_range_reduction, in_odd_range, x);
    const SDOUBLE from_behind = SUB_DOUBLE_S(pi_2, x);

    const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
    const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);

    /* obtaining bool vectors for the each quadrant */
    MASK8 m0 = CMP_MASK(quadrant, zero,  _CMP_EQ_OQ);
    MASK8 m1 = CMP_MASK(quadrant, one,   _CMP_EQ_OQ);
    MASK8 m2 = CMP_MASK(quadrant, two,   _CMP_EQ_OQ);
    MASK8 m3 = CMP_MASK(quadrant, three, _CMP_EQ_OQ);

    x = MASK_MUL_PD(x, m1, x, half);
    x = MASK_SUB_PD(x, m3, pi_2, x);

    x = MASK_MUL_PD(x, m2, x, half);
    x = MASK_SUB_PD(x, m2, pi_4, x);


    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);

    /* ---- Taylor Loop ---- */
    const SDOUBLE result_q0_t1 = FMADD_PD(taylor_coeff13, x_square, taylor_coeff12);
    const SDOUBLE result_q0_t2 = FMADD_PD(result_q0_t1, x_square, taylor_coeff11);
    const SDOUBLE result_q0_t3 = FMADD_PD(result_q0_t2, x_square, taylor_coeff10);
    const SDOUBLE result_q0_t4 = FMADD_PD(result_q0_t3, x_square, taylor_coeff9);
    const SDOUBLE result_q0_t5 = FMADD_PD(result_q0_t4, x_square, taylor_coeff8);
    const SDOUBLE result_q0_t6 = FMADD_PD(result_q0_t5, x_square, taylor_coeff7);
    const SDOUBLE result_q0_t7 = FMADD_PD(result_q0_t6, x_square, taylor_coeff6);
    const SDOUBLE result_q0_t8 = FMADD_PD(result_q0_t7, x_square, taylor_coeff5);
    const SDOUBLE result_q0_t9 = FMADD_PD(result_q0_t8, x_square, taylor_coeff4);
    const SDOUBLE result_q0_t10 = FMADD_PD(result_q0_t9, x_square, taylor_coeff3);
    const SDOUBLE result_q0_t11 = FMADD_PD(result_q0_t10, x_square, taylor_coeff2);
    const SDOUBLE result_q0_t12 = FMADD_PD(result_q0_t11, x_square, taylor_coeff1);

    /* ---- Correction Calculation ---- */
    const SDOUBLE one_over_from_behind = DIV_DOUBLE_S(one, from_behind);

    const SDOUBLE correction_sign_1 = SUB_DOUBLE_S(in_even_range, in_odd_range);
    const SDOUBLE correction_sign  = MUL_DOUBLE_S(correction_sign_1, x_negative);

    const SDOUBLE correction_term_1 = MUL_DOUBLE_S(b_correction, constants_away);
    const SDOUBLE correction_term_2 = MUL_DOUBLE_S(correction_term_1, correction_sign);
    const SDOUBLE correction_term = MUL_DOUBLE_S(correction_term_2, one_over_from_behind);

    const SDOUBLE coeff0 = MASK_ADD_PD(one, m3, correction_term, one);
    const SDOUBLE result_q0_1 = FMADD_PD(result_q0_t12, x_square, coeff0);

    const SDOUBLE result_q0 = MUL_DOUBLE_S(result_q0_1, x);

    /* ---- Readjusting for the second range ---- */
    const SDOUBLE nominator = MUL_DOUBLE_S(two, result_q0);
    const SDOUBLE result_q0_square = MUL_DOUBLE_S(result_q0, result_q0);
    const SDOUBLE denominator = SUB_DOUBLE_S(one, result_q0_square);

    /* Obtaining the interval results */
    const SDOUBLE result_q1 = DIV_DOUBLE_S(nominator, denominator);
    const SDOUBLE result_q2 = DIV_DOUBLE_S(denominator, nominator);
    const SDOUBLE result_q3 = DIV_DOUBLE_S(one, result_q0);

    /* Add quadrant results together */
    const SDOUBLE partial_result_0  = MASK_ADD_PD(zero, m0, result_q0, zero);
    const SDOUBLE partial_result_1  = MASK_ADD_PD(partial_result_0, m1, result_q1, partial_result_0);
    const SDOUBLE partial_result_2  = MASK_ADD_PD(partial_result_1, m2, result_q2, partial_result_1);
    const SDOUBLE partial_result_31 = MASK_ADD_PD(partial_result_2, m2, q2_bitshift, partial_result_2);
    const SDOUBLE partial_result_3  = MASK_ADD_PD(partial_result_31, m3, result_q3, partial_result_31);

    const SDOUBLE result = MUL_DOUBLE_S(sign_adjust, partial_result_3);

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
