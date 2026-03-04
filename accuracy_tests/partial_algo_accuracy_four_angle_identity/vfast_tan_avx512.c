#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdint.h>
#include <string.h>

#include <stdio.h>
#include "trig_simd.h"



void vfast_tan(double *input, double *res, int *lsb, size_t n) {
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
    const SDOUBLE four = LOAD_DOUBLE(4.0);
    const SDOUBLE six = LOAD_DOUBLE(6.0);

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

        x = MUL_DOUBLE_S(x, zero_two_five);

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

        const SDOUBLE tan1 = MUL_DOUBLE_S(result_q0_1, x);
        const SDOUBLE tan2 = MUL_DOUBLE_S(tan1, tan1);
        const SDOUBLE tan3 = MUL_DOUBLE_S(tan2, tan1);
        const SDOUBLE tan4 = MUL_DOUBLE_S(tan2, tan2);
        
        SDOUBLE nominator0 = MUL_DOUBLE_S(four, tan1);
        SDOUBLE nominator1 = MUL_DOUBLE_S(four, tan3);
        SDOUBLE nominator = SUB_DOUBLE_S(nominator0, nominator1);

        SDOUBLE denominator1 = MUL_DOUBLE_S(six, tan2);
        SDOUBLE denominator = SUB_DOUBLE_S(one, denominator1);
        denominator = ADD_DOUBLE_S(denominator, tan4);

        SDOUBLE result = DIV_DOUBLE_S(nominator, denominator);

        SIMD_TO_DOUBLE_VEC(&res[i], result);
    }

    /* Treatment of the left overs with glibc */
    int num_left_over = (n % simd_doubles);
    if (num_left_over != 0) {printf("[WARNING] Test got wrong input number");}
}

