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

    // ist richtig
    const SDOUBLE pi_2_lo = LOAD_DOUBLE(0x1.1a62633145c07p-54);

    const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);



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

        // Fehler Möglich erklärt aber nicht das pattern
        const SDOUBLE from_behind = SUB_DOUBLE_S(pi_2_hi, x);
        const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
        const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);

        // Masken haben keinen Fehler
        MASK8 m0 = CMP_MASK(quadrant, zero,  _CMP_EQ_OQ);
        MASK8 m1 = CMP_MASK(quadrant, one,   _CMP_EQ_OQ);
        MASK8 m2 = CMP_MASK(quadrant, two,   _CMP_EQ_OQ);
        MASK8 m3 = CMP_MASK(quadrant, three, _CMP_EQ_OQ);

        // Kann kein fehler haben da nur exponent geaendert wird
        x = MASK_MUL_PD(x, m1, x, half);

        // Die subtraction als solche is korrekt, da 
        // 
        // 1/2 a < b < 2 * a (Sterbenz)
        //
        // a = pi/2; b = x
        //
        // Pi ist als solches korrekt, weil (pi_hi + pi_lo) = pi
        x = MASK_SUB_PD(x, m2, pi_2_hi, x);
        x = MASK_ADD_PD(x, m2, x, pi_2_lo);

        // Kann kein fehler haben da nur exponent geaendert wird
        x = MASK_MUL_PD(x, m2, x, half);

        // analog zu davor bei m2
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

    int num_left_over = (n % simd_doubles);
    if (num_left_over != 0) {printf("[WARNING] Test got wrong input number");}
}
