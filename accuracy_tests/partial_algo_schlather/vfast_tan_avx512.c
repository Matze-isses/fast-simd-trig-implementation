#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdint.h>
#include <string.h>

#include <stdio.h>

//
#define DO_LONG_DOUBLE 1
#include "trig_simd.h"

//#define MATTHIS 1


void vfast_tan(double *input, double *res, int *lsb, size_t n) {
    int simd_doubles = SIMD_LENGTH / 64;

    const SDOUBLE pi_2_hi = SET1_DOUBLE(M_PI_2);

    // ist richtig
    const SDOUBLE pi_2_lo = SET1_DOUBLE(0x1.1a62633145c07p-54);

    const SDOUBLE one_over_pi_8 = SET1_DOUBLE(1/M_PI_8);

    // printf("%5.30Lf %5.30Lf\n", (long double) M_PI_8, (long double) one_over_pi_8);
     
    const SDOUBLE zero = SET_ZERO();
  
    const SDOUBLE half = SET1_DOUBLE(0.5);
    const SDOUBLE one = SET1_DOUBLE(1.0);
    const SDOUBLE two = SET1_DOUBLE(2.0);
    const SDOUBLE three = SET1_DOUBLE(3.0);

    const SDOUBLE taylor_coeff1  = SET1_DOUBLE(tan_tp1);
    const SDOUBLE taylor_coeff2  = SET1_DOUBLE(tan_tp2);
    const SDOUBLE taylor_coeff3  = SET1_DOUBLE(tan_tp3);
    const SDOUBLE taylor_coeff4  = SET1_DOUBLE(tan_tp4);
    const SDOUBLE taylor_coeff5  = SET1_DOUBLE(tan_tp5);
    const SDOUBLE taylor_coeff6  = SET1_DOUBLE(tan_tp6);
    const SDOUBLE taylor_coeff7  = SET1_DOUBLE(tan_tp7);
    const SDOUBLE taylor_coeff8  = SET1_DOUBLE(tan_tp8);
    const SDOUBLE taylor_coeff9  = SET1_DOUBLE(tan_tp9);
    const SDOUBLE taylor_coeff10 = SET1_DOUBLE(tan_tp10);
    const SDOUBLE taylor_coeff11 = SET1_DOUBLE(tan_tp11);
    const SDOUBLE taylor_coeff12 = SET1_DOUBLE(tan_tp12);
    const SDOUBLE taylor_coeff13 = SET1_DOUBLE(tan_tp13);
    
    for (int i = 0; i < (int) n; i += SIMD_DOUBLES) {
      SDOUBLE x   = LOAD_DOUBLE_VEC(input + i); // martin
  
        // Fehler Möglich erklärt aber nicht das pattern
        const SDOUBLE from_behind = SUB_DOUBLE_S(pi_2_hi, x);
        const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
        const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);

        // Masken haben keinen Fehler
        CMP_MASK(m0, quadrant, zero,  _CMP_EQ_OQ);
        CMP_MASK(m1, quadrant, one,   _CMP_EQ_OQ);
        CMP_MASK(m2, quadrant, two,   _CMP_EQ_OQ);
        CMP_MASK(m3, quadrant, three, _CMP_EQ_OQ);

        x = MASK_SUB_PD(x, m2, pi_2_hi, x); 
        x = MASK_MUL_PD(x, m2, x, half);

	
#if defined MATTHIS			   
	x = MASK_ADD_PD(x, m2, x, pi_2_lo);
#endif


#define F 0x1p+54
#define MF 0x1p-54

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

        // const DDOUBLE result_q0_1   = FMADD_PD(result_q0_t12, xsq_dd, one);
        // ab hier muss ddouble stehen
        SDOUBLE result_q0_1 = result_q0_t12 * x_square;
        SDOUBLE result_q0 = MUL_DOUBLE_S(result_q0_1, x);

        // (r1 + 1) * x

	if (i < 10 || i == 391)
	  printf("%d %f->%f->%1.20f %1.20f", m2, input[i], (double) x, (double) (pi_2_lo * (1.0 + 3.0 * taylor_coeff1 * x_square)), (double) result_q0);

#if !defined MATTHIS			   
        result_q0 = MASK_ADD_PD(result_q0, m2, result_q0,
                    (pi_2_lo) * (1.0 +  3.0 * taylor_coeff1 * x_square  /* hoehere Terme */ ) / 2 /* --> pi_lo + x_lo */
                ); 	
        result_q0 += x;
#endif
  
	

        SIMD_TO_DOUBLE_VEC(res + i, result_q0); // martin
	//
	if (i < 10 || i == 391)
	  printf(" %1.20f\n", (double) result_q0); 
    }

    
    //  printf("\n%3.30Lf\n%3.30Lf\n", ((long double) pi_2_hi ) * 2,
    //	   ((long double) pi_2_hi + pi_2_lo) * 2);exit(1);



	/* Treatment of the left overs with glibc */
    int num_left_over = (n % simd_doubles);
    if (num_left_over != 0) {printf("[WARNING] Test got wrong input number");}
}

/*
  
1 0.785398->0.392699->0.00000000000000007068 0.41421356237309503445 0.41421356237309508996 // long
1 0.785398->0.392699->0.00000000000000007068 0.41421356237309503445 0.41421356237309508996

  

1.108739->0.231029->0.00000000000000006450 0.23522857783480488614
1.109515->0.230641->0.00000000000000006449 0.23481946368979147999

  
  1.088579->0.241109->0.00000000000000006479 0.24589216122691867628

3.141592653589793 115997963468544
3.141592653589793 238512808959406 **
3.141592653589793 2384626433832795028841971 ****
3.141592653589793 17714696599672

*/
