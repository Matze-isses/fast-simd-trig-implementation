#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#pragma once
#include <immintrin.h>

//#define __AVX2__ 1

#define M_PI_8 (0x1.921fb54442d18p-2)
#define RANG_REDUCTION_CORRECTION (3.8981718325193755e-17)
#define MIN_POSITIVE_COS_VALUE (0x1.1a62633145c07p-54)


// Taylor Polynomaial Coefficiants
#define tan_tp1 (0.3333333333333333)
#define tan_tp2 (0.13333333333333333)
//#define tan_tp1 0.33333333333333333333333333333333333333333333333333333333
//#define tan_tp2 0.133333333333333333333333333333333333333333333333333333333
#define tan_tp3 (0.05396825396825397)
#define tan_tp4 (0.021869488536155203)
#define tan_tp5 (0.008863235529902197)
#define tan_tp6 (0.003592128036572481)
#define tan_tp7 (0.0014558343870513183)
#define tan_tp8 (0.000590027440945586)
#define tan_tp9 (0.00023912911424355248)
#define tan_tp10 (9.691537956929451e-05)
#define tan_tp11 (3.927832388331683e-05)
#define tan_tp12 (1.5918905069328964e-05)
#define tan_tp13 (6.451689215655431e-06)


#define sin_tp0 (1.0)
#define sin_tp1 (-0.16666666666666666)
#define sin_tp2 (0.0083333333333333332)
#define sin_tp3 (-0.00019841269841269841)
#define sin_tp4 (2.7557319223985893e-06)
#define sin_tp5 (-2.505210838544172e-08)
#define sin_tp6 (1.6059043836821613e-10)
#define sin_tp7 (-7.6471637318198164e-13)
#define sin_tp8 (2.8114572543455206e-15)
#define sin_tp9 (-8.2206352466243295e-18)
#define sin_tp10 (1.9572941063391263e-20)
#define sin_tp11 (-3.8681701706306835e-23)
#define sin_tp12 (6.4469502843844736e-26)
#define sin_tp13 (-9.183689863795546e-29)
#define sin_tp14 (1.1309962886447718e-31)
#define sin_tp15 (-1.2161250415535181e-34)
#define sin_tp16 (1.1516335620771951e-37)
#define sin_tp17 (-9.6775929586318907e-41)
#define sin_tp18 (7.2654601791530714e-44)
#define sin_tp19 (-4.9024697565135435e-47)
#define sin_tp20 (2.9893108271424046e-50)


#define DDOUBLE SDOUBLE

#if defined(__AVX512F__)

#define USE_AVX512 (true)
#define SIMD_LENGTH (512)
#define SIMD_DOUBLES (8)

// Double I/O
#define SDOUBLE __m512d
#define SET1_DOUBLE _mm512_set1_pd
#define LOAD_DOUBLE_VEC _mm512_loadu_pd
#define SIMD_TO_DOUBLE_VEC _mm512_storeu_pd

#define SET_ZERO _mm512_setzero_pd

#define MASK8 __mmask8
#define CMP_MASK(m,A,B,C) MASK8 m = _mm512_cmp_pd_mask(A,B,C)
#define MASK_ADD_PD _mm512_mask_add_pd
#define MASK_SUB_PD _mm512_mask_sub_pd
#define MASK_MUL_PD _mm512_mask_mul_pd

// Double Operations
#define MUL_DOUBLE_S _mm512_mul_pd
#define DIV_DOUBLE_S _mm512_div_pd

#define ADD_DOUBLE_S _mm512_add_pd
#define SUB_DOUBLE_S _mm512_sub_pd

// floor(a)
#define FLOOR_DOUBLE_S(a) _mm512_roundscale_pd((a), _MM_FROUND_TO_NEG_INF)

#define FMADD_PD _mm512_fmadd_pd
#define CMP_PD _mm512_cmp_pd_mask

#define ABS_PD(a) _mm512_abs_pd((a))

// Replace +/-inf (by absolute value) with 0.0
#define REMOVE_INF(a) \
    _mm512_mask_mov_pd( \
        (a), \
        _mm512_cmp_pd_mask( \
            _mm512_andnot_pd(_mm512_set1_pd(-0.0), (a)), \
            _mm512_set1_pd(INFINITY), \
            _CMP_EQ_OQ \
        ), \
        _mm512_setzero_pd() \
    )

#define PRINT_FULL_M512D(reg) do {                            \
    double vals[8];                                           \
    _mm512_storeu_pd(vals, (reg));                            \
    printf(#reg " = [%.19g, %.19g, %.19g, %.19g, %.19g, %.19g, %.19g, %.19g]\n", \
           vals[0], vals[1], vals[2], vals[3],                \
           vals[4], vals[5], vals[6], vals[7]);               \
} while (0)

#define SET1(z, a)


#else

#define USE_AVX512 (false)

#define MASK_ADD_PD(src, k, a, b) ((one - (k)) * (src) + (k) * ((a) + (b)))
#define MASK_SUB_PD(src, k, a, b) ((one - (k)) * (src) + (k) * ((a) - (b)))
#define MASK_MUL_PD(src, k, a, b) ((one - (k)) * (src) + (k) * ((a) * (b)))


// Double Operations
#define MUL_DOUBLE_S(a,b) ((a) * (b))
#define DIV_DOUBLE_S(a,b) ((a) / (b))

#define ADD_DOUBLE_S(a,b) ((a) + (b))
#define SUB_DOUBLE_S(a,b) ((a) - (b))

#define FMADD_PD(a,b,c) ((a) * (b) + (c)) // _mm256_fmadd_pd 

#define MASK8 SDOUBLE

#if defined(__AVX2__)

#define SIMD_LENGTH (256)
#define SIMD_DOUBLES (4)

// Double I/O
#define SDOUBLE __m256d
#define SET1_DOUBLE _mm256_set1_pd
#define LOAD_DOUBLE_VEC _mm256_loadu_pd
#define SIMD_TO_DOUBLE_VEC _mm256_storeu_pd

#define SET_ZERO _mm256_setzero_pd
#define SET1(z, a) const SDOUBLE z = _mm256_set1_pd(a)
#define CMP_MASK(m,A,B,C) MASK8 m = _mm256_cmp_pd(A,B,C); m = _mm256_and_pd(m,one);

// Double Operations

#define FLOOR_DOUBLE_S _mm256_floor_pd



// for (a, b, c) does return (a * b) + c
#define FMADD_PD(a,b,c) ((a) * (b) + (c)) // _mm256_fmadd_pd 
#define CEIL_PD _mm256_ceil_pd
#define CMP_PD _mm256_cmp_pd
#define BLEND_PD _mm256_blendv_pd


#define ABS_PD(a) \
    _mm256_max_pd((a), _mm256_mul_pd(_mm256_set1_pd(-1.0), (a)))

#define REMOVE_INF(a) \
    _mm256_blendv_pd( \
        (a), \
        _mm256_setzero_pd(), \
        _mm256_cmp_pd( \
            _mm256_andnot_pd(_mm256_set1_pd(-0.0), (a)), \
            _mm256_set1_pd(INFINITY), \
            _CMP_EQ_OQ \
        ) \
    )

#else


#define SDOUBLE double
#define SIMD_LENGTH (64)
#define SIMD_DOUBLES (1)

// Double I/O
#define SET1_DOUBLE 
#define LOAD_DOUBLE_VEC(z) (z)[0]
#define SIMD_TO_DOUBLE_VEC(a,b) (a)[0] = b

#define SET_ZERO() 0
#define SET1(z, a) const SDOUBLE z = (a)
#define CMP_MASK(m,A,B,C) bool m = (C) == _CMP_EQ_OQ ? (A) == (B) : true;



#if defined DO_LONG_DOUBLE
  #undef DDOUBLE
  #define DDOUBLE long double
  #define FLOOR_DOUBLE_S floorl
  #define TRUNC truncl
  #define CEIL_PD _ceill
  #warning long double
#else
  #define FLOOR_DOUBLE_S floor
  #define CEIL_PD _ceil
  #define TRUNC trunc
#endif


#endif // else AVX2
	 

#endif // else AVX512


#endif // end header
