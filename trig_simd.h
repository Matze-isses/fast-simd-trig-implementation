#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#pragma once
#include <immintrin.h>

//#define __AVX2__ 1

#define M_PI_8 (M_PI / 8)
#define RANG_REDUCTION_CORRECTION (3.8981718325193755e-17)
#define MIN_POSITIVE_COS_VALUE (0.00000000000000006123233995736766)


// Taylor Polynomaial Coefficiants
#define tan_tp1 (0x1.5555555555555p-2)
#define tan_tp2 (0x1.1111111111111p-3)
#define tan_tp3 (0x1.ba1ba1ba1ba1cp-5)
#define tan_tp4 (0x1.664f4882c10fap-6)
#define tan_tp5 (0x1.226e355e6c23dp-7)
#define tan_tp6 (0x1.d6d3d0e157dep-9)
#define tan_tp7 (0x1.7da36452b75e3p-10)
#define tan_tp8 (0x1.3558248036744p-11)
#define tan_tp9 (0x1.f57d7734d1664p-13)
#define tan_tp10 (0x1.967e18afcafadp-14)
#define tan_tp11 (0x1.497d8eea25259p-15)
#define tan_tp12 (0x1.0b132d39a605p-16)
#define tan_tp13 (0x1.b0f72d3ee24e9p-18)


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

#define __AVX512F__ 1

#if defined(__AVX512F__)

#define USE_AVX512 (true)
#define SIMD_LENGTH (512)
#define SIMD_DOUBLES (8)

// Double I/O
#define SDOUBLE __m512d
#define LOAD_DOUBLE _mm512_set1_pd
#define LOAD_DOUBLE_VEC _mm512_loadu_pd
#define SIMD_TO_DOUBLE_VEC _mm512_storeu_pd

#define SET_ZERO _mm512_setzero_pd

#define MASK8 __mmask8
#define CMP_MASK _mm512_cmp_pd_mask
#define MASK_ADD_PD _mm512_mask_add_pd
#define MASK_SUB_PD _mm512_mask_sub_pd
#define MASK_MUL_PD _mm512_mask_mul_pd

#define MASK_MOV_PD(dst, mask, src_false, src_true) \
    SDOUBLE dst = _mm512_mask_mov_pd((src_false), (mask), (src_true))

#define MASKZ_MOV_PD(dst, mask, vec) \
    const SDOUBLE (dst) = _mm512_maskz_mov_pd((mask), (vec))

#define GEN_MASK_IF_ODD(dst, vec)                    \
    MASK8 dst = _mm512_test_epi64_mask(              \
        _mm512_cvttpd_epi64((vec)),                   \
        _mm512_set1_epi64(1)                          \
    )

#define FLIP_SIGN_IF_MASK_PD(dst, mask, vec)         \
    SDOUBLE dst = _mm512_castsi512_pd(               \
        _mm512_mask_xor_epi64(                        \
            _mm512_castpd_si512((vec)),               \
            (mask),                                   \
            _mm512_castpd_si512((vec)),               \
            _mm512_set1_epi64(0x8000000000000000ULL)  \
        )                                             \
    )

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

#elif defined(__AVX2__)

#define USE_AVX512 (false)
#define SIMD_LENGTH (256)
#define SIMD_DOUBLES (4)

// Double I/O
#define SDOUBLE __m256d
#define LOAD_DOUBLE _mm256_set1_pd
#define LOAD_DOUBLE_VEC _mm256_loadu_pd
#define SIMD_TO_DOUBLE_VEC _mm256_storeu_pd

#define SET_ZERO _mm256_setzero_pd

// Double Operations
#define MUL_DOUBLE_S _mm256_mul_pd
#define DIV_DOUBLE_S _mm256_div_pd

#define ADD_DOUBLE_S _mm256_add_pd
#define SUB_DOUBLE_S _mm256_sub_pd

#define FLOOR_DOUBLE_S _mm256_floor_pd


// for (a, b, c) does return (a * b) + c
#define FMADD_PD _mm256_fmadd_pd 
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

#endif

void vfast_sin(double *input, double *res, size_t n);
void vfast_tan(double *input, double *res, size_t n);
void safe_vfast_tan(double *input, double *res, size_t n, double error_threshold);
#endif
