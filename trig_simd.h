#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#pragma once
#include <immintrin.h>

//#define __AVX2__ 1

#define M_PI_8 (M_PI / 8)
#define RANG_REDUCTION_CORRECTION (3.8981718325193755e-17)
#define MIN_POSITIVE_COS_VALUE (0.00000000000000006123233995736766)

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

void sin_simd(double *input, double *res, size_t n);
void tan_simd(double *input, double *res, size_t n);
void safe_tan_simd(double *input, double *res, size_t n, double error_threshold);
#endif
