#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#pragma once
#include <immintrin.h>

#define __AVX2__ 1
#define M_PI_8 (M_PI / 8)
#define TAN_CORRECTION (0.00000000000000006123233995736765)

#if defined(__AVX2__)

#define SIMD_LENGTH (256)

// Double I/O
#define SDOUBLE __m256d
#define LOAD_DOUBLE _mm256_set1_pd
#define LOAD_DOUBLE_VEC _mm256_loadu_pd
#define SIMD_TO_DOUBLE_VEC _mm256_storeu_pd

// Double Operations
#define MUL_DOUBLE_S _mm256_mul_pd
#define DIV_DOUBLE_S _mm256_div_pd

#define ADD_DOUBLE_S _mm256_add_pd
#define SUB_DOUBLE_S _mm256_sub_pd

#define FLOOR_DOUBLE_S _mm256_floor_pd

// for (a, b, c) does return (a * b) + c
#define FMADD_PD _mm256_fmadd_pd 
#define CEIL_PD _mm256_ceil_pd

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

#endif // __AVX2__ 

#if defined(__AVX512__) && defined(__AVX512F__)
#define SIMD_LENGTH (512)

// Double I/O
#define SDOUBLE __m512d
#define LOAD_DOUBLE         _mm512_set1_pd
#define LOAD_DOUBLE_VEC     _mm512_loadu_pd
#define SIMD_TO_DOUBLE_VEC  _mm512_storeu_pd

// Double Operations
#define MUL_DOUBLE_S _mm512_mul_pd
#define DIV_DOUBLE_S _mm512_div_pd

#define ADD_DOUBLE_S _mm512_add_pd
#define SUB_DOUBLE_S _mm512_sub_pd

// floor(a)
#define FLOOR_DOUBLE_S(a) \
    _mm512_roundscale_pd((a), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)

// for (a, b, c) returns (a * b) + c
#define FMADD_PD _mm512_fmadd_pd

// ceil(a)
#define CEIL_PD(a) \
    _mm512_roundscale_pd((a), _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC)

// |a|
#define ABS_PD(a) \
    _mm512_max_pd((a), _mm512_mul_pd(_mm512_set1_pd(-1.0), (a)))

// Replace +/-inf (by absolute value) with 0.0
#define REMOVE_INF(a) \
    _mm512_mask_mov_pd( \
        (a), \
        _mm512_cmp_pd_mask( \
            _mm512_andnot_pd(_mm512_set1_pd(-0.0), (a)), /* clear sign bit */ \
            _mm512_set1_pd(INFINITY), \
            _CMP_EQ_OQ \
        ), \
        _mm512_setzero_pd() \
    )

#endif

void sin_simd(double *input, double *res, size_t n);
void tan_simd(double *input, double *res, size_t n);
#endif
