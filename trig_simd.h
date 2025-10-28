#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#pragma once
#include <immintrin.h>

#define __AVX2__ 1
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

// Float I/O
#define SFLOAT __m256
#define LOAD_FLOAT _mm256_set1_ps
#define LOAD_FLOAT_VEC _mm256_loadu_ps
#define SIMD_TO_FLOAT_VEC _mm256_storeu_ps

// Float Operations
#define MUL_FLOAT_S _mm256_mul_ps
#define ADD_FLOAT_S _mm256_add_ps
#define SUB_FLOAT_S _mm256_sub_ps
#define FLOOR_FLOAT_S _mm256_floor_ps


#endif // __AVX2__ 

#if defined(__AVX512__) && defined(__AVX512F__)

#define SIMD_LENGTH (512)

// Double I/O
#define SDOUBLE __m512d
#define LOAD_DOUBLE _mm512_set1_pd
#define LOAD_DOUBLE_VEC _mm512_loadu_pd
#define SIMD_TO_DOUBLE_VEC _mm512_storeu_pd

// Double Operations
#define MUL_DOUBLE_S _mm512_mul_pd
#define ADD_DOUBLE_S _mm512_add_pd
#define SUB_DOUBLE_S _mm512_sub_pd

// this will cause trouble
#define FLOOR_DOUBLE_S _mm512_floor_pd

#endif

void sin_simd(double *input, double *res, size_t n, int prec);
void tan_simd(double *input, double *res, size_t n, int prec);
#endif
