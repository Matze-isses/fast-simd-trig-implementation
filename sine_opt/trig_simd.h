#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#pragma once
#include <immintrin.h>

#define M_PI_8 (M_PI / 8)
#define RANG_REDUCTION_CORRECTION (3.8981718325193755e-17)
#define PI_LO (0x1.1a62633145c07p-54)
#define COR_COEFF (0x1.1a62633145c07p-55)

#define TWO_PI (0x1.921fb54442d18p+2)
#define ONE_OVER_TWO_PI (0x1.45f306dc9c883p-3)
#define ONE_OVER_PI (0x1.45f306dc9c883p-2)
#define ONE_OVER_PI_2 (0x1.45f306dc9c883p-1)



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

#define cos_tp0 (1.0)
#define cos_tp1 (-0.5)
#define cos_tp2 (0.041666666666666664)
#define cos_tp3 (-0.001388888888888889)
#define cos_tp4 (2.48015873015873e-05)
#define cos_tp5 (-2.755731922398589e-07)
#define cos_tp6 (2.08767569878681e-09)
#define cos_tp7 (-1.1470745597729725e-11)
#define cos_tp8 (4.779477332387385e-14)
#define cos_tp9 (-1.5619206968586225e-16)
#define cos_tp10 (4.110317623312165e-19)
#define cos_tp11 (-8.896791392450574e-22)
#define cos_tp12 (1.6117375710961184e-24)
#define cos_tp13 (-2.4795962632247976e-27)
#define cos_tp14 (3.279889237069838e-30)

#define __AVX512F__ 1


#if defined(__AVX512F__)

#define USE_AVX512 (true)
#define SIMD_LENGTH (512)
#define SIMD_DOUBLES (8)

// Double I/O
#define MASK8 __mmask8
#define SDOUBLE __m512d

#define SET1_PD(dst, a)\
    const SDOUBLE (dst) = _mm512_set1_pd(a)

#define SIMD_TO_DOUBLE_VEC _mm512_storeu_pd

#define LOAD_DOUBLE_VEC(dst, src)\ 
    const SDOUBLE (dst) = _mm512_loadu_pd(src)


#define ADD_DOUBLE_S(dst, a, b)\
    const SDOUBLE (dst) = _mm512_add_pd(a, b)

#define SUB_DOUBLE_S(dst, a, b)\
    const SDOUBLE (dst) = _mm512_sub_pd(a, b)

// Double Operations
#define MUL_DOUBLE_S(dst, a, b)\
    const SDOUBLE (dst) = _mm512_mul_pd(a, b)

#define DIV_DOUBLE_S(dst, a, b)\
    const SDOUBLE (dst) = _mm512_div_pd(a, b)


#define FMADD_PD(dst, a, b, c)\
    const SDOUBLE (dst) = _mm512_fmadd_pd(a, b, c)


#define CMP_MASK(dst, vec, a, qualifier) \
    const MASK8 (dst) = _mm512_cmp_pd_mask(vec, a, qualifier)

#define MASK_ADD_PD(dst, src, mask, a, b) \
    const SDOUBLE (dst) = _mm512_mask_add_pd(src, mask, a, b)

#define MASK_SUB_PD(dst, src, mask, a, b) \
    const SDOUBLE (dst) = _mm512_mask_sub_pd(src, mask, a, b)

#define MASK_MUL_PD(dst, src, mask, a, b) \
    const SDOUBLE (dst) = _mm512_mask_mul_pd(src, mask, a, b)

#define MASK_MOV_PD(dst, mask, src_false, src_true) \
    const SDOUBLE (dst) = _mm512_mask_mov_pd((src_false), (mask), (src_true))

#define MASKZ_MOV_PD(dst, mask, vec) \
    const SDOUBLE (dst) = _mm512_maskz_mov_pd((mask), (vec))

#define GEN_MASK_IF_ODD(dst, vec)                    \
    const MASK8 dst = _mm512_test_epi64_mask(              \
        _mm512_cvttpd_epi64((vec)),                   \
        _mm512_set1_epi64(1)                          \
    )

#define FLIP_SIGN_IF_MASK_PD(dst, mask, vec)         \
    const SDOUBLE (dst) = _mm512_castsi512_pd(               \
        _mm512_mask_xor_epi64(                        \
            _mm512_castpd_si512((vec)),               \
            (mask),                                   \
            _mm512_castpd_si512((vec)),               \
            _mm512_set1_epi64(0x8000000000000000ULL)  \
        )                                             \
    )

#define HALF_PD_FAST(dst, vec)                                                \
    const SDOUBLE (dst) = _mm512_castsi512_pd(                                        \
        _mm512_sub_epi64(                                                      \
            _mm512_castpd_si512((vec)),                                        \
            _mm512_set1_epi64(0x0010000000000000ULL)                           \
        )                                                                      \
    )

#define DOUBLE_PD_FAST(dst, vec)                                              \
    const SDOUBLE dst = _mm512_castsi512_pd(                                        \
        _mm512_add_epi64(                                                      \
            _mm512_castpd_si512((vec)),                                        \
            _mm512_set1_epi64(0x0010000000000000ULL)                           \
        )                                                                      \
    )


// floor(a)
#define FLOOR_DOUBLE_S(dst, a)\
    const SDOUBLE (dst) = _mm512_roundscale_pd((a), _MM_FROUND_TO_NEG_INF)

#define ABS_PD(dst, a) \
    const SDOUBLE (dst) = _mm512_and_pd((a), _mm512_castsi512_pd(_mm512_set1_epi64(0x7FFFFFFFFFFFFFFF)))

#elif defined(__AVX2__)

#define USE_AVX512 (false)
#define SIMD_LENGTH (256)
#define SIMD_DOUBLES (4)

// Double I/O
#define MASK8   uint8_t
#define SDOUBLE __m256d
#define SINT    __m256i

// ---------- macros ----------
#define LOAD_DOUBLE_VEC(dst, src) \
    const SDOUBLE (dst) = _mm256_loadu_pd(src)

#define SIMD_TO_DOUBLE_VEC _mm256_storeu_pd

#define SET1_PD(dst, a) \
    const SDOUBLE (dst) = _mm256_set1_pd(a)

#define ADD_DOUBLE_S(dst, a, b) \
    const SDOUBLE (dst) = _mm256_add_pd(a, b)

#define SUB_DOUBLE_S(dst, a, b) \
    const SDOUBLE (dst) = _mm256_sub_pd(a, b)

// Double Operations
#define MUL_DOUBLE_S(dst, a, b) \
    const SDOUBLE (dst) = _mm256_mul_pd(a, b)

#define DIV_DOUBLE_S(dst, a, b) \
    const SDOUBLE (dst) = _mm256_div_pd(a, b)

// requires FMA support (-mfma)
#define FMADD_PD(dst, a, b, c) \
    const SDOUBLE (dst) = _mm256_fmadd_pd(a, b, c)

#define HALF_PD_FAST(dst, vec) \
    const SDOUBLE (dst) = _mm256_castsi256_pd( \
        _mm256_sub_epi64( \
            _mm256_castpd_si256((vec)), \
            _mm256_set1_epi64x(0x0010000000000000ULL) \
        ) \
    )

#define DOUBLE_PD_FAST(dst, vec) \
    const SDOUBLE (dst) = _mm256_castsi256_pd( \
        _mm256_add_epi64( \
            _mm256_castpd_si256((vec)), \
            _mm256_set1_epi64x(0x0010000000000000ULL) \
        ) \
    )

// floor(a)
#define FLOOR_DOUBLE_S(dst, a) \
    const SDOUBLE (dst) = _mm256_floor_pd((a))


// abs(a)
#define ABS_PD(dst, a) \
    const SDOUBLE (dst) = _mm256_and_pd((a), _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF)))

#else

#define SIMD_DOUBLES (1)

#endif

void vfast_sin(double *input, double *res, size_t n);
void vfast_tan(double *input, double *res, size_t n);
void safe_vfast_tan(double *input, double *res, size_t n, double error_threshold);
#endif
