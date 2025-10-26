#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#include <immintrin.h>

#define SIMD_LENGTH (256)


// Double I/O
#define SDOUBLE __m256d
#define LOAD_DOUBLE _mm256_set1_pd
#define LOAD_DOUBLE_VEC _mm256_loadu_pd
#define SIMD_TO_DOUBLE_VEC _mm256_storeu_pd

// Double Operations
#define MUL_DOUBLE_S _mm256_mul_pd
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




#define PRINT_M256D(reg) do {                                \
    double vals[4];                                          \
    _mm256_storeu_pd(vals, (reg));                           \
    printf(#reg " = [%.6f, %.6f, %.6f, %.6f]\n",             \
           vals[0], vals[1], vals[2], vals[3]);              \
} while (0)

#define PRINT_M256(reg) do {                                   \
    float vals[8];                                             \
    _mm256_storeu_ps(vals, (reg));                             \
    printf(#reg " = [%.6f, %.6f, %.6f, %.6f, "                 \
                  "%.6f, %.6f, %.6f, %.6f]\n",                 \
           vals[0], vals[1], vals[2], vals[3],                 \
           vals[4], vals[5], vals[6], vals[7]);               \
} while (0)

#define PRINT_ARRAY(arr, len) do {                                  \
    printf(#arr " = [");                                            \
    for (size_t _i = 0; _i < (len); _i++) {                          \
        printf("%.6f", (double)(arr)[_i]);                           \
        if (_i + 1 < (len)) printf(", ");                            \
    }                                                                \
    printf("]\n");                                                   \
} while (0)

#endif
