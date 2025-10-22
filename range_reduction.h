#ifndef RANGE_REDUCTION
#define RANGE_REDUCTION 1

#include <immintrin.h>

#define SIMD_LENGTH (256)

#define LOAD_DOU_VEC _mm256_loadu_pd

#define SDOUBLE __m256d
#define LOAD_DOU _mm256_set1_pd
#define MUL_DOU _mm256_mul_pd
#define SIMD_TO_DOUBLE_VEC _mm256_storeu_pd

#define FLOOR_S _mm256_floor_pd

#define MUL_DOUBLE_S _mm256_mul_pd
#define ADD_DOUBLE_S _mm256_add_pd
#define SUB_DOUBLE_S _mm256_sub_pd

#define PRINT_M256D(reg) do {                                \
    double vals[4];                                          \
    _mm256_storeu_pd(vals, (reg));                           \
    printf(#reg " = [%.6f, %.6f, %.6f, %.6f]\n",             \
           vals[0], vals[1], vals[2], vals[3]);              \
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
