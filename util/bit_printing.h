#ifndef BIT_PRINTING_H
#define BIT_PRINTING_H
#include <stdint.h> 

void print_double_bits(double value);
void print_bits_ulong(unsigned long x);
void print_bits_u8(uint8_t x);

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
