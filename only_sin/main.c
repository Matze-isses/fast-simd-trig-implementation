
#include <stdio.h>
#include <immintrin.h>
#include "trig_sin.h"

int main(void)
{
    __m512d x = _mm512_set_pd(
        0.0, 0.5, 1.0, 1.5,
        2.0, 2.5, 3.0, 3.5
    );

    __m512d y = SIN(x);

    alignas(64) double result[8];
    _mm512_store_pd(result, y);

    for (int i = 0; i < 8; ++i)
        printf("sin(%f) = %f\n", i * 0.5, result[7 - i]);

    return 0;
}
