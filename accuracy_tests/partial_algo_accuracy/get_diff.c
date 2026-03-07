#include <stdio.h>
#include <math.h>

int main(void) {
    // "True" π approximated at long-double precision (typically 80-bit or 128-bit)
    long double pi_true = acosl(-1.0L);

    // EXACT IEEE-754 double value of π (nearest double to mathematical π)
    double pi_d = 0x1.921fb54442d18p+1;

    // Compute difference in long double, then express it as a double
    long double diff_ld = pi_true - (long double)pi_d;
    double diff_d = (double)diff_ld;

    printf("pi_d (exact double pi) = %a\n", pi_d);
    printf("diff as hex double     = %a\n", diff_d);

    return 0;
}
