#include <math.h>
#include <stddef.h>
#include "trig_simd.h"


void vfast_tan(double *input, double *res, size_t n) {
    double one_over_pi_8 = 1 / M_PI_8;
    double one_over_pi_2 = 1 / M_PI_2;

    for (int i = 0; i < (int) n; i ++) {
        double ranges_away = input[i] * one_over_pi_2;
        double num_ranges_away = floor(ranges_away);
        double range_multiple = num_ranges_away * M_PI_2;
        double x_reduced_range = input[i]
    }
}


void safe_vfast_tan(double *input, double *res, size_t n, double error_threshold) {
}
