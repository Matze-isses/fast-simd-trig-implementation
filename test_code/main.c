#include <stdio.h>
#include <math.h>
#include "taylor_coefficiants.h"

const double bounds[2] = {0.0, M_PI_2};

void interval_centers(double *centers, int n) {
    double length = (bounds[1] - bounds[0]) / n;
    for (int i = 0; i < n; i++) {
        centers[i] = length / 2 + i * length;
    }
}


double calculate_taylor_first(double x, double center, double *coeff, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += coeff[i] * pow(x - center, i);
    }   
    return result;
}




int main(int argc, char* argv[]) {
    int order_taylor = 8;
    double coeffs[5];

    int num_centers = 5;
    double centers[5];

    interval_centers(centers, num_centers);
    double example_part = centers[0];

    get_taylor_coefficiants(coeffs, order_taylor, example_part);

    printf("Coeff:\n");
    for (int i = 0; i < order_taylor; i++){
        printf("Order: %d; Coeff: %.17g", i, coeffs[i]);
    }


    printf("First Center: %.17g\n", example_part);
    double test_var = (centers[0] + centers[1])/2;

    double test_result = calculate_taylor_first(test_var, example_part, coeffs, order_taylor);
    double true_result = sin(test_var);

    printf("True Result: %.17g\nTest Result: %.17g\nError: %.17g\n\n", true_result, test_result, true_result-test_result);


    return 0;
}

