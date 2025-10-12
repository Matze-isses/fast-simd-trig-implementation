#include <stdio.h>
#include <math.h>

double get_derivativ_sin(double dev_point, int order) {
    double result = 0.0;
    int derivative_option = order % 4;

    switch (derivative_option) {
        case 0: result = sin(dev_point); break;
        case 1: result = cos(dev_point); break;
        case 2: result = -sin(dev_point); break;
        case 3: result = -cos(dev_point); break;
    }

    return result;
}

unsigned long long factorial(unsigned int n) {
    if (n == 0 || n == 1) return 1;
    return n * factorial(n - 1);
}

void get_taylor_coefficiants(double *p, int n, double dev_point) {
    printf("\n");

    for (int i = 0; i < n; i++) {
        double f_prime = get_derivativ_sin(dev_point, i); 
        p[i] = f_prime / factorial(i);

        printf("Order: %d   |   Pre-Factor: %.17g\n", i, p[i]);
    }

    printf("\n");
}
