
#include <stdio.h>
#include <math.h>

int main() {

    const int N = 16;
    double x[N];

    double pi = M_PI;
    double start = pi / 4.0;
    double end = pi / 2.0;

    /* generate 16 x[i] between pi/4 and pi/2 */
    for (int i = 0; i < N; i++) {
        x[i] = start + (end - start) * i / (N - 1);
    }

    double sum = 0.0;
    double c = 0.0;   // Kahan correction

    for (int i = 0; i < N; i++) {

        double value = M_PI_2 - x[i];
        double y = x[i] + value - M_PI_2;
        printf("Kahan sum = %.17g\n", y);
    }
    return 0;
}
