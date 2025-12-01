#include <stdio.h>
#include <math.h>

double __tan(double x);  // from s_tan.c

int main(void) {
    double x = 0.7;

    double y_glibc = __tan(x);
    double y_sys   = tan(x);

    printf("x = %.17g\n", x);
    printf("__tan from s_tan.c : %.17g\n", y_glibc);
    printf("system tan()       : %.17g\n", y_sys);
    printf("difference         : %.17g\n", y_glibc - y_sys);
    return 0;
}
