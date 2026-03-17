#include <stdio.h>
#include <math.h>

static void split_multiple_of_pi(long double factor)
{
    long double x = factor * acosl(-1.0L);
    double a = (double)x;
    long double rest = x - (long double)a;
    double b = (double)rest;

    printf("factor           = %.21Lg\n", factor);
    printf("x (long double)  = %La\n", x);
    printf("a (double)       = %a\n", a);
    printf("b (double)       = %a   %25.17e\n", b, b);
    printf("a+b in long double = %La\n", (long double)a + (long double)b);
}

int main(void)
{
    split_multiple_of_pi(1.0L);
    return 0;
}
