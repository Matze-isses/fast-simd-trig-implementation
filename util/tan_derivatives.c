#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Computes the n-th derivative of tan(x) at a given x,
 * using the derivative polynomial recurrence:
 *   P_0(t) = t
 *   P_{k+1}(t) = (1 + t^2) * d/dt P_k(t)
 * where t = tan(x) and d^n/dx^n tan(x) = P_n(tan(x)).
 */
double tan_derivative(int n, double x) {
    if (n < 0) {
        return NAN;  // undefined
    }

    double t = tan(x);

    // n = 0: just tan(x)
    if (n == 0) {
        return t;
    }

    int max_deg = n + 1;  // degree of P_n is n+1
    // Allocate coefficient arrays for polynomials:
    // a = current P_k, dp = its derivative, b = next P_{k+1}
    double *a  = calloc(max_deg + 1, sizeof(double));
    double *dp = calloc(max_deg + 1, sizeof(double));
    double *b  = calloc(max_deg + 1, sizeof(double));

    if (!a || !dp || !b) {
        // allocation failed
        free(a); free(dp); free(b);
        return NAN;
    }

    // P_0(t) = t  --> a[1] = 1
    int deg = 1;
    a[0] = 0.0;
    a[1] = 1.0;

    // Build up to P_n
    for (int k = 0; k < n; ++k) {
        // Clear dp and b
        for (int i = 0; i <= max_deg; ++i) {
            dp[i] = 0.0;
            b[i]  = 0.0;
        }

        // dp(t) = d/dt P_k(t)
        // If P_k(t) = sum_{i=0}^{deg} a[i] t^i,
        // then dp(t) = sum_{i=1}^{deg} i * a[i] t^{i-1}.
        for (int i = 1; i <= deg; ++i) {
            dp[i - 1] += i * a[i];
        }

        // P_{k+1}(t) = (1 + t^2) * dp(t)
        // = dp(t) + t^2 * dp(t)
        // So in coefficients:
        //   b[j]        += dp[j]
        //   b[j + 2]    += dp[j]
        for (int j = 0; j <= deg - 1; ++j) {  // dp has degree deg-1
            b[j] += dp[j];
            if (j + 2 <= max_deg) {
                b[j + 2] += dp[j];
            }
        }

        // Now b is P_{k+1}, degree increases by 1
        deg = deg + 1;

        // Copy b into a
        for (int i = 0; i <= deg; ++i) {
            a[i] = b[i];
        }
    }

    // Now a[] holds the coefficients of P_n(t), degree = deg = n+1.
    // Evaluate P_n(t) at t = tan(x) using Horner's method.
    double result = 0.0;
    for (int i = deg; i >= 0; --i) {
        result = result * t + a[i];
    }

    free(a);
    free(dp);
    free(b);

    return result;
}

/*** Example usage ***/
int main(void) {
    double x = 0.0;
    double taylor_denominator = 1.0;
                     
    for (int n = 0; n <= 50; ++n) {
        double val = tan_derivative(n, x);
        taylor_denominator *= (n == 0.0) ? 1 : n;
        double taylor_eval = val / taylor_denominator;

        printf("Derivative: %d; Derivative Value: %.17g; Taylor Denominator: %.17g; Taylor Coeff: %.17g\n",
               n, val, taylor_denominator, taylor_eval);
    }
    return 0;
}
