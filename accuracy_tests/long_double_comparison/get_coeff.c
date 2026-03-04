// taylor_tan_ldhex.c
// Build:  gcc -std=c11 -O2 -Wall -Wextra taylor_tan_ldhex.c -o taylor_tan_ldhex
// Run:    ./taylor_tan_ldhex 27
//
// Prints Taylor coefficients a_n for tan(x) = sum a_n x^n around 0,
// as long double hex floats (C99 %La). Only odd n are nonzero.

#include <stdio.h>
#include <stdlib.h>

static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n, sz);
    if (!p) {
        perror("calloc");
        exit(1);
    }
    return p;
}

int main(int argc, char **argv) {
    int N = 27; // max degree to print
    if (argc >= 2) {
        N = atoi(argv[1]);
        if (N < 1) N = 1;
    }

    // P0(t) = t  => coefficients: P[0]=0, P[1]=1
    // We'll store P_n(t) = sum_{j=0..deg} P[j] * t^j
    int deg = 1;
    long double *P = (long double *)xcalloc((size_t)(N + 3), sizeof(long double));
    P[0] = 0.0L;
    P[1] = 1.0L;

    long double fact = 1.0L; // n!
    // We will iteratively build P_{n} from P_{n-1} and at each step output a_n = P_n(0)/n!
    for (int n = 0; n < N; ++n) {
        // Build derivative D(t) = d/dt P(t)
        // If P has degree 'deg', then D has degree deg-1.
        int ddeg = (deg > 0) ? (deg - 1) : 0;
        long double *D = (long double *)xcalloc((size_t)(N + 3), sizeof(long double));
        for (int j = 1; j <= deg; ++j) {
            D[j - 1] = (long double)j * P[j];
        }

        // Pnext(t) = (1 + t^2) * D(t) = D(t) + t^2*D(t)
        int next_deg = ddeg + 2;
        if (next_deg > N + 2) next_deg = N + 2;

        long double *Pnext = (long double *)xcalloc((size_t)(N + 3), sizeof(long double));
        for (int k = 0; k <= ddeg; ++k) {
            Pnext[k] += D[k];
            if (k + 2 <= N + 2) Pnext[k + 2] += D[k];
        }

        free(D);
        free(P);
        P = Pnext;
        deg = next_deg;

        // Now we have P_{n+1}. Its constant term is tan^{(n+1)}(0).
        int nn = n + 1;
        fact *= (long double)nn;               // (n+1)!
        long double deriv_at_0 = P[0];         // P_{n+1}(0)
        long double a_nn = deriv_at_0 / fact;  // Taylor coefficient

        if ((nn & 1) == 1) {
            // Print: n and the coefficient as hex long double.
            // Example format: n=3  a_n=0x1.1111111111111p-2
            printf("n=%d  a_n=%La\n", nn, a_nn);
        }
    }

    free(P);
    return 0;
}

