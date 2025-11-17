import numpy as np
import math
from scipy.optimize import least_squares

def taylor_sin_coeffs(deg):
    """Return [a0,...,adeg] for the Taylor poly of sin(x) up to given degree."""
    coeffs = np.zeros(deg + 1)
    for p in range(1, deg + 1, 2):
        k = (p - 1) // 2
        coeffs[p] = (-1)**k / math.factorial(p)
    return coeffs


def rational_eval(x, a, b):
    """
    Evaluate P(x)/Q(x).
    P(x) = sum_{i=0..m} a[i] x^i
    Q(x) = 1 + sum_{j=1..n} b[j-1] x^j   (so Q(0)=1)
    """
    x = np.asarray(x, dtype=float)
    # Horner for P
    px = np.zeros_like(x)
    for ai in reversed(a):
        px = px * x + ai
    # Direct for Q (small n, no need for Horner)
    qx = np.ones_like(x)
    for j, bj in enumerate(b, start=1):
        qx += bj * x**j
    return px / qx


def fit_rational_sin(m=6, n=5, interval=(-np.pi/2, np.pi/2), n_grid=2000, robust=True):
    """Fit R_{m,n}(x) = P_m(x)/Q_n(x) to sin(x) on the given interval."""
    assert m + n <= 16, "Parameter budget exceeded: need m+n <= 11 (since Q(0)=1)."

    x = np.linspace(interval[0], interval[1], n_grid)
    y = np.sin(x)

    # Initial guess
    a0 = taylor_sin_coeffs(m)
    b0 = np.zeros(n)

    def pack(a, b): return np.concatenate([a, b])
    def unpack(theta): return theta[:m+1], theta[m+1:]

    if robust:
        def residuals(theta):
            a, b = unpack(theta)
            r = rational_eval(x, a, b) - y
            return r / (0.5 + np.abs(r))
    else:
        def residuals(theta):
            a, b = unpack(theta)
            return rational_eval(x, a, b) - y

    theta0 = pack(a0, b0)
    res = least_squares(residuals, theta0, max_nfev=2000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    a_opt, b_opt = unpack(res.x)

    fx = rational_eval(x, a_opt, b_opt)
    err = fx - y
    print(f"Converged: {res.success}, message: {res.message}")
    print(f"L2 RMS error on grid: {np.sqrt(np.mean(err**2)):.3e}")
    print(f"L∞ (max) error on grid: {np.max(np.abs(err)):.3e}")

    return a_opt, b_opt

# Fit
a, b = fit_rational_sin(m=8, n=8, interval=(0, np.pi * 2))
print("Numerator coefficients a (low->high):", a)
print("Denominator coefficients b for x^1..x^n:", b)

# Clean scalar evaluation without warnings
def R_scalar(x):
    return rational_eval(np.array([x]), a, b).item()

print("R(0.1) vs sin(0.1):", R_scalar(0.1), math.sin(0.1))
print("R(1.0) vs sin(1.0):", R_scalar(1.0), math.sin(1.0))
