import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def taylor_tan_coeffs(n, a=0, numeric=True):
    x = sp.symbols('x')
    a_sym = sp.nsimplify(a)

    if sp.cos(a_sym) == 0:
        raise ValueError("tan(x) has a pole at x = a. Choose a such that cos(a) != 0.")

    series = sp.series(sp.tan(x), x, a_sym, n + 1).removeO()
    h = sp.symbols('h')
    series_h = sp.expand(series.subs(x, a_sym + h))

    coeffs = [sp.simplify(series_h.coeff(h, k)) for k in range(n + 1)]

    if numeric:
        coeffs = [float(c) for c in coeffs]

    return coeffs


def taylor_poly_from_coeffs(coeffs, a=0):
    """
    Given coefficients [c0, c1, ..., cn] for a Taylor polynomial in (x - a),
    return a Python function p(x) that evaluates:

        p(x) = Σ_{k=0}^n c_k * (x - a)^k

    Uses Horner's method for numerical stability.
    """
    def p(x):
        dx = x - a
        result = 0.0
        # Horner's method: ((((c_n)*dx + c_{n-1})*dx + ...) * dx + c_0)
        for c in reversed(coeffs):
            result = result * dx + c
        return result

    return p


lower = np.pi/8
upper = np.pi/4 

coeffs = taylor_tan_coeffs(n=21, a=3/16 * np.pi, numeric=True)
print("Coefficients:", coeffs)

approx_tan = taylor_poly_from_coeffs(coeffs, a=3/16 * np.pi)

x_vals = np.linspace(lower, upper, 100000)
y_true = np.tan(x_vals)

y_own = approx_tan(x_vals)




fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1]})

ax1.plot(x_vals, y_true, label="True")
ax1.plot(x_vals, np.array(y_own), label="OWN")
ax1.legend()

ax2.plot(x_vals, y_true - np.array(y_own), label="Error")
ax2.legend()

plt.show()
