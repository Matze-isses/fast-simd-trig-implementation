import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import brentq



def taylor(degree):
    """
    Returns a function that computes the Taylor polynomial of tan(x)
    of the given degree around point a (default 0).
    """

    x = sp.symbols('x')
    f = sp.tan(x)
    
    # Taylor expansion around a
    taylor_series = sp.series(f, x, 0, degree + 1)
    taylor_poly = taylor_series.removeO()
    
    # Convert to a numerical function
    f_taylor = sp.lambdify(x, taylor_poly, 'numpy')
    
    return f_taylor


import math

def nth_derivative_numeric(f, a, n, h=1e-5):
    """
    Approximate the n-th derivative of f at x = a
    using forward finite differences.
    """
    if n == 0:
        return f(a)
    
    # Sample f at a, a+h, ..., a+n*h
    values = [f(a + k*h) for k in range(n+1)]
    
    # Repeatedly apply forward differences n times
    for order in range(1, n+1):
        values = [values[i+1] - values[i] for i in range(len(values)-1)]
    
    # The first value is now Δ^n f(a)
    return values[0] / (h**n)


def taylor_from_callable(f, n, a=0.0, h=1e-5):
    """
    Build the n-th degree Taylor polynomial of a callable f
    around the point a, using numeric derivatives.

    Parameters
    ----------
    f : callable
        Python function f(x) that is sufficiently smooth near x = a.
    n : int
        Degree of the Taylor polynomial.
    a : float, default 0.0
        Expansion point.
    h : float, default 1e-5
        Step size for finite differences.

    Returns
    -------
    P : callable
        A Python function P(x) that evaluates the Taylor polynomial.
    coeffs : list of float
        Coefficients c_k such that
        P(x) = sum_{k=0}^n c_k * (x - a)^k
    """
    # Compute coefficients c_k = f^(k)(a) / k!
    coeffs = [
        nth_derivative_numeric(f, a, k, h) / math.factorial(k)
        for k in range(n+1)
    ]

    def P(x):
        # Evaluate polynomial using Horner's method in (x - a)
        dx = x - a
        result = 0.0
        for c in reversed(coeffs):
            result = result * dx + c
        return result

    return P, coeffs

func = taylor(15)


def adjust_x(x_vals):
    return x_vals + (x_vals * 0.001) ** 2


lower = np.pi/8
upper = np.pi/4

x_vals = np.linspace(lower, upper, 100000)
y_true = np.tan(x_vals)
x_addon = []

for i in range(x_vals.shape[0]):
    x_in = x_vals[i]
    y_true_res = np.tan(x_in)
    same_y_in_mine = 0

    lower_bound = x_in
    upper_bound = 100

    while upper_bound - lower_bound > 0.000000001:
        dist = (upper_bound - lower_bound)
        mid_point = lower_bound + dist/2
        if func(mid_point) > y_true_res:
            upper_bound = mid_point
        else:
            lower_bound = mid_point

    x_addon.append(lower_bound + (upper_bound - lower_bound)/2)

x_addon = np.array(x_addon) - x_vals

f = lambda x: np.tan(x) - func(x)
coeff, adjusted_func = taylor_from_callable(f, n=8, a=0.8)

y_own = func(x_vals)
y_adjusted = func(x_vals) + coeff(x_vals)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1]})

# ax1.plot(x_vals, x_addon, label="Addon")
# ax1.plot(x_vals, y_true, label="True Solution")
ax1.plot(x_vals, np.array(y_own), label="Plane Taylor Poly")
ax1.plot(x_vals, np.array(y_adjusted), label="Adjusted Taylor Poly")
ax1.legend()

error = y_true - np.array(y_own)
adjusted_error = y_true - y_adjusted

ax2.plot(x_vals, error, label="Error")
ax2.plot(x_vals, adjusted_error, label="Adjusted Error")

ax2.set_title("Error Plot")
ax2.legend()

plt.show()
