import numpy as np
import time
import sympy as sp
import matplotlib.pyplot as plt

CLOSE_TO_SINGULARITY = np.pi/2 - 0.0000000001


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


def lagrange(points):
    xs = [float(p[0]) for p in points]
    ys = [float(p[1]) for p in points]

    if len(set(xs)) != len(xs):
        raise ValueError("All x-values must be distinct.")

    def f(x):
        total = 0.0
        n = len(xs)
        for i in range(n):
            xi, yi = xs[i], ys[i]
            term = yi
            for j in range(n):
                if i != j:
                    term *= (x - xs[j]) / (xi - xs[j])
            total += term
        return total

    return f


def compare_from_left_to_true(func, start_x=0, max_error=1e-10):
    upper = CLOSE_TO_SINGULARITY
    lower = start_x
    error = np.inf

    x_test = 0
    
    while abs(error - max_error) > 1e-10:
        x_test = (upper + lower) / 2

        x_true = np.tan(x_test)
        x_approx = func(x_test)

        error = x_true - x_approx

        if abs(error) > max_error:
            upper = x_test
        else:
            lower = x_test

        # print(f"Test: {x_test}; Error: {error}")
    
    return x_test
    

def compare_from_right_to_true(func, start_x=0, max_error=1e-10):
    upper = CLOSE_TO_SINGULARITY
    lower = start_x
    error = 10000

    x_test = 0
    last_error = 100000
    
    while abs(error - max_error) > 1e-15 and last_error > error:
        last_error = error
        x_test = (upper + lower) / 2

        x_true = np.tan(x_test)
        x_approx = func(x_test)

        error = abs(x_true - x_approx)

        if abs(error) > max_error:
            lower = x_test
        else:
            upper = x_test

        # print(f"Test: {x_test}; Error: {error}")
    
    return x_test
    


taylor_degree = 11
taylor_poly = taylor(11)
fraction = lambda x: 1 / (np.pi/2 - x)

last_useable_x = compare_from_left_to_true(taylor_poly)
print(f"The end of the working range for taylor is: {last_useable_x}")

first_useable_x = compare_from_right_to_true(fraction)
print(f"The end of the working 1/x is:              {first_useable_x}")

x_vals = np.linspace(0, CLOSE_TO_SINGULARITY, 100000)

plt.plot(x_vals, np.tan(x_vals), label='True Tan', c='r')

plt.plot(x_vals, taylor_poly(x_vals), label='Taylor approx', c='g', alpha=0.7, linewidth=1)
plt.scatter([last_useable_x], [1], label='End of Taylor Range', c='g')

plt.plot(x_vals, fraction(x_vals), label='1/x approx', alpha=0.7, linewidth=1)
plt.scatter([first_useable_x], [1], label='Start of 1/x Range')

plt.ylim((0, 10))
plt.legend()
plt.show()

