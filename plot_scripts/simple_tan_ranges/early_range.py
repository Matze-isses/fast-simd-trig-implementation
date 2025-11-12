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
lagrange_dev_points = [0.19634954083686207,0.21695837360101924,0.2547110508438543,0.2975585889106345,0.35012499655863705,0.44374608248110736,0.49747151060090783,0.5501646861130133,0.6025439602473401,0.6430674208055671,0.7160011986194357,0.7499327302443118,0.774622280919591,0.8]

lagrange_points = [(x, np.tan(x)) for x in lagrange_dev_points]
lagrange_poly = lagrange(lagrange_points)
fraction = lambda x: 1 / (np.pi/2 - x)

end_taylor_range = compare_from_left_to_true(taylor_poly)
print(f"The End of Taylor range is:   {end_taylor_range}")

end_lagrange_range = compare_from_left_to_true(lagrange_poly, start_x=end_taylor_range)
print(f"The End of Lagrange range is: {end_lagrange_range}")

start_fraction_range = compare_from_right_to_true(fraction)
print(f"The Start of 1/x is:          {start_fraction_range}")



def combined_func(x_vals):
    y = []

    for x in x_vals:
        if x < end_taylor_range:
            y.append(taylor_poly(x))
        elif x < end_lagrange_range:
            y.append(lagrange_poly(x))
        elif x > start_fraction_range:
            y.append(fraction(x))
        else:
            y.append(0)

    return np.array(y)

x_vals = np.linspace(0, CLOSE_TO_SINGULARITY, 100000)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1, 1]})

ax1.plot(x_vals, np.tan(x_vals), label='True Tan', c='r')

ax1.plot(x_vals, taylor_poly(x_vals), label='Taylor approx', c='g', alpha=0.7, linewidth=1)
ax1.plot(x_vals, lagrange_poly(x_vals), label='Lagrange approx', c='b', alpha=0.7, linewidth=1)
ax1.plot(x_vals, fraction(x_vals), label='1/x approx', c='y', alpha=0.7, linewidth=1)

ax1.plot([0, end_taylor_range], [-0.2, -0.2], lw=3, label='Taylor Range (11 deg)', c='g')
ax1.plot([end_taylor_range, end_lagrange_range], [-0.2, -0.2], lw=3, c='b', label='Lagrange Range (13 deg)')
ax1.plot([start_fraction_range, CLOSE_TO_SINGULARITY], [-0.2, -0.2], lw=3, c='y', label='1/x Range')

ax1.set_ylim((-0.3, 10))
ax1.legend()

ax2.plot(x_vals, combined_func(x_vals), label="Combined approaches")
ax2.plot(x_vals, np.tan(x_vals), label="True Tan")

ax2.plot([0, end_taylor_range], [-0.2, -0.2], lw=3, label='Taylor Range (11 deg)', c='g')
ax2.plot([end_taylor_range, end_lagrange_range], [-0.2, -0.2], lw=3, c='b', label='Lagrange Range (13 deg)')
ax2.plot([start_fraction_range, CLOSE_TO_SINGULARITY], [-0.2, -0.2], lw=3, c='y', label='1/x Range')

ax2.set_ylim((-0.3, 10))
ax2.legend()

ax3.plot(x_vals, combined_func(x_vals) - np.tan(x_vals), label="Func - Tan")
ax3.set_yscale('log')
ax3.legend()

plt.savefig("./plots/problem_areas.png")
plt.show()

