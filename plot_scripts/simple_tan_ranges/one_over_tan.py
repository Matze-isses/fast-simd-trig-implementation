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


comparison_function = "tan"

lower = np.pi/8
upper = np.pi/4
# upper = 1.5

x_vals = np.linspace(lower, upper, 1000000)
y_true = np.tan(x_vals)

func = taylor(15)
y_own = func(x_vals)

y = np.linspace(lower, upper, 100000)
dx = []
for y0 in y:
    # Solve true(x) = y0
    f1 = lambda x: func(x) - y0
    x_true = brentq(f1, 0, 10000)

    # Solve approx(x) = y0
    f2 = lambda x: np.tan(x) - y0
    x_approx = brentq(f2, 0, upper)

    dx.append(x_true - x_approx)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1, 1]})

ax1.plot(x_vals, y_true, label="True Solution")

ax1.plot(x_vals, np.array(y_own), label="Taylor Poly")
ax1.legend()

error = y_true - np.array(y_own)

one_over_x = np.exp(23*(x_vals - np.pi/2))
adjusted_correcture = error[-1]/one_over_x[-1]*one_over_x

ax2.plot(x_vals, error, label="Error")
ax2.plot(x_vals, adjusted_correcture, label=comparison_function)
ax2.set_title("Error Plot")
ax2.legend()

ax3.plot(x_vals, error - adjusted_correcture, label="Left Error - Correcture")
ax2.set_title("Correcture Error Plot")

plt.savefig(f"./plots/error_fkt_{comparison_function}.png")

plt.show()
