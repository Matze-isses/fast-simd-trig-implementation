import numpy as np
import sympy as sp
import matplotlib.pyplot as plt



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



x_vals = np.linspace(0, 1.56, 100000)
y_true = np.tan(x_vals)

func = taylor(15)
cutoff = 0.4

y_own = []

for x in x_vals:
    if x < cutoff:
        y_own.append(func(x))
    elif x > np.pi/2 - cutoff:
        y_own.append(1/func(np.pi/2 - x))
    else:
        y_own.append(0)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1]})

ax1.plot(x_vals, y_true, label="True")
ax1.plot(x_vals, np.array(y_own), label="True")
ax1.legend()

ax2.plot(x_vals, y_true - np.array(y_own), label="True")
ax2.legend()

plt.show()
