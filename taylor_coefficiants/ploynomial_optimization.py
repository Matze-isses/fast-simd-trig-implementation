import numpy as np
from scipy.optimize import minimize


def mixed_poly_approx(x, m, *args):
    n = len(args)
    top_coeff = args[:(int)(n/2)]
    bot_coeff = args[(int)(n/2):]
    top = 0
    bot = 0

    for i, c, cl in zip(range(n), top_coeff, bot_coeff):
        top += (x - np.pi/2) ** i * c
        bot += (x - np.pi/2) ** i * cl

    return top/bot if bot > 0 else 0


def test_function(*args):
    x_vals = np.random.uniform(0, np.pi, 10000)
    y_vals = np.sin(x_vals)
    error = 0

    for x, y in zip(x_vals, y_vals):
        res = mixed_poly_approx(x, *args) 
        error += abs(res - y)

    print(f"\33[A\rCurrent Error: {error}") #]
    return error


def objective(x):
    return test_function(*x)
    
x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
print("")
result = minimize(objective, x0, method='Nelder-Mead')

print("Optimization successful:", result.success)
print("Minimum value:", result.fun)
print("Best parameters:", result.x)
