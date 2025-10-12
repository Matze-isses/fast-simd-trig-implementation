import ctypes
from itertools import product
import matplotlib.pyplot as plt
import numpy as np

from ctypes import c_int, c_double, POINTER

# Load the shared library (adjust path/name if needed)
lib = ctypes.CDLL("./libpolyfit.so")

# Declare signature: double compute_mae(int num_coeffs, const double *coeffs)
lib.compute_mae.argtypes = (c_int, POINTER(c_double))
lib.compute_mae.restype = c_double

def taylor_coefficients_sin_pi_over_4(n):
    coeffs = []
    sqrt2_over_2 = (2 ** 0.5) / 2.0  # no libraries used
    inv_fact = 1.0                   # 1/0!
    for k in range(n):
        sign = 1.0 if (k % 4) < 2 else -1.0
        coeffs.append(sign * sqrt2_over_2 * inv_fact)
        inv_fact /= (k + 1.0)        # update to 1/(k+1)!
    return coeffs

def run_c_mae(*coeffs):
    n = len(coeffs)
    arr = (c_double * n)(*coeffs)
    return lib.compute_mae(n, arr)


coeffs = taylor_coefficients_sin_pi_over_4(3)
print(" | ".join([str(c) for c in coeffs]))
points_to_evaluate = 100

x_vals = np.linspace(-0.5, 0.5, points_to_evaluate)
y_vals = np.linspace(-0.5, 0.5, points_to_evaluate)
X, Y = np.meshgrid(x_vals, y_vals)
solution = np.zeros_like(X)

total_runs = points_to_evaluate**2
current_run = 0

for i, j in product(range(X.shape[0]), range(Y.shape[0])):
        x = X[i, j]
        y = Y[i, j]
        solution[i, j] = run_c_mae(0.7071067811865476, 0.7071067811865476, x, y)
        current_run += 1
        print(f"Currently at {current_run}/{total_runs}")


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, solution, cmap='viridis')

ax.set_xlabel('x')
ax.set_ylabel('y+z')
ax.set_zlabel('MAE')
plt.show()
