import time
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import basinhopping


total_tests = 1000000
weights = [0.1, 0.2, 0.7]

np.random.seed(42)
x_vals = np.random.uniform(0, 1, (int)(weights[0] * total_tests))
np.append(x_vals, np.random.uniform(1, 1.5, (int)(weights[1] * total_tests)), axis=0)
np.append(x_vals, np.random.uniform(1.5, np.pi/2, (int)(weights[2] * total_tests)), axis=0)

x_vals = np.array([0.8, 1.0, 1.3, 1.4, 1.5, 1.51, 1.53, 1.54, 1.569, 1.57])




def mixed_poly_approx(x, *args):
    top_coeff = 9 

    y1 = np.zeros_like(x, dtype=float)
    y2 = np.zeros_like(x, dtype=float)
    for a in reversed(args[:top_coeff]):
        y1 = y1 * x + a

    for a in reversed(args[top_coeff:]):
        y2 = y2 * x + a

    return y1 / y2


def test_function(*args):
    y_vals = np.tan(x_vals)
    res = mixed_poly_approx(x_vals, *args) 

    error = np.sum(np.abs(res - y_vals))

    # print(f"\33[A\rCurrent Error: {error}") #]
    return error


def objective(x):
    return test_function(*x)


def objective_stable(x):
    return objective(x)


def plot_results(x0, result):
    x_vals = np.linspace(np.pi/4, 1.57, 1000000)
    tan_values = np.tan(x_vals)
    prior_results = mixed_poly_approx(x_vals, *x0)
    optimized_results = mixed_poly_approx(x_vals, *result.x)

    # Create two vertically stacked plots that share the same x-axis
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6), height_ratios=[2, 1])

# --- Top plot: function values ---
    ax1.plot(x_vals, tan_values, lw=2, alpha=0.6, linestyle="-", label="tan(x)", c="#aa0000")
    ax1.plot(x_vals, optimized_results, lw=1, linestyle="-", label="Optimized", c="#00aa00")
    ax1.plot(x_vals, prior_results, lw=1, linestyle="-", label="Prior Results", c="#5011a9")
    ax1.set_ylabel("Function Value")
    ax1.legend()
    ax1.grid(True)

# --- Bottom plot: absolute error ---
    error_optimized = np.abs(tan_values - optimized_results)
    error_prior = np.abs(tan_values - prior_results)

    ax2.plot(x_vals, error_optimized, lw=2, linestyle="-", label="Optimized", c="#00aa00")
    # ax2.plot(x_vals, error_prior, lw=2, linestyle="--", alpha=0.6, label="Prior Results", c="#5011a9")
    ax2.set_xlabel("x")
    ax2.set_ylabel("Absolute Error")
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.show()

    
def  run_optimization(x0):
    print("Run Optimization\n")
    minimizer_kwargs = dict(
        method="Powell",
        options={"maxiter": 100_000, "xtol": 1e-6, "ftol": 1e-8}
    )
    result = basinhopping(
        objective_stable, x0,
        niter=60, stepsize=0.25, minimizer_kwargs=minimizer_kwargs, disp=False
    )

    print("Optimization successful:", result.success)
    print("Minimum value:", result.fun)
    print("\n")

    for x in result.x[:8]:
        print(f"{x},")

    print(f"{result.x[9]}\n")

    for x in result.x[9:-1]:
        print(f"{x},")

    print(f"{result.x[-1]}")
    return result


x0 = [
    0.038633276436529584,
    0.9137390202439528,
    -0.07252515250382474,
    0.19142802900467548,
    -0.09060210206571162,
    0.08145495015662974,
    -0.12242428571428568,
    -0.004119999999999985,
    0.02891777777777782,
    1.01038,
    0.014700000000000157,
    -0.3602766666666666,
    -0.042852930416246715,
    0.0035719190611220863
]

print("")
print(objective(x0))
result = run_optimization(x0)

plot_results(x0, result)

