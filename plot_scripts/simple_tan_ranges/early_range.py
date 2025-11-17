import numpy as np
import time
import sympy as sp
import matplotlib.pyplot as plt

CLOSE_TO_SINGULARITY = np.pi/2 - 0.0000000001
START_POINT = 1.5


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


def compare_from_left_to_true(func, start_x=0, max_error=1e-15):
    upper = CLOSE_TO_SINGULARITY
    lower = start_x
    error = np.inf

    x_test = 0
    
    while abs(error - max_error) > max_error:
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
    

def compare_from_right_to_true(func, start_x=0, max_error=1e-15):
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
    


taylor_degree = 15
taylor_poly = taylor(taylor_degree)

lagrange_dev_points = [0.39269883169872416,0.3969635890885462,0.4054027821186821,0.42868615832886153,0.45763231552994654,0.4874267117698911,0.5407990116640636,0.5722270534603305,0.6115366496593349,0.6435851626043971,0.6768824719122251,0.7267832749992127,0.7490809660299766,0.7665916043368057,0.7795665200052884,0.7853981633974483]

deg_lagrange = len(lagrange_dev_points)-1

lagrange_points = [(x, np.tan(x)) for x in lagrange_dev_points]
lagrange_poly = lagrange(lagrange_points)

end_taylor_range = compare_from_left_to_true(taylor_poly)
print(f"The End of Taylor range is:   {end_taylor_range}")
end_taylor_range = np.pi/8

end_lagrange_range = compare_from_left_to_true(lagrange_poly, start_x=end_taylor_range)
print(f"The End of Lagrange range is: {end_lagrange_range}")


def combined_func(x_vals):
    y = []

    for x in x_vals:
        if x < end_taylor_range:
            y.append(taylor_poly(x))
        elif x < end_lagrange_range:
            y.append(lagrange_poly(x))
        elif np.pi/2 - end_taylor_range > x > end_lagrange_range:
            y.append(1 / lagrange_poly(np.pi/2 - x))
        elif x > np.pi/2 - end_taylor_range:
            y.append(1/taylor_poly(np.pi/2 - x))
        else:
            y.append(0)

    return np.array(y)



x_vals = np.linspace(START_POINT, CLOSE_TO_SINGULARITY, 1000000)

print("Max Difference: ", np.max(combined_func(x_vals) - np.tan(x_vals)))

# Create Full HD figure with better aspect ratio
fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, sharex=True, figsize=(9.6, 5.4), gridspec_kw={'height_ratios': [1, 1, 1]}
)

# --- AX1 ---
ax1.plot(x_vals, taylor_poly(x_vals), label='Taylor approx', c='g', alpha=0.7, linewidth=1)
ax1.plot(x_vals, lagrange_poly(x_vals), label='Lagrange approx', c='b', alpha=0.7, linewidth=1)
ax1.plot(x_vals, 1/taylor_poly(np.pi/2 - x_vals), label='1/Taylor approx', c='yellowgreen', alpha=0.7, linewidth=1)
ax1.plot(x_vals, 1/lagrange_poly(np.pi/2 - x_vals), label='1/Lagrange approx', c='lightblue', alpha=0.7, linewidth=1)
ax1.plot(x_vals, np.tan(x_vals), label='True Tan', linewidth=2.5, c='r')

ax1.set_ylim((-0.3, 10))
ax1.set_title("Approximation vs True Tan", fontsize=11)
ax1.legend(loc='upper left', fontsize=8, ncol=2, frameon=False)

# --- AX2 ---
ax2.plot(x_vals, combined_func(x_vals), label="Combined approaches", c='purple', linewidth=1.2)
ax2.plot(x_vals, np.tan(x_vals), label="True Tan", c='r', linewidth=2)

# Highlight ranges
ax2.plot([0, end_taylor_range], [-0.2, -0.2], lw=3, label=f'Taylor Range {taylor_degree}', c='g')
ax2.plot([np.pi/2-end_taylor_range, np.pi/2], [-0.2, -0.2], lw=3, label=f'1/Taylor Range {taylor_degree}', c='yellowgreen')
ax2.plot([end_taylor_range, np.pi/4], [-0.2, -0.2], lw=3, c='b', label=f'Lagrange Range {deg_lagrange}')
ax2.plot([np.pi/4, np.pi/2-end_taylor_range], [-0.2, -0.2], lw=3, c='lightblue', label=f'1/Lagrange Range {deg_lagrange}')

ax2.set_ylim((-0.3, 10))
ax2.set_title("Combined Approximation Behavior", fontsize=11)
ax2.legend(loc='upper left', fontsize=8, ncol=2, frameon=False)

# --- AX3 ---
ax3.plot(x_vals, combined_func(x_vals) - np.tan(x_vals), label="Func - Tan", c='k', linewidth=1)
ax3.set_yscale('symlog', linthresh=1e-15)
ax3.set_title("Error (Func - Tan)", fontsize=11)
ax3.legend(loc='upper left', fontsize=8, frameon=False)

# --- Layout & Save ---
fig.suptitle("Main Result", fontsize=14, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.97])  # adjust spacing to avoid overlap

plt.savefig("./plots/main_result.png", dpi=200, bbox_inches='tight')
plt.show()
