import ctypes
import matplotlib.pyplot as plt
import numpy as np

# Load the shared library
lib = ctypes.CDLL("./libtest_object.so")

lib.test_glibc_sin_time.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.test_glibc_sin_time.restype = ctypes.c_double

# Declare argument and return types
lib.test_sin_time.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.test_sin_time.restype = ctypes.c_double

lib.test_tan_time.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.test_tan_time.restype = ctypes.c_double

lib.test_sin_accuracy.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.test_sin_accuracy.restype = ctypes.c_double

lib.test_tan_accuracy.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.test_tan_accuracy.restype = ctypes.c_double


plausible_test_size = 10000000
plot_quality = (1920, 1080)

dpi = 200
width = plot_quality[0] / dpi
height = plot_quality[1] / dpi


def taylor2error(max_degree=10, block=False):
    lower_bound = 0
    upper_bound = 1000
    taylor_degrees = list(range(1, max_degree+1))
    errors = []

    # taylor to error comparison
    for i in taylor_degrees:
        print(f"Current Taylor Degree: {i}/{max_degree}", end=" ")
        error = lib.test_sin_accuracy(i, lower_bound, upper_bound, plausible_test_size)
        errors.append(error)
        print(f"Had Error: {error}")
    
    plt.figure(figsize=(width, height), dpi=dpi)

    # Plot the error on a semilog y-axis
    plt.semilogy(taylor_degrees, errors, marker='o', linestyle='-', label='Error')

    # Add labels on each point with the actual error values
    for x, y in zip(taylor_degrees, errors):
        plt.text(x, y, f"{y:.1e}", ha='left', va='bottom', fontsize=9)

    plt.xlabel("Taylor degree")
    plt.ylabel("Error")
    plt.title(f"Degree vs Error for n={plausible_test_size}; lower bound={lower_bound}; upper bound={upper_bound}")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()

    plt.savefig("./plots/plot_error.png", dpi=dpi, bbox_inches='tight')
    plt.show(block=block)


def taylor2time(max_degree=10, block=False):
    test_size = 1000000000

    lower_bound = 0
    upper_bound = 10000
    taylor_degrees = list(range(1, max_degree+1))
    times = []

    # taylor to error comparison
    for i in taylor_degrees:
        print(f"Current Taylor Degree: {i}/{max_degree}", end=" ")
        time = lib.test_sin_time(i, lower_bound, upper_bound, test_size)
        times.append(time)
        print(f"Was calculated in: {time}")

    plt.figure(figsize=(width, height), dpi=dpi)
    plt.plot(taylor_degrees, times, marker='x', linestyle='-')

    # Add labels on each point with the actual error values
    for x, t in zip(taylor_degrees, times):
        plt.text(x, t, f"{t:6.1f} ms", ha='left', va='bottom', fontsize=9)

    plt.xlabel("Taylor degree")
    plt.ylabel("Time")
    plt.title(f"Degree vs Time for n={test_size}; lower bound={lower_bound}; upper bound={upper_bound}")
    
    plt.savefig("./plots/plot_time.png", dpi=dpi, bbox_inches='tight')
    plt.show(block=block)


def bound2error(max_bound_exponent=6, block=False):
    taylor_degree = 20
    lower_bounds = [10 ** i for i in range(max_bound_exponent)]
    upper_bounds = [10 ** i for i in range(1, max_bound_exponent+1)]
    errors = []

    # taylor to error comparison
    for i, (lower_bound, upper_bound) in enumerate(zip(lower_bounds, upper_bounds)):
        print(f"Current Bound Exponent: {i}/{max_bound_exponent}", end=" ")
        error = lib.test_sin_accuracy(taylor_degree, lower_bound, upper_bound, plausible_test_size)
        errors.append(error)
        print(f"Had Error: {error}")
    
    plt.figure(figsize=(width, height), dpi=dpi)

    # Plot the error on a semilog y-axis
    plt.semilogy(list(range(1, max_bound_exponent+1)), errors, marker='o', linestyle='-', label='Error')

    # Add labels on each point with the actual error values
    for x, ub, lb, y in zip(list(range(1, max_bound_exponent+1)), upper_bounds, lower_bounds, errors):
        plt.text(x, y, f"[{lb:.1e}, {ub:.1e}]", ha='left', va='bottom', fontsize=9)

    plt.xlabel("Upper Bound Exponent (10^x)")
    plt.ylabel("Error")
    plt.title(f"Bounds vs Error for n={plausible_test_size}; taylor degree={taylor_degree}")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()

    plt.savefig("./plots/plot_bounds.png", dpi=dpi, bbox_inches='tight')
    plt.show(block=block)

def glibc_bounds2time(max_bound_exponent=6, block=False):
    test_size = 1000000000
    lower_bounds = [10 ** i for i in range(max_bound_exponent)]
    upper_bounds = [10 ** i for i in range(1, max_bound_exponent+1)]
    times = []

    # taylor to error comparison
    for i, (lower_bound, upper_bound) in enumerate(zip(lower_bounds, upper_bounds)):
        print(f"Current Bound Exponent: {i}/{max_bound_exponent}", end=" ")
        time = lib.test_glibc_sin_time(0, lower_bound, upper_bound, test_size)
        times.append(time)
        print(f"Had Time: {time}")
    
    plt.figure(figsize=(width, height), dpi=dpi)

    # Plot the error on a semilog y-axis
    plt.plot(list(range(1, max_bound_exponent+1)), times, marker='o', linestyle='-', label='Time')

    plt.xlabel("Upper Bound Exponent (10^x)")
    plt.ylabel("Time")
    plt.title(f"Bounds vs Time of glibc for n={test_size}")
    plt.legend()

    plt.savefig("./plots/plot_glibc_bounds.png", dpi=dpi, bbox_inches='tight')
    plt.show(block=block)

# taylor2error(30)
# taylor2time(30)
# bound2error(30)
glibc_bounds2time(20, True)
