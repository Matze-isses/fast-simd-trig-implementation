import ctypes
import matplotlib.pyplot as plt
import numpy as np

# Load the shared library
lib = ctypes.CDLL("./libtest_object.so")

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


def taylor2error(max_degree=10):
    lower_bound = 0
    upper_bound = 10000
    taylor_degrees = list(range(1, max_degree+1))
    errors = []

    # taylor to error comparison
    for i in taylor_degrees:
        print(f"Current Taylor Degree: {i}/{max_degree}", end=" ")
        error = lib.test_sin_accuracy(i, lower_bound, upper_bound, plausible_test_size)
        errors.append(error)
        print(f"Had Error: {error}")
    
    plt.plot(taylor_degrees, np.log(np.array(errors)), "x-") 
    plt.xlabel("Taylor Degree")
    plt.ylabel("Log Error")
    plt.title("Log Error for Increasing Taylor Degree")
    plt.savefig("./plots/deg2error")
    plt.show()


taylor2error()

# Example usage
taylor_degree = 10
lower_bound = 0.0
upper_bound = 3.14

test_size = 1000000

time_sin = lib.test_sin_time(taylor_degree, lower_bound, upper_bound, test_size)
print(f"sin time: {time_sin:.6f} s")

time_tan = lib.test_tan_time(taylor_degree, lower_bound, upper_bound, test_size)
print(f"tan time: {time_tan:.6f} s")

accuracy_sin = lib.test_sin_accuracy(taylor_degree, lower_bound, upper_bound, test_size)
print(f"sin mean abs error: {accuracy_sin:.6e}")

accuracy_tan = lib.test_tan_accuracy(taylor_degree, lower_bound, upper_bound, test_size)
print(f"tan mean abs error: {accuracy_tan:.6e}")
