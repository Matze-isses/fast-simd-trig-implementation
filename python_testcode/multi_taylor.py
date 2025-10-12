import math
from itertools import product
import numpy as np
import matplotlib.pyplot as plt


BASE_INTERVAL = [0, math.pi/2]


def taylor_sin(n, a):
    """
    Return a function P(x) giving the n-th order Taylor approximation
    of sin(x) around x = a.

    P(x) = sum_{k=0..n} sin^(k)(a) * (x-a)^k / k!

    Works with scalars and NumPy arrays for x.
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("n must be a non-negative integer")

    sa, ca = math.sin(a), math.cos(a)  # cache sin(a), cos(a)

    def kth_deriv_at_a(k):
        r = k & 3  # k % 4
        if r == 0: return sa
        if r == 1: return ca
        if r == 2: return -sa
        return -ca

    def P(x):
        dx = x - a
        res = 0.0
        term = 1.0  # (dx^0)/(0!)
        for k in range(n + 1):
            res += kth_deriv_at_a(k) * term
            # update for next k: (dx^(k+1))/(k+1)! = (dx^k/k!) * (dx/(k+1))
            if k != n:
                term = term * dx / (k + 1)
        return res

    return P


def own_sin(values, splits: int = 5, order: int = 3):
    taylor_functions = []

    a, b = BASE_INTERVAL
    step = (b - a) / splits
    center_points = [a + (i + 0.5) * step for i in range(splits)]

    for center in center_points: 
        taylor_functions.append(taylor_sin(order, center))

    y = []

    for x in values:
        func_index = math.floor(x / step)
        if func_index == splits:
            func_index -= 1
        y.append(taylor_functions[func_index](x))

    return y
        

def get_error_of(splits, order):
    x = np.linspace(BASE_INTERVAL[0], BASE_INTERVAL[1], 50000)
    y = np.array(own_sin(x, order=order, splits=splits))
    print(f"Done for: {splits} and {order}")
    return sum([abs(own - real) for own, real in zip(y, np.sin(x))])

split_tests = list(range(25))[5:]
order_tests = list(range(15))[5:]
coords = [(s, o, get_error_of(s, o)) for s, o in product(split_tests, order_tests)]

nx, ny = len(split_tests), len(order_tests)

X, Y = np.meshgrid(split_tests, order_tests, indexing='xy')  # shapes (ny, nx)
Z = np.array([z for _, _, z in coords]).reshape(nx, ny).T   # -> shape (ny, nx)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(X, Y, Z)  # default colormap; no explicit colors
ax.set_xlabel('split (s)')
ax.set_ylabel('order (o)')
ax.set_zlabel('error')
plt.show()


# x = np.linspace(BASE_INTERVAL[0], BASE_INTERVAL[1], 10000)
# y = np.array(own_sin(x, order=4, splits=4))

# plt.plot(x, y, c='b')
# plt.plot(x, np.sin(x), c='r')

# error = sum([(own - real) ** 2 for own, real in zip(y, np.sin(x))])
# print(error)
# plt.show()
 



