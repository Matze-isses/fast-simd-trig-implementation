import numpy as np
import math

def _poly_derivative(coeffs):
    """
    Given coeffs for P(t) = sum_{i=0}^k a_i t^i,
    return coeffs for P'(t).
    """
    if len(coeffs) <= 1:
        return [0.0]
    deriv = [0.0] * (len(coeffs) - 1)
    for i in range(1, len(coeffs)):
        deriv[i - 1] = coeffs[i] * i
    return deriv

def _poly_mul_1_plus_t2(coeffs):
    """
    Given coeffs for Q(t), return coeffs for (1 + t^2) * Q(t).
    """
    n = len(coeffs)
    res = [0.0] * (n + 2)
    for i, a in enumerate(coeffs):
        res[i]     += a      # 1 * Q(t)
        res[i + 2] += a      # t^2 * Q(t)

    # Optional: trim trailing tiny zeros
    while len(res) > 1 and abs(res[-1]) < 1e-15:
        res.pop()
    return res

def _poly_eval(coeffs, t):
    """
    Evaluate polynomial with given coeffs at t using Horner's method.
    coeffs[i] is coefficient of t^i.
    """
    result = 0.0
    for a in reversed(coeffs):
        result = result * t + a
    return result

def tan_nth_derivative(n, x):
    """
    Return the n-th derivative of tan(x) at the point x (in radians),
    using the iterative scheme:
        P0(t) = t
        P_{k+1}(t) = (1 + t^2) * d/dt P_k(t)
    and then tan^{(n)}(x) = P_n(tan(x)).
    """
    # P0(t) = t  -> coefficients [0, 1] (0 + 1*t)
    coeffs = [0.0, 1.0]

    # Build P_n iteratively:
    # P_{k+1}(t) = (1 + t^2) * (d/dt P_k(t))
    for _ in range(n):
        coeffs = _poly_derivative(coeffs)
        coeffs = _poly_mul_1_plus_t2(coeffs)

    t = math.tan(x)
    return _poly_eval(coeffs, t)

# --- small demo / sanity check ---
if __name__ == "__main__":
    error_dist = []
    error_der = []
    one_over_fak = []
    max_degree = 50

    fak = 1

    for i in range(1, max_degree):
        fak *= (i+1)
        x = (np.pi/4)**(i+1)
        error_dist.append(x)
        one_over_fak.append(1/fak)

    x = 0
    fak = 1
    for n in range(1, max_degree):
        fak *= (n)
        derivative = tan_nth_derivative(n, x)
        error_der.append(derivative)
        taylor_coeff = derivative / fak
        print(taylor_coeff)

    print("Degree; Distance Error; Derivative Error; One Over Faktorial; Total Error")
    for n, ind, der, fak in zip(range(1, max_degree), error_dist, error_der, one_over_fak):
        print(f"{n}; {ind}; {der}; {fak}; {ind * der * fak}")

