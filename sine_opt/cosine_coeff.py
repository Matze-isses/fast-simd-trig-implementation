import math

def cosine_taylor_coeffs(n_terms):
    coeffs = []
    for k in range(n_terms):
        if k % 2 == 0:
            n = k // 2
            coeffs.append(((-1)**n) / math.factorial(2*n))
        else:
            coeffs.append(0.0)
    return coeffs

coeffs = cosine_taylor_coeffs(30)

for i, c in enumerate(coeffs):
    print(f"x^{i}: {c}")
