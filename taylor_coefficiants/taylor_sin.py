import math
import sympy as sp

# --- Inputs ---
func_choice = input("Function (sin, cos, tan): ").strip().lower()
degree = int(input("Degree of the Taylor polynomial: "))

# --- Setup for symbolic poly / derivative printout ---
x = sp.Symbol('x')

def build_poly_from_coeffs(coeffs, a):
    """Return Sympy polynomial sum_{n=0}^{degree-1} coeffs[n]*(x-a)^n"""
    return sp.expand(sum(sp.Rational(1,1)*sp.nsimplify(coeffs[n]) * (x - a)**n for n in range(len(coeffs))))

coeffs = []
a = None  # expansion point

if func_choice in ("sin", "cos"):
    # Preserve your original behavior: expand around a = pi/2
    a = math.pi / 2

    # Use your derivative-cycle logic, selecting the right start (sin vs cos)
    if func_choice == "sin":
        # Your original "sin=True" branch
        for n in range(degree):
            if n % 4 == 0:
                deriv = math.sin(a)
            elif n % 4 == 1:
                deriv = math.cos(a)
            elif n % 4 == 2:
                deriv = -math.sin(a)
            else:  # n % 4 == 3
                deriv = -math.cos(a)
            coeffs.append(deriv / math.factorial(n))
    else:
        # Your original "sin=False" branch (cos first derivative pattern)
        for n in range(degree):
            if n % 4 == 0:
                deriv = math.cos(a)
            elif n % 4 == 1:
                deriv = -math.sin(a)
            elif n % 4 == 2:
                deriv = -math.cos(a)
            else:  # n % 4 == 3
                deriv = math.sin(a)
            coeffs.append(deriv / math.factorial(n))

    # Build + differentiate polynomial symbolically (for nice output)
    a_sym = sp.nsimplify(a)
    poly = build_poly_from_coeffs(coeffs, a_sym)
    dpoly = sp.diff(poly, x)

elif func_choice == "tan":
    # For tan, expand around a = 0 (safe)
    a = 0 

    xs = sp.Symbol('x')
    f = sp.tan(xs)
    coeffs = []
    for n in range(degree):
        if (n % 2 == 0 and a == 0):
            coeffs.append(0.0)
            print("{:.17g},".format(float(coeffs[-1])))
            continue
        deriv_n = sp.diff(f, xs, n)
        deriv_at_a = deriv_n.subs(xs, a)
        coeffs.append(sp.simplify(deriv_at_a / math.factorial(n)))
        print("{:.17g},".format(float(coeffs[-1])))

    poly = build_poly_from_coeffs(coeffs, 0)
    dpoly = sp.diff(poly, x)

else:
    raise ValueError("Unsupported function. Choose one of: sin, cos, tan")

# --- Output ---
print("\nExpansion point a =", a)

for c in coeffs:
    print("{:.17g},".format(float(c)))

