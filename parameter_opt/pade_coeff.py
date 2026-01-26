import sympy as sp

def pade_coeffs_tan(m: int, n: int):
    x = sp.Symbol('x')

    # Taylor coefficients a_0..a_{m+n} of tan(x)
    series = sp.series(sp.tan(x), x, 0, m + n + 1).removeO().expand()
    a = [series.coeff(x, k) for k in range(m + n + 1)]

    # Solve for q_1..q_n (with q_0 = 1) from:
    # a_k + sum_{j=1..n} q_j a_{k-j} = 0  for k = m+1..m+n
    q_syms = sp.symbols(f"q1:{n+1}")
    eqs = []
    for k in range(m + 1, m + n + 1):
        eqs.append(sp.Eq(a[k] + sum(q_syms[j-1] * a[k - j] for j in range(1, n + 1)), 0))

    sol = sp.solve(eqs, q_syms, dict=True)
    if not sol:
        raise RuntimeError("No solution found for Q coefficients (unexpected for tan at 0).")
    sol = sol[0]

    q = [sp.Integer(1)] + [sp.simplify(sol[s]) for s in q_syms]

    # Compute p_k = sum_{j=0..min(k,n)} q_j a_{k-j} for k=0..m
    p = []
    for k in range(m + 1):
        pk = sum(q[j] * a[k - j] for j in range(0, min(k, n) + 1))
        p.append(sp.simplify(pk))

    return p, q  # ascending powers


def print_c_array(name: str, coeffs, digits: int = 20):
    # Convert exact SymPy rationals to decimal strings for C.
    # If you prefer exact rationals in C, say so and I’ll output as fractions.
    def to_c_double(c):
        c = sp.nsimplify(c)  # keep exact if possible
        v = sp.N(c, digits)
        s = sp.sstr(v)  # e.g., '1.234e-5'
        if 'e' in s:
            s = s.replace('e', 'E')
        return s

    print(f"double {name}[] = {{")
    for c in coeffs:
        print(f"  {to_c_double(c)},")
    print("};\n")


if __name__ == "__main__":
    m = 30
    n = 30

    p_coeffs, q_coeffs = pade_coeffs_tan(m, n)

    # If you want the raw Python lists, you already have p_coeffs/q_coeffs here.
    # For C paste:
    print_c_array("PADE_COEFF_P", p_coeffs, digits=25)
    print_c_array("PADE_COEFF_Q", q_coeffs, digits=25)
