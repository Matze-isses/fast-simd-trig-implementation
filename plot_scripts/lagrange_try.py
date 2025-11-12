import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


used_params = {
    "lower_bound": np.pi/4,
    "upper_bound": 1.5,
    "max_degree": 20,
    "wanted_error": 1e-10, 
    "show_every_plot": True,
    "show_final_plot": True,
    "test_size": 100000
}


def taylor(points):
    """
    Returns a function that computes the Taylor polynomial of tan(x)
    of the given degree around point a (default 0).
    """
    degree = len(points) - 1 

    x = sp.symbols('x')
    f = sp.tan(x)
    
    # Taylor expansion around a
    taylor_series = sp.series(f, x, 0, degree + 1)
    taylor_poly = taylor_series.removeO()
    
    # Convert to a numerical function
    f_taylor = sp.lambdify(x, taylor_poly, 'numpy')
    
    return f_taylor


def newton(points):
    """
    Return a callable Newton interpolation function for given data points.

    Parameters
    ----------
    points : list of (x, y)
        Distinct x-values are required.

    Returns
    -------
    f : function
        Interpolating function f(x). Supports scalar or NumPy array inputs.
    """
    xs = np.array([float(p[0]) for p in points], dtype=float)
    ys = np.array([float(p[1]) for p in points], dtype=float)

    if len(set(xs.tolist())) != len(xs):
        raise ValueError("All x-values must be distinct for Newton interpolation.")

    n = len(xs)
    # Divided differences: c[k] = k-th order coefficient
    c = ys.astype(float).copy()
    for k in range(1, n):
        c[k:n] = (c[k:n] - c[k-1:n-1]) / (xs[k:n] - xs[0:n-k])

    # Horner-like evaluation for Newton form:
    # f(x) = c0 + (x-x0)[c1 + (x-x1)[c2 + ... ]]
    def f(x):
        x = np.array(x, dtype=float)
        val = np.full_like(x, c[-1], dtype=float)
        for k in range(n - 2, -1, -1):
            val = val * (x - xs[k]) + c[k]
        return float(val) if val.size == 1 else val

    return f



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


def hermite(points):
    # Unpack points
    xs = []
    ys = []
    dys = []

    for p in points:
        if len(p) == 2:
            x, y = p
            dy = None
        elif len(p) == 3:
            x, y, dy = p
        else:
            raise ValueError("Each point must be (x, y) or (x, y, y').")
        xs.append(float(x))
        ys.append(float(y))
        dys.append(None if dy is None else float(dy))

    # Distinct x-values required
    if len(set(xs)) != len(xs):
        raise ValueError("All x-values must be distinct for Hermite interpolation.")

    n = len(xs)
    z = np.zeros(2 * n)
    Q = np.zeros((2 * n, 2 * n))

    # Fill z and first column (f-values)
    for i in range(n):
        z[2 * i] = xs[i]
        z[2 * i + 1] = xs[i]
        Q[2 * i][0] = ys[i]
        Q[2 * i + 1][0] = ys[i]
        if dys[i] is not None:
            Q[2 * i + 1][1] = dys[i]
        else:
            Q[2 * i + 1][1] = 0.0  # assume 0 derivative if not given
        if i > 0:
            Q[2 * i][1] = (Q[2 * i][0] - Q[2 * i - 1][0]) / (z[2 * i] - z[2 * i - 1])

    # Build divided differences table
    for j in range(2, 2 * n):
        for i in range(2 * n - j):
            Q[i][j] = (Q[i + 1][j - 1] - Q[i][j - 1]) / (z[i + j] - z[i])

    coeffs = Q[0, :]  # Newton form coefficients

    def f(x):
        """Evaluate Hermite interpolant at x (supports scalar or array)."""
        x = np.array(x, dtype=float)
        res = np.zeros_like(x)
        for i in range(len(x)):
            val = coeffs[-1]
            for k in range(2 * n - 2, -1, -1):
                val = val * (x[i] - z[k]) + coeffs[k]
            res[i] = val
        return res if res.size > 1 else float(res)

    return f






def generate_poly_plot(
        poly_name, 
        func_gen, 
        lower_bound=np.pi/4,
        upper_bound=1.5,
        max_degree=20,
        wanted_error=1e-10, 
        show_every_plot=False,
        show_final_plot=True,
        test_size=100000
    ):

    dev_points = [lower_bound, upper_bound]

    current_degree = 0
    max_error = 1000 
    
    x_grid = np.linspace(lower_bound, upper_bound, test_size)
    y_true = np.tan(x_grid)

    y_lagrange = np.zeros(test_size)
    abs_error = np.abs(y_lagrange - y_true)

    func = None


    while current_degree < max_degree and wanted_error < max_error:
        current_degree += 1

        x_vals = np.array(dev_points)
        points = [(x, np.tan(x)) for x in x_vals]

        func = func_gen(points)
        y_lagrange = func(x_grid)

        abs_error = np.abs(y_lagrange - y_true)
        max_error = np.max(abs_error)

        problem_index = np.argmax(abs_error)
        problem_point = x_grid[problem_index]
        running = 0
        while problem_point in dev_points:
            problem_index += 1 if problem_index < x_grid.shape[0]/2 else -1
            problem_point = x_grid[problem_index]
            running += 1
            assert running < 100, "To long while iteration"

        dev_points.append(problem_point)

        
        if show_every_plot:
            print(problem_point)
            print(f"Mean Abs Error: {np.mean(abs_error)}")

            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1]})
            plt.title(f"{poly_name} Polynom with degree: {current_degree+1}")

            ax1.scatter(np.array(dev_points[:-1]), np.tan(np.array(dev_points[:-1])), c='b', label='Problem Points')
            ax1.scatter(np.array(dev_points[-1]), np.tan(np.array(dev_points[-1])), c='r', label='Problem Points')
            ax1.plot(x_grid, y_lagrange, label='Lagrange')
            ax1.plot(x_grid, y_true, label='Tan')
            ax1.set_ylabel('y')
            ax1.set_ylabel('x')
            ax1.legend()

            ax2.plot(x_grid, abs_error, label='Error')
            ax2.set_yscale('log')
            ax2.set_xlabel('x')
            ax2.set_ylabel('Error')
            ax2.legend()

            plt.show()

    if show_final_plot:
        fig, (ax1, ax2) = plt.subplots( 2, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1]})
        plt.title(f"{poly_name} Polynom with degree: {current_degree+1}")

        ax1.scatter(np.array(dev_points[:-1]), np.tan(np.array(dev_points[:-1])), c='b', label='Problem Points')
        ax1.scatter(np.array(dev_points[-1]), np.tan(np.array(dev_points[-1])), c='r', label='Problem Points')
        ax1.plot(x_grid, y_lagrange, label='Lagrange')
        ax1.plot(x_grid, y_true, label='Tan')
        ax1.set_ylabel('y')
        ax1.set_ylabel('x')
        ax1.legend()

        ax2.plot(x_grid, abs_error, label='Error')
        ax2.set_yscale('log')
        ax2.set_xlabel('x')
        ax2.set_ylabel('Error')
        ax2.legend()

        plt.show()

    print(f"Mean Abs Error: {np.mean(abs_error)}")
    print(f"Max  Abs Error: {np.max(abs_error)}")


def show_taylor(**kwargs):
    generate_poly_plot("Taylor", taylor, **kwargs)

def show_lagrange(**kwargs):
    generate_poly_plot("Lagrange", lagrange, **kwargs)

def show_hermit(**kwargs):
    generate_poly_plot("Hermit", hermite, **kwargs)

def show_newton(**kwargs):
    generate_poly_plot("Newton", newton, **kwargs)


if __name__ == "__main__":
    # poly_type = input("Show (a) lagrange (b) hermit -/ Polinomials")
    # show_taylor(**used_params)
    show_lagrange(**used_params)
    # show_hermit(**used_params)
    # show_newton(**used_params)



