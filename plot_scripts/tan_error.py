import numpy as np
import matplotlib.pyplot as plt


def get_data(path="tan_error_behavior.tsv"):
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    x = data[:, 0]
    err = data[:, 1]
    return x, err


def simple_error_plot(x, err):
    plt.figure(figsize=(14, 8))
    plt.title("tan_simd absolute error vs input")
    plt.xlabel("x")
    plt.ylabel("|tan(x) - tan_simd(x)|")
    plt.grid(True)
    plt.plot(x, err)
    plt.show()


def plot_range(x, err, bounds: tuple = (np.pi/2, 3/2 * np.pi)):
    plot_x, plot_y = [], []

    for _x, _y in zip(x, err):
        if bounds[0] <= _x <= bounds[1]:
            plot_x.append(_x)
            plot_y.append(_y)

    plt.figure(figsize=(14, 8))
    plt.title("tan_simd absolute error vs input")
    plt.xlabel("x")
    plt.ylabel("|tan(x) - tan_simd(x)|")
    plt.grid(True)
    plt.plot(plot_x, plot_y)
    plt.show()



def problem_area(x, err):
    pi = np.pi
    w = 1e-3   # window width

    xmin = x.min()
    xmax = x.max()

    n_min = int(np.floor((2 * xmin / pi - 1) / 2))
    n_max = int(np.ceil((2 * xmax / pi - 1) / 2))

    plt.figure(figsize=(14, 8))
    plt.title(r"tan_simd absolute error near poles $(1+2n)\pi/2$ (centered overlays)")
    plt.xlabel(r"offset from pole $x - (1+2n)\pi/2$")
    plt.ylabel(r"$|tan(x) - tan\_simd(x)|$")
    plt.grid(True)

    count = 0

    for n in range(n_min, n_max + 1):
        pole = (1 + 2*n) * pi / 2.0

        # window [pole, pole + w]
        mask = (x >= pole) & (x <= pole + w)
        if not np.any(mask):
            continue

        xw = x[mask]
        ew = err[mask]

        # center window so all poles overlay at 0
        t = xw - pole

        # give each curve its own legend entry
        label = rf"$x = \frac{{(1+2\cdot{n})\pi}}{{2}}$"
        plt.plot(t, ew, alpha=0.7, label=label)
        count += 1

    plt.xlim(0.0, w)
    plt.legend(loc="upper left", fontsize="small", ncol=2)
    print(f"Plotted {count} pole windows.")
    plt.show()


if __name__ == "__main__":
    x, err = get_data()
    plot_range(x, err)
    # problem_area(x, err)
