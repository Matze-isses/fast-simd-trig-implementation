import numpy as np
import matplotlib.pyplot as plt


def get_data(path="./tan_error_behavior.tsv"):
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    x = data[:, 0]
    err = data[:, 1]
    return x, err


def simple_error_plot(x, err):
    plt.figure(figsize=(14, 8))
    plt.title("tan_simd absolute error vs input")
    plt.xlabel("x")
    plt.ylabel("tan(x) - tan_simd(x)")
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
    plt.ylabel("tan(x) - tan_simd(x)")
    plt.grid(True)
    plt.plot(plot_x, plot_y)
    plt.show()



def problem_area_right(x, err):
    pi = np.pi
    w = 1e-4   # window width

    xmin = x.min()
    xmax = x.max()

    n_min = int(np.floor((2 * xmin / pi - 1) / 2))
    n_max = int(np.ceil((2 * xmax / pi - 1) / 2))

    plt.figure(figsize=(14, 8))
    plt.title(r"tan_simd absolute error near poles $(1+2n)\pi/2$ (centered overlays)")
    plt.xlabel(r"offset from pole $x - (1+2n)\pi/2$")
    plt.ylabel(r"$tan(x) - tan\_simd(x)$")
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


def problem_area_left(x, err):
    pi = np.pi
    w = 1e-5   # window width: pole - 0.0001 .. pole

    xmin = x.min()
    xmax = x.max()

    # Polstellen: (1 + 2n) * pi/2
    n_min = int(np.floor((2 * xmin / pi - 1) / 2))
    n_max = int(np.ceil((2 * xmax / pi - 1) / 2))

    plt.figure(figsize=(14, 8))
    plt.title(r"tan_simd absolute error left of poles $(1+2n)\pi/2$ (centered overlays)")
    plt.xlabel(r"offset from pole $x - (1+2n)\pi/2$")
    plt.ylabel(r"$tan(x) - tan\_simd(x)$")
    plt.grid(True)

    count = 0

    for n in range(n_min, n_max + 1):
        pole = (1 + 2*n) * pi / 2.0

        # window [pole - w, pole]
        mask = (x >= pole - w) & (x <= pole)
        if not np.any(mask):
            continue

        xw = x[mask]
        ew = err[mask]

        # Zentrieren: jetzt liegt der Bereich in [-w, 0]
        t = xw - pole

        label = rf"$x = \frac{{(1+2\cdot{n})\pi}}{{2}}$"
        plt.plot(t, ew, alpha=0.7, label=label)
        count += 1

    plt.xlim(-w, 0.0)
    plt.legend(loc="upper left", fontsize="small", ncol=2)
    print(f"Plotted {count} pole windows (left side).")
    plt.show()


def compare_correction(path1, path2):
    x1, y1 = get_data(path1)
    x2, y2 = get_data(path2)


    plt.figure(figsize=(14, 8))
    plt.title(fr"Compare {path1} and {path2}")
    plt.xlabel(r"offset from pole $x - (1+2n)\pi/2$")
    plt.ylabel(r"$tan(x) - tan\_simd(x)$")
    plt.grid(True)

    plt.plot(x1, y1, label=f'{path1}')
    plt.plot(x2, y2, label=f'{path2}')

    plt.legend()
    plt.show()



if __name__ == "__main__":
    print("Next: ", 2/3 * np.pi - 0.00001, 2/3 * np.pi + 0.00001)
    x, err = get_data()
    compare_correction('error_second_positive.tsv', 'tan_error_behavior.tsv')
    # simple_error_plot(x, err) 
    # plot_range(x, err)
    # problem_area_right(x, err)
    # problem_area_left(x, err)
