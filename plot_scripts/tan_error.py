import numpy as np
import matplotlib.pyplot as plt


def get_data(path="./tan_error_behavior.tsv"):
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    x = data[:, 0]
    err = data[:, 1]
    return np.array(x), np.array(err)


def simple_error_plot(x, err):
    plt.figure(figsize=(14, 8))
    plt.title("tan_simd absolute error vs input")
    plt.xlabel("x")
    plt.ylabel("tan(x) - tan_simd(x)")
    plt.grid(True)
    plt.plot(x, err)
    plt.show()

def simple_scatter_error_plot(x, err, title="Error Plot", error_metric="tan(x) - tan_simd(x)"):
    plt.figure(figsize=(14, 8))
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel(error_metric)
    plt.grid(True)
    plt.scatter(x, err)
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

        mask = (x >= pole) & (x <= pole + w)
        if not np.any(mask):
            continue

        xw = x[mask]
        ew = err[mask]

        # center window so all poles overlay at 0
        t = xw - pole

        # OPTIONAL but recommended: drop non-finite points (near pole you might get inf/nan)
        finite = np.isfinite(t) & np.isfinite(ew)
        t = t[finite]
        ew = ew[finite]
        if t.size == 0:
            continue

        # CRITICAL FIX: sort by t so the line follows the data monotonically in x
        idx = np.argsort(t, kind="mergesort")  # stable sort
        t = t[idx]
        ew = ew[idx]

        a, b = 0.0, 2e-7
        num = 7

        xs = np.linspace(a, b, num + 1)[1:]

        mask_small = (t > a) & (t <= b)
        ts = t[mask_small]
        es = ew[mask_small]

        print(f"\nn = {n}")

        if ts.size < 2:
            print("  not enough data points in (0, 2e-7]")
        else:
            ys = np.interp(xs, ts, es)
            for xi, yi in zip(xs, ys):
                print(f"  x = {xi:.12e}, err = {yi:.12e}")

        label = rf"$x = \frac{{(1+2\cdot{n})\pi}}{{2}}$"
        plt.plot(t, ew, alpha=0.7, label=label)
        count += 1

    plt.xlim(0.0, w)
    plt.legend(loc="upper left", fontsize="small", ncol=2)
    print(f"Plotted {count} pole windows.")
    plt.show()


def problem_area_left(x, err):
    pi = np.pi
    w = 1e-5   # window width: pole - w .. pole

    xmin = x.min()
    xmax = x.max()

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

        # LEFT window: [pole - w, pole]
        mask = (x >= pole - w) & (x <= pole)
        if not np.any(mask):
            continue

        xw = x[mask]
        ew = err[mask]

        # center window so all poles overlay at 0 (left side => t in [-w, 0])
        t = xw - pole

        # drop non-finite points (near pole you might get inf/nan)
        finite = np.isfinite(t) & np.isfinite(ew)
        t = t[finite]
        ew = ew[finite]
        if t.size == 0:
            continue

        # sort by t so the line follows the data monotonically in x
        idx = np.argsort(t, kind="mergesort")  # stable sort
        t = t[idx]
        ew = ew[idx]

        # mirror the "small interval near pole" reporting, but on the LEFT:
        # interval ( -b, 0 ) in terms of offsets
        b = 2e-7
        num = 7
        xs = -np.linspace(0.0, b, num + 1)[1:]  # [-b/num, ..., -b]

        mask_small = (t < 0.0) & (t >= -b)
        ts = t[mask_small]
        es = ew[mask_small]

        print(f"\nn = {n}")

        if ts.size < 2:
            print("  not enough data points in [-2e-7, 0)")
        else:
            # ensure increasing x for interp (ts is already sorted, but keep it explicit)
            ys = np.interp(xs, ts, es)
            for xi, yi in zip(xs, ys):
                print(f"  x = {xi:.12e}, err = {yi:.12e}")

        label = rf"$x = \frac{{(1+2\cdot{n})\pi}}{{2}}$"
        plt.plot(t, ew, alpha=0.7, label=label)
        count += 1

    plt.xlim(-w, 0.0)
    plt.legend(loc="upper left", fontsize="small", ncol=2)
    print(f"Plotted {count} pole windows (left side).")
    plt.show()


def problem_area_both(x, err, w_left=1e-5, w_right=1e-4):
    pi = np.pi

    xmin = np.min(x)
    xmax = np.max(x)

    n_min = int(np.floor((2 * xmin / pi - 1) / 2))
    n_max = int(np.ceil((2 * xmax / pi - 1) / 2))

    plt.figure(figsize=(14, 8))
    plt.title(r"tan_simd absolute error near poles $(1+2n)\pi/2$ (centered overlays)")
    plt.xlabel(r"offset from pole $x - (1+2n)\pi/2$")
    plt.ylabel(r"$tan(x) - tan\_simd(x)$")
    plt.grid(True)

    # get default matplotlib color cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, n in enumerate(range(n_min, n_max + 1)):
        pole = (1 + 2*n) * pi / 2.0
        label = rf"$x = \frac{{(1+2\cdot{n})\pi}}{{2}}$"
        color = colors[i % len(colors)]
        first_plot_for_n = True  # one legend entry per n

        # LEFT window: [pole - w_left, pole]
        maskL = (x >= pole - w_left) & (x <= pole)
        if np.any(maskL):
            xw = x[maskL]
            ew = err[maskL]
            t = xw - pole
            finite = np.isfinite(t) & np.isfinite(ew)
            t = t[finite]
            ew = ew[finite]
            if t.size:
                idx = np.argsort(t, kind="mergesort")
                plt.plot(
                    t[idx], ew[idx],
                    alpha=0.7,
                    color=color,
                    label=label if first_plot_for_n else None
                )
                first_plot_for_n = False

        # RIGHT window: [pole, pole + w_right]
        maskR = (x >= pole) & (x <= pole + w_right)
        if np.any(maskR):
            xw = x[maskR]
            ew = err[maskR]
            t = xw - pole
            finite = np.isfinite(t) & np.isfinite(ew)
            t = t[finite]
            ew = ew[finite]
            if t.size:
                idx = np.argsort(t, kind="mergesort")
                plt.plot(
                    t[idx], ew[idx],
                    alpha=0.7,
                    color=color,
                    label=label if first_plot_for_n else None
                )

    plt.xlim(-w_left, w_right)
    plt.legend(loc="upper left", fontsize="small", ncol=2)
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

def ulp_and_abs_error(path_ulp, path_abs):
    x1, y1 = get_data(path_ulp)
    x2, y2 = get_data(path_abs)
    y1 = y1 / 100000


    plt.figure(figsize=(14, 8))
    plt.title(fr"Compare {path_ulp} and {path_abs}")
    plt.xlabel(r"offset from pole $x - (1+2n)\pi/2$")
    plt.ylabel(r"$tan(x) - tan\_simd(x)$")
    plt.grid(True)

    plt.plot(x1, y1, label=f'{path_ulp}')
    plt.plot(x2, y2, label=f'{path_abs}')

    plt.legend()
    plt.show()


if __name__ == "__main__":
    print("Next: ", 3/2 * np.pi - 0.00001, 3/2 * np.pi + 0.00001)
    x, err = get_data('./tan_ulp_error_behavior.tsv')
    # x, err = get_data('./tan_error_behavior.tsv')
    # simple_error_plot(x, err) 
    simple_scatter_error_plot(x, err)
    # plot_range(x, err)
    # problem_area_right(x, err)
    # problem_area_left(x, err)
    # problem_area_both(x, err)
    # compare_correction('one_bit_smaller_correction.tsv', 'tan_error_behavior.tsv')
    # ulp_and_abs_error('./tan_ulp_error_behavior.tsv', './tan_error_behavior.tsv')
