import numpy as np
import matplotlib.pyplot as plt


def get_data(path="./tan_error_behavior.tsv"):
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    x = data[:, 0]
    err = data[:, 1]
    return np.array(x), np.array(err)


def simple_error_plot(x, err):
    plt.figure(figsize=(19.2, 10.8), dpi=100)
    plt.xlabel(r"$x$", fontsize=22)
    plt.ylabel(r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$", fontsize=22)
    plt.grid(True)

    plt.plot(x, err)

    ax = plt.gca()

    ax.ticklabel_format(axis="x", style="plain", useOffset=False)

    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)

    x_left  = x.min()
    x_right = x.max()

    ax.set_xticks([x_left, x_right])
    ax.set_xticklabels([rf"${x_left:.5f}$", r"$\pi/2$"])

    plt.ylim(-10, 0.1)

    plt.savefig("tan_difference_uncorrected_first_singularity.png", dpi=100, bbox_inches="tight")
    plt.show()

def simple_scatter_error_plot(x, err):
    plt.figure(figsize=(19.2, 10.8), dpi=100)

    plt.xlabel("x", fontsize=18)
    plt.ylabel(
            r"$\mathrm{ULP}(\tan_{\mathrm{impl}}(x),\tan_{\mathrm{ref}}(x), x)$",
            fontsize=18
            )

    plt.tick_params(axis="both", which="major", labelsize=14)
    plt.grid(True)

    plt.scatter(x, err, s=1.5, alpha=0.7)
    plt.savefig("ulp_error_thierd_quadrant_corrected.png", dpi=100, bbox_inches="tight")
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
    plt.ylabel(r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$")
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
    plt.ylabel(r"$\tan_{\mathrm{impl}}(x) - \tan_{\mathrm{ref}}(x)$")
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


def problem_area_both(x, err, w_left=1e-7, w_right=1e-7, scatter_plot: bool = False):
    pi = np.pi

    xmin = np.min(x)
    xmax = np.max(x)

    n_min = int(np.floor((2 * xmin / pi - 1) / 2))
    n_max = int(np.ceil((2 * xmax / pi - 1) / 2))

    print(n_min, n_max)

    plt.figure(figsize=(19.2, 10.8), dpi=100)
    plt.xlabel(r"offset from pole", fontsize=20)
    plt.ylabel(r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$", fontsize=20)
    plt.grid(True)

    # get default matplotlib color cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, n in enumerate(range(n_min, n_max + 1)):
        if n != 9:
            continue
        pole = (1 + 2*n) * pi / 2.0
        label = rf"$x = {1+2*n} \cdot pi/2 +" + r"\mathrm{offset}$"
        color = colors[0]
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
                if scatter_plot: 
                    plt.scatter(
                        t[idx], ew[idx],
                        alpha=0.9,
                        s=1.5,
                        color=color,
                        label=label if first_plot_for_n else None
                    )
                else:
                    plt.plot(
                        t[idx], ew[idx],
                        alpha=0.9,
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
                if scatter_plot: 
                    plt.scatter(
                        t[idx], ew[idx],
                        alpha=0.9,
                        s=1.5,
                        color=color,
                        label=label if first_plot_for_n else None
                    )
                else:
                    plt.plot(
                        t[idx], ew[idx],
                        alpha=0.9,
                        color=color,
                        label=label if first_plot_for_n else None
                    )

    plt.xlim(-w_left, w_right)
    plt.ylim(-1e-6, 1e-6)
    plt.legend(loc="lower right", fontsize="small", ncol=2)

    ax = plt.gca()
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)

    ax.xaxis.get_offset_text().set_fontsize(15)
    ax.yaxis.get_offset_text().set_fontsize(15)

    plt.savefig("tan_positive_singularities.png", dpi=100, bbox_inches="tight")
    plt.show()


def problem_area_both_with_ulp(
    x, err, x_ulp, err_ulp,
    w_left=1e-7, w_right=1e-7,
    n_only=9,
    outfile="tan_positive_singularities.png",
):
    pi = np.pi

    xmin = np.min(x)
    xmax = np.max(x)

    n_min = int(np.floor((2 * xmin / pi - 1) / 2))
    n_max = int(np.ceil((2 * xmax / pi - 1) / 2))

    # ---- FULL HD, FULL WIDTH ----
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True,
        figsize=(19.2, 10.8), dpi=100,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.05}
    )

    ax_top.set_ylabel(
        r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$", fontsize=20
    )
    ax_bot.set_xlabel(r"offset from pole", fontsize=20)

    ax_bot.set_ylabel(
        r"$\mathrm{ULP}(\tan_{\mathrm{HP}}, \tan_{\mathrm{approx}}, x)$", fontsize=20
    )

    ax_top.grid(True)
    ax_bot.grid(True)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    def _plot_window(ax, xdat, ydat, pole, *, color, label=None, scatter=False):
        mask = (xdat >= pole - w_left) & (xdat <= pole + w_right)
        t = xdat[mask] - pole
        y = ydat[mask]
        finite = np.isfinite(t) & np.isfinite(y)
        t = t[finite]
        y = y[finite]

        if t.size:
            idx = np.argsort(t, kind="mergesort")
            if scatter:
                ax.scatter(t[idx], y[idx], s=6, alpha=0.9, color=color)
            else:
                ax.plot(t[idx], y[idx], alpha=0.9, color=color, label=label)

    for i, n in enumerate(range(n_min, n_max + 1)):
        if (n_only is not None) and (n != n_only):
            continue

        pole = (1 + 2 * n) * pi / 2.0
        label = rf"$x = {1+2*n} \cdot \pi/2 + \mathrm{{offset}}$"
        color = colors[i % len(colors)]

        _plot_window(ax_top, x, err, pole, color=color, label=label, scatter=False)
        _plot_window(ax_bot, x_ulp, err_ulp, pole, color=color, scatter=True)

    # ---- limits ----
    ax_bot.set_xlim(-w_left, w_right)
    ax_top.set_ylim(-1e-6, 1e-6)

    # ---- ULP ticks: force ±3 ----
    ax_bot.set_yticks([-3, -2, -1, 0, 1, 2, 3])
    ax_bot.set_ylim(-3.2, 3.2)

    # ---- ticks + scientific offset font sizes ----
    for ax in (ax_top, ax_bot):
        ax.tick_params(axis="both", labelsize=15)
        ax.xaxis.get_offset_text().set_fontsize(15)
        ax.yaxis.get_offset_text().set_fontsize(15)

    ax_top.legend(
        loc="upper right",
        ncol=2,
        fontsize=16,      # ← increase this
    )

    fig.savefig(outfile, dpi=100, bbox_inches="tight")
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

def scatter_err_near_tan_poles(x, err, outfile="tan_err_near_poles.png", isulp=False):
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    if x.shape != err.shape:
        raise ValueError("x and err must have the same shape")

    m = np.isfinite(x) & np.isfinite(err)
    x = x[m]
    err = err[m]
    if x.size == 0:
        raise ValueError("No finite (x, err) pairs")

    pi = np.pi

    # distance to nearest tan pole: reduce x modulo π into [0, π)
    r = np.mod(x, pi)                 # in [0, π)
    d = np.abs(r - pi/2.0)            # distance to π/2 within the cell
    keep = d <= (pi/8.0)

    xk = x[keep]
    ek = err[keep]

    fig, ax = plt.subplots(figsize=(19.2, 10.8), dpi=100)
    ax.scatter(xk, ek, s=0.1, alpha=0.7)

    ax.set_xlabel("x", fontsize=18)

    if isulp:
        ax.set_ylabel(
            r"$\mathrm{ULP}(\tan_{\mathrm{HP}}, \tan_{\mathrm{approx}}, x)$", fontsize=18
        )
    else:
        ax.set_ylabel(r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$", fontsize=22)

    ax.grid(True)
    ax.tick_params(axis="both", labelsize=14)

    fig.savefig(outfile, dpi=100, bbox_inches="tight")
    plt.show()

def scatter_err_tan_poles_annulus(x, err, outfile="tan_err_pole_annulus.png", isulp=False):
    """
    Scatter-plot (x, err) only for points whose distance to the nearest tan pole is:
        π/8 <= d <= π/4
    where poles are at x = (2n+1)*π/2.

    Inputs: x, err, outfile
    """
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    if x.shape != err.shape:
        raise ValueError("x and err must have the same shape")

    m = np.isfinite(x) & np.isfinite(err)
    x = x[m]
    err = err[m]
    if x.size == 0:
        raise ValueError("No finite (x, err) pairs")

    pi = np.pi

    # distance to nearest tan pole using reduction modulo π
    r = np.mod(x, pi)                 # in [0, π)
    d = np.abs(r - pi/2.0)            # distance to π/2 within the cell

    keep = (d >= (pi/8.0)) & (d <= (pi/4.0))

    xk = x[keep]
    ek = err[keep]

    print(f"MAX EK: {np.max(np.abs(ek))}")

    fig, ax = plt.subplots(figsize=(19.2, 10.8), dpi=100)
    ax.scatter(xk, ek, s=0.1, alpha=0.7)

    ax.set_xlabel("x", fontsize=18)
    if isulp:
        ax.set_ylabel(
            r"$\mathrm{ULP}(\tan_{\mathrm{HP}}, \tan_{\mathrm{approx}}, x)$", fontsize=18
        )
    else:
        ax.set_ylabel(r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$", fontsize=22)

    ax.tick_params(axis="both", labelsize=15)
    ax.xaxis.get_offset_text().set_fontsize(15)
    ax.yaxis.get_offset_text().set_fontsize(15)

    ax.grid(True)
    ax.tick_params(axis="both", labelsize=14)

    fig.savefig(outfile, dpi=100, bbox_inches="tight")
    plt.show()


def scatter_err_tan_poles_outer_annulus(x, err, outfile="tan_err_pole_outer_annulus_small_range.png", isulp=False):
    """
    Scatter-plot (x, err) only for points whose distance to the nearest tan pole is:
        π/4 <= d <= 3π/8
    where poles are at x = (2n+1)*π/2.

    Inputs: x, err, outfile
    """
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    if x.shape != err.shape:
        raise ValueError("x and err must have the same shape")

    m = np.isfinite(x) & np.isfinite(err)
    x = x[m]
    err = err[m]
    if x.size == 0:
        raise ValueError("No finite (x, err) pairs")

    pi = np.pi

    # distance to nearest tan pole using reduction modulo π
    r = np.mod(x, pi)                 # in [0, π)
    d = np.abs(r - pi/2.0)            # distance to π/2 within the cell

    keep = (d >= (pi/4.0)) & (d <= (3.0*pi/8.0))

    xk = x[keep]
    ek = err[keep]

    fig, ax = plt.subplots(figsize=(19.2, 10.8), dpi=100)
    ax.scatter(xk, ek, s=0.1, alpha=0.7)

    ax.set_xlabel("x", fontsize=18)
    if isulp:
        ax.set_ylabel(
            r"$\mathrm{ULP}(\tan_{\mathrm{HP}}, \tan_{\mathrm{approx}}, x)$", fontsize=18
        )
    else:
        ax.set_ylabel(r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$", fontsize=22)

    ax.tick_params(axis="both", labelsize=15)
    ax.xaxis.get_offset_text().set_fontsize(15)
    ax.yaxis.get_offset_text().set_fontsize(15)

    ax.grid(True)
    ax.tick_params(axis="both", labelsize=14)

    fig.savefig(outfile, dpi=100, bbox_inches="tight")
    plt.show()


def scatter_err_tan_poles_far_region(x, err, outfile="tan_err_pole_far_region_small_range.png", isulp=False):
    """
    Scatter-plot (x, err) only for points whose distance to the nearest tan pole is:
        3π/8 < d <= π/2
    where poles are at x = (2n+1)*π/2.

    Inputs: x, err, outfile
    """
    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)
    if x.shape != err.shape:
        raise ValueError("x and err must have the same shape")

    m = np.isfinite(x) & np.isfinite(err)
    x = x[m]
    err = err[m]
    if x.size == 0:
        raise ValueError("No finite (x, err) pairs")

    pi = np.pi

    # distance to nearest tan pole using reduction modulo π
    r = np.mod(x, pi)                 # in [0, π)
    d = np.abs(r - pi/2.0)            # distance to π/2 within the cell

    keep = (d > (3.0*pi/8.0)) & (d <= (pi/2.0))

    xk = x[keep]
    ek = err[keep]

    fig, ax = plt.subplots(figsize=(19.2, 10.8), dpi=100)
    ax.scatter(xk, ek, s=0.1, alpha=0.7)

    ax.set_xlabel("x", fontsize=18)

    if isulp:
        ax.set_ylabel(
            r"$\mathrm{ULP}(\tan_{\mathrm{HP}}, \tan_{\mathrm{approx}}, x)$", fontsize=18
        )
    else:
        ax.set_ylabel(r"$\tan_{\mathrm{HP}}(x) - \tan_{\mathrm{approx}}(x)$", fontsize=22)

    ax.tick_params(axis="both", labelsize=15)
    ax.xaxis.get_offset_text().set_fontsize(15)
    ax.yaxis.get_offset_text().set_fontsize(15)

    ax.grid(True)
    ax.tick_params(axis="both", labelsize=14)

    fig.savefig(outfile, dpi=100, bbox_inches="tight")
    plt.show()

def fit_linear_constants_poleband(x, err, dmin, dmax):
    x = np.asarray(x, dtype=float).ravel()
    err = np.asarray(err, dtype=float).ravel()
    if x.shape != err.shape:
        raise ValueError("x and err must have the same shape")

    if not (np.isfinite(dmin) and np.isfinite(dmax) and dmin <= dmax):
        raise ValueError("Require finite dmin <= dmax")

    # keep finite pairs
    m = np.isfinite(x) & np.isfinite(err)
    x = x[m]
    err = err[m]
    if x.size == 0:
        raise ValueError("No finite (x, err) pairs")

    pi = np.pi

    # distance to nearest tan pole via reduction modulo pi
    r = np.mod(x, pi)          # in [0, pi)
    d = np.abs(r - pi/2.0)     # in [0, pi/2]

    # apply pole-distance band
    keep = (d >= dmin) & (d <= dmax)
    xk = x[keep]
    ek = err[keep]
    if xk.size < 2:
        raise ValueError("Not enough points after filtering (need >= 2)")

    # OLS fit: ek ≈ a*xk + b
    A = np.column_stack([xk, np.ones_like(xk)])
    (a, b), *_ = np.linalg.lstsq(A, ek, rcond=None)
    return a, b



if __name__ == "__main__":
    # x, err = get_data('./tan_ulp_error_behavior.tsv')
    # x, err = get_data('./tan_error_behavior.tsv')

    # simple_error_plot(x, err) 
    # simple_scatter_error_plot(x, err)

    # scatter_err_near_tan_poles(x, err)

#   x, err = get_data()
#   a, b = fit_linear_constants_poleband(x, err, 3*np.pi/8, np.pi/2)
#   print(a)


#   name = "large_uniform"
#   x, err = get_data('./tan_error_behavior.tsv')

#   scatter_err_tan_poles_far_region(x, err, f"error_first_range_{name}_range.png", False)
#   scatter_err_tan_poles_outer_annulus(x, err, f"error_second_range_{name}_range.png", False)
#   scatter_err_tan_poles_annulus(x, err, f"error_thierd_range_{name}_range.png", False)
    
    x, err = get_data('./tan_ulp_error_behavior.tsv')

    x = np.asarray(x, dtype=float)
    err = np.asarray(err, dtype=float)

    pi = np.pi
    r = np.mod(x, pi)
    d = np.abs(r - pi/2.0)

    keep = [
            ("fourth", (d > 0) & (d <= pi/8.0)),
            ("third", (d > pi/8.0) & (d <= pi/4.0)),
            ("second", (d > pi/4.0) & (d <= (3.0*pi/8.0))),
            ("first", (d > (3.0*pi/8.0)) & (d <= (pi/2.0)))
    ]

    for name, mask in reversed(keep):
        xk = x[mask]
        ek = err[mask]
        i = np.argmax(np.abs(ek))
        print(f"{np.abs(ek[i]/xk[i])} for {name} with {ek[i]}")
    
#   scatter_err_tan_poles_outer_annulus(x, err, f"error_second_range_{name}_range_ulp.png", True)
#   scatter_err_tan_poles_annulus(x, err, f"error_thierd_range_{name}_range_ulp.png", True)
#   scatter_err_near_tan_poles(x, err, f"error_fourth_range_{name}_range_ulp.png", True)

    # scatter_err_tan_poles_outer_annulus(x, err, "error_second_range_small_range_ulp.png")
    # scatter_err_tan_poles_outer_annulus(x, err, "error_second_range_large_range_ulp.png")

    # scatter_err_tan_poles_far_region(x, err, "error_first_range_small_range.png")

    # plot_range(x, err)
    # problem_area_right(x, err)
    # problem_area_left(x, err)
    # problem_area_both(x, err)

    # problem_area_both_with_ulp(x, err, *get_data('./tan_ulp_error_behavior.tsv'), n_only=9)
    # ulp_error_by_singularity(*get_data('./tan_ulp_error_behavior.tsv'), "test_plot.png")
    # compare_correction('one_bit_smaller_correction.tsv', 'tan_error_behavior.tsv')
    # ulp_and_abs_error('./tan_ulp_error_behavior.tsv', './tan_error_behavior.tsv')
