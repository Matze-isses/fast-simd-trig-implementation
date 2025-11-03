#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.optimize import basinhopping


class InteractiveApproxApp:
    # ---- Config ----
    X_MIN = np.pi / 4
    X_MAX = 1.57            # ~ just under π/2 to avoid tan blowup
    PLOT_SAMPLES = 4000
    EPS = 1e-4
    TOP_COEFF = 9           # first 9 coefficients are numerator
    BOUND_MIN = -10.0
    BOUND_MAX = 10.0

    def __init__(
        self,
        x_vals_init,
        coefs_init=None,
        figsize=(11.5, 8.5),
        guess_state_path="interactive_plot_state.json",
        result_state_path="optimized_params_state.json",
    ):
        self.result_state_path = result_state_path
        self._init_params(coefs_init, guess_state_path)

        # Store initial x-values (first & last fixed; internals adjustable)
        self.x_vals = np.array(x_vals_init, dtype=float)
        if self.x_vals.ndim != 1 or len(self.x_vals) < 3:
            raise ValueError("x_vals_init must be a 1D array with at least 3 points.")
        self.n_internal = len(self.x_vals) - 2  # number of sliders

        self.last_result = None
        self.running = False

        # Dense grid for plotting
        self.plot_x = np.linspace(self.X_MIN, self.X_MAX, self.PLOT_SAMPLES)
        self.tan_vals = np.tan(self.plot_x)

        # Build UI and render
        self._build_ui(figsize)
        self._update_scatter_and_current()
        plt.show()

    # ---- Parameter init / loading ----
    def _init_params(self, coefs_init, guess_state_path):
        """
        Build current parameter vector:
          [num_coeffs (TOP_COEFF), den_coeffs (m), num_centers (TOP_COEFF-1), den_centers (m-1)]
        """
        if coefs_init is not None:
            # Back-compat: coeffs-only provided; centers -> zeros
            num_len = self.TOP_COEFF
            den_len = max(1, len(coefs_init) - num_len)
            num_coeffs = list(coefs_init[:num_len]) + [0.0] * max(0, num_len - len(coefs_init))
            den_coeffs = list(coefs_init[num_len:num_len + den_len])
            if len(den_coeffs) == 0:
                den_coeffs = [1.0]
            num_centers = [0.0] * max(0, num_len - 1)
            den_centers = [0.0] * max(0, den_len - 1)
        else:
            # Load from guess app state
            try:
                with open(guess_state_path, "r", encoding="utf-8") as f:
                    state = json.load(f)
            except Exception:
                state = {}

            num_coeffs = list(state.get("sliders_top", []))
            den_coeffs = list(state.get("sliders_bot", []))

            # Ensure numerator length = TOP_COEFF
            if len(num_coeffs) < self.TOP_COEFF:
                num_coeffs = num_coeffs + [0.0] * (self.TOP_COEFF - len(num_coeffs))
            elif len(num_coeffs) > self.TOP_COEFF:
                num_coeffs = num_coeffs[:self.TOP_COEFF]

            if len(den_coeffs) == 0:
                den_coeffs = [1.0]

            num_centers = list(state.get("sliders_top_cent", []))
            den_centers = list(state.get("sliders_bot_cent", []))

            # Enforce centers length consistency
            ncl = max(0, self.TOP_COEFF - 1)
            dcl = max(0, len(den_coeffs) - 1)
            if len(num_centers) < ncl:
                num_centers += [0.0] * (ncl - len(num_centers))
            else:
                num_centers = num_centers[:ncl]
            if len(den_centers) < dcl:
                den_centers += [0.0] * (dcl - len(den_centers))
            else:
                den_centers = den_centers[:dcl]

        # Persist lengths
        self.num_len = self.TOP_COEFF
        self.den_len = len(den_coeffs)
        self.num_cent_len = max(0, self.num_len - 1)
        self.den_cent_len = max(0, self.den_len - 1)

        # Build current parameter vector
        self.current_params = np.array(
            list(num_coeffs) + list(den_coeffs) + list(num_centers) + list(den_centers),
            dtype=float
        )

    # ---- Model & objective ----
    def _split_params(self, params):
        """Split [num_coeffs, den_coeffs, num_centers, den_centers]."""
        n1 = self.num_len
        n2 = self.den_len
        c1 = self.num_cent_len
        c2 = self.den_cent_len

        num_coeffs = params[0:n1]
        den_coeffs = params[n1:n1+n2]
        num_centers = params[n1+n2:n1+n2+c1]
        den_centers = params[n1+n2+c1:n1+n2+c1+c2]
        return num_coeffs, den_coeffs, num_centers, den_centers

    @staticmethod
    def _poly_with_centers(x, coeffs, centers):
        """
        y = c0 + sum_{k=1}^{n-1} c_k * (x - centers[k-1])**k
        """
        n = len(coeffs)
        if n == 0:
            return np.zeros_like(x, dtype=float)
        y = np.full_like(x, coeffs[0], dtype=float)
        for k in range(1, n):
            a = centers[k-1] if (k-1) < len(centers) else 0.0
            y = y + coeffs[k] * (x - a) ** k
        return y

    def _mixed_poly_approx(self, x, *params, top_coeff=None):
        """
        Rational approximant with per-term centers:
            P(x) = c0 + Σ_{k=1}^{n-1} c_k (x - a_{k-1})^k
            Q(x) = d0 + Σ_{k=1}^{m-1} d_k (x - b_{k-1})^k
        Returns P(x)/Q(x). The split uses self.num_len / self.den_len.
        """
        if len(params) == 1 and isinstance(params[0], (list, np.ndarray)):
            params = params[0]
        num_coeffs, den_coeffs, num_centers, den_centers = self._split_params(np.array(params, dtype=float))
        P = self._poly_with_centers(x, num_coeffs, num_centers)
        Q = self._poly_with_centers(x, den_coeffs, den_centers)
        R = np.full_like(P, np.nan, dtype=float)
        mask = np.abs(Q) > self.EPS
        R[mask] = P[mask] / Q[mask]
        return R

    def _objective(self, params):
        """L1 error on the CURRENT (possibly slider-adjusted) x points."""
        y_true = np.tan(self.x_vals)
        y_hat = self._mixed_poly_approx(self.x_vals, params)
        return np.sum(np.abs(y_hat - y_true))


    def _get_bounds(self, n_params):
        return [(self.BOUND_MIN, self.BOUND_MAX)] * n_params

    def _run_optimization(self, start_params):
        # Use bounded local minimizer so result.x automatically respects [-10, 10].
        print("Starting New Optimization Process!")
        bounds = self._get_bounds(len(start_params))
        minimizer_kwargs = dict(
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 250_000, "ftol": 1e-11},
        )
        result = basinhopping(
            self._objective,
            np.array(start_params, dtype=float),
            niter=150,           # global hops
            stepsize=0.12,       # hop size before local refine
            T=0.4,               # acceptance temperature
            minimizer_kwargs=minimizer_kwargs,
            disp=False
        )
        # Bounded local search already enforces limits; clip final defensively
        if hasattr(result, "x"):
            result.x = np.clip(result.x, self.BOUND_MIN, self.BOUND_MAX)
       
        return result

    # ---- UI construction ----
    def _build_ui(self, figsize):
        plt.close('all')
        self.fig = plt.figure(figsize=figsize)

        # Two-column layout:
        # Left column: top plot (row 0), bottom plot (row 1)
        # Right column: vertical stack of sliders and the "Run Optimization" button
        gs = self.fig.add_gridspec(
            nrows=2, ncols=2,
            width_ratios=[3.0, 1.3],
            height_ratios=[2.0, 1.0],
            wspace=0.28, hspace=0.35
        )

        # Plots on the left
        self.ax_top = self.fig.add_subplot(gs[0, 0])
        self.ax_bot = self.fig.add_subplot(gs[1, 0], sharex=self.ax_top)

        # --- Top plot: tan, current params curve, optimized (after first run) ---
        (self.line_tan,) = self.ax_top.plot(self.plot_x, self.tan_vals, lw=2, alpha=0.7, label="tan(x)")

        y_cur = self._mixed_poly_approx(self.plot_x, self.current_params)
        (self.line_prior,) = self.ax_top.plot(self.plot_x, y_cur, lw=1.6, label="Current params")

        (self.line_opt,) = self.ax_top.plot(self.plot_x, np.nan * self.plot_x,
                                            lw=1.6, linestyle="--", label="Optimized (last)")

        # Scatter of the selected x points on tan(x)
        (self.scatter_pts,) = self.ax_top.plot(self.x_vals, np.tan(self.x_vals), "o", ms=5, label="x points")

        self.ax_top.set_ylabel("f(x)")
        self.ax_top.set_title("Interactive rational approx to tan(x)")
        self.ax_top.grid(True, which="both", alpha=0.3)
        self.ax_top.legend()

        # >>> Make ONLY the upper plot use a symmetric log scale for better readability near 0
        #     linthresh sets the half-width of the linear region around 0 (here ±1e-3),
        #     so behavior near zero is visible while still handling negatives and large values.
        self.ax_top.set_yscale('symlog', linthresh=1e-3, linscale=1.0, base=10)

        # --- Bottom plot: ONLY error of the latest optimization (linear scale, unchanged) ---
        (self.line_err_opt,) = self.ax_bot.plot(self.plot_x, np.nan * self.plot_x,
                                                lw=1.6, linestyle="--", label="|tan - optimized|")
        self.ax_bot.set_xlabel("x")
        self.ax_bot.set_ylabel("Absolute Error")
        self.ax_bot.grid(True)
        self.ax_bot.legend()

        # Right column container
        right_ax = self.fig.add_subplot(gs[:, 1])  # spans both rows
        right_ax.axis("off")

        # Stack sliders and the button inside the right panel with inset axes
        n = self.n_internal
        pad = 0.012
        slider_h = 0.03   # half-height sliders
        btn_h = 0.08
        top_margin = 0.96
        bottom_margin = 0.06
        total_slider_h = n * slider_h + max(0, n - 1) * pad
        remaining = top_margin - bottom_margin - total_slider_h - btn_h - pad
        if remaining < 0:
            scale = (top_margin - bottom_margin - btn_h - pad) / (total_slider_h + 1e-9)
            slider_h *= max(0.5, scale)
            total_slider_h = n * slider_h + max(0, n - 1) * pad

        y = top_margin
        self.slider_axes = []
        self.sliders = []
        self.internal_idx = list(range(1, len(self.x_vals) - 1))

        for k, i in enumerate(self.internal_idx):
            y0 = y - slider_h
            ax_s = right_ax.inset_axes([0.10, y0, 0.80, slider_h])
            s = Slider(ax_s, f"x[{i}]", valmin=self.X_MIN, valmax=self.X_MAX, valinit=self.x_vals[i])
            s.on_changed(self._on_slider_change)
            self.slider_axes.append(ax_s)
            self.sliders.append(s)
            y = y0 - pad

        # Button at the bottom of the right panel
        ax_btn = right_ax.inset_axes([0.10, bottom_margin, 0.80, btn_h])
        self.btn_opt = Button(ax_btn, "Run Optimization", hovercolor="0.95")
        self.btn_opt.on_clicked(self._on_click)

    # ---- UI callbacks & updates ----
    def _on_slider_change(self, _):
        # Read all slider values and sort them while keeping endpoints fixed
        for s, i in zip(self.sliders, self.internal_idx):
            self.x_vals[i] = s.val

        # Keep monotonic order: sort internal x values, keep endpoints fixed
        internal_sorted = np.sort(self.x_vals[1:-1])
        self.x_vals[1:-1] = internal_sorted

        # Update all sliders to reflect sorted order (without recursive events)
        for s in self.sliders:
            s.eventson = False
        for s, val in zip(self.sliders, internal_sorted):
            s.set_val(val)
        for s in self.sliders:
            s.eventson = True

        self._update_scatter_and_current()

    def _on_click(self, _event):
        if self.running:
            return
        self.running = True
        self.btn_opt.label.set_text("Optimizing...")
        self.fig.canvas.draw_idle()

        result = self._run_optimization(self.current_params)

        if result is not None and hasattr(result, "x"):
            self.last_result = result
            self.current_params = np.array(result.x, dtype=float)

            # ---- Required prints (unchanged) ----
            print("\nOptimization successful:", result.success)
            print("Minimum value:", result.fun, "\n")

            print("Optimized Coefficients:\n")
            for x in result.x[:8]:
                print(f"{x},")
            print(f"{result.x[8]}\n")
            for x in result.x[9:-1]:
                print(f"{x},")
            print(f"{result.x[-1]}")
            print("\n" + "-" * 60 + "\n")

            # Update optimized curve and error
            y_opt = self._mixed_poly_approx(self.plot_x, self.current_params)
            self.line_opt.set_ydata(y_opt)

            err_opt = np.abs(self.tan_vals - y_opt)
            self.line_err_opt.set_ydata(err_opt)

            max_err = np.nanmax(err_opt)
            if np.isfinite(max_err) and max_err > 0:
                self.ax_bot.set_ylim(0, max_err * 1.5)

            test_x = np.linspace(np.pi/4, 1.57, 1_000_000)
            test_y = self._mixed_poly_approx(test_x, self.current_params)
            print("Mean Abs Error: ", np.mean(np.abs(test_y - np.tan(test_x))))

            # Save results to a separate JSON (silent)
            self._save_results_json(self.current_params)

        self.btn_opt.label.set_text("Run Optimization")
        self.running = False
        self.fig.canvas.draw_idle()

    def _save_results_json(self, params):
        num_coeffs, den_coeffs, num_centers, den_centers = self._split_params(params)
        payload = {
            "top_degree": self.num_len,
            "bot_degree": self.den_len,
            "numerator_coeffs": list(map(float, num_coeffs)),
            "denominator_coeffs": list(map(float, den_coeffs)),
            "numerator_centers": list(map(float, num_centers)),
            "denominator_centers": list(map(float, den_centers)),
            "param_vector": list(map(float, params)),
        }
        try:
            with open(self.result_state_path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
        except Exception:
            # Silent by design to avoid altering stdout behavior
            pass

    def _update_scatter_and_current(self):
        # Update scatter of x points
        self.scatter_pts.set_xdata(self.x_vals)
        self.scatter_pts.set_ydata(np.tan(self.x_vals))
        # Update current (prior) curve
        y_prior = self._mixed_poly_approx(self.plot_x, self.current_params)
        self.line_prior.set_ydata(y_prior)
        self.fig.canvas.draw_idle()

    def show(self):
        """Re-show the figure in interactive sessions if needed."""
        plt.show()


# ---- Run as script ----
if __name__ == "__main__":
    # Initial x values (first & last fixed; sliders for the internals)
    x_vals_init = np.linspace(np.pi/4, np.pi/2-0.001, 20)

    # Auto-load from the guess app (coeffs + centers). Prints remain unchanged.
    app = InteractiveApproxApp(x_vals_init=x_vals_init, coefs_init=None)
    # app.show()
