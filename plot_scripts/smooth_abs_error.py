#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os  # <<< added
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import basinhopping


class AdaptiveApproxLoop:
    # ---- Config ----
    X_MIN = np.pi / 4
    X_MAX = np.pi / 2 - 0.0001                 # ~ just under π/2 to avoid tan blowup
    PLOT_SAMPLES = 7000          # dense grid for plotting
    COARSE_SAMPLES = 100_000     # grid to pick worst-error x
    EPS = 1e-4

    TOP_COEFF = 9

    BOUND_MIN = -10.0
    BOUND_MAX = 10.0
    TARGET_N = 10000             # stop when len(x_vals) hits this
    verbose = False

    def __init__(
        self,
        x_vals_init,
        coefs_init=None,
        guess_state_path="interactive_plot_state.json",
        result_state_path="adaptive_results_state.json",
        best_result_path="adaptive_approx_results.json",
        figsize=(12, 10),
    ):
        self.result_state_path = result_state_path
        self.best_result_path = best_result_path
        self.best_max_abs_err = np.inf
        self.mean_abs_error_history = []  # for third plot

        # Load/initialize parameters
        self._init_params(coefs_init, guess_state_path)

        # Initialize adaptive points
        self.x_vals = np.array(x_vals_init, dtype=float)
        self.x_vals = np.clip(np.unique(self.x_vals), self.X_MIN, self.X_MAX)

        # Dense plot grid
        self.plot_x = np.linspace(self.X_MIN, self.X_MAX, self.PLOT_SAMPLES)
        self.tan_vals = np.tan(self.plot_x)

        # Build figure
        self._build_figure(figsize)

        # <<< added
        self._plot_update_count = -1
        self._saved_plot_count = 0
        self._plot_dir = "./plots"
        # <<< end added

        # Run adaptive optimization loop
        self._adaptive_optimize()
        plt.show()

    # ---- Parameter init ----
    def _init_params(self, coefs_init, guess_state_path):
        if coefs_init is not None:
            num_len = self.TOP_COEFF
            den_len = max(1, len(coefs_init) - num_len)
            num_coeffs = list(coefs_init[:num_len]) + [0.0] * max(0, num_len - len(coefs_init))
            den_coeffs = list(coefs_init[num_len:num_len + den_len]) or [1.0]
            num_centers = [0.0] * max(0, num_len - 1)
            den_centers = [0.0] * max(0, den_len - 1)
        else:
            try:
                with open(guess_state_path, "r", encoding="utf-8") as f:
                    state = json.load(f)
            except Exception:
                state = {}
            num_coeffs = list(state.get("sliders_top", []))
            den_coeffs = list(state.get("sliders_bot", []))
            num_centers = list(state.get("sliders_top_cent", []))
            den_centers = list(state.get("sliders_bot_cent", []))

            if len(num_coeffs) < self.TOP_COEFF:
                num_coeffs += [0.0] * (self.TOP_COEFF - len(num_coeffs))
            elif len(num_coeffs) > self.TOP_COEFF:
                num_coeffs = num_coeffs[:self.TOP_COEFF]
            if len(den_coeffs) == 0:
                den_coeffs = [1.0]

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

        self.num_len = self.TOP_COEFF
        self.den_len = len(den_coeffs)
        self.num_cent_len = max(0, self.num_len - 1)
        self.den_cent_len = max(0, self.den_len - 1)

        self.current_params = np.array(
            list(num_coeffs) + list(den_coeffs) + list(num_centers) + list(den_centers),
            dtype=float
        )

    # ---- Model & objective ----
    def _split_params(self, params):
        n1, n2 = self.num_len, self.den_len
        c1, c2 = self.num_cent_len, self.den_cent_len
        num_coeffs = params[0:n1]
        den_coeffs = params[n1:n1 + n2]
        num_centers = params[n1 + n2:n1 + n2 + c1]
        den_centers = params[n1 + n2 + c1:n1 + n2 + c1 + c2]
        return num_coeffs, den_coeffs, num_centers, den_centers

    @staticmethod
    def _poly_with_centers(x, coeffs, centers):
        n = len(coeffs)
        if n == 0:
            return np.zeros_like(x)
        y = np.full_like(x, coeffs[0], dtype=float)
        for k in range(1, n):
            a = centers[k - 1] if (k - 1) < len(centers) else 0.0
            y += coeffs[k] * (x - a) ** k
        return y

    def _mixed_poly_approx(self, x, *params):
        if len(params) == 1 and isinstance(params[0], (list, np.ndarray)):
            params = params[0]
        num_coeffs, den_coeffs, num_centers, den_centers = self._split_params(np.array(params, float))
        P = self._poly_with_centers(x, num_coeffs, num_centers)
        Q = self._poly_with_centers(x, den_coeffs, den_centers)
        mask = np.abs(Q) > self.EPS
        R = np.full_like(P, np.nan)
        R[mask] = P[mask] / Q[mask]
        return R

    def _objective(self, params):
        y_true = np.tan(self.x_vals)
        y_hat = self._mixed_poly_approx(self.x_vals, params)
        return np.sum(np.abs(y_hat - y_true))

    # ---- Optimizer ----
    def _get_bounds(self, n_params):
        return [(self.BOUND_MIN, self.BOUND_MAX)] * n_params

    def _run_optimization(self, start_params):
        print("Starting New Optimization Process!")
        result = basinhopping(
            self._objective,
            np.array(start_params, dtype=float),
            niter=150,
            stepsize=0.12,
            T=0.4,
            minimizer_kwargs=dict(
                method="L-BFGS-B",
                bounds=self._get_bounds(len(start_params)),
                options={"maxiter": 250_000, "ftol": 1e-11},
            ),
            disp=False,
        )
        if hasattr(result, "x"):
            result.x = np.clip(result.x, self.BOUND_MIN, self.BOUND_MAX)
        return result

    # ---- Figure ----
    def _build_figure(self, figsize):
        plt.close('all')
        self.fig = plt.figure(figsize=figsize)
        gs = self.fig.add_gridspec(
            nrows=3, ncols=1,
            height_ratios=[2.0, 1.0, 0.6],
            hspace=0.35
        )

        # (1) Function plot
        self.ax_top = self.fig.add_subplot(gs[0, 0])
        (self.line_tan,) = self.ax_top.plot([], [], lw=2, alpha=0.7, label="tan(x)")
        (self.line_last,) = self.ax_top.plot([], [], lw=1.3, ls="--", label="Last")
        (self.line_curr,) = self.ax_top.plot([], [], lw=1.6, label="Current")
        (self.scatter_pts,) = self.ax_top.plot([], [], "o", ms=4, label="x points")
        self.ax_top.set_ylabel("f(x)")
        self.ax_top.set_yscale('symlog', linthresh=1e-3, base=10)
        self.ax_top.legend()
        self.ax_top.grid(True)

        # (2) Error plot
        self.ax_mid = self.fig.add_subplot(gs[1, 0])
        (self.line_err,) = self.ax_mid.plot([], [], lw=1.6, label="|tan - current|")
        self.ax_mid.set_ylabel("Abs Error")
        self.ax_mid.legend()
        self.ax_mid.grid(True)

        # (3) MAE evolution
        self.ax_bot = self.fig.add_subplot(gs[2, 0])
        (self.line_mae,) = self.ax_bot.plot([], [], lw=1.8, label="Mean Abs Error Evolution")
        self.ax_bot.set_xlabel("Iteration")
        self.ax_bot.set_yscale('symlog', linthresh=1e-3, base=10)
        self.ax_bot.legend()
        self.ax_bot.grid(True)

    def _update_plots(self, prev_params, curr_params):
        self.line_tan.set_data(self.plot_x, self.tan_vals)
        if prev_params is not None:
            self.line_last.set_data(self.plot_x, self._mixed_poly_approx(self.plot_x, prev_params))
        else:
            self.line_last.set_data(self.plot_x, np.nan * self.plot_x)
        y_curr = self._mixed_poly_approx(self.plot_x, curr_params)
        self.line_curr.set_data(self.plot_x, y_curr)

        err = np.abs(self.tan_vals - y_curr)
        self.line_err.set_data(self.plot_x, err)
        self.ax_mid.set_ylim(0, np.nanmax(err) * 1.5)

        self.scatter_pts.set_data(self.x_vals, np.tan(self.x_vals))

        # update mean abs error plot
        x_it = np.arange(1, len(self.mean_abs_error_history) + 1)
        self.line_mae.set_data(x_it, self.mean_abs_error_history)
        if self.mean_abs_error_history:
            self.ax_bot.set_xlim(1, len(self.mean_abs_error_history) + 1)
            self.ax_bot.set_ylim(0, max(self.mean_abs_error_history) * 1.1)

        self.ax_top.relim(); self.ax_top.autoscale_view()
        self.ax_mid.relim(); self.ax_mid.autoscale_view()
        self.fig.canvas.draw_idle()
        plt.pause(1)

        # --- Save every 10th plot ---  <<< added
        self._plot_update_count += 1
        if self._plot_update_count % 10 == 0:
            os.makedirs(self._plot_dir, exist_ok=True)
            self._saved_plot_count += 1
            fname = f"adaptive_approx_loop{self._saved_plot_count:02d}.png"
            fpath = os.path.join(self._plot_dir, fname)
            self.fig.savefig(fpath, dpi=200, bbox_inches="tight")
            if self.verbose:
                print(f"Saved plot #{self._saved_plot_count} to {fpath}")
        # <<< end added

    # ---- Adaptive loop ----
    def _adaptive_optimize(self):
        prev_params = None
        curr_params = np.array(self.current_params, float)
        coarse_x = np.linspace(self.X_MIN, self.X_MAX, self.COARSE_SAMPLES)

        while len(self.x_vals) < self.TARGET_N:
            current_iter = self.x_vals.shape[0]
            minor_change = np.random.uniform(-1/current_iter, 1/current_iter, curr_params.shape[0])
            opt_params = curr_params + minor_change
            result = self._run_optimization(opt_params)

            if result is None or not hasattr(result, "x"):
                print("Optimization failed.")
                break
            prev_params = curr_params.copy()
            curr_params = np.array(result.x, float)

            y_pred = self._mixed_poly_approx(self.plot_x, curr_params)
            mae = np.mean(np.abs(y_pred - self.tan_vals))
            self.mean_abs_error_history.append(mae)

            print(f"\nIteration {len(self.x_vals)} — MAE: {mae:.6e}")

            if self.verbose:
                print(f"Params {curr_params}")
                print(f"xvals  {self.x_vals}")

            tan_c = np.tan(coarse_x)
            y_c = self._mixed_poly_approx(coarse_x, curr_params)
            err_c = np.abs(tan_c - y_c)
            max_err = float(np.nanmax(err_c))
            idx = int(np.nanargmax(err_c))

            if max_err < self.best_max_abs_err:
                self.best_max_abs_err = max_err
                self._save_best_snapshot(curr_params, max_err)

            # >>> Choose the nearest unused grid index
            idx = self._nearest_unseen_index(coarse_x, idx, self.x_vals)
            if idx is None:
                print("All grid points used; stopping.")
                break

            worst_x = float(coarse_x[idx])
            self.x_vals = np.clip(np.unique(np.append(self.x_vals, worst_x)), self.X_MIN, self.X_MAX)

            self._update_plots(prev_params, curr_params)
            self._silent_save(curr_params)

        self.current_params = curr_params

    def _nearest_unseen_index(self, grid, start_idx, existing_points, round_decimals=12):
        if grid.size == 0:
            return None
        existing = set(np.round(existing_points, round_decimals))

        def used(i):
            return np.round(grid[i], round_decimals) in existing

        n = grid.size
        if not used(start_idx):
            return start_idx

        radius = 1
        while (start_idx - radius) >= 0 or (start_idx + radius) < n:
            r = start_idx + radius
            if r < n and not used(r):
                return r
            l = start_idx - radius
            if l >= 0 and not used(l):
                return l
            radius += 1
        return None

    def _silent_save(self, params):
        num_coeffs, den_coeffs, num_centers, den_centers = self._split_params(params)
        payload = {
            "top_degree": self.num_len,
            "bot_degree": self.den_len,
            "numerator_coeffs": list(map(float, num_coeffs)),
            "denominator_coeffs": list(map(float, den_coeffs)),
            "numerator_centers": list(map(float, num_centers)),
            "denominator_centers": list(map(float, den_centers)),
            "param_vector": list(map(float, params)),
            "x_vals": list(map(float, self.x_vals)),
        }
        try:
            with open(self.result_state_path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
        except Exception:
            pass

    def _save_best_snapshot(self, params, max_abs_error):
        num_coeffs, den_coeffs, num_centers, den_centers = self._split_params(params)
        payload = {
            "top_degree": self.num_len,
            "bot_degree": self.den_len,
            "numerator_coeffs": list(map(float, num_coeffs)),
            "denominator_coeffs": list(map(float, den_coeffs)),
            "numerator_centers": list(map(float, num_centers)),
            "denominator_centers": list(map(float, den_centers)),
            "param_vector": list(map(float, params)),
            "x_vals": list(map(float, self.x_vals)),
            "max_abs_error": float(max_abs_error),
        }
        try:
            with open(self.best_result_path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
        except Exception:
            pass


# ---- Run as script ----
if __name__ == "__main__":
    x_vals_init = np.linspace(np.pi/4, np.pi/2 - 0.0001, 5)
    app = AdaptiveApproxLoop(x_vals_init=x_vals_init, coefs_init=None)
