import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


class InteractivePlot:
    def __init__(self, top_degree, bot_degree, xlim=(-2*np.pi, 2*np.pi)) -> None:
        self.top_degree = int(top_degree)      # number of coeff sliders for P(x)
        self.bot_degree = int(bot_degree)      # number of coeff sliders for Q(x)
        self.x = np.linspace(xlim[0], xlim[1], 1200)
        self.xlim = xlim

        # Matplotlib holders
        self.fig = None
        self.ax_top = None
        self.ax_bot = None
        self.line_ratio = None   # R(x) = P/Q
        self.line_tan = None     # tan(x)
        self.line_err = None     # |R - tan|

        # Sliders
        self.sliders_top = []    # a_i for P(x)
        self.sliders_bot = []    # b_i for Q(x)

    # --- polynomial evaluator (Horner) for y = sum_i coeffs[i] * x**i ---
    def eval_sliders(self, x, *coeffs):
        # Horner's method with coeffs in ascending powers (c0 + c1 x + ...)
        y = np.zeros_like(x, dtype=float)
        for a in reversed(coeffs):
            y = y * x + a
        return y

    def plot(self):
        # Figure & axes
        n_sliders = self.top_degree + self.bot_degree
        self.fig = plt.figure(figsize=(10, 7 + 0.32 * max(1, n_sliders)))
        gs = self.fig.add_gridspec(
            nrows=3, ncols=1,
            height_ratios=[3, 1.2, 0.22 * max(1, n_sliders)],  # bigger top plot
            hspace=0.30
        )
        self.ax_top = self.fig.add_subplot(gs[0, 0])
        self.ax_bot = self.fig.add_subplot(gs[1, 0], sharex=self.ax_top)

        self.ax_top.set_xlim(*self.xlim)
        self.ax_top.set_title("Top: R(x) = P(x)/Q(x)  and  tan(x)")
        self.ax_bot.set_title("Bottom: Absolute error  |P(x)/Q(x) − tan(x)|")
        self.ax_bot.set_xlabel("x")
        self.ax_top.set_ylabel("y")
        self.ax_bot.set_ylabel("|error|")

        # --- Initial coefficients ---
        # Top: Taylor (Maclaurin) coefficients for tan(x), even degrees allowed (init 0)
        init_top = self._taylor_coeffs_tan(self.top_degree) if self.top_degree else []
        # Bottom: unchanged (identity denominator): b0=1, others 0
        init_bot = [1.0] + [0.0] * (self.bot_degree - 1) if self.bot_degree else []

        # Initial curves
        P = self.eval_sliders(self.x, *init_top) if self.top_degree else np.zeros_like(self.x)
        Q = self.eval_sliders(self.x, *init_bot) if self.bot_degree else np.ones_like(self.x)
        R = self._safe_ratio(P, Q)
        T = self._safe_tan(self.x)
        E = self._safe_abs_error(R, T)

        (self.line_ratio,) = self.ax_top.plot(self.x, R, lw=2, label="P(x)/Q(x)")
        (self.line_tan,) = self.ax_top.plot(self.x, T, lw=2, linestyle="--", label="tan(x)")
        (self.line_err,) = self.ax_bot.plot(self.x, E, lw=2, label="|R - tan|")

        self.ax_top.legend(loc="upper right", fontsize=9, frameon=False)
        self.ax_bot.legend(loc="upper right", fontsize=9, frameon=False)

        # Sliders container
        sliders_ax = self.fig.add_subplot(gs[2, 0])
        sliders_ax.axis("off")

        # Helper to add sliders (wider, ±5 range)
        def add_slider(y0, label, vinit=0.0, vmin=-5.0, vmax=5.0):
            # wider: move left to 0.06 and width to 0.88
            rect = [0.06, y0, 0.88, 0.065]  # (left, bottom, width, height) in axes fraction
            ax = sliders_ax.inset_axes(rect)
            s = Slider(ax=ax, label=label, valmin=vmin, valmax=vmax, valinit=vinit)
            return s

        pad = 0.082
        ycursor = 0.92

        # Top sliders a0..a_{top_degree-1}
        for i in range(self.top_degree):
            s = add_slider(ycursor, f"a{i}", vinit=init_top[i] if i < len(init_top) else 0.0)
            s.on_changed(self._on_any_change)
            self.sliders_top.append(s)
            ycursor -= pad

        # Bottom sliders b0..b_{bot_degree-1}
        for j in range(self.bot_degree):
            s = add_slider(ycursor, f"b{j}", vinit=init_bot[j] if j < len(init_bot) else 0.0)
            s.on_changed(self._on_any_change)
            self.sliders_bot.append(s)
            ycursor -= pad

        # Reset button
        reset_ax = sliders_ax.inset_axes([0.82, 0.02, 0.15, 0.10])
        reset_btn = Button(reset_ax, "Reset")
        reset_btn.on_clicked(self._on_reset)

        plt.show()

    # --- numerical helpers with safe masking around singularities ---
    def _safe_ratio(self, P, Q, eps=1e-10, clip=None):
        R = np.full_like(P, np.nan, dtype=float)
        mask = np.abs(Q) > eps
        R[mask] = P[mask] / Q[mask]
        if clip is not None:
            R = np.where(np.abs(R) > clip, np.nan, R)
        return R

    def _safe_tan(self, x, eps=1e-10, clip=None):
        c = np.cos(x)
        T = np.full_like(x, np.nan, dtype=float)
        mask = np.abs(c) > eps
        T[mask] = np.tan(x[mask])
        if clip is not None:
            T = np.where(np.abs(T) > clip, np.nan, T)
        return T

    def _safe_abs_error(self, R, T):
        return np.abs(R - T)

    # --- unified callback ---
    def _on_any_change(self, _):
        a = [s.val for s in self.sliders_top] if self.sliders_top else []
        b = [s.val for s in self.sliders_bot] if self.sliders_bot else []

        P = self.eval_sliders(self.x, *a) if a else np.zeros_like(self.x)
        Q = self.eval_sliders(self.x, *b) if b else np.ones_like(self.x)
        R = self._safe_ratio(P, Q)
        T = self._safe_tan(self.x)
        E = self._safe_abs_error(R, T)

        self.line_ratio.set_ydata(R)
        self.line_tan.set_ydata(T)
        self.line_err.set_ydata(E)
        self.fig.canvas.draw_idle()

    def _on_reset(self, _event):
        for s in self.sliders_top + self.sliders_bot:
            s.reset()

    # --- Taylor (Maclaurin) coefficients for tan(x) up to given degree ---
    def _taylor_coeffs_tan(self, max_degree):
        """
        Returns a list [c0, c1, ..., c_{max_degree-1}] such that
        P(x) = sum_i c_i x^i approximates tan(x) at x=0.
        Even degrees are included and set to 0 by default (but user can adjust).
        """
        coeffs = [0.0] * max(0, max_degree)

        # Odd-power coefficients from the series:
        # tan x = x + x^3/3 + 2x^5/15 + 17x^7/315 + 62x^9/2835
        #       + 1382x^11/155925 + 21844x^13/6081075 + 929569x^15/638512875
        #       + 6404582x^17/10854718875 + 443861162x^19/1856156927625 + ...
        odd_coeff = {
            1: 1.0,
            3: 1.0/3.0,
            5: 2.0/15.0,
            7: 17.0/315.0,
            9: 62.0/2835.0,
            11: 1382.0/155925.0,
            13: 21844.0/6081075.0,
            15: 929569.0/638512875.0,
            17: 6404582.0/10854718875.0,
            19: 443861162.0/1856156927625.0,
        }

        for k, v in odd_coeff.items():
            if k < max_degree:
                coeffs[k] = v

        # c0 (constant term) of tan at 0 is 0
        if max_degree >= 1:
            coeffs[0] = 0.0

        return coeffs


if __name__ == "__main__":
    interactive_plot = InteractivePlot(9, 5, (0, 1.5))  # e.g., allow even degrees on top
    interactive_plot.plot()
