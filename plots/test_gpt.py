
import json
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, SpanSelector  # <-- kept
# Note: No seaborn; pure matplotlib

class InteractivePlot:
    def __init__(self, top_degree, bot_degree, xlim=(-2*np.pi, 2*np.pi)) -> None:
        self.top_degree = int(top_degree)   # number of coefficients a0..a_{n-1}
        self.bot_degree = int(bot_degree)   # number of coefficients b0..b_{m-1}

        self.x = np.linspace(xlim[0], xlim[1], 10000)
        self.xlim = xlim

        # Matplotlib objects
        self.fig = None
        self.ax_top = None
        self.ax_bot = None
        self.line_ratio = None
        self.line_tan = None
        self.line_err = None

        # Sliders
        self.sliders_top = []          # coefficients for numerator P
        self.sliders_bot = []          # coefficients for denominator Q
        self.sliders_top_cent = []     # centers for numerator terms (length top_degree-1)
        self.sliders_bot_cent = []     # centers for denominator terms (length bot_degree-1)

        # State
        self.state_path = "interactive_plot_state.json"
        self._suspend_update = False
        self._auto_ylims = True

        # Zoom state
        self._span = None
        self._zoom_active = False
        self._orig_xlim = None
        self._orig_ylim_top = None
        self._orig_ylim_bot = None

    # --- polynomial evaluator with per-term centers ---
    # Implements: y = c[0] + sum_{k=1}^{n-1} c[k] * (x - centers[k-1])**k
    # centers length must be n-1; centers[k-1] used for power k
    def eval_with_centers(self, x, coeffs, centers):
        n = len(coeffs)
        y = np.zeros_like(x, dtype=float)
        if n == 0:
            return y
        y = y + coeffs[0]
        for k in range(1, n):
            # centers index k-1
            a = centers[k-1] if (k-1) < len(centers) else 0.0
            y = y + coeffs[k] * (x - a) ** k
        return y

    def plot(self):
        # Count sliders for layout
        n_coeff_top = self.top_degree
        n_coeff_bot = self.bot_degree
        n_cent_top  = max(0, self.top_degree - 1)
        n_cent_bot  = max(0, self.bot_degree - 1)

        # Height based on the max stack in any slider column
        max_stack_per_cell = max(n_coeff_top, n_coeff_bot, n_cent_top, n_cent_bot, 1)
        self.fig = plt.figure(figsize=(12, 7 + 0.32 * max_stack_per_cell))

        # 3 rows: plots, plots, controls; 3 cols: plots | coeff sliders | center sliders
        gs = self.fig.add_gridspec(
            nrows=3, ncols=3,
            height_ratios=[4, 3, 0.9],   # third row for controls
            width_ratios=[2.1, 1.25, 1.25],
            hspace=0.30, wspace=0.28
        )

        # Plots (col 0)
        self.ax_top = self.fig.add_subplot(gs[0, 0])
        self.ax_bot = self.fig.add_subplot(gs[1, 0], sharex=self.ax_top)

        self.ax_top.set_xlim(*self.xlim)
        self.ax_top.set_title("Top: R(x) = P(x)/Q(x)  and  tan(x)")
        self.ax_bot.set_title("Bottom: Absolute error  |P(x)/Q(x) − tan(x)|")
        self.ax_bot.set_xlabel("x")
        self.ax_top.set_ylabel("y")
        self.ax_bot.set_ylabel("|error|")

        # Initial coefficients
        init_top = self._taylor_coeffs_tan(self.top_degree) if self.top_degree else []
        init_bot = [1.0] + [0.0] * (self.bot_degree - 1) if self.bot_degree else []

        # Initial centers: zeros (length degree-1)
        init_top_cent = [0.0] * max(0, self.top_degree - 1)
        init_bot_cent = [0.0] * max(0, self.bot_degree - 1)

        # Initial curves using center-aware evaluator
        P = self.eval_with_centers(self.x, init_top, init_top_cent) if self.top_degree else np.zeros_like(self.x)
        Q = self.eval_with_centers(self.x, init_bot, init_bot_cent) if self.bot_degree else np.ones_like(self.x)
        R = self._safe_ratio(P, Q)
        T = self._safe_tan(self.x)
        E = self._safe_abs_error(R, T)

        (self.line_ratio,) = self.ax_top.plot(self.x, R, lw=1, label="P(x)/Q(x)")
        (self.line_tan,)  = self.ax_top.plot(self.x, T, lw=2, alpha=0.6, linestyle="-", label="tan(x)")
        (self.line_err,)  = self.ax_bot.plot(self.x, E, lw=2, label="|R - tan|")

        self.ax_top.legend(loc="upper right", fontsize=9, frameon=False)
        self.ax_bot.legend(loc="upper right", fontsize=9, frameon=False)

        # Save original limits for "Home View"
        self._orig_xlim = tuple(self.ax_top.get_xlim())
        self._orig_ylim_top = tuple(self.ax_top.get_ylim())
        self._orig_ylim_bot = tuple(self.ax_bot.get_ylim())

        # Autoscale bottom to error initially
        self._autoscale_error_axis(E, xlim=self._orig_xlim)

        # --- Slider parent axes ---
        # Column 1: coefficient sliders (row 0: numerator, row 1: denominator)
        coeff_top_ax_parent = self.fig.add_subplot(gs[0, 1])
        coeff_top_ax_parent.set_title("P Coefficients", fontsize=10)
        coeff_top_ax_parent.axis("off")

        coeff_bot_ax_parent = self.fig.add_subplot(gs[1, 1])
        coeff_bot_ax_parent.set_title("Q Coefficients", fontsize=10)
        coeff_bot_ax_parent.axis("off")

        # Column 2: center sliders (row 0: numerator, row 1: denominator)
        cent_top_ax_parent = self.fig.add_subplot(gs[0, 2])
        cent_top_ax_parent.set_title("P Centers", fontsize=10)
        cent_top_ax_parent.axis("off")

        cent_bot_ax_parent = self.fig.add_subplot(gs[1, 2])
        cent_bot_ax_parent.set_title("Q Centers", fontsize=10)
        cent_bot_ax_parent.axis("off")

        # Slider helpers
        def degree_range(i: int):
            # same ranges you used for coefficients; reuse for centers too
            r = 2.0 / (i + 1.0)
            return -r, r

        def add_slider(parent_ax, y0, label, vinit, vmin, vmax):
            rect = [0.06, y0, 0.88, 0.065]
            ax = parent_ax.inset_axes(rect)
            vinit = float(np.clip(vinit, vmin, vmax))
            s = Slider(ax=ax, label=label, valmin=vmin, valmax=vmax, valinit=vinit, valstep=1e-5)
            return s

        pad = 0.082

        # --- Top (P) coefficient sliders in (0,1) ---
        ycursor = 0.92
        for i in range(self.top_degree):
            vmin, vmax = degree_range(i)
            vinit = init_top[i] if i < len(init_top) else 0.0
            s = add_slider(coeff_top_ax_parent, ycursor, f"a{i}", vinit=vinit, vmin=vmin, vmax=vmax)
            s.on_changed(self._on_any_change)
            self.sliders_top.append(s)
            ycursor -= pad

        # --- Bottom (Q) coefficient sliders in (1,1) ---
        ycursor = 0.92
        for j in range(self.bot_degree):
            vmin, vmax = degree_range(j)
            vinit = init_bot[j] if j < len(init_bot) else 0.0
            s = add_slider(coeff_bot_ax_parent, ycursor, f"b{j}", vinit=vinit, vmin=vmin, vmax=vmax)
            s.on_changed(self._on_any_change)
            self.sliders_bot.append(s)
            ycursor -= pad

        # --- Top (P) center sliders in (0,2) ---
        # centers indices correspond to powers k=1..n-1, we label as a_c0 for k=1, a_c1 for k=2, ...
        ycursor = 0.92
        for i in range(max(0, self.top_degree - 1)):
            # reuse degree_range with index i for a practical bound
            vmin, vmax = degree_range(i)
            vinit = init_top_cent[i]
            s = add_slider(cent_top_ax_parent, ycursor, f"a_c{i}", vinit=vinit, vmin=vmin, vmax=vmax)
            s.on_changed(self._on_any_change)
            self.sliders_top_cent.append(s)
            ycursor -= pad

        # --- Bottom (Q) center sliders in (1,2) ---
        ycursor = 0.92
        for j in range(max(0, self.bot_degree - 1)):
            vmin, vmax = degree_range(j)
            vinit = init_bot_cent[j]
            s = add_slider(cent_bot_ax_parent, ycursor, f"b_c{j}", vinit=vinit, vmin=vmin, vmax=vmax)
            s.on_changed(self._on_any_change)
            self.sliders_bot_cent.append(s)
            ycursor -= pad

        # --- Controls row (span all columns) ---
        controls_ax = self.fig.add_subplot(gs[2, :])
        controls_ax.axis("off")

        load_ax = controls_ax.inset_axes([0.02, 0.15, 0.12, 0.70])
        save_ax = controls_ax.inset_axes([0.16, 0.15, 0.12, 0.70])
        reset_ax = controls_ax.inset_axes([0.30, 0.15, 0.12, 0.70])
        fix_ax   = controls_ax.inset_axes([0.44, 0.15, 0.16, 0.70])
        zoom_ax  = controls_ax.inset_axes([0.62, 0.15, 0.22, 0.70])  # kept
        home_ax  = controls_ax.inset_axes([0.86, 0.15, 0.12, 0.70])  # moved on same row for clarity

        load_btn = Button(load_ax, "Load")
        save_btn = Button(save_ax, "Save")
        reset_btn = Button(reset_ax, "Reset")
        self.fix_btn  = Button(fix_ax,   "Fix Y-Limits")
        self.zoom_btn = Button(zoom_ax,  "Select Zoom (OFF)")
        home_btn      = Button(home_ax,  "Home View")

        load_btn.on_clicked(self._on_load_clicked)
        save_btn.on_clicked(self._on_save_clicked)
        reset_btn.on_clicked(self._on_reset)
        self.fix_btn.on_clicked(self._on_fix_toggle)
        self.zoom_btn.on_clicked(self._on_zoom_toggle)
        home_btn.on_clicked(self._on_home_view)

        # SpanSelector for X-zoom (disabled until toggled)
        self._span = SpanSelector(
            self.ax_top, self._on_span_select,
            direction='horizontal', useblit=True,
            props=dict(alpha=0.2), interactive=False
        )
        self._span.set_active(False)

        plt.show()

    # --- numerical helpers ---
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

    # --- autoscale helpers (respect zoom window) ---
    def _autoscale_error_axis(self, E, xlim=None):
        """Auto-adjust bottom y-limits to [0, 1.5*max(E in view)] if autoscaling is enabled."""
        if not self._auto_ylims:
            return
        if xlim is None:
            xlim = self.ax_top.get_xlim()
        xmin, xmax = xlim
        mask = (self.x >= xmin) & (self.x <= xmax)
        Ein = E[mask]
        ymax = np.nanmax(Ein) if np.any(np.isfinite(Ein)) else 0.0
        if not np.isfinite(ymax) or ymax <= 0:
            ymax = 1.0
        self.ax_bot.set_ylim(0.0, 1.5 * ymax)

    def _autoscale_top_axis(self, R, T, xlim=None, pad=0.04):
        """Autoscale top y-limits to visible R & tan with a small padding."""
        if xlim is None:
            xlim = self.ax_top.get_xlim()
        xmin, xmax = xlim
        mask = (self.x >= xmin) & (self.x <= xmax)
        Rin = R[mask]
        Tin = T[mask]
        data = np.concatenate([Rin[np.isfinite(Rin)], Tin[np.isfinite(Tin)]])
        if data.size == 0:
            return
        ymin, ymax = np.min(data), np.max(data)
        if ymin == ymax:
            # Avoid zero-height view
            delta = 1.0 if ymax == 0 else abs(ymax) * 0.5
            ymin, ymax = ymin - delta, ymax + delta
        span = ymax - ymin
        self.ax_top.set_ylim(ymin - pad*span, ymax + pad*span)

    # --- callbacks ---
    def _on_any_change(self, _):
        if self._suspend_update:
            return
        # gather coefficients
        a = [s.val for s in self.sliders_top]
        b = [s.val for s in self.sliders_bot]
        # gather centers (length degree-1)
        ac = [s.val for s in self.sliders_top_cent]
        bc = [s.val for s in self.sliders_bot_cent]

        P = self.eval_with_centers(self.x, a, ac) if a else np.zeros_like(self.x)
        Q = self.eval_with_centers(self.x, b, bc) if b else np.ones_like(self.x)
        R = self._safe_ratio(P, Q)
        T = self._safe_tan(self.x)
        E = self._safe_abs_error(R, T)

        self.line_ratio.set_ydata(R)
        self.line_tan.set_ydata(T)
        self.line_err.set_ydata(E)
        print(np.sum(E))

        # Respect current zoom window
        current_xlim = self.ax_top.get_xlim()
        self._autoscale_top_axis(R, T, xlim=current_xlim)
        self._autoscale_error_axis(E, xlim=current_xlim)

        self.fig.canvas.draw_idle()

    def _on_reset(self, _event):
        for s in self.sliders_top + self.sliders_bot + self.sliders_top_cent + self.sliders_bot_cent:
            s.reset()

    # --- Zoom controls ---
    def _on_zoom_toggle(self, _event):
        self._zoom_active = not self._zoom_active
        self._span.set_active(self._zoom_active)
        self.zoom_btn.label.set_text("Select Zoom (ON)" if self._zoom_active else "Select Zoom (OFF)")

    def _on_span_select(self, xmin, xmax):
        if xmin == xmax:
            return
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        self.ax_top.set_xlim(xmin, xmax)
        R = self.line_ratio.get_ydata()
        T = self.line_tan.get_ydata()
        E = self.line_err.get_ydata()
        self._autoscale_top_axis(R, T, xlim=(xmin, xmax))
        self._autoscale_error_axis(E, xlim=(xmin, xmax))
        self.fig.canvas.draw_idle()

    def _on_home_view(self, _event):
        if self._orig_xlim is not None:
            self.ax_top.set_xlim(*self._orig_xlim)
        if self._orig_ylim_top is not None:
            self.ax_top.set_ylim(*self._orig_ylim_top)
        if self._orig_ylim_bot is not None and not self._auto_ylims:
            self.ax_bot.set_ylim(*self._orig_ylim_bot)
        R = self.line_ratio.get_ydata()
        T = self.line_tan.get_ydata()
        E = self.line_err.get_ydata()
        self._autoscale_top_axis(R, T, xlim=self._orig_xlim)
        self._autoscale_error_axis(E, xlim=self._orig_xlim)
        self.fig.canvas.draw_idle()

    # --- Fix/Auto toggle handler ---
    def _on_fix_toggle(self, _event):
        self._auto_ylims = not self._auto_ylims
        self.fix_btn.label.set_text("Fix Y-Limits" if self._auto_ylims else "Auto Y-Limits")
        if self._auto_ylims:
            E = self.line_err.get_ydata()
            self._autoscale_error_axis(E)
        self.fig.canvas.draw_idle()

    # --- Save/Load helpers (centers included) ---
    def _collect_state(self):
        return {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "xlim": list(self.xlim),
            "top_degree": self.top_degree,
            "bot_degree": self.bot_degree,
            "sliders_top": [s.val for s in self.sliders_top],
            "sliders_bot": [s.val for s in self.sliders_bot],
            "sliders_top_cent": [s.val for s in self.sliders_top_cent],
            "sliders_bot_cent": [s.val for s in self.sliders_bot_cent],
        }

    def _apply_state(self, state):
        self._suspend_update = True
        try:
            # coeffs
            for i, v in enumerate(state.get("sliders_top", [])):
                if i < len(self.sliders_top):
                    vmin, vmax = (-2.0 / (i + 1.0), 2.0 / (i + 1.0))
                    self.sliders_top[i].set_val(float(np.clip(v, vmin, vmax)))
            for j, v in enumerate(state.get("sliders_bot", [])):
                if j < len(self.sliders_bot):
                    vmin, vmax = (-2.0 / (j + 1.0), 2.0 / (j + 1.0))
                    self.sliders_bot[j].set_val(float(np.clip(v, vmin, vmax)))
            # centers
            for i, v in enumerate(state.get("sliders_top_cent", [])):
                if i < len(self.sliders_top_cent):
                    vmin, vmax = (-2.0 / (i + 1.0), 2.0 / (i + 1.0))
                    self.sliders_top_cent[i].set_val(float(np.clip(v, vmin, vmax)))
            for j, v in enumerate(state.get("sliders_bot_cent", [])):
                if j < len(self.sliders_bot_cent):
                    vmin, vmax = (-2.0 / (j + 1.0), 2.0 / (j + 1.0))
                    self.sliders_bot_cent[j].set_val(float(np.clip(v, vmin, vmax)))
        finally:
            self._suspend_update = False
        self._on_any_change(None)

    def _on_save_clicked(self, _event):
        state = self._collect_state()
        try:
            with open(self.state_path, "w", encoding="utf-8") as f:
                json.dump(state, f, indent=2)
            print(f"[saved] {self.state_path} at {state['timestamp']}")
        except Exception as e:
            print(f"[save error] {e}")

    def _on_load_clicked(self, _event):
        try:
            with open(self.state_path, "r", encoding="utf-8") as f:
                state = json.load(f)
            self._apply_state(state)
            print(f"[loaded] {self.state_path} from {state.get('timestamp','unknown time')}")
        except FileNotFoundError:
            print(f"[load] No saved state found at {self.state_path}")
        except Exception as e:
            print(f"[load error] {e}")

    # --- Taylor coefficients for tan(x) ---
    def _taylor_coeffs_tan(self, max_degree):
        coeffs = [0.0] * max(0, max_degree)
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
        if max_degree >= 1:
            coeffs[0] = 0.0
        return coeffs


if __name__ == "__main__":
    interactive_plot = InteractivePlot(9, 5, (0, np.pi/2-0.01))
    interactive_plot.plot()
