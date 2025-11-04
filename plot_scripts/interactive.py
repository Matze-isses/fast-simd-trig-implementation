import json
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, SpanSelector  # <-- added SpanSelector


class InteractivePlot:
    def __init__(self, top_degree, bot_degree, xlim=(-2*np.pi, 2*np.pi)) -> None:
        self.top_degree = int(top_degree)
        self.bot_degree = int(bot_degree)

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
        self.sliders_top = []
        self.sliders_bot = []

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

    # --- polynomial evaluator (Horner) ---
    def eval_sliders(self, x, *coeffs):
        y = np.zeros_like(x, dtype=float)
        for a in reversed(coeffs):
            y = y * x + a
        return y

    def plot(self):
        n_sliders = self.top_degree + self.bot_degree
        self.fig = plt.figure(figsize=(10, 7 + 0.32 * max(1, n_sliders)))
        gs = self.fig.add_gridspec(
            nrows=4, ncols=1,
            height_ratios=[4, 3, 0.22 * max(1, n_sliders), 0.2],
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

        # Initial coefficients
        init_top = self._taylor_coeffs_tan(self.top_degree) if self.top_degree else []
        init_bot = [1.0] + [0.0] * (self.bot_degree - 1) if self.bot_degree else []

        # Initial curves
        P = self.eval_sliders(self.x, *init_top) if self.top_degree else np.zeros_like(self.x)
        Q = self.eval_sliders(self.x, *init_bot) if self.bot_degree else np.ones_like(self.x)
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

        # Sliders area
        sliders_ax = self.fig.add_subplot(gs[2, 0])
        sliders_ax.axis("off")

        def degree_range(i: int):
            r = 2.0 / (i + 1.0)
            return -r, r

        def add_slider(y0, label, vinit, vmin, vmax):
            rect = [0.06, y0, 0.88, 0.065]
            ax = sliders_ax.inset_axes(rect)
            vinit = float(np.clip(vinit, vmin, vmax))
            s = Slider(ax=ax, label=label, valmin=vmin, valmax=vmax, valinit=vinit, valstep=1e-5)
            return s

        pad = 0.082
        ycursor = 0.92

        # Top sliders
        for i in range(self.top_degree):
            vmin, vmax = degree_range(i)
            vinit = init_top[i] if i < len(init_top) else 0.0
            s = add_slider(ycursor, f"a{i}", vinit=vinit, vmin=vmin, vmax=vmax)
            s.on_changed(self._on_any_change)
            self.sliders_top.append(s)
            ycursor -= pad

        # Bottom sliders
        for j in range(self.bot_degree):
            vmin, vmax = degree_range(j)
            vinit = init_bot[j] if j < len(init_bot) else 0.0
            s = add_slider(ycursor, f"b{j}", vinit=vinit, vmin=vmin, vmax=vmax)
            s.on_changed(self._on_any_change)
            self.sliders_bot.append(s)
            ycursor -= pad

        # --- Controls row ---
        controls_ax = self.fig.add_subplot(gs[3, 0])
        controls_ax.axis("off")

        load_ax = controls_ax.inset_axes([0.02, 0.15, 0.15, 0.70])
        save_ax = controls_ax.inset_axes([0.19, 0.15, 0.15, 0.70])
        reset_ax = controls_ax.inset_axes([0.36, 0.15, 0.15, 0.70])
        fix_ax   = controls_ax.inset_axes([0.53, 0.15, 0.18, 0.70])
        zoom_ax  = controls_ax.inset_axes([0.73, 0.15, 0.24, 0.70])  # new
        home_ax  = controls_ax.inset_axes([0.73, 0.02, 0.24, 0.10])  # tiny bar below zoom

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
        # add small proportional padding
        span = ymax - ymin
        self.ax_top.set_ylim(ymin - pad*span, ymax + pad*span)

    # --- callbacks ---
    def _on_any_change(self, _):
        if self._suspend_update:
            return
        a = [s.val for s in self.sliders_top]
        b = [s.val for s in self.sliders_bot]
        P = self.eval_sliders(self.x, *a)
        Q = self.eval_sliders(self.x, *b)
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
        for s in self.sliders_top + self.sliders_bot:
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
        # Set new x-limits (sharex syncs bottom)
        self.ax_top.set_xlim(xmin, xmax)
        # Re-run autoscale for y in this window
        R = self.line_ratio.get_ydata()
        T = self.line_tan.get_ydata()
        E = self.line_err.get_ydata()
        self._autoscale_top_axis(R, T, xlim=(xmin, xmax))
        self._autoscale_error_axis(E, xlim=(xmin, xmax))
        self.fig.canvas.draw_idle()

    def _on_home_view(self, _event):
        # Restore original limits
        if self._orig_xlim is not None:
            self.ax_top.set_xlim(*self._orig_xlim)
        if self._orig_ylim_top is not None:
            self.ax_top.set_ylim(*self._orig_ylim_top)
        if self._orig_ylim_bot is not None and not self._auto_ylims:
            # Only restore fixed limits if autoscaling is disabled
            self.ax_bot.set_ylim(*self._orig_ylim_bot)
        # If autoscaling is enabled, recompute for full range
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
        # If re-enabling autoscale, refresh immediately under current view
        if self._auto_ylims:
            E = self.line_err.get_ydata()
            self._autoscale_error_axis(E)
        self.fig.canvas.draw_idle()

    # --- Save/Load helpers (unchanged) ---
    def _collect_state(self):
        return {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "xlim": list(self.xlim),
            "top_degree": self.top_degree,
            "bot_degree": self.bot_degree,
            "sliders_top": [s.val for s in self.sliders_top],
            "sliders_bot": [s.val for s in self.sliders_bot],
        }

    def _apply_state(self, state):
        self._suspend_update = True
        try:
            for i, v in enumerate(state.get("sliders_top", [])):
                if i < len(self.sliders_top):
                    vmin, vmax = (-2.0 / (i + 1.0), 2.0 / (i + 1.0))
                    self.sliders_top[i].set_val(float(np.clip(v, vmin, vmax)))
            for j, v in enumerate(state.get("sliders_bot", [])):
                if j < len(self.sliders_bot):
                    vmin, vmax = (-2.0 / (j + 1.0), 2.0 / (j + 1.0))
                    self.sliders_bot[j].set_val(float(np.clip(v, vmin, vmax)))
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
