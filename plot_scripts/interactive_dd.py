import json
from datetime import datetime
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, SpanSelector  # <-- added SpanSelector


def load_json(path="results.json"):
    with open(path, "r+") as file:
        data = json.load(file)

    return data


def save_json(top_param, top_centers, bot_param, bot_centers, path="results.json"):
    data_dict = {
        "top_param": list(top_param),
        "bot_param": list(bot_param),
        "top_centers": list(top_centers),
        "bot_centers": list(bot_centers),
    }
    with open(path, "w+") as file:
        json.dump(data_dict, file)


class Interactive:

    def __init__(self, num_top_param=10, num_bot_param=10) -> None:
        self.top_param = np.zeros(num_top_param)
        self.top_centers = np.zeros(num_top_param-1)

        self.bot_param = np.ones(num_bot_param)
        self.bot_centers = np.zeros(num_bot_param-1)
        self.x = np.linspace(0, 1.5, 100000)

        if os.path.exists("results.json"):
            data = load_json()
            self.top_param = np.array(data["top_param"])
            self.bot_param = np.array(data["bot_param"])
            self.top_centers = np.array(data["top_centers"])
            self.bot_centers = np.array(data["bot_centers"])

        # Matplotlib objects
        self.fig = None
        self.ax_top_graph = None
        self.ax_bot_graph = None

        self.ax_top_sliders = None
        self.ax_bot_sliders = None

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

    def _mixed_poly_approx(self, x):
        y1 = np.zeros_like(x, dtype=float)
        for a, c in zip(reversed(self.top_param), reversed(self.top_centers)):
            y1 = y1 * (x-c) + a

        y2 = np.zeros_like(x, dtype=float)
        for a, c in zip(reversed(self.bot_param), reversed(self.bot_centers)):
            y2 = y2 * (x-c) + a

        return y1 / y2

    def _on_any_change(self, _):
        pass

    def plot(self):
        n_top_sliders = 2 * self.top_param.shape[0] - 1
        n_bot_sliders = 2 * self.bot_param.shape[0] - 1
        self.fig = plt.figure(figsize=(10, 10))

        gs = self.fig.add_gridspec(
            nrows=2, ncols=2,
            height_ratios=[1, 1],
            width_ratios=[1, 1],
            hspace=0.30
        )

        self.ax_top_graph = self.fig.add_subplot(gs[0, 0])
        self.ax_bot_graph = self.fig.add_subplot(gs[1, 0], sharex=self.ax_top_graph)

        self.ax_top_sliders = self.fig.add_subplot(gs[0, 1])
        self.ax_bot_sliders = self.fig.add_subplot(gs[1, 1])


        tan_vals = np.tan(self.x)
        own_approx = self._mixed_poly_approx(self.x)

        (self.line_ratio,) = self.ax_top_graph.plot(self.x, own_approx, lw=1, label="Tangens")
        (self.line_tan,) = self.ax_top_graph.plot(self.x, tan_vals, lw=2, alpha=0.6, linestyle="-", label="tan(x)")
        (self.line_err,) = self.ax_bot_graph.plot(self.x, np.abs(tan_vals - own_approx), lw=2, label="|R - tan|")


        def degree_range(i: int):
            r = 2.0 / (i + 1.0)
            return -r, r

        def add_slider(axes, y0, label, vinit, vmin, vmax):
            rect = [0.06, y0, 0.88, 0.065]
            ax = axes.inset_axes(rect)
            vinit = float(np.clip(vinit, vmin, vmax))
            s = Slider(ax=ax, label=label, valmin=vmin, valmax=vmax, valinit=vinit, valstep=1e-5)
            return s

        pad = 0.082
        ycursor = 0.92
        slider_inits_top = np.append(self.top_param, self.top_centers, axis=0)
        slider_inits_bot = np.append(self.bot_param, self.bot_centers, axis=0)

        # Top sliders
        for i in range(n_top_sliders):
            vmin, vmax = degree_range(i)
            vinit = slider_inits_top[i]
            s = add_slider(self.ax_top_sliders, ycursor, f"a{i}", vinit=vinit, vmin=vmin, vmax=vmax)

            s.on_changed(self._on_any_change)
            self.sliders_top.append(s)
            ycursor -= pad

        # Top sliders
        for i in range(n_bot_sliders):
            vmin, vmax = degree_range(i)
            vinit = slider_inits_bot[i]
            s = add_slider(self.ax_bot_sliders, ycursor, f"a{i}", vinit=vinit, vmin=vmin, vmax=vmax)

            s.on_changed(self._on_any_change)
            self.sliders_bot.append(s)
            ycursor -= pad

        plt.show()
        



if __name__ == "__main__":
    app = Interactive()
    app.plot()
