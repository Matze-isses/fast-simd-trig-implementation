import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import brentq



def taylor(degree):
    """
    Returns a function that computes the Taylor polynomial of tan(x)
    of the given degree around point a (default 0).
    """

    x = sp.symbols('x')
    f = sp.tan(x)
    
    # Taylor expansion around a
    taylor_series = sp.series(f, x, 0, degree + 1)
    taylor_poly = taylor_series.removeO()
    
    # Convert to a numerical function
    f_taylor = sp.lambdify(x, taylor_poly, 'numpy')
    
    return f_taylor


comparison_function = "tan"

lower = np.pi/8
upper = np.pi/4
# upper = 1.5

x_vals = np.linspace(lower, upper, 1000000)
y_true = np.tan(x_vals)

func = taylor(15)
y_own = func(x_vals)

y = np.linspace(lower, upper, 100000)
dx = []
for y0 in y:
    # Solve true(x) = y0
    f1 = lambda x: func(x) - y0
    x_true = brentq(f1, 0, 10000)

    # Solve approx(x) = y0
    f2 = lambda x: np.tan(x) - y0
    x_approx = brentq(f2, 0, upper)

    dx.append(x_true - x_approx)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 1, 1]})

ax1.plot(x_vals, y_true, label="True Solution")

ax1.plot(x_vals, np.array(y_own), label="Taylor Poly")
ax1.legend()

error = y_true - np.array(y_own)

one_over_x = 1/(np.pi/2-x_vals)
adjusted_correcture = error[-1]/one_over_x[-1]*one_over_x

ax2.plot(x_vals, error, label="Error")
ax2.plot(x_vals, adjusted_correcture, label=comparison_function)
ax2.set_title("Error Plot")
ax2.legend()

ax3.plot(x_vals, error - adjusted_correcture, label="Left Error - Correcture")
ax2.set_title("Correcture Error Plot")

plt.savefig(f"./plots/error_fkt_{comparison_function}.png")

plt.show()


dx = np.array(dx)


# -------------------------
# Bézier helpers
# -------------------------
def bezier_curve(P0, P1, P2, P3, n=400):
    t = np.linspace(0, 1, n)
    B0 = (1 - t) ** 3
    B1 = 3 * (1 - t) ** 2 * t
    B2 = 3 * (1 - t) * t ** 2
    B3 = t ** 3
    pts = (B0[:, None] * P0 +
           B1[:, None] * P1 +
           B2[:, None] * P2 +
           B3[:, None] * P3)
    return pts[:, 0], pts[:, 1]

# -------------------------
# Setup plot
# -------------------------
fig, ax = plt.subplots()
ax.plot(y, dx, label="x diff for same y values")

ax.set_xlabel("y")
ax.set_ylabel("Δx")
ax.grid(True)

# Endpoints fixed to start/end of curve
P0 = np.array([y[0],  dx[0]])
P3 = np.array([y[-1], dx[-1]])

# Initial guesses for control points (on chord)
P1 = P0 + (P3 - P0) / 3
P2 = P0 + 2 * (P3 - P0) / 3

# Initial Bézier
bx, by = bezier_curve(P0, P1, P2, P3)
bez_line, = ax.plot(bx, by, lw=2, label="cubic Bézier")

# Control polygon (optional visual)
ctrl_line, = ax.plot(
    [P0[0], P1[0], P2[0], P3[0]],
    [P0[1], P1[1], P2[1], P3[1]],
    "--", lw=1, color="gray"
)

# Draggable control points (green)
P1_artist, = ax.plot([P1[0]], [P1[1]], "o", ms=10, color="green")
P2_artist, = ax.plot([P2[0]], [P2[1]], "o", ms=10, color="green")

ax.legend()

# -------------------------
# Drag logic
# -------------------------
# We’ll store references in a dict so we can modify inside callbacks
state = {
    "dragging": None,  # None, "P1", or "P2"
    "P1": P1,
    "P2": P2,
}

def pick_control_point(event):
    """Return 'P1', 'P2', or None depending on click proximity."""
    if event.xdata is None or event.ydata is None:
        return None

    # Simple distance in data coords; adjust threshold as you like
    x, y_ = event.xdata, event.ydata
    d1 = np.hypot(x - state["P1"][0], y_ - state["P1"][1])
    d2 = np.hypot(x - state["P2"][0], y_ - state["P2"][1])
    thresh = 0.3  # data units

    if d1 < d2 and d1 < thresh:
        return "P1"
    elif d2 <= d1 and d2 < thresh:
        return "P2"
    else:
        return None

def on_press(event):
    if event.inaxes != ax:
        return
    cp = pick_control_point(event)
    state["dragging"] = cp

def on_release(event):
    if state["dragging"] is not None:
        # Print the coordinates of both control points
        P1 = state["P1"]
        P2 = state["P2"]
        print(f"P1 = ({P1[0]:.6f}, {P1[1]:.6f}),   P2 = ({P2[0]:.6f}, {P2[1]:.6f})")

    state["dragging"] = None

def on_motion(event):
    if state["dragging"] is None:
        return
    if event.inaxes != ax or event.xdata is None or event.ydata is None:
        return

    # Update the corresponding control point
    name = state["dragging"]
    state[name][0] = event.xdata
    state[name][1] = event.ydata

    # Recompute Bézier with new control points
    P1 = state["P1"]
    P2 = state["P2"]
    bx, by = bezier_curve(P0, P1, P2, P3)
    bez_line.set_data(bx, by)

    # Update control polygon & point artists
    ctrl_line.set_data(
        [P0[0], P1[0], P2[0], P3[0]],
        [P0[1], P1[1], P2[1], P3[1]]
    )
    P1_artist.set_data([P1[0]], [P1[1]])
    P2_artist.set_data([P2[0]], [P2[1]])

    fig.canvas.draw_idle()

# P1 = (0.649244, -0.000000),   P2 = (0.712166, -0.000000)
# Connect callbacks
fig.canvas.mpl_connect("button_press_event", on_press)
fig.canvas.mpl_connect("button_release_event", on_release)
fig.canvas.mpl_connect("motion_notify_event", on_motion)

plt.show()
