import tensorflow as tf
import matplotlib.pyplot as plt
import numpy as np
import json, os, tempfile

# ======================================================
# Configurable parameters
# ======================================================
TOP_NUM_PARAMS = 20
BOT_NUM_PARAMS = 20
TOP_NUM_CENTERS = TOP_NUM_PARAMS - 1
BOT_NUM_CENTERS = BOT_NUM_PARAMS - 1

X_MIN = 1.57
X_MAX = 1.5707933306386703
NUM_POINTS = 100_000

LEARNING_RATE = 1e-3
TARGET_MAIN_SUM = 1e-5
MAX_STEPS = 10_000_000
PLOT_EVERY = 25_000  # plot refresh interval

# Parameter projection bounds
PARAM_MIN = -8.0
PARAM_MAX = 8.0

# JSON saving
BEST_JSON_PATH = "best_params.json"
SAVE_AFTER_STEP = 100_000

# ======================================================
# Setup
# ======================================================

x = np.array([], dtype=np.float64)
weights = [1.0]
points = [X_MIN, X_MAX]

for i, w in enumerate(weights):
    x1 = np.linspace(points[i], points[i+1], (int)(w * NUM_POINTS), dtype=np.float64)
    x = np.concatenate([x1[:-1], x])

DTYPE = tf.float64
x = tf.linspace(tf.cast(X_MIN, DTYPE), tf.cast(X_MAX, DTYPE), NUM_POINTS)

def init_centers(num_centers: int):
    if num_centers <= 0:
        return tf.zeros([0], dtype=DTYPE)
    return tf.linspace(tf.cast(X_MIN, DTYPE), tf.cast(X_MAX, DTYPE), num_centers + 2)[1:-1]

# ---- Trainables ----
top_param   = tf.Variable(tf.zeros([TOP_NUM_PARAMS], dtype=DTYPE))
bot_param   = tf.Variable(tf.zeros([BOT_NUM_PARAMS], dtype=DTYPE))
top_centers = tf.Variable(init_centers(TOP_NUM_CENTERS))
bot_centers = tf.Variable(init_centers(BOT_NUM_CENTERS))

# ======================================================
# Fixed-size Horner evaluators (XLA-friendly)
# ======================================================
@tf.function(jit_compile=True)
def eval_shifted_horner_top(a, c, x_):
    y = tf.fill(tf.shape(x_), a[0])
    for i in range(TOP_NUM_CENTERS):
        y = y * (x_ - c[i]) + a[i + 1]
    return y

@tf.function(jit_compile=True)
def eval_shifted_horner_bot(a, c, x_):
    y = tf.fill(tf.shape(x_), a[0])
    for i in range(BOT_NUM_CENTERS):
        y = y * (x_ - c[i]) + a[i + 1]
    return y


# ======================================================
# Objective
# ======================================================
@tf.function(jit_compile=True)
def objective(top_param, bot_param, top_centers, bot_centers):
    y1 = eval_shifted_horner_top(top_param, top_centers, x)
    y2 = eval_shifted_horner_bot(bot_param, bot_centers, x)
    eps = tf.cast(1e-12, DTYPE)
    result_y = y1 / (y2 + eps)
    # Keep L1 objective for training stability
    main_loss = tf.reduce_sum(tf.abs(tf.tan(x) - result_y))
    reg = 1e-8 * (
        tf.reduce_sum(tf.square(top_param)) + tf.reduce_sum(tf.square(bot_param)) +
        tf.reduce_sum(tf.square(top_centers)) + tf.reduce_sum(tf.square(bot_centers))
    )
    return main_loss + reg, main_loss

optimizer = tf.keras.optimizers.Adam(learning_rate=LEARNING_RATE)

@tf.function(jit_compile=True)
def _project_bounds():
    top_param.assign(tf.clip_by_value(top_param, PARAM_MIN, PARAM_MAX))
    bot_param.assign(tf.clip_by_value(bot_param, PARAM_MIN, PARAM_MAX))
    top_centers.assign(tf.clip_by_value(top_centers, PARAM_MIN, PARAM_MAX))
    bot_centers.assign(tf.clip_by_value(bot_centers, PARAM_MIN, PARAM_MAX))

@tf.function(jit_compile=True)
def train_step():
    with tf.GradientTape() as tape:
        total_loss, main_loss = objective(top_param, bot_param, top_centers, bot_centers)

    grads = tape.gradient(total_loss, [top_param, bot_param, top_centers, bot_centers])
    optimizer.apply_gradients(zip(grads, [top_param, bot_param, top_centers, bot_centers]))
    _project_bounds()
    return total_loss, main_loss

# ======================================================
# Plot setup (interactive)
# ======================================================
plt.ion()
fig, axes = plt.subplots(3, 1, figsize=(19.2, 10.8), constrained_layout=True)  # Full HD @ dpi=100
fig.set_dpi(100)
ax_f, ax_err, ax_loss = axes

x_np = x.numpy()
tan_np = np.tan(x_np)

# function + error plots
approx_line, = ax_f.plot(x_np, np.zeros_like(tan_np), linewidth=1, label="Approx")
tan_line,    = ax_f.plot(x_np, tan_np, linewidth=1, alpha=0.7, label="tan(x)")
ax_f.set_title(f"Approx vs tan(x) — top deg {TOP_NUM_CENTERS}, bottom deg {BOT_NUM_CENTERS}")
ax_f.set_xlim([X_MIN, X_MAX])
ax_f.legend(loc="upper left")
ax_f.set_yscale('symlog', linthresh=1e-3, base=10)

# --- LEFTOVER (signed) error ---
err_line, = ax_err.plot(x_np, np.zeros_like(tan_np), linewidth=1, label="leftover = tan(x) - approx(x)")
ax_err.set_title("Leftover Error (signed)")
ax_err.set_xlim([X_MIN, X_MAX])
ax_err.set_yscale('symlog', linthresh=1e-3, base=10)
ax_err.legend(loc="upper left")

# total loss evolution plot
ax_loss.set_title("Total Loss Evolution")
ax_loss.set_xlabel("Step")
ax_loss.set_ylabel("Total Loss")
ax_loss.set_yscale('symlog', linthresh=1e-3, base=10)

loss_points, = ax_loss.plot([], [], "bo-", markersize=3, linewidth=1, label="total_loss")
ax_loss.legend(loc="upper right")

@tf.function(jit_compile=True)
def _current_approx():
    y1 = eval_shifted_horner_top(top_param, top_centers, x)
    y2 = eval_shifted_horner_bot(bot_param, bot_centers, x)
    eps = tf.cast(1e-12, DTYPE)
    return y1 / (y2 + eps)

def _save_current_figure(step:int):
    # Ensure directory exists
    out_dir = os.path.join(".", "plots", "tf_optimized")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"iter{step}.png")
    # Ensure Full HD output size
    prev_size = fig.get_size_inches()
    prev_dpi = fig.dpi
    try:
        fig.set_size_inches(19.2, 10.8)  # inches for 1920x1080 @ 100 dpi
        fig.set_dpi(100)
        fig.savefig(out_path, dpi=100, bbox_inches="tight")
    finally:
        # Restore previous settings for interactive session consistency
        fig.set_size_inches(prev_size)
        fig.set_dpi(prev_dpi)

def update_plot(step, main_loss_val, loss_hist_steps, loss_hist_vals):
    approx = _current_approx().numpy()
    leftover = tan_np - approx  # signed error
    approx_line.set_ydata(approx)
    err_line.set_ydata(leftover)

    # update loss evolution
    loss_points.set_data(loss_hist_steps, loss_hist_vals)
    ax_loss.relim()
    ax_loss.autoscale_view()

    ax_f.relim(); ax_f.autoscale_view()
    ax_err.relim(); ax_err.autoscale_view()
    fig.suptitle(
        f"Step {step} | sum L1 error = {main_loss_val:.6e} | "
        f"top deg {TOP_NUM_CENTERS}, bottom deg {BOT_NUM_CENTERS}",
        fontsize=12
    )
    fig.canvas.draw_idle()
    plt.pause(0.5)

    # --- Save figure each time we draw it ---
    _save_current_figure(step)

# Initial plot
update_plot(step=0, main_loss_val=np.inf, loss_hist_steps=[], loss_hist_vals=[])

# ======================================================
# Training loop with best-tracking, JSON save & loss tracking
# ======================================================
best_total = float("inf")
best_main  = float("inf")

loss_hist_steps = []
loss_hist_vals  = []

def _save_best_json(step:int, total_loss:float, main_loss:float):
    payload = {
        "step": step,
        "total_loss": float(total_loss),
        "main_loss": float(main_loss),
        "top_param": top_param.numpy().astype(float).tolist(),
        "bot_param": bot_param.numpy().astype(float).tolist(),
        "top_centers": top_centers.numpy().astype(float).tolist(),
        "bot_centers": bot_centers.numpy().astype(float).tolist(),
        "bounds": {"min": PARAM_MIN, "max": PARAM_MAX},
        "domain": {"X_MIN": X_MIN, "X_MAX": X_MAX, "NUM_POINTS": NUM_POINTS}
    }
    dir_name = os.path.dirname(BEST_JSON_PATH) or "."
    os.makedirs(dir_name, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=dir_name) as tmp:
        json.dump(payload, tmp, indent=2)
        tmp_path = tmp.name
    os.replace(tmp_path, BEST_JSON_PATH)

for step in range(1, MAX_STEPS + 1):
    total_loss_t, main_loss_t = train_step()
    total_loss = float(total_loss_t.numpy())
    main_loss  = float(main_loss_t.numpy())

    # track best loss and save if new best found after threshold
    if total_loss < best_total:
        best_total = total_loss
        best_main  = min(best_main, main_loss)
        if step >= SAVE_AFTER_STEP:
            _save_best_json(step, total_loss, main_loss)

    # PRINT STATEMENT (exactly as required)
    if step % 100 == 0:
        tf.print("step", step, f"total_loss: {total_loss_t:15.5f}", f"main_sum_loss {main_loss_t:15.5f}")

        if total_loss < 100_000_000:
            loss_hist_steps.append(step)
            loss_hist_vals.append(total_loss)

    # Plot update only every PLOT_EVERY steps
    if step % PLOT_EVERY == 0:
        update_plot(step, main_loss, loss_hist_steps, loss_hist_vals)

    if main_loss_t < tf.cast(TARGET_MAIN_SUM, DTYPE):
        tf.print("EARLY STOP at step", step, "main_sum_loss", main_loss_t)
        update_plot(step, main_loss, loss_hist_steps, loss_hist_vals)
        if step >= SAVE_AFTER_STEP and total_loss <= best_total:
            _save_best_json(step, total_loss, main_loss)
        break

plt.ioff()
plt.show()

# ---- Final report ----
tf.print("\nBest main sum loss (observed):", best_main)
tf.print("Best total loss (observed):", best_total)
tf.print("top_param:", top_param)
tf.print("bot_param:", bot_param)
tf.print("top_centers:", top_centers)
tf.print("bot_centers:", bot_centers)
print(f"\nBest parameters (when saving enabled) are written to: {BEST_JSON_PATH}")
