import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# 1) Your data
# -----------------------------
x = np.array([
    0.39269908169872414, 0.5157325170655326, 0.638765952432341,
    0.7617993877991494, 0.8848328231659579, 1.0078662585327662,
    1.1308996938995746, 1.2539331292663831, 1.3769665646331917, 1.5
])
y = np.array([
    7.905309740152688e-11, 8.54329873511972e-09, 3.4688503047775043e-07,
    7.561642495379495e-06, 0.00010796649218769439, 0.001145809455095792,
    0.009915949668291457, 0.07619521670525842, 0.5859987186389297,
    6.5978377206998085
])

# -----------------------------
# 2) Build weighted fit for 1/P(x)
# -----------------------------

# Reciprocal data
z = 1.0 / y

# Choose polynomial degree for P(x)
deg = 4     # you can play with 2,3,5,...

# Weights: make large y more important
# (square to reflect least-squares on original y)
w = y**2

# Fit polynomial P(x) ≈ z with weights
# coeffs = [a_n, ..., a1, a0] for P(x) = a_n x^n + ... + a0
coeffs = np.polyfit(x, z, deg, w=w)

print("Fitted coefficients for P(x) (highest degree first):")
print(coeffs)

# Define P(x) and f(x) = 1/P(x)
def P(x_val):
    return np.polyval(coeffs, x_val)

def f(x_val):
    return 1.0 / P(x_val)

# -----------------------------
# 3) Check fit quality at data points
# -----------------------------
y_hat = f(x)
residuals = y_hat - y
rmse = np.sqrt(np.mean(residuals**2))
print("RMSE on data:", rmse)

# Optionally, print each point’s approximation
for xi, yi, yh in zip(x, y, y_hat):
    print(f"x = {xi:.6f}, y = {yi:.6e}, y_hat = {yh:.6e}")

# -----------------------------
# 4) Plot data vs fitted function
# -----------------------------
xs = np.linspace(x.min(), x.max(), 400)
ys = f(xs)

plt.figure(figsize=(8,5))
plt.plot(xs, ys, label="Weighted 1 / P(x) fit")
plt.scatter(x, y, label="Data points")
plt.yscale("log")  # optional but very helpful with your huge range
plt.xlabel("x")
plt.ylabel("y")
plt.title("Data vs weighted reciprocal-polynomial fit")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.shoyw()
