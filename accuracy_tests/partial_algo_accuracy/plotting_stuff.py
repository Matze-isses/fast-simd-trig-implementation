import json
import matplotlib.pyplot as plt

# Load data
with open("data.json", "r") as f:
    data = json.load(f)

# Split into x and y
x = [point[0] for point in data if abs(point[1]) < 100]
y = [point[1] for point in data if abs(point[1]) < 100]

# Plot
plt.figure()
plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Simple Plot")
plt.show()
