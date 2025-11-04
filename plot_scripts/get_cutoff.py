import time
import numpy as np

cutoff = 0.1
error = 1

while error > 1/100:
    x = np.pi/2 - cutoff
    ls = np.linspace(x, np.pi/2, 10000)[:-1]

    error = np.sum(np.tan(ls) - (1/ls)) / 9999
    time.sleep(0.5)
    print(f"Error: {error}; Cutoff: {x}")
    cutoff = cutoff * 0.9
