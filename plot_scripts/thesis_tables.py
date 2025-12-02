import numpy as np

fak = 1

for i in range(1, 16):
    fak *= (i+1)
    x = (np.pi/16)**(i+1) / fak

    print(f"n: {i}; error: {x}")
