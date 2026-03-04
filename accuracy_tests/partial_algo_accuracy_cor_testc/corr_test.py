import json
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt

with open("corr_test.json") as f:
    data = np.array(json.load(f))

x = data[:,0]
ulp = data[:,1]
lsb = data[:,2]

# event: large error
large = np.abs(ulp) == 2


# contingency table
a = np.sum((lsb == 1) & large)
b = np.sum((lsb == 1) & (~large))
c = np.sum((lsb == 0) & large)
d = np.sum((lsb == 0) & (~large))

print("Total 2 Bit deviations: ", a + c)

table = [[a, b],
         [c, d]]

odds_ratio, p = fisher_exact(table)

print("contingency table:")
print(table)

print("odds ratio:", odds_ratio)
print("p-value:", p)

# also print probabilities
p1 = a / (a + b)
p0 = c / (c + d)

print("P(|ULP|=2 | LSB=1) =", p1)
print("P(|ULP|=2 | LSB=0) =", p0)

lsbone = lsb == 1

plt.scatter(x[lsbone], ulp[lsbone], label='ULP With LSB 1', s=0.2)
plt.scatter(x[~lsbone], ulp[~lsbone], label='ULP With LSB 0', s=0.2)
plt.legend()
plt.show()
