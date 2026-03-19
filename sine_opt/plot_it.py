import numpy as np
import matplotlib.pyplot as plt

def load_tsv(filename):
    data = np.loadtxt(filename, delimiter="\t", skiprows=1)
    x = data[:, 0]
    approx = data[:, 1]
    ref = data[:, 2]
    ulp = data[:, 3]
    return x, approx, ref, ulp


def plot_results(x, approx, ref, ulp):
    plt.figure()

    plt.plot(x**2, ulp)

    plt.xlabel("x")
    plt.ylabel("ULP error")
    plt.title("Signed ULP Error")
    plt.grid()

    plt.show()

def main():
    filename = "sin_results.tsv"
    x, approx, ref, ulp = load_tsv(filename)
    plot_results(x, approx, ref, ulp)

if __name__ == "__main__":
    main()
