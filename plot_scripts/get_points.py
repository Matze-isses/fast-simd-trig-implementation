import numpy as np
import matplotlib.pyplot as plt


def get_data(path="./tan_ulp_error_behavior.tsv"):
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    x = data[:, 0]
    err = data[:, 1]
    return np.array(x), np.array(err)


if __name__ == "__main__":
    x_values, error_values = get_data()

    positive_error = []
    negative_error = []

    for x, y in zip(x_values, error_values):
        if y == 3:
            positive_error.append(x)

        elif y == -3:
            negative_error.append(x)

    for e1, e2 in zip(positive_error, negative_error):
        print(e1.hex())
        print(e2.hex())

    print("\n")

    print(f"Min/Max Positive: {np.min(positive_error)}/{np.max(positive_error)}")
    print(f"Min/Max Negative: {np.min(negative_error)}/{np.max(negative_error)}")
