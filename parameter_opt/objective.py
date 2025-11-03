
import numpy as np


def objective_function(
        top_param: np.ndarray,
        bot_param: np.ndarray,
        top_centers: np.ndarray,
        bot_centers: np.ndarray
    ):

    x = np.linspace(0, 1.57, 100000)

    y1 = np.zeros_like(x, dtype=float)
    for a, c in zip(reversed(top_param), reversed(top_centers)):
        y1 = y1 * (x-c) + a

    y2 = np.zeros_like(x, dtype=float)
    for a, c in zip(reversed(bot_param), reversed(bot_centers)):
        y2 = y2 * (x-c) + a

    result_y = y1/y2

    return np.sum(np.abs(np.tan(x) - result_y))
