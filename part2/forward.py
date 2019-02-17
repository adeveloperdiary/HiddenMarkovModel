import pandas as pd
import numpy as np

data = pd.read_csv('data_python.csv')

V = data['Visible'].values

# Transition Probabilities
a = np.array(((0.54, 0.46), (0.49, 0.51)))

# Emission Probabilities
b = np.array(((0.16, 0.26, 0.58), (0.25, 0.28, 0.47)))

# Equal Probabilities for the initial distribution
initial_distribution = np.array((0.5, 0.5))


def forward(V, a, b, initial_distribution):
    alpha = np.zeros((V.shape[0], a.shape[0]))
    alpha[0, :] = initial_distribution * b[:, V[0]]

    for t in range(1, V.shape[0]):
        for j in range(a.shape[0]):
            # Matrix Computation Steps
            #                  ((1x2) . (1x2)) * (1)
            #                  (1) * (1)
            alpha[t, j] = alpha[t - 1].dot(a[j, :]) * b[j, V[t]]

    return alpha


alpha = forward(V, a, b, initial_distribution)
print(alpha)
