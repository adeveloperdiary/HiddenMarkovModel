import pandas as pd
import numpy as np

data = pd.read_csv('data_python.csv')

V = data['Visible'].values

a = np.array(((0.54, 0.46), (0.49, 0.51)))
b = np.array(((0.16, 0.26, 0.58), (0.25, 0.28, 0.47)))

initial_distribution = np.array((0.5, 0.5))


def forward(X, a, b, initial_distribution):
    alpha = np.zeros((X.shape[0], a.shape[0]))
    alpha[0, :] = initial_distribution * b[:, X[0] - 1]

    for t in range(1, X.shape[0]):
        for j in range(a.shape[0]):
            alpha[t, j] = alpha[t - 1].T.dot(a[j, :]) * b[j, V[t] - 1]

    return alpha


alpha = forward(X, a, b, initial_distribution)
print(alpha)
