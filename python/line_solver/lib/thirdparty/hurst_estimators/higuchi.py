"""
Higuchi's method for Hurst parameter estimation.

Ported from MATLAB implementation by Chu Chen (Version 1.0, 03/10/2008).
"""

import numpy as np


def higuchi(sequence, isplot=False):
    """
    Estimate the Hurst parameter of a given sequence with Higuchi's method.

    Parameters
    ----------
    sequence : array_like
        The input sequence for estimation.
    isplot : bool, optional
        Whether to display the plot (default False).

    Returns
    -------
    H : float
        The estimated Hurst parameter of the input sequence.
    """
    sequence = np.asarray(sequence, dtype=float).ravel()
    sequence = np.cumsum(sequence)
    N = len(sequence)
    mlarge = int(np.floor(N / 5))
    M = np.floor(np.logspace(0, np.log10(mlarge), 50)).astype(int)
    M = np.unique(M[M > 1])
    n = len(M)
    cut_min = int(np.ceil(n / 10)) - 1  # Convert to 0-based
    cut_max = int(np.floor(6 * n / 10))  # Exclusive end in Python

    curve_length = np.zeros(n)
    for h in range(n):
        m = M[h]
        k = int(np.floor((N - m) / m))
        if k == 0:
            continue
        temp_length = np.zeros((m, k))

        for i in range(m):
            for j in range(k):
                # MATLAB 1-based: sequence(i+j*m) - sequence(i+(j-1)*m)
                # Python 0-based: sequence[i+(j+1)*m] - sequence[i+j*m]
                temp_length[i, j] = abs(
                    sequence[i + (j + 1) * m] - sequence[i + j * m])

        curve_length[h] = np.sum(np.mean(temp_length, axis=1)) * ((N - 1) / m ** 3)

    x = np.log(M.astype(float))
    y = np.log(curve_length)
    X = x[cut_min:cut_max]
    Y = y[cut_min:cut_max]
    p1 = np.polyfit(X, Y, 1)
    Yfit = np.polyval(p1, X)
    yfit = np.polyval(p1, x)
    H = 2 + (Yfit[-1] - Yfit[0]) / (X[-1] - X[0])

    if isplot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x, y, 'b*')
        ax.plot(X, Yfit, 'r-', linewidth=2)
        ax.plot(x[:cut_min], yfit[:cut_min], 'r:', linewidth=2)
        ax.plot(x[cut_max:], yfit[cut_max:], 'r:', linewidth=2)
        ax.set_xlabel('Log(Aggregate Level)')
        ax.set_ylabel('Log(Curve Length)')
        ax.set_title('Higuchi Method')
        plt.show()

    return H
