"""
Absolute Moment method for Hurst parameter estimation.

Ported from MATLAB implementation by Chu Chen (Version 1.0, 03/10/2008).
"""

import numpy as np


def absval(sequence, isplot=False, moment=1):
    """
    Estimate the Hurst parameter of a given sequence with the absolute
    moment method.

    Parameters
    ----------
    sequence : array_like
        The input sequence for estimation.
    isplot : bool, optional
        Whether to display the plot (default False).
    moment : int, optional
        A positive integer used to calculate the absolute moment (default 1).

    Returns
    -------
    H : float
        The estimated Hurst parameter of the input sequence.
    """
    sequence = np.asarray(sequence, dtype=float).ravel()
    sequence = sequence - np.mean(sequence)
    N = len(sequence)
    mlarge = int(np.floor(N / 5))
    M = np.floor(np.logspace(0, np.log10(mlarge), 50)).astype(int)
    M = np.unique(M[M > 1])
    n = len(M)
    cut_min = int(np.ceil(n / 10)) - 1  # Convert to 0-based
    cut_max = int(np.floor(6 * n / 10))  # Exclusive end in Python

    A = np.zeros(n)
    for i in range(n):
        m = M[i]
        k = int(np.floor(N / m))
        if k == 0:
            continue
        matrix_sequence = sequence[:m * k].reshape((m, k), order='F')
        A[i] = np.sum(np.abs(np.mean(matrix_sequence, axis=0)) ** moment) / k

    x = np.log10(M.astype(float))
    y = np.log10(A)
    y1 = -x + y[0] + x[0]
    X = x[cut_min:cut_max]
    Y = y[cut_min:cut_max]

    p1 = np.polyfit(X, Y, 1)
    Yfit = np.polyval(p1, X)
    yfit = np.polyval(p1, x)
    alpha = (Yfit[-1] - Yfit[0]) / (X[-1] - X[0])
    H = 1 + alpha / moment

    if isplot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x, y, 'b*')
        ax.plot(x, y1)
        ax.plot(X, Yfit, 'r-', linewidth=2)
        ax.plot(x[:cut_min], yfit[:cut_min], 'r:', linewidth=2)
        ax.plot(x[cut_max:], yfit[cut_max:], 'r:', linewidth=2)
        ax.set_xlabel('Log of Aggregate Level')
        ax.set_ylabel('Log of Absolute Moment')
        ax.set_title('Absolute Moment Method')
        plt.show()

    return H
