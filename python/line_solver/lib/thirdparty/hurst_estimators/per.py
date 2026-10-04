"""
Periodogram method for Hurst parameter estimation.

Ported from MATLAB implementation by Chu Chen (Version 1.0, 03/10/2008).
"""

import numpy as np


def per(sequence, isplot=False):
    """
    Estimate the Hurst parameter of a given sequence with the periodogram
    method.

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
    n = len(sequence)
    Xk = np.fft.fft(sequence)
    P_origin = np.abs(Xk) ** 2 / (2 * np.pi * n)
    P = P_origin[:int(np.floor(n / 2)) + 1]

    # MATLAB: x = log10((pi/n)*[2:floor(0.5*n)])
    # MATLAB: y = log10(P(2:floor(0.5*n)))
    # MATLAB 2:floor(0.5*n) is 1-based, giving indices 2,3,...,floor(0.5*n)
    # In Python (0-based), P[1:floor(0.5*n)] gives the same elements
    half_n = int(np.floor(0.5 * n))
    freq_indices = np.arange(2, half_n + 1)  # MATLAB [2:floor(0.5*n)]
    x = np.log10((np.pi / n) * freq_indices)
    y = np.log10(P[1:half_n])  # P(2:floor(0.5*n)) in MATLAB = P[1:half_n] in Python

    # Use the lowest 20% part of periodogram to estimate
    n_fit = int(np.floor(len(x) / 5))
    X = x[:n_fit]
    Y = y[:n_fit]
    p1 = np.polyfit(X, Y, 1)
    Yfit = np.polyval(p1, X)
    H = (1 - (Yfit[-1] - Yfit[0]) / (X[-1] - X[0])) / 2

    if isplot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x, y, 'b.', markersize=2)
        ax.plot(X, Yfit, 'r-', linewidth=3)
        ax.set_xlabel('log10(Frequency)')
        ax.set_ylabel('log10(Periodogram)')
        ax.set_title('Periodogram Method')
        plt.show()

    return H
