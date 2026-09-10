"""
Boxed Periodogram method for Hurst parameter estimation.

Ported from MATLAB implementation by Chu Chen (Version 1.0, 03/10/2008).
"""

import numpy as np


def boxper(sequence, isplot=False, boxnumber=50):
    """
    Estimate the Hurst parameter of a given sequence with the modified
    (boxed) periodogram method.

    Parameters
    ----------
    sequence : array_like
        The input sequence for estimation.
    isplot : bool, optional
        Whether to display the plot (default False).
    boxnumber : int, optional
        Number of boxes for the periodogram (default 50, must be in [30, 100]).

    Returns
    -------
    H : float
        The estimated Hurst parameter of the input sequence.
    """
    if boxnumber < 30 or boxnumber > 100:
        raise ValueError(
            'The input argument boxnumber must be an integer between [30, 100]')

    sequence = np.asarray(sequence, dtype=float).ravel()
    n = len(sequence)
    Xk = np.fft.fft(sequence)
    P_origin = np.abs(Xk) ** 2 / (2 * np.pi * n)
    P = P_origin[:int(np.floor(n / 2)) + 1]

    cut_min = int(np.ceil(0.001 * n / 2))
    M = np.floor(np.logspace(
        np.log10(cut_min),
        np.log10(0.1 * n - cut_min),
        boxnumber + 1
    )).astype(int)
    M = np.unique(M)
    N_boxes = len(M) - 1

    x = np.zeros(N_boxes)
    y = np.zeros(N_boxes)
    for i in range(N_boxes):
        m1 = M[i] + cut_min
        m2 = M[i + 1] + cut_min
        x[i] = np.log10((np.pi * (m2 - m1)) / n)
        # MATLAB P(m1:m2) is 1-based inclusive; Python equivalent: P[m1-1:m2]
        p_slice = P[(m1 - 1):m2]
        y[i] = np.log10(np.sum(p_slice) / len(p_slice))

    X = x
    Y = y
    p1 = np.polyfit(X, Y, 1)
    Yfit = np.polyval(p1, X)
    H = (1 - (Yfit[-1] - Yfit[0]) / (X[-1] - X[0])) / 2

    if isplot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x, y, 'b.', markersize=4)
        ax.plot(X, Yfit, 'r-', linewidth=2)
        ax.set_xlabel('Log10(Frequency)')
        ax.set_ylabel('Log10(Periodogram)')
        ax.set_title('Boxed Periodogram Method')
        plt.show()

    return H
