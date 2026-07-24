"""
Rescaled Range (R/S) method for Hurst parameter estimation.

Ported from MATLAB implementation by Chu Chen (Version 1.0, 03/10/2008).
"""

import numpy as np


def rs(sequence, isplot=False):
    """
    Estimate the Hurst parameter of a given sequence with the R/S method.

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
    N = len(sequence)
    dlarge = int(np.floor(N / 5))
    dsmall = max(10, np.log10(N) ** 2)
    D = np.floor(np.logspace(np.log10(dsmall), np.log10(dlarge), 50)).astype(int)
    D = np.unique(D)
    n = len(D)
    x = np.zeros(n)
    y = np.zeros(n)

    R_all = [None] * n
    S_all = [None] * n

    for i in range(n):
        d = D[i]
        m = int(np.floor(N / d))
        if m == 0:
            continue
        # Reshape sequence into d x m matrix (column-major like MATLAB)
        matrix_sequence = sequence[:d * m].reshape((d, m), order='F')

        Z1 = np.cumsum(matrix_sequence, axis=0)
        col_means = np.mean(matrix_sequence, axis=0)
        Z2 = np.cumsum(np.tile(col_means, (d, 1)), axis=0)
        diff = Z1 - Z2
        R_i = np.max(diff, axis=0) - np.min(diff, axis=0)
        S_i = np.std(matrix_sequence, axis=0, ddof=1)

        R_all[i] = R_i
        S_all[i] = S_i

        if np.min(R_i) == 0 or np.min(S_i) == 0:
            continue

        x[i] = np.log10(d)
        y[i] = np.mean(np.log10(R_i / S_i))

    # Fit a line with middle part of sequence
    index = x != 0
    x = x[index]
    y = y[index]
    n2 = len(x)
    cut_min = int(np.ceil(3 * n2 / 10)) - 1  # Convert to 0-based
    cut_max = int(np.floor(9 * n2 / 10))  # Exclusive end in Python

    X = x[cut_min:cut_max]
    Y = y[cut_min:cut_max]
    p1 = np.polyfit(X, Y, 1)
    Yfit = np.polyval(p1, X)
    H = (Yfit[-1] - Yfit[0]) / (X[-1] - X[0])

    if isplot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        bound = int(np.ceil(np.log10(N)))
        ax.set_xlim(0, bound)
        ax.set_ylim(0, 0.75 * bound)

        # Build the index mapping for valid entries
        temp = np.arange(1, len(index) + 1) * index
        valid_indices = temp[index].astype(int) - 1  # Convert to 0-based

        for i in range(n2):
            idx = valid_indices[i]
            if R_all[idx] is not None and S_all[idx] is not None:
                ratio = R_all[idx] / S_all[idx]
                # Filter out zero or negative ratios
                valid = ratio > 0
                if np.any(valid):
                    ax.plot(np.full(np.sum(valid), x[i]),
                            np.log10(ratio[valid]), 'b.', markersize=2)

        xline = np.linspace(0, bound, 10)
        y1 = 0.5 * xline
        y2 = xline
        ax.plot(xline, y1, 'b--', linewidth=2, label='slope 1/2')
        ax.plot(xline, y2, 'b-.', linewidth=2, label='slope 1')
        ax.plot(X, Yfit, 'r-', linewidth=3)
        ax.set_xlabel('log10(blocks of size m)')
        ax.set_ylabel('log10(R/S)')
        ax.set_title('R/S Method')
        plt.show()

    return H
