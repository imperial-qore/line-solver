"""
Peng (Residuals of Regression) method for Hurst parameter estimation.

Ported from MATLAB implementation by Chu Chen (Version 1.0, 03/10/2008).
"""

import numpy as np


def peng(sequence, isplot=False):
    """
    Estimate the Hurst parameter of a given sequence with the residuals
    of regression (Peng) method.

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
    FBM = np.cumsum(sequence)
    mlarge = int(np.floor(N / 5))
    msmall = max(10, np.log10(N) ** 2)
    M = np.floor(np.logspace(np.log10(msmall), np.log10(mlarge), 50)).astype(int)
    M = np.unique(M)
    n = len(M)
    cut_min = int(np.ceil(n / 10)) - 1  # Convert to 0-based
    cut_max = int(np.floor(7 * n / 10))  # Exclusive end in Python

    # Calculate residuals under different aggregate levels
    Goble_residuals = np.zeros(n)
    for i in range(n):
        m = M[i]
        k = int(np.floor(N / m))
        if k == 0:
            continue
        matrix_FBM = FBM[:m * k].reshape((m, k), order='F')
        x_reg = np.arange(1, m + 1, dtype=float)
        Local_residual = np.zeros(k)

        for j in range(k):
            y_col = matrix_FBM[:, j]
            # MATLAB: vv = [x' ones(length(x),1)]; p = vv\y;
            vv = np.column_stack([x_reg, np.ones(m)])
            p, _, _, _ = np.linalg.lstsq(vv, y_col, rcond=None)
            norm_xx = np.linalg.norm(y_col - vv @ p)
            Local_residual[j] = norm_xx ** 2 / m

        Goble_residuals[i] = np.mean(Local_residual)

    # Fit and calculate H
    x = np.log10(M.astype(float))
    y = np.log10(Goble_residuals)
    X = x[cut_min:cut_max]
    Y = y[cut_min:cut_max]
    p1 = np.polyfit(X, Y, 1)
    Yfit = np.polyval(p1, X)
    yfit = np.polyval(p1, x)
    H = 0.5 * (Yfit[-1] - Yfit[0]) / (X[-1] - X[0])

    if isplot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x, y, 'b*')
        ax.plot(X, Yfit, 'r-', linewidth=2)
        ax.plot(x[:cut_min], yfit[:cut_min], 'r:', linewidth=2)
        ax.plot(x[cut_max:], yfit[cut_max:], 'r:', linewidth=2)
        ax.set_xlabel('Log of Aggregate Level')
        ax.set_ylabel('Log of Residual Variance')
        ax.set_title('Peng Method')
        plt.show()

    return H
