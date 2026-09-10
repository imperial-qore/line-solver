"""
Numerical inverse Laplace transform using Abate-Whitt framework.

Implements three methods:
  - CME (Concentrated Matrix Exponential / Talbot contour): uses pre-computed
    optimal parameters from iltcme.json to select abscissae and weights on a
    Talbot-style contour, yielding very accurate inversions with few function
    evaluations.
  - Euler: uses the Euler summation acceleration of the Bromwich integral with
    binomial-coefficient weights and geometric abscissae, following the
    Abate-Whitt (1995) formulation.
  - Gaver-Stehfest: uses the Gaver functional with Stehfest's factorial-based
    weights and logarithmic abscissae.

All three methods share a common final computation step that evaluates the
Laplace-domain function on a meshgrid of abscissae/time values and performs a
weighted real-part summation.

This is a Python port of the MATLAB ``matlab_ilt.m`` function from the iltcme
library.
"""

import json
import os
import math
import numpy as np

__all__ = ["matlab_ilt", "matlab_ilt_matrix", "cme_parameters"]

# ---------------------------------------------------------------------------
# Module-level cache for the CME parameter table (loaded once on first use)
# ---------------------------------------------------------------------------
_cme_params = None


def _load_cme_params():
    """Load the pre-computed CME/Talbot contour parameters from iltcme.json."""
    global _cme_params
    if _cme_params is None:
        json_path = os.path.join(os.path.dirname(__file__), "iltcme.json")
        with open(json_path, "r") as fh:
            _cme_params = json.load(fh)
    return _cme_params


def cme_parameters():
    """Return the pre-computed CME parameter table from ``iltcme.json``.

    The table is a list of dictionaries with keys ``n`` (number of harmonic
    terms), ``optim``, ``a``, ``b``, ``c``, ``omega``, ``mu1`` and ``cv2``. It
    is loaded once and cached, and is shared by the inverse Laplace transform
    above and by the CME distribution class, which reads the same entries as a
    matrix-exponential representation of order ``2*n+1``.

    The returned list is the cached object, so callers must not modify it.
    """
    return _load_cme_params()


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------


def matlab_ilt(fun, T, maxFnEvals, method="cme"):
    """Numerically invert a Laplace transform.

    Parameters
    ----------
    fun : callable
        A function ``fun(s)`` that evaluates the Laplace-domain transform at a
        (possibly complex) point *s*.  Must accept and return scalar or
        array-like values.
    T : array_like
        One-dimensional array of positive real time points at which the
        inverse transform is evaluated.
    maxFnEvals : int
        Budget of Laplace-domain function evaluations.  For the CME method
        this determines which pre-computed parameter set is selected (the one
        with smallest ``cv2`` whose ``n + 1 <= maxFnEvals``).  For Euler /
        Gaver-Stehfest it controls the number of summation terms.
    method : str, optional
        Inversion method: ``'cme'`` (default), ``'euler'``, or ``'gaver'``.

    Returns
    -------
    ilt : numpy.ndarray
        Real-valued array of the same length as *T* containing the
        approximate inverse Laplace transform values.
    """

    T = np.atleast_1d(np.asarray(T, dtype=float))

    if method == "cme":
        eta, beta = _cme_weights(maxFnEvals)
    elif method == "euler":
        eta, beta = _euler_weights(maxFnEvals)
    elif method == "gaver":
        eta, beta = _gaver_weights(maxFnEvals)
    else:
        raise ValueError(
            "Unknown inverse Laplace transform method '{}'. "
            "Supported: cme, euler, gaver".format(method)
        )

    # ------------------------------------------------------------------
    # Common Abate-Whitt summation (shared by all three methods)
    #
    #   f(t) ≈ (1/t) * sum_k  Re[ eta_k * F(beta_k / t) ]
    #
    # Vectorised over all time points T simultaneously.
    # ------------------------------------------------------------------
    eta = np.asarray(eta, dtype=complex)
    beta = np.asarray(beta, dtype=complex)

    # eta_mesh has shape (len(T), len(eta)); rows replicate eta
    # beta_mesh has the same shape; rows replicate beta
    eta_mesh, T_mesh = np.meshgrid(eta, T)
    beta_mesh = np.meshgrid(beta, T)[0]

    # Evaluate the Laplace-domain function at beta/T for every combination
    s_vals = beta_mesh / T_mesh
    F_vals = np.vectorize(fun)(s_vals)

    ilt = (1.0 / T) * np.sum(np.real(eta_mesh * F_vals), axis=1)
    return ilt


def matlab_ilt_matrix(fun, T, maxFnEvals, method="cme"):
    """Numerically invert a matrix-valued Laplace transform.

    Like :func:`matlab_ilt`, but ``fun(s)`` returns a 2-D array (matrix) rather
    than a scalar. ``fun`` is evaluated once per Abate-Whitt node and the full
    complex matrix is accumulated, sharing the eta/beta weights with the scalar
    variant.

    Parameters
    ----------
    fun : callable
        ``fun(s)`` returning an ``(nr, nc)`` complex array at complex point *s*.
    T : array_like
        One-dimensional array of positive real time points.
    maxFnEvals : int
        Budget of Laplace-domain function evaluations.
    method : str, optional
        ``'cme'`` (default), ``'euler'``, or ``'gaver'``.

    Returns
    -------
    vals : numpy.ndarray
        Real-valued array of shape ``(len(T), nr, nc)``.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))

    if method == "cme":
        eta, beta = _cme_weights(maxFnEvals)
    elif method == "euler":
        eta, beta = _euler_weights(maxFnEvals)
    elif method == "gaver":
        eta, beta = _gaver_weights(maxFnEvals)
    else:
        raise ValueError(
            "Unknown inverse Laplace transform method '{}'. "
            "Supported: cme, euler, gaver".format(method)
        )

    eta = np.asarray(eta, dtype=complex)
    beta = np.asarray(beta, dtype=complex)

    # Probe output dimensions with one cheap evaluation.
    probe = np.asarray(fun(beta[0] / T[0]), dtype=complex)
    nr, nc = probe.shape

    vals = np.zeros((len(T), nr, nc), dtype=float)
    for i in range(len(T)):
        t = T[i]
        acc = np.zeros((nr, nc), dtype=complex)
        for k in range(len(eta)):
            acc += eta[k] * np.asarray(fun(beta[k] / t), dtype=complex)
        vals[i, :, :] = np.real(acc) / t
    return vals


# ---------------------------------------------------------------------------
# CME / Talbot method
# ---------------------------------------------------------------------------


def _cme_weights(maxFnEvals):
    """Compute eta and beta vectors for the CME (Talbot contour) method.

    Selects the parameter set from ``iltcme.json`` with the smallest ``cv2``
    (i.e. steepest contour / highest accuracy) whose number of function
    evaluations ``n + 1`` does not exceed *maxFnEvals*.

    Returns
    -------
    eta : numpy.ndarray (complex)
    beta : numpy.ndarray (complex)
    """
    cme_params = _load_cme_params()

    # Find the best-fitting parameter set ----------------------------------
    # Start with the first entry; then scan for one with lower cv2 that fits
    # within the evaluation budget.
    best = cme_params[0]
    for entry in cme_params[1:]:
        if entry["cv2"] < best["cv2"] and entry["n"] + 1 <= maxFnEvals:
            best = entry

    params = best

    a = np.asarray(params["a"], dtype=float)
    b = np.asarray(params["b"], dtype=float)
    c = float(params["c"])
    mu1 = float(params["mu1"])
    omega = float(params["omega"])
    n = int(params["n"])

    # eta = [c*mu1,  (a + i*b)*mu1]          length n+1
    eta_0 = c * mu1
    eta_rest = (a + 1j * b) * mu1
    eta = np.concatenate(([eta_0], eta_rest))

    # beta = [mu1,  (1 + i*k*omega)*mu1]    k = 1..n      length n+1
    k = np.arange(1, n + 1)
    beta_0 = mu1
    beta_rest = (1.0 + 1j * k * omega) * mu1
    beta = np.concatenate(([beta_0], beta_rest))

    return eta, beta


# ---------------------------------------------------------------------------
# Euler method
# ---------------------------------------------------------------------------


def _euler_weights(maxFnEvals):
    """Compute eta and beta vectors for the Euler summation method.

    Uses binomial-coefficient weights computed in log-space to avoid overflow,
    with a scaling factor of ``10^(n_euler/3)`` and geometric abscissae along
    the imaginary axis.

    Returns
    -------
    eta : numpy.ndarray (complex)
    beta : numpy.ndarray (complex)
    """
    n_euler = int(math.floor((maxFnEvals - 1) / 2))

    # Build the raw eta vector of length 2*n_euler + 1
    # eta = [0.5,  1, 1, ..., 1 (n_euler ones),  0, 0, ..., 0 (n_euler-1 zeros),  2^{-n_euler}]
    eta = np.zeros(2 * n_euler + 1)
    eta[0] = 0.5
    eta[1 : n_euler + 1] = 1.0
    # eta[n_euler+1 : 2*n_euler] are already 0.0
    eta[2 * n_euler] = 2.0 ** (-n_euler)

    # Accumulate binomial coefficients (in reverse) using log-space arithmetic
    # Matches MATLAB:
    #   for k = 1:n_euler-1
    #       eta(2*n_euler - k + 1) = eta(2*n_euler - k + 2)
    #           + exp(sum(log(1:n_euler)) - n_euler*log(2)
    #                 - sum(log(1:k)) - sum(log(1:(n_euler-k))))
    #   end
    log_n_fact = _log_factorial(n_euler)
    for k in range(1, n_euler):
        binom_log = (
            log_n_fact
            - n_euler * math.log(2)
            - _log_factorial(k)
            - _log_factorial(n_euler - k)
        )
        # MATLAB index: 2*n_euler - k + 1  (1-based) → Python: 2*n_euler - k
        # MATLAB index: 2*n_euler - k + 2  (1-based) → Python: 2*n_euler - k + 1
        # But note MATLAB's eta has 1-based indexing; translating:
        #   MATLAB eta(2*n_euler - k + 1) ↔ Python eta[2*n_euler - k]
        #   MATLAB eta(2*n_euler - k + 2) ↔ Python eta[2*n_euler - k + 1]
        eta[2 * n_euler - k] = eta[2 * n_euler - k + 1] + math.exp(binom_log)

    # k = 0 : 2*n_euler
    k = np.arange(0, 2 * n_euler + 1)

    # beta = n_euler * log(10)/3 + i*pi*k
    beta = n_euler * math.log(10) / 3.0 + 1j * math.pi * k

    # eta = 10^(n_euler/3) * (-1)^k .* eta
    sign = 1.0 - 2.0 * np.mod(k, 2)  # +1 for even k, -1 for odd k
    eta = (10.0 ** (n_euler / 3.0)) * sign * eta

    return eta, beta


# ---------------------------------------------------------------------------
# Gaver-Stehfest method
# ---------------------------------------------------------------------------


def _gaver_weights(maxFnEvals):
    """Compute eta and beta vectors for the Gaver-Stehfest method.

    Uses factorial-based weights computed in log-space and purely real
    logarithmic abscissae.

    Returns
    -------
    eta : numpy.ndarray (complex)
    beta : numpy.ndarray (complex)
    """
    # Ensure maxFnEvals is even
    if maxFnEvals % 2 == 1:
        maxFnEvals = maxFnEvals - 1

    ndiv2 = maxFnEvals // 2

    eta = np.zeros(maxFnEvals, dtype=complex)
    beta = np.zeros(maxFnEvals, dtype=complex)

    ln2 = math.log(2.0)

    for k in range(1, maxFnEvals + 1):  # k = 1 .. maxFnEvals
        inside_sum = 0.0
        j_start = int(math.floor((k + 1) / 2))
        j_end = min(k, ndiv2)

        for j in range(j_start, j_end + 1):
            # Compute the summand in log-space to avoid overflow:
            #   j^(ndiv2+1) * C(ndiv2,j) * C(2j,j) * C(j,k-j)  / factorial(ndiv2)
            #
            # MATLAB uses:
            #   exp((ndiv2+1)*log(j)
            #       - sum(log(1:(ndiv2-j)))
            #       + sum(log(1:2*j))
            #       - 2*sum(log(1:j))
            #       - sum(log(1:(k-j)))
            #       - sum(log(1:(2*j-k))))
            #
            # Note: MATLAB's formula already incorporates the factorial(ndiv2)
            # denominator into the log-sum decomposition.
            log_term = (
                (ndiv2 + 1) * math.log(j)
                - _log_factorial(ndiv2 - j)
                + _log_factorial(2 * j)
                - 2.0 * _log_factorial(j)
                - _log_factorial(k - j)
                - _log_factorial(2 * j - k)
            )
            inside_sum += math.exp(log_term)

        eta[k - 1] = ln2 * ((-1.0) ** (k + ndiv2)) * inside_sum
        beta[k - 1] = k * ln2

    return eta, beta


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def _log_factorial(n):
    """Compute log(n!) = sum(log(1), log(2), ..., log(n)).

    Returns 0.0 for n <= 0, consistent with MATLAB's ``sum(log(1:0)) == 0``.
    """
    if n <= 0:
        return 0.0
    return sum(math.log(i) for i in range(1, n + 1))
