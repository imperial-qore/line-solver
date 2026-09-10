"""
Deterministic (Semi-Markov) Process functions.

Native Python implementations of DET functions from the KPC-Toolbox,
ported from MATLAB: matlab/lib/kpctoolbox/smp/det/

A DET process is represented as a list [D0, D1] of two numpy arrays,
following the MAP (Markovian Arrival Process) convention where:
  - D0: hidden transition rate matrix (diagonal entries are negative rates)
  - D1: visible (arrival) transition rate matrix
"""

import numpy as np
from numpy.linalg import inv, matrix_power
from math import factorial
from typing import List, Tuple, Optional, Union


def det_embedded(DET: list) -> np.ndarray:
    """
    Compute the embedded DTMC of a deterministic MAP.

    P = inv(-D0) * D1

    Args:
        DET: Deterministic MAP as [D0, D1], each an (n, n) numpy array.

    Returns:
        P: Embedded discrete-time transition probability matrix (n, n).
    """
    D0 = np.asarray(DET[0], dtype=np.float64)
    D1 = np.asarray(DET[1], dtype=np.float64)
    P = -inv(D0) @ D1
    return P


def det_moment(DET: list, kset: Union[int, List[int], np.ndarray]) -> np.ndarray:
    """
    Compute moments of a deterministic process.

    For each order k in kset, the k-th moment is:
        M_k = pi * inv(-D0)^k * e(n)

    where pi is the stationary distribution of the embedded DTMC
    and e(n) is a column vector of ones.

    Args:
        DET: Deterministic MAP as [D0, D1].
        kset: Moment order(s) to compute. Can be a single int or
              a list/array of ints.

    Returns:
        Array of moment values, one per element of kset.
        If kset is a scalar int, returns a 1-element array.
    """
    from .mc import dtmc_solve
    from .basic import e

    D0 = np.asarray(DET[0], dtype=np.float64)
    n = D0.shape[0]

    # Ensure kset is iterable
    if np.isscalar(kset):
        kset = [kset]

    P = det_embedded(DET)
    pi = dtmc_solve(P)  # 1-d vector
    inv_neg_D0 = inv(-D0)
    ones_col = e(n)  # (n, 1) column vector

    M = []
    for k in kset:
        # pi * inv(-D0)^k * e(n)
        inv_neg_D0_k = matrix_power(inv_neg_D0, k)
        # pi is 1-d (n,), inv_neg_D0_k is (n,n), ones_col is (n,1)
        val = pi @ inv_neg_D0_k @ ones_col
        M.append(float(val))

    return np.array(M)


def det_scv(DET: list) -> float:
    """
    Compute the squared coefficient of variation of a deterministic process.

    SCV = E[X^2] / E[X]^2 - 1

    Args:
        DET: Deterministic MAP as [D0, D1].

    Returns:
        Squared coefficient of variation (scalar).
    """
    E1 = det_moment(DET, 1)[0]
    E2 = det_moment(DET, 2)[0]
    return E2 / (E1 ** 2) - 1


def det_acf(DET: list, kset: Union[int, List[int], np.ndarray]) -> np.ndarray:
    """
    Compute the autocorrelation function of a deterministic process.

    For each lag k in kset:
        rho(k) = (pi * K * P^(k-1) * K * e(n) - E1^2) / (E2 - E1^2)

    where K(i,j) = (-D0(i,i))^{-1} * P(i,j), P is the embedded DTMC,
    pi its stationary distribution, E1 the first moment, E2 the second moment.

    Args:
        DET: Deterministic MAP as [D0, D1].
        kset: Lag(s) at which to compute ACF. Can be a single int
              or a list/array of ints.

    Returns:
        Array of ACF values, one per element of kset.
    """
    from .mc import dtmc_solve
    from .basic import e

    D0 = np.asarray(DET[0], dtype=np.float64)
    n = D0.shape[0]

    # Ensure kset is iterable
    if np.isscalar(kset):
        kset = [kset]

    P = det_embedded(DET)
    pi = dtmc_solve(P)  # 1-d vector
    ones_col = e(n)  # (n, 1)

    # Build K matrix: K(i,j) = (-D0(i,i))^{-1} * P(i,j)
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            K[i, j] = (1.0 / (-D0[i, i])) * P[i, j]

    E1 = det_moment(DET, 1)[0]
    E2 = det_moment(DET, 2)[0]
    variance = E2 - E1 ** 2

    M = []
    for k in kset:
        # pi * K * P^(k-1) * K * e(n)
        Pk_minus_1 = matrix_power(P, k - 1)
        val = pi @ K @ Pk_minus_1 @ K @ ones_col
        rho = (float(val) - E1 ** 2) / variance if variance > 0 else 0.0
        M.append(rho)

    return np.array(M)


def det_sample(DET: list, nSamples: int,
               initState: Optional[Union[int, np.ndarray]] = None,
               seed: Optional[int] = None
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate samples from a deterministic MAP.

    Simulates the MAP and returns inter-arrival times. The sampling
    is done in batches of 10000 for efficiency (matching MATLAB).

    Args:
        DET: Deterministic MAP as [D0, D1].
        nSamples: Number of inter-arrival time samples to generate.
        initState: Initial state (1-based index matching MATLAB convention),
                   or a probability vector to sample from, or None to use
                   the interval-stationary distribution (map_pie).
        seed: Random seed (optional, currently unused -- included for
              signature compatibility).

    Returns:
        Tuple of (SAMPLES, LAST, FIRST):
          - SAMPLES: (nSamples,) array of inter-arrival times
          - LAST: (nSamples,) array of last visited states (1-based)
          - FIRST: (nSamples,) array of first visited states (1-based)
    """
    from line_solver.api.mam.map_analysis import map_pie

    D0 = np.asarray(DET[0], dtype=np.float64)
    D1 = np.asarray(DET[1], dtype=np.float64)

    # Determine initial state
    if initState is None:
        # Interval-stationary initialization
        pi = map_pie(D0, D1)
        x = np.cumsum(pi)
        r = np.random.rand()
        initState = int(np.min(np.where(r <= x)[0])) + 1  # 1-based
    elif hasattr(initState, '__len__') and len(initState) > 1:
        # initState is a probability vector
        pi = np.asarray(initState, dtype=np.float64)
        x = np.cumsum(pi)
        r = np.random.rand()
        initState = int(np.min(np.where(r <= x)[0])) + 1  # 1-based

    # Run in batches of 10000 (matching MATLAB)
    RUNS = nSamples // 10000
    LS = initState  # last state carried between batches
    SAMPLES = []
    LAST = []
    FIRST = []

    for _ in range(RUNS):
        S, L, F, LS = _sub_map_sample(DET, 10000, LS)
        SAMPLES.append(S)
        LAST.append(L)
        FIRST.append(F)

    remainder = nSamples % 10000
    if remainder > 0:
        S, L, F, LS = _sub_map_sample(DET, remainder, LS)
        SAMPLES.append(S)
        LAST.append(L)
        FIRST.append(F)

    if len(SAMPLES) > 0:
        SAMPLES = np.concatenate(SAMPLES)
        LAST = np.concatenate(LAST)
        FIRST = np.concatenate(FIRST)
    else:
        SAMPLES = np.array([])
        LAST = np.array([], dtype=int)
        FIRST = np.array([], dtype=int)

    return SAMPLES, LAST, FIRST


def _sub_map_sample(DET: list, nSamples: int, initState: int
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Internal sub-routine for sampling from a MAP in a single batch.

    This mirrors the MATLAB sub_map_sample nested function.

    Args:
        DET: Deterministic MAP as [D0, D1].
        nSamples: Number of samples in this batch.
        initState: Current state (1-based).

    Returns:
        Tuple of (SAMPLES, LAST, FIRST, lastState):
          - SAMPLES: (nSamples,) inter-arrival times
          - LAST: (nSamples,) last state per sample (1-based)
          - FIRST: (nSamples,) first state per sample (1-based)
          - lastState: final state for chaining to next batch (1-based)
    """
    if nSamples <= 0:
        return (np.array([]), np.array([], dtype=int),
                np.array([], dtype=int), initState)

    D0 = np.asarray(DET[0], dtype=np.float64)
    D1 = np.asarray(DET[1], dtype=np.float64)
    nStates = D0.shape[0]

    # Build transition probability matrix p (nStates x 2*nStates)
    # Columns 0..nStates-1: hidden transitions (from D0)
    # Columns nStates..2*nStates-1: arrival transitions (from D1)
    p = np.zeros((nStates, 2 * nStates))
    for b in range(2):
        for i in range(nStates):
            for j in range(nStates):
                p[i, b * nStates + j] = DET[b][i, j] / abs(D0[i, i])

    # Zero out diagonal of hidden transitions (no self-loop without arrival)
    for i in range(nStates):
        p[i, i] = 0.0

    # Compute cumulative distribution for each row
    cdf = np.zeros_like(p)
    for i in range(nStates):
        cdf[i, :] = np.cumsum(p[i, :])
    cdf = np.abs(cdf)

    # Holding times: -1/diag(D0)
    holdTimes = -1.0 / np.diag(D0)

    curState = initState  # 1-based
    # visits tracks the state path for each sample; dynamically grown
    maxpathlen = 20
    visits = np.zeros((nSamples, maxpathlen), dtype=int)

    for i in range(nSamples):
        arrival = False
        last = 1  # next write position (0-based index into visits[i,:])
        visits[i, 0] = curState  # 1-based state
        while not arrival:
            r = np.random.rand()
            # Find first column where cdf >= r (1-based state indexing)
            destState = int(np.argmax(cdf[curState - 1, :] >= r)) + 1  # 1-based
            if destState > nStates:
                # Arrival transition
                arrival = True
                destState = destState - nStates
                curState = destState
            else:
                # Hidden transition
                if last >= maxpathlen:
                    maxpathlen = maxpathlen + 5
                    visits = np.pad(visits,
                                    ((0, 0), (0, 5)),
                                    mode='constant', constant_values=0)
                visits[i, last] = destState
                curState = destState
                last += 1

    # Compute LAST: last non-zero state in each row of visits
    LAST = np.zeros(nSamples, dtype=int)
    for i in range(nSamples):
        nonzero_idx = np.where(visits[i, :] > 0)[0]
        LAST[i] = visits[i, nonzero_idx[-1]]

    # Replace state indices in visits with holding times
    H = visits.astype(np.float64)
    for s in range(1, nStates + 1):
        H[visits == s] = holdTimes[s - 1]

    # Sum holding times across path for each sample
    SAMPLES = np.sum(H, axis=1)

    FIRST = visits[:, 0].copy()

    return SAMPLES, LAST, FIRST, curState


def det_sum(DETs: list) -> list:
    """
    Compute the sum of independent deterministic processes.

    Delegates to map_sumind which constructs the MAP representing
    the sum (concatenation) of independent random variables.

    Args:
        DETs: List of DET processes, each as [D0, D1].
              For example: [[D0_1, D1_1], [D0_2, D1_2], ...]

    Returns:
        Combined DET as [D0, D1].
    """
    from line_solver.api.mam.map_analysis import map_sumind

    # Convert each DET [D0, D1] to (D0, D1) tuple as expected by map_sumind
    maps = []
    for det in DETs:
        D0 = np.asarray(det[0], dtype=np.float64)
        D1 = np.asarray(det[1], dtype=np.float64)
        maps.append((D0, D1))

    D0_new, D1_new = map_sumind(maps)
    return [D0_new, D1_new]
