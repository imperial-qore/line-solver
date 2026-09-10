"""
MAP and PH Queue Analysis.

Native Python implementations for analyzing queues with Markovian Arrival
Processes (MAP) and Phase-Type (PH) service distributions using BuTools.

Key functions:
    qsys_phph1: PH/PH/1 queue
    qsys_mapph1: MAP/PH/1 queue
    qsys_mapm1: MAP/M/1 queue
    qsys_mapmc: MAP/M/c queue
    qsys_mapmap1: MAP/MAP/1 queue

References:
    Original MATLAB: matlab/src/api/qsys/qsys_*.m
    BuTools library for matrix-analytic queue analysis
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple
from dataclasses import dataclass


@dataclass
class QueueResult:
    """Result structure for queue analysis."""
    meanQueueLength: float
    meanWaitingTime: float
    meanSojournTime: float
    utilization: float
    queueLengthDist: Optional[np.ndarray] = None
    queueLengthMoments: Optional[np.ndarray] = None
    sojournTimeMoments: Optional[np.ndarray] = None
    analyzer: str = "native"


def ph_to_map(alpha: np.ndarray, T: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a PH distribution to its equivalent MAP representation.

    For a PH renewal process, the MAP has:
        D0 = T (transitions within the PH, no arrival)
        D1 = t * alpha where t = -T*e (exit rates times restart distribution)

    Args:
        alpha: Initial probability vector (1 x n)
        T: Sub-generator matrix (n x n)

    Returns:
        Tuple of (D0, D1):
            D0: MAP hidden transition matrix
            D1: MAP observable transition matrix
    """
    alpha = np.asarray(alpha, dtype=float).flatten()
    T = np.atleast_2d(np.asarray(T, dtype=float))

    n = T.shape[0]

    # D0 = T (hidden transitions)
    D0 = T.copy()

    # Exit rate vector: t = -T * e
    exit_rates = -T @ np.ones(n)

    # D1 = t * alpha (restart to initial distribution)
    D1 = np.outer(exit_rates, alpha)

    return D0, D1


def qsys_phph1(alpha: np.ndarray, T: np.ndarray,
               beta: np.ndarray, S: np.ndarray,
               numQLMoms: int = 3,
               numQLProbs: int = 100,
               numSTMoms: int = 3) -> QueueResult:
    """
    Analyze a PH/PH/1 queue using matrix-analytic methods.

    Converts the arrival PH to MAP representation and uses the
    MMAPPH1FCFS solver from BuTools.

    Args:
        alpha: Arrival PH initial probability vector (1 x n)
        T: Arrival PH generator matrix (n x n)
        beta: Service PH initial probability vector (1 x m)
        S: Service PH generator matrix (m x m)
        numQLMoms: Number of queue length moments to compute (default: 3)
        numQLProbs: Number of queue length probabilities (default: 100)
        numSTMoms: Number of sojourn time moments (default: 3)

    Returns:
        QueueResult with queue performance metrics

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_phph1.m
    """
    alpha = np.asarray(alpha, dtype=float).flatten()
    T = np.atleast_2d(np.asarray(T, dtype=float))
    beta = np.asarray(beta, dtype=float).flatten()
    S = np.atleast_2d(np.asarray(S, dtype=float))

    # Convert arrival PH to MAP
    D0, D1 = ph_to_map(alpha, T)

    # Compute arrival rate
    negTinv = np.linalg.inv(-T)
    mean_interarrival = alpha @ negTinv @ np.ones(T.shape[0])
    lambda_val = 1.0 / mean_interarrival

    # Compute service rate
    negSinv = np.linalg.inv(-S)
    mean_service = beta @ negSinv @ np.ones(S.shape[0])
    mu = 1.0 / mean_service

    rho = lambda_val / mu

    # Check stability
    if rho >= 1:
        return QueueResult(
            meanQueueLength=np.inf,
            meanWaitingTime=np.inf,
            meanSojournTime=np.inf,
            utilization=rho,
            analyzer="native:unstable"
        )

    try:
        # see _kb/03-api-layer.md for rationale
        from line_solver.lib.thirdparty.butools.queues import MMAPPH1FCFS

        D = [np.matrix(D0), np.matrix(D1)]
        ncMoms, ncDistr, stMoms = MMAPPH1FCFS(
            D, [np.matrix(beta)], [np.matrix(S)],
            'ncMoms', numQLMoms, 'ncDistr', numQLProbs, 'stMoms', numSTMoms
        )

        meanQL = ncMoms[0] if len(ncMoms) > 0 else 0.0
        meanST = stMoms[0] if len(stMoms) > 0 else 0.0
        meanWT = max(0, meanST - mean_service)

        return QueueResult(
            meanQueueLength=meanQL,
            meanWaitingTime=meanWT,
            meanSojournTime=meanST,
            utilization=rho,
            queueLengthDist=ncDistr,
            queueLengthMoments=ncMoms,
            sojournTimeMoments=stMoms,
            analyzer="BuTools:MMAPPH1FCFS"
        )

    except ImportError:
        # Fallback to approximation
        # Use M/M/1 approximation with adjusted parameters
        meanQL = rho / (1 - rho)
        meanWT = rho / (lambda_val * (1 - rho))
        meanST = meanWT + mean_service

        return QueueResult(
            meanQueueLength=meanQL,
            meanWaitingTime=meanWT,
            meanSojournTime=meanST,
            utilization=rho,
            analyzer="native:MM1_approx"
        )


def qsys_mapph1(D0: np.ndarray, D1: np.ndarray,
                beta: np.ndarray, S: np.ndarray,
                numQLMoms: int = 3,
                numQLProbs: int = 100,
                numSTMoms: int = 3) -> QueueResult:
    """
    Analyze a MAP/PH/1 queue.

    Args:
        D0: MAP hidden transition matrix
        D1: MAP observable transition matrix
        beta: Service PH initial probability vector
        S: Service PH generator matrix
        numQLMoms: Number of queue length moments
        numQLProbs: Number of queue length probabilities
        numSTMoms: Number of sojourn time moments

    Returns:
        QueueResult with queue performance metrics

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mapph1.m
    """
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    beta = np.asarray(beta, dtype=float).flatten()
    S = np.atleast_2d(np.asarray(S, dtype=float))

    # Compute arrival rate from MAP
    D = D0 + D1
    n = D.shape[0]

    # Stationary distribution of MAP
    pi = np.ones(n) / n
    for _ in range(1000):
        pi_new = pi @ np.linalg.matrix_power(np.eye(n) + D / 100, 100)
        pi_new /= np.sum(pi_new)
        if np.linalg.norm(pi_new - pi) < 1e-12:
            break
        pi = pi_new

    lambda_val = pi @ D1 @ np.ones(n)

    # Compute service rate
    negSinv = np.linalg.inv(-S)
    mean_service = beta @ negSinv @ np.ones(S.shape[0])
    mu = 1.0 / mean_service

    rho = lambda_val / mu

    if rho >= 1:
        return QueueResult(
            meanQueueLength=np.inf,
            meanWaitingTime=np.inf,
            meanSojournTime=np.inf,
            utilization=rho,
            analyzer="native:unstable"
        )

    try:
        # See qsys_phph1 above: mtfcfs (re-exported by the package, not cfcfs),
        # positional option pairs, np.matrix arguments.
        from line_solver.lib.thirdparty.butools.queues import MMAPPH1FCFS

        D = [np.matrix(D0), np.matrix(D1)]
        ncMoms, ncDistr, stMoms = MMAPPH1FCFS(
            D, [np.matrix(beta)], [np.matrix(S)],
            'ncMoms', numQLMoms, 'ncDistr', numQLProbs, 'stMoms', numSTMoms
        )

        meanQL = ncMoms[0] if len(ncMoms) > 0 else 0.0
        meanST = stMoms[0] if len(stMoms) > 0 else 0.0
        meanWT = max(0, meanST - mean_service)

        return QueueResult(
            meanQueueLength=meanQL,
            meanWaitingTime=meanWT,
            meanSojournTime=meanST,
            utilization=rho,
            queueLengthDist=ncDistr,
            queueLengthMoments=ncMoms,
            sojournTimeMoments=stMoms,
            analyzer="BuTools:MMAPPH1FCFS"
        )

    except ImportError:
        # Approximation
        meanQL = rho / (1 - rho)
        meanWT = rho / (lambda_val * (1 - rho))
        meanST = meanWT + mean_service

        return QueueResult(
            meanQueueLength=meanQL,
            meanWaitingTime=meanWT,
            meanSojournTime=meanST,
            utilization=rho,
            analyzer="native:MM1_approx"
        )


def qsys_mapm1(D0: np.ndarray, D1: np.ndarray,
               mu: float) -> QueueResult:
    """
    Analyze a MAP/M/1 queue.

    Args:
        D0: MAP hidden transition matrix
        D1: MAP observable transition matrix
        mu: Service rate

    Returns:
        QueueResult with queue performance metrics

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mapm1.m
    """
    # Service is exponential, so use PH with single phase
    beta = np.array([1.0])
    S = np.array([[-mu]])

    return qsys_mapph1(D0, D1, beta, S)


def _stat_gen(G: np.ndarray) -> np.ndarray:
    """
    Left null vector of G normalized to sum one, i.e. x G = 0 with x e = 1.

    This is SMCSolver's stat(A) with A = G + I: the augmented system [G, e] is
    solved directly, so G is NOT required to have zero row sums. The boundary
    generators of the Gaver-Jacobs-Latouche recursion below do not have them.
    """
    G = np.atleast_2d(np.asarray(G, dtype=float))
    n = G.shape[0]
    B = np.hstack([G, np.ones((n, 1))])
    y = np.zeros(n + 1)
    y[n] = 1.0
    x, _, _, _ = np.linalg.lstsq(B.T, y, rcond=None)
    return x


def qsys_mapmc(D0: np.ndarray, D1: np.ndarray,
               mu: float, c: int, max_num_comp: int = 1000) -> QueueResult:
    """
    Analyze a MAP/M/c queue.

    The chain is a level-dependent QBD in the number in system: above level c
    the c servers are all busy, so the blocks repeat as A0 = c*mu*I (down),
    A1 = D0 - c*mu*I (local), A2 = D1 (up) and the tail is matrix-geometric in
    R. Below level c the departure rate is level dependent and the boundary
    vector comes from the Gaver, Jacobs and Latouche backward recursion. The
    waiting time is phase type; its generator is the fixed point of a Sylvester
    equation, and the representation is the time reversal with respect to the
    arrival-epoch vector.

    Args:
        D0: MAP hidden transition matrix
        D1: MAP observable transition matrix
        mu: Service rate per server
        c: Number of servers
        max_num_comp: Cap on the number of queue length probabilities

    Returns:
        QueueResult with queue performance metrics

    References:
        Perez, Van Velthoven, Van Houdt, Q-MAM, ValueTools 2008 (Q_CT_MAP_M_C)
        Gaver, Jacobs, Latouche, Adv. Appl. Probab. 16:715-731, 1984
        Asmussen, Moller, Queueing Systems 37(1):9-29, 2001
        Original MATLAB: matlab/src/api/qsys/qsys_mapmc.m
    """
    from scipy.linalg import solve_sylvester
    from ..mam.qbd import qbd_R_logred

    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    c = int(c)
    m = D0.shape[0]
    eye_m = np.eye(m)

    theta = _stat_gen(D0 + D1)
    lambda_val = float(theta @ D1 @ np.ones(m))
    rho = lambda_val / (c * mu)

    if rho >= 1:
        return QueueResult(
            meanQueueLength=np.inf,
            meanWaitingTime=np.inf,
            meanSojournTime=np.inf,
            utilization=rho,
            analyzer="native:unstable"
        )

    A0 = c * mu * eye_m
    A1 = D0 - c * mu * eye_m
    A2 = D1
    R = qbd_R_logred(A0, A1, A2)

    # Gaver-Jacobs-Latouche boundary levels 0..c-1
    piGJL = np.zeros(c * m)
    if c > 1:
        invC = {}
        invC[1] = np.linalg.inv(-D0)
        for i in range(2, c):
            invC[i] = np.linalg.inv(-D0 + (i - 1) * mu * eye_m
                                    - (i - 1) * mu * invC[i - 1] @ D1)
        piGJL[(c - 1) * m:c * m] = _stat_gen(
            D0 - (c - 1) * mu * eye_m + R @ A0 + (c - 1) * mu * invC[c - 1] @ D1)
        for i in range(c - 1, 0, -1):
            piGJL[(i - 1) * m:i * m] = piGJL[i * m:(i + 1) * m] @ (i * mu * invC[i])
    else:
        piGJL[0:m] = _stat_gen(D0 + R @ A0)

    ImR_inv = np.linalg.inv(eye_m - R)
    norm_k = (np.sum(piGJL[0:(c - 1) * m])
              + float(piGJL[(c - 1) * m:c * m] @ ImR_inv @ np.ones(m)))
    piGJL = piGJL / norm_k

    # Matrix-geometric tail from level c upwards
    piC = [piGJL[(c - 1) * m:c * m].copy()]
    sumpi = float(np.sum(piGJL))
    while sumpi < 1 - 1e-10 and len(piC) < 1 + max_num_comp - c:
        piC.append(piC[-1] @ R)
        sumpi += float(np.sum(piC[-1]))
    piC = np.asarray(piC)

    ql = np.concatenate([
        np.array([np.sum(piGJL[i * m:(i + 1) * m]) for i in range(c - 1)]),
        np.sum(piC, axis=1)])
    meanQL = float(np.sum(np.arange(len(ql)) * ql))

    # Waiting time: phase-type, via the Sylvester fixed point
    piT = np.concatenate([piGJL[0:(c - 1) * m], piC.reshape(-1)])
    d1e = D1 @ np.ones(m)
    denom = float(np.sum(piT.reshape(-1, m) @ d1e))
    prob_zero = float(np.sum(piGJL.reshape(-1, m) @ d1e)) / denom

    temp = piGJL[(c - 1) * m:c * m] @ ImR_inv @ D1
    alpha_vec = temp / np.sum(temp)

    Tnew = -A0.copy()
    for _ in range(1000):
        Told = Tnew
        # Tnew * L + L * D0 = -I
        Lmat = solve_sylvester(Tnew, D0, -eye_m)
        Tnew = -A0 + Lmat @ D1 * (mu * c)
        if np.max(np.abs(Told - Tnew)) <= 1e-10:
            break

    rho_vec = np.sum(Tnew + A0, axis=1)
    nonz = np.flatnonzero(alpha_vec > 0)
    theta_red = alpha_vec[nonz]
    Tred = Tnew[np.ix_(nonz, nonz)]
    Smat = np.diag(1.0 / theta_red) @ Tred.T @ np.diag(theta_red)
    wait_alpha = alpha_vec * rho_vec / float(alpha_vec @ rho_vec)
    wait_alpha = (1 - prob_zero) * wait_alpha[nonz]

    # E[W] = alpha(-S)^-1 e. wait_alpha is DEFECTIVE by design: its deficit is
    # the atom at zero, i.e. the probability of finding a free server.
    meanWT = float(wait_alpha @ np.linalg.solve(-Smat, np.ones(len(nonz))))
    meanST = meanWT + 1.0 / mu

    return QueueResult(
        meanQueueLength=meanQL,
        meanWaitingTime=meanWT,
        meanSojournTime=meanST,
        utilization=rho,
        queueLengthDist=ql,
        analyzer="Q-MAM:MAP/M/%d" % c
    )


def qsys_mapmap1(D0_arr: np.ndarray, D1_arr: np.ndarray,
                 D0_srv: np.ndarray, D1_srv: np.ndarray) -> QueueResult:
    """
    Analyze a MAP/MAP/1 queue.

    Both arrival and service processes are Markovian Arrival Processes.

    Args:
        D0_arr: Arrival MAP hidden transition matrix
        D1_arr: Arrival MAP observable transition matrix
        D0_srv: Service MAP hidden transition matrix
        D1_srv: Service MAP observable transition matrix

    Returns:
        QueueResult with queue performance metrics

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mapmap1.m
    """
    # Convert service MAP to equivalent PH representation
    # The service MAP can be viewed as PH with the MAP structure

    D0_srv = np.atleast_2d(np.asarray(D0_srv, dtype=float))
    D1_srv = np.atleast_2d(np.asarray(D1_srv, dtype=float))

    m = D0_srv.shape[0]

    # Service PH initial distribution (stationary distribution of embedded chain)
    D_srv = D0_srv + D1_srv
    pi = np.ones(m) / m
    for _ in range(1000):
        pi_new = pi @ np.linalg.matrix_power(np.eye(m) + D_srv / 100, 100)
        pi_new /= np.sum(pi_new)
        if np.linalg.norm(pi_new - pi) < 1e-12:
            break
        pi = pi_new

    # Use pi as initial distribution and D0_srv as sub-generator
    return qsys_mapph1(D0_arr, D1_arr, pi, D0_srv)


def qsys_mapg1(
    D0: np.ndarray,
    D1: np.ndarray,
    service_moments: np.ndarray,
    num_ql_moms: int = 3,
    num_ql_probs: int = 100,
    num_st_moms: int = 3
) -> QueueResult:
    """
    Analyze a MAP/G/1 queue using BuTools MMAPPH1FCFS.

    The general service time distribution is fitted to a Phase-Type (PH)
    distribution using moment matching before analysis.

    Args:
        D0: MAP hidden transition matrix (n x n)
        D1: MAP arrival transition matrix (n x n)
        service_moments: First k raw moments of service time [E[S], E[S^2], ...]
                        (k = 2 or 3 for best accuracy)
        num_ql_moms: Number of queue length moments to compute (default: 3)
        num_ql_probs: Number of queue length probabilities (default: 100)
        num_st_moms: Number of sojourn time moments to compute (default: 3)

    Returns:
        QueueResult with queue performance metrics

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mapg1.m

    Note:
        Uses the MMAPPH1FCFS solver from BuTools after fitting the general
        service distribution to a PH distribution.
    """
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    service_moments = np.asarray(service_moments, dtype=float).flatten()

    # Fit service distribution to PH using moment matching
    sigma, S = _fit_service_to_ph(service_moments)

    # Build arrival MMAP structure for BuTools (single class)
    D = [D0, D1]

    # Service parameters as cell arrays
    sigma_list = [sigma]
    S_list = [S]

    # Call BuTools solver
    try:
        from line_solver.lib.thirdparty.butools.queues import MMAPPH1FCFS
        nc_moms, nc_distr, st_moms = MMAPPH1FCFS(
            D, sigma_list, S_list,
            'ncMoms', num_ql_moms,
            'ncDistr', num_ql_probs,
            'stMoms', num_st_moms
        )
    except Exception as e:
        # Fallback to simpler approximation
        from ..mam import map_lambda
        lambda_val = map_lambda(D0, D1)
        mean_service = service_moments[0]
        mu = 1.0 / mean_service
        rho = lambda_val * mean_service

        if rho >= 1:
            return QueueResult(
                meanQueueLength=np.inf,
                meanWaitingTime=np.inf,
                meanSojournTime=np.inf,
                utilization=rho,
                analyzer=f"native:fallback_failed:{e}"
            )

        # Simple M/G/1 approximation
        if len(service_moments) >= 2:
            cv2 = (service_moments[1] - service_moments[0]**2) / service_moments[0]**2
        else:
            cv2 = 1.0

        # Pollaczek-Khintchine formula
        mean_ql = rho + (rho**2 * (1 + cv2)) / (2 * (1 - rho))
        mean_st = mean_ql / lambda_val
        mean_wt = max(0, mean_st - mean_service)

        return QueueResult(
            meanQueueLength=mean_ql,
            meanWaitingTime=mean_wt,
            meanSojournTime=mean_st,
            utilization=rho,
            analyzer="native:MG1_approx"
        )

    # Compute utilization from arrival and service rates
    from ..mam import map_lambda
    lambda_val = map_lambda(D0, D1)
    mean_service = service_moments[0]
    rho = lambda_val * mean_service

    # Extract results
    mean_ql = nc_moms[0] if len(nc_moms) > 0 else 0
    mean_st = st_moms[0] if len(st_moms) > 0 else 0
    mean_wt = max(0, mean_st - mean_service)

    return QueueResult(
        meanQueueLength=mean_ql,
        meanWaitingTime=mean_wt,
        meanSojournTime=mean_st,
        utilization=rho,
        queueLengthDist=np.array(nc_distr) if nc_distr is not None else None,
        queueLengthMoments=np.array(nc_moms) if nc_moms is not None else None,
        sojournTimeMoments=np.array(st_moms) if st_moms is not None else None,
        analyzer="BuTools:MMAPPH1FCFS"
    )


def _fit_service_to_ph(moments: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a general service time distribution to a PH distribution.

    Uses moment matching with APH2From3Moments when 3 moments are available,
    otherwise falls back to simpler approximations.

    Args:
        moments: Raw moments [E[S], E[S^2], ...] (at least mean required)

    Returns:
        Tuple of (sigma, S) where sigma is initial distribution and S is generator
    """
    m1 = moments[0]

    if len(moments) >= 3:
        # Try to use BuTools APH from 3 moments
        try:
            from line_solver.lib.thirdparty.butools.ph import APH2From3Moments
            sigma, S = APH2From3Moments(moments[:3])
            return sigma, S
        except:
            pass

    if len(moments) >= 2:
        # Use 2 moments - create PH(2) with matching mean and variance
        m2 = moments[1]
        cv2 = m2 / (m1 * m1) - 1.0  # Squared coefficient of variation

        if cv2 <= 0.001:
            # Near-deterministic: use Erlang with many phases
            k = max(1, min(100, int(round(1.0 / max(cv2, 0.001)))))
            return _create_erlang_ph(m1, k)
        elif cv2 < 1:
            # Hypoexponential: use Erlang-k approximation
            k = max(2, int(round(1.0 / cv2)))
            return _create_erlang_ph(m1, k)
        elif abs(cv2 - 1.0) < 0.001:
            # Exponential
            return _create_exponential_ph(m1)
        else:
            # Hyperexponential: use 2-phase hyperexponential
            return _create_hyperexp2_ph(m1, cv2)
    else:
        # Only mean provided - use exponential
        return _create_exponential_ph(m1)


def _create_exponential_ph(mean: float) -> Tuple[np.ndarray, np.ndarray]:
    """Create an exponential PH distribution."""
    sigma = np.array([1.0])
    S = np.array([[-1.0 / mean]])
    return sigma, S


def _create_erlang_ph(mean: float, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Create an Erlang-k PH distribution."""
    mu = k / mean
    sigma = np.zeros(k)
    sigma[0] = 1.0

    S = np.zeros((k, k))
    for i in range(k):
        S[i, i] = -mu
        if i < k - 1:
            S[i, i + 1] = mu

    return sigma, S


def _create_hyperexp2_ph(mean: float, cv2: float) -> Tuple[np.ndarray, np.ndarray]:
    """Create a 2-phase hyperexponential with matched mean and cv2."""
    # Balanced means approach
    p = 0.5 * (1.0 + np.sqrt((cv2 - 1.0) / (cv2 + 1.0)))
    lambda1 = 2.0 * p / mean
    lambda2 = 2.0 * (1.0 - p) / mean

    sigma = np.array([p, 1.0 - p])
    S = np.diag([-lambda1, -lambda2])

    return sigma, S


__all__ = [
    'QueueResult',
    'ph_to_map',
    'qsys_phph1',
    'qsys_mapph1',
    'qsys_mapm1',
    'qsys_mapmc',
    'qsys_mapmap1',
    'qsys_mapg1',
]
