"""
Analytical Age of Information (AoI) Formulas for Standard Queues.

Native Python implementations of closed-form and semi-analytical formulas
for computing AoI metrics in standard queueing systems.

Key functions:
    FCFS queues: aoi_fcfs_mm1, aoi_fcfs_md1, aoi_fcfs_dm1, aoi_fcfs_mgi1, aoi_fcfs_gim1
    LCFS-PR queues: aoi_lcfspr_mm1, aoi_lcfspr_md1, aoi_lcfspr_dm1, aoi_lcfspr_mgi1, aoi_lcfspr_gim1
    LCFS-S/D queues: aoi_lcfss_mgi1, aoi_lcfss_gim1, aoi_lcfsd_mgi1, aoi_lcfsd_gim1

References:
    Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
    the Stationary Distribution of the Age of Information and Its
    Application to Single-Server Queues," IEEE Trans. Information Theory,
    vol. 65, no. 12, pp. 8305-8324, 2019.

Original MATLAB: matlab/src/api/aoi/aoi_*.m
"""

import numpy as np
from scipy.optimize import brentq
from typing import Tuple, Callable, Optional


def aoi_fcfs_mm1(lambd: float, mu: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for M/M/1 FCFS queue.

    Args:
        lambd: Arrival rate (Poisson arrivals)
        mu: Service rate (exponential service)

    Returns:
        Tuple of (meanAoI, varAoI, peakAoI)

    Raises:
        ValueError: If parameters invalid or system unstable

    References:
        Inoue et al., IEEE Trans. IT, 2019, Section III-A
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")

    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Mean AoI
    meanAoI = (1 / mu) * (1 + 1 / rho + rho ** 2 / (1 - rho))

    # Peak AoI
    peakAoI = (1 / mu) * (1 + 1 / rho + rho / (1 - rho))

    # Variance (second moment from paper)
    E_A = meanAoI
    E_A2 = (2 / mu ** 2) * (1 / rho ** 2 + 1 / rho + rho ** 2 / (1 - rho) ** 2 +
                            1 / (1 - rho) + rho / (1 - rho))
    varAoI = max(0, E_A2 - E_A ** 2)

    return meanAoI, varAoI, peakAoI


def aoi_fcfs_md1(lambd: float, d: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for M/D/1 FCFS queue.

    Args:
        lambd: Arrival rate (Poisson arrivals)
        d: Deterministic service time

    Returns:
        Tuple of (meanAoI, varAoI, peakAoI)

    Raises:
        ValueError: If parameters invalid or system unstable
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if d <= 0:
        raise ValueError("Service time d must be positive")

    rho = lambd * d
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Mean interarrival time
    E_Y = 1 / lambd

    # Mean waiting time (P-K formula with zero variance service)
    E_W = lambd * d ** 2 / (2 * (1 - rho))

    # Mean system time
    E_T = E_W + d

    # Mean AoI
    meanAoI = E_Y + E_T + lambd * d ** 2 / (2 * (1 - rho))

    # Peak AoI
    peakAoI = E_T + E_Y

    # Variance approximation
    varAoI = E_Y ** 2 + (lambd * d ** 2 / (2 * (1 - rho))) ** 2

    return meanAoI, varAoI, peakAoI


def aoi_fcfs_dm1(d: float, mu: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for D/M/1 FCFS queue.

    Args:
        d: Deterministic interarrival time
        mu: Service rate (exponential service)

    Returns:
        Tuple of (meanAoI, varAoI, peakAoI)

    Raises:
        ValueError: If parameters invalid or system unstable
    """
    if d <= 0:
        raise ValueError("Interarrival time d must be positive")
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")

    lambd = 1 / d
    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Find sigma: root of exp(-mu*d*(1-sigma)) = sigma
    def sigma_eq(sig):
        return np.exp(-mu * d * (1 - sig)) - sig

    try:
        sigma = brentq(sigma_eq, 0.001, 0.999)
    except:
        sigma = rho  # Fallback approximation

    # Mean delay
    E_D = 1 / (mu * (1 - sigma))

    # Mean AoI
    meanAoI = d + E_D + sigma / (mu * (1 - sigma))

    # Peak AoI
    peakAoI = d + E_D

    # Variance approximation
    varAoI = 1 / mu ** 2 + (sigma / (mu * (1 - sigma))) ** 2

    return meanAoI, varAoI, peakAoI


def aoi_lcfspr_mm1(lambd: float, mu: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for M/M/1 preemptive LCFS queue.

    In LCFS-PR, a new arrival preempts the current update in service,
    ensuring the freshest update is always served.

    Args:
        lambd: Arrival rate (Poisson arrivals)
        mu: Service rate (exponential service)

    Returns:
        Tuple of (meanAoI, varAoI, peakAoI)

    Raises:
        ValueError: If parameters invalid or system unstable

    References:
        Inoue et al., IEEE Trans. IT, 2019, Section IV
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")

    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Mean AoI for preemptive LCFS
    meanAoI = (1 / mu) * (1 + 1 / rho)

    # Peak AoI
    peakAoI = 1 / mu + 1 / lambd

    # Variance
    E_A = meanAoI
    E_A2 = 2 * (1 / lambd ** 2 + 1 / (lambd * mu) + 1 / mu ** 2)
    varAoI = max(0, E_A2 - E_A ** 2)

    return meanAoI, varAoI, peakAoI


def aoi_lcfspr_md1(lambd: float, d: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for M/D/1 preemptive LCFS queue.

    Args:
        lambd: Arrival rate (Poisson arrivals)
        d: Deterministic service time

    Returns:
        Tuple of (meanAoI, varAoI, peakAoI)
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if d <= 0:
        raise ValueError("Service time d must be positive")

    # For preemptive LCFS with deterministic service
    # Mean AoI = 1/lambda + d
    meanAoI = 1 / lambd + d

    # Peak AoI = 1/lambda + d
    peakAoI = 1 / lambd + d

    # Variance
    varAoI = 1 / lambd ** 2

    return meanAoI, varAoI, peakAoI


def aoi_lcfspr_dm1(d: float, mu: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for D/M/1 preemptive LCFS queue.

    Args:
        d: Deterministic interarrival time
        mu: Service rate (exponential service)

    Returns:
        Tuple of (meanAoI, varAoI, peakAoI)
    """
    if d <= 0:
        raise ValueError("Interarrival time d must be positive")
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")

    # For preemptive LCFS
    meanAoI = d + 1 / mu

    # Peak AoI
    peakAoI = d + 1 / mu

    # Variance
    varAoI = 1 / mu ** 2

    return meanAoI, varAoI, peakAoI


def aoi_fcfs_mgi1(lambd: float, H_lst: Callable, E_H: float,
                  E_H2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for M/GI/1 FCFS queue.

    Args:
        lambd: Arrival rate (Poisson arrivals)
        H_lst: LST of service time distribution, callable(s) -> complex
        E_H: Mean service time (first moment)
        E_H2: Second moment of service time

    Returns:
        Tuple of (meanAoI, lstAoI, peakAoI)

    Raises:
        ValueError: If parameters invalid or system unstable

    References:
        Inoue et al., IEEE Trans. IT, 2019, Theorem 2
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if E_H <= 0:
        raise ValueError("Mean service time E_H must be positive")
    if E_H2 < E_H ** 2:
        raise ValueError("Second moment E_H2 must be >= E_H^2")

    rho = lambd * E_H
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Mean interarrival time
    E_Y = 1 / lambd

    # Mean waiting time (Pollaczek-Khinchine)
    E_W = lambd * E_H2 / (2 * (1 - rho))

    # Mean system time
    E_T = E_W + E_H

    # Mean AoI
    meanAoI = E_Y + E_T + lambd * E_H2 / (2 * (1 - rho))

    # Peak AoI
    peakAoI = E_T + E_Y

    # LST of AoI
    def lstAoI(s):
        s = np.asarray(s, dtype=complex)
        H_s = H_lst(s)
        # Pollaczek-Khinchine LST for waiting time
        W_s = (1 - rho) * s / (s - lambd + lambd * H_s)
        # AoI LST
        return (lambd * H_s) / (s + lambd - lambd * H_s) * W_s

    return meanAoI, lstAoI, peakAoI


def aoi_fcfs_gim1(Y_lst: Callable, mu: float, E_Y: float,
                  E_Y2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for GI/M/1 FCFS queue.

    Args:
        Y_lst: LST of interarrival time distribution, callable(s) -> complex
        mu: Service rate (exponential service)
        E_Y: Mean interarrival time (first moment)
        E_Y2: Second moment of interarrival time

    Returns:
        Tuple of (meanAoI, lstAoI, peakAoI)

    Raises:
        ValueError: If parameters invalid or system unstable

    References:
        Inoue et al., IEEE Trans. IT, 2019, Theorem 3
    """
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")
    if E_Y <= 0:
        raise ValueError("Mean interarrival time E_Y must be positive")
    if E_Y2 < E_Y ** 2:
        raise ValueError("Second moment E_Y2 must be >= E_Y^2")

    lambd = 1 / E_Y
    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Find sigma: root of Y*(mu - mu*sigma) = sigma
    def sigma_eq(sig):
        return float(np.real(Y_lst(mu - mu * sig))) - sig

    try:
        sigma = brentq(sigma_eq, 0.001, 0.999)
    except:
        sigma = 0.5
        for _ in range(100):
            sigma_new = float(np.real(Y_lst(mu - mu * sigma)))
            if abs(sigma_new - sigma) < 1e-10:
                break
            sigma = sigma_new

    # Mean delay
    E_D = 1 / (mu * (1 - sigma))

    # Mean AoI
    meanAoI = E_Y + E_D + sigma / (mu * (1 - sigma))

    # Peak AoI
    peakAoI = E_Y + E_D

    # LST of AoI
    def lstAoI(s):
        s_arr = np.atleast_1d(np.asarray(s, dtype=complex))
        result = np.zeros_like(s_arr)

        for i, si in enumerate(s_arr):
            # Find sigma(s)
            def sig_eq(sig):
                return float(np.real(Y_lst(si + mu - mu * sig))) - sig

            try:
                sigma_s = brentq(sig_eq, 0.001, 0.999)
            except:
                sigma_s = sigma

            D_s = (1 - sigma) * mu / (si + mu - mu * sigma_s)
            result[i] = (mu * sigma_s) / (si + mu - mu * sigma_s) * D_s

        if np.isscalar(s):
            return result[0]
        return result

    return meanAoI, lstAoI, peakAoI


def aoi_lcfspr_mgi1(lambd: float, H_lst: Callable, E_H: float,
                    E_H2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for M/GI/1 preemptive LCFS queue.

    Args:
        lambd: Arrival rate (Poisson arrivals)
        H_lst: LST of service time distribution
        E_H: Mean service time
        E_H2: Second moment of service time

    Returns:
        Tuple of (meanAoI, lstAoI, peakAoI)
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if E_H <= 0:
        raise ValueError("Mean service time E_H must be positive")

    # For preemptive LCFS with M/GI/1
    meanAoI = 1 / lambd + E_H

    # Peak AoI
    peakAoI = 1 / lambd + E_H

    # LST (simplified)
    def lstAoI(s):
        s = np.asarray(s, dtype=complex)
        return H_lst(s) * lambd / (lambd + s)

    return meanAoI, lstAoI, peakAoI


def aoi_lcfspr_gim1(Y_lst: Callable, mu: float, E_Y: float,
                    E_Y2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for GI/M/1 preemptive LCFS queue.

    Args:
        Y_lst: LST of interarrival time distribution
        mu: Service rate (exponential service)
        E_Y: Mean interarrival time
        E_Y2: Second moment of interarrival time

    Returns:
        Tuple of (meanAoI, lstAoI, peakAoI)
    """
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")
    if E_Y <= 0:
        raise ValueError("Mean interarrival time E_Y must be positive")

    # For preemptive LCFS
    meanAoI = E_Y + 1 / mu

    # Peak AoI
    peakAoI = E_Y + 1 / mu

    # LST (simplified)
    def lstAoI(s):
        s = np.asarray(s, dtype=complex)
        return Y_lst(s) * mu / (mu + s)

    return meanAoI, lstAoI, peakAoI


def aoi_lcfss_mgi1(lambd: float, H_lst: Callable, E_H: float,
                   E_H2: float) -> Tuple[float, float]:
    """
    Mean AoI for M/GI/1 LCFS with service discarding.

    In LCFS-S, if a new update arrives while another is being served,
    the arriving update is discarded (not preemptive, but skips queue).

    Args:
        lambd: Arrival rate
        H_lst: LST of service time distribution
        E_H: Mean service time
        E_H2: Second moment of service time

    Returns:
        Tuple of (meanAoI, peakAoI)
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if E_H <= 0:
        raise ValueError("Mean service time E_H must be positive")

    rho = lambd * E_H
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Approximate mean AoI for LCFS-S
    meanAoI = 1 / lambd + E_H / (1 - rho)
    peakAoI = 1 / lambd + E_H

    return meanAoI, peakAoI


def aoi_lcfss_gim1(Y_lst: Callable, mu: float, E_Y: float,
                   E_Y2: float) -> Tuple[float, float]:
    """
    Mean AoI for GI/M/1 LCFS with service discarding.

    Args:
        Y_lst: LST of interarrival time distribution
        mu: Service rate
        E_Y: Mean interarrival time
        E_Y2: Second moment of interarrival time

    Returns:
        Tuple of (meanAoI, peakAoI)
    """
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")
    if E_Y <= 0:
        raise ValueError("Mean interarrival time E_Y must be positive")

    lambd = 1 / E_Y
    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    meanAoI = E_Y + 1 / (mu * (1 - rho))
    peakAoI = E_Y + 1 / mu

    return meanAoI, peakAoI


def aoi_lcfsd_mgi1(lambd: float, H_lst: Callable, E_H: float,
                   E_H2: float) -> Tuple[float, float]:
    """
    Mean AoI for M/GI/1 LCFS with departure discarding.

    In LCFS-D, if multiple updates accumulate, only the most recent
    is kept when service completes.

    Args:
        lambd: Arrival rate
        H_lst: LST of service time distribution
        E_H: Mean service time
        E_H2: Second moment of service time

    Returns:
        Tuple of (meanAoI, peakAoI)
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if E_H <= 0:
        raise ValueError("Mean service time E_H must be positive")

    rho = lambd * E_H
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # For LCFS-D, updates are effectively sampled at departures
    meanAoI = 1 / lambd + E_H + lambd * E_H2 / (2 * (1 - rho))
    peakAoI = 1 / lambd + E_H

    return meanAoI, peakAoI


def aoi_lcfsd_gim1(Y_lst: Callable, mu: float, E_Y: float,
                   E_Y2: float) -> Tuple[float, float]:
    """
    Mean AoI for GI/M/1 LCFS with departure discarding.

    Args:
        Y_lst: LST of interarrival time distribution
        mu: Service rate
        E_Y: Mean interarrival time
        E_Y2: Second moment of interarrival time

    Returns:
        Tuple of (meanAoI, peakAoI)
    """
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")
    if E_Y <= 0:
        raise ValueError("Mean interarrival time E_Y must be positive")

    lambd = 1 / E_Y
    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    # Find sigma
    def sigma_eq(sig):
        return float(np.real(Y_lst(mu - mu * sig))) - sig

    try:
        sigma = brentq(sigma_eq, 0.001, 0.999)
    except:
        sigma = rho

    E_D = 1 / (mu * (1 - sigma))

    meanAoI = E_Y + E_D
    peakAoI = E_Y + E_D

    return meanAoI, peakAoI


__all__ = [
    # FCFS queues
    'aoi_fcfs_mm1',
    'aoi_fcfs_md1',
    'aoi_fcfs_dm1',
    'aoi_fcfs_mgi1',
    'aoi_fcfs_gim1',
    # LCFS preemptive queues
    'aoi_lcfspr_mm1',
    'aoi_lcfspr_md1',
    'aoi_lcfspr_dm1',
    'aoi_lcfspr_mgi1',
    'aoi_lcfspr_gim1',
    # LCFS with discarding
    'aoi_lcfss_mgi1',
    'aoi_lcfss_gim1',
    'aoi_lcfsd_mgi1',
    'aoi_lcfsd_gim1',
]
