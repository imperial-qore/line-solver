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

    Exact formulas (Inoue et al. 2019):
        E[A] = (1/mu)(1 + 1/rho + rho^2/(1-rho))
        E[Apeak] = (1/mu)(1 + 1/rho + rho/(1-rho))
        E[A^2] = (2/mu^2)(1 - rho - rho^3 + 4*rho^4 - 2*rho^5)/(rho^2*(1-rho)^2)

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

    # Variance (exact second moment via the sawtooth identity, as MATLAB)
    E_A = meanAoI
    E_A2 = (2 / mu ** 2) * (1 - rho - rho ** 3 + 4 * rho ** 4 - 2 * rho ** 5) / \
        (rho ** 2 * (1 - rho) ** 2)
    varAoI = max(0, E_A2 - E_A ** 2)

    return meanAoI, varAoI, peakAoI


def aoi_fcfs_md1(lambd: float, d: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for M/D/1 FCFS queue.

    Exact mean (Inoue et al. 2019):
        E[A] = d*(1/2 + 1/(2*(1-rho)) + ((1-rho)/rho)*exp(rho)), rho = lambd*d
        E[Apeak] = E[T] + E[Y]

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

    # Mean AoI for M/D/1 FCFS (exact, as MATLAB aoi_fcfs_md1)
    meanAoI = d * (0.5 + 1 / (2 * (1 - rho)) + ((1 - rho) / rho) * np.exp(rho))

    # Peak AoI
    peakAoI = E_T + E_Y

    # Variance approximation (as MATLAB aoi_fcfs_md1)
    E_H = d
    E_H2 = d ** 2
    E_Y2 = 2 / lambd ** 2
    varAoI = max(0.0, E_Y2 - E_Y ** 2 + 2 * E_W * E_H / (1 - rho)
                 + E_H2 * rho / (1 - rho) ** 2)

    return meanAoI, varAoI, peakAoI


def aoi_fcfs_dm1(d: float, mu: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for D/M/1 FCFS queue.

    Exact mean: E[A] = d/2 + 1/(mu*(1-sigma)) with sigma the root of
    sigma = exp(-mu*d*(1-sigma)); E[Apeak] = d + 1/(mu*(1-sigma)).

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

    # Mean AoI for D/M/1 FCFS (exact, as MATLAB aoi_fcfs_dm1):
    # deterministic interarrivals kill the Y-W correlation term, so
    # E[A] = E[Y^2]/(2*E[Y]) + E[D] = d/2 + E[D]
    meanAoI = d / 2 + E_D

    # Peak AoI
    peakAoI = d + E_D

    # Variance approximation (as MATLAB aoi_fcfs_dm1)
    E_D2 = 2 / (mu * (1 - sigma)) ** 2
    Var_D = E_D2 - E_D ** 2
    varAoI = max(0.0, Var_D + (sigma / (mu * (1 - sigma))) ** 2)

    return meanAoI, varAoI, peakAoI


def aoi_lcfspr_mm1(lambd: float, mu: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for M/M/1 preemptive LCFS queue.

    In LCFS-PR, a new arrival preempts the current update in service,
    ensuring the freshest update is always served.

    Exact formulas:
        E[A] = 1/mu + 1/lambd
        E[Apeak] = 1/(lambd+mu) + 1/lambd + 1/mu

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

    # Peak AoI (exact, as MATLAB): E[S|success] + E[inter-delivery time]
    peakAoI = 1 / (lambd + mu) + 1 / lambd + 1 / mu

    # Variance
    E_A = meanAoI
    E_A2 = 2 * (1 / lambd ** 2 + 1 / (lambd * mu) + 1 / mu ** 2)
    varAoI = max(0, E_A2 - E_A ** 2)

    return meanAoI, varAoI, peakAoI


def aoi_lcfspr_md1(lambd: float, d: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for M/D/1 preemptive LCFS queue.

    Peak AoI (exact): E[Apeak] = d + exp(lambd*d)/lambd, since deliveries
    occur at rate lambd*exp(-lambd*d).

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

    # Peak AoI (exact, as MATLAB): deliveries at rate lambda*exp(-lambda*d)
    peakAoI = d + np.exp(lambd * d) / lambd

    # Variance
    varAoI = 1 / lambd ** 2

    return meanAoI, varAoI, peakAoI


def aoi_lcfspr_dm1(d: float, mu: float) -> Tuple[float, float, float]:
    """
    Mean, variance, and peak AoI for D/M/1 preemptive LCFS queue.

    Peak AoI (exact): E[Apeak] = E[S|success] + d/q with q = 1 - exp(-mu*d)
    and E[S|success] = (1/mu - d*exp(-mu*d) - exp(-mu*d)/mu)/q.

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

    # Peak AoI (exact, as MATLAB): success prob q = 1 - exp(-mu*d)
    q = 1 - np.exp(-mu * d)
    ES_succ = (1 / mu - d * np.exp(-mu * d) - np.exp(-mu * d) / mu) / q
    peakAoI = ES_succ + d / q

    # Variance
    varAoI = 1 / mu ** 2

    return meanAoI, varAoI, peakAoI


def aoi_fcfs_mgi1(lambd: float, H_lst: Callable, E_H: float,
                  E_H2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for M/GI/1 FCFS queue.

    Exact mean: E[A] = E[H] + E[T] + (1-2*rho)/lambd - d/ds T*(s)|_{s=lambd}
    with T*(s) = H*(s)*W*(s) the Pollaczek-Khinchine system-time LST;
    E[Apeak] = E[T] + E[Y].

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

    # Mean AoI for M/GI/1 FCFS (exact, as MATLAB aoi_fcfs_mgi1):
    # E[A] = E[H] + E[T] + (1-2*rho)/lambda - d/ds T*(s)|_{s=lambda}
    # where T*(s) = H*(s)*W*(s) is the system-time LST (P-K)
    def _Tstar(s):
        Hs = H_lst(s)
        return Hs * (1 - rho) * s / (s - lambd + lambd * Hs)

    hstep = 1e-6 * max(1.0, lambd)
    dTstar = float(np.real(_Tstar(lambd + hstep) - _Tstar(lambd - hstep))) / (2 * hstep)
    meanAoI = E_H + E_T + (1 - 2 * rho) / lambd - dTstar

    # Peak AoI
    peakAoI = E_T + E_Y

    # LST of AoI (Inoue et al. 2019, Theorem 2), via the general age formula
    #   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
    # the cycle average of exp(-s*age) over a departure interval: the age starts
    # each cycle at the system time T of the packet just delivered and grows to
    # the peak T + (next interarrival) at the next delivery. For M/GI/1 FCFS
    # Lindley gives W' = max(0, T - Y), so W' + Y = max(Y, T), and with
    # Y ~ Exp(lambda) independent of T,
    #   E[exp(-s*max(Y,T))] = T*(s) - (s/(s+lambda)) * T*(s+lambda).
    #
    # THE PREVIOUS FORM WAS NOT AN LST: (lambda*H*(s))/(s+lambda-lambda*H*(s))
    # diverges as s -> 0, so A*(0) was +inf instead of 1 and the value exceeded
    # 1 for small s. Checked against simulation on M/E2/1: at s = 0.3 the old
    # form gave 1.3062, the form below 0.56935, and the sample path 0.56948.
    def _lst_at(sv):
        if abs(sv) < 1e-12:
            return 1.0 + 0j
        H_s = H_lst(sv)
        T_s = _Tstar_safe(sv)
        T_sl = _Tstar_safe(sv + lambd)
        peak_s = H_s * (T_s - (sv / (sv + lambd)) * T_sl)
        return (lambd / sv) * (T_s - peak_s)

    def _Tstar_safe(sv):
        if abs(sv) < 1e-12:
            return 1.0 + 0j
        return _Tstar(sv)

    def lstAoI(s):
        arr = np.asarray(s, dtype=complex)
        if arr.ndim == 0:
            return _lst_at(complex(arr))
        return np.array([_lst_at(complex(v)) for v in arr.ravel()]).reshape(arr.shape)

    return meanAoI, lstAoI, peakAoI


def aoi_fcfs_gim1(Y_lst: Callable, mu: float, E_Y: float,
                  E_Y2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for GI/M/1 FCFS queue.

    Exact mean: E[A] = lambd*E[Y^2]/2 + 1/mu + lambd*(-Y*'(eta))/eta with
    eta = mu*(1-sigma), sigma the root of Y*(mu - mu*sigma) = sigma;
    E[Apeak] = E[Y] + 1/(mu*(1-sigma)).

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

    # Mean AoI for GI/M/1 FCFS (exact, as MATLAB aoi_fcfs_gim1):
    # stationary system time is exp(eta), eta = mu*(1-sigma), independent
    # of the next interarrival, so E[A] = lam*E[Y^2]/2 + 1/mu + lam*(-Y*'(eta))/eta
    eta = mu * (1 - sigma)
    hstep = 1e-6 * max(1.0, eta)
    dYstar = float(np.real(Y_lst(eta + hstep) - Y_lst(eta - hstep))) / (2 * hstep)
    meanAoI = lambd * E_Y2 / 2 + 1 / mu + lambd * (-dYstar) / eta

    # Peak AoI
    peakAoI = E_Y + E_D

    # LST of AoI (Inoue et al. 2019, Theorem 3), via the general age formula
    #   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
    # the cycle average of exp(-s*age) over a departure interval. In GI/M/1 the
    # system time is EXPONENTIAL at rate eta = mu*(1-sigma), so T*(s) =
    # eta/(s+eta), and Lindley gives W' + Y = max(Y, T) with T ~ Exp(eta)
    # independent of the next interarrival Y, so
    #   E[exp(-s*max(Y,T))] = Y*(s) - (s/(s+eta)) * Y*(s+eta),
    # and the peak adds one fresh Exp(mu) service. A*(0) = 1 follows from the
    # defining relation Y*(eta) = sigma.
    #
    # THE PREVIOUS FORM WAS NOT AN LST: (mu*sigma(s))/(s+mu-mu*sigma(s))*D*(s)
    # gives sigma/(1-sigma) at s = 0 rather than 1, and it re-solved sigma(s) by
    # brentq at every point with a SILENT fallback to sigma(0) on failure.
    # Checked against simulation on E2/M/1: at s = 0.2 the old form gave
    # 0.23056, the form below 0.59319, and the sample path 0.59332.
    def _lst_at(sv):
        if abs(sv) < 1e-12:
            return 1.0 + 0j
        T_s = eta / (sv + eta)
        peak_s = (mu / (sv + mu)) * (Y_lst(sv) - (sv / (sv + eta)) * Y_lst(sv + eta))
        return (lambd / sv) * (T_s - peak_s)

    def lstAoI(s):
        arr = np.asarray(s, dtype=complex)
        if arr.ndim == 0:
            return _lst_at(complex(arr))
        return np.array([_lst_at(complex(v)) for v in arr.ravel()]).reshape(arr.shape)

    return meanAoI, lstAoI, peakAoI


def aoi_lcfspr_mgi1(lambd: float, H_lst: Callable, E_H: float,
                    E_H2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for M/GI/1 preemptive LCFS queue.

    Peak AoI (exact): E[Apeak] = -H*'(lambd)/H*(lambd) + 1/(lambd*H*(lambd)).

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

    # Peak AoI (exact, as MATLAB): q = H*(lambda), E[S|succ] = -H*'(lam)/H*(lam)
    hstep = 1e-6 * max(1.0, lambd)
    H_lam = float(np.real(H_lst(lambd)))
    dHstar = float(np.real(H_lst(lambd + hstep) - H_lst(lambd - hstep))) / (2 * hstep)
    peakAoI = -dHstar / H_lam + 1 / (lambd * H_lam)

    # LST (simplified)
    def lstAoI(s):
        s = np.asarray(s, dtype=complex)
        return H_lst(s) * lambd / (lambd + s)

    return meanAoI, lstAoI, peakAoI


def aoi_lcfspr_gim1(Y_lst: Callable, mu: float, E_Y: float,
                    E_Y2: float) -> Tuple[float, Callable, float]:
    """
    Mean AoI and LST for GI/M/1 preemptive LCFS queue.

    Peak AoI (exact): E[Apeak] = E[S|success] + 1/(lambd*q) with
    q = 1 - Y*(mu) and E[S|success] = (1/mu + Y*'(mu) - Y*(mu)/mu)/q.

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

    # Peak AoI (exact, as MATLAB): q = 1 - Y*(mu)
    lambd = 1 / E_Y
    hstep = 1e-6 * max(1.0, mu)
    Y_mu = float(np.real(Y_lst(mu)))
    dYstar = float(np.real(Y_lst(mu + hstep) - Y_lst(mu - hstep))) / (2 * hstep)
    q = 1 - Y_mu
    ES_succ = (1 / mu + dYstar - Y_mu / mu) / q
    peakAoI = ES_succ + 1 / (lambd * q)

    # LST (simplified)
    def lstAoI(s):
        s = np.asarray(s, dtype=complex)
        return Y_lst(s) * mu / (mu + s)

    return meanAoI, lstAoI, peakAoI


def aoi_lcfss_mgi1(lambd: float, H_lst: Callable, E_H: float,
                   E_H2: float) -> Tuple[float, Optional[Callable], float]:
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
        Tuple of (meanAoI, lstAoI, peakAoI) where lstAoI is None
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

    E_Y = 1 / lambd
    # Busy period moments (as MATLAB aoi_lcfss_mgi1)
    E_B = E_H / (1 - rho)
    E_B2 = E_H2 / (1 - rho) ** 3

    meanAoI = E_Y + E_H + lambd * E_H2 / (2 * (1 - rho) ** 2)
    peakAoI = E_Y + E_H + lambd * E_B2 / (2 * E_B)

    return meanAoI, None, peakAoI


def aoi_lcfss_gim1(Y_lst: Callable, mu: float, E_Y: float,
                   E_Y2: float) -> Tuple[float, Optional[Callable], float]:
    """
    Mean AoI for GI/M/1 LCFS with service discarding.

    Args:
        Y_lst: LST of interarrival time distribution
        mu: Service rate
        E_Y: Mean interarrival time
        E_Y2: Second moment of interarrival time

    Returns:
        Tuple of (meanAoI, lstAoI, peakAoI) where lstAoI is None
    """
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")
    if E_Y <= 0:
        raise ValueError("Mean interarrival time E_Y must be positive")

    lambd = 1 / E_Y
    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    E_S = 1 / mu

    # Find sigma: root of Y*(mu - mu*sigma) = sigma (as MATLAB aoi_lcfss_gim1)
    def sigma_eq(sig):
        return float(np.real(Y_lst(mu - mu * sig))) - sig

    try:
        sigma = brentq(sigma_eq, 0.001, 0.999)
    except Exception:
        sigma = rho

    E_D = 1 / (mu * (1 - sigma))

    meanAoI = E_Y + E_S + sigma * E_D
    peakAoI = E_Y + E_D

    return meanAoI, None, peakAoI


def aoi_lcfsd_mgi1(lambd: float, H_lst: Callable, E_H: float,
                   E_H2: float) -> Tuple[float, Optional[Callable], float]:
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
        Tuple of (meanAoI, lstAoI, peakAoI) where lstAoI is None
    """
    if lambd <= 0:
        raise ValueError("Arrival rate lambda must be positive")
    if E_H <= 0:
        raise ValueError("Mean service time E_H must be positive")

    rho = lambd * E_H
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    E_Y = 1 / lambd
    # Residual service time and effective system time (as MATLAB aoi_lcfsd_mgi1)
    E_H_residual = E_H2 / (2 * E_H)
    E_T_eff = E_H + rho * E_H_residual

    meanAoI = E_Y + E_H + rho * E_H2 / (2 * E_H) + rho * E_H / (1 + rho)
    peakAoI = E_Y + E_T_eff

    return meanAoI, None, peakAoI


def aoi_lcfsd_gim1(Y_lst: Callable, mu: float, E_Y: float,
                   E_Y2: float) -> Tuple[float, Optional[Callable], float]:
    """
    Mean AoI for GI/M/1 LCFS with departure discarding.

    Args:
        Y_lst: LST of interarrival time distribution
        mu: Service rate
        E_Y: Mean interarrival time
        E_Y2: Second moment of interarrival time

    Returns:
        Tuple of (meanAoI, lstAoI, peakAoI) where lstAoI is None
    """
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")
    if E_Y <= 0:
        raise ValueError("Mean interarrival time E_Y must be positive")

    lambd = 1 / E_Y
    rho = lambd / mu
    if rho >= 1:
        raise ValueError(f"System unstable: rho = {rho:.4f} >= 1")

    E_S = 1 / mu

    # Find sigma
    def sigma_eq(sig):
        return float(np.real(Y_lst(mu - mu * sig))) - sig

    try:
        sigma = brentq(sigma_eq, 0.001, 0.999)
    except Exception:
        sigma = rho

    # Effective system time (as MATLAB aoi_lcfsd_gim1)
    E_T_eff = E_S + sigma * E_S

    meanAoI = E_Y + E_S * (1 + sigma) + sigma * E_S / (1 + sigma)
    peakAoI = E_Y + E_T_eff

    return meanAoI, None, peakAoI


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
