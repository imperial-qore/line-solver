"""
Markov Modulated Poisson Process (MMPP) functions for KPC-Toolbox.

Native Python implementations of MMPP fitting and analysis.
Matches MATLAB: matlab/lib/kpctoolbox/mmpp/
"""

import numpy as np
from typing import Tuple, Optional
from scipy.optimize import fsolve


def mmpp2_fit3(E1: float, E2: float, E3: float, G2: float
               ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Core MMPP2 fitter matching MATLAB mmpp2_fit3(E1, E2, E3, G2).

    Fits a 2-state MMPP to match first three moments and G2 parameter.
    G2 specifies the ratio of consecutive autocorrelations: rho(i)/rho(i-1).

    Args:
        E1: First moment (mean)
        E2: Second moment
        E3: Third moment
        G2: Autocorrelation decay ratio

    Returns:
        Tuple of (D0, D1) matrices
    """
    SCV = (E2 - E1**2) / E1**2

    if (G2 < 1e-6) or (G2 == 0):
        # If G2 is very close to zero, fit with MAP(1)
        mu00 = 2 * (6 * E1**3 * SCV - E3) / E1 / (6 * E1**3 * SCV + 3 * E1**3 * SCV**2 + 3 * E1**3 - 2 * E3)
        mu11 = 0.0
        q01 = 9 * E1**5 * (SCV - 1) * (SCV**2 - 2 * SCV + 1) / (6 * E1**3 * SCV - E3) / (6 * E1**3 * SCV + 3 * E1**3 * SCV**2 + 3 * E1**3 - 2 * E3)
        q10 = -3 * (SCV - 1) * E1**2 / (6 * E1**3 * SCV - E3)
    else:
        # Full MMPP(2) fitting - use exact MATLAB symbolic expressions
        # Compute the discriminant that appears throughout
        disc_val = (E3**2 - 12 * E1**3 * SCV * E3 + 6 * E1**3 * G2 * E3 -
                   6 * G2 * SCV * E1**3 * E3 + 18 * G2 * SCV**3 * E1**6 -
                   18 * E1**6 * G2 * SCV**2 + 9 * E1**6 * G2**2 +
                   36 * E1**6 * SCV**2 + 18 * E1**6 * G2 * SCV -
                   18 * E1**6 * SCV * G2**2 + 9 * E1**6 * SCV**2 * G2**2 -
                   18 * E1**6 * G2)

        disc = np.sqrt(disc_val) if disc_val >= 0 else np.sqrt(complex(disc_val))

        denom = -3 * E1**3 * SCV**2 - 6 * E1**3 * SCV - 3 * E1**3 + 2 * E3
        inner = -3 * E1**3 * G2 + 3 * E1**3 * G2 * SCV - 6 * E1**3 * SCV + E3 + disc

        # mu11
        mu11 = inner / E1 / denom

        # mu00 - full MATLAB symbolic expression
        t = inner / denom  # = term used in mu00 formula
        mu00 = G2 * (-4 * E3 * G2 + 4 * t * E3 * G2 - 18 * E1**3 * t * G2 - 18 * E1**3 * t * G2 * SCV**2 - 12 * E1**3 * G2**2 - 12 * E1**3 * t * G2**2 * SCV + 12 * E1**3 * t * G2 * SCV + 12 * E1**3 * G2 * SCV**2 - 9 * E1**3 * t * SCV + 3 * E1**3 * t + 12 * E1**3 * G2**2 * SCV + 9 * E1**3 * t * SCV**2 + 12 * E1**3 * G2 + 12 * E1**3 * t * G2**2 - 3 * E1**3 * t * SCV**3 - 4 * E3 * G2**2) / (12 * E1**3 * G2**3 * SCV + 3 * E1**3 * SCV**3 * G2 - 12 * E1**3 * G2**3 + 18 * E1**3 * G2**2 * SCV**2 - 3 * E1**3 * G2 + 27 * E1**3 * t * G2 * SCV**2 - 9 * E1**3 * G2 * SCV**2 + 18 * E1**3 * G2**2 - 12 * E1**3 * G2**2 * SCV + 9 * E1**3 * G2 * SCV - 12 * E1**3 * t * G2**3 * SCV - 9 * E1**3 * t * SCV**3 * G2 - 24 * E1**3 * t * G2**2 * SCV**2 - t * E3 * SCV**2 + 4 * t * E3 * G2**2 + 12 * E1**3 * t * G2**3 - t * E3 + 2 * t * E3 * SCV + 9 * E1**3 * t * G2 + 24 * E1**3 * t * G2**2 * SCV - 27 * E1**3 * t * G2 * SCV + 6 * E1**3 * t * SCV - 12 * E1**3 * t * SCV**2 - 24 * E1**3 * t * G2**2 + 6 * E1**3 * t * SCV**3 - 4 * E3 * G2**2) / E1

        # q01 - full MATLAB symbolic expression
        q01 = -3 * E1**2 * (-6 * t * E1**2 * SCV + 12 * t * E1**2 * G2 * SCV - 6 * G2 * SCV * E1**2 - 3 * t * E1**2 * G2 + t / E1 * E3 + 3 * E1**2 * G2 + 6 * t * E1**2 * SCV**2 - 9 * t * E1**2 * SCV**2 * G2 + 3 * E1**2 * G2 * SCV**2 - t / E1 * E3 * SCV - 6 * t * E1**2 * G2**2 * SCV + 6 * E1**2 * G2**2 * SCV + 3 * t * E1**2 * G2**2 - G2 * t / E1 * E3 - 3 * E1**2 * G2**2 + 3 * t * E1**2 * SCV**2 * G2**2 - 3 * E1**2 * SCV**2 * G2**2 + G2 * SCV * t / E1 * E3) / (-45 * t * E1**5 * G2 * SCV**2 + 18 * G2**2 * E1**5 * SCV + 18 * E1**5 * G2**3 - 27 * E1**5 * G2**2 * SCV**2 + 6 * E1**2 * G2**2 * E3 - 27 * E1**5 * G2**2 - 18 * E1**5 * G2**3 * SCV - 18 * E1**5 * G2 * SCV + 18 * E1**5 * G2 * SCV**2 + 3 * E1**2 * G2 * E3 - 3 * E1**2 * G2 * E3 * SCV + t / E1 * E3**2 + 3 * t * E1**2 * G2 * SCV * E3 - 36 * t * E1**5 * G2**2 * SCV + 36 * t * E1**5 * G2**2 + 36 * t * E1**5 * SCV**2 + 45 * t * E1**5 * G2 * SCV - 12 * t * E1**2 * SCV * E3 - 3 * t * E1**2 * G2 * E3 + 9 * t * E1**5 * G2 * SCV**3 + 36 * t * E1**5 * G2**2 * SCV**2 - 6 * t * E1**2 * G2**2 * E3 + 18 * t * E1**5 * G2**3 * SCV - 18 * t * E1**5 * G2**3 - 9 * t * E1**5 * G2)

        # q10 - full MATLAB symbolic expression
        q10 = 3 * (-3 * E1**3 * t * SCV**3 - 3 * E1**3 * t * G2 * SCV**2 + 6 * E1**3 * SCV**2 + 3 * E1**3 * G2 * SCV**2 + 3 * E1**3 * t * SCV**2 + 6 * E1**3 * t * G2 * SCV - E3 * SCV - 6 * E1**3 * SCV + t * E3 * SCV - 6 * E1**3 * G2 * SCV - 3 * E1**3 * t * SCV - t * E3 + 3 * E1**3 * G2 - 3 * E1**3 * t * G2 + 3 * E1**3 * t + E3) * E1**2 * (-1 + G2) / (E3**2 - 12 * E1**3 * SCV * E3 + 6 * E1**3 * G2 * E3 - 6 * G2 * SCV * E1**3 * E3 + 18 * G2 * SCV**3 * E1**6 - 18 * E1**6 * G2 * SCV**2 + 9 * E1**6 * G2**2 + 36 * E1**6 * SCV**2 + 18 * E1**6 * G2 * SCV - 18 * E1**6 * SCV * G2**2 + 9 * E1**6 * SCV**2 * G2**2 - 18 * E1**6 * G2)

    D0 = np.array([[-mu00 - q01, q01],
                   [q10, -mu11 - q10]])
    D1 = np.array([[mu00, 0.0],
                   [0.0, mu11]])

    # Take real part if complex
    D0 = np.real(D0)
    D1 = np.real(D1)

    return (D0, D1)


def mmpp2_fit(E1: float, E2: float, E3: float, ACFLAG1: float
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit MMPP2 from moments and lag-1 autocorrelation.

    Convenience wrapper: converts ACFLAG1 to G2, then calls mmpp2_fit3.

    Args:
        E1: First moment (mean)
        E2: Second moment
        E3: Third moment
        ACFLAG1: Lag-1 autocorrelation

    Returns:
        Tuple of (D0, D1) matrices
    """
    SCV = (E2 - E1**2) / (E1**2)
    if abs(SCV - 1) < 1e-10 or SCV == 0:
        G2 = 0.0
    else:
        G2 = ACFLAG1 / ((1 - 1 / SCV) / 2)
    return mmpp2_fit3(E1, E2, E3, G2)


def mmpp2_fit1(mean: float, scv: float, skew: float, idc: float
               ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit MMPP2 from mean, SCV, skewness, and index of dispersion.

    Matches MATLAB: mmpp2_fit1(mean, scv, skew, idc)

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness (-1 for automatic)
        idc: Index of dispersion for counts

    Returns:
        Tuple of (D0, D1) matrices
    """
    from ..mam import map2_fit
    E1 = mean
    E2 = (1 + scv) * E1**2
    g2 = -(scv - idc) / (-1 + idc)
    if skew == -1:
        E3 = -1.0
    else:
        E3 = -(2 * E1**3 - 3 * E1 * E2 - skew * (E2 - E1**2)**1.5)
    MAP, err = map2_fit(E1, E2, E3, g2)
    return MAP


def mmpp2_fit2(mean: float, scv: float, skew: float, g2: float
               ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit MMPP2 from mean, SCV, skewness, and G2 parameter.

    Matches MATLAB: mmpp2_fit2(mean, scv, skew, g2)

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness
        g2: Autocorrelation decay ratio

    Returns:
        Tuple of (D0, D1) matrices
    """
    if scv == 1:
        from ..mam import map_exponential
        return map_exponential(mean)
    E1 = mean
    E2 = (1 + scv) * E1**2
    E3 = -(2 * E1**3 - 3 * E1 * E2 - skew * (E2 - E1**2)**1.5)
    return mmpp2_fit3(E1, E2, E3, g2)


def mmpp2_fit4(mean: float, scv: float, skew: float, acf1: float
               ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit MMPP2 from mean, SCV, skewness, and lag-1 autocorrelation.

    Matches MATLAB: mmpp2_fit4(mean, scv, skew, acf1)

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness (-1 for automatic)
        acf1: Lag-1 autocorrelation

    Returns:
        Tuple of (D0, D1) matrices
    """
    E1 = mean
    E2 = (1 + scv) * E1**2
    if skew == -1:
        E3 = -1.0
    else:
        E3 = -(2 * E1**3 - 3 * E1 * E2 - skew * (E2 - E1**2)**1.5)
    rho0 = (1 - 1 / scv) / 2
    g2 = acf1 / rho0
    return mmpp2_fit3(E1, E2, E3, g2)


def mmpp2_fitc(mu: float, bt1: float, bt2: float, binf: float,
               m3t2: float, t1: float, t2: float
               ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit MMPP2 from counting process statistics (Heffes-Lucantoni method).

    Matches MATLAB: mmpp2_fitc(mu, bt1, bt2, binf, m3t2, t1, t2)

    Args:
        mu: Arrival rate
        bt1: IDC at scale t1
        bt2: IDC at scale t2
        binf: IDC for t->inf
        m3t2: Third central moment at scale t2
        t1: First time scale
        t2: Second time scale

    Returns:
        Tuple of (D0, D1) matrices
    """
    # Degenerate case
    if abs(binf - 1) < 1e-8 and abs(binf - bt1) < 1e-8:
        return (np.array([[-mu]]), np.array([[mu]]))

    if not (binf > bt1 and bt1 > 1):
        return (np.array([[-mu]]), np.array([[mu]]))

    # d = r1 + r2
    c = (binf - 1) / (binf - bt1)
    # Solve ProductLog: z = w * exp(w) where z = -c * exp(-c)
    z = -c * np.exp(-c)
    w = fsolve(lambda w: z - w * np.exp(w), 1.0, full_output=False)[0]
    x = (w + c) / t1

    k1 = mu**3 * t2**3
    k2 = 3 * mu**2 * (binf - 1) * t2**2
    k3 = 3 * mu * (binf - 1) / x * t2
    k4 = 3 * mu / x**2 * (binf - 1) * t2 * np.exp(-x * t2)
    k5 = 6 * mu / x**3 * (binf - 1) * (1 - np.exp(-x * t2))
    g1t2 = m3t2 + 3 * mu * t2 * (mu * t2 - 1) * bt2 + mu * t2 * (mu * t2 - 1) * (mu * t2 - 2)
    h = (g1t2 - k1 - k2 - k3 * (-mu) - k4 * mu * x) / ((k3 / x) + k4 - k5)

    if abs(h) < 1e-4:
        r1 = x / 2
        r2 = x / 2
        l2 = mu - 0.5 * np.sqrt(2 * (binf - 1) * mu * x)
        l1 = mu + 0.5 * np.sqrt(2 * (binf - 1) * mu * x)
    else:
        y = (binf - 1) * mu * x**3 / (2 * h**2)
        r1 = x / 2 * (1 + 1 / np.sqrt(4 * y + 1))
        r2 = x - r1
        if r1 < r2:
            r1, r2 = r2, r1
        w_val = h / (r1 - r2)
        w_min = -mu / r1 * (r1 + r2)
        w_max = mu / r2 * (r1 + r2)
        if w_val < w_min or w_val > w_max:
            z_val = (binf - 1) * x**3 * mu
            u = x * z_val / (2 * mu**2 * x**2 + z_val)
            r1 = u + (x - u) / 2
            r2 = x - r1
            delta = np.sqrt(z_val / (2 * r1 * r2))
            l2 = mu - r2 / x * delta
            l1 = l2 + delta
        else:
            l2 = mu - h / (r1 - r2) * (r2 / (r1 + r2))
            l1 = h / (r1 - r2) + l2

    D0 = np.array([[-(r1 + l1), r1],
                   [r2, -(r2 + l2)]])
    D1 = np.array([[l1, 0.0],
                   [0.0, l2]])
    return (D0, D1)


def mmpp2_fitc_approx(a: float, bt1: float, bt2: float, binf: float,
                      m3t2: float, t1: float, t2: float
                      ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit MMPP2 from counting process statistics using optimization.

    Matches MATLAB: mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)

    Args:
        a: Arrival rate
        bt1: IDC at scale t1
        bt2: IDC at scale t2
        binf: IDC for t->inf
        m3t2: Third central moment at scale t2
        t1: First time scale
        t2: Second time scale

    Returns:
        Tuple of (D0, D1) matrices
    """
    from ..mam import map_scale
    from scipy.optimize import minimize

    def compute_obj(params):
        l1, l2, r1, r2 = params
        if l1 <= 0 or l2 <= 0 or r1 < 0 or r2 < 0:
            return 1e10

        xa = (l1 * r2 + l2 * r1) / (r1 + r2)
        if xa <= 0:
            return 1e10
        factor = a / xa

        exp_val_t1 = np.exp(-(r1 * t1 * factor + r2 * t1 * factor))
        xbt1 = (r1 * (2 * l1**2 * r2**2 * t1 * factor - 2 * l2**2 * r2 - 2 * l1**2 * r2 +
                      2 * l2**2 * r2**2 * t1 * factor + 4 * l1 * l2 * r2 +
                      2 * l1**2 * r2 * exp_val_t1 + 2 * l2**2 * r2 * exp_val_t1 -
                      4 * l1 * l2 * r2**2 * t1 * factor - 4 * l1 * l2 * r2 * exp_val_t1) +
                r1**2 * (2 * r2 * t1 * factor * l1**2 - 4 * r2 * t1 * factor * l1 * l2 +
                        2 * r2 * t1 * factor * l2**2)) / \
               (t1 * factor * (r1 + r2)**3 * (l1 * r2 + l2 * r1)) + 1

        if abs(t1 - t2) > 1e-10:
            exp_val_t2 = np.exp(-(r1 * t2 * factor + r2 * t2 * factor))
            xbt2 = (r1 * (2 * l1**2 * r2**2 * t2 * factor - 2 * l2**2 * r2 - 2 * l1**2 * r2 +
                          2 * l2**2 * r2**2 * t2 * factor + 4 * l1 * l2 * r2 +
                          2 * l1**2 * r2 * exp_val_t2 + 2 * l2**2 * r2 * exp_val_t2 -
                          4 * l1 * l2 * r2**2 * t2 * factor - 4 * l1 * l2 * r2 * exp_val_t2) +
                    r1**2 * (2 * r2 * t2 * factor * l1**2 - 4 * r2 * t2 * factor * l1 * l2 +
                            2 * r2 * t2 * factor * l2**2)) / \
                   (t2 * factor * (r1 + r2)**3 * (l1 * r2 + l2 * r1)) + 1
        else:
            xbt2 = xbt1

        xbinf = ((2 * r2 * l1**2 - 4 * r2 * l1 * l2 + 2 * r2 * l2**2) * r1**2 +
                 (2 * l1**2 * r2**2 - 4 * l1 * l2 * r2**2 + 2 * l2**2 * r2**2) * r1) / \
                ((r1 + r2)**3 * (l1 * r2 + l2 * r1)) + 1

        t = t2 * factor
        d = r1 + r2
        p = (l1 - l2) * (r1 - r2)
        xg3t = (xa**3 * t**3 +
                3 * xa**2 * (xbinf - 1) * t**2 +
                3 * xa * (xbinf - 1) / d * (p / d - xa) * t +
                3 * xa / d**2 * (xbinf - 1) * (p + xa * d) * t * np.exp(-t * d) -
                6 * xa / d**3 * (xbinf - 1) * p * (1 - np.exp(-t * d)))
        xm3t2 = xg3t - 3 * xa * t * (xa * t - 1) * xbt2 - xa * t * (xa * t - 1) * (xa * t - 2)

        obj = 0.0
        obj += (xa / a - 1)**2
        obj += (xbt1 / bt1 - 1)**2
        if abs(t1 - t2) > 1e-10:
            obj += (xbt2 / bt2 - 1)**2
        obj += (xbinf / binf - 1)**2
        if abs(m3t2) > 1e-10:
            obj += (xm3t2 / m3t2 - 1)**2

        return obj

    x0 = np.array([a * 0.75, a * 1.5, 1.0 / 3.0, 2.0 / 3.0])
    bounds = [(1e-6, None), (1e-6, None), (0, None), (0, None)]
    result = minimize(compute_obj, x0, method='L-BFGS-B', bounds=bounds,
                     options={'maxiter': 10000, 'ftol': 1e-12})

    l1, l2, r1, r2 = result.x
    D0 = np.array([[-(l1 + r1), r1],
                   [r2, -(l2 + r2)]])
    D1 = np.array([[l1, 0.0],
                   [0.0, l2]])

    # Scale to match exact rate
    D0, D1 = map_scale(D0, D1, 1.0 / a)
    return (D0, D1)


def mmpp2_fitc_theoretical(MAP, t1: float = 1.0, t2: float = 10.0, tinf: float = 1e8
                           ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit theoretical characteristics of a MAP(n) with a MMPP(2).

    Matches MATLAB: mmpp2_fitc_theoretical(map, t1, t2, tinf)

    Args:
        MAP: Input MAP as (D0, D1) tuple
        t1: First time scale (default 1)
        t2: Second time scale (default 10)
        tinf: Large time scale for asymptotic IDC (default 1e8)

    Returns:
        Tuple of (D0, D1) matrices for the fitted MMPP2
    """
    from ..mam import map_count_mean, map_count_var, map_count_moment

    D0, D1 = MAP[0], MAP[1]

    a_val = map_count_mean(D0, D1, np.array([t1]))[0] / t1
    bt1_val = map_count_var(D0, D1, np.array([t1]))[0] / (a_val * t1)
    bt2_val = map_count_var(D0, D1, np.array([t2]))[0] / (a_val * t2)
    binf_val = map_count_var(D0, D1, np.array([tinf]))[0] / (a_val * tinf)
    mt2 = map_count_moment(D0, D1, t2, [1, 2, 3])
    m3t2_val = mt2[2] - 3 * mt2[1] * mt2[0] + 2 * mt2[0]**3

    return mmpp2_fitc(a_val, bt1_val, bt2_val, binf_val, m3t2_val, t1, t2)


def mmpp_rand(K: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a random MMPP with K states.

    Matches MATLAB: mmpp_rand(K) - generates random D0, D1 (diagonal), normalizes.

    Args:
        K: Number of states

    Returns:
        Tuple of (D0, D1) matrices
    """
    from ..mam import map_normalize

    D0 = np.random.rand(K, K)
    D1 = np.diag(np.diag(np.random.rand(K, K)))
    return map_normalize(D0, D1)


__all__ = [
    'mmpp2_fit',
    'mmpp2_fit1',
    'mmpp2_fit2',
    'mmpp2_fit3',
    'mmpp2_fit4',
    'mmpp2_fitc',
    'mmpp2_fitc_approx',
    'mmpp2_fitc_theoretical',
    'mmpp_rand',
]
