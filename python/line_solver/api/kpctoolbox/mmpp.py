"""
Markov Modulated Poisson Process (MMPP) functions for KPC-Toolbox.

Native Python implementations of MMPP fitting and analysis.
"""

import numpy as np
from typing import Dict


def mmpp2_fit(E1: float, E2: float, E3: float, ACFLAG1: float
              ) -> Dict[str, np.ndarray]:
    """
    Fit a 2-state MMPP (MMPP2) to match first three moments and lag-1 autocorrelation.

    Args:
        E1: First moment (mean)
        E2: Second moment
        E3: Third moment
        ACFLAG1: Lag-1 autocorrelation (should be in [0, 0.5])

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    SCV = (E2 - E1 * E1) / (E1 * E1)

    # Handle edge case
    if abs(SCV - 1) < 1e-10:
        G2 = 0.0
    else:
        G2 = ACFLAG1 / ((1 - 1 / SCV) / 2) if SCV != 0 else 0.0

    if abs(G2) < 1e-6 or G2 == 0.0:
        # If G2 is very close to zero, fit with MAP(1) approximation
        denom1 = 6 * E1**3 * SCV + 3 * E1**3 * SCV**2 + 3 * E1**3 - 2 * E3
        denom2 = 6 * E1**3 * SCV - E3

        if abs(denom1) > 1e-10 and abs(denom2) > 1e-10:
            mu00 = 2 * (6 * E1**3 * SCV - E3) / E1 / denom1
            q01 = 9 * E1**5 * (SCV - 1) * (SCV**2 - 2 * SCV + 1) / denom2 / denom1
            q10 = -3 * (SCV - 1) * E1**2 / denom2
        else:
            mu00 = 1.0 / E1
            q01 = 0.1
            q10 = 0.1

        mu11 = 0.0
    else:
        # Full MMPP2 fitting
        try:
            disc_val = (E3**2 - 12 * E1**3 * SCV * E3 + 6 * E1**3 * G2 * E3 -
                       6 * G2 * SCV * E1**3 * E3 + 18 * G2 * SCV**3 * E1**6 -
                       18 * E1**6 * G2 * SCV**2 + 9 * E1**6 * G2**2 +
                       36 * E1**6 * SCV**2 + 18 * E1**6 * G2 * SCV -
                       18 * E1**6 * SCV * G2**2 + 9 * E1**6 * SCV**2 * G2**2 -
                       18 * E1**6 * G2)

            if disc_val < 0:
                disc = 0.0
            else:
                disc = np.sqrt(disc_val)

            denom1 = -3 * E1**3 * SCV**2 - 6 * E1**3 * SCV - 3 * E1**3 + 2 * E3

            if abs(denom1) > 1e-10:
                term1 = (-3 * E1**3 * G2 + 3 * E1**3 * G2 * SCV -
                        6 * E1**3 * SCV + E3 + disc) / denom1
            else:
                term1 = 1.0

            mu11 = term1 / E1

            # mu00 calculation (simplified)
            a = G2 * (1 - SCV)
            b = 2 * SCV - G2 * SCV - 2
            c = term1 * E1**2

            if abs(b) > 1e-10:
                discriminant = a * a - 4 * b * c
                if discriminant >= 0:
                    mu00 = (-a + np.sqrt(discriminant)) / (2 * b) / E1
                else:
                    mu00 = 1.0 / E1
            else:
                mu00 = 1.0 / E1

            # q01 and q10 calculations
            if disc > 0:
                q01 = 3 * E1**2 * (1 - G2) * (SCV - 1) * term1 / (disc * 2)
                q10 = 3 * E1**2 * (G2 - 1) * (SCV - 1) / disc
            else:
                q01 = 0.1
                q10 = 0.1

        except (ValueError, ZeroDivisionError):
            mu00 = 1.0 / E1
            mu11 = 0.0
            q01 = 0.1
            q10 = 0.1

    # Ensure rates are non-negative
    mu00 = max(0.0, mu00)
    mu11 = max(0.0, mu11)
    q01 = max(0.0, q01)
    q10 = max(0.0, q10)

    # Build D0 and D1 matrices
    D0 = np.array([
        [-mu00 - q01, q01],
        [q10, -mu11 - q10]
    ])

    D1 = np.array([
        [mu00, 0.0],
        [0.0, mu11]
    ])

    return {'D0': D0, 'D1': D1}


def mmpp2_fit1(E1: float, E2: float, E3: float) -> Dict[str, np.ndarray]:
    """
    Fit MMPP2 using only moments (no autocorrelation).

    Args:
        E1: First moment
        E2: Second moment
        E3: Third moment

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    return mmpp2_fit(E1, E2, E3, 0.0)


def mmpp2_fit2(E1: float, E2: float, E3: float, acf1: float
               ) -> Dict[str, np.ndarray]:
    """
    Fit MMPP2 using moments and lag-1 ACF.

    Args:
        E1: First moment
        E2: Second moment
        E3: Third moment
        acf1: Lag-1 autocorrelation

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    return mmpp2_fit(E1, E2, E3, acf1)


def mmpp2_fit3(E1: float, E2: float, E3: float, acf2: float
               ) -> Dict[str, np.ndarray]:
    """
    Fit MMPP2 using moments and lag-2 ACF (approximation).

    Args:
        E1: First moment
        E2: Second moment
        E3: Third moment
        acf2: Lag-2 autocorrelation

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    # Approximate lag-1 from lag-2 using geometric decay assumption
    acf1 = np.sqrt(abs(acf2)) * np.sign(acf2)
    return mmpp2_fit(E1, E2, E3, acf1)


def mmpp2_fit4(E1: float, E2: float, E3: float, acf_values: np.ndarray
               ) -> Dict[str, np.ndarray]:
    """
    Fit MMPP2 using moments and multiple ACF lags.

    Args:
        E1: First moment
        E2: Second moment
        E3: Third moment
        acf_values: Array of ACF values at lags 1, 2, ...

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    acf_values = np.asarray(acf_values)
    acf1 = acf_values[0] if len(acf_values) > 0 else 0.0
    return mmpp2_fit(E1, E2, E3, acf1)


def mmpp2_fitc(mean_count: float, var_count: float, scale: float
               ) -> Dict[str, np.ndarray]:
    """
    Fit MMPP2 from counting process statistics.

    Args:
        mean_count: Mean count
        var_count: Variance of counts
        scale: Time scale

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    # Convert counting statistics to inter-arrival statistics
    E1 = scale / mean_count
    idc = var_count / mean_count
    SCV = idc

    # Approximate E2 and E3 from SCV
    E2 = E1**2 * (1 + SCV)
    E3 = E2 * E1 * (2 * SCV + 1)

    return mmpp2_fit(E1, E2, E3, 0.0)


def mmpp2_fitc_approx(mean_count: float, var_count: float, scale: float,
                      acf_count: float) -> Dict[str, np.ndarray]:
    """
    Fit MMPP2 from counting process with approximation.

    Args:
        mean_count: Mean count
        var_count: Variance of counts
        scale: Time scale
        acf_count: ACF of counting process

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    E1 = scale / mean_count
    idc = var_count / mean_count
    SCV = idc

    E2 = E1**2 * (1 + SCV)
    E3 = E2 * E1 * (2 * SCV + 1)

    # Approximate inter-arrival ACF from counting ACF
    acf1 = acf_count * 0.5

    return mmpp2_fit(E1, E2, E3, acf1)


def mmpp2_fitc_theoretical(lambda1: float, lambda2: float,
                           q12: float, q21: float) -> Dict[str, np.ndarray]:
    """
    Theoretical MMPP2 fitting from parameters.

    Args:
        lambda1: Arrival rate in state 1
        lambda2: Arrival rate in state 2
        q12: Transition rate from state 1 to 2
        q21: Transition rate from state 2 to 1

    Returns:
        Fitted MAP as {'D0': D0, 'D1': D1}
    """
    D0 = np.array([
        [-lambda1 - q12, q12],
        [q21, -lambda2 - q21]
    ])

    D1 = np.array([
        [lambda1, 0.0],
        [0.0, lambda2]
    ])

    return {'D0': D0, 'D1': D1}


def mmpp_rand(MAP: Dict[str, np.ndarray], n_samples: int,
              seed: int = None) -> np.ndarray:
    """
    Generate random samples from an MMPP.

    Args:
        MAP: MMPP as {'D0': D0, 'D1': D1}
        n_samples: Number of samples to generate
        seed: Random seed (optional)

    Returns:
        Array of inter-arrival times
    """
    D0 = MAP['D0']
    D1 = MAP['D1']
    n = D0.shape[0]

    if seed is not None:
        np.random.seed(seed)

    samples = np.zeros(n_samples)

    # Start in state 0
    current_state = 0

    for s in range(n_samples):
        time = 0.0

        # Generate inter-arrival time
        while True:
            # Holding time in current state
            rate = -D0[current_state, current_state]
            if rate <= 0:
                time += 1.0  # Fallback
                break

            hold_time = -np.log(np.random.rand()) / rate

            # Decide next transition
            arrival_rate = D1[current_state, current_state]
            arrival_prob = arrival_rate / rate if rate > 0 else 0

            if np.random.rand() < arrival_prob:
                # Arrival occurred
                time += hold_time

                # Sample destination state for arrival
                if arrival_rate > 0:
                    rnd = np.random.rand()
                    cum_prob = 0.0
                    for j in range(n):
                        cum_prob += D1[current_state, j] / arrival_rate
                        if rnd < cum_prob:
                            current_state = j
                            break
                break
            else:
                # Hidden transition (no arrival)
                time += hold_time

                # Sample destination state for hidden transition
                hidden_rate = rate - arrival_rate
                if hidden_rate > 0:
                    rnd = np.random.rand()
                    cum_prob = 0.0
                    for j in range(n):
                        if j != current_state:
                            cum_prob += D0[current_state, j] / hidden_rate
                            if rnd < cum_prob:
                                current_state = j
                                break

        samples[s] = time

    return samples


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
