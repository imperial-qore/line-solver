"""
Markov-Modulated Deterministic Process (MMDP) for fluid queue modeling.

Native Python implementation of MMDP, a process characterized by deterministic
fluid flow rates modulated by a background continuous-time Markov chain (CTMC).

MMDP is the deterministic analogue of MMPP:
- MMPP: Poisson arrival rates modulated by a Markov chain
- MMDP: Deterministic flow rates modulated by a Markov chain

The (Q, R) parameterization follows BUTools conventions:
- Q: Generator matrix of the modulating CTMC (row sums = 0)
- R: Diagonal matrix of deterministic rates per state

Key Classes:
    MMDP: General n-state Markov-Modulated Deterministic Process
    MMDP2: Specialized 2-state MMDP with convenient parameterization

Key Functions:
    mmdp_isfeasible: Check if (Q, R) defines a valid MMDP
    mmdp_mean_rate: Compute stationary mean rate
    mmdp_scv: Compute squared coefficient of variation
    mmdp_from_map: Convert MAP to MMDP representation

Usage:
    from line_solver.api.mmdp import MMDP, MMDP2, mmdp_isfeasible

    # Create 2-state MMDP
    Q = np.array([[-0.5, 0.5], [0.3, -0.3]])
    R = np.array([[2.0, 0], [0, 5.0]])
    mmdp = MMDP(Q, R)
    print(f"Mean rate: {mmdp.get_mean_rate()}")

    # Or use MMDP2 with convenient parameterization
    mmdp2 = MMDP2(r0=2.0, r1=5.0, sigma0=0.5, sigma1=0.3)
    print(f"Mean rate: {mmdp2.get_mean_rate()}")

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg
from typing import Tuple, Optional, Union

# Import CTMC solver from mc module
from ..mc import ctmc_solve


def mmdp_isfeasible(
    Q: np.ndarray,
    R: np.ndarray,
    tol: float = 1e-10
) -> bool:
    """
    Check if (Q, R) defines a valid MMDP.

    Requirements:
    - Q must be a valid generator (square, row sums = 0, proper signs)
    - R must be diagonal with non-negative entries
    - Q and R must have compatible dimensions

    Args:
        Q: Generator matrix (n x n)
        R: Rate matrix (n x n diagonal) or rate vector (n,)
        tol: Numerical tolerance for validation

    Returns:
        True if (Q, R) defines a valid MMDP, False otherwise

    Examples:
        >>> Q = np.array([[-0.5, 0.5], [0.3, -0.3]])
        >>> R = np.array([[2.0, 0], [0, 5.0]])
        >>> mmdp_isfeasible(Q, R)
        True
    """
    Q = np.asarray(Q, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)

    # Handle scalar case (single state)
    if Q.ndim == 0:
        Q = np.array([[Q]])
    if R.ndim == 0:
        R = np.array([[R]])

    # Convert vector R to diagonal matrix
    if R.ndim == 1:
        R = np.diag(R)

    n = Q.shape[0]

    # Q must be square
    if Q.shape[0] != Q.shape[1]:
        return False

    # R must be n x n
    if R.shape[0] != n or R.shape[1] != n:
        return False

    # Q must be a valid generator
    for i in range(n):
        # Diagonal elements must be non-positive
        if Q[i, i] > tol:
            return False
        # Off-diagonal elements must be non-negative
        for j in range(n):
            if i != j and Q[i, j] < -tol:
                return False
        # Row sums must be zero
        if abs(np.sum(Q[i, :])) > tol:
            return False

    # R must be diagonal with non-negative entries
    for i in range(n):
        # Diagonal entries must be non-negative
        if R[i, i] < -tol:
            return False
        # Off-diagonal entries must be zero
        for j in range(n):
            if i != j and abs(R[i, j]) > tol:
                return False

    return True


def mmdp_mean_rate(
    Q: np.ndarray,
    R: np.ndarray
) -> float:
    """
    Compute the stationary mean rate of an MMDP.

    The mean rate is E[r] = pi * diag(R), where pi is the stationary
    distribution of the modulating CTMC with generator Q.

    Args:
        Q: Generator matrix (n x n)
        R: Rate matrix (n x n diagonal) or rate vector (n,)

    Returns:
        Stationary mean deterministic rate

    Examples:
        >>> Q = np.array([[-0.5, 0.5], [0.3, -0.3]])
        >>> R = np.array([[2.0, 0], [0, 5.0]])
        >>> mmdp_mean_rate(Q, R)
        3.125
    """
    Q = np.asarray(Q, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)

    # Handle scalar case
    if Q.ndim == 0 or (Q.ndim == 2 and Q.shape[0] == 1):
        if R.ndim == 0:
            return float(R)
        elif R.ndim == 1:
            return float(R[0])
        else:
            return float(R[0, 0])

    # Convert vector R to diagonal matrix
    if R.ndim == 1:
        r = R
    else:
        r = np.diag(R)

    # Compute stationary distribution
    pi = ctmc_solve(Q)

    # Mean rate = pi * r
    return float(np.dot(pi, r))


def mmdp_scv(
    Q: np.ndarray,
    R: np.ndarray
) -> float:
    """
    Compute the squared coefficient of variation of MMDP rates.

    Computes Var[r] / E[r]^2 where expectation is over the
    stationary distribution of the modulating CTMC.

    Args:
        Q: Generator matrix (n x n)
        R: Rate matrix (n x n diagonal) or rate vector (n,)

    Returns:
        Squared coefficient of variation

    Examples:
        >>> Q = np.array([[-0.5, 0.5], [0.3, -0.3]])
        >>> R = np.array([[2.0, 0], [0, 5.0]])
        >>> mmdp_scv(Q, R)  # Returns SCV of rates
    """
    Q = np.asarray(Q, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)

    # Handle scalar case (no variability)
    if Q.ndim == 0 or (Q.ndim == 2 and Q.shape[0] == 1):
        return 0.0

    # Convert vector R to rate vector
    if R.ndim == 1:
        r = R
    else:
        r = np.diag(R)

    # Compute stationary distribution
    pi = ctmc_solve(Q)

    # Mean rate
    mean_rate = np.dot(pi, r)

    # E[r^2]
    mean_rate_sq = np.dot(pi, r * r)

    # Variance
    var_rate = mean_rate_sq - mean_rate * mean_rate

    # SCV
    if mean_rate > 0:
        return float(var_rate / (mean_rate * mean_rate))
    else:
        return float('inf')


def mmdp_from_map(
    D0: np.ndarray,
    D1: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a MAP to MMDP representation.

    Extracts the full generator Q = D0 + D1 and uses row sums of D1
    as the deterministic rates.

    Args:
        D0: MAP's D0 matrix (transitions without arrivals)
        D1: MAP's D1 matrix (transitions with arrivals)

    Returns:
        Tuple (Q, R) where:
        - Q: Generator matrix (D0 + D1)
        - R: Diagonal rate matrix (row sums of D1)

    Examples:
        >>> D0 = np.array([[-2.0, 0.5], [0.3, -1.5]])
        >>> D1 = np.array([[1.0, 0.5], [0.4, 0.8]])
        >>> Q, R = mmdp_from_map(D0, D1)
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    # Q = D0 + D1
    Q = D0 + D1

    # R = diag(row sums of D1)
    n = D1.shape[0]
    R = np.diag(np.sum(D1, axis=1))

    return Q, R


class MMDP:
    """
    Markov-Modulated Deterministic Process for fluid queue modeling.

    Models fluid flow with deterministic rates modulated by a background
    Markov chain. Suitable for arrival and service processes in Markovian
    fluid queues.

    The (Q, R) parameterization follows BUTools conventions:
    - Q: Generator matrix of the modulating CTMC (row sums = 0)
    - R: Diagonal matrix of deterministic rates per state

    Attributes:
        Q: Generator matrix (n x n)
        R: Rate matrix (n x n diagonal)
        n_phases: Number of phases in the modulating CTMC

    Examples:
        >>> Q = np.array([[-0.5, 0.5], [0.3, -0.3]])
        >>> R = np.array([[2.0, 0], [0, 5.0]])
        >>> mmdp = MMDP(Q, R)
        >>> mmdp.get_mean_rate()
        3.125
    """

    def __init__(
        self,
        Q: np.ndarray,
        R: np.ndarray,
        validate: bool = True
    ):
        """
        Create an MMDP with specified generator Q and rate matrix R.

        Args:
            Q: n x n generator matrix (row sums = 0)
            R: n x n diagonal matrix of deterministic rates, or n-vector
            validate: If True, validate that (Q, R) is feasible

        Raises:
            ValueError: If validate=True and (Q, R) is not a valid MMDP
        """
        self._Q = np.asarray(Q, dtype=np.float64)
        R_arr = np.asarray(R, dtype=np.float64)

        # Handle scalar case
        if self._Q.ndim == 0:
            self._Q = np.array([[self._Q]])
        if R_arr.ndim == 0:
            R_arr = np.array([[R_arr]])

        # Convert vector to diagonal matrix if needed
        if R_arr.ndim == 1:
            self._R = np.diag(R_arr)
        else:
            self._R = R_arr

        self.n_phases = self._Q.shape[0]

        # Cache for stationary distribution
        self._pi = None

        if validate and not mmdp_isfeasible(self._Q, self._R):
            raise ValueError("Invalid MMDP: (Q, R) does not satisfy feasibility constraints")

    @property
    def Q(self) -> np.ndarray:
        """Get the generator matrix Q."""
        return self._Q.copy()

    @property
    def R(self) -> np.ndarray:
        """Get the rate matrix R (diagonal)."""
        return self._R.copy()

    def r(self) -> np.ndarray:
        """Get the rate vector (diagonal of R)."""
        return np.diag(self._R)

    def get_number_of_phases(self) -> int:
        """Get the number of phases in the modulating CTMC."""
        return self.n_phases

    def _get_stationary(self) -> np.ndarray:
        """Get cached stationary distribution."""
        if self._pi is None:
            self._pi = ctmc_solve(self._Q)
        return self._pi

    def get_mean_rate(self) -> float:
        """
        Compute the stationary mean rate.

        For MMDP: E[r] = pi * diag(R), where pi is the stationary distribution
        of the modulating CTMC with generator Q.

        Returns:
            Stationary mean deterministic rate
        """
        pi = self._get_stationary()
        r = np.diag(self._R)
        return float(np.dot(pi, r))

    def get_mean(self) -> float:
        """
        Compute the mean inter-event time (inverse of mean rate).

        Returns:
            Mean inter-event time (1 / mean_rate), or inf if mean_rate = 0
        """
        rate = self.get_mean_rate()
        if rate > 0:
            return 1.0 / rate
        return float('inf')

    def get_rate(self) -> float:
        """
        Get the rate (alias for get_mean_rate).

        Returns:
            Stationary mean rate
        """
        return self.get_mean_rate()

    def get_scv(self) -> float:
        """
        Compute the squared coefficient of variation of rates.

        Computes Var[r] / E[r]^2 where expectation is over the
        stationary distribution of the modulating CTMC.

        Returns:
            Squared coefficient of variation
        """
        pi = self._get_stationary()
        r = np.diag(self._R)

        mean_rate = np.dot(pi, r)
        mean_rate_sq = np.dot(pi, r * r)
        var_rate = mean_rate_sq - mean_rate * mean_rate

        if mean_rate > 0:
            return float(var_rate / (mean_rate * mean_rate))
        return float('inf')

    def get_process(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the (Q, R) representation.

        Returns:
            Tuple (Q, R) of generator and rate matrices
        """
        return (self.Q, self.R)

    def is_feasible(self) -> bool:
        """
        Check if this MMDP is feasible.

        Returns:
            True if (Q, R) satisfies all MMDP constraints
        """
        return mmdp_isfeasible(self._Q, self._R)

    @staticmethod
    def from_map(D0: np.ndarray, D1: np.ndarray) -> 'MMDP':
        """
        Convert a MAP to MMDP (deterministic representation).

        Converts a Markovian Arrival Process to a Markov-Modulated
        Deterministic Process by extracting the full generator and
        using row sums of D1 as the deterministic rates.

        Args:
            D0: MAP's D0 matrix (transitions without arrivals)
            D1: MAP's D1 matrix (transitions with arrivals)

        Returns:
            MMDP representation of the MAP
        """
        Q, R = mmdp_from_map(D0, D1)
        return MMDP(Q, R, validate=False)

    @staticmethod
    def from_mmpp2(
        lambda0: float,
        lambda1: float,
        sigma0: float,
        sigma1: float
    ) -> 'MMDP':
        """
        Create MMDP from MMPP2 parameters.

        Creates a 2-state MMDP using the same parameterization as MMPP2.

        Args:
            lambda0: Rate in state 0
            lambda1: Rate in state 1
            sigma0: Transition rate from state 0 to state 1
            sigma1: Transition rate from state 1 to state 0

        Returns:
            2-state MMDP
        """
        Q = np.array([[-sigma0, sigma0], [sigma1, -sigma1]])
        R = np.diag([lambda0, lambda1])
        return MMDP(Q, R, validate=False)

    def __repr__(self) -> str:
        mean_rate = self.get_mean_rate()
        scv = self.get_scv()
        return f"MMDP({self.n_phases}x{self.n_phases}, mean_rate={mean_rate:.6f}, scv={scv:.6f})"


class MMDP2(MMDP):
    """
    2-state Markov-Modulated Deterministic Process.

    A specialized MMDP with exactly 2 phases, using a convenient
    parameterization analogous to MMPP2.

    Parameterization:
        r0, r1: Deterministic rates in states 0 and 1
        sigma0: Transition rate from state 0 to state 1
        sigma1: Transition rate from state 1 to state 0

    The generator matrix is:
        Q = [[-sigma0, sigma0], [sigma1, -sigma1]]

    The rate matrix is:
        R = diag([r0, r1])

    Attributes:
        r0: Rate in state 0
        r1: Rate in state 1
        sigma0: Transition rate from state 0 to state 1
        sigma1: Transition rate from state 1 to state 0

    Examples:
        >>> mmdp2 = MMDP2(r0=2.0, r1=5.0, sigma0=0.5, sigma1=0.3)
        >>> mmdp2.get_mean_rate()
        3.125
    """

    def __init__(
        self,
        r0: float,
        r1: float,
        sigma0: float,
        sigma1: float
    ):
        """
        Create a 2-state MMDP.

        Args:
            r0: Deterministic rate in state 0
            r1: Deterministic rate in state 1
            sigma0: Transition rate from state 0 to state 1
            sigma1: Transition rate from state 1 to state 0
        """
        self.r0 = float(r0)
        self.r1 = float(r1)
        self.sigma0 = float(sigma0)
        self.sigma1 = float(sigma1)

        Q = np.array([[-sigma0, sigma0], [sigma1, -sigma1]])
        R = np.diag([r0, r1])

        super().__init__(Q, R, validate=False)

    def get_mean_rate(self) -> float:
        """
        Compute stationary mean rate (closed-form).

        For a 2-state MMDP, the mean rate has the closed form:
            E[r] = (r0 * sigma1 + r1 * sigma0) / (sigma0 + sigma1)

        Returns:
            Stationary mean deterministic rate
        """
        return (self.r0 * self.sigma1 + self.r1 * self.sigma0) / (self.sigma0 + self.sigma1)

    def get_scv(self) -> float:
        """
        Compute squared coefficient of variation (closed-form).

        Returns:
            Squared coefficient of variation
        """
        # Stationary probabilities
        pi0 = self.sigma1 / (self.sigma0 + self.sigma1)
        pi1 = self.sigma0 / (self.sigma0 + self.sigma1)

        # Mean and variance
        mean_rate = pi0 * self.r0 + pi1 * self.r1
        var_rate = pi0 * self.r0**2 + pi1 * self.r1**2 - mean_rate**2

        if mean_rate > 0:
            return float(var_rate / (mean_rate * mean_rate))
        return float('inf')

    def get_number_of_phases(self) -> int:
        """Always returns 2 for MMDP2."""
        return 2

    def __repr__(self) -> str:
        mean_rate = self.get_mean_rate()
        return f"MMDP2(r0={self.r0}, r1={self.r1}, sigma0={self.sigma0}, sigma1={self.sigma1}, mean_rate={mean_rate:.6f})"


# Module exports
__all__ = [
    'MMDP',
    'MMDP2',
    'mmdp_isfeasible',
    'mmdp_mean_rate',
    'mmdp_scv',
    'mmdp_from_map',
]
