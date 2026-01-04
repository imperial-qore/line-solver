"""
Markovian distributions for LINE (pure Python).

This module provides Markovian (phase-type) distribution implementations
including PH, APH, Coxian, and MAP (Markovian Arrival Process).
"""

from typing import Optional, Tuple, Union
import numpy as np
from scipy import linalg

from .base import ContinuousDistribution, Markovian


class PH(ContinuousDistribution, Markovian):
    """
    Phase-type distribution.

    A phase-type distribution is defined by an initial probability
    vector alpha and a sub-generator matrix T. The distribution
    represents the time until absorption in a continuous-time
    Markov chain.

    Args:
        alpha: Initial probability vector (1 x n).
        T: Sub-generator matrix (n x n). Must have negative diagonal
           and non-negative off-diagonal elements.
    """

    def __init__(self, alpha: Union[list, np.ndarray], T: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'PH'
        self._alpha = np.atleast_1d(np.array(alpha, dtype=float))
        self._T = np.atleast_2d(np.array(T, dtype=float))

        n = len(self._alpha)
        if self._T.shape != (n, n):
            raise ValueError(f"T must be {n}x{n} to match alpha of length {n}")

        # Compute exit rate vector
        self._t = -self._T.sum(axis=1)

    @property
    def alpha(self) -> np.ndarray:
        """Get the initial probability vector."""
        return self._alpha.copy()

    @property
    def T(self) -> np.ndarray:
        """Get the sub-generator matrix."""
        return self._T.copy()

    @property
    def t(self) -> np.ndarray:
        """Get the exit rate vector."""
        return self._t.copy()

    def getMean(self) -> float:
        """Get the mean."""
        # E[X] = -alpha * T^(-1) * e
        try:
            T_inv = linalg.inv(self._T)
            e = np.ones(len(self._alpha))
            return float(-self._alpha @ T_inv @ e)
        except linalg.LinAlgError:
            return float('nan')

    def getVar(self) -> float:
        """Get the variance."""
        try:
            T_inv = linalg.inv(self._T)
            e = np.ones(len(self._alpha))
            mean = float(-self._alpha @ T_inv @ e)
            second_moment = float(2 * self._alpha @ T_inv @ T_inv @ e)
            return second_moment - mean ** 2
        except linalg.LinAlgError:
            return float('nan')

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return len(self._alpha)

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix (equals T)."""
        return self._T.copy()

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix."""
        return np.outer(self._t, self._alpha)

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        return -np.diag(self._T)

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        mu = self.getMu()
        return np.where(mu > 0, self._t / mu, 0)

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        return self._alpha.copy()

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        e = np.ones(len(self._alpha))
        expm = linalg.expm(self._T * x)
        return 1.0 - float(self._alpha @ expm @ e)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        expm = linalg.expm(self._T * x)
        return float(self._alpha @ expm @ self._t)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples using simulation."""
        if rng is None:
            rng = np.random.default_rng()

        samples = np.zeros(n)
        num_phases = len(self._alpha)

        for i in range(n):
            # Start in initial phase
            phase = rng.choice(num_phases, p=self._alpha)
            time = 0.0

            while True:
                # Time in current phase
                rate = -self._T[phase, phase]
                time += rng.exponential(1.0 / rate)

                # Transition probabilities
                probs = np.zeros(num_phases + 1)
                probs[:num_phases] = self._T[phase, :].copy()
                probs[phase] = 0  # No self-loops counted
                probs[num_phases] = self._t[phase]  # Absorption
                probs = np.maximum(probs, 0)

                total = probs.sum()
                if total <= 0:
                    break
                probs /= total

                # Next phase or absorption
                next_state = rng.choice(num_phases + 1, p=probs)
                if next_state == num_phases:  # Absorbed
                    break
                phase = next_state

            samples[i] = time

        return samples



class APH(PH):
    """
    Acyclic Phase-Type distribution.

    An APH distribution is a phase-type distribution where the
    underlying Markov chain has an acyclic structure (upper triangular T).

    Can be constructed in two ways:
    1. From matrices: APH(alpha, T) - like PH
    2. From moments: APH(mean, scv=1.0, skew=None) - moment matching

    Args (matrix form):
        alpha: Initial probability vector (1 x n).
        T: Sub-generator matrix (n x n).

    Args (moment form):
        mean: Target mean (must be a positive scalar).
        scv: Target squared coefficient of variation (default: 1.0).
        skew: Target skewness (optional, for 3-moment matching).
    """

    def __init__(self, alpha_or_mean: Union[list, np.ndarray, float] = None,
                 T_or_scv: Union[list, np.ndarray, float] = 1.0,
                 skew: Optional[float] = None,
                 mean: Optional[float] = None,
                 scv: Optional[float] = None):
        """Initialize APH from matrices or moments."""
        # Support for keyword-only moment matching: APH(mean=..., scv=...)
        if mean is not None:
            alpha_or_mean = mean
            T_or_scv = scv if scv is not None else 1.0

        # Detect which constructor form is being used
        is_matrix_form = (isinstance(alpha_or_mean, (list, np.ndarray)) and
                          isinstance(T_or_scv, (list, np.ndarray)))

        if is_matrix_form:
            # Matrix-based construction: APH(alpha, T)
            super().__init__(alpha_or_mean, T_or_scv)
            self._name = 'APH'
            self._target_mean = None
            self._target_scv = None
            self._target_skew = None
        else:
            # Moment-based construction: APH(mean, scv, skew)
            mean = float(alpha_or_mean)
            scv = float(T_or_scv)

            if mean <= 0:
                raise ValueError("Mean must be positive")
            if scv < 0:
                raise ValueError("SCV must be non-negative")

            self._target_mean = mean
            self._target_scv = scv
            self._target_skew = skew

            # Moment matching to construct alpha and T
            alpha, T = self._moment_match(mean, scv, skew)

            super().__init__(alpha, T)
            self._name = 'APH'

    def _moment_match(self, mean: float, scv: float, skew: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Match moments to construct APH representation.

        Uses standard 2-moment or 3-moment matching algorithms.
        """
        if scv <= 0:
            # Deterministic (approximated by high-order Erlang)
            k = 1000
            rate = k / mean
            alpha = np.zeros(k)
            alpha[0] = 1.0
            T = np.diag([-rate] * k) + np.diag([rate] * (k - 1), 1)
            return alpha, T

        if scv < 1:
            # SCV < 1: Use Erlang-like structure
            k = max(1, int(np.ceil(1.0 / scv)))
            rate = k / mean
            alpha = np.zeros(k)
            alpha[0] = 1.0
            T = np.diag([-rate] * k) + np.diag([rate] * (k - 1), 1)

            # Adjust first phase for exact moment matching
            if k > 1:
                var = scv * mean ** 2
                # Solve for rates that match mean and variance
                alpha, T = self._erlang_adjust(mean, scv, k)
            return alpha, T

        elif scv == 1.0:
            # SCV = 1: Exponential
            rate = 1.0 / mean
            return np.array([1.0]), np.array([[-rate]])

        else:
            # SCV > 1: Use hyperexponential structure
            return self._hyperexp_match(mean, scv)

    def _erlang_adjust(self, mean: float, scv: float, k: int) -> Tuple[np.ndarray, np.ndarray]:
        """Adjust Erlang to match mean and SCV < 1."""
        # Use mixture of Erlang-k and Erlang-(k-1)
        # This gives SCV in range [1/k, 1/(k-1)]
        if k <= 1:
            rate = 1.0 / mean
            return np.array([1.0]), np.array([[-rate]])

        # Probability of using k phases
        k_float = float(k)
        p = k_float * (scv * k_float - 1) / (k_float - 1) if k > 1 else 1.0
        p = max(0, min(1, p))

        # Rate to match mean
        rate = (p * k + (1 - p) * (k - 1)) / mean

        # Construct mixed Erlang representation
        alpha = np.zeros(k)
        alpha[0] = p
        if k > 1:
            alpha[1] = 1 - p
        else:
            alpha[0] = 1.0

        T = np.diag([-rate] * k) + np.diag([rate] * (k - 1), 1)
        return alpha, T

    def _hyperexp_match(self, mean: float, scv: float) -> Tuple[np.ndarray, np.ndarray]:
        """Match moments using 2-phase hyperexponential."""
        # Balanced means hyperexponential
        cv = np.sqrt(scv)
        # Two phases with rates mu1, mu2 and probabilities p, 1-p
        # Use balanced means approach
        p = 0.5 * (1 + np.sqrt((scv - 1) / (scv + 1)))
        p = max(0.01, min(0.99, p))

        mu1 = 2 * p / mean
        mu2 = 2 * (1 - p) / mean

        alpha = np.array([p, 1 - p])
        T = np.diag([-mu1, -mu2])
        return alpha, T

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'APH':
        """
        Create an APH distribution from mean and SCV.

        Uses moment matching to construct an acyclic phase-type
        distribution with the specified mean and squared coefficient
        of variation.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.

        Returns:
            APH distribution with given mean and SCV.
        """
        return cls(mean, scv)  # Positional args: APH(alpha_or_mean, T_or_scv)

    # CamelCase alias
    fitMeanAndScv = fit_mean_and_scv



class Coxian(PH):
    """
    Coxian distribution.

    A Coxian distribution is a special case of phase-type distributions
    where transitions can only go to the next phase or to absorption.

    Args:
        means_or_rates: Service rates (mu) for each phase. Rates are converted
                        to means internally.
        probs: Transition probabilities (phi). If length equals number of phases,
               the last element (which should be 1.0 for absorption) is dropped.
               If length is n-1, used as-is.
    """

    def __init__(self, means_or_rates: Union[list, np.ndarray],
                 probs: Optional[Union[list, np.ndarray]] = None):
        # Flatten 2D column vectors to 1D
        rates = np.array(means_or_rates, dtype=float).flatten()
        n = len(rates)

        # Convert rates to means (1/rate)
        self._means = 1.0 / rates

        if probs is None:
            # Default: all continuation probabilities = 1 (except last)
            self._probs = np.ones(n - 1) if n > 1 else np.array([])
        else:
            probs_arr = np.array(probs, dtype=float).flatten()
            # If probs has n elements, drop the last one (implicit 1.0 for absorption)
            if len(probs_arr) == n:
                self._probs = probs_arr[:-1]
            else:
                self._probs = probs_arr

        if len(self._probs) != n - 1:
            raise ValueError(f"probs must have length {n - 1} or {n}, got {len(self._probs) + (n - (n-1))}")

        if np.any(self._probs < 0) or np.any(self._probs > 1):
            raise ValueError("All continuation probabilities must be in [0, 1]")

        if np.any(self._means <= 0):
            raise ValueError("All means must be positive")

        # Construct alpha and T
        alpha, T = self._build_representation()

        # Call parent constructor
        ContinuousDistribution.__init__(self)
        Markovian.__init__(self)
        self._name = 'Coxian'
        self._alpha = alpha
        self._T = T
        self._t = -self._T.sum(axis=1)

    def _build_representation(self) -> Tuple[np.ndarray, np.ndarray]:
        """Build the phase-type representation."""
        n = len(self._means)
        rates = 1.0 / self._means

        # Initial probability: start in phase 1
        alpha = np.zeros(n)
        alpha[0] = 1.0

        # Sub-generator matrix
        T = np.zeros((n, n))
        for i in range(n):
            T[i, i] = -rates[i]
            if i < n - 1:
                T[i, i + 1] = rates[i] * self._probs[i]

        return alpha, T

    @property
    def means(self) -> np.ndarray:
        """Get the phase means."""
        return self._means.copy()

    @property
    def probs(self) -> np.ndarray:
        """Get the continuation probabilities."""
        return self._probs.copy()



class MAP(ContinuousDistribution, Markovian):
    """
    Markovian Arrival Process.

    A MAP is a generalization of the Poisson process that can capture
    correlation between inter-arrival times. It is defined by two
    matrices D0 and D1 where:
    - D0 contains transition rates without arrivals
    - D1 contains transition rates with arrivals
    - D0 + D1 is a valid generator matrix

    Args:
        D0: Matrix of transition rates without arrivals.
        D1: Matrix of transition rates with arrivals.
    """

    def __init__(self, D0: Union[list, np.ndarray], D1: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'MAP'
        self._D0 = np.atleast_2d(np.array(D0, dtype=float))
        self._D1 = np.atleast_2d(np.array(D1, dtype=float))

        n = self._D0.shape[0]
        if self._D0.shape != (n, n) or self._D1.shape != (n, n):
            raise ValueError("D0 and D1 must be square matrices of the same size")

        # Verify D0 + D1 is a valid generator
        Q = self._D0 + self._D1
        if not np.allclose(Q.sum(axis=1), 0):
            raise ValueError("D0 + D1 must have zero row sums (generator matrix)")

        # Compute stationary distribution
        self._pi = self._compute_stationary()

    def _compute_stationary(self) -> np.ndarray:
        """Compute the stationary distribution of the underlying CTMC."""
        Q = self._D0 + self._D1
        n = Q.shape[0]

        # Solve pi * Q = 0 with sum(pi) = 1
        A = np.vstack([Q.T, np.ones(n)])
        b = np.zeros(n + 1)
        b[-1] = 1.0

        try:
            pi, _, _, _ = linalg.lstsq(A, b)
            pi = np.maximum(pi, 0)
            pi /= pi.sum()
            return pi
        except linalg.LinAlgError:
            return np.ones(n) / n

    @property
    def D0(self) -> np.ndarray:
        """Get the D0 matrix."""
        return self._D0.copy()

    @property
    def D1(self) -> np.ndarray:
        """Get the D1 matrix."""
        return self._D1.copy()

    @property
    def pi(self) -> np.ndarray:
        """Get the stationary distribution."""
        return self._pi.copy()

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix."""
        return self._D0.copy()

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix."""
        return self._D1.copy()

    def getMean(self) -> float:
        """Get the mean inter-arrival time."""
        # Mean = pi * (-D0)^(-1) * e
        try:
            D0_inv = linalg.inv(-self._D0)
            e = np.ones(len(self._pi))
            return float(self._pi @ D0_inv @ e)
        except linalg.LinAlgError:
            return float('nan')

    def getVar(self) -> float:
        """Get the variance of inter-arrival times."""
        try:
            D0_inv = linalg.inv(-self._D0)
            e = np.ones(len(self._pi))
            mean = float(self._pi @ D0_inv @ e)
            second_moment = float(2 * self._pi @ D0_inv @ D0_inv @ e)
            return second_moment - mean ** 2
        except linalg.LinAlgError:
            return float('nan')

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return self._D0.shape[0]

    def getMu(self) -> np.ndarray:
        """Get the arrival rates from each phase."""
        return self._D1.sum(axis=1)

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities (all 1 for MAP arrivals)."""
        return np.ones(self.getNumberOfPhases())

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector (stationary distribution)."""
        return self._pi.copy()

    def getRate(self) -> float:
        """Get the arrival rate (1/mean)."""
        mean = self.getMean()
        if mean == 0 or np.isnan(mean):
            return float('inf')
        return 1.0 / mean

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random inter-arrival times."""
        if rng is None:
            rng = np.random.default_rng()

        samples = np.zeros(n)
        num_phases = len(self._pi)

        # Start from stationary distribution
        phase = rng.choice(num_phases, p=self._pi)

        for i in range(n):
            time = 0.0

            while True:
                # Total rate out of current phase
                rate_out = -self._D0[phase, phase]
                time += rng.exponential(1.0 / rate_out)

                # Transition probabilities
                probs = np.zeros(2 * num_phases)
                probs[:num_phases] = np.maximum(self._D0[phase, :], 0)
                probs[phase] = 0  # No self-loops
                probs[num_phases:] = np.maximum(self._D1[phase, :], 0)

                total = probs.sum()
                if total <= 0:
                    break
                probs /= total

                # Next state
                next_state = rng.choice(2 * num_phases, p=probs)

                if next_state >= num_phases:
                    # Arrival occurred, record time
                    samples[i] = time
                    phase = next_state - num_phases
                    break
                else:
                    # No arrival, continue
                    phase = next_state

        return samples



class MMPP(MAP):
    """
    Markov-Modulated Poisson Process.

    A special case of MAP where arrivals occur according to Poisson
    processes with rates that depend on an underlying Markov chain state.

    Args:
        Q: Generator matrix of the modulating Markov chain.
        lambda_: Vector of Poisson rates for each state.
    """

    def __init__(self, Q: Union[list, np.ndarray], lambda_: Union[list, np.ndarray]):
        self._Q = np.atleast_2d(np.array(Q, dtype=float))
        self._lambda = np.atleast_1d(np.array(lambda_, dtype=float))

        n = self._Q.shape[0]
        if self._Q.shape != (n, n):
            raise ValueError("Q must be a square matrix")
        if len(self._lambda) != n:
            raise ValueError(f"lambda must have length {n}")

        # Construct D0 and D1 from Q and lambda
        D0 = self._Q - np.diag(self._lambda)
        D1 = np.diag(self._lambda)

        super().__init__(D0, D1)
        self._name = 'MMPP'

    @property
    def Q(self) -> np.ndarray:
        """Get the modulating Markov chain generator."""
        return self._Q.copy()

    @property
    def lambda_(self) -> np.ndarray:
        """Get the Poisson rates for each state."""
        return self._lambda.copy()



class Cox2(Coxian):
    """
    2-phase Coxian distribution.

    Convenience class for the common 2-phase case.

    Args:
        mean1: Mean time in phase 1.
        mean2: Mean time in phase 2.
        p: Probability of going to phase 2 (vs absorbing after phase 1).
    """

    def __init__(self, mean1: float, mean2: float, p: float):
        if p < 0 or p > 1:
            raise ValueError("Continuation probability p must be in [0, 1]")
        super().__init__([mean1, mean2], [p])
        self._name = 'Cox2'


class MMPP2(MAP):
    """
    Markov-Modulated Poisson Process with 2 states.

    A special case of MAP with two states, where arrivals occur according
    to Poisson processes with rates lambda0 and lambda1, and the modulating
    chain switches between states with rates sigma0 and sigma1.

    Args:
        lambda0: Arrival rate in state 0.
        lambda1: Arrival rate in state 1.
        sigma0: Transition rate from state 0 to state 1.
        sigma1: Transition rate from state 1 to state 0.
    """

    def __init__(self, lambda0: float, lambda1: float, sigma0: float, sigma1: float):
        self._lambda0 = lambda0
        self._lambda1 = lambda1
        self._sigma0 = sigma0
        self._sigma1 = sigma1

        # Build D0 and D1 matrices
        D0 = np.array([
            [-lambda0 - sigma0, sigma0],
            [sigma1, -lambda1 - sigma1]
        ])
        D1 = np.array([
            [lambda0, 0],
            [0, lambda1]
        ])

        super().__init__(D0, D1)
        self._name = 'MMPP2'

    @property
    def lambda0(self) -> float:
        """Get arrival rate in state 0."""
        return self._lambda0

    @property
    def lambda1(self) -> float:
        """Get arrival rate in state 1."""
        return self._lambda1

    @property
    def sigma0(self) -> float:
        """Get transition rate from state 0 to 1."""
        return self._sigma0

    @property
    def sigma1(self) -> float:
        """Get transition rate from state 1 to 0."""
        return self._sigma1

    def getIDC(self) -> float:
        """Get the Index of Dispersion for Counts."""
        numerator = 2 * (self._lambda0 - self._lambda1) ** 2 * self._sigma0 * self._sigma1
        denominator = ((self._sigma0 + self._sigma1) ** 2 *
                       (self._lambda0 * self._sigma1 + self._lambda1 * self._sigma0))
        return 1 + numerator / denominator if denominator > 0 else 1.0

    @staticmethod
    def fit_mean_scv_acf(mean: float, scv: float, acf_decay: float) -> 'MMPP2':
        """
        Fit MMPP2 to match mean, SCV, and ACF decay.

        Args:
            mean: Target mean inter-arrival time.
            scv: Target squared coefficient of variation.
            acf_decay: ACF decay parameter (gamma2).

        Returns:
            Fitted MMPP2 distribution.
        """
        # Simplified fitting based on moments
        rate = 1.0 / mean
        if scv <= 1.0:
            # Low variability: use similar rates
            lambda0 = rate * (1 + 0.1)
            lambda1 = rate * (1 - 0.1)
        else:
            # High variability: use different rates
            diff = np.sqrt(scv - 1) * rate
            lambda0 = rate + diff
            lambda1 = max(0.01 * rate, rate - diff)

        # Switching rates based on ACF decay
        sigma = rate * (1 - acf_decay) / 2 if acf_decay < 1 else rate
        return MMPP2(lambda0, lambda1, sigma, sigma)



class ME(ContinuousDistribution, Markovian):
    """
    Matrix Exponential distribution.

    ME distributions generalize Phase-Type distributions by allowing
    the initial vector alpha to have entries outside [0,1] and the
    matrix A to have arbitrary structure (not necessarily a sub-generator).

    Args:
        alpha: Initial vector (may have negative entries).
        A: Matrix parameter (eigenvalues must have negative real parts).
    """

    def __init__(self, alpha: Union[list, np.ndarray], A: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'ME'
        self._alpha = np.atleast_1d(np.array(alpha, dtype=float)).flatten()
        self._A = np.atleast_2d(np.array(A, dtype=float))

        n = len(self._alpha)
        if self._A.shape != (n, n):
            raise ValueError(f"A must be {n}x{n} to match alpha of length {n}")

        # Verify eigenvalues have negative real parts
        eigenvalues = linalg.eigvals(self._A)
        if not all(np.real(eigenvalues) < 0):
            raise ValueError("All eigenvalues of A must have negative real parts")

    @property
    def alpha(self) -> np.ndarray:
        """Get the initial vector."""
        return self._alpha.copy()

    @property
    def A(self) -> np.ndarray:
        """Get the matrix parameter."""
        return self._A.copy()

    def getAlpha(self) -> np.ndarray:
        """Get the initial vector (getter method for API compatibility)."""
        return self._alpha.copy()

    def getA(self) -> np.ndarray:
        """Get the matrix parameter (getter method for API compatibility)."""
        return self._A.copy()

    def getMean(self) -> float:
        """Get the mean."""
        try:
            A_inv = linalg.inv(self._A)
            e = np.ones(len(self._alpha))
            return float(-self._alpha @ A_inv @ e)
        except linalg.LinAlgError:
            return float('nan')

    def getVar(self) -> float:
        """Get the variance."""
        try:
            A_inv = linalg.inv(self._A)
            e = np.ones(len(self._alpha))
            mean = float(-self._alpha @ A_inv @ e)
            second_moment = float(2 * self._alpha @ A_inv @ A_inv @ e)
            return second_moment - mean ** 2
        except linalg.LinAlgError:
            return float('nan')

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return len(self._alpha)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        e = np.ones(len(self._alpha))
        expm = linalg.expm(self._A * x)
        return 1.0 - float(self._alpha @ expm @ e)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        expm = linalg.expm(self._A * x)
        t = -self._A.sum(axis=1)
        return float(self._alpha @ expm @ t)

    @staticmethod
    def from_exp(rate: float) -> 'ME':
        """Create ME from exponential distribution."""
        return ME(np.array([1.0]), np.array([[-rate]]))

    @staticmethod
    def from_erlang(k: int, rate: float) -> 'ME':
        """Create ME from Erlang distribution."""
        alpha = np.zeros(k)
        alpha[0] = 1.0
        A = np.diag([-rate] * k) + np.diag([rate] * (k - 1), 1)
        return ME(alpha, A)

    @staticmethod
    def from_hyper_exp(p: Union[float, list, np.ndarray],
                       rates: Union[float, list, np.ndarray],
                       rate2: Optional[float] = None) -> 'ME':
        """Create ME from HyperExp distribution.

        Can be called as:
        - from_hyper_exp(p, rate1, rate2) for 2-phase
        - from_hyper_exp([p1, p2, ...], [r1, r2, ...]) for n-phase
        """
        if rate2 is not None:
            # 2-phase form: from_hyper_exp(p, rate1, rate2)
            alpha = np.array([float(p), 1.0 - float(p)])
            A = np.array([[-float(rates), 0.0], [0.0, -float(rate2)]])
        else:
            # n-phase form: from_hyper_exp(probs, rates)
            probs = np.atleast_1d(np.array(p, dtype=float))
            rate_arr = np.atleast_1d(np.array(rates, dtype=float))
            if len(probs) != len(rate_arr):
                raise ValueError("probs and rates must have the same length")
            alpha = probs
            A = np.diag(-rate_arr)
        return ME(alpha, A)

    @staticmethod
    def fit_moments(moments: Union[list, np.ndarray]) -> 'ME':
        """Create ME distribution matching given moments.

        Uses the van de Liefvoort algorithm to construct a matrix-exponential
        distribution that has the specified moments.

        Args:
            moments: List of moments. To obtain an ME of order M, provide
                    2*M-1 moments (e.g., 3 moments for order 2, 5 for order 3).

        Returns:
            ME distribution matching the given moments.

        Raises:
            ValueError: If moments cannot be matched (e.g., invalid structure).
        """
        import numpy.matlib as ml

        def _reduced_moms_from_moms(m):
            """Convert raw moments to reduced moments (m_i / i!)."""
            rm = m[:]
            f = 1.0
            for i in range(len(m)):
                f = f / (i + 1)
                rm[i] *= f
            return rm

        def _appie_algorithm(rmom):
            """Internal implementation of van de Liefvoort's algorithm."""
            m = len(rmom)
            if m % 2 == 0:
                rm = rmom[0:m-1]
                m = int(m / 2)
            else:
                rm = rmom
                m = int((m+1) / 2)

            rm = np.insert(np.array(rm), 0, 1.0)
            f = np.zeros((2*m, 1))
            f[0] = 1.0
            y = np.zeros((2*m, 1))
            n = 0
            k = 0
            q = 1
            d = [0] * m
            alpha_mat = np.zeros((m, m))
            beta = np.zeros((m, 1))

            def shift(arr):
                sh = np.roll(arr, 1)
                sh[0] = 0
                return sh

            for i in range(2*m):
                ro = q * np.dot(rm, f)
                nold = n
                n = nold + 1
                yold = y
                if n > 0 and ro != 0:
                    if k > 0:
                        beta[k-1] = ro / np.power(rm[1], d[k-1]+n-1)
                    k = k + 1
                    d[k-1] = n
                    n = -n
                    q = q / ro
                    y = shift(f)
                elif n <= 0:
                    j = nold + d[k-1]
                    alpha_mat[k-1, j] = ro / rm[1]**j
                f = shift(f) - ro * yold

            if sum(d) != m:
                raise Exception("Insufficient matrix order!")

            K = ml.zeros((m, m))
            K[0, 0] = rm[1]
            for i in range(m-1):
                K[i, i+1] = rm[1]
            ind = d[0]
            for i in range(1, m):
                if ind < m:
                    inc = d[i]
                    ind = ind + inc
                    if ind <= m:
                        K[ind-1, ind-inc-d[i-1]] = beta[i-1]
                        for j in range(1, inc+1):
                            K[ind-1, ind-j] = alpha_mat[i, j-1]
            return K

        try:
            import math
            moms = list(moments)
            rmoms = _reduced_moms_from_moms(moms)
            K = _appie_algorithm(rmoms)
            N = int(math.ceil(len(moms) / 2.0))

            T = ml.zeros((N, N))
            for i in range(N):
                for j in range(i+1):
                    T[i, j] = 1.0
            U = ml.zeros((N, N))
            for i in range(N):
                for j in range(i, N):
                    U[i, j] = 1.0 / (N-i)

            alpha = ml.zeros((1, N))
            alpha[0, 0] = 1.0
            alpha = alpha * np.linalg.inv(T) * U
            A = np.linalg.inv(U) * T * K * np.linalg.inv(T) * U
            A = np.linalg.inv(-A)

            alpha_arr = np.asarray(alpha).flatten()
            A_arr = np.asarray(A)
            return ME(alpha_arr, A_arr)
        except Exception as e:
            raise ValueError(f"Could not fit ME to moments: {e}")

    # CamelCase aliases for API compatibility
    fromExp = from_exp
    fromErlang = from_erlang
    fromHyperExp = from_hyper_exp
    fitMoments = fit_moments

    # Snake_case aliases for getters
    def get_alpha(self) -> np.ndarray:
        """Get the initial vector."""
        return self._alpha.copy()

    def get_a(self) -> np.ndarray:
        """Get the matrix parameter."""
        return self._A.copy()

    def get_mean(self) -> float:
        """Get the mean."""
        return self.getMean()

    def get_var(self) -> float:
        """Get the variance."""
        return self.getVar()

    def get_number_of_phases(self) -> int:
        """Get the number of phases."""
        return self.getNumberOfPhases()

    def eval_cdf(self, x: float) -> float:
        """Evaluate CDF at point x."""
        return self.evalCDF(x)

    def eval_pdf(self, x: float) -> float:
        """Evaluate PDF at point x."""
        return self.evalPDF(x)



class MarkedMAP(ContinuousDistribution, Markovian):
    """
    Marked Markovian Arrival Process.

    A MarkedMAP extends MAP to support multiple arrival types (marks/classes).
    It is defined by matrices D0, D1, D2, ..., Dk where:
    - D0: transitions without arrivals
    - Dk (k >= 1): transitions with arrivals of type k

    Args:
        process: List or array of matrices [D0, D1, D2, ..., Dk].
    """

    def __init__(self, process: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'MarkedMAP'

        # Convert to list of numpy arrays
        self._process = [np.atleast_2d(np.array(m, dtype=float)) for m in process]

        if len(self._process) < 2:
            raise ValueError("MarkedMAP requires at least D0 and D1")

        n = self._process[0].shape[0]
        for i, m in enumerate(self._process):
            if m.shape != (n, n):
                raise ValueError(f"All matrices must be {n}x{n}")

        # Verify D0 + sum(Dk) is a valid generator
        generator = sum(self._process)
        if not np.allclose(generator.sum(axis=1), 0, atol=1e-6):
            raise ValueError("D0 + sum(Dk) must have zero row sums")

        # Compute stationary distribution
        self._pi = self._compute_stationary()

    def _compute_stationary(self) -> np.ndarray:
        """Compute the stationary distribution."""
        Q = sum(self._process)
        n = Q.shape[0]
        A = np.vstack([Q.T, np.ones(n)])
        b = np.zeros(n + 1)
        b[-1] = 1.0
        try:
            pi, _, _, _ = linalg.lstsq(A, b)
            pi = np.maximum(pi, 0)
            pi /= pi.sum()
            return pi
        except linalg.LinAlgError:
            return np.ones(n) / n

    def D(self, k: int) -> np.ndarray:
        """Get the k-th matrix (D0, D1, ..., Dk)."""
        if k < 0 or k >= len(self._process):
            raise IndexError(f"Matrix index {k} out of range")
        return self._process[k].copy()

    @property
    def num_types(self) -> int:
        """Get the number of arrival types."""
        return len(self._process) - 1

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return self._process[0].shape[0]

    def getMean(self) -> float:
        """Get the mean inter-arrival time."""
        try:
            D0_inv = linalg.inv(-self._process[0])
            e = np.ones(len(self._pi))
            return float(self._pi @ D0_inv @ e)
        except linalg.LinAlgError:
            return float('nan')

    def getVar(self) -> float:
        """Get the variance."""
        try:
            D0_inv = linalg.inv(-self._process[0])
            e = np.ones(len(self._pi))
            mean = float(self._pi @ D0_inv @ e)
            second_moment = float(2 * self._pi @ D0_inv @ D0_inv @ e)
            return second_moment - mean ** 2
        except linalg.LinAlgError:
            return float('nan')

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None
               ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate random inter-arrival times and their types.

        Returns:
            Tuple of (times, types) arrays.
        """
        if rng is None:
            rng = np.random.default_rng()

        times = np.zeros(n)
        types = np.zeros(n, dtype=int)
        num_phases = len(self._pi)

        # Aggregate D1 for total arrival rate per transition
        D1_total = sum(self._process[1:])

        phase = rng.choice(num_phases, p=self._pi)

        for i in range(n):
            time = 0.0

            while True:
                rate_out = -self._process[0][phase, phase]
                time += rng.exponential(1.0 / rate_out)

                # Build transition probabilities
                probs = []
                # D0 transitions (no arrival)
                for j in range(num_phases):
                    if j != phase:
                        probs.append(max(0, self._process[0][phase, j]))
                    else:
                        probs.append(0)
                # Arrival transitions (per type)
                for k in range(1, len(self._process)):
                    for j in range(num_phases):
                        probs.append(max(0, self._process[k][phase, j]))

                probs = np.array(probs)
                total = probs.sum()
                if total <= 0:
                    break
                probs /= total

                next_state = rng.choice(len(probs), p=probs)

                if next_state < num_phases:
                    # D0 transition, no arrival
                    phase = next_state
                else:
                    # Arrival occurred
                    idx = next_state - num_phases
                    arrival_type = idx // num_phases + 1
                    new_phase = idx % num_phases
                    times[i] = time
                    types[i] = arrival_type
                    phase = new_phase
                    break

        return times, types



class BMAP(MarkedMAP):
    """
    Batch Markovian Arrival Process.

    BMAP is a point process where arrivals occur in batches.
    Each matrix Dk represents transitions generating k arrivals.

    Args:
        process: List of matrices [D0, D1, D2, ..., Dk] where
                 Dk is the rate matrix for batch size k.
    """

    def __init__(self, process: Union[list, np.ndarray]):
        super().__init__(process)
        self._name = 'BMAP'

    @property
    def max_batch_size(self) -> int:
        """Get the maximum batch size."""
        return len(self._process) - 1

    def getMeanBatchSize(self) -> float:
        """Get the mean batch size."""
        total_rate = 0.0
        weighted_sum = 0.0

        for k in range(1, len(self._process)):
            rate_k = float(self._pi @ self._process[k] @ np.ones(len(self._pi)))
            total_rate += rate_k
            weighted_sum += k * rate_k

        return weighted_sum / total_rate if total_rate > 0 else 0.0

    def getBatchRates(self) -> np.ndarray:
        """Get arrival rates for each batch size."""
        rates = np.zeros(self.max_batch_size)
        for k in range(1, len(self._process)):
            rates[k - 1] = float(self._pi @ self._process[k] @ np.ones(len(self._pi)))
        return rates

    @staticmethod
    def from_map_with_batch_pmf(D0: np.ndarray, D1: np.ndarray,
                                 batch_sizes: list, pmf: list) -> 'BMAP':
        """
        Create BMAP from base MAP and batch size distribution.

        Args:
            D0: Base MAP D0 matrix.
            D1: Base MAP D1 matrix.
            batch_sizes: List of possible batch sizes.
            pmf: Probability mass function for batch sizes.

        Returns:
            BMAP with specified batch distribution.
        """
        D0 = np.atleast_2d(np.array(D0, dtype=float))
        D1 = np.atleast_2d(np.array(D1, dtype=float))

        # Normalize PMF
        pmf = np.array(pmf) / sum(pmf)

        max_batch = max(batch_sizes)
        Dk = [np.zeros_like(D0) for _ in range(max_batch)]

        for size, prob in zip(batch_sizes, pmf):
            Dk[size - 1] += D1 * prob

        return BMAP([D0] + Dk)



class RAP(ContinuousDistribution, Markovian):
    """
    Rational Arrival Process.

    RAP generalizes MAP by allowing the representation matrices
    to have entries that may lead to complex eigenvalues, as long as
    the resulting distribution is valid (non-negative PDF and CDF in [0,1]).

    Similar to MAP, RAP is defined by two matrices H0 and H1 where:
    - H0: Rate transitions without arrivals (like D0 in MAP)
    - H1: Rate transitions with arrivals (like D1 in MAP)

    Args:
        H0: Square matrix for transitions without arrivals.
        H1: Square matrix for transitions with arrivals.
    """

    def __init__(self, H0: Union[list, np.ndarray], H1: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'RAP'
        self._H0 = np.atleast_2d(np.array(H0, dtype=float))
        self._H1 = np.atleast_2d(np.array(H1, dtype=float))

        n = self._H0.shape[0]
        if self._H0.shape != (n, n) or self._H1.shape != (n, n):
            raise ValueError("H0 and H1 must be square matrices of the same size")

        # Compute stationary distribution
        self._pi = self._compute_stationary()

    def _compute_stationary(self) -> np.ndarray:
        """Compute the stationary distribution of the underlying CTMC."""
        Q = self._H0 + self._H1
        n = Q.shape[0]

        # Solve pi * Q = 0 with sum(pi) = 1
        A = np.vstack([Q.T, np.ones(n)])
        b = np.zeros(n + 1)
        b[-1] = 1.0

        try:
            pi, _, _, _ = linalg.lstsq(A, b)
            pi = np.maximum(pi, 0)
            if pi.sum() > 0:
                pi /= pi.sum()
            else:
                pi = np.ones(n) / n
            return pi
        except linalg.LinAlgError:
            return np.ones(n) / n

    @property
    def H0(self) -> np.ndarray:
        """Get the H0 matrix."""
        return self._H0.copy()

    @property
    def H1(self) -> np.ndarray:
        """Get the H1 matrix."""
        return self._H1.copy()

    @property
    def pi(self) -> np.ndarray:
        """Get the stationary distribution."""
        return self._pi.copy()

    def getH0(self) -> np.ndarray:
        """Get the H0 matrix."""
        return self._H0.copy()

    def getH1(self) -> np.ndarray:
        """Get the H1 matrix."""
        return self._H1.copy()

    def getMean(self) -> float:
        """Get the mean inter-arrival time."""
        try:
            H0_inv = linalg.inv(-self._H0)
            n = self._H0.shape[0]
            e = np.ones(n)
            # Use embedded stationary distribution (at arrivals)
            # beta = pi * H1 / (pi * H1 * e)
            arrival_rate = self._pi @ self._H1 @ e
            if arrival_rate > 0:
                beta = (self._pi @ self._H1) / arrival_rate
            else:
                beta = self._pi
            return float(beta @ H0_inv @ e)
        except linalg.LinAlgError:
            return float('nan')

    def getVar(self) -> float:
        """Get the variance."""
        try:
            H0_inv = linalg.inv(-self._H0)
            n = self._H0.shape[0]
            e = np.ones(n)
            # Use embedded stationary distribution (at arrivals)
            arrival_rate = self._pi @ self._H1 @ e
            if arrival_rate > 0:
                beta = (self._pi @ self._H1) / arrival_rate
            else:
                beta = self._pi
            mean = float(beta @ H0_inv @ e)
            second_moment = float(2 * beta @ H0_inv @ H0_inv @ e)
            return second_moment - mean ** 2
        except linalg.LinAlgError:
            return float('nan')

    def getRate(self) -> float:
        """Get the arrival rate (1/mean)."""
        mean = self.getMean()
        if mean > 0:
            return 1.0 / mean
        return float('nan')

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return self._H0.shape[0]

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        n = self._H0.shape[0]
        e = np.ones(n)
        expm = linalg.expm(self._H0 * x)
        return 1.0 - float(self._pi @ expm @ e)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        expm = linalg.expm(self._H0 * x)
        h = -self._H0.sum(axis=1)
        return float(self._pi @ expm @ h)

    @staticmethod
    def from_exp(rate: float) -> 'RAP':
        """Create RAP from exponential distribution."""
        H0 = np.array([[-rate]])
        H1 = np.array([[rate]])
        return RAP(H0, H1)

    @staticmethod
    def from_erlang(k: int, rate: float) -> 'RAP':
        """Create RAP from Erlang distribution."""
        H0 = np.diag([-rate] * k) + np.diag([rate] * (k - 1), 1)
        H1 = np.zeros((k, k))
        H1[k - 1, 0] = rate
        return RAP(H0, H1)

    @staticmethod
    def from_map(map_dist: 'MAP') -> 'RAP':
        """Create RAP from MAP distribution."""
        return RAP(map_dist.D0, map_dist.D1)

    @staticmethod
    def from_poisson(rate: float) -> 'RAP':
        """Create RAP from Poisson process (same as exponential)."""
        return RAP.from_exp(rate)

    # CamelCase aliases for API compatibility
    fromExp = from_exp
    fromErlang = from_erlang
    fromMAP = from_map
    fromPoisson = from_poisson

    # Snake_case aliases for getters
    def get_h0(self) -> np.ndarray:
        """Get the H0 matrix."""
        return self._H0.copy()

    def get_h1(self) -> np.ndarray:
        """Get the H1 matrix."""
        return self._H1.copy()

    def get_mean(self) -> float:
        """Get the mean."""
        return self.getMean()

    def get_var(self) -> float:
        """Get the variance."""
        return self.getVar()

    def get_rate(self) -> float:
        """Get the arrival rate."""
        return self.getRate()

    def get_number_of_phases(self) -> int:
        """Get the number of phases."""
        return self.getNumberOfPhases()

    def eval_cdf(self, x: float) -> float:
        """Evaluate CDF at point x."""
        return self.evalCDF(x)

    def eval_pdf(self, x: float) -> float:
        """Evaluate PDF at point x."""
        return self.evalPDF(x)



class MMDP(ContinuousDistribution, Markovian):
    """
    Markov-Modulated Deterministic Process.

    MMDP models a process where deterministic inter-arrival times
    are modulated by an underlying Markov chain. In each state,
    arrivals occur at deterministic intervals.

    Args:
        Q: Generator matrix of the modulating Markov chain.
        d: Deterministic inter-arrival times for each state.
    """

    def __init__(self, Q: Union[list, np.ndarray], d: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'MMDP'
        self._Q = np.atleast_2d(np.array(Q, dtype=float))
        self._d = np.atleast_1d(np.array(d, dtype=float))

        n = self._Q.shape[0]
        if self._Q.shape != (n, n):
            raise ValueError("Q must be a square matrix")
        if len(self._d) != n:
            raise ValueError(f"d must have length {n}")

        # Compute stationary distribution
        A = np.vstack([self._Q.T, np.ones(n)])
        b = np.zeros(n + 1)
        b[-1] = 1.0
        try:
            pi, _, _, _ = linalg.lstsq(A, b)
            pi = np.maximum(pi, 0)
            self._pi = pi / pi.sum()
        except linalg.LinAlgError:
            self._pi = np.ones(n) / n

    @property
    def Q(self) -> np.ndarray:
        """Get the generator matrix."""
        return self._Q.copy()

    @property
    def d(self) -> np.ndarray:
        """Get the deterministic times."""
        return self._d.copy()

    def getMean(self) -> float:
        """Get the mean inter-arrival time."""
        return float(np.dot(self._pi, self._d))

    def getVar(self) -> float:
        """Get the variance."""
        mean = self.getMean()
        return float(np.dot(self._pi, self._d ** 2)) - mean ** 2

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return len(self._d)



class MMDP2(MMDP):
    """
    Markov-Modulated Deterministic Process with 2 states.

    Convenience class for the common 2-state case.

    Args:
        d0: Deterministic time in state 0.
        d1: Deterministic time in state 1.
        sigma0: Transition rate from state 0 to state 1.
        sigma1: Transition rate from state 1 to state 0.
    """

    def __init__(self, d0: float, d1: float, sigma0: float, sigma1: float):
        Q = np.array([
            [-sigma0, sigma0],
            [sigma1, -sigma1]
        ])
        d = np.array([d0, d1])
        super().__init__(Q, d)
        self._name = 'MMDP2'

        self._d0 = d0
        self._d1 = d1
        self._sigma0 = sigma0
        self._sigma1 = sigma1

    @property
    def d0(self) -> float:
        """Get deterministic time in state 0."""
        return self._d0

    @property
    def d1(self) -> float:
        """Get deterministic time in state 1."""
        return self._d1

    @property
    def sigma0(self) -> float:
        """Get transition rate from state 0 to 1."""
        return self._sigma0

    @property
    def sigma1(self) -> float:
        """Get transition rate from state 1 to 0."""
        return self._sigma1



class MarkedMMPP(ContinuousDistribution, Markovian):
    """
    Marked Markov-Modulated Poisson Process (M3PP).

    A MarkedMMPP extends MMPP to support multiple arrival types (marks).
    In each state, arrivals of different types occur according to Poisson
    processes. The D1k matrices must be diagonal (MMPP constraint).

    Uses the M3A representation format: D = {D0, D1, D11, D12, ..., D1K}
    where K is the number of marking types.

    Args:
        process: List of matrices [D0, D1, D11, D12, ..., D1K] or
                 [D0, D11, D12, ..., D1K] (D1 computed as sum).
        num_types: Number of marking types K.
    """

    def __init__(self, process: Union[list, np.ndarray], num_types: int):
        super().__init__()
        self._name = 'MarkedMMPP'

        # Convert to list of numpy arrays
        self._process = [np.atleast_2d(np.array(m, dtype=float)) for m in process]
        self._num_types = num_types

        # Validate format
        if len(self._process) < 2:
            raise ValueError("MarkedMMPP requires at least D0 and one D1k matrix")

        n = self._process[0].shape[0]
        for i, m in enumerate(self._process):
            if m.shape != (n, n):
                raise ValueError(f"All matrices must be {n}x{n}")

        # Handle different input formats
        if num_types == len(self._process) - 1:
            # Format: [D0, D11, D12, ..., D1K] - compute D1
            D0 = self._process[0]
            D1k_list = self._process[1:]
            D1 = sum(D1k_list)
            self._D0 = D0
            self._D1 = D1
            self._D1k = D1k_list
        elif num_types == len(self._process) - 2:
            # Format: [D0, D1, D11, D12, ..., D1K]
            self._D0 = self._process[0]
            self._D1 = self._process[1]
            self._D1k = self._process[2:]
        else:
            raise ValueError("Inconsistency between num_types and number of matrices")

        # Verify D1 and D1k matrices are diagonal (MMPP constraint)
        if not np.allclose(self._D1 - np.diag(np.diag(self._D1)), 0, atol=1e-10):
            raise ValueError("D1 must be diagonal for MarkedMMPP")
        for k, D1k in enumerate(self._D1k):
            if not np.allclose(D1k - np.diag(np.diag(D1k)), 0, atol=1e-10):
                raise ValueError(f"D1{k+1} must be diagonal for MarkedMMPP")

        # Compute stationary distribution
        Q = self._D0 + self._D1
        A = np.vstack([Q.T, np.ones(n)])
        b = np.zeros(n + 1)
        b[-1] = 1.0
        try:
            pi, _, _, _ = linalg.lstsq(A, b)
            pi = np.maximum(pi, 0)
            self._pi = pi / pi.sum()
        except linalg.LinAlgError:
            self._pi = np.ones(n) / n

    def D(self, i: int, j: int = 0) -> np.ndarray:
        """
        Get representation matrix.

        Args:
            i: Primary index (0 for D0, 1 for D1)
            j: Secondary index for D1j (0 for aggregate D1, k for D1k)

        Returns:
            The requested matrix.
        """
        if i == 0:
            return self._D0.copy()
        elif i == 1:
            if j == 0:
                return self._D1.copy()
            elif 1 <= j <= self._num_types:
                return self._D1k[j - 1].copy()
            else:
                raise IndexError(f"D1{j} out of range")
        else:
            raise IndexError(f"D{i} not valid")

    @property
    def num_types(self) -> int:
        """Get the number of marking types."""
        return self._num_types

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return self._D0.shape[0]

    def getNumberOfTypes(self) -> int:
        """Get the number of marking types."""
        return self._num_types

    def getMean(self) -> float:
        """Get the mean inter-arrival time."""
        try:
            D0_inv = linalg.inv(-self._D0)
            e = np.ones(len(self._pi))
            return float(self._pi @ D0_inv @ e)
        except linalg.LinAlgError:
            return float('nan')

    def getVar(self) -> float:
        """Get the variance."""
        try:
            D0_inv = linalg.inv(-self._D0)
            e = np.ones(len(self._pi))
            mean = float(self._pi @ D0_inv @ e)
            second_moment = float(2 * self._pi @ D0_inv @ D0_inv @ e)
            return second_moment - mean ** 2
        except linalg.LinAlgError:
            return float('nan')

    def getMarkedMeans(self) -> np.ndarray:
        """Get mean inter-arrival times for each marking type."""
        means = np.zeros(self._num_types)
        for k in range(self._num_types):
            rate_k = float(self._pi @ np.diag(self._D1k[k]))
            means[k] = 1.0 / rate_k if rate_k > 0 else float('inf')
        return means

    def getRate(self) -> float:
        """Get the aggregate arrival rate."""
        return float(self._pi @ np.diag(self._D1))

    def toMAP(self) -> MAP:
        """Convert to aggregate MAP (ignoring marks)."""
        return MAP(self._D0, self._D1)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None
               ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate random inter-arrival times and their types.

        Returns:
            Tuple of (times, types) arrays.
        """
        if rng is None:
            rng = np.random.default_rng()

        times = np.zeros(n)
        types = np.zeros(n, dtype=int)
        num_phases = len(self._pi)

        phase = rng.choice(num_phases, p=self._pi)

        for i in range(n):
            time = 0.0

            while True:
                rate_out = -self._D0[phase, phase]
                if rate_out <= 0:
                    break
                time += rng.exponential(1.0 / rate_out)

                # Build transition probabilities
                probs = []
                # D0 transitions (no arrival)
                for j in range(num_phases):
                    if j != phase:
                        probs.append(max(0, self._D0[phase, j]))
                    else:
                        probs.append(0)
                # Arrival transitions (per type) - diagonal entries only
                for k in range(self._num_types):
                    probs.append(max(0, self._D1k[k][phase, phase]))

                probs = np.array(probs)
                total = probs.sum()
                if total <= 0:
                    break
                probs /= total

                next_state = rng.choice(len(probs), p=probs)

                if next_state < num_phases:
                    # D0 transition, no arrival
                    phase = next_state
                else:
                    # Arrival occurred
                    arrival_type = next_state - num_phases + 1
                    times[i] = time
                    types[i] = arrival_type
                    # Phase doesn't change on diagonal arrival
                    break

        return times, types

    @staticmethod
    def rand(order: int = 2, num_classes: int = 2) -> 'MarkedMMPP':
        """Generate a random MarkedMMPP."""
        # Random D0 (off-diagonal transitions)
        D0 = np.random.rand(order, order)
        np.fill_diagonal(D0, 0)

        # Random diagonal D1k matrices
        D1k_list = []
        for _ in range(num_classes):
            D1k = np.diag(np.random.rand(order))
            D1k_list.append(D1k)

        # Set D0 diagonal for valid generator
        D1_sum = sum(D1k_list)
        for i in range(order):
            D0[i, i] = -(D0[i, :].sum() - D0[i, i] + D1_sum[i, i])

        return MarkedMMPP([D0] + D1k_list, num_classes)

