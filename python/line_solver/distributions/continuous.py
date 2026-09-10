"""
Continuous probability distributions for LINE (pure Python).

This module provides continuous distribution implementations including
exponential, deterministic, Erlang, hyperexponential, and other common
service time distributions.
"""

from typing import Optional, Tuple, Union, List
import numpy as np
from scipy import stats, linalg

from .base import ContinuousDistribution, Markovian


class Exp(ContinuousDistribution, Markovian):
    """
    Exponential distribution.

    The exponential distribution is the simplest continuous distribution
    for modeling service times in queueing systems. It has the memoryless
    property and SCV = 1.

    Args:
        rate: The rate parameter (lambda = 1/mean).
    """

    def __init__(self, rate: float):
        super().__init__()
        self._name = 'Exp'
        if rate <= 0:
            raise ValueError("Rate must be positive")
        self._rate = rate

    @classmethod
    def fit(cls, mean: float, scv: float = 1.0, skew: float = None) -> 'Exp':
        """
        Fit an exponential to the given moments (MATLAB Exp.fit).

        The exponential has SCV = 1 and skewness 2, so only the mean is used;
        the other moments are accepted for signature compatibility.
        """
        return cls.fit_mean_and_scv(mean, scv)

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float = 1.0) -> 'Exp':
        """
        Fit an exponential to a mean and SCV (MATLAB Exp.fitMeanAndSCV).

        An exponential cannot represent SCV != 1; MATLAB warns and uses SCV = 1,
        which is what happens here.
        """
        import warnings
        if abs(scv - 1.0) > 1e-3:
            warnings.warn('The exponential distribution cannot fit SCV != 1, '
                          'changing SCV to 1.', RuntimeWarning)
        return cls.fit_mean(mean)

    @classmethod
    def fitMeanAndSCV(cls, mean: float, scv: float = 1.0) -> 'Exp':
        """camelCase alias of fit_mean_and_scv (MATLAB/JAR spelling)."""
        return cls.fit_mean_and_scv(mean, scv)

    @classmethod
    def fit_mean(cls, mean: float) -> 'Exp':
        """
        Create an exponential distribution with the given mean.

        THE RATE IS CLAMPED to [GlobalConstants.Zero, GlobalConstants.Immediate],
        which is what MATLAB `Exp.fitMean` and the JAR twin both do
        (`min(Immediate, max(Zero, 1/MEAN))`) and what `fitRate` below already
        did here. Without it a mean BELOW FineTol (1e-8) built a different model
        in each codebase from the same script: `Exp.fit_mean(5e-10)` gave a rate
        of 2e9 where MATLAB gave 1e8, so `lqn_sockshop`, whose bookkeeping
        activities are written `Exp.fitMean(0.0000000005)`, round-tripped
        Python -> MATLAB into a model with 20x smaller demands at those
        activities. It disagreed only where the answer is near zero, which is
        exactly where a RELATIVE comparison is most severe: the JSON parity row
        reported maxrel=1 on QLen.

        Args:
            mean: Target mean.

        Returns:
            Exp distribution with the clamped rate.
        """
        if mean < 0:
            raise ValueError("Mean must be positive")
        from ..constants import GlobalConstants
        # A zero mean is the immediate activity, not a malformed one; MATLAB
        # takes 1/0 = Inf through the same min() and lands on Immediate.
        rate = GlobalConstants.Immediate if mean == 0 else 1.0 / mean
        return cls.fitRate(rate)

    # CamelCase alias
    fitMean = fit_mean

    @property
    def rate(self) -> float:
        """Get the rate parameter."""
        return self._rate

    @rate.setter
    def rate(self, value: float):
        """Set the rate parameter."""
        if value <= 0:
            raise ValueError("Rate must be positive")
        self._rate = value
        

    def getMean(self) -> float:
        """Get the mean (1/rate)."""
        return 1.0 / self._rate

    def getVar(self) -> float:
        """Get the variance (1/rate^2)."""
        return 1.0 / (self._rate ** 2)

    def getSCV(self) -> float:
        """Get the squared coefficient of variation (always 1 for exponential)."""
        return 1.0

    def getSkew(self) -> float:
        """Get the skewness (always 2 for exponential)."""
        return 2.0

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        return 1.0 - np.exp(-self._rate * x)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return self._rate * np.exp(-self._rate * x)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.exponential(scale=1.0/self._rate, size=n)

    def getNumberOfPhases(self) -> int:
        """Get the number of phases (1 for exponential)."""
        return 1

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix for MAP representation."""
        return np.array([[-self._rate]])

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix for MAP representation."""
        return np.array([[self._rate]])

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        return np.array([self._rate])

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        return np.array([1.0])

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        return np.array([1.0])


    @classmethod
    def fitRate(cls, rate: float) -> 'Exp':
        """
        Create an exponential distribution with the given rate.

        THE RATE IS CLAMPED TO [GlobalConstants.Zero, GlobalConstants.Immediate],
        exactly as MATLAB Exp.fitRate, the JAR twin and the C++ `exp_rate` do. A
        fitter is fed a COMPUTED rate -- an iterate of SolverLN, a refreshed
        arrival rate -- and a rate of zero is a starved element rather than a
        malformed model: raising here aborted the whole layered solve of
        lqn_ofbiz on an activity whose throughput was still zero. A rate the
        caller writes itself still goes through the constructor, which refuses a
        non-positive one.

        Args:
            rate: The rate parameter (lambda).

        Returns:
            Exp distribution with the clamped rate.
        """
        from ..constants import GlobalConstants
        return cls(rate=min(GlobalConstants.Immediate,
                            max(GlobalConstants.Zero, float(rate))))

    # Snake_case alias
    fit_rate = fitRate

    # Aliases
    get_rate = lambda self: self._rate


class Det(ContinuousDistribution):
    """
    Deterministic (constant) distribution.

    All service times are exactly equal to the specified value.
    Has SCV = 0 (no variability).

    Args:
        value: The constant service time value.
    """

    def __init__(self, value: float):
        super().__init__()
        self._name = 'Det'
        if value < 0:
            raise ValueError("Value must be non-negative")
        self._value = value

    @property
    def value(self) -> float:
        """Get the constant value."""
        return self._value

    @value.setter
    def value(self, val: float):
        """Set the constant value."""
        if val < 0:
            raise ValueError("Value must be non-negative")
        self._value = val
        

    def getMean(self) -> float:
        """Get the mean (equals the constant value)."""
        return self._value

    def getVar(self) -> float:
        """Get the variance (always 0)."""
        return 0.0

    def getSCV(self) -> float:
        """Get the SCV (always 0)."""
        return 0.0

    def getSkew(self) -> float:
        """Get the skewness (undefined, return 0)."""
        return 0.0

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        return 1.0 if x >= self._value else 0.0

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x (delta function, return inf at value)."""
        return float('inf') if x == self._value else 0.0

    def evalLST(self, s):
        """LST of a deterministic time: exp(-s*t). Matches MATLAB Det.evalLST.
        numpy rather than math, so a COMPLEX argument is admissible: transform
        inversion and root location both need one."""
        val = np.exp(-s * self._value)
        return complex(val) if isinstance(s, complex) else float(np.real(val))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples (all equal to value)."""
        return np.full(n, self._value)

    def isImmediate(self) -> bool:
        """Check if this is an immediate (zero) service."""
        return self._value == 0.0

    @classmethod
    def fit_mean(cls, mean: float) -> 'Det':
        """
        Create a deterministic distribution with the given mean.

        Since Det has zero variance, the mean equals the constant value.

        Args:
            mean: The mean (and constant value) of the distribution.

        Returns:
            Det distribution with the specified mean.
        """
        return cls(mean)

    # MATLAB-compatible alias
    fitMean = fit_mean


class Immediate(Det):
    """
    Immediate (zero delay) distribution.

    Represents instantaneous service with zero delay.
    """

    _instance = None

    def __init__(self):
        super().__init__(0.0)
        self._name = 'Immediate'

    @classmethod
    def getInstance(cls) -> 'Immediate':
        """Get singleton instance of Immediate distribution."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    # snake_case alias
    get_instance = getInstance

    def getSCV(self) -> float:
        """SCV of an immediate service, 1 as in MATLAB and the JAR. The variance
        over a zero mean is undefined, and the deterministic 0 that Det returns
        made an Immediate look like a Det service to any SCV-driven fit."""
        return 1.0

    def isImmediate(self) -> bool:
        """Check if this is immediate service."""
        return True


class Disabled(ContinuousDistribution):
    """
    Disabled distribution.

    Represents a disabled service (no service at all).
    Used for nodes that don't serve a particular job class.
    """

    _instance = None

    def __init__(self):
        super().__init__()
        self._name = 'Disabled'

    @classmethod
    def getInstance(cls) -> 'Disabled':
        """Get singleton instance of Disabled distribution."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    # snake_case alias
    get_instance = getInstance

    def getMean(self) -> float:
        """Get the mean (NaN: a disabled class has no service law at all).

        NaN and not infinity, matching MATLAB `Disabled.getMean` and the JAR's
        `Disabled.getMean`. The difference is load bearing wherever a caller
        selects the served classes with a `getMean() > tol` test: NaN fails that
        test, infinity passes it and admits every disabled class."""
        return float('nan')

    def getVar(self) -> float:
        """Get the variance (NaN, as in MATLAB and the JAR)."""
        return float('nan')

    def getSCV(self) -> float:
        """Get the SCV (NaN, as in MATLAB and the JAR)."""
        return float('nan')

    def getRate(self) -> float:
        """Get the rate (NaN, as in the JAR; not 1/inf = 0)."""
        return float('nan')

    def getSkew(self) -> float:
        """Get the skewness (NaN, as in the JAR)."""
        return float('nan')

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF (NaN, as in MATLAB and the JAR)."""
        return float('nan')

    def evalLST(self, s: float) -> float:
        """Evaluate the Laplace-Stieltjes transform (NaN, as in the JAR)."""
        return float('nan')

    def sample(self, n: int = 1, rng=None) -> np.ndarray:
        """Draw n samples, all NaN, as in MATLAB and the JAR."""
        return np.full(int(n), float('nan'))

    def isDisabled(self) -> bool:
        """Check if this distribution is disabled."""
        return True



class Erlang(ContinuousDistribution, Markovian):
    """
    Erlang distribution (sum of k exponentials).

    The Erlang distribution is the distribution of the sum of k
    independent exponential random variables with the same rate.
    It has SCV = 1/k.

    Args:
        phase_rate: Rate parameter for each exponential phase (alpha).
        nphases: Number of sequential exponential phases (r).

    Mean = nphases / phase_rate = r / alpha
    """

    def __init__(self, phase_rate: float, nphases: int):
        super().__init__()
        self._name = 'Erlang'
        if phase_rate <= 0:
            raise ValueError("Phase rate must be positive")
        if nphases < 1:
            raise ValueError("Number of phases must be at least 1")
        self._phase_rate = phase_rate
        self._phases = int(round(nphases))
        # Mean = nphases / phase_rate
        self._mean = self._phases / self._phase_rate

    @classmethod
    def fit(cls, mean: float, scv: float, skew: float = None) -> 'Erlang':
        """Fit an Erlang to the given moments (MATLAB Erlang.fit).

        The Erlang has one shape degree of freedom, so the skewness cannot be
        set independently and is ignored, as in MATLAB."""
        return cls.fit_mean_and_scv(mean, scv)

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Erlang':
        """
        Create an Erlang distribution from mean and SCV.

        For Erlang, SCV = 1/k where k is the number of phases, so the order is
        k = ceil(1/SCV): the achievable SCVs are 1, 1/2, 1/3, ... and the fit
        takes the first one AT OR BELOW the request. Rounding instead would
        return a different law -- at SCV=0.4, ceil gives 3 phases and round
        gives 2 -- and every solver downstream would answer a different model
        with no error raised. MATLAB `Erlang.fitMeanAndSCV`, the JAR and the
        C++ port all use ceil.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation, which must be <= 1.

        Returns:
            Erlang distribution with the given mean and the closest achievable
            SCV at or below the requested one.
        """
        if scv <= 0:
            raise ValueError("SCV must be positive")
        if mean <= 0:
            raise ValueError("Mean must be positive")
        if scv > 1:
            raise ValueError(
                "The Erlang distribution requires squared coefficient of variation <= 1")
        import math
        phases = int(math.ceil(1.0 / scv))
        phase_rate = phases / mean
        return cls(phase_rate=phase_rate, nphases=phases)

    @classmethod
    def fit_mean_and_order(cls, mean: float, phases: int) -> 'Erlang':
        """
        Create an Erlang distribution from mean and number of phases.

        Args:
            mean: Target mean.
            phases: Number of phases (order).

        Returns:
            Erlang distribution with given mean and phases.
        """
        if mean <= 0:
            raise ValueError("Mean must be positive")
        phase_rate = phases / mean
        return cls(phase_rate=phase_rate, nphases=phases)

    # CamelCase aliases
    fitMeanAndScv = fit_mean_and_scv
    fitMeanAndSCV = fit_mean_and_scv
    fitMeanAndOrder = fit_mean_and_order

    @property
    def phases(self) -> int:
        """Get the number of phases."""
        return self._phases

    def getMean(self) -> float:
        """Get the mean."""
        return self._mean

    def getVar(self) -> float:
        """Get the variance."""
        return self._mean ** 2 / self._phases

    def getSCV(self) -> float:
        """Get the SCV (1/phases)."""
        return 1.0 / self._phases

    def getSkew(self) -> float:
        """Get the skewness."""
        return 2.0 / np.sqrt(self._phases)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        return stats.gamma.cdf(x, a=self._phases, scale=1.0/self._phase_rate)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return stats.gamma.pdf(x, a=self._phases, scale=1.0/self._phase_rate)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.gamma(shape=self._phases, scale=1.0/self._phase_rate, size=n)

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return self._phases

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix for MAP representation."""
        k = self._phases
        D0 = np.zeros((k, k))
        for i in range(k):
            D0[i, i] = -self._phase_rate
            if i < k - 1:
                D0[i, i + 1] = self._phase_rate
        return D0

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix for MAP representation."""
        k = self._phases
        D1 = np.zeros((k, k))
        D1[k - 1, 0] = self._phase_rate
        return D1

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        return np.full(self._phases, self._phase_rate)

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        phi = np.zeros(self._phases)
        phi[-1] = 1.0
        return phi

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        alpha = np.zeros(self._phases)
        alpha[0] = 1.0
        return alpha



class HyperExp(ContinuousDistribution, Markovian):
    """
    Hyperexponential distribution (mixture of exponentials).

    The hyperexponential distribution is a mixture of exponential
    distributions. It has SCV >= 1.

    Supports two calling conventions (matching MATLAB API):
        - HyperExp(p, rate1, rate2): 2-phase with probability p of rate1,
          probability (1-p) of rate2
        - HyperExp(probs, rates): n-phase with vectors of probabilities and rates

    Args:
        p_or_probs: For 2-phase: probability of first component (scalar).
                    For n-phase: list/array of probabilities.
        rate1_or_rates: For 2-phase: rate of first component (scalar).
                        For n-phase: list/array of rates.
        rate2: For 2-phase only: rate of second component.
    """

    def __init__(self, p_or_probs: Union[float, list, np.ndarray],
                 rate1_or_rates: Union[float, list, np.ndarray],
                 rate2: Optional[float] = None):
        super().__init__()
        self._name = 'HyperExp'

        # Determine if this is 2-phase (3 scalar args) or n-phase (2 vector args)
        if rate2 is not None:
            # 2-phase form: HyperExp(p, rate1, rate2)
            p = float(p_or_probs)
            r1 = float(rate1_or_rates)
            r2 = float(rate2)
            if not (0 <= p <= 1):
                raise ValueError("Probability p must be in [0, 1]")
            if r1 <= 0 or r2 <= 0:
                raise ValueError("All rates must be positive")
            self._probs = np.array([p, 1.0 - p], dtype=float)
            self._rates = np.array([r1, r2], dtype=float)
        else:
            # n-phase form: HyperExp(probs, rates)
            self._probs = np.array(p_or_probs, dtype=float)
            self._rates = np.array(rate1_or_rates, dtype=float)

            if len(self._probs) != len(self._rates):
                raise ValueError("Probabilities and rates must have same length")
            if not np.allclose(np.sum(self._probs), 1.0):
                raise ValueError("Probabilities must sum to 1")
            if np.any(self._rates <= 0):
                raise ValueError("All rates must be positive")
            if np.any(self._probs < 0):
                raise ValueError("All probabilities must be non-negative")

        self._means = 1.0 / self._rates

    @classmethod
    def fit(cls, mean, scv: float = None, skew: float = None, **kwargs) -> 'HyperExp':
        """Fit a two-phase hyperexponential to three moments (MATLAB
        HyperExp.fit).

        MATLAB tries a Prony fit of the moment triple first and falls back to
        the two-moment fit when it is infeasible. The fallback is used here,
        which is what MATLAB itself returns whenever the triple is not
        hyperexponential-feasible.

        ``HyperExp.fit(dist, method='feldmannwhitt', ...)`` instead fits the
        ccdf of the distribution ``dist`` ITSELF at points spread over decades
        of time scale, rather than matching moments (hyperexp_fit_longtail,
        Feldmann and Whitt 1998). That is the only form available for a
        long-tail law: a Pareto with tail index below 2 has no finite variance,
        so the moment fit above does not exist at all, and even where the
        moments are finite they say nothing about the orders of magnitude over
        which such a law acts. Any further keyword arguments (``k``, ``c1``,
        ``b``, ``decade``, ``points``) are passed through.
        """
        if hasattr(mean, 'evalCDF') or hasattr(mean, 'eval_cdf'):
            dist = mean
            method = str(kwargs.pop('method', 'feldmannwhitt')).lower()
            if method != 'feldmannwhitt':
                raise ValueError("HyperExp.fit on a distribution supports method "
                                 "'feldmannwhitt' only; '%s' was requested." % method)
            from ..api.mam.hyperexp_longtail import hyperexp_fit_longtail
            cdf = getattr(dist, 'evalCDF', None) or getattr(dist, 'eval_cdf')
            res = hyperexp_fit_longtail(lambda t: 1.0 - float(cdf(t)), **kwargs)
            return cls(np.asarray(res['p'], dtype=float).ravel(),
                       np.asarray(res['lambda'], dtype=float).ravel())
        return cls.fit_mean_and_scv(mean, scv)

    @classmethod
    def fit_mean(cls, mean: float) -> 'HyperExp':
        """Two-phase hyperexponential with both rates 1/mean (MATLAB
        HyperExp.fitMean). Both phases share the rate, so the mixing
        probability is immaterial and is taken as 0.5."""
        return cls(0.5, 1.0 / mean, 1.0 / mean)

    @classmethod
    def fit_rate(cls, rate: float) -> 'HyperExp':
        """Two-phase hyperexponential with both rates equal to rate (MATLAB
        HyperExp.fitRate)."""
        return cls(0.5, rate, rate)

    @classmethod
    def fitMean(cls, mean: float) -> 'HyperExp':
        """camelCase alias of fit_mean (MATLAB/JAR spelling)."""
        return cls.fit_mean(mean)

    @classmethod
    def fitRate(cls, rate: float) -> 'HyperExp':
        """camelCase alias of fit_rate (MATLAB/JAR spelling)."""
        return cls.fit_rate(rate)

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float, p: float = 0.99) -> 'HyperExp':
        """
        Create a 2-phase hyperexponential distribution from mean and SCV.

        Uses the same algorithm as MATLAB's map_hyperexp function.

        Args:
            mean: Target mean (MEAN).
            scv: Target squared coefficient of variation (must be >= 1).
            p: Probability of being served in phase 1 (default: 0.99).

        Returns:
            HyperExp distribution with given mean and SCV.
        """
        if scv < 1.0:
            raise ValueError("HyperExp requires SCV >= 1")
        if mean <= 0:
            raise ValueError("Mean must be positive")

        # Port of MATLAB's map_hyperexp algorithm
        # E2 = (1 + SCV) * MEAN^2
        E2 = (1.0 + scv) * mean * mean

        # Delta = -4*p*MEAN^2 + 4*p^2*MEAN^2 + 2*E2*p - 2*E2*p^2
        Delta = -4.0 * p * mean * mean + 4.0 * p * p * mean * mean + 2.0 * E2 * p - 2.0 * E2 * p * p

        if Delta < 0:
            # Try decreasing p if solution not feasible
            if p > 1e-6:
                return cls.fit_mean_and_scv(mean, scv, p / 10.0)
            else:
                raise ValueError(f"Cannot fit HyperExp with mean={mean}, scv={scv}")

        # Try first root
        denom = E2 * p - 2.0 * mean * mean
        if abs(denom) < 1e-12:
            # Avoid division by zero
            if p > 1e-6:
                return cls.fit_mean_and_scv(mean, scv, p / 10.0)
            else:
                raise ValueError(f"Cannot fit HyperExp with mean={mean}, scv={scv}")

        mu2 = (-2.0 * mean + 2.0 * p * mean + np.sqrt(Delta)) / denom
        denom2 = p - 1.0 + mean * mu2
        if abs(denom2) < 1e-12:
            # Try second root
            mu2 = (-2.0 * mean + 2.0 * p * mean - np.sqrt(Delta)) / denom
            denom2 = p - 1.0 + mean * mu2

        mu1 = mu2 * p / denom2

        # Check feasibility (all rates must be positive)
        if mu1 <= 0 or mu2 <= 0 or p < 0 or p > 1:
            # Try second root
            mu2 = (-2.0 * mean + 2.0 * p * mean - np.sqrt(Delta)) / denom
            denom2 = p - 1.0 + mean * mu2
            if abs(denom2) > 1e-12:
                mu1 = mu2 * p / denom2

            # Still not feasible? Try decreasing p
            if mu1 <= 0 or mu2 <= 0:
                if p > 1e-6:
                    return cls.fit_mean_and_scv(mean, scv, p / 10.0)
                else:
                    raise ValueError(f"Cannot fit HyperExp with mean={mean}, scv={scv}")

        # Return HyperExp with p, mu1, mu2
        return cls(p, mu1, mu2)

    @classmethod
    def fit_mean_and_scv_balanced(cls, mean: float, scv: float) -> 'HyperExp':
        """
        Create a 2-phase hyperexponential distribution with balanced means.

        Uses balanced means representation where p/mu1 = (1-p)/mu2.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation (must be >= 1).

        Returns:
            HyperExp distribution with given mean and SCV.
        """
        if scv < 1.0:
            raise ValueError("HyperExp requires SCV >= 1")
        if mean <= 0:
            raise ValueError("Mean must be positive")

        # Port of MATLAB's fitMeanAndSCVBalanced
        mu1 = -(2.0 * (np.sqrt((scv - 1.0) / (scv + 1.0)) / 2.0 - 0.5)) / mean
        p = 0.5 - np.sqrt((scv - 1.0) / (scv + 1.0)) / 2.0

        if mu1 < 0 or p < 0 or p > 1:
            p = np.sqrt((scv - 1.0) / (scv + 1.0)) / 2.0 + 0.5
            mu1 = (2.0 * (np.sqrt((scv - 1.0) / (scv + 1.0)) / 2.0 + 0.5)) / mean

        mu2 = (1.0 - p) / p * mu1

        return cls(float(np.real(p)), float(np.real(mu1)), float(np.real(mu2)))

    # CamelCase alias
    fitMeanAndScvBalanced = fit_mean_and_scv_balanced
    fitMeanAndSCVBalanced = fit_mean_and_scv_balanced

    # CamelCase alias
    fitMeanAndScv = fit_mean_and_scv
    fitMeanAndSCV = fit_mean_and_scv

    @property
    def means(self) -> np.ndarray:
        """Get the means of each component."""
        return self._means

    @property
    def probs(self) -> np.ndarray:
        """Get the probabilities of each component."""
        return self._probs

    @property
    def rates(self) -> np.ndarray:
        """Get the rates of each component."""
        return self._rates

    def getMean(self) -> float:
        """Get the mean."""
        return float(np.sum(self._probs * self._means))

    def getVar(self) -> float:
        """Get the variance."""
        mean = self.getMean()
        second_moment = 2 * np.sum(self._probs * self._means ** 2)
        return second_moment - mean ** 2

    def getSkew(self) -> float:
        """Get the skewness.

        A phase-type mixture has raw moments E[S^k] = k! * sum_i p_i * m_i^k,
        so the third central moment follows in closed form. Without this the
        base-class default returned 0, i.e. a symmetric law, which silently
        understates E[S^3] for every consumer that reconstructs it from the
        skewness (the ForkTail branch variance among them: it read E[S^3] = 13
        instead of 141.55 for the SCV = 4 fit).
        """
        mean = self.getMean()
        var = self.getVar()
        if var <= 0:
            return 0.0
        third_moment = 6 * float(np.sum(self._probs * self._means ** 3))
        return (third_moment - 3 * mean * var - mean ** 3) / var ** 1.5

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        cdf = 0.0
        for p, r in zip(self._probs, self._rates):
            cdf += p * (1.0 - np.exp(-r * x))
        return cdf

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        pdf = 0.0
        for p, r in zip(self._probs, self._rates):
            pdf += p * r * np.exp(-r * x)
        return pdf

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        components = rng.choice(len(self._means), size=n, p=self._probs)
        samples = np.zeros(n)
        for i, c in enumerate(components):
            samples[i] = rng.exponential(scale=self._means[c])
        return samples

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return len(self._means)

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix for MAP representation."""
        k = len(self._rates)
        return np.diag(-self._rates)

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix for MAP representation."""
        k = len(self._rates)
        D1 = np.zeros((k, k))
        for i in range(k):
            for j in range(k):
                D1[i, j] = self._rates[i] * self._probs[j]
        return D1

    def getMu(self) -> np.ndarray:
        """Get the service rates in each phase."""
        return self._rates.copy()

    def getPhi(self) -> np.ndarray:
        """Get the completion probabilities from each phase."""
        return np.ones(len(self._rates))

    def getInitProb(self) -> np.ndarray:
        """Get the initial probability vector."""
        return self._probs.copy()



class Gamma(ContinuousDistribution):
    """
    Gamma distribution.

    The gamma distribution is a two-parameter continuous distribution
    that generalizes the exponential and Erlang distributions.

    Args:
        shape: Shape parameter (k or alpha).
        scale: Scale parameter (theta).
    """

    def __init__(self, shape: float, scale: float):
        super().__init__()
        self._name = 'Gamma'
        if shape <= 0:
            raise ValueError("Shape must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self._shape = shape
        self._scale = scale

    @property
    def shape(self) -> float:
        """Get the shape parameter."""
        return self._shape

    @property
    def scale(self) -> float:
        """Get the scale parameter."""
        return self._scale

    def getMean(self) -> float:
        """Get the mean (shape * scale)."""
        return self._shape * self._scale

    def getVar(self) -> float:
        """Get the variance (shape * scale^2)."""
        return self._shape * self._scale ** 2

    def getSCV(self) -> float:
        """Get the SCV (1/shape)."""
        return 1.0 / self._shape

    def getSkew(self) -> float:
        """Get the skewness."""
        return 2.0 / np.sqrt(self._shape)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        return stats.gamma.cdf(x, a=self._shape, scale=self._scale)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return stats.gamma.pdf(x, a=self._shape, scale=self._scale)

    def evalLST(self, s):
        """
        LST of the Gamma law, (beta/(s+beta))^shape with beta = 1/scale.

        MATLAB, the JAR and the cpp port all carry this closed form; without it
        the base rectangle rule answered here, which is both approximate (3.3e-4
        relative at s = 0.5+1i) and real-only. Being analytic it also serves the
        complex arguments transform inversion needs.
        """
        shape = self._shape
        scale = self._scale
        beta = 1.0 / scale
        val = (beta / (s + beta)) ** shape
        return complex(val) if isinstance(s, complex) else float(np.real(val))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.gamma(shape=self._shape, scale=self._scale, size=n)

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Gamma':
        """
        Create a Gamma distribution from mean and SCV.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.

        Returns:
            Gamma distribution with given mean and SCV.
        """
        # For Gamma: SCV = 1/shape, mean = shape * scale
        # So: shape = 1/scv, scale = mean * scv
        if scv <= 0:
            raise ValueError("SCV must be positive")
        shape = 1.0 / scv
        scale = mean * scv
        return cls(shape, scale)

    # CamelCase alias
    fitMeanAndSCV = fit_mean_and_scv
    fitMeanAndScv = fit_mean_and_scv


class Lognormal(ContinuousDistribution):
    """
    Lognormal distribution.

    A random variable X has a lognormal distribution if log(X) is
    normally distributed.

    Args:
        mu: Mean of the underlying normal distribution.
        sigma: Standard deviation of the underlying normal distribution.
    """

    def __init__(self, mu: float, sigma: float):
        super().__init__()
        self._name = 'Lognormal'
        if sigma <= 0:
            raise ValueError("Sigma must be positive")
        self._mu = mu
        self._sigma = sigma

    @property
    def mu(self) -> float:
        """Get the mu parameter."""
        return self._mu

    @property
    def sigma(self) -> float:
        """Get the sigma parameter."""
        return self._sigma

    def getMean(self) -> float:
        """Get the mean."""
        return np.exp(self._mu + self._sigma ** 2 / 2)

    def getVar(self) -> float:
        """Get the variance."""
        return (np.exp(self._sigma ** 2) - 1) * np.exp(2 * self._mu + self._sigma ** 2)

    def getSkew(self) -> float:
        """Get the skewness."""
        es2 = np.exp(self._sigma ** 2)
        return (es2 + 2) * np.sqrt(es2 - 1)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x <= 0:
            return 0.0
        return stats.lognorm.cdf(x, s=self._sigma, scale=np.exp(self._mu))

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x <= 0:
            return 0.0
        return stats.lognorm.pdf(x, s=self._sigma, scale=np.exp(self._mu))

    def evalLST(self, s):
        """Numerical LST (rectangle rule, n=1000) matching MATLAB Lognormal.evalLST.

        numpy rather than math for the kernel, so a COMPLEX argument is
        admissible: transform inversion and root location both need one, and
        MATLAB's own quadrature extends to the complex plane unchanged.
        """
        import math
        mu = self._mu
        sigma = self._sigma
        upper = math.exp(mu + 5.0 * sigma)
        n = 1000
        dx = upper / n
        cplx = isinstance(s, complex)
        total = 0.0 + 0.0j if cplx else 0.0
        for i in range(1, n + 1):
            x = i * dx
            logx = math.log(x)
            pdf = math.exp(-(logx - mu) ** 2 / (2.0 * sigma ** 2)) / (x * sigma * math.sqrt(2.0 * math.pi))
            total += np.exp(-s * x) * pdf
        return complex(total * dx) if cplx else float(np.real(total * dx))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.lognormal(mean=self._mu, sigma=self._sigma, size=n)

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Lognormal':
        """
        Construct a Lognormal from a target mean and squared coefficient of
        variation (SCV = variance/mean^2), converting to log-space (mu, sigma).

        Port of MATLAB Lognormal.fitMeanAndSCV.
        """
        if mean <= 0:
            raise ValueError("Mean must be positive")
        if scv <= 0:
            raise ValueError("SCV must be positive")
        c2 = scv  # c = sqrt(scv); c*c == scv
        mu = np.log(mean / np.sqrt(c2 + 1.0))
        sigma = np.sqrt(np.log(c2 + 1.0))
        return cls(float(mu), float(sigma))

    # CamelCase aliases
    fitMeanAndSCV = fit_mean_and_scv
    fitMeanAndScv = fit_mean_and_scv



class Pareto(ContinuousDistribution):
    """
    Pareto distribution.

    The Pareto distribution is a power-law distribution often used
    to model heavy-tailed phenomena.

    Args:
        alpha: Shape parameter (tail index).
        scale: Scale parameter (minimum value).
    """

    def __init__(self, alpha: float, scale: float):
        super().__init__()
        self._name = 'Pareto'
        if alpha <= 0:
            raise ValueError("Alpha must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self._alpha = alpha
        self._scale = scale

    @property
    def alpha(self) -> float:
        """Get the alpha parameter."""
        return self._alpha

    @property
    def scale(self) -> float:
        """Get the scale parameter."""
        return self._scale

    def getMean(self) -> float:
        """Get the mean."""
        if self._alpha <= 1:
            return float('inf')
        return self._alpha * self._scale / (self._alpha - 1)

    def getVar(self) -> float:
        """Get the variance."""
        if self._alpha <= 2:
            return float('inf')
        return (self._scale ** 2 * self._alpha) / ((self._alpha - 1) ** 2 * (self._alpha - 2))

    def getSkew(self) -> float:
        """Get the skewness.

        For a Pareto law with shape alpha the third moment exists only when
        alpha > 3, and the skewness is 2*(1+alpha)/(alpha-3)*sqrt((alpha-2)/alpha).
        Without this the base-class default returned 0, i.e. a symmetric law,
        which silently understates E[S^3] for every consumer that reconstructs
        it from the skewness.
        """
        if self._alpha <= 3:
            return float('inf')
        a = self._alpha
        return 2.0 * (1.0 + a) / (a - 3.0) * np.sqrt((a - 2.0) / a)

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [scale, inf)."""
        return (self._scale, float('inf'))

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < self._scale:
            return 0.0
        return 1.0 - (self._scale / x) ** self._alpha

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < self._scale:
            return 0.0
        return self._alpha * self._scale ** self._alpha / x ** (self._alpha + 1)

    def evalLST(self, s: float) -> float:
        """Laplace-Stieltjes transform E[e^{-sX}] of the Pareto distribution.

        ``A*(s) = int_k^inf e^{-sx} alpha k^alpha x^{-(alpha+1)} dx``. Substituting
        ``x = k/u`` maps the infinite tail onto a unit interval and cancels the
        scale exactly::

            A*(s) = alpha * int_0^1 u^(alpha-1) exp(-s*k/u) du

        This is the same transform as the closed form of Nadarajah & Kotz,
        ``A*(s) = alpha*(s*k)^alpha*Gamma(-alpha, s*k) = alpha*E_{alpha+1}(s*k)``
        (Queueing Syst (2006) 54:243-244, DOI 10.1007/s11134-006-0299-1), but in a
        form that stays accurate as ``s -> 0``, where the incomplete-gamma product
        underflows to 0/inf. Here ``s = 0`` gives ``alpha*int_0^1 u^(alpha-1) du =
        1`` exactly, and the integrand is bounded and C^inf on a FINITE interval
        for ``alpha >= 2`` (the shape floor the constructor enforces).

        Accuracy: adaptive Gauss-Kronrod at 1e-12 relative, matching the MATLAB and
        JAR implementations, verified against mpmath to 1e-15. The previous
        implementation was a 1000-point right-endpoint rectangle sum truncated at
        ``k*1000**(1/alpha)``; it lost the mass beyond the truncation point and
        biased the transform low by ~3.1% at alpha=2.0078 (it returned
        ``A*(0)=0.96914``, not 1).
        """
        from scipy.integrate import quad
        alpha = self._alpha
        k = self._scale
        if s == 0.0:
            return 1.0  # A*(0) = 1 exactly; skip the quadrature

        def integrand(u):
            # u=0 is an essential zero of the integrand (exp(-s*k/u) and all its
            # derivatives vanish there), so the guard only avoids 0/0.
            if u <= 0.0:
                return 0.0
            return (u ** (alpha - 1.0)) * np.exp(-s * k / u)

        # see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
        if isinstance(s, complex):
            # quad integrates a REAL integrand, so the two parts are taken
            # separately; the interval and the rule are otherwise unchanged.
            re, _ = quad(lambda u: float(np.real(integrand(u))), 0.0, 1.0,
                         epsabs=0.0, epsrel=1e-12, limit=200)
            im, _ = quad(lambda u: float(np.imag(integrand(u))), 0.0, 1.0,
                         epsabs=0.0, epsrel=1e-12, limit=200)
            return alpha * complex(re, im)
        val, _ = quad(integrand, 0.0, 1.0, epsabs=0.0, epsrel=1e-12, limit=200)
        return alpha * val

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return (rng.pareto(a=self._alpha, size=n) + 1) * self._scale

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Pareto':
        """
        Create a Pareto distribution from mean and SCV.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.

        Returns:
            Pareto distribution with given mean and SCV.

        Note:
            For Pareto with alpha > 2:
            mean = alpha * scale / (alpha - 1)
            var = scale^2 * alpha / ((alpha - 1)^2 * (alpha - 2))
            scv = var / mean^2 = 1 / (alpha * (alpha - 2))

            Solving for alpha: alpha = (1 + sqrt(1 + 4*scv)) / (2*scv)
            Then: scale = mean * (alpha - 1) / alpha
        """
        if scv <= 0:
            raise ValueError("SCV must be positive")

        # see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
        discriminant = 1.0 + 1.0 / scv
        alpha = 1.0 + np.sqrt(discriminant)

        if alpha <= 2:
            # For very high SCV, use minimum alpha = 2.01 to ensure finite variance
            alpha = 2.01

        scale = mean * (alpha - 1) / alpha
        return cls(alpha, scale)

    # CamelCase alias
    fitMeanAndSCV = fit_mean_and_scv
    fitMeanAndScv = fit_mean_and_scv


class Uniform(ContinuousDistribution):
    """
    Uniform distribution on [min, max].

    Args:
        min_val: Minimum value.
        max_val: Maximum value.
    """

    def __init__(self, min_val: float, max_val: float):
        super().__init__()
        self._name = 'Uniform'
        if max_val < min_val:
            raise ValueError("max_val must be >= min_val")
        self._min = min_val
        self._max = max_val

    @property
    def min_val(self) -> float:
        """Get the minimum value."""
        return self._min

    @property
    def max_val(self) -> float:
        """Get the maximum value."""
        return self._max

    def getMean(self) -> float:
        """Get the mean."""
        return (self._min + self._max) / 2

    def getVar(self) -> float:
        """Get the variance."""
        return (self._max - self._min) ** 2 / 12

    def getSkew(self) -> float:
        """Get the skewness (always 0 for uniform)."""
        return 0.0

    def getSupport(self) -> Tuple[float, float]:
        """Get the support [min, max]."""
        return (self._min, self._max)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < self._min:
            return 0.0
        if x > self._max:
            return 1.0
        return (x - self._min) / (self._max - self._min)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < self._min or x > self._max:
            return 0.0
        return 1.0 / (self._max - self._min)

    def evalLST(self, s: float) -> float:
        """LST of Uniform[min,max]: (e^{-s*min}-e^{-s*max})/(s*(max-min)). Matches
        MATLAB. numpy rather than math, so a COMPLEX argument is admissible."""
        if abs(s) < 1e-14:
            return 1.0
        val = (np.exp(-s * self._min) - np.exp(-s * self._max)) / (s * (self._max - self._min))
        return complex(val) if isinstance(s, complex) else float(np.real(val))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.uniform(self._min, self._max, size=n)



class Weibull(ContinuousDistribution):
    """
    Weibull distribution.

    The Weibull distribution is commonly used in reliability engineering
    to model time to failure.

    Args:
        shape: Shape parameter (k).
        scale: Scale parameter (lambda).
    """

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Weibull':
        """
        Fit a Weibull to a mean and squared coefficient of variation.

        Port of MATLAB Weibull.fitMeanAndSCV and jline.lang.processes.Weibull.
        The shape comes from the Justus et al. (1976) approximation
        k = CV^(-1.086) with CV = sqrt(scv); the scale then makes the MEAN
        exact, scale = mean / Gamma(1 + 1/k). Only the SCV is approximate: the
        error is below 3% inside the range the approximation was published for
        (k in [1,10], i.e. scv <= 1) and grows quickly outside it (12% at
        scv = 2, 48% at scv = 4).

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.
        """
        from scipy.special import gamma as gamma_fn
        c = np.sqrt(scv)
        shape = c ** (-1.086)          # Justus approximation (1976)
        scale = mean / gamma_fn(1 + 1.0 / shape)
        return cls(shape, scale)

    @classmethod
    def fitMeanAndSCV(cls, mean: float, scv: float) -> 'Weibull':
        """camelCase alias of fit_mean_and_scv (MATLAB/JAR spelling)."""
        return cls.fit_mean_and_scv(mean, scv)

    def __init__(self, shape: float, scale: float):
        super().__init__()
        self._name = 'Weibull'
        if shape <= 0:
            raise ValueError("Shape must be positive")
        if scale <= 0:
            raise ValueError("Scale must be positive")
        self._shape = shape
        self._scale = scale

    @property
    def shape(self) -> float:
        """Get the shape parameter."""
        return self._shape

    @property
    def scale(self) -> float:
        """Get the scale parameter."""
        return self._scale

    def getMean(self) -> float:
        """Get the mean."""
        from scipy.special import gamma
        return self._scale * gamma(1 + 1 / self._shape)

    def getVar(self) -> float:
        """Get the variance."""
        from scipy.special import gamma
        g1 = gamma(1 + 1 / self._shape)
        g2 = gamma(1 + 2 / self._shape)
        return self._scale ** 2 * (g2 - g1 ** 2)

    def getSkew(self) -> float:
        """Get the skewness.

        With g_k = Gamma(1 + k/shape), the Weibull third central moment gives
        skew = (g3 - 3*g1*g2 + 2*g1^3) / (g2 - g1^2)^1.5, scale-free. Without
        this the base-class default returned 0, i.e. a symmetric law.
        """
        from scipy.special import gamma
        g1 = gamma(1 + 1 / self._shape)
        g2 = gamma(1 + 2 / self._shape)
        g3 = gamma(1 + 3 / self._shape)
        var_std = g2 - g1 ** 2
        if var_std <= 0:
            return 0.0
        return (g3 - 3 * g1 * g2 + 2 * g1 ** 3) / var_std ** 1.5

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        if x < 0:
            return 0.0
        return 1.0 - np.exp(-(x / self._scale) ** self._shape)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        if x < 0:
            return 0.0
        return (self._shape / self._scale) * (x / self._scale) ** (self._shape - 1) * \
               np.exp(-(x / self._scale) ** self._shape)

    def evalLST(self, s):
        """Numerical LST (rectangle rule, n=1000) matching MATLAB Weibull.evalLST.

        numpy rather than math for the kernel, so a COMPLEX argument is
        admissible, as in MATLAB.
        """
        import math
        alpha = self._scale  # MATLAB scale param
        r = self._shape      # MATLAB shape param
        upper = alpha * ((-math.log(1e-10)) ** (1.0 / r))
        n = 1000
        dx = upper / n
        cplx = isinstance(s, complex)
        total = 0.0 + 0.0j if cplx else 0.0
        for i in range(1, n + 1):
            x = i * dx
            pdf = (r / alpha) * ((x / alpha) ** (r - 1.0)) * math.exp(-((x / alpha) ** r))
            total += np.exp(-s * x) * pdf
        return complex(total * dx) if cplx else float(np.real(total * dx))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return self._scale * rng.weibull(a=self._shape, size=n)



class Normal(ContinuousDistribution):
    """
    Normal (Gaussian) distribution.

    Note: For queueing applications, a truncated or shifted version
    may be needed since normal distributions can take negative values.

    Args:
        mean: Mean of the distribution.
        std: Standard deviation.
    """

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Normal':
        """
        Fit a normal to a mean and squared coefficient of variation.

        Port of MATLAB Normal.fitMeanAndSCV: var = scv*mean^2, std = sqrt(var),
        floored at the fine tolerance as MATLAB does.
        """
        var = scv * mean ** 2
        return cls(mean, max(1e-8, np.sqrt(var)))

    @classmethod
    def fitMeanAndSCV(cls, mean: float, scv: float) -> 'Normal':
        """camelCase alias of fit_mean_and_scv (MATLAB/JAR spelling)."""
        return cls.fit_mean_and_scv(mean, scv)

    def __init__(self, mean: float, std: float):
        super().__init__()
        self._name = 'Normal'
        if std <= 0:
            raise ValueError("Standard deviation must be positive")
        self._mean_val = mean
        self._std = std

    @property
    def std(self) -> float:
        """Get the standard deviation."""
        return self._std

    def getMean(self) -> float:
        """Get the mean."""
        return self._mean_val

    def getVar(self) -> float:
        """Get the variance."""
        return self._std ** 2

    def getSCV(self) -> float:
        """Get the squared coefficient of variation (var/mean^2)."""
        if self._mean_val == 0:
            return float('inf')
        return (self._std ** 2) / (self._mean_val ** 2)

    def getStd(self) -> float:
        """Get the standard deviation."""
        return self._std

    def getSkew(self) -> float:
        """Get the skewness (always 0 for normal)."""
        return 0.0

    def getSupport(self) -> Tuple[float, float]:
        """Get the support (-inf, inf)."""
        return (float('-inf'), float('inf'))

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x."""
        return stats.norm.cdf(x, loc=self._mean_val, scale=self._std)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF at point x."""
        return stats.norm.pdf(x, loc=self._mean_val, scale=self._std)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random samples."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.normal(self._mean_val, self._std, size=n)

    @classmethod
    def fitMean(cls, mean: float, std: float = 1.0) -> 'Normal':
        """Create a Normal distribution with given mean and std."""
        return cls(mean, std)

    @classmethod
    def fitMeanAndStd(cls, mean: float, std: float) -> 'Normal':
        """Create a Normal distribution with given mean and std."""
        return cls(mean, std)

    @classmethod
    def fitMeanAndVar(cls, mean: float, var: float) -> 'Normal':
        """Create a Normal distribution with given mean and variance."""
        return cls(mean, np.sqrt(var))

    def getSkewness(self) -> float:
        """Get the skewness (always 0 for normal)."""
        return 0.0

    # Snake_case aliases
    fit_mean = fitMean
    fit_mean_and_std = fitMeanAndStd
    fit_mean_and_var = fitMeanAndVar
    get_mean = lambda self: self._mean_val
    get_std = getStd
    get_var = getVar
    get_scv = getSCV
    get_skewness = getSkewness
    eval_cdf = evalCDF
    eval_pdf = evalPDF


class MultivariateNormal(ContinuousDistribution):
    """
    Multivariate Normal (Gaussian) distribution.

    Represents a d-dimensional normal distribution with mean vector mu
    and covariance matrix Sigma.

    Args:
        mu: d-dimensional mean vector.
        Sigma: d x d positive definite covariance matrix.
    """

    def __init__(self, mu: Union[list, np.ndarray], Sigma: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'MultivariateNormal'

        self._mu = np.atleast_1d(np.array(mu, dtype=float)).flatten()
        self._Sigma = np.atleast_2d(np.array(Sigma, dtype=float))

        d = len(self._mu)
        if self._Sigma.shape != (d, d):
            raise ValueError(f"Sigma must be {d}x{d} to match mu of length {d}")

        # Check positive definite via Cholesky
        try:
            self._L = linalg.cholesky(self._Sigma, lower=True)
        except linalg.LinAlgError:
            raise ValueError("Sigma must be positive definite")

        self._dimension = d

    @property
    def dimension(self) -> int:
        """Get the dimensionality."""
        return self._dimension

    def getMeanVector(self) -> np.ndarray:
        """Get the mean vector."""
        return self._mu.copy()

    def getCovariance(self) -> np.ndarray:
        """Get the covariance matrix."""
        return self._Sigma.copy()

    def getCorrelation(self) -> np.ndarray:
        """Get the correlation matrix."""
        d = self._dimension
        R = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                std_i = np.sqrt(self._Sigma[i, i])
                std_j = np.sqrt(self._Sigma[j, j])
                if std_i > 1e-10 and std_j > 1e-10:
                    R[i, j] = self._Sigma[i, j] / (std_i * std_j)
                else:
                    R[i, j] = float(i == j)
        return R

    def getMean(self) -> float:
        """Get the mean of the first component (for compatibility)."""
        return float(self._mu[0])

    def getVar(self) -> float:
        """Get the variance of the first component."""
        return float(self._Sigma[0, 0])

    def getSkew(self) -> float:
        """Get skewness (0 for normal)."""
        return 0.0

    def evalPDF(self, x: Union[list, np.ndarray]) -> Union[float, np.ndarray]:
        """Evaluate the multivariate normal PDF at point(s) x.

        Args:
            x: Single point (1D array of length d) or multiple points (2D array of shape n x d)

        Returns:
            Single float for one point, or numpy array for multiple points.
        """
        x_arr = np.atleast_1d(np.array(x, dtype=float))

        # Handle 2D array (multiple points)
        if x_arr.ndim == 2:
            n_points = x_arr.shape[0]
            if x_arr.shape[1] != self._dimension:
                raise ValueError(f"Each point must have dimension {self._dimension}")

            inv_Sigma = linalg.inv(self._Sigma)
            det_Sigma = linalg.det(self._Sigma)
            norm_const = 1.0 / np.sqrt((2 * np.pi) ** self._dimension * det_Sigma)

            results = np.zeros(n_points)
            for i in range(n_points):
                diff = x_arr[i] - self._mu
                exponent = -0.5 * diff @ inv_Sigma @ diff
                results[i] = norm_const * np.exp(exponent)
            return results

        # Handle 1D array (single point)
        x_arr = x_arr.flatten()
        if len(x_arr) != self._dimension:
            raise ValueError(f"x must have dimension {self._dimension}")

        diff = x_arr - self._mu
        inv_Sigma = linalg.inv(self._Sigma)
        det_Sigma = linalg.det(self._Sigma)

        norm_const = 1.0 / np.sqrt((2 * np.pi) ** self._dimension * det_Sigma)
        exponent = -0.5 * diff @ inv_Sigma @ diff

        return float(norm_const * np.exp(exponent))

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """
        Generate n samples from the multivariate normal.

        Returns:
            n x d matrix of samples.
        """
        if rng is None:
            rng = np.random.default_rng()

        # X = mu + L @ Z where L = chol(Sigma), Z ~ N(0, I)
        Z = rng.standard_normal((self._dimension, n))
        X = (self._mu[:, np.newaxis] + self._L @ Z).T

        return X

    def getMarginal(self, indices: Union[list, np.ndarray]) -> 'MultivariateNormal':
        """Extract a marginal distribution for a subset of dimensions."""
        indices = np.atleast_1d(np.array(indices, dtype=int))
        mu_marg = self._mu[indices]
        Sigma_marg = self._Sigma[np.ix_(indices, indices)]
        return MultivariateNormal(mu_marg, Sigma_marg)

    def getMarginalUniv(self, index: int) -> 'Normal':
        """Extract a univariate marginal distribution."""
        mean_marg = self._mu[index]
        std_marg = np.sqrt(self._Sigma[index, index])
        return Normal(mean_marg, std_marg)

    def getDimension(self) -> int:
        """Get the dimensionality."""
        return self._dimension

    # Snake_case aliases
    def get_dimension(self) -> int:
        """Get the dimensionality."""
        return self._dimension

    def get_mean_vector(self) -> np.ndarray:
        """Get the mean vector."""
        return self._mu.copy()

    def get_covariance(self) -> np.ndarray:
        """Get the covariance matrix."""
        return self._Sigma.copy()

    def get_correlation(self) -> np.ndarray:
        """Get the correlation matrix."""
        return self.getCorrelation()

    def get_marginal(self, indices: Union[list, np.ndarray]) -> 'MultivariateNormal':
        """Extract a marginal distribution for a subset of dimensions."""
        return self.getMarginal(indices)

    def get_marginal_univ(self, index: int) -> 'Normal':
        """Extract a univariate marginal distribution."""
        return self.getMarginalUniv(index)

    def eval_pdf(self, x: Union[list, np.ndarray]) -> Union[float, np.ndarray]:
        """Evaluate the multivariate normal PDF at point(s) x."""
        x_arr = np.atleast_1d(np.array(x, dtype=float))

        # Handle multiple points: if x is 2D (n_points x d)
        if x_arr.ndim == 2:
            results = np.zeros(x_arr.shape[0])
            for i in range(x_arr.shape[0]):
                results[i] = self.evalPDF(x_arr[i])
            return results
        else:
            return self.evalPDF(x_arr)

    @classmethod
    def fitMeanAndCovariance(cls, mu: Union[list, np.ndarray],
                              Sigma: Union[list, np.ndarray]) -> 'MultivariateNormal':
        """Create a MultivariateNormal distribution with given mean and covariance."""
        return cls(mu, Sigma)

    # CamelCase to snake_case alias for fit method
    fit_mean_and_covariance = fitMeanAndCovariance


class Prior(ContinuousDistribution):
    """
    Discrete prior distribution over alternative distributions.

    Prior represents parameter uncertainty by specifying a discrete set of
    alternative distributions with associated probabilities. Used with the
    UQ solver for Bayesian-style analysis.

    This is NOT a mixture distribution - each alternative represents a
    separate model realization.

    Args:
        distributions: List of Distribution objects.
        probabilities: List of probabilities (must sum to 1).
    """

    def __init__(self, distributions: list, probabilities: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'Prior'

        if not isinstance(distributions, list) or len(distributions) == 0:
            raise ValueError("distributions must be a non-empty list")

        self._distributions = distributions
        self._probabilities = np.array(probabilities, dtype=float)

        if len(self._distributions) != len(self._probabilities):
            raise ValueError("Number of distributions must match number of probabilities")

        if np.any(self._probabilities < 0):
            raise ValueError("Probabilities must be non-negative")

        if not np.isclose(self._probabilities.sum(), 1.0, atol=1e-6):
            raise ValueError(f"Probabilities must sum to 1 (got {self._probabilities.sum()})")

    @property
    def distributions(self) -> list:
        """Get the alternative distributions."""
        return self._distributions

    @property
    def probabilities(self) -> np.ndarray:
        """Get the probabilities."""
        return self._probabilities.copy()

    def getNumAlternatives(self) -> int:
        """Get the number of alternative distributions."""
        return len(self._distributions)

    def getAlternative(self, idx: int):
        """Get the distribution at index idx."""
        if idx < 0 or idx >= len(self._distributions):
            raise IndexError("Index out of bounds")
        return self._distributions[idx]

    def getProbability(self, idx: int) -> float:
        """Get the probability of alternative idx."""
        if idx < 0 or idx >= len(self._probabilities):
            raise IndexError("Index out of bounds")
        return float(self._probabilities[idx])

    def discretize(self, n: int = 11, method: str = 'quadrature', rng=None):
        """Reduce the prior to n weighted alternatives.

        The method is honoured, which is what makes SolverUQ's own
        'quadrature'/'montecarlo' methods mean anything:

          'quadrature'  the alternatives and their probabilities unchanged, and
                        n is ignored: a discrete set is already exact.
          'montecarlo'  n i.i.d. draws of the ALTERNATIVE INDEX against its
                        probabilities, weights 1/n. Returning the alternatives
                        unweighted here would silently drop the prior.

        Mirrors Prior.discretize in MATLAB. The continuous form of the prior --
        a parameter density plus a distribution factory -- is not representable
        by this class, which holds an explicit alternative list, so a continuous
        request has nothing to discretize and is refused rather than answered
        from the discrete branch.

        Args:
            n: number of alternatives (ignored by 'quadrature')
            method: 'quadrature' (default) or 'montecarlo'
            rng: optional numpy Generator, for a reproducible draw

        Returns:
            (dists, weights) with weights summing to 1
        """
        if method not in ('quadrature', 'montecarlo'):
            raise ValueError("Unknown discretization method: %s" % method)
        if method == 'quadrature':
            return list(self._distributions), self._probabilities.copy()
        n = int(n) if n else 11
        if n < 1:
            raise ValueError("montecarlo needs at least one draw")
        draws = rng if rng is not None else np.random.default_rng()
        cumprob = np.cumsum(self._probabilities)
        dists = []
        for _ in range(n):
            u = float(draws.random())
            idx = int(np.searchsorted(cumprob, u, side='left'))
            idx = min(idx, len(self._distributions) - 1)
            dists.append(self._distributions[idx])
        return dists, np.full(n, 1.0 / n)

    def getMean(self) -> float:
        """Get prior-weighted mean (expected mean over alternatives)."""
        mean = 0.0
        for i in range(len(self._distributions)):
            mean += self._probabilities[i] * self._distributions[i].getMean()
        return mean

    def getVar(self) -> float:
        """Get prior-weighted variance using law of total variance."""
        E_mean = 0.0       # E[E[X|D]]
        E_var = 0.0        # E[Var(X|D)]
        E_mean_sq = 0.0    # E[E[X|D]^2]

        for i in range(len(self._distributions)):
            m = self._distributions[i].getMean()
            v = self._distributions[i].getVar()
            E_mean += self._probabilities[i] * m
            E_var += self._probabilities[i] * v
            E_mean_sq += self._probabilities[i] * m ** 2

        # Var(X) = E[Var(X|D)] + Var(E[X|D])
        return E_var + (E_mean_sq - E_mean ** 2)

    def getSCV(self) -> float:
        """Get prior-weighted SCV."""
        mean = self.getMean()
        var = self.getVar()
        return var / mean ** 2 if mean > 0 else 0.0

    def evalCDF(self, t: float) -> float:
        """Evaluate mixture CDF at t."""
        cdf = 0.0
        for i in range(len(self._distributions)):
            cdf += self._probabilities[i] * self._distributions[i].evalCDF(t)
        return cdf

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Sample from prior (mixture sampling)."""
        if rng is None:
            rng = np.random.default_rng()

        samples = np.zeros(n)
        cumprob = np.cumsum(self._probabilities)

        for i in range(n):
            r = rng.random()
            idx = np.searchsorted(cumprob, r)
            idx = min(idx, len(self._distributions) - 1)
            samples[i] = self._distributions[idx].sample(1, rng)[0]

        return samples

    def isPrior(self) -> bool:
        """Return True (used for detection by UQ solver)."""
        return True

    def isPriorDistribution(self) -> bool:
        """Alias for isPrior (used for detection by UQ solver)."""
        return True

    def getProbabilities(self) -> np.ndarray:
        """Get all probabilities."""
        return self._probabilities.copy()

    # Snake_case aliases
    get_num_alternatives = getNumAlternatives
    get_alternative = getAlternative
    get_probability = getProbability
    get_probabilities = getProbabilities
    is_prior = isPrior
    is_prior_distribution = isPriorDistribution


class Expolynomial(ContinuousDistribution):
    """
    Expolynomial distribution with density f(x) = sum ci * x^ai * exp(-li*x).

    Represents an expolynomial density over a bounded domain [eft, lft],
    matching the GEN expolynomial format of external stochastic Petri net tools.

    Args:
        density: Density expression string in expolynomial (GEN) format.
        eft: Earliest firing time (lower bound of support).
        lft: Latest firing time (upper bound of support, use math.inf for unbounded).
    """

    def __init__(self, density: str, eft: float, lft: float):
        super().__init__()
        self._name = 'Expolynomial'
        self._density = density
        self._eft = float(eft)
        self._lft = float(lft)

    @property
    def density(self) -> str:
        """Get the density expression string."""
        return self._density

    @property
    def eft(self) -> float:
        """Get the earliest firing time."""
        return self._eft

    @property
    def lft(self) -> float:
        """Get the latest firing time."""
        return self._lft

    def getMean(self) -> float:
        """Get the mean (returns NaN - numerical integration not supported in Python)."""
        return float('nan')

    def getVar(self) -> float:
        """Get the variance (returns NaN)."""
        return float('nan')

    def getSCV(self) -> float:
        """Get the squared coefficient of variation (returns NaN)."""
        return float('nan')

    def getRate(self) -> float:
        """Get the rate 1/mean (returns NaN)."""
        return float('nan')

    def getSupport(self):
        """Get the support [eft, lft]."""
        return (self._eft, self._lft)

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF at point x (returns NaN - not supported)."""
        return float('nan')

    def sample(self, n: int = 1, rng=None) -> np.ndarray:
        """Generate random samples (returns NaN - not supported)."""
        return np.full(n, float('nan'))


class NHPP(ContinuousDistribution):
    """
    Non-homogeneous Poisson process (NHPP) with a piecewise-constant intensity.

    The intensity is a step function of the wall clock: segment i covers
    [breakpoints[i], breakpoints[i+1]) and carries rate rates[i], so breakpoints
    has one more entry than rates. With cyclic=True the schedule repeats,
    giving a cyclic Poisson process.

    Two horizon conventions:
      cyclic     : the schedule repeats with period
                   T = breakpoints[-1] - breakpoints[0]; the active segment at
                   time t follows from (t - breakpoints[0]) % T.
      non-cyclic : the intensity is zero outside
                   [breakpoints[0], breakpoints[-1]), so the process emits
                   nothing once the schedule is exhausted. A non-cyclic NHPP is
                   therefore a transient construct: run to steady state it
                   converges to the empty system, so callers should use a time
                   span within the horizon.

    This is NOT a renewal process. Successive intervals are dependent, because
    the position within the schedule carries over from one event to the next.
    Accordingly the scalar summaries that presuppose an i.i.d. interval
    distribution -- getSCV, getSkew, getVar, evalCDF -- are undefined and return
    NaN rather than a representative exponential value, which would silently
    misreport the process as Poisson. The schedule is the parameterisation: read
    it with getRateSchedule. getMean is well defined and returns the
    arrival-stationary (Palm) mean interval 1/timeAverageRate.

    Solver support: the LDES simulation engine honours the exact schedule in
    both steady state (cyclic only) and transient analysis. SolverFLD honours it
    in getTranAvg, by injecting the intensity as a time-varying rate multiplier
    on the closing ODE; SolverFLD.getAvg uses the time-average rate, which is
    the steady state of a cyclic schedule. Every other solver rejects a model
    using it via the standard unsupported-feature check.

    Args:
        breakpoints: strictly increasing segment boundaries, length n+1.
        rates: non-negative rate on each segment, length n.
        cyclic: whether the schedule repeats with the horizon as period.
    """

    def __init__(self, breakpoints, rates, cyclic: bool = True):
        super().__init__()
        self._name = 'NHPP'
        breakpoints = np.asarray(breakpoints, dtype=float).ravel()
        rates = np.asarray(rates, dtype=float).ravel()
        if rates.size == 0 or breakpoints.size != rates.size + 1:
            raise ValueError(
                "NHPP: breakpoints must be non-empty with one more entry than rates")
        if np.any(np.diff(breakpoints) <= 0):
            raise ValueError("NHPP: breakpoints must be strictly increasing")
        if np.any(rates < 0) or np.any(np.isinf(rates)):
            raise ValueError("NHPP: rates must be finite and non-negative")
        if float(np.sum(rates * np.diff(breakpoints))) <= 0:
            raise ValueError(
                "NHPP: the schedule has zero total intensity, so no event can ever occur")
        self._breakpoints = breakpoints
        self._rates = rates
        self._cyclic = bool(cyclic)
        # Wall-clock position of the next sample; see sample().
        self._sample_clock = float(breakpoints[0])

    @property
    def breakpoints(self) -> np.ndarray:
        """Segment boundaries, length n+1."""
        return self._breakpoints

    @property
    def rates(self) -> np.ndarray:
        """Per-segment rates, length n."""
        return self._rates

    @property
    def cyclic(self) -> bool:
        """Whether the schedule repeats."""
        return self._cyclic

    def getBreakpoints(self) -> np.ndarray:
        """Segment boundaries, length n+1 (MATLAB/JAR accessor name)."""
        return self._breakpoints

    def getRates(self) -> np.ndarray:
        """Per-segment rates, length n (MATLAB/JAR accessor name)."""
        return self._rates

    def isCyclic(self) -> bool:
        """Whether the schedule repeats (MATLAB/JAR accessor name)."""
        return self._cyclic

    def getNumSegments(self) -> int:
        return int(self._rates.size)

    def getPeriod(self) -> float:
        """Horizon length, which is the period when cyclic."""
        return float(self._breakpoints[-1] - self._breakpoints[0])

    def getTimeAverageRate(self) -> float:
        """sum(rates*widths)/sum(widths) over the horizon."""
        return float(np.sum(self._rates * np.diff(self._breakpoints))) / self.getPeriod()

    def getRateAt(self, t: float) -> float:
        """Rate in force at t; zero past a non-cyclic horizon."""
        period = self.getPeriod()
        offset = float(t) - self._breakpoints[0]
        if self._cyclic:
            offset = offset % period
        elif offset < 0.0 or offset >= period:
            return 0.0
        pos = self._breakpoints[0] + offset
        idx = int(np.searchsorted(self._breakpoints[1:], pos, side='right'))
        idx = min(idx, self._rates.size - 1)
        return float(self._rates[idx])

    def getRateSchedule(self) -> dict:
        """The parameterisation of the process; the scalar summaries are not.

        Model compilation recognises a schedule-bearing process by this method
        rather than by class name.
        """
        return {'breakpoints': self._breakpoints, 'rates': self._rates,
                'cyclic': self._cyclic}

    def getMean(self) -> float:
        """Arrival-stationary (Palm) mean interval."""
        return 1.0 / self.getTimeAverageRate()

    def getRate(self) -> float:
        return self.getTimeAverageRate()

    def getVar(self) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def getSCV(self) -> float:
        """NaN: an NHPP is not a renewal process, so there is no i.i.d. interval
        distribution for an SCV to summarise. Returning a representative value
        would report a time-varying process as an exponential one to every
        consumer of sn.scv."""
        return float('nan')

    def getSkew(self) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def getSkewness(self) -> float:
        """NaN; see getSCV (MATLAB/JAR accessor name)."""
        return float('nan')

    def evalCDF(self, x: float) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def evalLST(self, s: float) -> float:
        """NaN: no i.i.d. interval distribution, so no Laplace-Stieltjes
        transform. Overrides the base numerical quadrature, which would
        integrate against an undefined CDF."""
        return float('nan')

    def __repr__(self) -> str:
        return "line_solver.NHPP(%d segments, %s, avgRate=%f)" % (
            self.getNumSegments(),
            "cyclic" if self._cyclic else "non-cyclic",
            self.getTimeAverageRate())

    def getProcess(self):
        return [self._breakpoints, self._rates, self._cyclic]

    def resetSampleClock(self) -> None:
        """Restart the sample path at the schedule start."""
        self._sample_clock = float(self._breakpoints[0])

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Draw n successive interarrival times along ONE sample path.

        The intensity depends on absolute time, so this advances an internal
        clock across calls: consecutive samples form a realisation of the process
        starting at breakpoints[0], not independent draws from a marginal. Use
        resetSampleClock() to restart. A non-cyclic schedule that runs out
        returns 0 for every remaining sample, the intensity there being zero.
        """
        if rng is None:
            rng = np.random.default_rng()
        out = np.zeros(n, dtype=float)
        for i in range(n):
            residual = float(rng.exponential(1.0))
            interval = self.nextInterval(self._sample_clock, residual)
            out[i] = interval
            if interval <= 0.0:
                break  # horizon exhausted: no further event can occur
            self._sample_clock += interval
        return out

    def nextInterval(self, frm: float, residual: float) -> float:
        """Solve int_{frm}^{frm+x} lambda(u) du = residual for x.

        Walks the schedule forward, consuming the budget segment by segment.
        Returns 0 when a non-cyclic horizon is exhausted first, which callers
        read as "no further event".

        Exact for an NHPP: conditional on no event since the last one, the
        residual is governed by the intensity from the current instant onward, so
        a holding time drawn under a rate that has since changed is not a sample
        from this process.
        """
        period = self.getPeriod()
        offset = float(frm) - self._breakpoints[0]
        if self._cyclic:
            offset = offset % period
        elif offset >= period:
            return 0.0
        elif offset < 0.0:
            offset = 0.0
        pos = self._breakpoints[0] + offset
        idx = 0
        while idx < self._rates.size - 1 and pos >= self._breakpoints[idx + 1]:
            idx += 1
        elapsed = 0.0
        while True:
            remaining = self._breakpoints[idx + 1] - pos
            mass = self._rates[idx] * remaining
            # The rate guard also keeps a zero-rate segment from dividing 0/0 on
            # the measure-zero draw residual == 0.
            if self._rates[idx] > 0.0 and mass >= residual:
                return float(elapsed + residual / self._rates[idx])
            residual -= mass
            elapsed += remaining
            idx += 1
            if idx >= self._rates.size:
                if not self._cyclic:
                    return 0.0
                idx = 0
            pos = self._breakpoints[idx]


def _schedule_support(M: np.ndarray, ignore_diagonal: bool) -> np.ndarray:
    """Boolean support pattern of M, optionally excluding the diagonal."""
    pattern = M != 0.0
    if ignore_diagonal:
        pattern = pattern.copy()
        np.fill_diagonal(pattern, False)
    return pattern


def _check_common_support(mats: List[np.ndarray], ignore_diagonal: bool,
                          cls: str, label: str) -> None:
    """Reject a schedule whose matrices do not share one sparsity pattern.

    The fluid solver expresses a time-varying process as a per-entry multiplier
    on a nominal (time-averaged) matrix, and that multiplier is undefined where
    the nominal entry is zero. Requiring one support pattern across segments is
    what makes the nominal nonzero wherever any segment is. A process whose
    phase-transition topology changes in time is therefore refused outright
    rather than silently losing the transitions absent from the nominal.
    """
    ref = _schedule_support(mats[0], ignore_diagonal)
    for k in range(1, len(mats)):
        if not np.array_equal(_schedule_support(mats[k], ignore_diagonal), ref):
            raise ValueError(
                "%s: the %s sparsity pattern must be identical across segments; "
                "segment %d differs from segment 1. A schedule that switches a "
                "transition on or off cannot be expressed as a per-entry "
                "multiplier on the time-averaged process. Keep the entry present "
                "with a small positive rate instead." % (cls, label, k + 1))


class MAPt(ContinuousDistribution):
    """
    Time-inhomogeneous Markovian arrival process (MAP_t).

    Following Ko and Pender (Oper. Res. Lett. 45, 2017), a MAP_t is an ordinary
    MAP whose two matrices are functions of the wall clock, D0(t) and D1(t),
    required only to be locally integrable. This class realises that definition
    with a piecewise-constant schedule, which is dense in L1_loc and is the form
    that serialises: segment k covers [breakpoints[k], breakpoints[k+1]) and
    carries the pair (D0[k], D1[k]), so breakpoints has one more entry than the
    matrix lists. D0 holds transition rates without an arrival, D1 the rates
    that generate one, and D0+D1 is a generator in every segment.

    Two horizon conventions, as for NHPP:
      cyclic     : the schedule repeats with period
                   T = breakpoints[-1] - breakpoints[0].
      non-cyclic : outside [breakpoints[0], breakpoints[-1]) the process is
                   frozen in its last phase and emits nothing, so a non-cyclic
                   MAP_t is a transient construct.

    Setting h = 1 with D0 = [[-lambda_k]], D1 = [[lambda_k]] recovers exactly
    the NHPP with the same breakpoints and rates.

    This is neither a renewal process nor a time-homogeneous one, so the scalar
    summaries that presuppose an i.i.d. interval distribution -- getSCV, getVar,
    getSkew, evalCDF, evalLST -- are undefined and return NaN rather than a
    representative value that would misreport the process as stationary. The
    schedule is the parameterisation: read it with getRateSchedule.

    The class deliberately does NOT extend Markovian. Code gated on
    isMarkovian() reads getProcess() as a single stationary (D0, D1) pair and
    would silently drop the schedule; NHPP avoids the base class for the same
    reason.

    Nominal process. getTimeAverageProcess returns the width-weighted average
    pair (D0bar, D1bar), which is what the fluid solver substitutes as the
    stationary carrier of the phase structure and what sn.rates summarises via
    its MAP arrival rate. For h = 1 that rate coincides with the NHPP
    time-average rate.

    Args:
        breakpoints: strictly increasing segment boundaries, length n+1.
        D0: list of n square matrices of no-arrival transition rates.
        D1: list of n square matrices of arrival-generating rates.
        cyclic: whether the schedule repeats with the horizon as period.
    """

    def __init__(self, breakpoints, D0, D1, cyclic: bool = True):
        super().__init__()
        self._name = 'MAPt'
        breakpoints = np.asarray(breakpoints, dtype=float).ravel()
        if isinstance(D0, np.ndarray) and D0.ndim == 2:
            D0 = [D0]
        if isinstance(D1, np.ndarray) and D1.ndim == 2:
            D1 = [D1]
        D0 = [np.atleast_2d(np.asarray(M, dtype=float)) for M in D0]
        D1 = [np.atleast_2d(np.asarray(M, dtype=float)) for M in D1]
        n = len(D0)
        if n == 0 or len(D1) != n:
            raise ValueError("MAPt: D0 and D1 must be non-empty lists of equal length")
        if breakpoints.size != n + 1:
            raise ValueError(
                "MAPt: breakpoints must have one more entry than the number of segments")
        if np.any(np.diff(breakpoints) <= 0):
            raise ValueError("MAPt: breakpoints must be strictly increasing")
        h = D0[0].shape[0]
        for k in range(n):
            if D0[k].shape != (h, h) or D1[k].shape != (h, h):
                raise ValueError(
                    "MAPt: every D0 and D1 must be square of the same order; "
                    "segment %d has shapes %s and %s against order %d"
                    % (k + 1, D0[k].shape, D1[k].shape, h))
            if np.any(D1[k] < 0.0):
                raise ValueError("MAPt: D1 must be non-negative in segment %d" % (k + 1))
            off = D0[k] - np.diag(np.diag(D0[k]))
            if np.any(off < 0.0):
                raise ValueError(
                    "MAPt: off-diagonal D0 entries must be non-negative in segment %d"
                    % (k + 1))
            if not np.allclose((D0[k] + D1[k]).sum(axis=1), 0.0, atol=1e-10):
                raise ValueError(
                    "MAPt: D0+D1 must have zero row sums (generator) in segment %d"
                    % (k + 1))
        _check_common_support(D0, True, 'MAPt', 'off-diagonal D0')
        _check_common_support(D1, False, 'MAPt', 'D1')
        if all(float(np.sum(D1[k])) <= 0.0 for k in range(n)):
            raise ValueError(
                "MAPt: every segment has zero arrival intensity, so no event can ever occur")
        self._breakpoints = breakpoints
        self._D0 = D0
        self._D1 = D1
        self._cyclic = bool(cyclic)
        # Wall-clock position and phase of the next sample; see sample().
        self._sample_clock = float(breakpoints[0])
        self._sample_phase = 0

    @property
    def breakpoints(self) -> np.ndarray:
        """Segment boundaries, length n+1."""
        return self._breakpoints

    @property
    def D0(self) -> List[np.ndarray]:
        """Per-segment no-arrival rate matrices."""
        return [M.copy() for M in self._D0]

    @property
    def D1(self) -> List[np.ndarray]:
        """Per-segment arrival-generating rate matrices."""
        return [M.copy() for M in self._D1]

    @property
    def cyclic(self) -> bool:
        """Whether the schedule repeats."""
        return self._cyclic

    def getBreakpoints(self) -> np.ndarray:
        """Segment boundaries, length n+1 (MATLAB/JAR accessor name)."""
        return self._breakpoints

    def getD0Segments(self) -> List[np.ndarray]:
        """Per-segment D0 matrices (MATLAB/JAR accessor name)."""
        return [M.copy() for M in self._D0]

    def getD1Segments(self) -> List[np.ndarray]:
        """Per-segment D1 matrices (MATLAB/JAR accessor name)."""
        return [M.copy() for M in self._D1]

    def isCyclic(self) -> bool:
        """Whether the schedule repeats (MATLAB/JAR accessor name)."""
        return self._cyclic

    def getNumSegments(self) -> int:
        return len(self._D0)

    def getNumberOfPhases(self) -> int:
        return int(self._D0[0].shape[0])

    def getPeriod(self) -> float:
        """Horizon length, which is the period when cyclic."""
        return float(self._breakpoints[-1] - self._breakpoints[0])

    def getSegmentIndexAt(self, t: float) -> int:
        """Index of the segment in force at t, or -1 past a non-cyclic horizon."""
        period = self.getPeriod()
        offset = float(t) - self._breakpoints[0]
        if self._cyclic:
            offset = offset % period
        elif offset < 0.0 or offset >= period:
            return -1
        pos = self._breakpoints[0] + offset
        idx = int(np.searchsorted(self._breakpoints[1:], pos, side='right'))
        return min(idx, len(self._D0) - 1)

    def getD0At(self, t: float) -> np.ndarray:
        """D0 in force at t; the zero matrix past a non-cyclic horizon."""
        idx = self.getSegmentIndexAt(t)
        if idx < 0:
            return np.zeros_like(self._D0[0])
        return self._D0[idx].copy()

    def getD1At(self, t: float) -> np.ndarray:
        """D1 in force at t; the zero matrix past a non-cyclic horizon."""
        idx = self.getSegmentIndexAt(t)
        if idx < 0:
            return np.zeros_like(self._D1[0])
        return self._D1[idx].copy()

    def getTimeAverageProcess(self) -> Tuple[np.ndarray, np.ndarray]:
        """Width-weighted average (D0bar, D1bar) over the horizon.

        This is the nominal stationary MAP that carries the phase structure
        where a solver needs a time-homogeneous carrier. It is a valid MAP: a
        convex combination of generators is a generator, and non-negativity is
        preserved entrywise.
        """
        widths = np.diff(self._breakpoints)
        total = float(np.sum(widths))
        D0bar = sum(w * M for w, M in zip(widths, self._D0)) / total
        D1bar = sum(w * M for w, M in zip(widths, self._D1)) / total
        return D0bar, D1bar

    def getTimeAverageRate(self) -> float:
        """Arrival rate of the time-averaged MAP, i.e. pi*D1bar*e.

        For h = 1 this is exactly the NHPP width-weighted average intensity.
        """
        from ..api.mam import map_lambda
        D0bar, D1bar = self.getTimeAverageProcess()
        return float(map_lambda(D0bar, D1bar))

    def getRateAt(self, t: float) -> float:
        """Arrival rate of the MAP in force at t; zero past a non-cyclic horizon.

        This is the stationary rate of that segment's MAP, not the instantaneous
        conditional intensity, which depends on the current phase.
        """
        from ..api.mam import map_lambda
        idx = self.getSegmentIndexAt(t)
        if idx < 0:
            return 0.0
        return float(map_lambda(self._D0[idx], self._D1[idx]))

    def getRateSchedule(self) -> dict:
        """The parameterisation of the process; the scalar summaries are not.

        Model compilation recognises a schedule-bearing process by this method
        rather than by class name.
        """
        return {'breakpoints': self._breakpoints,
                'D0': [M.copy() for M in self._D0],
                'D1': [M.copy() for M in self._D1],
                'cyclic': self._cyclic}

    def getMean(self) -> float:
        """Arrival-stationary (Palm) mean interval of the time-averaged MAP."""
        return 1.0 / self.getTimeAverageRate()

    def getRate(self) -> float:
        return self.getTimeAverageRate()

    def getVar(self) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def getSCV(self) -> float:
        """NaN: a MAP_t is neither renewal nor time-homogeneous, so there is no
        i.i.d. interval distribution for an SCV to summarise. Returning the SCV
        of the time-averaged MAP would report a time-varying process as a
        stationary one to every consumer of sn.scv."""
        return float('nan')

    def getSkew(self) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def getSkewness(self) -> float:
        """NaN; see getSCV (MATLAB/JAR accessor name)."""
        return float('nan')

    def evalCDF(self, x: float) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def evalLST(self, s: float) -> float:
        """NaN: no i.i.d. interval distribution, so no Laplace-Stieltjes
        transform. Overrides the base numerical quadrature, which would
        integrate against an undefined CDF."""
        return float('nan')

    def __repr__(self) -> str:
        return "line_solver.MAPt(%d segments, %d phases, %s, avgRate=%f)" % (
            self.getNumSegments(), self.getNumberOfPhases(),
            "cyclic" if self._cyclic else "non-cyclic",
            self.getTimeAverageRate())

    def getProcess(self):
        return [self._breakpoints, [M.copy() for M in self._D0],
                [M.copy() for M in self._D1], self._cyclic]

    def resetSampleClock(self) -> None:
        """Restart the sample path at the schedule start, in phase 1."""
        self._sample_clock = float(self._breakpoints[0])
        self._sample_phase = 0

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Draw n successive interarrival times along ONE sample path.

        Both the intensity and the phase depend on absolute time, so this
        advances an internal clock and phase across calls: consecutive samples
        form a realisation of the process starting at breakpoints[0] in phase 1,
        not independent draws from a marginal. Use resetSampleClock() to
        restart. A non-cyclic schedule that runs out returns 0 for every
        remaining sample.
        """
        if rng is None:
            rng = np.random.default_rng()
        out = np.zeros(n, dtype=float)
        for i in range(n):
            interval, phase = self.nextArrival(self._sample_clock, self._sample_phase, rng)
            out[i] = interval
            if interval <= 0.0:
                break  # horizon exhausted: no further event can occur
            self._sample_clock += interval
            self._sample_phase = phase
        return out

    def nextArrival(self, frm: float, phase: int,
                    rng: Optional[np.random.Generator] = None) -> Tuple[float, int]:
        """Time to the next arrival from wall clock frm in the given phase.

        Returns (interval, phase after the arrival). Exact: within a segment the
        phase process is a homogeneous CTMC, and by the memoryless property the
        residual holding time may be redrawn at a breakpoint without biasing the
        path, so the boundary is crossed by advancing the clock and resampling
        under the new matrices. Returns (0, phase) when a non-cyclic horizon is
        exhausted, which callers read as "no further arrival".
        """
        if rng is None:
            rng = np.random.default_rng()
        elapsed = 0.0
        pos = float(frm)
        while True:
            idx = self.getSegmentIndexAt(pos)
            if idx < 0:
                return 0.0, phase
            # Time left in the active segment, unrolling a cyclic schedule.
            period = self.getPeriod()
            offset = pos - self._breakpoints[0]
            if self._cyclic:
                offset = offset % period
            seg_end_offset = self._breakpoints[idx + 1] - self._breakpoints[0]
            to_boundary = seg_end_offset - offset
            D0 = self._D0[idx]
            D1 = self._D1[idx]
            total = -float(D0[phase, phase])
            if total <= 0.0:
                # Absorbing phase in this segment: only a boundary can free it.
                if not self._cyclic and idx == len(self._D0) - 1:
                    return 0.0, phase
                elapsed += to_boundary
                pos += to_boundary
                continue
            holding = float(rng.exponential(1.0 / total))
            if holding >= to_boundary:
                if not self._cyclic and idx == len(self._D0) - 1:
                    return 0.0, phase
                elapsed += to_boundary
                pos += to_boundary
                continue
            elapsed += holding
            pos += holding
            # Competing transitions out of the current phase, arrivals first.
            h = D0.shape[0]
            weights = np.concatenate([D1[phase, :], D0[phase, :].copy()])
            weights[h + phase] = 0.0
            u = float(rng.random()) * total
            cum = 0.0
            for j, wgt in enumerate(weights):
                cum += wgt
                if u < cum:
                    if j < h:
                        return elapsed, j
                    phase = j - h
                    break
            else:
                # Rounding shortfall: attribute the draw to the last positive entry.
                j = int(np.max(np.nonzero(weights)[0]))
                if j < h:
                    return elapsed, j
                phase = j - h


class PHt(ContinuousDistribution):
    """
    Time-inhomogeneous phase-type distribution (Ph_t).

    Following Ko and Pender (Oper. Res. Lett. 45, 2017), a Ph_t is an ordinary
    phase-type distribution whose initial vector and sub-generator are functions
    of the wall clock, alpha(t) and S(t), required only to be locally
    integrable. This class realises that definition with a piecewise-constant
    schedule: segment k covers [breakpoints[k], breakpoints[k+1]) and carries
    the pair (alpha[k], S[k]), so breakpoints has one more entry than the lists.
    The exit vector is s(t) = -S(t)e.

    Because both the phase and the elapsed service depend on absolute time, a
    Ph_t service time is a function of the epoch at which service starts:
    sampleFrom(t0) is the operative sampler, and sample() walks one path.

    Setting h = 1 with S = [[-mu_k]] recovers a time-varying exponential, whose
    completion stream at a saturated server is the NHPP with rates mu_k.

    Like MAPt this does NOT extend Markovian, so that isMarkovian()-gated code
    cannot read it as a single stationary (alpha, S) pair; and the scalar
    summaries getSCV, getVar, getSkew, evalCDF, evalLST return NaN, the
    distribution of a service time being different at every start epoch.

    Args:
        breakpoints: strictly increasing segment boundaries, length n+1.
        alpha: list of n initial probability row vectors, each summing to 1.
        S: list of n sub-generator matrices with non-positive row sums.
        cyclic: whether the schedule repeats with the horizon as period.
    """

    def __init__(self, breakpoints, alpha, S, cyclic: bool = True):
        super().__init__()
        self._name = 'PHt'
        breakpoints = np.asarray(breakpoints, dtype=float).ravel()
        if isinstance(S, np.ndarray) and S.ndim == 2:
            S = [S]
        if isinstance(alpha, np.ndarray) and alpha.ndim == 1:
            alpha = [alpha]
        alpha = [np.asarray(a, dtype=float).ravel() for a in alpha]
        S = [np.atleast_2d(np.asarray(M, dtype=float)) for M in S]
        n = len(S)
        if n == 0 or len(alpha) != n:
            raise ValueError("PHt: alpha and S must be non-empty lists of equal length")
        if breakpoints.size != n + 1:
            raise ValueError(
                "PHt: breakpoints must have one more entry than the number of segments")
        if np.any(np.diff(breakpoints) <= 0):
            raise ValueError("PHt: breakpoints must be strictly increasing")
        h = S[0].shape[0]
        for k in range(n):
            if S[k].shape != (h, h) or alpha[k].size != h:
                raise ValueError(
                    "PHt: every S must be square of order %d with a matching alpha; "
                    "segment %d has shapes %s and %s" % (h, k + 1, S[k].shape, alpha[k].shape))
            if np.any(alpha[k] < 0.0) or not np.isclose(float(np.sum(alpha[k])), 1.0, atol=1e-10):
                raise ValueError(
                    "PHt: alpha must be a probability vector in segment %d" % (k + 1))
            off = S[k] - np.diag(np.diag(S[k]))
            if np.any(off < 0.0):
                raise ValueError(
                    "PHt: off-diagonal S entries must be non-negative in segment %d" % (k + 1))
            exit_rates = -S[k].sum(axis=1)
            if np.any(exit_rates < -1e-10):
                raise ValueError(
                    "PHt: S must have non-positive row sums in segment %d" % (k + 1))
        _check_common_support(S, True, 'PHt', 'off-diagonal S')
        _check_common_support([(-M.sum(axis=1)).reshape(-1, 1) for M in S], False,
                              'PHt', 'exit vector')
        _check_common_support([a.reshape(1, -1) for a in alpha], False, 'PHt', 'alpha')
        if all(float(np.sum(-S[k].sum(axis=1))) <= 0.0 for k in range(n)):
            raise ValueError(
                "PHt: every segment has zero exit rate, so service can never complete")
        self._breakpoints = breakpoints
        self._alpha = alpha
        self._S = S
        self._cyclic = bool(cyclic)
        self._sample_clock = float(breakpoints[0])

    @property
    def breakpoints(self) -> np.ndarray:
        """Segment boundaries, length n+1."""
        return self._breakpoints

    @property
    def alpha(self) -> List[np.ndarray]:
        """Per-segment initial probability vectors."""
        return [a.copy() for a in self._alpha]

    @property
    def S(self) -> List[np.ndarray]:
        """Per-segment sub-generators."""
        return [M.copy() for M in self._S]

    @property
    def cyclic(self) -> bool:
        """Whether the schedule repeats."""
        return self._cyclic

    def getBreakpoints(self) -> np.ndarray:
        """Segment boundaries, length n+1 (MATLAB/JAR accessor name)."""
        return self._breakpoints

    def getAlphaSegments(self) -> List[np.ndarray]:
        """Per-segment initial vectors (MATLAB/JAR accessor name)."""
        return [a.copy() for a in self._alpha]

    def getSSegments(self) -> List[np.ndarray]:
        """Per-segment sub-generators (MATLAB/JAR accessor name)."""
        return [M.copy() for M in self._S]

    def isCyclic(self) -> bool:
        """Whether the schedule repeats (MATLAB/JAR accessor name)."""
        return self._cyclic

    def getNumSegments(self) -> int:
        return len(self._S)

    def getNumberOfPhases(self) -> int:
        return int(self._S[0].shape[0])

    def getPeriod(self) -> float:
        """Horizon length, which is the period when cyclic."""
        return float(self._breakpoints[-1] - self._breakpoints[0])

    def getSegmentIndexAt(self, t: float) -> int:
        """Index of the segment in force at t, or -1 past a non-cyclic horizon."""
        period = self.getPeriod()
        offset = float(t) - self._breakpoints[0]
        if self._cyclic:
            offset = offset % period
        elif offset < 0.0 or offset >= period:
            return -1
        pos = self._breakpoints[0] + offset
        idx = int(np.searchsorted(self._breakpoints[1:], pos, side='right'))
        return min(idx, len(self._S) - 1)

    def getAlphaAt(self, t: float) -> np.ndarray:
        """alpha in force at t; the last segment's vector past a non-cyclic
        horizon, where getSAt is zero so no service can complete anyway."""
        idx = self.getSegmentIndexAt(t)
        if idx < 0:
            return self._alpha[-1].copy()
        return self._alpha[idx].copy()

    def getSAt(self, t: float) -> np.ndarray:
        """S in force at t; the zero matrix past a non-cyclic horizon."""
        idx = self.getSegmentIndexAt(t)
        if idx < 0:
            return np.zeros_like(self._S[0])
        return self._S[idx].copy()

    def getTimeAverageProcess(self) -> Tuple[np.ndarray, np.ndarray]:
        """Width-weighted average (alphabar, Sbar) over the horizon.

        A convex combination of sub-generators is a sub-generator and of
        probability vectors a probability vector, so the nominal is a valid
        phase-type representation.
        """
        widths = np.diff(self._breakpoints)
        total = float(np.sum(widths))
        abar = sum(w * a for w, a in zip(widths, self._alpha)) / total
        Sbar = sum(w * M for w, M in zip(widths, self._S)) / total
        return abar, Sbar

    def getTimeAverageProcessMAP(self) -> Tuple[np.ndarray, np.ndarray]:
        """The nominal as a (D0, D1) pair, D1 = s*alpha, for the fluid carrier."""
        abar, Sbar = self.getTimeAverageProcess()
        sbar = -Sbar.sum(axis=1).reshape(-1, 1)
        return Sbar, sbar @ abar.reshape(1, -1)

    def getTimeAverageRate(self) -> float:
        """Completion rate of the time-averaged phase-type, 1/(-alphabar*Sbar^-1*e)."""
        abar, Sbar = self.getTimeAverageProcess()
        return 1.0 / float(-abar @ np.linalg.solve(Sbar, np.ones(Sbar.shape[0])))

    def getRateAt(self, t: float) -> float:
        """Completion rate of the phase-type in force at t; zero past a
        non-cyclic horizon."""
        idx = self.getSegmentIndexAt(t)
        if idx < 0:
            return 0.0
        S = self._S[idx]
        return 1.0 / float(-self._alpha[idx] @ np.linalg.solve(S, np.ones(S.shape[0])))

    def getRateSchedule(self) -> dict:
        """The parameterisation of the process; the scalar summaries are not.

        Model compilation recognises a schedule-bearing process by this method
        rather than by class name.
        """
        return {'breakpoints': self._breakpoints,
                'alpha': [a.copy() for a in self._alpha],
                'S': [M.copy() for M in self._S],
                'cyclic': self._cyclic}

    def getMean(self) -> float:
        """Mean of the time-averaged phase-type."""
        return 1.0 / self.getTimeAverageRate()

    def getRate(self) -> float:
        return self.getTimeAverageRate()

    def getVar(self) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def getSCV(self) -> float:
        """NaN: the service-time distribution differs at every start epoch, so
        there is no single i.i.d. law for an SCV to summarise. Returning the SCV
        of the time-averaged representation would report a time-varying process
        as a stationary one to every consumer of sn.scv."""
        return float('nan')

    def getSkew(self) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def getSkewness(self) -> float:
        """NaN; see getSCV (MATLAB/JAR accessor name)."""
        return float('nan')

    def evalCDF(self, x: float) -> float:
        """NaN; see getSCV."""
        return float('nan')

    def evalLST(self, s: float) -> float:
        """NaN: no single interval distribution, so no Laplace-Stieltjes
        transform. Overrides the base numerical quadrature, which would
        integrate against an undefined CDF."""
        return float('nan')

    def __repr__(self) -> str:
        return "line_solver.PHt(%d segments, %d phases, %s, avgRate=%f)" % (
            self.getNumSegments(), self.getNumberOfPhases(),
            "cyclic" if self._cyclic else "non-cyclic",
            self.getTimeAverageRate())

    def getProcess(self):
        return [self._breakpoints, [a.copy() for a in self._alpha],
                [M.copy() for M in self._S], self._cyclic]

    def resetSampleClock(self) -> None:
        """Restart the sample path at the schedule start."""
        self._sample_clock = float(self._breakpoints[0])

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Draw n successive service times along ONE sample path.

        The law depends on absolute time, so this advances an internal clock
        across calls: sample i starts where sample i-1 completed, not at a fixed
        epoch. Use resetSampleClock() to restart, or sampleFrom() to draw a
        service time starting at a chosen epoch. A non-cyclic schedule that runs
        out returns 0 for every remaining sample.
        """
        if rng is None:
            rng = np.random.default_rng()
        out = np.zeros(n, dtype=float)
        for i in range(n):
            interval = self.sampleFrom(self._sample_clock, rng)
            out[i] = interval
            if interval <= 0.0:
                break  # horizon exhausted: service can never complete
            self._sample_clock += interval
        return out

    def sampleFrom(self, t0: float, rng: Optional[np.random.Generator] = None) -> float:
        """Service time for a job whose service starts at wall clock t0.

        Exact: within a segment the phase process is a homogeneous absorbing
        CTMC, and by the memoryless property the residual holding time may be
        redrawn at a breakpoint, so the boundary is crossed by advancing the
        clock and resampling under the new sub-generator. The initial phase is
        drawn from alpha in force at t0. Returns 0 when a non-cyclic horizon is
        exhausted before absorption.
        """
        if rng is None:
            rng = np.random.default_rng()
        idx = self.getSegmentIndexAt(t0)
        if idx < 0:
            return 0.0
        alpha = self._alpha[idx]
        phase = int(np.searchsorted(np.cumsum(alpha), float(rng.random()) * float(np.sum(alpha))))
        phase = min(phase, alpha.size - 1)
        elapsed = 0.0
        pos = float(t0)
        while True:
            idx = self.getSegmentIndexAt(pos)
            if idx < 0:
                return 0.0
            period = self.getPeriod()
            offset = pos - self._breakpoints[0]
            if self._cyclic:
                offset = offset % period
            to_boundary = (self._breakpoints[idx + 1] - self._breakpoints[0]) - offset
            S = self._S[idx]
            total = -float(S[phase, phase])
            if total <= 0.0:
                if not self._cyclic and idx == len(self._S) - 1:
                    return 0.0
                elapsed += to_boundary
                pos += to_boundary
                continue
            holding = float(rng.exponential(1.0 / total))
            if holding >= to_boundary:
                if not self._cyclic and idx == len(self._S) - 1:
                    return 0.0
                elapsed += to_boundary
                pos += to_boundary
                continue
            elapsed += holding
            pos += holding
            # Competing transitions: absorption first, then phase changes.
            exit_rate = float(-S[phase, :].sum())
            weights = np.concatenate([[exit_rate], S[phase, :].copy()])
            weights[1 + phase] = 0.0
            u = float(rng.random()) * total
            cum = 0.0
            for j, wgt in enumerate(weights):
                cum += wgt
                if u < cum:
                    if j == 0:
                        return elapsed
                    phase = j - 1
                    break
            else:
                j = int(np.max(np.nonzero(weights)[0]))
                if j == 0:
                    return elapsed
                phase = j - 1
