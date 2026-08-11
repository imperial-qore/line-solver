"""
Markovian distributions for LINE (pure Python).

This module provides Markovian (phase-type) distribution implementations
including PH, APH, Coxian, and MAP (Markovian Arrival Process).
"""

from typing import Optional, Tuple, Union
import math
import numpy as np
from scipy import linalg

from .base import ContinuousDistribution, Markovian
from ..api.io.logging import line_warning

# see _kb/01-model-classes.md (ME/RAP section) for rationale
_ME_CHECK_PREC = 1e-14

# see _kb/01-model-classes.md (ME/RAP section) for rationale
_ME_SCAN_TAIL = 1e-12       # residual mass left beyond the scan horizon
_ME_SCAN_HORIZON_CAP = 1e4  # hard cap on the horizon for near-degenerate A
_ME_SCAN_PER_PERIOD = 20    # samples per period of the fastest oscillation
_ME_SCAN_MIN_PTS = 2001
_ME_SCAN_MAX_PTS = 200001
_ME_SCAN_RELTOL = 1e-10     # negative only if below -reltol * max|f|

# Inverse-CDF sampling parameters for ME.
_ME_TAIL_TOL = 1e-12        # survival value at which the CDF table is complete
_ME_GRID_STEPS = 2000       # grid points over the initial [0, mean + 10*sigma]
_ME_MAX_STEPS = 400000      # hard cap on table length
_ME_NEWTON_ITERS = 12
_ME_NEWTON_TOL = 1e-13

# Conditional-survival inversion parameters for RAP sampling.
_RAP_BRACKET_ITERS = 200
_RAP_ROOT_ITERS = 200
_RAP_ROOT_TOL = 1e-13


def _butools_settings():
    """Return the BuTools settings module (carries checkInput/verbose flags)."""
    from ..lib.thirdparty import butools
    return butools


def _check_me_representation(alpha_row, A):
    """BuTools CheckMERepresentation on a (1,n) initial vector and (n,n) matrix."""
    from ..lib.thirdparty.butools.ph.check import CheckMERepresentation
    return bool(CheckMERepresentation(np.matrix(alpha_row), np.matrix(A),
                                      _ME_CHECK_PREC))


def _scan_negative_density(alpha_row, A):
    """Search the density f(t) = -alpha*expm(A*t)*A*e for a negative value.

    Returns ``(is_negative, f_min, t_min)``. A negative value found here is a
    witness: it proves the representation is not a distribution. Finding none
    proves nothing, so the caller must not report the converse.

    This replaces BuTools CheckMEPositiveDensity as the trigger for the
    construction-time warning. That routine searches for a Markovian monocyclic
    equivalent, which is a sufficient condition only, and its verdict depends on
    the representation rather than on the distribution: for
    ``alpha = [1,0,0]``, ``A = [[-0.5,0,0],[0,-1,w],[0,-w,-1]]`` the
    distribution is Exp(0.5) for every w, yet the search fails once w >= 2*pi.
    It also costs of the order of a second per call at search order 1000, which
    is far too slow for a constructor.

    The horizon covers all but _ME_SCAN_TAIL of the mass, using the dominant
    (least negative) eigenvalue of A; the sampling rate resolves the fastest
    oscillation present, taken from the largest imaginary part.
    """
    alpha = np.asarray(alpha_row, dtype=float).reshape(1, -1)
    A = np.asarray(A, dtype=float)
    n = A.shape[0]

    lam = np.linalg.eigvals(A)
    decay = np.max(lam.real)
    if not np.isfinite(decay) or decay >= 0.0:
        # Not a valid ME; CheckMERepresentation has already rejected it.
        return False, 0.0, 0.0
    horizon = min(-math.log(_ME_SCAN_TAIL) / abs(decay), _ME_SCAN_HORIZON_CAP)

    npts = _ME_SCAN_MIN_PTS
    wmax = float(np.max(np.abs(lam.imag)))
    if wmax > 0.0:
        npts = max(npts, int(math.ceil(horizon * wmax * _ME_SCAN_PER_PERIOD
                                       / (2.0 * math.pi))) + 1)
    npts = min(npts, _ME_SCAN_MAX_PTS)

    step = horizon / (npts - 1)
    # One matrix exponential, then propagate: v_k = alpha*expm(A*k*step).
    E = linalg.expm(A * step)
    ve = -(A @ np.ones((n, 1)))

    v = alpha
    f_min = math.inf
    t_min = 0.0
    f_absmax = 0.0
    for k in range(npts):
        f = float(v @ ve)
        if f < f_min:
            f_min = f
            t_min = k * step
        if abs(f) > f_absmax:
            f_absmax = abs(f)
        v = v @ E

    return bool(f_min < -_ME_SCAN_RELTOL * f_absmax), f_min, t_min


def _check_rap_representation(H0, H1):
    """BuTools CheckRAPRepresentation on the (H0,H1) pair."""
    from ..lib.thirdparty.butools.map.check import CheckRAPRepresentation
    return bool(CheckRAPRepresentation(np.matrix(H0), np.matrix(H1),
                                       _ME_CHECK_PREC))


def _spectral_expm_data(A):
    """Diagonalize A so that v*exp(A*t)*w can be evaluated in closed form.

    Returns ``(lam, V, Vinv)`` with ``A = V diag(lam) Vinv``, or ``None`` when A
    is defective or the eigenvector basis is too ill-conditioned to be used, in
    which case callers fall back to ``scipy.linalg.expm``.
    """
    try:
        lam, V = linalg.eig(np.asarray(A, dtype=float))
        cond = np.linalg.cond(V)
        if not np.isfinite(cond) or cond > 1e10:
            return None
        Vinv = linalg.inv(V)
    except (linalg.LinAlgError, ValueError):
        return None
    if not np.allclose(V @ np.diag(lam) @ Vinv, A, atol=1e-9, rtol=1e-9):
        return None
    return lam, V, Vinv


def _rap_row_expm(v, x, H0, spectral):
    """Return the row vector v*exp(H0*x)."""
    if spectral is not None:
        lam, V, Vinv = spectral
        return np.real((v @ V) * np.exp(lam * x) @ Vinv)
    return v @ linalg.expm(H0 * x)


def _rap_survival_and_derivative(v, x, H0, spectral, e, coeff):
    """Return (v*exp(H0*x)*e, v*exp(H0*x)*H0*e) at the scalar point x."""
    if spectral is not None:
        lam, _, _ = spectral
        ex = np.exp(lam * x)
        return float(np.real(coeff @ ex)), float(np.real((coeff * lam) @ ex))
    w = v @ linalg.expm(H0 * x)
    return float(w @ e), float(w @ (H0 @ e))


def _rap_invert_survival(v, target, H0, spectral, e, x_guess):
    """Solve v*exp(H0*x)*e = target for x >= 0.

    The survival function decreases monotonically from v*e = 1, so the root is
    bracketed by doubling and then located by Newton steps that fall back to
    bisection whenever they leave the bracket.
    """
    coeff = None
    if spectral is not None:
        lam, V, Vinv = spectral
        coeff = (v @ V) * (Vinv @ e)

    g0, _ = _rap_survival_and_derivative(v, 0.0, H0, spectral, e, coeff)
    if target >= g0:
        return 0.0

    lo = 0.0
    hi = x_guess if x_guess > 0 else 1.0
    for _ in range(_RAP_BRACKET_ITERS):
        ghi, _ = _rap_survival_and_derivative(v, hi, H0, spectral, e, coeff)
        if ghi <= target:
            break
        lo = hi
        hi *= 2.0
    else:
        raise ValueError(
            "RAP sampling could not bracket the inter-arrival time: the "
            "survival function stays above %g out to x = %g" % (target, hi))

    x = 0.5 * (lo + hi)
    for _ in range(_RAP_ROOT_ITERS):
        g, dg = _rap_survival_and_derivative(v, x, H0, spectral, e, coeff)
        resid = g - target
        if abs(resid) < _RAP_ROOT_TOL:
            break
        if resid > 0:
            lo = x
        else:
            hi = x
        x_new = x - resid / dg if dg < 0 else x
        if not np.isfinite(x_new) or x_new <= lo or x_new >= hi:
            x_new = 0.5 * (lo + hi)
        if abs(x_new - x) < _RAP_ROOT_TOL * max(1.0, abs(x)):
            x = x_new
            break
        x = x_new
    return x


def _dominant_eigenvalue(A):
    """Return the eigenvalue of A with the largest (least negative) real part."""
    ev = linalg.eigvals(np.asarray(A, dtype=float))
    return float(np.max(np.real(ev)))


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

    def getSkew(self) -> float:
        """Get the skewness, (E3 - 3 E1 E2 + 2 E1^3) / (E2 - E1^2)^(3/2).

        The base Distribution returns 0.0, which is silently wrong for a
        phase-type; the moments are available in closed form here.
        """
        try:
            T_inv = linalg.inv(self._T)
            e = np.ones(len(self._alpha))
            m1 = float(-self._alpha @ T_inv @ e)
            m2 = float(2 * self._alpha @ T_inv @ T_inv @ e)
            m3 = float(-6 * self._alpha @ T_inv @ T_inv @ T_inv @ e)
            var = m2 - m1 ** 2
            if var <= 0:
                return 0.0
            return (m3 - 3 * m1 * m2 + 2 * m1 ** 3) / var ** 1.5
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



    @classmethod
    def fit(cls, mean: float, scv: float, skew: float = None) -> 'PH':
        """Fit a PH to (mean, SCV, skewness).

        Mirrors MATLAB PH.fit, whose two-state search falls back to the acyclic
        fitter when it finds nothing; APH is a PH, so the acyclic fit is a
        valid return value here.
        """
        from .markovian import APH as _APH
        return _APH.fit(mean, scv, skew)

    @classmethod
    def fit_central(cls, mean: float, var: float, skew: float = None) -> 'PH':
        """Fit a PH to (mean, variance, skewness), as MATLAB PH.fitCentral."""
        from .markovian import APH as _APH
        return _APH.fit_central(mean, var / (mean * mean), skew)

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'PH':
        """Fit a PH to (mean, SCV), as MATLAB PH.fitMeanAndSCV.

        Two moments only, so the minimum-skewness third moment of the acyclic
        class is used. The MATLAB method referenced an undefined SKEW until
        2026-07-22 and could not be called at all.
        """
        from .markovian import APH as _APH
        return _APH.fit_mean_and_scv(mean, scv)

    @classmethod
    def fit_raw_moments(cls, m1: float, m2: float, m3: float) -> 'PH':
        """Fit a PH to the first three raw moments, as MATLAB
        PH.fitRawMoments."""
        from .markovian import APH as _APH
        return _APH.fit_raw_moments(m1, m2, m3)


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
        # see _kb/01-model-classes.md (ME/RAP section) for rationale
        if skew is not None:
            try:
                from ..lib.thirdparty.butools.ph.canonical import APHFrom3Moments
                e1 = float(mean)
                e2 = (1.0 + float(scv)) * e1 ** 2
                e3 = -(2.0 * e1 ** 3 - 3.0 * e1 * e2 - float(skew) * (e2 - e1 ** 2) ** 1.5)
                alpha, T = APHFrom3Moments([e1, e2, e3])
                return np.asarray(alpha, dtype=float).ravel(), np.asarray(T, dtype=float)
            except Exception:
                pass  # infeasible moments -> fall back to 2-moment matching

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
        # Use mixture of Erlang-k and Erlang-(k-1) via starting phase selection
        # Starting in phase 0 gives Erlang-k, starting in phase 1 gives Erlang-(k-1)
        # SCV of Erlang-k is 1/k, SCV of Erlang-(k-1) is 1/(k-1)
        if k <= 1:
            rate = 1.0 / mean
            return np.array([1.0]), np.array([[-rate]])

        # Check if we can use pure Erlang-k (when scv = 1/k)
        scv_k = 1.0 / k
        if abs(scv - scv_k) < 1e-10:
            # Pure Erlang-k
            rate = k / mean
            alpha = np.zeros(k)
            alpha[0] = 1.0
            T = np.diag([-rate] * k) + np.diag([rate] * (k - 1), 1)
            return alpha, T

        # For SCV between 1/k and 1/(k-1), use mixture
        # p = probability of starting in phase 0 (Erlang-k)
        # 1-p = probability of starting in phase 1 (Erlang-(k-1))
        #
        # Using moment matching:
        # Mean = (p*k + (1-p)*(k-1)) / r = mean  =>  r = (k - 1 + p) / mean
        # For the mixture SCV, solving gives:
        # p = (k - 1) * (1 - scv*(k-1)) / (scv*(k-1) - 1 + 1/k * (k - scv*(k-1)*(k-1)))
        # Simplified: p such that target SCV is achieved

        # Simpler approach: linear interpolation based on SCV
        # SCV = 1/k when p=1, SCV = 1/(k-1) when p=0
        # p = (1/(k-1) - scv) / (1/(k-1) - 1/k) = (1/(k-1) - scv) * k * (k-1) / 1
        scv_k_minus_1 = 1.0 / (k - 1) if k > 1 else 1.0
        p = (scv_k_minus_1 - scv) / (scv_k_minus_1 - scv_k)
        p = max(0, min(1, p))

        # Rate to match mean: mean = (p*k + (1-p)*(k-1)) / rate
        effective_phases = p * k + (1 - p) * (k - 1)
        rate = effective_phases / mean

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
        """Match moments using BUTools APHFrom2Moments algorithm.

        This produces a 2-phase APH with the same structure as MATLAB's
        BUTools APHFrom2Moments, ensuring identical simulation results
        when using the same seed in JMT.

        The structure is:
        - Phase 1 transitions to phase 2 with rate lambda*p*N
        - Phase 2 absorbs with rate lambda*N
        - Initial probability: alpha = [p, 1-p]
        """
        # BUTools APHFrom2Moments algorithm
        # cv2 = moms(2)/moms(1)^2 - 1.0 = scv
        cv2 = scv
        lam = 1.0 / mean  # base rate (lambda)
        N = max(int(np.ceil(1.0 / cv2)), 2)  # number of phases, at least 2

        # Probability of starting in phase 1
        p = 1.0 / (cv2 + 1.0 + (cv2 - 1.0) / (N - 1))

        # Build the generator matrix A
        A = -lam * p * N * np.eye(N)
        for i in range(N - 1):
            A[i, i + 1] = -A[i, i]  # transition from phase i to i+1
        A[N - 1, N - 1] = -lam * N  # last phase rate

        # Initial probability vector
        alpha = np.zeros(N)
        alpha[0] = p
        alpha[N - 1] = 1.0 - p

        return alpha, A

    @classmethod
    def fit(cls, mean: float, scv: float, skew: float = None) -> 'APH':
        """Fit an APH to (mean, SCV, skewness), as MATLAB APH.fit does."""
        return cls.fit_central(mean, scv, skew)

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
    fitMeanAndSCV = fit_mean_and_scv

    @classmethod
    def fit_central(cls, mean: float, scv: float, skew: float = None) -> 'APH':
        """
        Create an APH distribution from central moments.

        Uses moment matching to construct an acyclic phase-type
        distribution with the specified mean, SCV, and optionally skewness.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.
            skew: Target skewness (optional).

        Returns:
            APH distribution with given moments.
        """
        if skew is None:
            return cls(mean, scv)
        # see _kb/01-model-classes.md (ME/RAP section) for rationale
        try:
            from ..lib.thirdparty.butools.ph.canonical import APHFrom3Moments
            e1 = float(mean)
            e2 = (1.0 + float(scv)) * e1 ** 2
            e3 = -(2.0 * e1 ** 3 - 3.0 * e1 * e2 - float(skew) * (e2 - e1 ** 2) ** 1.5)
            alpha, T = APHFrom3Moments([e1, e2, e3])
            return cls(np.asarray(alpha, dtype=float).ravel(), np.asarray(T, dtype=float))
        except Exception:
            return cls(mean, scv, skew)

    # CamelCase alias
    fitCentral = fit_central

    @classmethod
    def fit_raw_moments(cls, m1: float, m2: float, m3: float) -> 'APH':
        """
        Create an APH distribution from the first three raw moments.

        Args:
            m1: First raw moment E[X].
            m2: Second raw moment E[X^2].
            m3: Third raw moment E[X^3].

        Returns:
            APH distribution matching the given raw moments.

        References:
            MATLAB: matlab/src/lang/processes/APH.m (fitRawMoments)
        """
        try:
            from ..lib.thirdparty.butools.ph.canonical import APHFrom3Moments
            alpha, T = APHFrom3Moments([float(m1), float(m2), float(m3)])
            return cls(np.asarray(alpha, dtype=float).ravel(), np.asarray(T, dtype=float))
        except Exception:
            mean = float(m1)
            var = float(m2) - mean ** 2
            scv = var / (mean ** 2) if mean != 0 else 0.0
            # Standardized skewness from the third central moment.
            mu3 = float(m3) - 3.0 * mean * float(m2) + 2.0 * mean ** 3
            std = np.sqrt(var) if var > 0 else 0.0
            skew = (mu3 / std ** 3) if std > 0 else None
            return cls(mean, scv, skew)

    # CamelCase alias
    fitRawMoments = fit_raw_moments


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
        # phi[i] is the COMPLETION probability (absorbing from phase i)
        # so (1 - phi[i]) is the probability of continuing to phase i+1
        T = np.zeros((n, n))
        for i in range(n):
            T[i, i] = -rates[i]
            if i < n - 1:
                T[i, i + 1] = rates[i] * (1 - self._probs[i])

        return alpha, T

    @property
    def means(self) -> np.ndarray:
        """Get the phase means."""
        return self._means.copy()

    @property
    def probs(self) -> np.ndarray:
        """Get the continuation probabilities."""
        return self._probs.copy()

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Coxian':
        """
        Create a Coxian distribution from mean and SCV.

        Uses moment matching to construct a Coxian distribution.
        Matches MATLAB Coxian.fitMeanAndSCV.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.

        Returns:
            Coxian distribution with given mean and SCV.
        """
        tol = 1e-3  # CoarseTol
        if scv >= 1.0 - tol and scv <= 1.0 + tol:
            # Exponential
            mu = [1.0 / mean]
            phi = [1.0]
            return cls(mu, phi)
        elif scv > 0.5 + tol and scv < 1.0 - tol:
            # 2-phase Coxian with phi1=0 (all jobs pass through both phases)
            sq = np.sqrt(1.0 + 2.0 * (scv - 1.0))
            mu1 = 2.0 / mean / (1.0 + sq)
            mu2 = 2.0 / mean / (1.0 - sq)
            return cls([mu1, mu2], [0.0, 1.0])
        elif scv <= 0.5 + tol:
            # Erlang-like
            k = max(2, int(np.ceil(1.0 / scv)))
            rate = k / mean
            rates = [rate] * k
            phi = [0.0] * (k - 1) + [1.0]
            return cls(rates, phi)
        else:
            # SCV > 1+tol: HyperExp-like structure
            mu1 = 2.0 / mean
            mu2 = mu1 / (2.0 * scv)
            phi1 = 1.0 - mu2 / mu1
            return cls([mu1, mu2], [phi1, 1.0])

    @classmethod
    def fit_central(cls, mean: float, var: float, skew: float = None) -> 'Coxian':
        """
        Create a Coxian distribution from the first three central moments.

        Port of MATLAB Coxian.fitCentral and jline.lang.processes.Coxian
        .fitCentral, which agree line for line: fit the third moment exactly
        with a 2-phase Coxian, and fall back to the two-moment fit only when
        that solution misses the target SCV by more than 1%.

        Args:
            mean: Target mean.
            var: Target variance (NOT SCV -- matches the MATLAB/JAR signature).
            skew: Target skewness. When omitted, reduces to a two-moment fit.

        Returns:
            Coxian distribution with given moments.
        """
        scv = var / mean / mean
        if skew is None:
            return cls.fit_mean_and_scv(mean, scv)
        cx2 = Cox2.fit_central(mean, var, skew)
        cx = cls(1.0 / cx2.means, np.append(cx2.probs, 1.0))
        if abs(1.0 - cx.getSCV() / scv) > 0.01:
            cx = cls.fit_mean_and_scv(mean, scv)
        return cx

    # CamelCase aliases
    fitMeanAndSCV = fit_mean_and_scv
    fitMeanAndScv = fit_mean_and_scv
    fitCentral = fit_central


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

    def toTimeReversed(self) -> 'MAP':
        """Return the time-reversed MAP.

        Useful in departure-process and reversed-time arguments.

        References:
            api/mam/map_analysis.py (map_timereverse)
        """
        from ..api.mam import map_timereverse
        D0r, D1r = map_timereverse(self._D0, self._D1)
        return MAP(D0r, D1r)

    def to_time_reversed(self) -> 'MAP':
        """snake_case alias for :meth:`toTimeReversed`."""
        return self.toTimeReversed()

    def getACF(self, lags=1) -> np.ndarray:
        """Return the autocorrelation coefficients of the inter-arrival times
        at the requested lag(s).

        References:
            api/mam/map_analysis.py (map_acf)
        """
        from ..api.mam import map_acf
        return map_acf(self._D0, self._D1, lags)

    def get_acf(self, lags=1) -> np.ndarray:
        """snake_case alias for :meth:`getACF`."""
        return self.getACF(lags)

    def getIDC(self) -> float:
        """Return the asymptotic index of dispersion for counts (IDC).

        References:
            api/mam/map_analysis.py (map_idc)
        """
        from ..api.mam import map_idc
        return map_idc(self._D0, self._D1)

    def get_idc(self) -> float:
        """snake_case alias for :meth:`getIDC`."""
        return self.getIDC()

    def getMean(self) -> float:
        """Get the mean inter-arrival time."""
        # see _kb/01-model-classes.md (ME/RAP section) for rationale
        from ..api.mam import map_mean
        return float(map_mean(self._D0, self._D1))

    def getVar(self) -> float:
        """Get the variance of inter-arrival times."""
        from ..api.mam import map_var
        return float(map_var(self._D0, self._D1))

    def getSkew(self) -> float:
        """Get the skewness of the inter-arrival times."""
        from ..api.mam import map_moment
        m1 = float(map_moment(self._D0, self._D1, 1))
        m2 = float(map_moment(self._D0, self._D1, 2))
        m3 = float(map_moment(self._D0, self._D1, 3))
        var = m2 - m1 ** 2
        if var <= 0:
            return 0.0
        return (m3 - 3 * m1 * m2 + 2 * m1 ** 3) / var ** 1.5

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return self._D0.shape[0]

    def getMu(self) -> np.ndarray:
        """Get the total outgoing rate from each phase.

        mu_i = -D0(i,i), matching MATLAB Markovian.getMu. This is the rate at
        which phase i is left by ANY transition, not only by an arrival. Summing
        D1 instead counts only the arrival transitions, so the two agree only
        when D0 carries no off-diagonal mass -- which is why every renewal
        PH-like MAP in the test suite hid the difference.
        """
        return -np.diag(self._D0)

    def getPhi(self) -> np.ndarray:
        """Get the probability that a transition out of a phase is an arrival.

        phi_i = (D1*e)_i / -D0(i,i), matching MATLAB Markovian.getPhi. It equals
        1 only when every exit from phase i emits an arrival; a MAP with hidden
        phase changes has phi < 1. The D0(0,0) == 0 guard mirrors the MATLAB
        special case for an immediate process.
        """
        if self._D0[0, 0] == 0:
            return np.ones(self.getNumberOfPhases())
        return -(self._D1.sum(axis=1)) / np.diag(self._D0)

    def getInitProb(self) -> np.ndarray:
        """Get the phase distribution embedded at arrival instants.

        This is map_pie, the equilibrium of the embedded DTMC
        P = (-D0)^(-1) D1, and it is the vector that governs the inter-arrival
        time marginal. It is NOT the time-stationary phase distribution
        self._pi, which is a different vector whenever the phase process is not
        symmetric in its arrival rates. Matches MATLAB Markovian.getInitProb,
        which calls map_pie, and agrees with RAP.getInitProb on a (D0,D1) pair
        that is both a MAP and a RAP.
        """
        from ..api.mam import map_pie
        return np.asarray(map_pie(self._D0, self._D1), dtype=float).ravel()

    def getRate(self) -> float:
        """Get the arrival rate (1/mean)."""
        mean = self.getMean()
        if mean == 0 or np.isnan(mean):
            return float('inf')
        return 1.0 / mean

    def set_mean(self, target_mean: float) -> 'MAP':
        """
        Create a new MAP with the specified mean by scaling rates.

        The structure of the MAP is preserved, only the rates are scaled
        to achieve the target mean.

        Args:
            target_mean: Target mean inter-arrival/service time.

        Returns:
            New MAP with the specified mean.
        """
        current_mean = self.getMean()
        if current_mean <= 0 or np.isnan(current_mean):
            raise ValueError("Cannot scale MAP with invalid mean")

        scale = current_mean / target_mean
        new_D0 = self._D0 * scale
        new_D1 = self._D1 * scale

        return MAP(new_D0, new_D1)

    # CamelCase alias
    setMean = set_mean

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

    @classmethod
    def rand(cls, n: int = 2, seed: int = None) -> 'MAP':
        """
        Create a random MAP with n phases.

        Generates a random MAP with the specified number of phases.
        The D0 and D1 matrices are randomly generated to form a valid MAP.

        Args:
            n: Number of phases (default: 2).
            seed: Random seed for reproducibility (optional).

        Returns:
            MAP distribution with random parameters.
        """
        if seed is not None:
            np.random.seed(seed)

        # Generate random D0 (sub-generator with negative diagonal)
        D0 = np.random.rand(n, n) * 0.5
        np.fill_diagonal(D0, 0)

        # Generate random D1 (arrival transitions)
        D1 = np.random.rand(n, n) * 0.5

        # Ensure row sums are zero (D0 + D1 is a generator)
        for i in range(n):
            row_sum = D0[i, :].sum() + D1[i, :].sum()
            D0[i, i] = -row_sum

        return cls(D0, D1)


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

    Convenience class for the common 2-phase case. Mirrors
    jline.lang.processes.Cox2 and MATLAB Cox2: the phases are given by their
    RATES, and phi1 is the COMPLETION probability of phase 1 (the job absorbs
    after phase 1 with probability phi1 and continues to phase 2 with
    probability 1 - phi1), matching Coxian's phi convention.

    Args:
        mu1: Rate of phase 1.
        mu2: Rate of phase 2.
        phi1: Completion probability after phase 1.
    """

    def __init__(self, mu1: float, mu2: float, phi1: float):
        if phi1 < 0 or phi1 > 1:
            raise ValueError("Completion probability phi1 must be in [0, 1]")
        super().__init__([mu1, mu2], [phi1, 1.0])
        self._name = 'Cox2'

    @classmethod
    def fit_mean_and_scv(cls, mean: float, scv: float) -> 'Cox2':
        """
        Fit a 2-phase Coxian to a mean and SCV.

        Port of jline.lang.processes.Cox2.fitMeanAndSCV and MATLAB
        Cox2.fitMeanAndSCV, which agree branch for branch.

        This override is REQUIRED: Coxian.fit_mean_and_scv builds its result as
        cls(mu_list, phi_list), a signature Cox2 does not have, so the
        inherited classmethod raises TypeError on every call. It also may
        return an Erlang-like fit with more than two phases, which a Cox2
        cannot represent.

        A 2-phase Coxian cannot achieve SCV < 0.5; there phi1 comes out
        negative and the constructor rejects it. MATLAB and the JAR build the
        infeasible object silently instead.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation.

        Returns:
            Cox2 with the given mean and SCV.
        """
        if scv < 1.0 and scv >= 0.5:
            l0 = 2.0 / mean / (1.0 + np.sqrt(1.0 + 2.0 * (scv - 1.0)))
            l1 = 2.0 / mean / (1.0 - np.sqrt(1.0 + 2.0 * (scv - 1.0)))
            phi = 0.0
        elif scv == 1.0:
            l0 = 1.0 / mean
            l1 = 1.0 / mean
            phi = 1.0
        else:
            l0 = 2.0 / mean
            l1 = l0 / 2.0 / scv
            phi = 1.0 - l1 / l0
        return cls(l0, l1, phi)

    @classmethod
    def fit_mean(cls, mean: float) -> 'Cox2':
        """
        Fit a 2-phase Coxian to a mean alone.

        Port of jline.lang.processes.Cox2.fitMean and MATLAB Cox2.fitMean.
        """
        from ..constants import GlobalConstants  # local: avoids an import cycle
        phi = 1.0 - GlobalConstants.CoarseTol
        l0 = 1.0 / mean
        l1 = 1.0 / mean
        return cls(l0, l1, phi)

    @classmethod
    def fit_central(cls, mean: float, var: float, skew: float = None) -> 'Cox2':
        """
        Fit a 2-phase Coxian to the first three central moments.

        Port of jline.lang.processes.Cox2.fitCentral and MATLAB
        Cox2.fitCentral. This override is REQUIRED for a different reason than
        fit_mean_and_scv: it is the exact three-moment solver that
        Coxian.fit_central delegates to, and it returns a Cox2 rather than the
        general Coxian its caller then rebuilds.

        Falls back to a two-moment fit when the three-moment solution is
        infeasible, and to a mean-only fit when SCV < 0.5, as the JAR does.

        Args:
            mean: Target mean.
            var: Target variance (NOT SCV -- matches the JAR/MATLAB signature).
            skew: Target skewness. When omitted, reduces to a two-moment fit.
        """
        scv = var / mean / mean
        if skew is None:
            return cls.fit_mean_and_scv(mean, scv)

        e1 = mean
        e2 = (1.0 + scv) * e1 ** 2
        e3 = -(2.0 * e1 ** 3 - 3.0 * e1 * e2 - skew * (e2 - e1 ** 2) ** 1.5)

        # Consider the two possible solutions.
        disc = (24.0 * e1 ** 3 * e3 - 27.0 * e1 ** 2 * e2 ** 2
                - 18.0 * e1 * e2 * e3 + 18.0 * e2 ** 3 + e3 ** 2)
        den = -3.0 * e2 ** 2 + 2.0 * e1 * e3
        root = np.sqrt(disc) if disc >= 0 else np.nan

        mu11 = (2.0 * (e3 - 3.0 * e1 * e2)) / den + (3.0 * e1 * e2 - e3 + root) / den
        mu12 = (2.0 * (e3 - 3.0 * e1 * e2)) / den - (e3 - 3.0 * e1 * e2 + root) / den
        mu21 = -(3.0 * e1 * e2 - e3 + root) / den
        mu22 = (e3 - 3.0 * e1 * e2 + root) / den

        # Completion probability of phase 1. With the two rates fixed by the
        # second and third moments, the mean determines it:
        # e1 = 1/mu1 + (1-phi)/mu2  =>  phi = 1 - mu2*e1 + mu2/mu1. Verified in
        # sage/proofs/distribution_fitters.py.
        phi1 = 1.0 - mu21 * e1 + mu21 / mu11
        phi2 = 1.0 - mu22 * e1 + mu22 / mu12

        if 0 <= phi1 <= 1 and mu11 >= 0 and mu21 >= 0:
            return cls(mu11, mu21, phi1)
        elif 0 <= phi2 <= 1 and mu12 >= 0 and mu22 >= 0:
            return cls(mu12, mu22, phi2)
        else:
            line_warning('Cox2.fitCentral',
                         'Third moment could not be fitted exactly.')
            if scv >= 0.5:
                return cls.fit_mean_and_scv(mean, scv)
            return cls.fit_mean(mean)

    # CamelCase aliases. Re-bound here so they resolve to the Cox2 overrides
    # rather than to the Coxian aliases, which are bound to Coxian's functions.
    fitMeanAndSCV = fit_mean_and_scv
    fitMeanAndScv = fit_mean_and_scv
    fitMean = fit_mean
    fitCentral = fit_central


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

    def get_idc(self) -> float:
        """Index of Dispersion for Counts (snake_case alias for getIDC)."""
        return self.getIDC()

    @staticmethod
    def fitRawMomentsAndACFDecay(m1: float, m2: float, m3: float, gamma2: float) -> 'MMPP2':
        """MMPP(2) matching three raw moments and the ACF decay rate.

        Port of MATLAB MMPP2.fitRawMomentsAndACFDecay: mmpp2_fit3 on the raw
        moments, then the MMPP2 parameters read off (D1 diagonal, D0
        off-diagonal).
        """
        from line_solver.lib.kpctoolbox.mmpp import mmpp2_fit3
        if m1 <= 1e-8:
            return MMPP2(float('inf'), float('inf'), 1.0, 1.0)
        D0, D1 = mmpp2_fit3(m1, m2, m3, gamma2)
        return MMPP2(D1[0, 0], D1[1, 1], D0[0, 1], D0[1, 0])

    @staticmethod
    def fitRawMomentsAndACFLag1(m1: float, m2: float, m3: float, rho1: float) -> 'MMPP2':
        """MMPP(2) matching three raw moments and the lag-1 autocorrelation.

        rho1 = gamma2 * (1 - 1/SCV)/2 for every MMPP(2) (proved in
        sage/proofs/mmpp2_fit3.py), which is the conversion used here.
        """
        scv = (m2 - m1 * m1) / (m1 * m1)
        rho0 = (1 - 1 / scv) / 2
        return MMPP2.fitRawMomentsAndACFDecay(m1, m2, m3, rho1 / rho0)

    @staticmethod
    def fitRawMomentsAndIDC(m1: float, m2: float, m3: float, idc: float) -> 'MMPP2':
        """MMPP(2) matching three raw moments and the asymptotic index of
        dispersion, using gamma2 = (IDC - SCV)/(IDC - 1)."""
        scv = (m2 - m1 * m1) / (m1 * m1)
        gamma2 = -(scv - idc) / (-1 + idc)
        return MMPP2.fitRawMomentsAndACFDecay(m1, m2, m3, gamma2)

    @staticmethod
    def fitCentralAndACFDecay(mean: float, var: float, skew: float, gamma2: float) -> 'MMPP2':
        """MMPP(2) from central moments and the ACF decay rate (MATLAB
        MMPP2.fitCentralAndACFDecay -> mmpp2_fit2)."""
        from line_solver.lib.kpctoolbox.mmpp import mmpp2_fit2
        scv = var / (mean * mean)
        D0, D1 = mmpp2_fit2(mean, scv, skew, gamma2)
        return MMPP2(D1[0, 0], D1[1, 1], D0[0, 1], D0[1, 0])

    @staticmethod
    def fitCentralAndACFLag1(mean: float, var: float, skew: float, rho1: float) -> 'MMPP2':
        """MMPP(2) from central moments and the lag-1 autocorrelation (MATLAB
        MMPP2.fitCentralAndACFLag1 -> mmpp2_fit4)."""
        from line_solver.lib.kpctoolbox.mmpp import mmpp2_fit4
        scv = var / (mean * mean)
        D0, D1 = mmpp2_fit4(mean, scv, skew, rho1)
        return MMPP2(D1[0, 0], D1[1, 1], D0[0, 1], D0[1, 0])

    @staticmethod
    def fitCentralAndIDC(mean: float, var: float, skew: float, idc: float) -> 'MMPP2':
        """MMPP(2) from central moments and the asymptotic index of dispersion
        (MATLAB MMPP2.fitCentralAndIDC -> mmpp2_fit1)."""
        from line_solver.lib.kpctoolbox.mmpp import mmpp2_fit1
        scv = var / (mean * mean)
        D0, D1 = mmpp2_fit1(mean, scv, skew, idc)
        return MMPP2(D1[0, 0], D1[1, 1], D0[0, 1], D0[1, 0])

    @staticmethod
    def fit_mean_scv_acf(mean: float, scv: float, acf_decay: float) -> 'MMPP2':
        """
        Fit MMPP2 to match mean, SCV and the autocorrelation decay rate.

        Uses the exact closed form of mmpp2_fit3 at the minimum-skewness third
        moment, E3 = (3/2 + 1e-3) * E2^2 / E1, which is the SKEW=-1 convention
        of map_mmpp2. The fitted process matches E1, E2 and
        acf_decay = rho(k+1)/rho(k) identically; the lag-1 autocorrelation it
        realizes is acf_decay * (1 - 1/scv) / 2.

        Args:
            mean: Target mean inter-arrival time.
            scv: Target squared coefficient of variation (>= 1).
            acf_decay: ACF decay rate gamma2, in [0, 1).

        Returns:
            Fitted MMPP2 distribution.

        Raises:
            ValueError: if the request is outside the MMPP(2) feasible set.
        """
        from line_solver.lib.kpctoolbox.mmpp import mmpp2_fit3

        if mean <= 0:
            raise ValueError("MMPP2.fit_mean_scv_acf: mean must be positive")
        if scv < 1 - 1e-8:
            raise ValueError(
                "MMPP2.fit_mean_scv_acf: scv=%g is infeasible, the "
                "inter-arrival times of an MMPP(2) are over-dispersed "
                "(scv>=1)." % scv)
        if not (0 <= acf_decay < 1):
            raise ValueError(
                "MMPP2.fit_mean_scv_acf: acf_decay=%g is infeasible, the "
                "decay rate of an MMPP(2) lies in [0,1)." % acf_decay)
        if abs(scv - 1) <= 1e-8:
            # Poisson boundary: no autocorrelation is representable
            rate = 1.0 / mean
            return MMPP2(rate, rate, rate, rate)

        E1 = mean
        E2 = (1 + scv) * E1**2
        E3 = (1.5 + 1e-3) * E2**2 / E1
        D0, D1 = mmpp2_fit3(E1, E2, E3, acf_decay)
        return MMPP2(D1[0, 0], D1[1, 1], D0[0, 1], D0[1, 0])


class DMAP(ContinuousDistribution, Markovian):
    """
    Discrete-time Markovian Arrival Process.

    A DMAP models discrete-time arrival streams with correlation.
    D0 + D1 is a stochastic matrix (row sums = 1).

    Args:
        D0: Sub-stochastic transition matrix (no arrivals).
        D1: Transition matrix with arrivals.
    """

    def __init__(self, D0: Union[list, np.ndarray], D1: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'DMAP'
        self._D0 = np.atleast_2d(np.array(D0, dtype=float))
        self._D1 = np.atleast_2d(np.array(D1, dtype=float))

        n = self._D0.shape[0]
        if self._D0.shape != (n, n) or self._D1.shape != (n, n):
            raise ValueError("D0 and D1 must be square matrices of the same size")

        P = self._D0 + self._D1
        if not np.allclose(P.sum(axis=1), 1.0):
            raise ValueError("D0 + D1 must have row sums equal to 1 (stochastic matrix)")

        self._pi = self._compute_stationary()

    def _compute_stationary(self) -> np.ndarray:
        """Compute stationary distribution of the embedded DTMC."""
        n = self._D0.shape[0]
        I = np.eye(n)
        D0inv = linalg.inv(I - self._D0)
        P = D0inv @ self._D1

        A = np.vstack([(P.T - np.eye(n)), np.ones(n)])
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
        return self._D0.copy()

    @property
    def D1(self) -> np.ndarray:
        return self._D1.copy()

    @property
    def pi(self) -> np.ndarray:
        return self._pi.copy()

    def getD0(self) -> np.ndarray:
        return self._D0.copy()

    def getD1(self) -> np.ndarray:
        return self._D1.copy()

    def getMean(self) -> float:
        """Mean inter-arrival time: pi * (I-D0)^{-1} * e"""
        n = self._D0.shape[0]
        D0inv = linalg.inv(np.eye(n) - self._D0)
        e = np.ones(n)
        return float(self._pi @ D0inv @ e)

    def getVar(self) -> float:
        """Variance of inter-arrival times."""
        n = self._D0.shape[0]
        D0inv = linalg.inv(np.eye(n) - self._D0)
        e = np.ones(n)
        mean = float(self._pi @ D0inv @ e)
        second_moment = float(2 * self._pi @ D0inv @ D0inv @ e) - mean
        return second_moment - mean ** 2

    def getNumberOfPhases(self) -> int:
        return self._D0.shape[0]

    def getMu(self) -> np.ndarray:
        return (np.eye(self._D0.shape[0]) - self._D0).sum(axis=1)

    def getPhi(self) -> np.ndarray:
        return np.ones(self.getNumberOfPhases())

    def getInitProb(self) -> np.ndarray:
        return self._pi.copy()

    def getRate(self) -> float:
        mean = self.getMean()
        if mean == 0 or np.isnan(mean):
            return float('inf')
        return 1.0 / mean

    def set_mean(self, target_mean: float) -> 'DMAP':
        """Scale DMAP to match target mean."""
        n = self._D0.shape[0]
        I = np.eye(n)
        current_mean = self.getMean()
        if current_mean <= 0 or np.isnan(current_mean):
            raise ValueError("Cannot scale DMAP with invalid mean")
        scale = current_mean / target_mean
        new_ImD0 = (I - self._D0) * scale
        new_D0 = I - new_ImD0
        new_D1 = self._D1 * scale
        P = new_D0 + new_D1
        rs = P.sum(axis=1, keepdims=True)
        new_D0 = new_D0 / rs
        new_D1 = new_D1 / rs
        return DMAP(new_D0, new_D1)

    setMean = set_mean

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate random inter-arrival times (integer-valued)."""
        if rng is None:
            rng = np.random.default_rng()

        num_phases = len(self._pi)
        samples = np.zeros(n)
        phase = rng.choice(num_phases, p=self._pi)

        for i in range(n):
            t = 0
            while True:
                t += 1
                probs = np.concatenate([self._D0[phase, :], self._D1[phase, :]])
                probs = np.maximum(probs, 0)
                total = probs.sum()
                if total <= 0:
                    break
                probs /= total
                next_state = rng.choice(2 * num_phases, p=probs)
                if next_state >= num_phases:
                    samples[i] = t
                    phase = next_state - num_phases
                    break
                else:
                    phase = next_state

        return samples

    @classmethod
    def rand(cls, n: int = 2, seed: int = None) -> 'DMAP':
        """Create a random DMAP with n phases."""
        rng = np.random.default_rng(seed)
        D0 = rng.random((n, n))
        for i in range(n):
            D0[i, :] /= (D0[i, :].sum() + 0.5 + rng.random())
        remaining = 1.0 - D0.sum(axis=1)
        D1 = rng.random((n, n))
        for i in range(n):
            D1[i, :] = D1[i, :] / D1[i, :].sum() * remaining[i]
        return cls(D0, D1)


def _cme_table_entry(order: int) -> dict:
    """Select the CME table entry realizing the given number of phases.

    The table is keyed by the number of harmonic terms n, so an order of 2*n+1
    phases maps to the entries with that n. Several entries can share an n (the
    'full' and 'approx' optimizations), and the most concentrated one is taken,
    matching the selection rule of the CME inverse Laplace transform.
    """
    from ..lib.thirdparty.iltcme import cme_parameters

    order = int(order)
    if order < 3 or order % 2 == 0:
        raise ValueError(
            "CME order must be an odd integer of the form 2*n+1 with n >= 1, "
            "got %d" % order)
    n = (order - 1) // 2
    candidates = [entry for entry in cme_parameters() if int(entry['n']) == n]
    if not candidates:
        available = sorted({2 * int(e['n']) + 1 for e in cme_parameters()})
        nearest = min(available, key=lambda o: abs(o - order))
        raise ValueError(
            "No tabulated CME of order %d; the nearest available order is %d. "
            "Use CME.getSupportedOrders() for the full list." % (order, nearest))
    return min(candidates, key=lambda entry: float(entry['cv2']))


def _cme_representation(order: int) -> Tuple[np.ndarray, np.ndarray, float]:
    """Build the unit-mean (alpha, A) matrix-exponential form of a CME.

    The tabulated density is
    f(x) = mu1*exp(-mu1*x)*(c + sum_k [a_k*cos(k*w*mu1*x) + b_k*sin(k*w*mu1*x)]),
    and A = blkdiag(-mu1, mu1*[[-1,-k*w],[k*w,-1]] for k=1..n) reproduces the
    exponential envelope in its first phase and the k-th harmonic in its k-th
    2x2 rotation block, since expm(mu1*x*[[-1,-kw],[kw,-1]]) is exp(-mu1*x)
    times the rotation by k*w*mu1*x. The entries of alpha follow by matching
    -alpha*expm(A*x)*A*e term by term: with wk = k*w and d = 2*(1+wk^2),

        alpha_0     = c
        alpha_{2k-1} = ((1+wk)*a_k - (1-wk)*b_k)/d
        alpha_{2k}   = ((1-wk)*a_k + (1+wk)*b_k)/d

    The result has unit mean and sums to one, as any (alpha, A) whose density
    integrates to one must.

    Returns:
        Tuple of (alpha, A, scv), with scv the tabulated cv2 of the entry.
    """
    entry = _cme_table_entry(order)
    a = np.asarray(entry['a'], dtype=float)
    b = np.asarray(entry['b'], dtype=float)
    c = float(entry['c'])
    mu1 = float(entry['mu1'])
    w = float(entry['omega'])
    n = int(entry['n'])

    size = 2 * n + 1
    A = np.zeros((size, size))
    alpha = np.zeros(size)
    A[0, 0] = -mu1
    alpha[0] = c
    for k in range(1, n + 1):
        i = 2 * k - 1
        wk = k * w
        A[i, i] = -mu1
        A[i, i + 1] = -wk * mu1
        A[i + 1, i] = wk * mu1
        A[i + 1, i + 1] = -mu1
        d = 2.0 * (1.0 + wk * wk)
        alpha[i] = ((1.0 + wk) * a[k - 1] - (1.0 - wk) * b[k - 1]) / d
        alpha[i + 1] = ((1.0 - wk) * a[k - 1] + (1.0 + wk) * b[k - 1]) / d

    # see _kb/01-model-classes.md (Concentrated matrix exponential (CME)) for rationale
    alpha /= alpha.sum()

    return alpha, A, float(entry['cv2'])


class ME(ContinuousDistribution, Markovian):
    """
    Matrix Exponential distribution.

    ME distributions generalize Phase-Type distributions by allowing
    the initial vector alpha to have entries outside [0,1] and the
    matrix A to have arbitrary structure (not necessarily a sub-generator).

    Args:
        alpha: Initial vector (may have negative entries).
        A: Matrix parameter (eigenvalues must have negative real parts).
        checkDensity: Scan the density for a negative value (default True). It
            is set to False only by subclasses whose representation is a
            density by construction, such as CME, where the scan would cost
            O(1e5) propagations of a large matrix and can never fire.
    """

    def __init__(self, alpha: Union[list, np.ndarray], A: Union[list, np.ndarray],
                 checkDensity: bool = True):
        super().__init__()
        self._name = 'ME'
        self._alpha = np.atleast_1d(np.array(alpha, dtype=float)).flatten()
        self._A = np.atleast_2d(np.array(A, dtype=float))

        n = len(self._alpha)
        if self._A.shape != (n, n):
            raise ValueError(f"A must be {n}x{n} to match alpha of length {n}")

        # Strict validation via BuTools, matching MATLAB ME.m and
        # jline.lang.processes.ME. Checking only that every eigenvalue has a
        # negative real part is weaker than CheckMERepresentation, which also
        # requires the sum of alpha to lie in [0,1] and the dominant eigenvalue
        # of A to be real; the weaker test accepted representations that MATLAB
        # and the JAR reject.
        alpha_row = self._alpha.reshape(1, n)
        if not _check_me_representation(alpha_row, self._A):
            raise ValueError(
                "Invalid ME representation: Check that A is square, alpha and A "
                "have compatible dimensions, all eigenvalues of A have negative "
                "real parts, and the dominant eigenvalue is real.")

        # see _kb/01-model-classes.md (ME/RAP section) for rationale
        if checkDensity and _butools_settings().checkInput:
            is_neg, f_min, t_min = _scan_negative_density(alpha_row, self._A)
            if is_neg:
                line_warning('ME',
                             'The ME representation has a negative density: '
                             'f(t) = -alpha*expm(A*t)*A*e reaches %g at '
                             't = %g. Moments and transforms remain well '
                             'defined, but evalPDF returns negative values and '
                             'sample() will not reproduce a proper '
                             'distribution.' % (f_min, t_min))

        # Inverse-CDF table for sample(), built lazily and cached so that
        # repeated single-sample draws do not rebuild it.
        self._cdf_x = None
        self._cdf_F = None
        self._cdf_tail_surv = None
        self._spectral = None

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

    def getSCV(self) -> float:
        """Get the squared coefficient of variation, var / mean^2."""
        mean = self.getMean()
        if mean == 0 or np.isnan(mean):
            return float('nan')
        return self.getVar() / (mean ** 2)

    def evalLST(self, s: float) -> float:
        """Evaluate the Laplace-Stieltjes transform at s.

        LST(s) = alpha * (s*I - A)^(-1) * (-A) * e, mirroring MATLAB ME.evalLST.
        """
        n = len(self._alpha)
        e = np.ones(n)
        rhs = -(self._A @ e)
        return float(self._alpha @ linalg.solve(s * np.eye(n) - self._A, rhs))

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix of the process representation (equals A)."""
        return self._A.copy()

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix of the process representation, -A*e*alpha."""
        e = np.ones(len(self._alpha))
        return np.outer(-(self._A @ e), self._alpha)

    def getProcess(self) -> list:
        """Get the process representation as ``[D0, D1] = [A, -A*e*alpha]``.

        Same convention as MATLAB ME.getProcess (a {D0,D1} cell) and
        jline.lang.processes.ME.getProcess (a two-entry MatrixCell).
        """
        return [self.getD0(), self.getD1()]

    def get_process(self) -> list:
        """snake_case alias for :meth:`getProcess`."""
        return self.getProcess()

    def getMu(self) -> np.ndarray:
        """Get the rate out of each phase."""
        return -np.diag(self._A)

    def getPhi(self) -> np.ndarray:
        """Get the completion probability out of each phase."""
        mu = self.getMu()
        t = -self._A.sum(axis=1)
        return np.where(mu > 0, t / mu, 0.0)

    def getInitProb(self) -> np.ndarray:
        """Get the initial vector alpha (entries may be negative for an ME)."""
        return self._alpha.copy()

    def _build_cdf_table(self) -> None:
        """Build and cache the monotone inverse-CDF table used by sample().

        The table holds F(x) = 1 - alpha*exp(A*x)*e on a uniform grid that is
        extended until the survival function drops below ``_ME_TAIL_TOL``, so the
        tail is covered rather than truncated. Successive rows are obtained by
        propagating the row vector with a single ``expm(A*dt)``, which keeps the
        cost linear in the number of grid points.
        """
        if self._cdf_x is not None:
            return

        n = len(self._alpha)
        e = np.ones(n)
        mean = self.getMean()
        var = self.getVar()
        sigma = np.sqrt(var) if np.isfinite(var) and var > 0 else 0.0
        span = mean + 10.0 * sigma
        if not np.isfinite(span) or span <= 0:
            # Degenerate moment estimates: fall back to the slowest decay mode.
            span = -10.0 / _dominant_eigenvalue(self._A)
        dt = span / _ME_GRID_STEPS

        E = linalg.expm(self._A * dt)
        xs = [0.0]
        surv = [float(self._alpha @ e)]
        v = self._alpha.copy()
        k = 0
        while surv[-1] > _ME_TAIL_TOL and k < _ME_MAX_STEPS:
            v = v @ E
            k += 1
            xs.append(k * dt)
            surv.append(float(v @ e))

        x_arr = np.array(xs)
        F = 1.0 - np.array(surv)
        # An ME density may dip negative, which makes F non-monotone; force a
        # monotone table so the inversion below is well posed.
        F = np.maximum.accumulate(F)
        self._cdf_x = x_arr
        self._cdf_F = F
        self._cdf_tail_surv = max(1.0 - F[-1], 0.0)
        self._spectral = _spectral_expm_data(self._A)

    def _survival_and_density(self, x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Return (alpha*exp(A*x)*e, -alpha*exp(A*x)*A*e) evaluated at each x."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        n = len(self._alpha)
        e = np.ones(n)
        if self._spectral is not None:
            lam, V, Vinv = self._spectral
            c = (self._alpha @ V) * (Vinv @ e)
            E = np.exp(np.outer(x, lam))
            surv = np.real(E @ c)
            dens = np.real(E @ (-lam * c))
            return surv, dens
        surv = np.empty(x.shape)
        dens = np.empty(x.shape)
        t = -(self._A @ e)
        for i, xi in enumerate(x):
            expm = linalg.expm(self._A * xi)
            av = self._alpha @ expm
            surv[i] = float(av @ e)
            dens[i] = float(av @ t)
        return surv, dens

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate n samples by inverting the CDF.

        A PH-style walk over the phase process is not applicable: alpha may have
        negative entries and A need not be a sub-generator, so there is no
        Markov chain to walk. Instead the CDF F(x) = 1 - alpha*exp(A*x)*e is
        tabulated on a grid covering the tail, inverted by binary search on the
        table and refined by safeguarded Newton steps against the exact density
        -alpha*exp(A*x)*A*e. Uniforms beyond the last tabulated point are mapped
        through the exponential tail of the dominant eigenvalue of A rather than
        clamped to the grid endpoint, which would truncate the tail and bias the
        mean low.
        """
        if rng is None:
            rng = np.random.default_rng()
        self._build_cdf_table()

        u = np.asarray(rng.random(n), dtype=float)
        x_grid = self._cdf_x
        F = self._cdf_F
        out = np.zeros(n)

        idx = np.searchsorted(F, u, side='left')
        below = idx <= 0
        beyond = idx >= len(F)
        interior = ~(below | beyond)

        out[below] = x_grid[0]

        if np.any(beyond):
            # Exponential tail: S(x) ~ S(x_last) * exp(eta * (x - x_last)) with
            # eta the dominant (real, negative) eigenvalue of A.
            eta = _dominant_eigenvalue(self._A)
            s_last = max(self._cdf_tail_surv, np.finfo(float).tiny)
            target = np.maximum(1.0 - u[beyond], np.finfo(float).tiny)
            out[beyond] = x_grid[-1] + np.log(target / s_last) / eta

        if np.any(interior):
            hi_i = idx[interior]
            lo_i = hi_i - 1
            lo = x_grid[lo_i]
            hi = x_grid[hi_i]
            Flo = F[lo_i]
            Fhi = F[hi_i]
            span = Fhi - Flo
            frac = np.where(span > 0, (u[interior] - Flo) / np.where(span > 0, span, 1.0), 0.5)
            x = lo + frac * (hi - lo)
            ui = u[interior]

            for _ in range(_ME_NEWTON_ITERS):
                surv, dens = self._survival_and_density(x)
                resid = (1.0 - surv) - ui
                if np.max(np.abs(resid)) < _ME_NEWTON_TOL:
                    break
                # Tighten the bracket with the freshly evaluated residual.
                lo = np.where(resid < 0, x, lo)
                hi = np.where(resid > 0, x, hi)
                step = np.where(dens > 0, resid / np.where(dens > 0, dens, 1.0), 0.0)
                x_new = x - step
                # Fall back to bisection whenever Newton leaves the bracket.
                outside = (x_new <= lo) | (x_new >= hi) | ~np.isfinite(x_new)
                x = np.where(outside, 0.5 * (lo + hi), x_new)

            out[interior] = x

        return out

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

    def get_scv(self) -> float:
        """Get the squared coefficient of variation."""
        return self.getSCV()

    def eval_lst(self, s: float) -> float:
        """Evaluate the Laplace-Stieltjes transform at s."""
        return self.evalLST(s)



def fit_me_mean_scv(mean: float, scv: float, max_phases: int = 0) -> 'ME':
    """Fit a matrix exponential to a given mean and squared coefficient of variation.

    For ``scv < 1`` the fit is the convolution ``X = c*Y + Z`` of a scaled
    concentrated matrix exponential ``Y`` (unit mean, minimal SCV ``sY`` for its
    order) with an independent exponential ``Z``. Writing ``c + d = mean`` and
    ``c^2*sY + d^2 = scv*mean^2`` gives

        c = mean*(1 - sqrt(1 - (1+sY)*(1-scv)))/(1 + sY),   d = mean - c,

    so every target in ``[sY/(1+sY), 1]`` is matched EXACTLY in ``2n+2`` phases.
    The exponential tail is what makes the convolution reach up to SCV 1; the
    concentrated part is what makes it reach far below the Erlang bound
    ``1/order`` at the same order.

    The order is the smallest tabulated one that reaches the target, capped by
    ``max_phases`` when given: with a phase budget an Erlang can only reach
    ``1/max_phases``, while this construction reaches ``O(1/max_phases^2)``, and
    the residual SCV is then the closest achievable from below.

    ``scv >= 1`` is outside the range of a concentrated ME (its SCV never
    exceeds 0.34), and the caller keeps its own hyperexponential fit there.

    Args:
        mean: Target mean, positive.
        scv: Target squared coefficient of variation, in (0, 1).
        max_phases: Optional cap on the number of phases, 0 for no cap.

    Returns:
        An :class:`ME` with the requested mean and, budget permitting, SCV.

    Raises:
        ValueError: if mean is not positive, or scv is not in (0, 1).
    """
    mean = float(mean)
    scv = float(scv)
    if not np.isfinite(mean) or mean <= 0:
        raise ValueError("fit_me_mean_scv mean must be a positive finite number")
    if not np.isfinite(scv) or scv <= 0 or scv >= 1:
        raise ValueError("fit_me_mean_scv requires 0 < scv < 1; use a "
                         "hyperexponential for scv >= 1 and a CME for scv = 0")

    # Smallest tabulated order whose convolution range covers the target, subject
    # to the phase budget. reach(order) = sY/(1+sY) is the minimum SCV of the
    # convolution, which sits just below the sY of the CME alone.
    best_order = None
    for order in CME.getSupportedOrders():
        if max_phases and order + 1 > max_phases:
            continue
        s_y = CME.getMinSCV(order)
        if s_y / (1.0 + s_y) <= scv:
            best_order = order
            break
        best_order = order  # budget-limited: keep the most concentrated one that fits
    if best_order is None:
        raise ValueError("no CME order fits a budget of %d phases; the smallest "
                         "is 3 phases plus one exponential" % max_phases)

    alpha_y, A_y, s_y = _cme_representation(best_order)
    reach = s_y / (1.0 + s_y)
    if scv < reach:
        # Budget-limited: the target is below what this order can reach, so the
        # most concentrated member of the family is returned and the caller gets
        # the closest achievable SCV rather than a silent Erlang truncation.
        c = mean / (1.0 + s_y)
    else:
        c = mean * (1.0 - math.sqrt(1.0 - (1.0 + s_y) * (1.0 - scv))) / (1.0 + s_y)
    d = mean - c

    n = len(alpha_y)
    if d <= mean * 1e-12:
        return CME(mean, best_order)
    if c <= mean * 1e-12:
        return ME(np.array([1.0]), np.array([[-1.0 / mean]]))

    # Convolution of two matrix exponentials: the exit flow of the first block
    # feeds the entry of the second, exactly as for a phase-type.
    size = n + 1
    alpha = np.zeros(size)
    alpha[:n] = alpha_y
    A = np.zeros((size, size))
    A[:n, :n] = A_y / c
    A[:n, n] = -(A_y / c) @ np.ones(n)
    A[n, n] = -1.0 / d
    return ME(alpha, A, checkDensity=False)


class CME(ME):
    """
    Concentrated Matrix Exponential distribution.

    A CME is the matrix-exponential distribution of odd order ``2*n+1`` whose
    squared coefficient of variation is (numerically) minimal for that order,
    from the tables of Horvath, Horvath and Telek. Its SCV decays as O(1/n^2),
    so it goes far below the Erlang bound ``1/order`` reachable by a phase-type
    distribution of the same order: order 101 gives SCV 3.9e-4, where Erlang-101
    gives 9.9e-3.

    The density of the unit-mean CME with ``n`` harmonic terms is

        f(x) = mu1 * exp(-mu1*x) * (c + sum_k [a_k*cos(k*w*mu1*x) + b_k*sin(k*w*mu1*x)])

    with ``w = omega``, which is exactly ``alpha*expm(A*x)*(-A*e)`` for the
    block-diagonal ``A = blkdiag(-mu1, mu1*[[-1,-k*w],[k*w,-1]] for k=1..n)``.
    The parameters ``a``, ``b``, ``c``, ``omega`` and ``mu1`` are read from the
    same ``iltcme.json`` table used by the CME inverse Laplace transform.

    Args:
        mean: Mean of the distribution (the tabulated CME has unit mean and is
            rescaled by ``A/mean``).
        order: Number of phases, an odd integer ``2*n+1`` with ``n`` present in
            the table.
    """

    def __init__(self, mean: float, order: int):
        alpha, A, scv = _cme_representation(order)
        mean = float(mean)
        if not np.isfinite(mean) or mean <= 0:
            raise ValueError("CME mean must be a positive finite number")
        # see _kb/01-model-classes.md (Concentrated matrix exponential (CME)) for rationale
        super().__init__(alpha, A / mean, checkDensity=False)
        # see _kb/01-model-classes.md (Concentrated matrix exponential (CME)) for rationale
        self._name = 'ME'
        self._cme_order = int(order)
        self._cme_scv = scv

    def getOrder(self) -> int:
        """Get the CME order, i.e. the number of phases."""
        return self._cme_order

    def get_order(self) -> int:
        """snake_case alias for :meth:`getOrder`."""
        return self.getOrder()

    @staticmethod
    def getSupportedOrders() -> list:
        """Get the sorted list of CME orders (phase counts) in the table."""
        from ..lib.thirdparty.iltcme import cme_parameters
        return sorted({2 * int(entry['n']) + 1 for entry in cme_parameters()})

    @staticmethod
    def getMinSCV(order: int) -> float:
        """Get the tabulated minimal SCV attained by a CME of the given order."""
        return _cme_representation(order)[2]

    @staticmethod
    def fitMeanAndSCV(mean: float, scv: float) -> 'CME':
        """Create the lowest-order CME with the given mean and SCV at most scv.

        Args:
            mean: Target mean.
            scv: Target squared coefficient of variation, an upper bound. The
                lowest tabulated order whose minimal SCV does not exceed it is
                selected, so the returned distribution is at least as
                concentrated as requested.

        Returns:
            The CME of that order, rescaled to the requested mean.

        Raises:
            ValueError: if no tabulated order reaches the requested SCV.
        """
        from ..lib.thirdparty.iltcme import cme_parameters
        best_order = None
        best_scv = None
        for entry in cme_parameters():
            if float(entry['cv2']) <= scv:
                order = 2 * int(entry['n']) + 1
                if best_order is None or order < best_order:
                    best_order = order
                    best_scv = float(entry['cv2'])
        if best_order is None:
            raise ValueError(
                "No tabulated CME reaches SCV %g; the most concentrated entry "
                "has SCV %g at order %d" % (
                    scv,
                    min(float(e['cv2']) for e in cme_parameters()),
                    2 * max(int(e['n']) for e in cme_parameters()) + 1))
        _ = best_scv
        return CME(mean, best_order)

    @staticmethod
    def fit_mean_and_scv(mean: float, scv: float) -> 'CME':
        """snake_case alias for :meth:`fitMeanAndSCV`."""
        return CME.fitMeanAndSCV(mean, scv)


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

    def getNumberOfTypes(self) -> int:
        """Get the number of arrival types (MATLAB/JAR API parity)."""
        return self.num_types

    def to_m3a(self) -> list:
        """
        Return the process in the canonical M3A layout
        [D0, D1_agg, D11, ..., D1K] used by the mam API (mmap_* functions)
        and by the NetworkStruct proc entries, where D1_agg = sum of the
        per-type matrices. The internal _process stays in the constructor
        layout [D0, D11, ..., D1K].
        """
        D1_agg = sum(self._process[1:])
        return [self._process[0].copy(), D1_agg] + [Dk.copy() for Dk in self._process[1:]]

    def to_map(self) -> 'MAP':
        """Aggregate MAP over all marks: MAP(D0, sum(D1k)). Mirrors MATLAB toMAP."""
        return MAP(self._process[0].copy(), sum(self._process[1:]))

    def to_marginal_map(self, k: int) -> 'MAP':
        """
        Marginal MAP of mark k (1-based): arrivals fire only on D1k, with the
        other marks' transitions folded into the hidden part, i.e.
        MAP(D0 + D1_agg - D1k, D1k). Mirrors MATLAB MarkedMAP.toMAPs for a
        single type.
        """
        if k < 1 or k > self.num_types:
            raise IndexError(f"Mark index {k} out of range")
        D1k = self._process[k]
        Df = self._process[0] + sum(self._process[1:]) - D1k
        return MAP(Df.copy(), D1k.copy())

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

        # Strict validation via BuTools, matching MATLAB RAP.m and
        # jline.lang.processes.RAP. Nothing used to be checked here beyond the
        # shapes, so representations whose H0+H1 has nonzero row sums, or whose
        # dominant eigenvalue of H0 is complex, were accepted in Python only.
        if not _check_rap_representation(self._H0, self._H1):
            raise ValueError(
                "Invalid RAP representation: Check that H0 and H1 are square "
                "matrices of the same size, H0 + H1 forms a valid infinitesimal "
                "generator (row sums = 0), all eigenvalues of H0 have negative "
                "real parts, and the dominant eigenvalue of H0 is real.")

        # Compute stationary distribution
        self._pi = self._compute_stationary()
        self._spectral = None

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

    def getPie(self) -> np.ndarray:
        """Get the phase vector embedded at arrival epochs.

        This is the equilibrium distribution of the embedded DTMC
        P = (-H0)^(-1) H1, i.e. pie = pi*H1 / (pi*H1*e). It, and not the
        time-stationary pi, is the vector that governs the inter-arrival time
        marginal, matching MATLAB map_pie.
        """
        e = np.ones(self._H0.shape[0])
        num = self._pi @ self._H1
        denom = float(num @ e)
        if denom == 0:
            raise ValueError("RAP has zero arrival rate: pi*H1*e = 0")
        return num / denom

    def get_pie(self) -> np.ndarray:
        """snake_case alias for :meth:`getPie`."""
        return self.getPie()

    def evalCDF(self, x: float) -> float:
        """Evaluate the CDF of the inter-arrival time at point x."""
        if x <= 0:
            return 0.0
        n = self._H0.shape[0]
        e = np.ones(n)
        expm = linalg.expm(self._H0 * x)
        return 1.0 - float(self.getPie() @ expm @ e)

    def evalPDF(self, x: float) -> float:
        """Evaluate the PDF of the inter-arrival time at point x."""
        if x < 0:
            return 0.0
        expm = linalg.expm(self._H0 * x)
        h = -self._H0.sum(axis=1)
        return float(self.getPie() @ expm @ h)

    def evalLST(self, s: float) -> float:
        """Evaluate the Laplace-Stieltjes transform of the inter-arrival time.

        LST(s) = pie * (s*I - H0)^(-1) * (-H0) * e, mirroring MATLAB RAP.evalLST.
        """
        n = self._H0.shape[0]
        e = np.ones(n)
        rhs = -(self._H0 @ e)
        return float(self.getPie() @ linalg.solve(s * np.eye(n) - self._H0, rhs))

    def getSCV(self) -> float:
        """Get the squared coefficient of variation of the inter-arrival time."""
        mean = self.getMean()
        if mean == 0 or np.isnan(mean):
            return float('nan')
        return self.getVar() / (mean ** 2)

    def getD0(self) -> np.ndarray:
        """Get the D0 matrix of the process representation (equals H0)."""
        return self._H0.copy()

    def getD1(self) -> np.ndarray:
        """Get the D1 matrix of the process representation (equals H1)."""
        return self._H1.copy()

    def getProcess(self) -> list:
        """Get the process representation ``[H0, H1]``, as MATLAB RAP.getProcess."""
        return [self._H0.copy(), self._H1.copy()]

    def get_process(self) -> list:
        """snake_case alias for :meth:`getProcess`."""
        return self.getProcess()

    def getACF(self, lags=1) -> np.ndarray:
        """Get the autocorrelation of the inter-arrival times at the given lags."""
        from ..api.mam import map_acf
        return map_acf(self._H0, self._H1, lags)

    def get_acf(self, lags=1) -> np.ndarray:
        """snake_case alias for :meth:`getACF`."""
        return self.getACF(lags)

    def getIDC(self) -> float:
        """Get the asymptotic index of dispersion for counts."""
        from ..api.mam import map_idc
        return float(map_idc(self._H0, self._H1))

    def get_idc(self) -> float:
        """snake_case alias for :meth:`getIDC`."""
        return self.getIDC()

    def getInitProb(self) -> np.ndarray:
        """Get the phase vector embedded at arrival epochs."""
        return self.getPie()

    def getMu(self) -> np.ndarray:
        """Get the total outgoing rate from each phase.

        mu_i = -H0(i,i), the same MATLAB Markovian.getMu formula MAP uses. A RAP
        whose matrices happen to be nonnegative IS a MAP, so the two classes have
        to answer identically on such a pair; summing H1 instead counts only the
        arrival transitions and contradicted MAP on exactly that degenerate case.
        """
        return -np.diag(self._H0)

    def getPhi(self) -> np.ndarray:
        """Get the probability that a transition out of a phase is an arrival.

        phi_i = (H1*e)_i / -H0(i,i), matching MATLAB Markovian.getPhi and MAP.
        The H0(0,0) == 0 guard mirrors the MATLAB special case for an immediate
        process.
        """
        if self._H0[0, 0] == 0:
            return np.ones(self._H0.shape[0])
        return -(self._H1.sum(axis=1)) / np.diag(self._H0)

    def sample(self, n: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Generate n successive inter-arrival times along one sample path.

        A RAP has no underlying Markov chain over phases to walk, so sampling
        carries the conditional phase vector v (normalized so that v*e = 1)
        instead. Starting from the arrival-embedded vector pie, each draw
        inverts the conditional survival function v*exp(H0*x)*e at a uniform
        variate and then updates
        ``v <- v*exp(H0*x)*H1 / (v*exp(H0*x)*H1*e)``.
        This reproduces the autocorrelation of the process; drawing from the ME
        marginal alone would produce independent inter-arrival times.
        """
        if rng is None:
            rng = np.random.default_rng()

        if self._spectral is None:
            self._spectral = _spectral_expm_data(self._H0)
        spectral = self._spectral

        e = np.ones(self._H0.shape[0])
        v = self.getPie()
        mean = self.getMean()
        x_guess = mean if np.isfinite(mean) and mean > 0 else 1.0

        out = np.zeros(n)
        for i in range(n):
            target = 1.0 - float(rng.random())
            x = _rap_invert_survival(v, target, self._H0, spectral, e, x_guess)
            out[i] = x
            w = _rap_row_expm(v, x, self._H0, spectral)
            vv = w @ self._H1
            s = float(vv @ e)
            if not np.isfinite(s) or s <= 0:
                raise ValueError(
                    "RAP sampling reached a phase vector with nonpositive "
                    "arrival mass (v*exp(H0*x)*H1*e = %g); the (H0,H1) pair "
                    "does not define a proper arrival process" % s)
            v = vv / s
        return out

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

    MMDP models a process where a deterministic RATE is modulated by an
    underlying Markov chain: in state i arrivals occur deterministically at
    rate r_i, i.e. one every 1/r_i time units.

    Parameterized by rates, mirroring MATLAB MMDP.m and the JAR MMDP.java
    (MMDP(Q, R), R = diag(r)). This class previously took inter-arrival TIMES
    d, which made it a different process sharing a name: it reported
    getMean() = pi.d (the state-averaged time), whereas averaging rates gives
    1/(pi.r). The two differ whenever the rates are not all equal, so the same
    model meant different things per codebase. MATLAB is the reference.

    Args:
        Q: Generator matrix of the modulating Markov chain (row sums 0).
        R: n x n diagonal matrix of deterministic rates, or an n-vector of them.
    """

    def __init__(self, Q: Union[list, np.ndarray], R: Union[list, np.ndarray]):
        super().__init__()
        self._name = 'MMDP'
        self._Q = np.atleast_2d(np.array(Q, dtype=float))
        # Accept the diagonal matrix or the bare rate vector, as MATLAB does.
        R_arr = np.array(R, dtype=float)
        if R_arr.ndim == 2 and R_arr.shape[0] == R_arr.shape[1] and R_arr.shape[0] > 1:
            self._r = np.diag(R_arr).astype(float).copy()
        else:
            self._r = np.atleast_1d(R_arr).ravel().astype(float)

        n = self._Q.shape[0]
        if self._Q.shape != (n, n):
            raise ValueError("Q must be a square matrix")
        if len(self._r) != n:
            raise ValueError(f"R must be an {n}x{n} diagonal matrix or a length-{n} vector")

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

    def r(self) -> np.ndarray:
        """Vector of deterministic rates, one per state. Mirrors MATLAB MMDP.r."""
        return self._r.copy()

    def R(self) -> np.ndarray:
        """Diagonal rate matrix. Mirrors MATLAB MMDP.R."""
        return np.diag(self._r)

    def getMeanRate(self) -> float:
        """Stationary mean rate pi.r. Mirrors MATLAB MMDP.getMeanRate."""
        return float(np.dot(self._pi, self._r))

    def getRate(self) -> float:
        """Alias of getMeanRate, as in MATLAB MMDP.getRate."""
        return self.getMeanRate()

    def getMean(self) -> float:
        """Mean inter-arrival time, the inverse of the mean rate.

        Mirrors MATLAB MMDP.getMean. This is 1/(pi.r), not the state-averaged
        inter-arrival time pi.(1/r): the rate is what the chain modulates.
        """
        rate = self.getMeanRate()
        return 1.0 / rate if rate > 0 else float('inf')

    def getSCV(self) -> float:
        """SCV of the modulated rate. Mirrors MATLAB MMDP.getSCV."""
        mean_rate = self.getMeanRate()
        var_rate = float(np.dot(self._pi, self._r ** 2)) - mean_rate ** 2
        if mean_rate > 0:
            return var_rate / (mean_rate ** 2)
        return float('nan')

    def getVar(self) -> float:
        """Variance, as MATLAB's Distribution base derives it: SCV * mean^2.

        MMDP.m does not override getVar, so this reproduces the inherited
        relation rather than re-deriving a variance from the rates.
        """
        return self.getSCV() * self.getMean() ** 2

    def getNumberOfPhases(self) -> int:
        """Get the number of phases."""
        return len(self._r)



class MMDP2(MMDP):
    """
    Markov-Modulated Deterministic Process with 2 states.

    Convenience class for the common 2-state case.

    Parameterized by deterministic RATES, mirroring MATLAB MMDP2.m and the JAR
    MMDP2.java (which build R = diag([r0, r1])). This class previously took
    deterministic TIMES d0/d1, which made it a different stochastic object
    sharing a name: averaging times as pi.d yields the arithmetic mean of the
    per-state times, whereas the rate parameterization averages rates and
    reports 1/(pi.r). The two disagree whenever r0 ~= r1, so a model written by
    one codebase and read by another silently changed meaning. MATLAB is the
    reference, so the rate form is canonical here too; the JSON wire keys are
    r0/r1 accordingly.

    Args:
        r0: Deterministic rate in state 0.
        r1: Deterministic rate in state 1.
        sigma0: Transition rate from state 0 to state 1.
        sigma1: Transition rate from state 1 to state 0.
    """

    def __init__(self, r0: float, r1: float, sigma0: float, sigma1: float):
        Q = np.array([
            [-sigma0, sigma0],
            [sigma1, -sigma1]
        ])
        super().__init__(Q, np.array([r0, r1]))
        self._name = 'MMDP2'

        self._r0 = r0
        self._r1 = r1
        self._sigma0 = sigma0
        self._sigma1 = sigma1

    @property
    def r0(self) -> float:
        """Get deterministic rate in state 0."""
        return self._r0

    @property
    def r1(self) -> float:
        """Get deterministic rate in state 1."""
        return self._r1

    def getMeanRate(self) -> float:
        """Stationary mean rate (closed form). Mirrors MATLAB MMDP2.getMeanRate."""
        return ((self._r0 * self._sigma1 + self._r1 * self._sigma0)
                / (self._sigma0 + self._sigma1))

    # getMean is inherited from MMDP (1/getMeanRate), as in MATLAB MMDP2.m,
    # which overrides only the two closed forms below.

    def getSCV(self) -> float:
        """SCV of the modulated rate. Mirrors MATLAB MMDP2.getSCV."""
        pi0 = self._sigma1 / (self._sigma0 + self._sigma1)
        pi1 = self._sigma0 / (self._sigma0 + self._sigma1)
        mean_rate = pi0 * self._r0 + pi1 * self._r1
        var_rate = pi0 * self._r0 ** 2 + pi1 * self._r1 ** 2 - mean_rate ** 2
        if mean_rate > 0:
            return var_rate / (mean_rate ** 2)
        return float('nan')

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

