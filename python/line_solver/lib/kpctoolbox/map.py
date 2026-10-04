"""
Markovian Arrival Process (MAP) functions for KPC-Toolbox.

Native Python implementations of MAP distribution conversion, fitting,
approximation, and sampling functions.
"""

import numpy as np
from scipy import linalg as la
from typing import Tuple, List, Union, Callable, Optional

from line_solver.api.mam.map_analysis import (
    map_pie,
    map_mean,
    map_var,
    map_cdf,
    map_feastol,
    map_erlang,
)

__all__ = [
    'map2ph',
    'map_mmpp2',
    'map_bernstein',
    'me_sample',
    'rap_sample',
]


def map2ph(
    MAP: Union[List[np.ndarray], Tuple[np.ndarray, np.ndarray]]
) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray]]:
    """
    Convert MAP {D0, D1} to PH distribution (alpha, T).

    Extracts the Phase-Type representation from a Markovian Arrival Process.
    The subgenerator T is D0, and the initial probability vector alpha
    is the equilibrium distribution of the embedded DTMC.

    Args:
        MAP: MAP representation as a list/tuple [D0, D1]

    Returns:
        Tuple of (alpha, T, PHR) where:
            alpha: Initial probability vector of the PH distribution
            T: Subgenerator matrix (= D0)
            PHR: PH-renewal process [D0, D1_renewal] in MAP notation
    """
    D0 = np.asarray(MAP[0], dtype=np.float64)
    D1 = np.asarray(MAP[1], dtype=np.float64)

    T = D0
    alpha = map_pie(D0, D1)

    n = D1.shape[0]
    ones_vec = np.ones((n, 1))
    # PHR{2} = D1 * ones(n,1) * alpha
    # D1 * ones = column vector of row sums; then outer product with alpha
    PHR_D1 = (D1 @ ones_vec) @ alpha.reshape(1, -1)
    PHR = [D0, PHR_D1]

    return alpha, T, PHR


def map_mmpp2(
    MEAN: float, SCV: float, SKEW: float, ACF1: float
) -> List[np.ndarray]:
    """
    Fit an MMPP(2) as a MAP given moments and autocorrelation.

    Args:
        MEAN: Mean inter-arrival time of the process
        SCV: Squared coefficient of variation of inter-arrival times
        SKEW: Skewness of inter-arrival times (-1 => automatic minimization,
              applies only to SCV > 1)
        ACF1: Lag-1 autocorrelation coefficient (-1 => maximum feasible
              autocorrelation)

    Returns:
        MAP as a list [D0, D1]

    Examples:
        >>> MAP = map_mmpp2(1, 2, -1, 0.2)  # MMPP(2) with minimal skewness
        >>> MAP = map_mmpp2(1, 2, -1, -1)   # MMPP(2) with minimal skewness
                                             # and maximal autocorrelation
    """
    SCV_REQ = SCV  # keep the request for diagnostics, SCV is recomputed below
    FEASTOL = 10**(-map_feastol())

    E1 = MEAN
    E2 = (1 + SCV) * E1**2
    E3 = -(2 * E1**3 - 3 * E1 * E2 - SKEW * (E2 - E1**2)**(3.0 / 2.0))

    # The closed form below solves the moment-matching equations, but nothing
    # in it constrains the solution to be a MAP: outside the MMPP(2) feasible
    # set it returns negative rates, i.e. a D1 with negative entries and a D0
    # with a positive diagonal. Rowsums stay zero, so the usual generator check
    # does not catch it. Reject the request instead of returning a non-MAP.
    if SCV < 1 - FEASTOL:
        raise ValueError(
            "map_mmpp2: SCV=%g is infeasible, the inter-arrival times of an "
            "MMPP(2) are over-dispersed (SCV>=1)." % SCV_REQ)
    if abs(SCV - 1) <= FEASTOL:
        raise ValueError(
            "map_mmpp2: SCV=1 is the Poisson boundary, where the MMPP(2) fit "
            "is degenerate: the decay rate G2=ACF1/(1-1/SCV)/0.5 divides by "
            "zero and every rate comes back NaN. Use map_exponential(%g) for "
            "a Poisson process." % MEAN)

    # ACF1=RHO0MAX is attained in the limit of a decay rate G2->1; no MMPP(2)
    # exceeds it, and none is negatively autocorrelated
    RHO0MAX = 0.5 * (1 - 1 / SCV)
    if ACF1 != -1:
        if ACF1 < -FEASTOL:
            raise ValueError(
                "map_mmpp2: ACF1=%g is infeasible, an MMPP(2) cannot be "
                "negatively autocorrelated. Pass ACF1=-1 to request the "
                "maximum feasible autocorrelation." % ACF1)
        if ACF1 > RHO0MAX + FEASTOL:
            raise ValueError(
                "map_mmpp2: ACF1=%g exceeds the maximum lag-1 autocorrelation "
                "%g feasible at SCV=%g. Pass ACF1=-1 to request it."
                % (ACF1, RHO0MAX, SCV_REQ))

    if ACF1 == -1:
        G2 = 1 - 10 * 10**(-map_feastol())  # autocorrelation decay rate
    else:
        G2 = ACF1 / (1 - 1 / SCV) / 0.5  # autocorrelation decay rate

    if SKEW == -1 and SCV > 1:
        # determine MAP with nearly minimum third moment
        E3 = (3.0 / 2.0 + 0.001) * E2**2 / E1

    # the SKEW==-1 branch above sits just above E3MIN, which is the infimum of
    # the third moment over the class
    E3MIN = (3.0 / 2.0) * E2**2 / E1
    if E3 < E3MIN - FEASTOL:
        raise ValueError(
            "map_mmpp2: SKEW=%g gives E3=%g, below the minimum third moment %g "
            "feasible at SCV=%g. Pass SKEW=-1 for the minimum-skewness fit."
            % (SKEW, E3, E3MIN, SCV_REQ))

    SCV = (E2 - E1**2) / E1**2

    # Precompute powers for readability
    E1_3 = E1**3
    E1_6 = E1**6
    SCV2 = SCV**2
    SCV3 = SCV**3
    G2_2 = G2**2
    G2_3 = G2**3

    if G2 < 1e-6:
        mu00 = (2 * (6 * E1_3 * SCV - E3) / E1
                / (6 * E1_3 * SCV + 3 * E1_3 * SCV2 + 3 * E1_3 - 2 * E3))
        mu11 = 0
        q01 = (9 * E1**5 * (SCV - 1) * (SCV2 - 2 * SCV + 1)
               / (6 * E1_3 * SCV - E3)
               / (6 * E1_3 * SCV + 3 * E1_3 * SCV2 + 3 * E1_3 - 2 * E3))
        q10 = -3 * (SCV - 1) * E1**2 / (6 * E1_3 * SCV - E3)
    else:
        # Discriminant under the square root
        DISC = (E3**2
                - 12 * E1_3 * SCV * E3
                + 6 * E1_3 * G2 * E3
                - 6 * G2 * SCV * E1_3 * E3
                + 18 * G2 * SCV3 * E1_6
                - 18 * E1_6 * G2 * SCV2
                + 9 * E1_6 * G2_2
                + 36 * E1_6 * SCV2
                + 18 * E1_6 * G2 * SCV
                - 18 * E1_6 * SCV * G2_2
                + 9 * E1_6 * SCV2 * G2_2
                - 18 * E1_6 * G2)
        sqD = np.sqrt(DISC)

        # The repeated subexpression F = A / B where:
        #   A = (-3*E1^3*G2 + 3*E1^3*G2*SCV - 6*E1^3*SCV + E3 + sqD)
        #   B = (-3*E1^3*SCV^2 - 6*E1^3*SCV - 3*E1^3 + 2*E3)
        A = (-3 * E1_3 * G2
             + 3 * E1_3 * G2 * SCV
             - 6 * E1_3 * SCV
             + E3
             + sqD)
        B = (-3 * E1_3 * SCV2
             - 6 * E1_3 * SCV
             - 3 * E1_3
             + 2 * E3)
        F = A / B

        # mu11 = F / E1
        mu11 = F / E1

        # mu00: line 39 of MATLAB, with F replacing every occurrence of A/B
        # mu00 = G2 * NUMER / DENOM / E1
        # where NUMER and DENOM are polynomials in F, E1, E3, SCV, G2
        mu00_numer = (
            - 4 * E3 * G2
            + 4 * F * E3 * G2
            - 18 * E1_3 * F * G2
            - 18 * E1_3 * F * G2 * SCV2
            - 12 * E1_3 * G2_2
            - 12 * E1_3 * F * G2_2 * SCV
            + 12 * E1_3 * F * G2 * SCV
            + 12 * E1_3 * G2 * SCV2
            - 9 * E1_3 * F * SCV
            + 3 * E1_3 * F
            + 12 * E1_3 * G2_2 * SCV
            + 9 * E1_3 * F * SCV2
            + 12 * E1_3 * G2
            + 12 * E1_3 * F * G2_2
            - 3 * E1_3 * F * SCV3
        )
        mu00_denom = (
            12 * E1_3 * G2_3 * SCV
            + 3 * E1_3 * SCV3 * G2
            - 12 * E1_3 * G2_3
            + 18 * E1_3 * G2_2 * SCV2
            - 3 * E1_3 * G2
            + 27 * E1_3 * F * G2 * SCV2
            - 9 * E1_3 * G2 * SCV2
            + 18 * E1_3 * G2_2
            - 12 * E1_3 * G2_2 * SCV
            + 9 * E1_3 * G2 * SCV
            - 12 * E1_3 * F * G2_3 * SCV
            - 9 * E1_3 * F * SCV3 * G2
            - 24 * E1_3 * F * G2_2 * SCV2
            - F * E3 * SCV2
            + 4 * F * E3 * G2_2
            + 12 * E1_3 * F * G2_3
            - F * E3
            + 2 * F * E3 * SCV
            + 9 * E1_3 * F * G2
            + 24 * E1_3 * F * G2_2 * SCV
            - 27 * E1_3 * F * G2 * SCV
            + 6 * E1_3 * F * SCV
            - 12 * E1_3 * F * SCV2
            - 24 * E1_3 * F * G2_2
            + 6 * E1_3 * F * SCV3
            # MATLAB closes the denominator with this term (...*SCV^3-4*E3*G2^2)/E1;
            # it used to sit in mu00_numer here, which biased mu00 by ~0.7%
            - 4 * E3 * G2_2
        )
        mu00 = G2 * mu00_numer / mu00_denom / E1

        # q01: line 41 of MATLAB, with F replacing A/B
        # q01 = -3 * E1^2 * q01_inner / q01_denom
        # First build the inner expression (everything inside the outer
        # parentheses that multiplies -3*E1^2):
        q01_inner = (
            - 6 * F * E1**2 * SCV
            + 12 * F * E1**2 * G2 * SCV
            - 6 * G2 * SCV * E1**2
            - 3 * F * E1**2 * G2
            + F / E1 * E3
            + 3 * E1**2 * G2
            + 6 * F * E1**2 * SCV2
            - 9 * F * E1**2 * SCV2 * G2
            + 3 * E1**2 * G2 * SCV2
            - E3 * F / E1 * SCV
            - 6 * F * E1**2 * G2_2 * SCV
            + 6 * E1**2 * G2_2 * SCV
            + 3 * F * E1**2 * G2_2
            - G2 * F / E1 * E3
            - 3 * E1**2 * G2_2
            + 3 * F * E1**2 * SCV2 * G2_2
            - 3 * E1**2 * SCV2 * G2_2
            + G2 * SCV * F / E1 * E3
        )

        # q01 denominator from MATLAB line 41
        q01_denom = (
            - 45 * F * E1**5 * G2 * SCV2
            + 18 * G2_2 * E1**5 * SCV
            + 18 * E1**5 * G2_3
            - 27 * E1**5 * G2_2 * SCV2
            + 6 * E1**2 * G2_2 * E3
            - 27 * E1**5 * G2_2
            - 18 * E1**5 * G2_3 * SCV
            - 18 * E1**5 * G2 * SCV
            + 18 * E1**5 * G2 * SCV2
            + 3 * E1**2 * G2 * E3
            - 3 * E1**2 * G2 * E3 * SCV
            + F / E1 * E3**2
            + 3 * F * E1**2 * G2 * SCV * E3
            - 36 * F * E1**5 * G2_2 * SCV
            + 36 * F * E1**5 * G2_2
            + 36 * F * E1**5 * SCV2
            + 45 * F * E1**5 * G2 * SCV
            - 12 * F * E1**2 * SCV * E3
            - 3 * F * E1**2 * G2 * E3
            + 9 * F * E1**5 * G2 * SCV3
            + 36 * F * E1**5 * G2_2 * SCV2
            - 6 * F * E1**2 * G2_2 * E3
            + 18 * F * E1**5 * G2_3 * SCV
            - 18 * F * E1**5 * G2_3
            - 9 * F * E1**5 * G2
        )

        q01 = -3 * E1**2 * q01_inner / q01_denom

        # q10: line 42 of MATLAB, with F replacing A/B
        # q10 = 3 * q10_inner * E1^2 * (-1 + G2) / DISC
        q10_inner = (
            - 3 * E1_3 * F * SCV3
            - 3 * E1_3 * F * G2 * SCV2
            + 6 * E1_3 * SCV2
            + 3 * E1_3 * G2 * SCV2
            + 3 * E1_3 * F * SCV2
            + 6 * E1_3 * F * G2 * SCV
            - E3 * SCV
            - 6 * E1_3 * SCV
            + F * E3 * SCV
            - 6 * E1_3 * G2 * SCV
            - 3 * E1_3 * F * SCV
            - F * E3
            + 3 * E1_3 * G2
            - 3 * E1_3 * F * G2
            + 3 * E1_3 * F
            + E3
        )

        q10 = 3 * q10_inner * E1**2 * (-1 + G2) / DISC

    # Catch-all: the checks above cover the known infeasible directions, but
    # the authoritative test is the solution itself. An MMPP(2) has
    # non-negative rates by construction, so anything else is not a MAP and
    # must not be returned.
    rates = np.array([mu00, mu11, q01, q10])
    if np.iscomplexobj(rates) or np.any(np.isnan(rates)) or np.any(np.real(rates) < -FEASTOL):
        raise ValueError(
            "map_mmpp2: (MEAN=%g,SCV=%g,SKEW=%g,ACF1=%g) is not "
            "MMPP(2)-feasible: the fit gives [mu00 mu11 q01 q10]=%s, which is "
            "not a MAP." % (MEAN, SCV_REQ, SKEW, ACF1, np.array2string(rates)))
    # a request on the feasibility boundary lands on a zero rate up to
    # roundoff: clear that sign flip, but leave small POSITIVE rates alone --
    # ACF1=-1 asks for G2=1-1e-7, whose near-uncoupled chain has legitimate
    # rates around 1e-9
    rates[rates < 0] = 0
    mu00, mu11, q01, q10 = rates

    D0 = np.array([[-mu00 - q01, q01],
                    [q10, -mu11 - q10]])
    D1 = np.array([[mu00, 0],
                    [0, mu11]])
    return [D0, D1]


def map_bernstein(
    f: Callable[[float], float], n: int = 20
) -> List[np.ndarray]:
    """
    Bernstein polynomial approximation to convert PDF to MAP.

    Approximates a continuous distribution (specified by its PDF) as a
    Markovian Arrival Process using Bernstein polynomial basis.

    Args:
        f: PDF function handle f(x) - probability density function
        n: Number of phases for the approximation (default: 20)

    Returns:
        MAP as a list [D0, D1]

    Note:
        Caller must rescale to target mean using map_scale.

    Examples:
        >>> from scipy.stats import gamma
        >>> pdf_func = lambda x: gamma.pdf(x, a=2, scale=1)
        >>> MAP = map_bernstein(pdf_func, 20)
    """
    from line_solver.api.mam.map_analysis import map_bernstein as _map_bernstein
    D0, D1 = _map_bernstein(f, n)
    return [D0, D1]


def me_sample(
    ME: Union[List[np.ndarray], Tuple[np.ndarray, np.ndarray]],
    n: int = 1,
    xs: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Generate random samples from a Matrix Exponential (ME) distribution.

    Uses inverse CDF interpolation: computes CDF on a fine grid and uses
    linear interpolation to invert uniform random variables.

    Args:
        ME: ME distribution as a list/tuple [D0, D1]
        n: Number of samples to generate (default: 1)
        xs: Optional pre-computed grid for CDF evaluation.
            If not provided, auto-generates grid based on mean and variance.

    Returns:
        Column vector (n,) of samples from the ME distribution
    """
    D0 = np.asarray(ME[0], dtype=np.float64)
    D1 = np.asarray(ME[1], dtype=np.float64)

    # Auto-generate grid if not provided
    if xs is None:
        mean_val = map_mean(D0, D1)
        var_val = map_var(D0, D1)
        std_val = np.sqrt(max(0, var_val))

        # Create grid from 0 to mean + 10*sigma with 1000 points
        xs = np.linspace(0, mean_val + 10 * std_val, 1000)

    xs = np.asarray(xs, dtype=np.float64).ravel()

    # Compute CDF at grid points
    Fxs = map_cdf(D0, D1, xs)

    # Ensure CDF is strictly increasing for interpolation
    for i in range(1, len(Fxs)):
        if Fxs[i] <= Fxs[i - 1]:
            Fxs[i] = Fxs[i - 1] + np.finfo(float).eps

    # Generate samples via inverse CDF interpolation
    sample = np.zeros(n)
    for i in range(n):
        u = np.random.rand()

        if u <= Fxs[0]:
            sample[i] = xs[0]
        elif u >= Fxs[-1]:
            # Extrapolate beyond grid using exponential tail approximation
            sample[i] = xs[-1] + np.log(1.0 / (1.0 - u + np.finfo(float).eps))
        else:
            # Linear interpolation between grid points
            sample[i] = np.interp(u, Fxs, xs)

    return sample


def rap_sample(
    RAP: Union[List[np.ndarray], Tuple[np.ndarray, np.ndarray]],
    n: int = 1,
    xs: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Generate random samples from a RAP (Rational Arrival Process) distribution.

    The marginal distribution of RAP inter-arrival times is a Matrix
    Exponential (ME) distribution, so this function delegates to me_sample.

    Args:
        RAP: RAP distribution as a list/tuple [H0, H1]
        n: Number of samples to generate (default: 1)
        xs: Optional pre-computed grid for CDF evaluation.
            If not provided, auto-generates grid based on mean and variance.

    Returns:
        Column vector (n,) of samples from the RAP marginal distribution

    Note:
        This function generates samples from the marginal distribution only.
        It does not preserve the correlation structure of the RAP.
    """
    # Delegate to me_sample (marginal of RAP is ME with same {H0, H1})
    if xs is not None:
        return me_sample(RAP, n, xs)
    else:
        return me_sample(RAP, n)
