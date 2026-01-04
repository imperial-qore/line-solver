"""
BuTools library functions - comprehensive toolkit for:
- Phase-type (PH) and matrix-exponential (ME) distributions
- Markov chains (CTMC/DTMC) analysis
- Moment transformations (raw, factorial, Hankel, normalized, reduced)
- MAP/MMAP validation and analysis
- Fluid/MAM queue solvers
- Representation transformations
- Matrix analysis utilities
- APH/DPH analysis and bounds
- Advanced queue analysis (QMAM)
"""

import numpy as np
from line_solver import native_to_array


# ========== PHASE-TYPE (PH) FUNCTIONS ==========

def lib_butools_ph_from_moments(moments):
    """Create matrix-exponential from moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        D0 = native_to_array(result.first)
        D1 = native_to_array(result.second)
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"meFromMoments failed: {str(e)}")


def lib_butools_ph_moments(alpha, A):
    """Compute moments from PH distribution."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
            native_to_array(alpha), native_to_array(A)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"momentsFromPH failed: {str(e)}")


def lib_butools_ph_me_moments(D0, D1):
    """Compute moments from ME distribution."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
            native_to_array(D0), native_to_array(D1)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"momentsFromME failed: {str(e)}")


def lib_butools_ph2_from_moments(m1, m2, m3):
    """Fit 2-phase APH from 3 moments."""
    try:
            float(m1), float(m2), float(m3)
        )
        alpha = native_to_array(result.first)
        A = native_to_array(result.second)
        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"ph2From3Moments failed: {str(e)}")


def lib_butools_ph3_from_moments(moments):
    """Fit 3-phase APH from 5 moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        alpha = native_to_array(result.first)
        A = native_to_array(result.second)
        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"ph3From5Moments failed: {str(e)}")


def lib_butools_aph_bounds(m1, m2, m3, order=2):
    """Compute moment bounds for APH of given order."""
    try:
        if order == 2:
            return {'m2_lower': float(lb2), 'm3_lower': float(lb3), 'm3_upper': float(ub3)}
        else:
            raise ValueError("Bounds implemented for order 2 only")
    except Exception as e:
        raise RuntimeError(f"APH bounds computation failed: {str(e)}")


def lib_butools_ph_cdf(alpha, A, x):
    """Compute CDF of PH distribution."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    x = np.atleast_1d(x)
    try:
            native_to_array(alpha), native_to_array(A), native_to_array(x)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"cdfFromPH failed: {str(e)}")


def lib_butools_ph_pdf(alpha, A, x):
    """Compute PDF of PH distribution."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    x = np.atleast_1d(x)
    try:
            native_to_array(alpha), native_to_array(A), native_to_array(x)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"pdfFromPH failed: {str(e)}")


def lib_butools_ph_canonical2(alpha, A):
    """Convert 2-phase PH to canonical form."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
            native_to_array(alpha), native_to_array(A)
        )
        alpha_c = native_to_array(result.first)
        A_c = native_to_array(result.second)
        return {'alpha': alpha_c, 'A': A_c}
    except Exception as e:
        raise RuntimeError(f"canonicalFromPH2 failed: {str(e)}")


def lib_butools_ph_check(alpha, A):
    """Validate PH representation."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
            native_to_array(alpha), native_to_array(A)
        )
        return bool(result)
    except Exception as e:
        return False


# ========== MOMENT TRANSFORMATIONS ==========

def lib_butools_moment_factorial_from_raw(moments):
    """Convert raw moments to factorial moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"factorialMomsFromMoms failed: {str(e)}")


def lib_butools_moment_raw_from_factorial(moments):
    """Convert factorial moments to raw moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"momsFromFactorialMoms failed: {str(e)}")


def lib_butools_moment_hankel_from_raw(moments):
    """Convert raw moments to Hankel moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"hankelMomsFromMoms failed: {str(e)}")


def lib_butools_moment_normalized_from_raw(moments):
    """Convert raw moments to normalized moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"normMomsFromMoms failed: {str(e)}")


def lib_butools_moment_reduced_from_raw(moments):
    """Convert raw moments to reduced moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"reducedMomsFromMoms failed: {str(e)}")


def lib_butools_moment_joint_factorial(moments):
    """Convert joint moments to joint factorial moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"jFactorialMomsFromJMoms failed: {str(e)}")


# ========== MAP/MMAP FUNCTIONS ==========

def lib_butools_map_check(D0, D1):
    """Validate MAP representation."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
            native_to_array(D0), native_to_array(D1)
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_map_marginal_moments(D0, D1):
    """Compute marginal moments from MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
            native_to_array(D0), native_to_array(D1)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"marginalMomentsFromMAP failed: {str(e)}")


def lib_butools_map_lag_correlations(D0, D1, lags):
    """Compute lag-k autocorrelations from MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    lags = np.atleast_1d(lags)
    try:
            native_to_array(D0), native_to_array(D1), native_to_array(lags)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"lagCorrelationsFromMAP failed: {str(e)}")


# ========== MARKOV CHAIN FUNCTIONS ==========

def lib_butools_markov_check_generator(Q):
    """Validate CTMC generator matrix."""
    Q = np.asarray(Q)
    try:
            native_to_array(Q)
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_markov_check_probability_matrix(P):
    """Validate DTMC transition matrix."""
    P = np.asarray(P)
    try:
            native_to_array(P)
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_markov_ctmc_solve(Q):
    """Solve CTMC for stationary distribution."""
    Q = np.asarray(Q)
    try:
            native_to_array(Q)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"ctmcSolve failed: {str(e)}")


def lib_butools_markov_dtmc_solve(P):
    """Solve DTMC for stationary distribution."""
    P = np.asarray(P)
    try:
            native_to_array(P)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"dtmcSolve failed: {str(e)}")


def lib_butools_markov_crp_solve(D0, D1):
    """Solve continuous renewal process."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
            native_to_array(D0), native_to_array(D1)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"crpSolve failed: {str(e)}")


def lib_butools_markov_drp_solve(D0, D1):
    """Solve discrete renewal process."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
            native_to_array(D0), native_to_array(D1)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"drpSolve failed: {str(e)}")


# ========== UTILITY FUNCTIONS ==========

def lib_butools_check_moments(moments):
    """Validate moment sequence feasibility."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        return bool(result)
    except Exception as e:
        return False


# ========== DISCRETE PHASE-TYPE (DPH) FUNCTIONS ==========

def lib_butools_dph2_from_moments(m1, m2, m3):
    """Fit 2-phase discrete PH from 3 moments."""
    try:
            float(m1), float(m2), float(m3)
        )
        sigma = native_to_array(result.first)
        D = native_to_array(result.second)
        return {'sigma': sigma, 'D': D}
    except Exception as e:
        raise RuntimeError(f"dph2From3Moments failed: {str(e)}")


def lib_butools_dph3_from_moments(moments):
    """Fit 3-phase discrete PH from 5 moments."""
    moments = np.asarray(moments).flatten()
    try:
            native_to_array(moments)
        )
        sigma = native_to_array(result.first)
        D = native_to_array(result.second)
        return {'sigma': sigma, 'D': D}
    except Exception as e:
        raise RuntimeError(f"dph3From5Moments failed: {str(e)}")


def lib_butools_dph_canonical_form(sigma, D):
    """Convert discrete PH to canonical form."""
    sigma = np.asarray(sigma)
    D = np.asarray(D)
    try:
            native_to_array(sigma), native_to_array(D)
        )
        sigma_c = native_to_array(result.first)
        D_c = native_to_array(result.second)
        return {'sigma': sigma_c, 'D': D_c}
    except Exception as e:
        raise RuntimeError(f"canonicalFromDPH failed: {str(e)}")


def lib_butools_dph_pmf(sigma, D, n):
    """Compute probability mass function for discrete PH."""
    sigma = np.asarray(sigma)
    D = np.asarray(D)
    try:
            native_to_array(sigma), native_to_array(D), int(n)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"pmfFromDPH failed: {str(e)}")


def lib_butools_dph_moments(sigma, D, num_moments=3):
    """Compute moments from discrete PH representation."""
    sigma = np.asarray(sigma)
    D = np.asarray(D)
    try:
            native_to_array(sigma), native_to_array(D), int(num_moments)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"momentsFromDPH failed: {str(e)}")


# ========== REPRESENTATION TRANSFORMATIONS ==========

def lib_butools_transform_ph_to_me(alpha, A):
    """Transform PH representation to matrix-exponential form."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
            native_to_array(alpha), native_to_array(A)
        )
        D0 = native_to_array(result.first)
        D1 = native_to_array(result.second)
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"phToMe transformation failed: {str(e)}")


def lib_butools_transform_me_to_ph(D0, D1):
    """Transform matrix-exponential representation to PH form."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
            native_to_array(D0), native_to_array(D1)
        )
        alpha = native_to_array(result.first)
        A = native_to_array(result.second)
        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"meToPh transformation failed: {str(e)}")


def lib_butools_transform_ph_to_dph(alpha, A):
    """Transform PH to discrete PH representation."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
            native_to_array(alpha), native_to_array(A)
        )
        sigma = native_to_array(result.first)
        D = native_to_array(result.second)
        return {'sigma': sigma, 'D': D}
    except Exception as e:
        raise RuntimeError(f"phToDph transformation failed: {str(e)}")


def lib_butools_transform_dph_to_ph(sigma, D):
    """Transform discrete PH to continuous PH representation."""
    sigma = np.asarray(sigma)
    D = np.asarray(D)
    try:
            native_to_array(sigma), native_to_array(D)
        )
        alpha = native_to_array(result.first)
        A = native_to_array(result.second)
        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"dphToPh transformation failed: {str(e)}")


# ========== DISCRETE MAP (DMAP) FUNCTIONS ==========

def lib_butools_dmap_check(D0, D1, prec=1e-14):
    """Validate discrete MAP representation."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = dmap.CheckDMAPRepresentationKt.checkDMAPRepresentation(
            native_to_array(D0), native_to_array(D1), float(prec)
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_drap_check(H0, H1, prec=1e-14):
    """Validate discrete rational arrival process representation."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        result = dmap.CheckDRAPRepresentationKt.checkDRAPRepresentation(
            native_to_array(H0), native_to_array(H1), float(prec)
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_dmap_marginal_moments(D0, D1, K=0, prec=1e-14):
    """Compute marginal moments from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = dmap.MarginalMomentsFromDMAPKt.marginalMomentsFromDMAP(
            native_to_array(D0), native_to_array(D1), int(K), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"marginalMomentsFromDMAP failed: {str(e)}")


def lib_butools_drap_marginal_moments(H0, H1, K=0, prec=1e-14):
    """Compute marginal moments from discrete rational arrival process."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        result = dmap.MarginalMomentsFromDRAPKt.marginalMomentsFromDRAP(
            native_to_array(H0), native_to_array(H1), int(K), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"marginalMomentsFromDRAP failed: {str(e)}")


def lib_butools_dmap_lag_correlations(D0, D1, L=1, prec=1e-14):
    """Compute lag autocorrelations from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = dmap.LagCorrelationsFromDMAPKt.lagCorrelationsFromDMAP(
            native_to_array(D0), native_to_array(D1), int(L), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"lagCorrelationsFromDMAP failed: {str(e)}")


def lib_butools_drap_lag_correlations(H0, H1, L=1, prec=1e-14):
    """Compute lag autocorrelations from discrete rational arrival process."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        result = dmap.LagCorrelationsFromDRAPKt.lagCorrelationsFromDRAP(
            native_to_array(H0), native_to_array(H1), int(L), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"lagCorrelationsFromDRAP failed: {str(e)}")


def lib_butools_dmap_marginal_distribution(D0, D1, prec=1e-14):
    """Get marginal distribution from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = dmap.MarginalDistributionFromDMAPKt.marginalDistributionFromDMAP(
            native_to_array(D0), native_to_array(D1), float(prec)
        )
        alpha = native_to_array(result.getAlpha())
        A = native_to_array(result.getA())
        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"marginalDistributionFromDMAP failed: {str(e)}")


def lib_butools_dmap_samples(D0, D1, K, prec=1e-14):
    """Generate random samples from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = dmap.SamplesFromDMAPKt.samplesFromDMAP(
            native_to_array(D0), native_to_array(D1), int(K), None, float(prec), None
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"samplesFromDMAP failed: {str(e)}")


def lib_butools_dmap_from_drap(H0, H1, prec=1e-14):
    """Convert discrete rational arrival process to discrete MAP."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        result = dmap.DMAPFromDRAPKt.dmapFromDRAP(
            native_to_array(H0), native_to_array(H1), float(prec)
        )
        D0 = native_to_array(result.getFirst())
        D1 = native_to_array(result.getSecond())
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"dmapFromDRAP failed: {str(e)}")


def lib_butools_dmap2_from_moments(moms, corr1):
    """Create discrete MAP(2) from 3 moments and lag-1 autocorrelation."""
    moms = np.asarray(moms).flatten()
    try:
        result = dmap.DMAP2FromMomentsKt.dmap2FromMoments(
            moms.tolist(), float(corr1)
        )
        D0 = native_to_array(result.getFirst())
        D1 = native_to_array(result.getSecond())
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"dmap2FromMoments failed: {str(e)}")


def lib_butools_dmap_canonical2(D0, D1, prec=1e-14):
    """Convert order-2 discrete MAP to canonical form."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = dmap.CanonicalFromDMAP2Kt.canonicalFromDMAP2(
            native_to_array(D0), native_to_array(D1), float(prec)
        )
        G0 = native_to_array(result.getFirst())
        G1 = native_to_array(result.getSecond())
        return {'G0': G0, 'G1': G1}
    except Exception as e:
        raise RuntimeError(f"canonicalFromDMAP2 failed: {str(e)}")


def lib_butools_random_dmap(order, mean=10.0, zero_entries=0, max_trials=1000, prec=1e-7):
    """Generate a random discrete MAP."""
    try:
        result = dmap.RandomDMAPKt.randomDMAP(
            int(order), float(mean), int(zero_entries), int(max_trials), float(prec), None
        )
        D0 = native_to_array(result.getFirst())
        D1 = native_to_array(result.getSecond())
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"randomDMAP failed: {str(e)}")
