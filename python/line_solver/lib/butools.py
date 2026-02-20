"""
BuTools library functions - comprehensive toolkit for:
- Phase-type (PH) and matrix-exponential (ME) distributions
- Markov chains (CTMC/DTMC) analysis
- Moment transformations (raw, factorial, Hankel, normalized, reduced)
- MAP/MMAP validation and analysis
- Representation transformations
- APH/DPH analysis and bounds
- Discrete MAP (DMAP) functions
"""

import numpy as np


# ========== PHASE-TYPE (PH) FUNCTIONS ==========

def lib_butools_ph_from_moments(moments):
    """Create matrix-exponential representation from moments."""
    from line_solver.api.butools.ph.appie import MEFromMoments
    moments = list(np.asarray(moments).flatten())
    alpha, A = MEFromMoments(moments)
    return {'D0': np.asarray(A), 'D1': np.asarray(alpha)}


def lib_butools_ph_moments(alpha, A):
    """Compute moments from PH distribution."""
    from line_solver.api.butools.ph.baseph import MomentsFromPH
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    A = np.asarray(A, dtype=float)
    result = MomentsFromPH(alpha, A)
    return np.asarray(result)


def lib_butools_ph_me_moments(D0, D1):
    """Compute moments from ME distribution (alpha=D1, A=D0)."""
    from line_solver.api.butools.ph.baseph import MomentsFromME
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    result = MomentsFromME(D1, D0)
    return np.asarray(result)


def lib_butools_ph2_from_moments(m1, m2, m3):
    """Fit 2-phase APH from 3 moments."""
    from line_solver.api.butools.ph.canonical import PH2From3Moments
    alpha, A = PH2From3Moments([float(m1), float(m2), float(m3)])
    return {'alpha': np.asarray(alpha), 'A': np.asarray(A)}


def lib_butools_ph3_from_moments(moments):
    """Fit 3-phase APH from 5 moments."""
    from line_solver.api.butools.ph.canonical import PH3From5Moments
    moments = list(np.asarray(moments).flatten())
    alpha, A = PH3From5Moments(moments)
    return {'alpha': np.asarray(alpha), 'A': np.asarray(A)}


def lib_butools_aph_bounds(m1, m2, m3, order=2):
    """Compute moment bounds for APH of given order."""
    from line_solver.api.butools.ph.canonical import (
        APH2ndMomentLowerBound, APH3rdMomentLowerBound, APH3rdMomentUpperBound
    )
    if order == 2:
        lb2 = APH2ndMomentLowerBound(float(m1), 2)
        lb3 = APH3rdMomentLowerBound(float(m1), float(m2), 2)
        ub3 = APH3rdMomentUpperBound(float(m1), float(m2), 2)
        return {'m2_lower': float(lb2), 'm3_lower': float(lb3), 'm3_upper': float(ub3)}
    else:
        raise ValueError("Bounds implemented for order 2 only")


def lib_butools_ph_cdf(alpha, A, x):
    """Compute CDF of PH distribution."""
    from line_solver.api.butools.ph.baseph import CdfFromPH
    alpha = np.matrix(np.atleast_2d(np.asarray(alpha, dtype=float)))
    A = np.matrix(np.asarray(A, dtype=float))
    x = np.atleast_1d(np.asarray(x, dtype=float))
    result = CdfFromPH(alpha, A, x)
    return np.asarray(result)


def lib_butools_ph_pdf(alpha, A, x):
    """Compute PDF of PH distribution."""
    from line_solver.api.butools.ph.baseph import PdfFromPH
    alpha = np.matrix(np.atleast_2d(np.asarray(alpha, dtype=float)))
    A = np.matrix(np.asarray(A, dtype=float))
    x = np.atleast_1d(np.asarray(x, dtype=float))
    result = PdfFromPH(alpha, A, x)
    return np.asarray(result)


def lib_butools_ph_canonical2(alpha, A):
    """Convert 2-phase PH to canonical form."""
    from line_solver.api.butools.ph.canonical import CanonicalFromPH2
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    A = np.asarray(A, dtype=float)
    alpha_c, A_c = CanonicalFromPH2(alpha, A)
    return {'alpha': np.asarray(alpha_c), 'A': np.asarray(A_c)}


def lib_butools_ph_check(alpha, A):
    """Validate PH representation."""
    from line_solver.api.butools.ph.check import CheckPHRepresentation
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    A = np.asarray(A, dtype=float)
    return bool(CheckPHRepresentation(alpha, A))


# ========== MOMENT TRANSFORMATIONS ==========

def lib_butools_moment_factorial_from_raw(moments):
    """Convert raw moments to factorial moments."""
    from line_solver.api.butools.moments.conv import FactorialMomsFromMoms
    moments = list(np.asarray(moments).flatten())
    result = FactorialMomsFromMoms(moments)
    return np.asarray(result)


def lib_butools_moment_raw_from_factorial(moments):
    """Convert factorial moments to raw moments."""
    from line_solver.api.butools.moments.conv import MomsFromFactorialMoms
    moments = list(np.asarray(moments).flatten())
    result = MomsFromFactorialMoms(moments)
    return np.asarray(result)


def lib_butools_moment_hankel_from_raw(moments):
    """Convert raw moments to Hankel moments."""
    from line_solver.api.butools.moments.conv import HankelMomsFromMoms
    moments = list(np.asarray(moments).flatten())
    result = HankelMomsFromMoms(moments)
    return np.asarray(result)


def lib_butools_moment_normalized_from_raw(moments):
    """Convert raw moments to normalized moments."""
    from line_solver.api.butools.moments.conv import NormMomsFromMoms
    moments = list(np.asarray(moments).flatten())
    result = NormMomsFromMoms(moments)
    return np.asarray(result)


def lib_butools_moment_reduced_from_raw(moments):
    """Convert raw moments to reduced moments."""
    from line_solver.api.butools.moments.conv import ReducedMomsFromMoms
    moments = list(np.asarray(moments).flatten())
    result = ReducedMomsFromMoms(moments)
    return np.asarray(result)


def lib_butools_moment_joint_factorial(moments):
    """Convert joint moments to joint factorial moments."""
    from line_solver.api.butools.moments.conv import JFactorialMomsFromJMoms
    moments = list(np.asarray(moments).flatten())
    result = JFactorialMomsFromJMoms(moments)
    return np.asarray(result)


def lib_butools_moments_from_moms(moments):
    """Convert raw moments to factorial moments (alias)."""
    return lib_butools_moment_factorial_from_raw(moments)


def lib_butools_me_from_moments(moments):
    """Create ME representation from moments (returns alpha, A)."""
    from line_solver.api.butools.ph.appie import MEFromMoments
    moments = list(np.asarray(moments).flatten())
    alpha, A = MEFromMoments(moments)
    return np.asarray(alpha), np.asarray(A)


# ========== MAP/MMAP FUNCTIONS ==========

def lib_butools_map_check(D0, D1):
    """Validate MAP representation."""
    from line_solver.api.butools.map.check import CheckMAPRepresentation
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    return bool(CheckMAPRepresentation(D0, D1))


def lib_butools_map_marginal_moments(D0, D1):
    """Compute marginal moments from MAP."""
    from line_solver.api.butools.map.basemap import MarginalMomentsFromMAP
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    result = MarginalMomentsFromMAP(D0, D1)
    return np.asarray(result)


def lib_butools_map_lag_correlations(D0, D1, L=1):
    """Compute lag-k autocorrelations from MAP."""
    from line_solver.api.butools.map.basemap import LagCorrelationsFromMAP
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    result = LagCorrelationsFromMAP(D0, D1, int(L))
    return np.asarray(result)


# ========== MARKOV CHAIN FUNCTIONS ==========

def lib_butools_markov_check_generator(Q):
    """Validate CTMC generator matrix."""
    from line_solver.api.butools.mc.check import CheckGenerator
    Q = np.asarray(Q, dtype=float)
    return bool(CheckGenerator(Q))


def lib_butools_markov_check_probability_matrix(P):
    """Validate DTMC transition matrix."""
    from line_solver.api.butools.mc.check import CheckProbMatrix
    P = np.asarray(P, dtype=float)
    return bool(CheckProbMatrix(P))


def lib_butools_markov_ctmc_solve(Q):
    """Solve CTMC for stationary distribution."""
    from line_solver.api.butools.mc.stst import CTMCSolve
    Q = np.asarray(Q, dtype=float)
    result = CTMCSolve(Q)
    return np.asarray(result)


def lib_butools_markov_dtmc_solve(P):
    """Solve DTMC for stationary distribution."""
    from line_solver.api.butools.mc.stst import DTMCSolve
    P = np.asarray(P, dtype=float)
    result = DTMCSolve(P)
    return np.asarray(result)


def lib_butools_markov_crp_solve(Q):
    """Solve continuous renewal process."""
    from line_solver.api.butools.mc.stst import CRPSolve
    Q = np.asarray(Q, dtype=float)
    result = CRPSolve(Q)
    return np.asarray(result)


def lib_butools_markov_drp_solve(P):
    """Solve discrete renewal process."""
    from line_solver.api.butools.mc.stst import DRPSolve
    P = np.asarray(P, dtype=float)
    result = DRPSolve(P)
    return np.asarray(result)


# ========== UTILITY FUNCTIONS ==========

def lib_butools_check_moments(moments):
    """Validate moment sequence feasibility."""
    from line_solver.api.butools.moments.check import CheckMoments
    moments = list(np.asarray(moments).flatten())
    return bool(CheckMoments(moments))


# ========== DISCRETE PHASE-TYPE (DPH) FUNCTIONS ==========

def lib_butools_dph2_from_moments(m1, m2, m3):
    """Fit 2-phase discrete PH from 3 moments."""
    from line_solver.api.butools.dph.canonical import DPH2From3Moments
    sigma, D = DPH2From3Moments([float(m1), float(m2), float(m3)])
    return {'sigma': np.asarray(sigma), 'D': np.asarray(D)}


def lib_butools_dph3_from_moments(moments):
    """Fit 3-phase discrete PH from 5 moments."""
    from line_solver.api.butools.dph.canonical import DPH3From5Moments
    moments = list(np.asarray(moments).flatten())
    sigma, D = DPH3From5Moments(moments)
    return {'sigma': np.asarray(sigma), 'D': np.asarray(D)}


def lib_butools_dph_canonical_form(sigma, D):
    """Convert discrete PH to canonical form."""
    from line_solver.api.butools.dph.canonical import CanonicalFromDPH2
    sigma = np.atleast_2d(np.asarray(sigma, dtype=float))
    D = np.asarray(D, dtype=float)
    sigma_c, D_c = CanonicalFromDPH2(sigma, D)
    return {'sigma': np.asarray(sigma_c), 'D': np.asarray(D_c)}


def lib_butools_dph_pmf(sigma, D, n):
    """Compute probability mass function for discrete PH."""
    from line_solver.api.butools.dph.basedph import PmfFromDPH
    sigma = np.atleast_2d(np.asarray(sigma, dtype=float))
    D = np.asarray(D, dtype=float)
    x = np.arange(1, int(n) + 1)
    result = PmfFromDPH(sigma, D, x)
    return np.asarray(result)


def lib_butools_dph_moments(sigma, D, num_moments=3):
    """Compute moments from discrete PH representation."""
    from line_solver.api.butools.dph.basedph import MomentsFromDPH
    sigma = np.atleast_2d(np.asarray(sigma, dtype=float))
    D = np.asarray(D, dtype=float)
    result = MomentsFromDPH(sigma, D, int(num_moments))
    return np.asarray(result)


# ========== REPRESENTATION TRANSFORMATIONS ==========

def lib_butools_transform_ph_to_me(alpha, A):
    """Transform PH representation to matrix-exponential form.
    PH is already a valid ME, so this is an identity transform."""
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    A = np.asarray(A, dtype=float)
    return {'D0': np.array(A), 'D1': np.array(alpha)}


def lib_butools_transform_me_to_ph(alpha, A):
    """Transform matrix-exponential representation to PH form."""
    from line_solver.api.butools.ph.monocyclic import MonocyclicPHFromME
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    A = np.asarray(A, dtype=float)
    beta, B = MonocyclicPHFromME(alpha, A)
    return {'alpha': np.asarray(beta), 'A': np.asarray(B)}


def lib_butools_transform_ph_to_dph(alpha, A):
    """Transform continuous PH to discrete PH representation via uniformization."""
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    A = np.asarray(A, dtype=float)
    max_rate = np.max(-np.diag(A))
    D = np.eye(A.shape[0]) + A / max_rate
    return {'sigma': np.array(alpha), 'D': np.array(D)}


def lib_butools_transform_dph_to_ph(sigma, D):
    """Transform discrete PH to continuous PH representation."""
    sigma = np.atleast_2d(np.asarray(sigma, dtype=float))
    D = np.asarray(D, dtype=float)
    import scipy.linalg as la
    A = la.logm(D)
    return {'alpha': np.array(sigma), 'A': np.real(np.array(A))}


# ========== DISCRETE MAP (DMAP) FUNCTIONS ==========

def lib_butools_dmap_check(D0, D1, prec=1e-14):
    """Validate discrete MAP representation."""
    from line_solver.api.butools.dmap.check import CheckDMAPRepresentation
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    return bool(CheckDMAPRepresentation(D0, D1, prec))


def lib_butools_drap_check(H0, H1, prec=1e-14):
    """Validate discrete rational arrival process representation."""
    from line_solver.api.butools.dmap.check import CheckDRAPRepresentation
    H0 = np.asarray(H0, dtype=float)
    H1 = np.asarray(H1, dtype=float)
    return bool(CheckDRAPRepresentation(H0, H1, prec))


def lib_butools_dmap_marginal_moments(D0, D1, K=0):
    """Compute marginal moments from discrete MAP."""
    from line_solver.api.butools.dmap.basedmap import MarginalMomentsFromDMAP
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    result = MarginalMomentsFromDMAP(D0, D1, int(K))
    return np.asarray(result)


def lib_butools_drap_marginal_moments(H0, H1, K=0):
    """Compute marginal moments from discrete rational arrival process."""
    from line_solver.api.butools.dmap.basedmap import MarginalMomentsFromDRAP
    H0 = np.asarray(H0, dtype=float)
    H1 = np.asarray(H1, dtype=float)
    result = MarginalMomentsFromDRAP(H0, H1, int(K))
    return np.asarray(result)


def lib_butools_dmap_lag_correlations(D0, D1, L=1):
    """Compute lag autocorrelations from discrete MAP."""
    from line_solver.api.butools.dmap.basedmap import LagCorrelationsFromDMAP
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    result = LagCorrelationsFromDMAP(D0, D1, int(L))
    return np.asarray(result)


def lib_butools_drap_lag_correlations(H0, H1, L=1):
    """Compute lag autocorrelations from discrete rational arrival process."""
    from line_solver.api.butools.dmap.basedmap import LagCorrelationsFromDRAP
    H0 = np.asarray(H0, dtype=float)
    H1 = np.asarray(H1, dtype=float)
    result = LagCorrelationsFromDRAP(H0, H1, int(L))
    return np.asarray(result)


def lib_butools_dmap_marginal_distribution(D0, D1):
    """Get marginal distribution from discrete MAP."""
    from line_solver.api.butools.dmap.basedmap import MarginalDistributionFromDMAP
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    alpha, A = MarginalDistributionFromDMAP(D0, D1)
    return {'alpha': np.asarray(alpha), 'A': np.asarray(A)}


def lib_butools_dmap_samples(D0, D1, K):
    """Generate random samples from discrete MAP."""
    from line_solver.api.butools.dmap.misc import SamplesFromDMAP
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    result = SamplesFromDMAP(D0, D1, int(K))
    return np.asarray(result)


def lib_butools_dmap_from_drap(H0, H1, prec=1e-14):
    """Convert discrete rational arrival process to discrete MAP."""
    from line_solver.api.butools.dmap.dmapfromdrap import DMAPFromDRAP
    H0 = np.asarray(H0, dtype=float)
    H1 = np.asarray(H1, dtype=float)
    D0, D1 = DMAPFromDRAP(H0, H1, prec)
    return {'D0': np.asarray(D0), 'D1': np.asarray(D1)}


def lib_butools_dmap2_from_moments(moms, corr1):
    """Create discrete MAP(2) from 3 moments and lag-1 autocorrelation."""
    from line_solver.api.butools.dmap.matching import DMAP2FromMoments
    moms = list(np.asarray(moms).flatten())
    D0, D1 = DMAP2FromMoments(moms, float(corr1))
    return {'D0': np.asarray(D0), 'D1': np.asarray(D1)}


def lib_butools_dmap_canonical2(D0, D1, prec=1e-14):
    """Convert order-2 discrete MAP to canonical form."""
    from line_solver.api.butools.dmap.matching import CanonicalFromDMAP2
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    G0, G1 = CanonicalFromDMAP2(D0, D1)
    return {'G0': np.asarray(G0), 'G1': np.asarray(G1)}


def lib_butools_random_dmap(order, mean=10.0, zero_entries=0, max_trials=1000, prec=1e-7):
    """Generate a random discrete MAP."""
    from line_solver.api.butools.dmap.basedmap import RandomDMAP
    D0, D1 = RandomDMAP(int(order), float(mean), int(zero_entries), int(max_trials), float(prec))
    return {'D0': np.asarray(D0), 'D1': np.asarray(D1)}
