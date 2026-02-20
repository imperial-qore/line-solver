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

import jpype
import numpy as np
from line_solver import jlineMatrixFromArray, jlineMatrixToArray


def _to_double_array(arr):
    """Convert numpy array/list to Java DoubleArray."""
    arr = np.asarray(arr, dtype=float).flatten()
    return jpype.JArray(jpype.JDouble)(arr.tolist())


def _to_int_array(arr):
    """Convert numpy array/list to Java IntArray."""
    arr = np.asarray(arr, dtype=int).flatten()
    return jpype.JArray(jpype.JInt)(arr.tolist())


# ========== PHASE-TYPE (PH) FUNCTIONS ==========

def lib_butools_ph_from_moments(moments):
    """Create matrix-exponential from moments.

    Returns (alpha, A) tuple where alpha is the initial probability vector
    and A is the transient generator matrix.
    """
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.MEFromMomentsKt.meFromMoments(_to_double_array(moments))
        alpha = jlineMatrixToArray(result.getAlpha())
        A = jlineMatrixToArray(result.getA())
        return alpha, A
    except Exception as e:
        raise RuntimeError(f"meFromMoments failed: {str(e)}")


def lib_butools_ph_moments(alpha, A, K=5):
    """Compute moments from PH distribution."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.MomentsFromMEKt.momentsFromPH(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A), int(K)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"momentsFromPH failed: {str(e)}")


def lib_butools_ph_me_moments(alpha, A, K=5):
    """Compute moments from ME distribution."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.MomentsFromMEKt.momentsFromME(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A), int(K)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"momentsFromME failed: {str(e)}")


def lib_butools_ph2_from_moments(m1, m2, m3):
    """Fit 2-phase APH from 3 moments.

    Returns (alpha, A) tuple.
    """
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.PH2From3MomentsKt.ph2From3Moments(
            _to_double_array([float(m1), float(m2), float(m3)]), 1e-14
        )
        alpha = jlineMatrixToArray(result.getAlpha())
        A = jlineMatrixToArray(result.getA())
        return alpha, A
    except Exception as e:
        raise RuntimeError(f"ph2From3Moments failed: {str(e)}")


def lib_butools_ph3_from_moments(moments):
    """Fit 3-phase APH from 5 moments.

    Returns (alpha, A) tuple.
    """
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.PH3From5MomentsKt.ph3From5Moments(_to_double_array(moments), 1e-10)
        alpha = jlineMatrixToArray(result.getAlpha())
        A = jlineMatrixToArray(result.getA())
        return alpha, A
    except Exception as e:
        raise RuntimeError(f"ph3From5Moments failed: {str(e)}")


def lib_butools_aph_bounds(m1, m2, m3=None, order=2):
    """Compute moment bounds for APH of given order."""
    try:
        bt = jpype.JPackage('jline').lib.butools
        lb2 = float(bt.APH2ndMomentLowerBoundKt.APH2ndMomentLowerBound(
            float(m1), int(order)
        ))
        if m2 is not None:
            lb3 = float(bt.APH3rdMomentLowerBoundKt.APH3rdMomentLowerBound(
                float(m1), float(m2), int(order)
            ))
            ub3 = float(bt.APH3rdMomentUpperBoundKt.APH3rdMomentUpperBound(
                float(m1), float(m2), int(order)
            ))
            return {'m2_lower': lb2, 'm3_lower': lb3, 'm3_upper': ub3}
        return {'m2_lower': lb2}
    except Exception as e:
        raise RuntimeError(f"APH bounds computation failed: {str(e)}")


def lib_butools_ph_cdf(alpha, A, x):
    """Compute CDF of PH distribution."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    x = np.atleast_1d(np.asarray(x, dtype=float))
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.CdfFromMEKt.cdfFromPH(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A),
            _to_double_array(x)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"cdfFromPH failed: {str(e)}")


def lib_butools_ph_pdf(alpha, A, x):
    """Compute PDF of PH distribution."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    x = np.atleast_1d(np.asarray(x, dtype=float))
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.PdfFromMEKt.pdfFromPH(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A),
            _to_double_array(x)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"pdfFromPH failed: {str(e)}")


def lib_butools_ph_canonical2(alpha, A):
    """Convert 2-phase PH to canonical form.

    Returns (alpha, A) tuple.
    """
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.CanonicalFromPH2Kt.canonicalFromPH2(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A), 1e-14
        )
        alpha_c = jlineMatrixToArray(result.getAlpha())
        A_c = jlineMatrixToArray(result.getA())
        return alpha_c, A_c
    except Exception as e:
        raise RuntimeError(f"canonicalFromPH2 failed: {str(e)}")


def lib_butools_ph_check(alpha, A):
    """Validate PH representation."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.CheckPHRepresentationKt.checkPHRepresentation(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A), 1e-14
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_me_check(alpha, A):
    """Validate ME representation."""
    alpha = np.asarray(alpha)
    A = np.asarray(A)
    try:
        ph = jpype.JPackage('jline').lib.butools.ph
        result = ph.CheckMERepresentationKt.checkMERepresentation(
            jlineMatrixFromArray(alpha), jlineMatrixFromArray(A), 1e-14
        )
        return bool(result)
    except Exception as e:
        return False


# ========== MOMENT TRANSFORMATIONS ==========

def lib_butools_moment_factorial_from_raw(moments):
    """Convert raw moments to factorial moments."""
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        bt = jpype.JPackage('jline').lib.butools
        result = bt.FactorialMomsFromMomsKt.factorialMomsFromMoms(
            jlineMatrixFromArray(moments)
        )
        return jlineMatrixToArray(result).flatten()
    except Exception as e:
        raise RuntimeError(f"factorialMomsFromMoms failed: {str(e)}")


def lib_butools_moment_raw_from_factorial(moments):
    """Convert factorial moments to raw moments."""
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        bt = jpype.JPackage('jline').lib.butools
        result = bt.MomsFromFactorialMomsKt.MomsFromFactorialMoms(
            jlineMatrixFromArray(moments)
        )
        return jlineMatrixToArray(result).flatten()
    except Exception as e:
        raise RuntimeError(f"MomsFromFactorialMoms failed: {str(e)}")


def lib_butools_moment_hankel_from_raw(moments):
    """Convert raw moments to Hankel moments."""
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        bt = jpype.JPackage('jline').lib.butools
        result = bt.HankelMomsFromMomsKt.hankelMomsFromMoms(
            jlineMatrixFromArray(moments)
        )
        return jlineMatrixToArray(result).flatten()
    except Exception as e:
        raise RuntimeError(f"hankelMomsFromMoms failed: {str(e)}")


def lib_butools_moment_normalized_from_raw(moments):
    """Convert raw moments to normalized moments."""
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        bt = jpype.JPackage('jline').lib.butools
        result = bt.NormMomsFromMomsKt.NormMomsFromMoms(
            _to_double_array(moments)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"NormMomsFromMoms failed: {str(e)}")


def lib_butools_moment_reduced_from_raw(moments):
    """Convert raw moments to reduced moments."""
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        bt = jpype.JPackage('jline').lib.butools
        result = bt.ReducedMomsFromMomsKt.reducedMomsFromMoms(
            jlineMatrixFromArray(moments)
        )
        return jlineMatrixToArray(result).flatten()
    except Exception as e:
        raise RuntimeError(f"reducedMomsFromMoms failed: {str(e)}")


def lib_butools_moment_joint_factorial(moments):
    """Convert joint moments to joint factorial moments."""
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        bt = jpype.JPackage('jline').lib.butools
        result = bt.JFactorialMomsFromJMomsKt.jFactorialMomsFromJMoms(
            jlineMatrixFromArray(moments)
        )
        return jlineMatrixToArray(result).flatten()
    except Exception as e:
        raise RuntimeError(f"jFactorialMomsFromJMoms failed: {str(e)}")


# ========== MAP/MMAP FUNCTIONS ==========

def lib_butools_map_check(D0, D1):
    """Validate MAP representation."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        m = jpype.JPackage('jline').lib.butools.map
        result = m.CheckMAPRepresentationKt.checkMAPRepresentation(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), 1e-12
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_map_marginal_moments(D0, D1, K=5):
    """Compute marginal moments from MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        m = jpype.JPackage('jline').lib.butools.map
        result = m.MarginalMomentsFromMAPKt.marginalMomentsFromMAP(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), int(K)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"marginalMomentsFromMAP failed: {str(e)}")


def lib_butools_map_lag_correlations(D0, D1, L=1):
    """Compute lag-k autocorrelations from MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        m = jpype.JPackage('jline').lib.butools.map
        result = m.LagCorrelationsFromMAPKt.lagCorrelationsFromMAP(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), int(L)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"lagCorrelationsFromMAP failed: {str(e)}")


# ========== MARKOV CHAIN FUNCTIONS ==========

def lib_butools_markov_check_generator(Q):
    """Validate CTMC generator matrix."""
    Q = np.asarray(Q)
    try:
        mc = jpype.JPackage('jline').lib.butools.mc
        result = mc.CheckGeneratorKt.checkGenerator(jlineMatrixFromArray(Q), False, 1e-14)
        return bool(result)
    except Exception as e:
        return False


def lib_butools_markov_check_probability_matrix(P):
    """Validate DTMC transition matrix."""
    P = np.asarray(P)
    try:
        mc = jpype.JPackage('jline').lib.butools.mc
        result = mc.CheckProbMatrixKt.checkProbMatrix(jlineMatrixFromArray(P), False, 1e-14)
        return bool(result)
    except Exception as e:
        return False


def lib_butools_markov_ctmc_solve(Q):
    """Solve CTMC for stationary distribution."""
    Q = np.asarray(Q)
    try:
        mc = jpype.JPackage('jline').lib.butools.mc
        result = mc.CTMCSolveKt.ctmcSolve(jlineMatrixFromArray(Q), 1e-14)
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"ctmcSolve failed: {str(e)}")


def lib_butools_markov_dtmc_solve(P):
    """Solve DTMC for stationary distribution."""
    P = np.asarray(P)
    try:
        mc = jpype.JPackage('jline').lib.butools.mc
        result = mc.DTMCSolveKt.dtmcSolve(jlineMatrixFromArray(P), 1e-14)
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"dtmcSolve failed: {str(e)}")


def lib_butools_markov_crp_solve(Q):
    """Solve continuous renewal process."""
    Q = np.asarray(Q)
    try:
        mc = jpype.JPackage('jline').lib.butools.mc
        result = mc.CRPSolveKt.crpSolve(jlineMatrixFromArray(Q), 1e-14)
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"crpSolve failed: {str(e)}")


def lib_butools_markov_drp_solve(P):
    """Solve discrete renewal process."""
    P = np.asarray(P)
    try:
        mc = jpype.JPackage('jline').lib.butools.mc
        result = mc.CRPSolveKt.drpSolve(jlineMatrixFromArray(P))
        return jlineMatrixToArray(result)
    except Exception as e:
        raise RuntimeError(f"drpSolve failed: {str(e)}")


# ========== UTILITY FUNCTIONS ==========

def lib_butools_check_moments(moments):
    """Validate moment sequence feasibility."""
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        bt = jpype.JPackage('jline').lib.butools
        result = bt.CheckMomentsKt.checkMoments(jlineMatrixFromArray(moments), 1e-14)
        return bool(result)
    except Exception as e:
        return False


# ========== DISCRETE PHASE-TYPE (DPH) FUNCTIONS ==========

def lib_butools_dph2_from_moments(m1, m2, m3):
    """Fit 2-phase discrete PH from 3 moments.

    Returns (beta, B) tuple.
    """
    try:
        dph = jpype.JPackage('jline').lib.butools.dph
        result = dph.DPH2From3MomentsKt.dph2From3Moments(
            _to_double_array([float(m1), float(m2), float(m3)])
        )
        beta = jlineMatrixToArray(result.getBeta())
        B = jlineMatrixToArray(result.getB())
        return beta, B
    except Exception as e:
        raise RuntimeError(f"dph2From3Moments failed: {str(e)}")


def lib_butools_dph3_from_moments(moments):
    """Fit 3-phase discrete PH from 5 moments.

    Returns (beta, B) tuple.
    """
    moments = np.asarray(moments, dtype=float).flatten()
    try:
        dph = jpype.JPackage('jline').lib.butools.dph
        result = dph.DPH3From5MomentsKt.dph3From5Moments(
            _to_double_array(moments), 1e-14
        )
        beta = jlineMatrixToArray(result.getBeta())
        B = jlineMatrixToArray(result.getB())
        return beta, B
    except Exception as e:
        raise RuntimeError(f"dph3From5Moments failed: {str(e)}")


def lib_butools_dph_canonical_form(sigma, D):
    """Convert discrete PH to canonical form.

    Returns (beta, B) tuple.
    """
    sigma = np.asarray(sigma)
    D = np.asarray(D)
    try:
        dph = jpype.JPackage('jline').lib.butools.dph
        result = dph.CanonicalFromDPH2Kt.canonicalFromDPH2(
            jlineMatrixFromArray(sigma), jlineMatrixFromArray(D), 1e-14
        )
        beta = jlineMatrixToArray(result.getBeta())
        B = jlineMatrixToArray(result.getB())
        return beta, B
    except Exception as e:
        raise RuntimeError(f"canonicalFromDPH2 failed: {str(e)}")


def lib_butools_dph_pmf(sigma, D, x):
    """Compute probability mass function for discrete PH."""
    sigma = np.asarray(sigma)
    D = np.asarray(D)
    x = np.atleast_1d(np.asarray(x, dtype=int))
    try:
        dph = jpype.JPackage('jline').lib.butools.dph
        result = dph.PmfFromDPHKt.pmfFromDPH(
            jlineMatrixFromArray(sigma), jlineMatrixFromArray(D),
            _to_int_array(x)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"pmfFromDPH failed: {str(e)}")


def lib_butools_dph_moments(sigma, D, K=3):
    """Compute moments from discrete PH representation."""
    sigma = np.asarray(sigma)
    D = np.asarray(D)
    try:
        dph = jpype.JPackage('jline').lib.butools.dph
        result = dph.MomentsFromDPHKt.momentsFromDPH(
            jlineMatrixFromArray(sigma), jlineMatrixFromArray(D), int(K)
        )
        return np.array(list(result))
    except Exception as e:
        raise RuntimeError(f"momentsFromDPH failed: {str(e)}")


# ========== REPRESENTATION TRANSFORMATIONS ==========
# These functions are not implemented in the JAR backend.

def lib_butools_transform_ph_to_me(alpha, A):
    """Transform PH representation to matrix-exponential form."""
    raise NotImplementedError("phToMe is not implemented in the JAR backend")


def lib_butools_transform_me_to_ph(D0, D1):
    """Transform matrix-exponential representation to PH form."""
    raise NotImplementedError("meToPh is not implemented in the JAR backend")


def lib_butools_transform_ph_to_dph(alpha, A):
    """Transform PH to discrete PH representation."""
    raise NotImplementedError("phToDph is not implemented in the JAR backend")


def lib_butools_transform_dph_to_ph(sigma, D):
    """Transform discrete PH to continuous PH representation."""
    raise NotImplementedError("dphToPh is not implemented in the JAR backend")


# ========== DISCRETE MAP (DMAP) FUNCTIONS ==========

def lib_butools_dmap_check(D0, D1, prec=1e-14):
    """Validate discrete MAP representation."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.CheckDMAPRepresentationKt.checkDMAPRepresentation(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), float(prec)
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_drap_check(H0, H1, prec=1e-14):
    """Validate discrete rational arrival process representation."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.CheckDRAPRepresentationKt.checkDRAPRepresentation(
            jlineMatrixFromArray(H0), jlineMatrixFromArray(H1), float(prec)
        )
        return bool(result)
    except Exception as e:
        return False


def lib_butools_dmap_marginal_moments(D0, D1, K=0, prec=1e-14):
    """Compute marginal moments from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.MarginalMomentsFromDMAPKt.marginalMomentsFromDMAP(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), int(K), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"marginalMomentsFromDMAP failed: {str(e)}")


def lib_butools_drap_marginal_moments(H0, H1, K=0, prec=1e-14):
    """Compute marginal moments from discrete rational arrival process."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.MarginalMomentsFromDRAPKt.marginalMomentsFromDRAP(
            jlineMatrixFromArray(H0), jlineMatrixFromArray(H1), int(K), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"marginalMomentsFromDRAP failed: {str(e)}")


def lib_butools_dmap_lag_correlations(D0, D1, L=1, prec=1e-14):
    """Compute lag autocorrelations from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.LagCorrelationsFromDMAPKt.lagCorrelationsFromDMAP(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), int(L), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"lagCorrelationsFromDMAP failed: {str(e)}")


def lib_butools_drap_lag_correlations(H0, H1, L=1, prec=1e-14):
    """Compute lag autocorrelations from discrete rational arrival process."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.LagCorrelationsFromDRAPKt.lagCorrelationsFromDRAP(
            jlineMatrixFromArray(H0), jlineMatrixFromArray(H1), int(L), float(prec)
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"lagCorrelationsFromDRAP failed: {str(e)}")


def lib_butools_dmap_marginal_distribution(D0, D1, prec=1e-14):
    """Get marginal distribution from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.MarginalDistributionFromDMAPKt.marginalDistributionFromDMAP(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), float(prec)
        )
        alpha = jlineMatrixToArray(result.getAlpha())
        A = jlineMatrixToArray(result.getA())
        return {'alpha': alpha, 'A': A}
    except Exception as e:
        raise RuntimeError(f"marginalDistributionFromDMAP failed: {str(e)}")


def lib_butools_dmap_samples(D0, D1, K, prec=1e-14):
    """Generate random samples from discrete MAP."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.SamplesFromDMAPKt.samplesFromDMAP(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), int(K), None, float(prec), None
        )
        return np.array(result)
    except Exception as e:
        raise RuntimeError(f"samplesFromDMAP failed: {str(e)}")


def lib_butools_dmap_from_drap(H0, H1, prec=1e-14):
    """Convert discrete rational arrival process to discrete MAP."""
    H0 = np.asarray(H0)
    H1 = np.asarray(H1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.DMAPFromDRAPKt.dmapFromDRAP(
            jlineMatrixFromArray(H0), jlineMatrixFromArray(H1), float(prec)
        )
        D0 = jlineMatrixToArray(result.getFirst())
        D1 = jlineMatrixToArray(result.getSecond())
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"dmapFromDRAP failed: {str(e)}")


def lib_butools_dmap2_from_moments(moms, corr1):
    """Create discrete MAP(2) from 3 moments and lag-1 autocorrelation."""
    moms = np.asarray(moms).flatten()
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.DMAP2FromMomentsKt.dmap2FromMoments(
            moms.tolist(), float(corr1)
        )
        D0 = jlineMatrixToArray(result.getFirst())
        D1 = jlineMatrixToArray(result.getSecond())
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"dmap2FromMoments failed: {str(e)}")


def lib_butools_dmap_canonical2(D0, D1, prec=1e-14):
    """Convert order-2 discrete MAP to canonical form."""
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.CanonicalFromDMAP2Kt.canonicalFromDMAP2(
            jlineMatrixFromArray(D0), jlineMatrixFromArray(D1), float(prec)
        )
        G0 = jlineMatrixToArray(result.getFirst())
        G1 = jlineMatrixToArray(result.getSecond())
        return {'G0': G0, 'G1': G1}
    except Exception as e:
        raise RuntimeError(f"canonicalFromDMAP2 failed: {str(e)}")


def lib_butools_random_dmap(order, mean=10.0, zero_entries=0, max_trials=1000, prec=1e-7):
    """Generate a random discrete MAP."""
    try:
        dmap = jpype.JPackage('jline').lib.butools.dmap
        result = dmap.RandomDMAPKt.randomDMAP(
            int(order), float(mean), int(zero_entries), int(max_trials), float(prec), None
        )
        D0 = jlineMatrixToArray(result.getFirst())
        D1 = jlineMatrixToArray(result.getSecond())
        return {'D0': D0, 'D1': D1}
    except Exception as e:
        raise RuntimeError(f"randomDMAP failed: {str(e)}")
