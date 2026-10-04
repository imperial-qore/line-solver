"""
Markov-modulated fluid queue (MFQ) API wrappers.

Thin native-Python entry points that delegate to the BUTools-family fluid tools
in ``line_solver.lib.thirdparty.butools``. These mirror the MATLAB
``matlab/src/api/mam/mfq_*.m`` functions:

  * mfq_prio_queue     - continuous-time fluid priority queue (Horvath 2015)
  * mfq_ld_solve       - first/second-order level-dependent (multi-regime) solve
  * mfq_ld_distr       - stationary pdf/cdf of a level-dependent fluid queue
  * mfq_ld_mean        - stationary mean fluid level
  * mfq_multiregime    - multi-regime feedback fluid queue (Kankaya-Akar)
  * mfq_sojourn        - fluid-queue sojourn-time distribution as ME/PH
  * mfq_fluflu_sojourn - fluid/fluid-queue sojourn-time distribution as ME/PH
"""
from line_solver.lib.thirdparty.butools.queues.fluidprio import FluidPrioQueue
from line_solver.lib.thirdparty.butools.queues.fluidstd import FluidQueueSTD, FluFluSTD
from line_solver.lib.thirdparty.butools.mam.ldfluid import (
    SecondOrderLevelDependentFluidSolve,
    LevelDependentFluidStationaryDistr,
    LevelDependentFluidStationaryMean,
    multiregime,
)

__all__ = [
    "mfq_prio_queue", "mfq_ld_solve", "mfq_ld_distr", "mfq_ld_mean",
    "mfq_multiregime", "mfq_sojourn", "mfq_fluflu_sojourn",
]


def mfq_prio_queue(Q, R, d, *args, prec=1e-14, erlMaxOrder=200, classes=None):
    """Performance measures of a continuous-time fluid priority queue."""
    return FluidPrioQueue(Q, R, d, *args, prec=prec, erlMaxOrder=erlMaxOrder, classes=classes)


def mfq_ld_solve(Q, R, S, T, boundaryL=None, boundaryU=None, Qt=None, prec=1e-14):
    """Matrix-exponential solution of a multi-regime first/second-order fluid queue."""
    return SecondOrderLevelDependentFluidSolve(Q, R, S, T, boundaryL, boundaryU, Qt, prec)


def mfq_ld_distr(masses, iniF, KF, cloF, iniB, KB, cloB, T, what, points):
    """Stationary fluid-level distribution ('pdf','pdfd','cdf','cdfm')."""
    return LevelDependentFluidStationaryDistr(masses, iniF, KF, cloF, iniB, KB, cloB, T, what, points)


def mfq_ld_mean(masses, iniF, KF, cloF, iniB, KB, cloB, T):
    """Stationary mean fluid level E[X]."""
    return LevelDependentFluidStationaryMean(masses, iniF, KF, cloF, iniB, KB, cloB, T)


def mfq_multiregime(Q, R, Qt, Rt, T, pdfpoints, cdfpoints):
    """Multi-regime feedback fluid queue; returns (pdf, pdfd, cdf, cdfm)."""
    return multiregime(Q, R, Qt, Rt, T, pdfpoints, cdfpoints)


def mfq_sojourn(Q, Rin, Rout, Q0=None, transToPH=False):
    """Sojourn-time distribution of a fluid queue as an ME/PH representation."""
    return FluidQueueSTD(Q, Rin, Rout, Q0, transToPH)


def mfq_fluflu_sojourn(Qin, Rin, Qout, Rout, srv0stop, transToPH=False):
    """Sojourn-time distribution of a fluid/fluid queue as an ME/PH representation."""
    return FluFluSTD(Qin, Rin, Qout, Rout, srv0stop, transToPH)
