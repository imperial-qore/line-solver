"""
MAP distance library functions.

Provides distance measures between continuous and discrete-time
Markovian Arrival Processes (MAPs and D-MAPs).

Continuous-time distances based on:
    G. Horvath, "Measuring the distance between MAPs and some
    applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
    https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8

Discrete-time extensions by QORE Lab (https://qore.doc.ic.ac.uk/)
"""


# ========== CONTINUOUS MAP DISTANCES ==========

def lib_mapdist_exp_mul_int(A0, A1, B0, B1, L, alA=None, alB=None):
    """Joint density inner product of two MAPs via recursive Lyapunov."""
    from line_solver.api.mapdist.continuous import map_exp_mul_int
    return map_exp_mul_int(A0, A1, B0, B1, L, alA, alB)


def lib_mapdist_dist(A0, A1, B0, B1, L, alA=None, alB=None):
    """Squared L2 distance between lag-L joint densities of two MAPs."""
    from line_solver.api.mapdist.continuous import map_dist
    return map_dist(A0, A1, B0, B1, L, alA, alB)


def lib_mapdist_dist_lag1(A0, A1, B0, B1, alA=None, alB=None):
    """Lag-1 joint density L2 distance via Kronecker/Lyapunov."""
    from line_solver.api.mapdist.continuous import map_dist_lag1
    return map_dist_lag1(A0, A1, B0, B1, alA, alB)


def lib_mapdist_geo_mul_sum(A0, A1, B0, B1, alA=None, alB=None):
    """Geometric sum for autocorrelation distance computation."""
    from line_solver.api.mapdist.continuous import map_geo_mul_sum
    return map_geo_mul_sum(A0, A1, B0, B1, alA, alB)


def lib_mapdist_dist_acf(A0, A1, B0, B1, alA=None, alB=None):
    """Squared L2 distance between autocorrelation functions of two MAPs."""
    from line_solver.api.mapdist.continuous import map_dist_acf
    return map_dist_acf(A0, A1, B0, B1, alA, alB)


def lib_mapdist_optim_dist(A0, A1, alA, B0, alB, L):
    """Find B1 minimizing lag-L joint density distance given fixed B0."""
    from line_solver.api.mapdist.continuous import map_optim_dist
    return map_optim_dist(A0, A1, alA, B0, alB, L)


def lib_mapdist_optim_dist_acf(A0, A1, alA, B0, alB):
    """Find B1 minimizing autocorrelation distance given fixed B0."""
    from line_solver.api.mapdist.continuous import map_optim_dist_acf
    return map_optim_dist_acf(A0, A1, alA, B0, alB)


# ========== DISCRETE MAP (D-MAP) DISTANCES ==========

def lib_dmapdist_geo_mul_sum(D0A, D1A, D0B, D1B, L, alA=None, alB=None):
    """Joint PMF inner product of two D-MAPs via discrete Lyapunov."""
    from line_solver.api.mapdist.discrete import dmap_geo_mul_sum
    return dmap_geo_mul_sum(D0A, D1A, D0B, D1B, L, alA, alB)


def lib_dmapdist_dist(D0A, D1A, D0B, D1B, L, alA=None, alB=None):
    """Squared L2 distance between lag-L joint PMFs of two D-MAPs."""
    from line_solver.api.mapdist.discrete import dmap_dist
    return dmap_dist(D0A, D1A, D0B, D1B, L, alA, alB)


def lib_dmapdist_dist_lag1(D0A, D1A, D0B, D1B, alA=None, alB=None):
    """Lag-1 joint PMF L2 distance via Kronecker/discrete-Lyapunov."""
    from line_solver.api.mapdist.discrete import dmap_dist_lag1
    return dmap_dist_lag1(D0A, D1A, D0B, D1B, alA, alB)


def lib_dmapdist_geo_mul_sum_acf(D0A, D1A, D0B, D1B, alA=None, alB=None):
    """Geometric sum for discrete autocorrelation distance."""
    from line_solver.api.mapdist.discrete import dmap_geo_mul_sum_acf
    return dmap_geo_mul_sum_acf(D0A, D1A, D0B, D1B, alA, alB)


def lib_dmapdist_dist_acf(D0A, D1A, D0B, D1B, alA=None, alB=None):
    """Squared L2 distance between autocorrelation functions of two D-MAPs."""
    from line_solver.api.mapdist.discrete import dmap_dist_acf
    return dmap_dist_acf(D0A, D1A, D0B, D1B, alA, alB)


def lib_dmapdist_optim_dist(D0A, D1A, alA, D0B, alB, L):
    """Find D1B minimizing lag-L joint PMF distance given fixed D0B."""
    from line_solver.api.mapdist.discrete import dmap_optim_dist
    return dmap_optim_dist(D0A, D1A, alA, D0B, alB, L)


def lib_dmapdist_optim_dist_acf(D0A, D1A, alA, D0B, alB):
    """Find D1B minimizing autocorrelation distance given fixed D0B."""
    from line_solver.api.mapdist.discrete import dmap_optim_dist_acf
    return dmap_optim_dist_acf(D0A, D1A, alA, D0B, alB)
