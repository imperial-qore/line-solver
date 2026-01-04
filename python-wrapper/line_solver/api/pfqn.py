
"""
Product-form queueing network (PFQN) algorithms.

This module implements analytical algorithms for product-form queueing
networks, including Mean Value Analysis (MVA), normalizing constant
methods, and various approximation techniques.

Key algorithms:
- pfqn_mva: Standard Mean Value Analysis
- pfqn_ca: Convolution Algorithm  
- pfqn_nc: Normalizing Constant methods
- pfqn_bs: Balanced System analysis
- pfqn_aql: Approximate queueing lengths
- Various load-dependent and multi-class extensions

These functions provide the computational core for product-form
network analysis in LINE.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def pfqn_ca(N, L, Z):
    """
    Convolution Algorithm for product-form networks.

    Computes normalizing constants using the convolution algorithm.

    Args:
        N: Population vector.
        L: Demand matrix (service demands = mean visits × mean service times).
        Z: Think time vector.

    Returns:
        tuple: (G, lG) where G is normalizing constant and lG is log(G).
    """
    # Note: Kotlin function expects (L, N, Z) order
    ret = jpype.JPackage('jline').api.pfqn.nc.Pfqn_caKt.pfqn_ca(
        jlineMatrixFromArray(np.asarray(L)),
        jlineMatrixFromArray(np.asarray(N)),
        jlineMatrixFromArray(np.asarray(Z))
    )

    return float(ret.G), float(ret.lG)


def pfqn_panacea(N, L, Z):
    """
    PANACEA algorithm for product-form networks.

    A hybrid algorithm that combines convolution and MVA techniques
    for efficient computation of normalizing constants.

    Args:
        N: Population vector
        L: Demand matrix (service demands)
        Z: Think time vector

    Returns:
        tuple: (G, lG) where G is normalizing constant and lG is log(G)
    """
    ret = jpype.JPackage('jline').api.pfqn.nc.Pfqn_panaceaKt.pfqn_panacea(
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(Z)
    )

    return ret.G, ret.lG


def pfqn_bs(N, L, Z):
    """
    Balanced System analysis for product-form networks.

    Computes performance measures assuming balanced system conditions.

    Args:
        N: Population vector
        L: Demand matrix (service demands)
        Z: Think time vector

    Returns:
        tuple: (XN, CN, QN, UN, RN, TN, AN) where:
            - XN: Throughputs
            - CN: Response times
            - QN: Queue lengths
            - UN: Utilizations
            - RN: Residence times
            - TN: Node throughputs
            - AN: Arrival rates
    """
    ret = jpype.JPackage('jline').api.pfqn.mva.Pfqn_bsKt.pfqn_bs(
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(Z)
    )

    XN = jlineMatrixToArray(ret.XN)
    QN = jlineMatrixToArray(ret.QN)
    UN = jlineMatrixToArray(ret.UN)
    RN = jlineMatrixToArray(ret.RN)
    TN = jlineMatrixToArray(ret.TN)
    AN = jlineMatrixToArray(ret.AN)

    CN = np.sum(RN, axis=0) + np.diag(Z)

    return XN, CN, QN, UN, RN, TN, AN


def pfqn_mva(N, L, Z, mi=None):
    """
    Mean Value Analysis for product-form networks.

    Computes performance measures using the MVA algorithm.

    Args:
        N: Population vector.
        L: Demand matrix (service demands = mean visits × mean service times).
        Z: Think time vector.
        mi: Server multiplicity vector (optional, defaults to 1 for each station).

    Returns:
        Performance measures including throughputs, response times,
        queue lengths, and utilizations.
    """
    L_mat = jlineMatrixFromArray(np.asarray(L))
    N_mat = jlineMatrixFromArray(np.asarray(N))
    Z_mat = jlineMatrixFromArray(np.asarray(Z))

    # Create mi matrix (server multiplicities) - defaults to all 1s
    M = np.asarray(L).shape[0] if np.asarray(L).ndim > 1 else len(np.asarray(L))
    if mi is None:
        mi_mat = jlineMatrixFromArray(np.ones((1, M)))
    else:
        mi_mat = jlineMatrixFromArray(np.asarray(mi))

    ret = jpype.JPackage('jline').api.pfqn.mva.Pfqn_mvaKt.pfqn_mva(
        L_mat, N_mat, Z_mat, mi_mat
    )

    XN = jlineMatrixToArray(ret.X)  # Throughputs
    QN = jlineMatrixToArray(ret.Q)  # Queue lengths
    UN = jlineMatrixToArray(ret.U)  # Utilizations
    RN = jlineMatrixToArray(ret.R)  # Residence times

    # Compute derived metrics
    Z_arr = np.asarray(Z)
    if Z_arr.ndim == 1:
        CN = np.sum(RN, axis=0) + Z_arr
    else:
        CN = np.sum(RN, axis=0) + np.diag(Z_arr)

    # TN (node throughputs) and AN (arrival rates) - compute from XN and visits
    # For now, use XN as placeholder since visit ratios aren't passed
    TN = XN.copy()
    AN = XN.copy()

    return XN, CN, QN, UN, RN, TN, AN


def pfqn_aql(N, L, Z):
    """
    Approximate Queue Length algorithm for product-form networks.

    Provides approximations for queue lengths in product-form networks
    when exact algorithms are computationally expensive.

    Args:
        N: Population vector
        L: Demand matrix (service demands)
        Z: Think time vector

    Returns:
        tuple: (XN, CN, QN, UN, RN, TN, AN) - Performance measures
    """
    ret = jpype.JPackage('jline').api.pfqn.mva.Pfqn_aqlKt.pfqn_aql(
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(Z),
        jpype.JInt(1000)
    )

    XN = jlineMatrixToArray(ret.XN)
    QN = jlineMatrixToArray(ret.QN)
    UN = jlineMatrixToArray(ret.UN)
    RN = jlineMatrixToArray(ret.RN)
    TN = jlineMatrixToArray(ret.TN)
    AN = jlineMatrixToArray(ret.AN)

    CN = np.sum(RN, axis=0) + np.diag(Z)

    return XN, CN, QN, UN, RN, TN, AN


def pfqn_mvald(L, N, Z, mu, stabilize=True):
    """
    MVA with Load-Dependent service rates.

    Mean Value Analysis for networks with load-dependent service rates.

    Args:
        L: Demand matrix
        N: Population vector  
        Z: Think time vector
        mu: Load-dependent service rate functions
        stabilize: Whether to use numerical stabilization (default: True)

    Returns:
        tuple: (XN, QN, UN, CN, lGN, isNumStable, newpi) where:
            - Performance measures and numerical stability indicators
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_mvaldKt.pfqn_mvald(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jlineMatrixFromArray(mu),
        jpype.JBoolean(stabilize)
    )

    XN = jlineMatrixToArray(result.XN)
    QN = jlineMatrixToArray(result.QN)
    UN = jlineMatrixToArray(result.UN)
    CN = jlineMatrixToArray(result.CN)
    lGN = jlineMatrixToArray(result.lGN)

    return XN, QN, UN, CN, lGN, result.isNumStable, jlineMatrixToArray(result.newpi)


def pfqn_mvaldms(lambda_rates, D, N, Z, S):
    """
    MVA with Load-Dependent Multi-Server stations.

    Mean Value Analysis for networks with load-dependent multi-server stations.

    Args:
        lambda_rates: Arrival rates
        D: Service demands
        N: Population vector
        Z: Think times
        S: Server configurations

    Returns:
        Performance measures for load-dependent multi-server network
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_mvaldmsKt.pfqn_mvaldms(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(D),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jlineMatrixFromArray(S)
    )

    XN = jlineMatrixToArray(result.XN)
    QN = jlineMatrixToArray(result.QN)
    UN = jlineMatrixToArray(result.UN)
    CN = jlineMatrixToArray(result.CN)
    lGN = jlineMatrixToArray(result.lGN)

    return XN, QN, UN, CN, lGN


def pfqn_mvaldmx(lambda_rates, D, N, Z, mu=None, S=None):
    """
    MVA with Load-Dependent mixed service stations.

    Mean Value Analysis for networks with mixed load-dependent service stations
    that can have both load-dependent rates and multiple servers.

    Args:
        lambda_rates: Arrival rates
        D: Service demands
        N: Population vector
        Z: Think times
        mu: Load-dependent service rates (optional)
        S: Server configurations (optional)

    Returns:
        tuple: (XN, QN, UN, CN, lGN, newPc) - Performance measures with corrected parameters
    """
    mu_matrix = jlineMatrixFromArray(mu) if mu is not None else None
    S_matrix = jlineMatrixFromArray(S) if S is not None else None

    result = jpype.JPackage("jline").api.pfqn.Pfqn_mvaldmxKt.pfqn_mvaldmx(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(D),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        mu_matrix,
        S_matrix
    )

    XN = jlineMatrixToArray(result.XN)
    QN = jlineMatrixToArray(result.QN)
    UN = jlineMatrixToArray(result.UN)
    CN = jlineMatrixToArray(result.CN)
    lGN = jlineMatrixToArray(result.lGN)
    newPc = jlineMatrixToArray(result.newPc)

    return XN, QN, UN, CN, lGN, newPc


def pfqn_mvams(lambda_rates, L, N, Z, mi=None, S=None):
    """
    MVA for Multi-Server stations.

    Mean Value Analysis algorithm extended for multi-server stations.

    Args:
        lambda_rates: Arrival rates
        L: Service demands
        N: Population vector
        Z: Think times
        mi: Service rates (optional)
        S: Number of servers per station (optional)

    Returns:
        Performance measures for multi-server network
    """
    mi_matrix = jlineMatrixFromArray(mi) if mi is not None else None
    S_matrix = jlineMatrixFromArray(S) if S is not None else None

    result = jpype.JPackage("jline").api.pfqn.Pfqn_mvamsKt.pfqn_mvams(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        mi_matrix,
        S_matrix
    )

    XN = jlineMatrixToArray(result.XN)
    QN = jlineMatrixToArray(result.QN)
    UN = jlineMatrixToArray(result.UN)
    CN = jlineMatrixToArray(result.CN)
    lGN = jlineMatrixToArray(result.lGN)

    return XN, QN, UN, CN, lGN


def pfqn_mvamx(lambda_rates, D, N, Z, mi=None):
    """
    MVA for mixed open/closed networks.

    Mean Value Analysis for mixed networks containing both open and closed classes.

    Args:
        lambda_rates: Arrival rates for open classes
        D: Service demands
        N: Population vector for closed classes
        Z: Think times
        mi: Service rates (optional)

    Returns:
        tuple: (XN, QN, UN, CN, lGN) - Performance measures
    """
    mi_matrix = jlineMatrixFromArray(mi) if mi is not None else None

    result = jpype.JPackage("jline").api.pfqn.Pfqn_mvamxKt.pfqn_mvamx(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(D),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        mi_matrix
    )

    XN = jlineMatrixToArray(result.XN)
    QN = jlineMatrixToArray(result.QN)
    UN = jlineMatrixToArray(result.UN)
    CN = jlineMatrixToArray(result.CN)
    lGN = jlineMatrixToArray(result.lGN)

    return XN, QN, UN, CN, lGN


def pfqn_nc(lambda_rates, L, N, Z, options):
    """
    Normalizing Constant algorithm for product-form networks.

    Computes performance measures using normalizing constant methods.

    Args:
        lambda_rates: Arrival rates
        L: Service demands
        N: Population vector
        Z: Think times
        options: Algorithm options

    Returns:
        tuple: (lG, X, Q, method) where:
            - lG: Log normalizing constant
            - X: Throughputs
            - Q: Queue lengths
            - method: Method used
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_ncKt.pfqn_nc(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        options
    )

    lG = result.lG
    X = jlineMatrixToArray(result.X)
    Q = jlineMatrixToArray(result.Q)
    method = result.method

    return lG, X, Q, method


def pfqn_gld(L, N, mu=None, options=None):
    """
    Generalized Load-Dependent algorithm.

    Computes normalizing constants for load-dependent service stations.

    Args:
        L: Service demands
        N: Population vector
        mu: Load-dependent service rates (optional)
        options: Algorithm options (optional)

    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm
    """
    mu_matrix = jlineMatrixFromArray(mu) if mu is not None else None

    result = jpype.JPackage("jline").api.pfqn.Pfqn_gldKt.pfqn_gld(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        mu_matrix,
        options
    )

    return result.G, result.lG


def pfqn_gldsingle(L, N, mu, options=None):
    """
    GLD algorithm for single load-dependent station.

    Specialized version for networks with a single load-dependent station.

    Args:
        L: Service demands
        N: Population vector
        mu: Load-dependent service rates
        options: Algorithm options (optional)

    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_gldsingleKt.pfqn_gldsingle(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(mu),
        options
    )

    return result.G, result.lG


def pfqn_comomrm(L, N, Z, m=None, atol=1e-8):
    """
    Co-moment matching algorithm for product-form networks.
    
    Uses moment matching techniques to approximate performance metrics
    in product-form queueing networks with high accuracy.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        m: Number of moments to match (optional)
        atol: Absolute tolerance for convergence
        
    Returns:
        Performance metrics computed via co-moment matching
    """
    m_val = jpype.JInt(m) if m is not None else None

    result = jpype.JPackage("jline").api.pfqn.Pfqn_comomrmKt.pfqn_comomrm(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        m_val,
        jpype.JDouble(atol)
    )

    return result.lG, jlineMatrixToArray(result.lGbasis)



def pfqn_linearizer(L, N, Z, schedule_types, tol=1e-8, maxiter=1000):
    """
    Linearizer algorithm for non-product-form networks.

    Uses iterative approximation to analyze networks with non-product-form
    scheduling disciplines (e.g., priority, LCFS-PR).

    Args:
        L: Service demands
        N: Population vector
        Z: Think times
        schedule_types: Scheduling disciplines for each station
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)

    Returns:
        tuple: (Q, U, R, C, X, totiter) where:
            - Q: Queue lengths
            - U: Utilizations
            - R: Response times
            - C: Class response times
            - X: Throughputs
            - totiter: Total iterations performed
    """
    java_schedule_types = jpype.JArray(jpype.JPackage("jline").lang.constant.SchedStrategy)(len(schedule_types))
    for i, sched in enumerate(schedule_types):
        java_schedule_types[i] = sched

    result = jpype.JPackage("jline").api.pfqn.Pfqn_linearizerKt.pfqn_linearizer(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        java_schedule_types,
        jpype.JDouble(tol),
        jpype.JInt(maxiter)
    )

    Q = jlineMatrixToArray(result.Q)
    U = jlineMatrixToArray(result.U)
    R = jlineMatrixToArray(result.R)
    C = jlineMatrixToArray(result.C)
    X = jlineMatrixToArray(result.X)

    return Q, U, R, C, X, result.totiter


def pfqn_linearizerms(L, N, Z, nservers, schedule_types, tol=1e-8, maxiter=1000):
    """
    Linearizer algorithm for multi-server non-product-form networks.

    Extended linearizer for networks with multi-server stations and
    non-product-form scheduling disciplines.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times
        nservers: Number of servers per station
        schedule_types: Scheduling disciplines
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)

    Returns:
        tuple: (Q, U, W, C, X, totiter) - Performance measures and iterations
    """
    java_schedule_types = jpype.java.util.ArrayList()
    for sched in schedule_types:
        java_schedule_types.add(sched)

    result = jpype.JPackage("jline").api.pfqn.Pfqn_linearizermxKt.pfqn_linearizerms(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jlineMatrixFromArray(nservers),
        java_schedule_types,
        jpype.JDouble(tol),
        jpype.JInt(maxiter)
    )

    Q = jlineMatrixToArray(result.Q)
    U = jlineMatrixToArray(result.U)
    W = jlineMatrixToArray(result.W)
    C = jlineMatrixToArray(result.C)
    X = jlineMatrixToArray(result.X)

    return Q, U, W, C, X, result.totiter


def pfqn_linearizerpp(L, N, Z, level=2, tol=1e-8, maxiter=1000, flag=0):
    """
    Priority Preemptive linearizer algorithm.

    Specialized linearizer for networks with priority preemptive scheduling.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times
        level: Priority level (default: 2)
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)
        flag: Algorithm flag (default: 0)

    Returns:
        Performance measures for priority preemptive network
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_linearizerppKt.pfqn_linearizerpp(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jpype.JInt(level),
        jpype.JDouble(tol),
        jpype.JInt(maxiter),
        jpype.JInt(flag)
    )

    Q = jlineMatrixToArray(result.Q)
    U = jlineMatrixToArray(result.U)
    R = jlineMatrixToArray(result.R)
    C = jlineMatrixToArray(result.C)
    X = jlineMatrixToArray(result.X)

    return Q, U, R, C, X, result.totiter


def pfqn_linearizermx(lambda_rates, L, N, Z, nservers, schedule_types, tol=1e-8, maxiter=1000, method="lin"):
    """
    Mixed linearizer algorithm for open/closed networks with multi-server stations.

    Extended linearizer for mixed networks containing both open and closed
    classes with multi-server stations and non-product-form scheduling.

    Args:
        lambda_rates: Arrival rates for open classes
        L: Service demands
        N: Population vector for closed classes
        Z: Think times
        nservers: Number of servers per station
        schedule_types: Scheduling disciplines
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)
        method: Solution method (default: "lin")

    Returns:
        tuple: (Q, U, R, C, X, totiter) - Performance measures and iterations
    """
    java_schedule_types = jpype.JArray(jpype.JPackage("jline").lang.constant.SchedStrategy)(len(schedule_types))
    for i, sched in enumerate(schedule_types):
        java_schedule_types[i] = sched

    result = jpype.JPackage("jline").api.pfqn.Pfqn_linearizermxKt.pfqn_linearizermx(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jlineMatrixFromArray(nservers),
        java_schedule_types,
        jpype.JDouble(tol),
        jpype.JInt(maxiter),
        jpype.JString(method)
    )

    Q = jlineMatrixToArray(result.Q)
    U = jlineMatrixToArray(result.U)
    R = jlineMatrixToArray(result.R)
    C = jlineMatrixToArray(result.C)
    X = jlineMatrixToArray(result.X)

    return Q, U, R, C, X, result.totiter


def pfqn_kt(L, N, Z):
    """
    Korkmazoglu-Tucci (KT) algorithm for normalizing constants.

    Efficient algorithm for computing normalizing constants in
    product-form networks.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times

    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_ktKt.pfqn_kt(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z)
    )

    return result.G, result.lG


def pfqn_recal(L, N, Z=None, m0=None):
    """
    Recurrence algorithm for normalizing constants.

    Uses recurrence relations to compute normalizing constants.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times (optional)
        m0: Initial multiplicities (optional)

    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm
    """
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None
    m0_matrix = jlineMatrixFromArray(m0) if m0 is not None else None

    result = jpype.JPackage("jline").api.pfqn.Pfqn_recalKt.pfqn_recal(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        Z_matrix,
        m0_matrix
    )

    return result.G, result.lG


def pfqn_cub(L, N, Z, order=3, atol=1e-8):
    """
    Cubic spline approximation for normalizing constants.

    Uses cubic spline interpolation for approximate computation
    of normalizing constants.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times
        order: Spline order (default: 3)
        atol: Absolute tolerance (default: 1e-8)

    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_cubKt.pfqn_cub(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jpype.JInt(order),
        jpype.JDouble(atol)
    )

    return result.G, result.lG


def pfqn_mmint2(L, N, Z):
    """
    Two-moment interpolation algorithm for normalizing constants.

    Uses moment interpolation to approximate normalizing constants.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times

    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_mmint2Kt.pfqn_mmint2(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z)
    )

    return result.G, result.lG


def pfqn_ls(L, N, Z=None, I=10000, seed=12345):
    """
    Large-Scale approximation algorithm for normalizing constants.

    Monte Carlo based approximation for large-scale product-form networks.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times (optional)
        I: Number of Monte Carlo samples (default: 10000)
        seed: Random seed (default: 12345)

    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm
    """
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None

    result = jpype.JPackage("jline").api.pfqn.Pfqn_lsKt.pfqn_ls(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        Z_matrix,
        jpype.JLong(I),
        jpype.JLong(seed)
    )

    return result.G, result.lG


def pfqn_rd(L, N, Z, mu, options=None):
    """
    Recursive Doubling algorithm for load-dependent networks.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times
        mu: Load-dependent service rates
        options: Algorithm options (optional)

    Returns:
        tuple: (lG, Cgamma) - Log normalizing constant and performance measures
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_rdKt.pfqn_rd(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jlineMatrixFromArray(mu),
        options
    )

    return result.lG, jlineMatrixToArray(result.Cgamma)


def pfqn_fnc(alpha, c):
    """
    Flow-equivalent Node Centralization algorithm.

    Args:
        alpha: Alpha parameters
        c: Capacity parameters

    Returns:
        tuple: (mu, c) - Service rates and updated capacities
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_fncKt.pfqn_fnc(
        jlineMatrixFromArray(alpha),
        jlineMatrixFromArray(c)
    )

    return jlineMatrixToArray(result.mu), jlineMatrixToArray(result.c)


def pfqn_propfair(L, N, Z):
    """
    Proportional fair algorithm for product-form networks.

    Args:
        L: Service demands
        N: Population vector
        Z: Think times

    Returns:
        tuple: (G, lG, X, Q, method) - Performance measures and method used
    """
    result = jpype.JPackage("jline").api.pfqn.Pfqn_propfairKt.pfqn_propfair(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z)
    )

    return result.G, result.lG, jlineMatrixToArray(result.X), jlineMatrixToArray(result.Q), result.method


def pfqn_xia(L, N, s, options=None):
    """
    Xia's algorithm for multi-server networks.

    Args:
        L: Service demands
        N: Population count
        s: Server configurations
        options: Algorithm options (optional)

    Returns:
        Algorithm results
    """
    return jpype.JPackage("jline").api.pfqn.Pfqn_xiaKt.pfqn_xia(
        jlineMatrixFromArray(L),
        jpype.JInt(N),
        jlineMatrixFromArray(s),
        options
    )


def pfqn_xzabalow(L, N, Z):
    """
    Lower bound algorithm by Zahorjan and Abadir.

    Args:
        L: Service demands
        N: Population count
        Z: Think time

    Returns:
        Lower bound estimates
    """
    return jpype.JPackage("jline").api.pfqn.Pfqn_xzabalowKt.pfqn_xzabalow(
        jlineMatrixFromArray(L),
        jpype.JDouble(N),
        jpype.JDouble(Z)
    )


def pfqn_xzabaup(L, N, Z):
    """
    Upper bound algorithm by Zahorjan and Abadir.

    Args:
        L: Service demands
        N: Population count
        Z: Think time

    Returns:
        Upper bound estimates
    """
    return jpype.JPackage("jline").api.pfqn.Pfqn_xzabaupKt.pfqn_xzabaup(
        jlineMatrixFromArray(L),
        jpype.JDouble(N),
        jpype.JDouble(Z)
    )


def pfqn_xzgsblow(L, N, Z):
    """
    Lower bound algorithm by Zahorjan, Gesbert, and Sevcik.

    Args:
        L: Service demands
        N: Population count
        Z: Think time

    Returns:
        Lower bound estimates
    """
    return jpype.JPackage("jline").api.pfqn.Pfqn_xzgsblowKt.pfqn_xzgsblow(
        jlineMatrixFromArray(L),
        jpype.JDouble(N),
        jpype.JDouble(Z)
    )


def pfqn_xzgsbup(L, N, Z):
    """
    Upper bound algorithm by Zahorjan, Gesbert, and Sevcik.

    Args:
        L: Service demands
        N: Population count
        Z: Think time

    Returns:
        Upper bound estimates
    """
    return jpype.JPackage("jline").api.pfqn.Pfqn_xzgsbupKt.pfqn_xzgsbup(
        jlineMatrixFromArray(L),
        jpype.JDouble(N),
        jpype.JDouble(Z)
    )




def pfqn_conwayms(lambda_rates, L, N, Z, nservers, options=None):
    """
    Conway Multi-Server algorithm for product-form networks.
    
    Extends the standard algorithms to handle multi-server stations
    using Conway's method for efficient computation.
    
    Args:
        lambda_rates: Arrival rate matrix
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        nservers: Number of servers per station
        options: Optional solver options
        
    Returns:
        Performance metrics for multi-server network
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_conwaymsKt.pfqn_conwayms(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jlineMatrixFromArray(nservers),
        options
    )

    XN = jlineMatrixToArray(result.XN)
    QN = jlineMatrixToArray(result.QN)
    UN = jlineMatrixToArray(result.UN)
    RN = jlineMatrixToArray(result.RN)

    return XN, QN, UN, RN


def pfqn_egflinearizer(L, N, Z, tol=1e-8, maxiter=1000):
    """
    Enhanced Gauss-Friedman linearizer for product-form networks.
    
    Advanced linearization algorithm that provides improved convergence
    and accuracy compared to standard linearizer methods.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        tol: Convergence tolerance
        maxiter: Maximum number of iterations
        
    Returns:
        Performance metrics computed via enhanced linearizer
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_egflinearizerKt.pfqn_egflinearizer(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jpype.JDouble(tol),
        jpype.JInt(maxiter)
    )

    Q = jlineMatrixToArray(result.Q)
    U = jlineMatrixToArray(result.U)
    R = jlineMatrixToArray(result.R)
    X = jlineMatrixToArray(result.X)

    return Q, U, R, X, result.totiter


def pfqn_gflinearizer(L, N, Z, tol=1e-8, maxiter=1000):
    """
    Gauss-Friedman linearizer for product-form networks.
    
    Uses the Gauss-Friedman linearization method to approximate
    performance metrics in product-form queueing networks.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        tol: Convergence tolerance
        maxiter: Maximum number of iterations
        
    Returns:
        Performance metrics computed via Gauss-Friedman linearizer
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_gflinearizerKt.pfqn_gflinearizer(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jpype.JDouble(tol),
        jpype.JInt(maxiter)
    )

    Q = jlineMatrixToArray(result.Q)
    U = jlineMatrixToArray(result.U)
    R = jlineMatrixToArray(result.R)
    X = jlineMatrixToArray(result.X)

    return Q, U, R, X, result.totiter


def pfqn_gld_complex(L, N, mu=None, options=None):
    """
    Complex variant of Generalized Load-Dependent algorithm.
    
    Enhanced GLD method for complex load-dependent service configurations
    with support for advanced service rate dependencies.
    
    Args:
        L: Service demand matrix
        N: Population vector
        mu: Load-dependent service rates (optional)
        options: Optional solver options
        
    Returns:
        Performance metrics for complex load-dependent systems
    """
    mu_matrix = jlineMatrixFromArray(mu) if mu is not None else None

    result = jpype.JPackage('jline').api.pfqn.Pfqn_gld_complexKt.pfqn_gld_complex(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        mu_matrix,
        options
    )

    return result.G, result.lG


def pfqn_gldsingle_complex(L, N, mu, options=None):
    """
    Complex GLD algorithm for single load-dependent station.
    
    Specialized version for networks with one complex load-dependent station
    with sophisticated service rate modeling.
    
    Args:
        L: Service demand matrix
        N: Population vector
        mu: Complex load-dependent service rates
        options: Optional solver options
        
    Returns:
        Performance metrics for single complex load-dependent station
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_gldsingle_complexKt.pfqn_gldsingle_complex(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(mu),
        options
    )

    return result.G, result.lG


def pfqn_le_hessian(L, N, Z, mu=None, options=None):
    """
    Load Evaluation algorithm with Hessian computation.
    
    Computes performance metrics using Load Evaluation method
    with Hessian matrix computation for sensitivity analysis.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        mu: Service rates (optional).
        options: Algorithm options (optional).
        
    Returns:
        numpy.ndarray: Hessian matrix for sensitivity analysis.
    """
    mu_matrix = jlineMatrixFromArray(mu) if mu is not None else None

    return jlineMatrixToArray(
        jpype.JPackage('jline').api.pfqn.Pfqn_le_hessianKt.pfqn_le_hessian(
            jlineMatrixFromArray(L),
            jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z),
            mu_matrix,
            options
        )
    )


def pfqn_le_hessianZ(L, N, Z, mu=None, options=None):
    """
    Load Evaluation algorithm with Z-parameter Hessian computation.
    
    Computes performance metrics using Load Evaluation method
    with Hessian matrix computation for think time parameters.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        mu: Service rates (optional).
        options: Algorithm options (optional).
        
    Returns:
        numpy.ndarray: Hessian matrix with respect to Z parameters.
    """
    mu_matrix = jlineMatrixFromArray(mu) if mu is not None else None

    return jlineMatrixToArray(
        jpype.JPackage('jline').api.pfqn.Pfqn_le_hessianZKt.pfqn_le_hessianZ(
            jlineMatrixFromArray(L),
            jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z),
            mu_matrix,
            options
        )
    )


def pfqn_lldfun(L, N, Z, mu, options=None):
    """
    Log-likelihood derivative function for parameter estimation.
    
    Computes the derivative of the log-likelihood function for
    parameter estimation in product-form queueing networks.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        mu: Service rate parameters.
        options: Algorithm options (optional).
        
    Returns:
        numpy.ndarray: Log-likelihood derivative values.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.pfqn.Pfqn_lldfunKt.pfqn_lldfun(
            jlineMatrixFromArray(L),
            jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z),
            jlineMatrixFromArray(mu),
            options
        )
    )


def pfqn_mci(L, N, Z, nsample=10000, seed=12345):
    """
    Monte Carlo Integration algorithm for product-form networks.
    
    Uses Monte Carlo sampling to estimate performance metrics,
    particularly useful for large or complex networks.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        nsample: Number of Monte Carlo samples
        seed: Random seed for reproducibility
        
    Returns:
        Estimated performance metrics with statistical confidence
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mciKt.pfqn_mci(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jpype.JInt(nsample),
        jpype.JInt(seed)
    )

    return result.G, result.lG, result.variance


def pfqn_mmint2_gausslegendre(L, N, Z, order=5):
    """
    Two-moment interpolation with Gauss-Legendre quadrature.
    
    Computes normalizing constants using moment interpolation
    with Gauss-Legendre quadrature for numerical integration.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        order: Quadrature order (default: 5).
        
    Returns:
        tuple: (G, lG) - Normalizing constant and its logarithm.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mmint2_gausslegendreKt.pfqn_mmint2_gausslegendre(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jpype.JInt(order)
    )

    return result.G, result.lG


def pfqn_mmsample2(L, N, Z, nsample=10000, seed=12345):
    """
    Two-stage moment sampling algorithm.
    
    Computes normalizing constants using a two-stage Monte Carlo
    sampling approach with variance reduction techniques.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        nsample: Number of samples (default: 10000).
        seed: Random seed (default: 12345).
        
    Returns:
        tuple: (G, lG, samples) - Constants and sample data.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mmsample2Kt.pfqn_mmsample2(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jpype.JInt(nsample),
        jpype.JInt(seed)
    )

    return result.G, result.lG, jlineMatrixToArray(result.samples)


def pfqn_mushift(L, N, Z, mu, shift_factor=1.0):
    """
    Service rate shifting algorithm for load-dependent networks.
    
    Applies a shifting transformation to service rates for improved
    numerical stability in load-dependent network analysis.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        mu: Load-dependent service rates.
        shift_factor: Shifting parameter (default: 1.0).
        
    Returns:
        tuple: (mu_shifted, G, lG) - Shifted rates and normalizing constants.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mushiftKt.pfqn_mushift(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        jlineMatrixFromArray(mu),
        jpype.JDouble(shift_factor)
    )

    return jlineMatrixToArray(result.mu_shifted), result.G, result.lG


def pfqn_le(L, N, Z, options=None):
    """
    Load Evaluation algorithm for product-form networks.
    
    Computes performance metrics using the Load Evaluation method,
    an efficient approximation technique for large networks.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        options: Optional solver options
        
    Returns:
        Performance metrics (throughput, utilization, response time)
    """
    java_options = None
    if options is not None:
        pass

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_leKt.pfqn_le(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_leKt.pfqn_le(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z)
        )

    return {
        'X': jlineMatrixToArray(result.getX()) if hasattr(result, 'getX') else None,
        'R': jlineMatrixToArray(result.getR()) if hasattr(result, 'getR') else None,
        'Q': jlineMatrixToArray(result.getQ()) if hasattr(result, 'getQ') else None
    }


def pfqn_le_fpi(L, N, Z, max_iter=1000, tol=1e-6):
    """
    Load Evaluation with Fixed Point Iteration.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        max_iter: Maximum iterations (default: 1000)
        tol: Convergence tolerance (default: 1e-6)

    Returns:
        dict: Performance metrics with convergence information
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_le_fpiKt.pfqn_le_fpi(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jpype.JInt(max_iter), jpype.JDouble(tol)
    )

    return {
        'X': jlineMatrixToArray(result.getX()),
        'R': jlineMatrixToArray(result.getR()),
        'Q': jlineMatrixToArray(result.getQ()),
        'converged': bool(result.getConverged()) if hasattr(result, 'getConverged') else None,
        'iterations': int(result.getIterations()) if hasattr(result, 'getIterations') else None
    }


def pfqn_le_fpiZ(L, N, Z, max_iter=1000, tol=1e-6):
    """
    Load Evaluation FPI with Z-parameter computation.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        max_iter: Maximum iterations (default: 1000)
        tol: Convergence tolerance (default: 1e-6)

    Returns:
        dict: Performance metrics with Z coefficients
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_le_fpiZKt.pfqn_le_fpiZ(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jpype.JInt(max_iter), jpype.JDouble(tol)
    )

    return {
        'X': jlineMatrixToArray(result.getX()),
        'R': jlineMatrixToArray(result.getR()),
        'Q': jlineMatrixToArray(result.getQ()),
        'Z_coeffs': jlineMatrixToArray(result.getZCoeffs()) if hasattr(result, 'getZCoeffs') else None
    }


def pfqn_le_hessian(L, N, Z, step_size=1e-6):
    result = jpype.JPackage('jline').api.pfqn.Pfqn_le_hessianKt.pfqn_le_hessian(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jpype.JDouble(step_size)
    )

    return {
        'X': jlineMatrixToArray(result.getX()),
        'R': jlineMatrixToArray(result.getR()),
        'Q': jlineMatrixToArray(result.getQ()),
        'hessian': jlineMatrixToArray(result.getHessian()) if hasattr(result, 'getHessian') else None
    }


def pfqn_le_hessianZ(L, N, Z, step_size=1e-6):
    result = jpype.JPackage('jline').api.pfqn.Pfqn_le_hessianZKt.pfqn_le_hessianZ(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jpype.JDouble(step_size)
    )

    return {
        'X': jlineMatrixToArray(result.getX()),
        'R': jlineMatrixToArray(result.getR()),
        'Q': jlineMatrixToArray(result.getQ()),
        'hessian_Z': jlineMatrixToArray(result.getHessianZ()) if hasattr(result, 'getHessianZ') else None
    }


def pfqn_lldfun(L, N, Z, lambda_vec):
    result = jpype.JPackage('jline').api.pfqn.Pfqn_lldfunKt.pfqn_lldfun(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jlineMatrixFromArray(lambda_vec)
    )

    return {
        'lld': jlineMatrixToArray(result.getLLD()),
        'gradient': jlineMatrixToArray(result.getGradient()) if hasattr(result, 'getGradient') else None,
        'likelihood': float(result.getLikelihood()) if hasattr(result, 'getLikelihood') else None
    }


def pfqn_mci(L, N, Z, num_samples=10000, confidence=0.95):
    """
    Monte Carlo Integration with confidence intervals.
    
    Extended MCI method that provides confidence interval estimates
    for the computed performance metrics.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        num_samples: Number of Monte Carlo samples
        confidence: Confidence level for intervals (e.g., 0.95)
        
    Returns:
        Performance metrics with confidence intervals
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mciKt.pfqn_mci(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jpype.JInt(num_samples), jpype.JDouble(confidence)
    )

    return {
        'X_mean': jlineMatrixToArray(result.getXMean()),
        'X_ci': jlineMatrixToArray(result.getXCI()) if hasattr(result, 'getXCI') else None,
        'R_mean': jlineMatrixToArray(result.getRMean()),
        'R_ci': jlineMatrixToArray(result.getRCI()) if hasattr(result, 'getRCI') else None,
        'samples_used': int(result.getSamplesUsed()) if hasattr(result, 'getSamplesUsed') else num_samples
    }


def pfqn_mmint2_gausslegendre(L, N, Z, order=10):
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mmint2_gausslegendreKt.pfqn_mmint2_gausslegendre(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jpype.JInt(order)
    )

    return {
        'integral': jlineMatrixToArray(result.getIntegral()),
        'nodes': jlineMatrixToArray(result.getNodes()) if hasattr(result, 'getNodes') else None,
        'weights': jlineMatrixToArray(result.getWeights()) if hasattr(result, 'getWeights') else None
    }


def pfqn_mmsample2(L, N, Z, num_samples=1000):
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mmsample2Kt.pfqn_mmsample2(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jpype.JInt(num_samples)
    )

    return {
        'stage1_X': jlineMatrixToArray(result.getStage1X()),
        'stage2_X': jlineMatrixToArray(result.getStage2X()),
        'combined_X': jlineMatrixToArray(result.getCombinedX()),
        'variance_reduction': float(result.getVarianceReduction()) if hasattr(result, 'getVarianceReduction') else None
    }


def pfqn_mu_ms(L, N, Z, server_counts):
    """
    Service rate computation for multi-server stations.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        server_counts: Number of servers per station

    Returns:
        dict: Performance metrics with multi-server utilization
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mu_msKt.pfqn_mu_ms(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jlineMatrixFromArray(server_counts)
    )

    return {
        'X': jlineMatrixToArray(result.getX()),
        'R': jlineMatrixToArray(result.getR()),
        'Q': jlineMatrixToArray(result.getQ()),
        'U': jlineMatrixToArray(result.getU()) if hasattr(result, 'getU') else None
    }


def pfqn_nrl(L, N, Z, options=None):
    """
    Newton-Raphson Linearization for product-form networks.
    
    Uses Newton-Raphson optimization to solve the linearized system
    of equations for performance metrics computation.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        options: Optional solver options
        
    Returns:
        Performance metrics computed via Newton-Raphson linearization
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_nrlKt.pfqn_nrl(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_nrlKt.pfqn_nrl(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z)
        )

    return {
        'G': result.G if hasattr(result, 'G') else None,
        'logG': result.logG if hasattr(result, 'logG') else None,
        'converged': result.converged if hasattr(result, 'converged') else None
    }


def pfqn_nrp(L, N, Z, options=None):
    """
    Newton-Raphson method for performance measures.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        options: Algorithm options (optional)

    Returns:
        dict: Performance measures with convergence status
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_nrpKt.pfqn_nrp(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_nrpKt.pfqn_nrp(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z)
        )

    return {
        'G': result.G if hasattr(result, 'G') else None,
        'prob': jlineMatrixToArray(result.prob) if hasattr(result, 'prob') else None,
        'converged': result.converged if hasattr(result, 'converged') else None
    }


def pfqn_stdf(L, N, Z, S, fcfs_nodes, rates, tset):
    """
    Service Time Distribution Function computation.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        S: Server configuration
        fcfs_nodes: FCFS node indicators
        rates: Service rates
        tset: Time set for evaluation

    Returns:
        dict: Distribution function values (CDF, PDF, mean)
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_stdfKt.pfqn_stdf(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jlineMatrixFromArray(S),
        jlineMatrixFromArray(fcfs_nodes), jlineMatrixFromArray(rates),
        jlineMatrixFromArray(tset)
    )

    return {
        'cdf': jlineMatrixToArray(result.cdf) if hasattr(result, 'cdf') else None,
        'pdf': jlineMatrixToArray(result.pdf) if hasattr(result, 'pdf') else None,
        'mean': jlineMatrixToArray(result.mean) if hasattr(result, 'mean') else None
    }


def pfqn_stdf_heur(L, N, Z, S, fcfs_nodes, rates, tset, options=None):
    """
    Service Time Distribution Function with heuristic approximation.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        S: Server configuration
        fcfs_nodes: FCFS node indicators
        rates: Service rates
        tset: Time set for evaluation
        options: Algorithm options (optional)

    Returns:
        dict: Approximate distribution function values
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_stdf_heurKt.pfqn_stdf_heur(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S),
            jlineMatrixFromArray(fcfs_nodes), jlineMatrixFromArray(rates),
            jlineMatrixFromArray(tset), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_stdf_heurKt.pfqn_stdf_heur(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S),
            jlineMatrixFromArray(fcfs_nodes), jlineMatrixFromArray(rates),
            jlineMatrixFromArray(tset)
        )

    return {
        'cdf': jlineMatrixToArray(result.cdf) if hasattr(result, 'cdf') else None,
        'pdf': jlineMatrixToArray(result.pdf) if hasattr(result, 'pdf') else None,
        'mean': jlineMatrixToArray(result.mean) if hasattr(result, 'mean') else None
    }


def pfqn_conwayms_core(L, N, Z, S, options=None):
    """
    Conway multi-server algorithm core computation.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        S: Server configuration
        options: Algorithm options (optional)

    Returns:
        dict: Core performance metrics (XN, QN, RN, UN)
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwayms_coreKt.pfqn_conwayms_core(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwayms_coreKt.pfqn_conwayms_core(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S)
        )

    return {
        'XN': jlineMatrixToArray(result.XN) if hasattr(result, 'XN') else None,
        'QN': jlineMatrixToArray(result.QN) if hasattr(result, 'QN') else None,
        'RN': jlineMatrixToArray(result.RN) if hasattr(result, 'RN') else None,
        'UN': jlineMatrixToArray(result.UN) if hasattr(result, 'UN') else None
    }


def pfqn_conwayms_estimate(L, N, Z, S, options=None):
    """
    Conway multi-server algorithm with estimation.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        S: Server configuration
        options: Algorithm options (optional)

    Returns:
        dict: Estimated performance metrics
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwayms_estimateKt.pfqn_conwayms_estimate(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwayms_estimateKt.pfqn_conwayms_estimate(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S)
        )

    return {
        'XN': jlineMatrixToArray(result.XN) if hasattr(result, 'XN') else None,
        'QN': jlineMatrixToArray(result.QN) if hasattr(result, 'QN') else None,
        'RN': jlineMatrixToArray(result.RN) if hasattr(result, 'RN') else None,
        'UN': jlineMatrixToArray(result.UN) if hasattr(result, 'UN') else None
    }


def pfqn_conwayms_forwardmva(L, N, Z, S, options=None):
    """
    Conway multi-server algorithm with forward MVA.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        S: Server configuration
        options: Algorithm options (optional)

    Returns:
        dict: Performance metrics via forward MVA
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwayms_forwardmvaKt.pfqn_conwayms_forwardmva(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwayms_forwardmvaKt.pfqn_conwayms_forwardmva(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S)
        )

    return {
        'XN': jlineMatrixToArray(result.XN) if hasattr(result, 'XN') else None,
        'QN': jlineMatrixToArray(result.QN) if hasattr(result, 'QN') else None,
        'RN': jlineMatrixToArray(result.RN) if hasattr(result, 'RN') else None,
        'UN': jlineMatrixToArray(result.UN) if hasattr(result, 'UN') else None
    }


def pfqn_mu_ms_gnaux(L, N, Z, S, mu, options=None):
    """
    Multi-server service rate computation with auxiliary variables.

    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think time vector
        S: Server configuration
        mu: Service rates
        options: Algorithm options (optional)

    Returns:
        dict: Scaling factors and convergence information
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_mu_ms_gnauxKt.pfqn_mu_ms_gnaux(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S),
            jlineMatrixFromArray(mu), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_mu_ms_gnauxKt.pfqn_mu_ms_gnaux(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S),
            jlineMatrixFromArray(mu)
        )

    return {
        'scaling': jlineMatrixToArray(result.scaling) if hasattr(result, 'scaling') else None,
        'G': result.G if hasattr(result, 'G') else None,
        'converged': result.converged if hasattr(result, 'converged') else None
    }


def pfqn_nc(N, L, Z, options=None):
    """
    Normalizing Constant algorithm variant for product-form networks.
    
    Alternative interface to the normalizing constant method that
    computes performance measures using exact or approximate NC techniques.
    
    Args:
        N: Population vector
        L: Service demand matrix
        Z: Think time vector
        options: Optional solver options
        
    Returns:
        Performance metrics computed via normalizing constants
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_ncKt.pfqn_nc(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_ncKt.pfqn_nc(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN),
        'CN': jlineMatrixToArray(result.CN) if hasattr(result, 'CN') else None
    }


def pfqn_gld(L, N, mu, options=None):
    """
    Generalized Load-Dependent algorithm variant.
    
    Alternative interface to the GLD method with explicit
    service rate parameters for load-dependent stations.
    
    Args:
        L: Service demand matrix
        N: Population vector
        mu: Load-dependent service rate matrix
        options: Optional solver options
        
    Returns:
        Performance metrics for load-dependent network
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_gldKt.pfqn_gld(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(mu), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_gldKt.pfqn_gld(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(mu)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN)
    }


def pfqn_le(N, L, Z, mu, options=None):
    """
    Load Evaluation algorithm with service rates for product-form networks.
    
    Extended Load Evaluation method that includes explicit service rate
    parameters for more accurate modeling of multi-server stations.
    
    Args:
        N: Population vector
        L: Service demand matrix
        Z: Think time vector
        mu: Service rate matrix
        options: Optional solver options
        
    Returns:
        Performance metrics including detailed service rate effects
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_leKt.pfqn_le(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(mu), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_leKt.pfqn_le(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(mu)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN)
    }


def pfqn_conwayms(N, L, Z, S, options=None):
    """
    Conway Multi-Server algorithm with server specifications.
    
    Alternative interface for Conway's multi-server method that
    takes explicit server configuration parameters.
    
    Args:
        N: Population vector
        L: Service demand matrix
        Z: Think time vector
        S: Server configuration matrix
        options: Optional solver options
        
    Returns:
        Performance metrics with detailed server utilization
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwaymsKt.pfqn_conwayms(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_conwaymsKt.pfqn_conwayms(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(S)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN)
    }


def pfqn_cdfun(nvec, cdscaling, M):
    """
    Evaluate class-dependent (CD) scaling function for load-dependent analysis.

    .. note::
        This is a low-level internal function that's refactored to use index-based
        parameters instead of object-based parameters. The JAR function now takes
        cdscaling as a list indexed by station index rather than a map.

    Args:
        nvec: Per-class queue-length values matrix
        cdscaling: List of CD scaling functions indexed by station index
        M: Number of stations

    Returns:
        Matrix of scaling factors
    """
    result = jpype.JPackage('jline').api.pfqn.ld.Pfqn_cdfunKt.pfqn_cdfun(
        nvec, cdscaling, M
    )
    return result


def pfqn_nca(N, L, Z, options=None):
    """
    Normalizing Constant Approximation algorithm.
    
    Computes performance measures using approximation techniques
    for the normalizing constant in product-form networks.
    
    Args:
        N: Population vector.
        L: Service demand matrix.
        Z: Think time vector.
        options: Algorithm options (optional).
        
    Returns:
        dict: Performance measures (QN, UN, RN, TN).
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_ncaKt.pfqn_nca(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_ncaKt.pfqn_nca(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN)
    }


def pfqn_ncld(N, L, Z, mu, options=None):
    """
    Normalizing Constant for Load-Dependent networks.
    
    Computes performance measures using normalizing constant methods
    for networks with load-dependent service rates.
    
    Args:
        N: Population vector.
        L: Service demand matrix.
        Z: Think time vector.
        mu: Load-dependent service rates.
        options: Algorithm options (optional).
        
    Returns:
        dict: Performance measures (QN, UN, RN, TN).
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_ncldKt.pfqn_ncld(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(mu), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_ncldKt.pfqn_ncld(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), jlineMatrixFromArray(mu)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN)
    }


def pfqn_pff_delay(N, L, Z, options=None):
    """
    Product-form approximation for delay networks.

    Args:
        N: Population vector
        L: Service demand matrix
        Z: Think time vector
        options: Algorithm options (optional)

    Returns:
        dict: Performance metrics for delay networks
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_pff_delayKt.pfqn_pff_delay(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_pff_delayKt.pfqn_pff_delay(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN)
    }


def pfqn_sqni(N, L, Z, options=None):
    """
    Single Queue Network Iteration algorithm.

    Args:
        N: Population vector
        L: Service demand matrix
        Z: Think time vector
        options: Algorithm options (optional)

    Returns:
        dict: Performance metrics via SQNI method
    """
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_sqniKt.pfqn_sqni(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_sqniKt.pfqn_sqni(
            jlineMatrixFromArray(N), jlineMatrixFromArray(L),
            jlineMatrixFromArray(Z)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN)
    }


def pfqn_nc_sanitize(L, N, Z, options=None):
    java_options = jpype.JObject(options) if options is not None else None

    if java_options is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_nc_sanitizeKt.pfqn_nc_sanitize(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z), java_options
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_nc_sanitizeKt.pfqn_nc_sanitize(
            jlineMatrixFromArray(L), jlineMatrixFromArray(N),
            jlineMatrixFromArray(Z)
        )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN),
        'G': float(result.G) if hasattr(result, 'G') else None
    }


def pfqn_qzgblow(M, N):
    """
    Zahorjan-Gesbert lower bound algorithm.
    
    Computes lower bounds for queue lengths using the
    Zahorjan-Gesbert approximation method.
    
    Args:
        M: Service demand matrix.
        N: Population vector.
        
    Returns:
        float: Lower bound estimate.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_qzgblowKt.pfqn_qzgblow(
        jlineMatrixFromArray(M), jlineMatrixFromArray(N)
    )
    return float(result)


def pfqn_qzgbup(M, N):
    """
    Zahorjan-Gesbert upper bound algorithm.
    
    Computes upper bounds for queue lengths using the
    Zahorjan-Gesbert approximation method.
    
    Args:
        M: Service demand matrix.
        N: Population vector.
        
    Returns:
        float: Upper bound estimate.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_qzgbupKt.pfqn_qzgbup(
        jlineMatrixFromArray(M), jlineMatrixFromArray(N)
    )
    return float(result)


def pfqn_nc_sanitize(L, N, Z):
    """
    Sanitize parameters for normalizing constant computation.
    
    Preprocesses and validates input parameters to ensure numerical
    stability in normalizing constant algorithms.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        
    Returns:
        dict: Sanitized parameters with scaling information.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_nc_sanitizeKt.pfqn_nc_sanitize(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z)
    )

    return {
        'L_sanitized': jlineMatrixToArray(result.L),
        'N_sanitized': jlineMatrixToArray(result.N),
        'Z_sanitized': jlineMatrixToArray(result.Z),
        'scaling_factor': float(result.scalingFactor) if hasattr(result, 'scalingFactor') else None
    }


def pfqn_comomrm_ld(L, N, Z, S):
    """
    Co-moment matching algorithm for load-dependent networks.
    
    Applies moment matching techniques to product-form networks
    with load-dependent service stations.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        S: Server configuration matrix.
        
    Returns:
        dict: Performance measures including normalizing constant G.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_comomrm_ldKt.pfqn_comomrm_ld(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jlineMatrixFromArray(S)
    )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN),
        'G': float(result.G) if hasattr(result, 'G') else None
    }


def pfqn_mvaldmx_ec(L, N, Z, S):
    """
    MVA for Load-Dependent Mixed networks with Error Control.
    
    Mean Value Analysis for mixed networks with load-dependent
    service rates and enhanced error control mechanisms.
    
    Args:
        L: Service demand matrix.
        N: Population vector.
        Z: Think time vector.
        S: Server configuration matrix.
        
    Returns:
        dict: Performance measures with error-controlled computation.
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mvaldmx_ecKt.pfqn_mvaldmx_ec(
        jlineMatrixFromArray(L), jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z), jlineMatrixFromArray(S)
    )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN),
        'G': float(result.G) if hasattr(result, 'G') else None
    }


# Additional PFQN functions for complete API coverage

def pfqn_ab(L, N, Z=None, type='exact'):
    """
    Asymptotic bounds for product-form queueing networks.
    
    Args:
        L: Service demand matrix (stations x classes)
        N: Population vector
        Z: Think times (optional)
        type: Type of bounds ('exact' or 'approx')
        
    Returns:
        dict: Asymptotic bounds on throughput and response time
    """
    L_matrix = jlineMatrixFromArray(L)
    N_matrix = jlineMatrixFromArray(N)
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None
    
    if Z_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_abKt.pfqn_ab(
            L_matrix, N_matrix, Z_matrix, str(type)
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_abKt.pfqn_ab(
            L_matrix, N_matrix, str(type)
        )
    
    return {
        'throughput_lower': jlineMatrixToArray(result.getThroughputLower()),
        'throughput_upper': jlineMatrixToArray(result.getThroughputUpper()),
        'response_time_lower': jlineMatrixToArray(result.getResponseTimeLower()),
        'response_time_upper': jlineMatrixToArray(result.getResponseTimeUpper())
    }


def pfqn_le(L, N, Z=None):
    """
    Little's law equations for product-form networks.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think times (optional)
        
    Returns:
        dict: Performance metrics from Little's law
    """
    L_matrix = jlineMatrixFromArray(L)
    N_matrix = jlineMatrixFromArray(N)
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None
    
    if Z_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_leKt.pfqn_le(L_matrix, N_matrix, Z_matrix)
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_leKt.pfqn_le(L_matrix, N_matrix)
    
    return {
        'throughput': jlineMatrixToArray(result.getThroughput()),
        'queue_length': jlineMatrixToArray(result.getQueueLength())
    }


def pfqn_le_fpi(L, N, Z=None, max_iter=100):
    """
    Fixed-point iteration for Little's law equations.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think times (optional)
        max_iter: Maximum iterations
        
    Returns:
        dict: Converged performance metrics
    """
    L_matrix = jlineMatrixFromArray(L)
    N_matrix = jlineMatrixFromArray(N)
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None
    
    if Z_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_le_fpiKt.pfqn_le_fpi(
            L_matrix, N_matrix, Z_matrix, int(max_iter)
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_le_fpiKt.pfqn_le_fpi(
            L_matrix, N_matrix, int(max_iter)
        )
    
    return {
        'throughput': jlineMatrixToArray(result.getThroughput()),
        'queue_length': jlineMatrixToArray(result.getQueueLength()),
        'iterations': result.getIterations()
    }


def pfqn_le_fpiz(L, N, Z, max_iter=100):
    """
    Fixed-point iteration with think time adjustment.
    
    Args:
        L: Service demand matrix
        N: Population vector  
        Z: Think times
        max_iter: Maximum iterations
        
    Returns:
        dict: Converged performance metrics
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_le_fpizKt.pfqn_le_fpiz(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z),
        int(max_iter)
    )
    
    return {
        'throughput': jlineMatrixToArray(result.getThroughput()),
        'queue_length': jlineMatrixToArray(result.getQueueLength()),
        'iterations': result.getIterations()
    }


def pfqn_le_hessianz(L, N, Z):
    """
    Hessian computation for product-form networks with think times.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think times
        
    Returns:
        numpy.ndarray: Hessian matrix
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_le_hessianzKt.pfqn_le_hessianz(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(Z)
    )
    return jlineMatrixToArray(result)


def pfqn_mom(L, N, Z=None):
    """
    Method of moments for product-form networks.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think times (optional)
        
    Returns:
        dict: First and second moments
    """
    L_matrix = jlineMatrixFromArray(L)
    N_matrix = jlineMatrixFromArray(N)
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None
    
    if Z_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_momKt.pfqn_mom(L_matrix, N_matrix, Z_matrix)
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_momKt.pfqn_mom(L_matrix, N_matrix)
    
    return {
        'mean': jlineMatrixToArray(result.getMean()),
        'variance': jlineMatrixToArray(result.getVariance())
    }


def pfqn_mu_ms(L, N, m):
    """
    Service rates for multi-server product-form networks.
    
    Args:
        L: Service demand matrix
        N: Population vector
        m: Number of servers at each station
        
    Returns:
        numpy.ndarray: Effective service rates
    """
    result = jpype.JPackage('jline').api.pfqn.Pfqn_mu_msKt.pfqn_mu_ms(
        jlineMatrixFromArray(L),
        jlineMatrixFromArray(N),
        jlineMatrixFromArray(m)
    )
    return jlineMatrixToArray(result)


def pfqn_procomom2(L, N, Z=None):
    """
    Second-order product-form complementary moments.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think times (optional)
        
    Returns:
        dict: Second-order complementary moments
    """
    L_matrix = jlineMatrixFromArray(L)
    N_matrix = jlineMatrixFromArray(N)
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None
    
    if Z_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_procomom2Kt.pfqn_procomom2(
            L_matrix, N_matrix, Z_matrix
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_procomom2Kt.pfqn_procomom2(
            L_matrix, N_matrix
        )
    
    return {
        'moments': jlineMatrixToArray(result.getMoments()),
        'cross_moments': jlineMatrixToArray(result.getCrossMoments())
    }


def pfqn_schmidt(L, N, Z=None):
    """
    Schmidt's approximation for product-form networks.
    
    Args:
        L: Service demand matrix
        N: Population vector
        Z: Think times (optional)
        
    Returns:
        dict: Performance metrics using Schmidt's method
    """
    L_matrix = jlineMatrixFromArray(L)
    N_matrix = jlineMatrixFromArray(N)
    Z_matrix = jlineMatrixFromArray(Z) if Z is not None else None
    
    if Z_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_schmidtKt.pfqn_schmidt(
            L_matrix, N_matrix, Z_matrix
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.Pfqn_schmidtKt.pfqn_schmidt(
            L_matrix, N_matrix
        )
    
    return {
        'throughput': jlineMatrixToArray(result.getThroughput()),
        'queue_length': jlineMatrixToArray(result.getQueueLength()),
        'utilization': jlineMatrixToArray(result.getUtilization())
    }


# LCFS Queueing Network Functions

def pfqn_lcfsqn_ca(alpha, beta, N=None):
    """
    Convolution algorithm for multiclass LCFS queueing networks.

    Computes the normalizing constant for a 2-station closed queueing network with:
      - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
      - Station 2: LCFS-PR (LCFS with Preemption-Resume)

    Args:
        alpha: Vector of inverse service rates at station 1 (LCFS).
               alpha[r] = 1/mu(1,r) for class r
        beta: Vector of inverse service rates at station 2 (LCFS-PR).
              beta[r] = 1/mu(2,r) for class r
        N: Population vector, N[r] = number of jobs of class r.
           Default: ones(1,R) - one job per class

    Returns:
        tuple: (G, V) where G is the normalizing constant and V is the auxiliary term.

    Reference:
        G. Casale, "A family of multiclass LCFS queueing networks with
        order-dependent product-form solutions", QUESTA 2026.
    """
    alpha_matrix = jlineMatrixFromArray(alpha)
    beta_matrix = jlineMatrixFromArray(beta)
    N_matrix = jlineMatrixFromArray(N) if N is not None else None

    if N_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.lcfs.Pfqn_lcfsqn_caKt.pfqn_lcfsqn_ca(
            alpha_matrix, beta_matrix, N_matrix
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.lcfs.Pfqn_lcfsqn_caKt.pfqn_lcfsqn_ca(
            alpha_matrix, beta_matrix
        )

    return result.getG(), result.getV()


def pfqn_lcfsqn_mva(alpha, beta, N=None):
    """
    Mean Value Analysis for multiclass LCFS queueing networks.

    Computes performance metrics for a 2-station closed queueing network with:
      - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
      - Station 2: LCFS-PR (LCFS with Preemption-Resume)

    This implementation uses log-space arithmetic to prevent numerical
    underflow. The results are mathematically exact (up to floating-point
    precision) - no approximations are made.

    Args:
        alpha: Vector of inverse service rates at station 1 (LCFS).
               alpha[r] = 1/mu(1,r) for class r
        beta: Vector of inverse service rates at station 2 (LCFS-PR).
              beta[r] = 1/mu(2,r) for class r
        N: Population vector, N[r] = number of jobs of class r.
           Default: ones(1,R) - one job per class

    Returns:
        dict: Dictionary containing:
            - 'T': Throughput vector (1 x R)
            - 'Q': Queue length matrix (2 x R)
            - 'U': Utilization matrix (2 x R)
            - 'B': Back probability matrix (2 x R)

    Reference:
        G. Casale, "A family of multiclass LCFS queueing networks with
        order-dependent product-form solutions", QUESTA 2026.
    """
    alpha_matrix = jlineMatrixFromArray(alpha)
    beta_matrix = jlineMatrixFromArray(beta)
    N_matrix = jlineMatrixFromArray(N) if N is not None else None

    if N_matrix is not None:
        result = jpype.JPackage('jline').api.pfqn.lcfs.Pfqn_lcfsqn_mvaKt.pfqn_lcfsqn_mva(
            alpha_matrix, beta_matrix, N_matrix
        )
    else:
        result = jpype.JPackage('jline').api.pfqn.lcfs.Pfqn_lcfsqn_mvaKt.pfqn_lcfsqn_mva(
            alpha_matrix, beta_matrix
        )

    return {
        'T': jlineMatrixToArray(result.getT()),
        'Q': jlineMatrixToArray(result.getQ()),
        'U': jlineMatrixToArray(result.getU()),
        'B': jlineMatrixToArray(result.getB())
    }


def pfqn_lcfsqn_nc(alpha, beta, N):
    """
    Normalizing constant for multiclass LCFS queueing networks.

    Computes the normalizing constant using matrix permanent calculations
    for a 2-station closed queueing network with Station 1 using LCFS
    (Last-Come-First-Served, non-preemptive) and Station 2 using LCFS-PR
    (LCFS with Preemption-Resume).

    Args:
        alpha: Vector of inverse service rates at station 1 (LCFS).
               alpha[r] = 1/mu(1,r) for class r
        beta: Vector of inverse service rates at station 2 (LCFS-PR).
              beta[r] = 1/mu(2,r) for class r
        N: Population vector, N[r] = number of jobs of class r.

    Returns:
        tuple: (G, Ax) where G is the normalizing constant and Ax is the
               array of A matrices for each state.

    Reference:
        G. Casale, "A family of multiclass LCFS queueing networks with
        order-dependent product-form solutions", QUESTA 2026.
    """
    alpha_matrix = jlineMatrixFromArray(alpha)
    beta_matrix = jlineMatrixFromArray(beta)
    N_matrix = jlineMatrixFromArray(N)

    result = jpype.JPackage('jline').api.pfqn.lcfs.Pfqn_lcfsqn_ncKt.pfqn_lcfsqn_nc(
        alpha_matrix, beta_matrix, N_matrix
    )

    # Convert Ax array
    Ax = []
    java_Ax = result.getAx()
    for i in range(len(java_Ax)):
        Ax.append(jlineMatrixToArray(java_Ax[i]))

    return result.getG(), Ax


# =============================================================================
# Replicated Stations Support
# =============================================================================

def pfqn_unique(L, mu=None, gamma=None):
    """
    Identify and consolidate replicated (identical) stations.

    Analyzes the load matrix to detect stations with identical workloads
    across all classes. These replicated stations are consolidated into
    unique representatives, enabling more efficient MVA/RECAL analysis.

    The algorithm groups stations by their normalized load vectors:
    - Stations with identical L[i,:]/||L[i,:]|| are considered replicas
    - Returns a mapping from original to unique stations
    - Optionally processes mu (service rates) and gamma (visit ratios)

    Args:
        L: Load matrix (M x K) where M = stations, K = classes.
            L[i,j] = load at station i from class j.
        mu: Optional service rate matrix (M x K). If provided, returns
            consolidated service rates for unique stations.
        gamma: Optional visit ratio matrix (M x K). If provided, returns
            consolidated visit ratios for unique stations.

    Returns:
        dict: Dictionary containing:
            - 'L_unique': Consolidated load matrix for unique stations
            - 'mu_unique': Consolidated service rates (if mu provided)
            - 'gamma_unique': Consolidated visit ratios (if gamma provided)
            - 'mi': Multiplicity vector - count of replicas for each unique station
            - 'mapping': Index mapping from original to unique stations

    Examples:
        >>> # 4 stations with 2 pairs of replicas
        >>> L = [[1.0, 2.0], [1.0, 2.0], [3.0, 4.0], [3.0, 4.0]]
        >>> result = pfqn_unique(L)
        >>> print(result['L_unique'])  # 2 unique stations
        >>> print(result['mi'])        # [2, 2] - each appears twice
        >>> print(result['mapping'])   # [0, 0, 1, 1]
    """
    L_matrix = jlineMatrixFromArray(np.array(L, dtype=float))

    mu_matrix = None if mu is None else jlineMatrixFromArray(np.array(mu, dtype=float))
    gamma_matrix = None if gamma is None else jlineMatrixFromArray(np.array(gamma, dtype=float))

    result = jpype.JPackage('jline').api.pfqn.Pfqn_replicasKt.pfqn_unique(
        L_matrix, mu_matrix, gamma_matrix
    )

    response = {
        'L_unique': jlineMatrixToArray(result.getL_unique()),
        'mi': jlineMatrixToArray(result.getMi()).flatten().astype(int),
        'mapping': np.array(result.getMapping(), dtype=int)
    }

    if result.getMu_unique() is not None:
        response['mu_unique'] = jlineMatrixToArray(result.getMu_unique())
    if result.getGamma_unique() is not None:
        response['gamma_unique'] = jlineMatrixToArray(result.getGamma_unique())

    return response


def pfqn_expand(QN, UN, CN, mapping):
    """
    Expand performance metrics from unique stations back to original dimensions.

    After solving a consolidated model with unique stations (via pfqn_unique),
    this function expands the results back to the original station dimensions
    using the mapping array.

    Args:
        QN: Queue lengths for unique stations (M_unique x K).
        UN: Utilizations for unique stations (M_unique x K).
        CN: Response times for unique stations (M_unique x K).
        mapping: Index mapping from original to unique stations.
            mapping[i] = index of unique station that original station i maps to.

    Returns:
        tuple: (QN_expanded, UN_expanded, CN_expanded) - metrics at original dimensions.

    Examples:
        >>> # Expand from 2 unique stations back to 4 original
        >>> mapping = [0, 0, 1, 1]  # Stations 0,1 map to unique 0; 2,3 to unique 1
        >>> QN_unique = [[0.5, 0.3], [0.8, 0.4]]  # 2 unique x 2 classes
        >>> UN_unique = [[0.2, 0.1], [0.3, 0.15]]
        >>> CN_unique = [[1.0, 0.8], [1.5, 1.2]]
        >>> QN, UN, CN = pfqn_expand(QN_unique, UN_unique, CN_unique, mapping)
        >>> print(QN.shape)  # (4, 2) - back to original dimensions
    """
    QN_matrix = jlineMatrixFromArray(np.array(QN, dtype=float))
    UN_matrix = jlineMatrixFromArray(np.array(UN, dtype=float))
    CN_matrix = jlineMatrixFromArray(np.array(CN, dtype=float))
    mapping_array = jpype.JArray(jpype.JInt)(list(mapping))

    result = jpype.JPackage('jline').api.pfqn.Pfqn_replicasKt.pfqn_expand(
        QN_matrix, UN_matrix, CN_matrix, mapping_array
    )

    return (
        jlineMatrixToArray(result.getFirst()),
        jlineMatrixToArray(result.getSecond()),
        jlineMatrixToArray(result.getThird())
    )


def pfqn_combine_mi(mi, mapping, M_unique):
    """
    Combine user-provided multiplicities with detected replica mappings.

    This utility function merges explicit server counts (multiplicities)
    with the replica detection from pfqn_unique, creating a combined
    multiplicity vector for the unique stations.

    Args:
        mi: Original multiplicity vector (M x 1). mi[i] = number of servers
            at original station i.
        mapping: Index mapping from original to unique stations.
        M_unique: Number of unique stations after consolidation.

    Returns:
        ndarray: Combined multiplicity vector (M_unique x 1) where each
            element is the sum of multiplicities for that unique station.

    Examples:
        >>> # 4 stations with multiplicities, mapped to 2 unique
        >>> mi = [2, 2, 3, 3]  # Each station has these server counts
        >>> mapping = [0, 0, 1, 1]  # Pairs map to same unique station
        >>> M_unique = 2
        >>> mi_combined = pfqn_combine_mi(mi, mapping, M_unique)
        >>> print(mi_combined)  # [4, 6] - summed for each unique station
    """
    mi_matrix = jlineMatrixFromArray(np.array(mi, dtype=float).flatten())
    mapping_array = jpype.JArray(jpype.JInt)(list(mapping))

    result = jpype.JPackage('jline').api.pfqn.Pfqn_replicasKt.pfqn_combine_mi(
        mi_matrix, mapping_array, int(M_unique)
    )

    return jlineMatrixToArray(result).flatten()
