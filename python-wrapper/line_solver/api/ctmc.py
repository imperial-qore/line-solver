
"""
Continuous-Time Markov Chain (CTMC) analysis algorithms.

This module provides low-level functions for analyzing continuous-time
Markov chains, including transient and steady-state analysis, uniformization,
simulation, and various computational methods.

These functions support the CTMC solver and other analytical methods.
"""

import jpype

from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def ctmc_uniformization(pi0, Q, t):
    """
    Compute transient probabilities using uniformization method.
    
    Args:
        pi0: Initial probability distribution.
        Q: Infinitesimal generator matrix.
        t: Time points for transient analysis.
        
    Returns:
        numpy.ndarray: Transient probability distributions.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Ctmc_uniformizationKt.ctmc_uniformization(
            jlineMatrixFromArray(pi0),
            jlineMatrixFromArray(Q),
            jlineMatrixFromArray(t)
        )
    )


def ctmc_timereverse(matrix):
    """
    Compute time-reversed CTMC.
    
    Constructs the time-reversed continuous-time Markov chain using
    the detailed balance equations and steady-state probabilities.
    
    Args:
        matrix: Original infinitesimal generator matrix.
        
    Returns:
        numpy.ndarray: Time-reversed generator matrix.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Ctmc_timereverseKt.ctmc_timereverse(
            jlineMatrixFromArray(matrix)
        )
    )


def ctmc_makeinfgen(matrix):
    """
    Construct infinitesimal generator matrix.
    
    Converts a rate matrix to a proper infinitesimal generator by
    setting diagonal elements to negative row sums.
    
    Args:
        matrix: Rate matrix with non-negative off-diagonal elements.
        
    Returns:
        numpy.ndarray: Infinitesimal generator matrix.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Ctmc_makeinfgenKt.ctmc_makeinfgen(
            jlineMatrixFromArray(matrix)
        )
    )


def ctmc_solve(matrix):
    """
    Solve for steady-state probabilities of a CTMC.
    
    Computes the stationary distribution of a continuous-time Markov chain
    by solving the system πQ = 0 with the normalization constraint.
    
    Args:
        matrix: Infinitesimal generator matrix Q.
        
    Returns:
        numpy.ndarray: Steady-state probability distribution π.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Ctmc_solveKt.ctmc_solve(
            jlineMatrixFromArray(matrix)
        )
    )


def ctmc_transient(Q, pi0=None, t0=None, t1=None):
    """
    Compute transient probabilities of a CTMC.

    Args:
        Q: Infinitesimal generator matrix.
        pi0: Initial probability distribution (optional).
        t0: Start time (optional).
        t1: End time (required).

    Returns:
        tuple: (times, probabilities) arrays.
    """
    Q_matrix = jlineMatrixFromArray(Q)

    if pi0 is not None and t0 is not None and t1 is not None:
        result = jpype.JPackage('jline').api.mc.Ctmc_transientKt.ctmc_transient(
            Q_matrix, jlineMatrixFromArray(pi0), t0, t1
        )
    elif pi0 is not None and t1 is not None:
        result = jpype.JPackage('jline').api.mc.Ctmc_transientKt.ctmc_transient(
            Q_matrix, jlineMatrixFromArray(pi0), t1
        )
    elif t1 is not None:
        result = jpype.JPackage('jline').api.mc.Ctmc_transientKt.ctmc_transient(
            Q_matrix, t1
        )
    else:
        raise ValueError("t1 (end time) must be provided")

    times = list(result.first)
    probabilities = [list(prob_array) for prob_array in result.second]

    return times, probabilities


def ctmc_simulate(Q, pi0=None, n=1000):
    """
    Simulate a CTMC using Monte Carlo methods.

    Args:
        Q: Infinitesimal generator matrix.
        pi0: Initial distribution (optional).
        n: Number of simulation steps (default: 1000).

    Returns:
        tuple: (states, sojourn_times) arrays.
    """
    Q_matrix = jlineMatrixFromArray(Q)

    if pi0 is not None:
        pi0_array = jpype.JArray(jpype.JDouble)(len(pi0))
        for i, val in enumerate(pi0):
            pi0_array[i] = float(val)
    else:
        pi0_array = None

    result = jpype.JPackage('jline').api.mc.Ctmc_simulateKt.ctmc_simulate(
        Q_matrix, pi0_array, jpype.JInt(n)
    )

    states = list(result.states)
    sojourn_times = list(result.sojournTimes)

    return states, sojourn_times


def ctmc_rand(length):
    """
    Generate random CTMC generator matrix.
    
    Creates a random infinitesimal generator matrix suitable for use
    as a continuous-time Markov chain generator.
    
    Args:
        length: Size of the square matrix (number of states).
        
    Returns:
        numpy.ndarray: Random infinitesimal generator matrix.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Ctmc_randKt.ctmc_rand(jpype.JInt(length))
    )


def ctmc_ssg(sn, options):
    """
    Generate state space for CTMC from service network.

    Args:
        sn: Service network model.
        options: Generation options.

    Returns:
        dict: State space information including matrices and mappings.
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_ssgKt.ctmc_ssg(sn, options)

    return {
        'state_space': jlineMatrixToArray(result.stateSpace),
        'state_space_aggr': jlineMatrixToArray(result.stateSpaceAggr),
        'state_space_hashed': jlineMatrixToArray(result.stateSpaceHashed),
        'node_state_space': {node: jlineMatrixToArray(space)
                           for node, space in result.nodeStateSpace.items()},
        'sn': result.sn
    }


def ctmc_stochcomp(Q, I_list=None):
    """
    Compute stochastic complementarity decomposition of CTMC.

    Args:
        Q: Infinitesimal generator matrix.
        I_list: Index list for decomposition (optional).

    Returns:
        dict: Decomposed matrices (S, Q11, Q12, Q21, Q22, T).
    """
    Q_matrix = jlineMatrixFromArray(Q)

    if I_list is None:
        I_list = []

    java_list = jpype.java.util.ArrayList()
    for idx in I_list:
        java_list.add(jpype.JDouble(float(idx)) if idx is not None else None)

    result = jpype.JPackage('jline').api.mc.Ctmc_stochcompKt.ctmc_stochcomp(
        Q_matrix, java_list
    )

    return {
        'S': jlineMatrixToArray(result.S),
        'Q11': jlineMatrixToArray(result.Q11),
        'Q12': jlineMatrixToArray(result.Q12),
        'Q21': jlineMatrixToArray(result.Q21),
        'Q22': jlineMatrixToArray(result.Q22),
        'T': jlineMatrixToArray(result.T)
    }

def ctmc_ssg_reachability(sn, options=None):
    """
    Generate reachable state space for CTMC from service network.

    Args:
        sn: Service network model.
        options: Reachability analysis options (optional).

    Returns:
        dict: Reachable state space information.
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_ssg_reachabilityKt.ctmc_ssg_reachability(
        sn, options
    )

    node_state_space = {}
    if result.nodeStateSpace:
        for node_entry in result.nodeStateSpace.entrySet():
            node_key = node_entry.getKey()
            node_value = jlineMatrixToArray(node_entry.getValue())
            node_state_space[node_key] = node_value

    return {
        'state_space': jlineMatrixToArray(result.stateSpace),
        'state_space_aggr': jlineMatrixToArray(result.stateSpaceAggr),
        'state_space_hashed': jlineMatrixToArray(result.stateSpaceHashed),
        'node_state_space': node_state_space,
        'sn': result.sn
    }


def ctmc_randomization(Q, q=None):
    """
    Apply randomization (uniformization) to CTMC.

    Args:
        Q: Infinitesimal generator matrix.
        q: Randomization rate (optional, auto-computed if None).

    Returns:
        tuple: (P_matrix, q_rate) - transition matrix and rate.
    """
    Q_matrix = jlineMatrixFromArray(Q)

    if q is not None:
        result = jpype.JPackage('jline').api.mc.Ctmc_randomizationKt.ctmc_randomization(
            Q_matrix, jpype.JDouble(q)
        )
    else:
        result = jpype.JPackage('jline').api.mc.Ctmc_randomizationKt.ctmc_randomization(
            Q_matrix, None
        )

    P_matrix = jlineMatrixToArray(result.getFirst())
    q_rate = result.getSecond()

    return P_matrix, q_rate


def ctmc_krylov(Q, pi0, t, options=None):
    """
    Compute CTMC transient probabilities using Krylov subspace methods.

    Args:
        Q: Infinitesimal generator matrix.
        pi0: Initial probability distribution.
        t: Time point.
        options: Krylov method options (optional).

    Returns:
        numpy.ndarray: Transient probability distribution at time t.
    """
    Q_matrix = jlineMatrixFromArray(Q)
    pi0_matrix = jlineMatrixFromArray(pi0)

    if options is not None:
        result = jpype.JPackage('jline').api.mc.Ctmc_krylovKt.ctmc_krylov(
            Q_matrix, pi0_matrix, jpype.JDouble(t), jpype.JObject(options)
        )
    else:
        result = jpype.JPackage('jline').api.mc.Ctmc_krylovKt.ctmc_krylov(
            Q_matrix, pi0_matrix, jpype.JDouble(t)
        )

    return jlineMatrixToArray(result)

