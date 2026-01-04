
"""
Markov Chain analysis functions.

This module provides unified access to both continuous-time (CTMC) and 
discrete-time (DTMC) Markov chain analysis functions. It serves as a
consolidated interface for Markov chain computations used throughout LINE.

Functions include both CTMC and DTMC methods:
- ctmc_makeinfgen: Construct infinitesimal generator
- ctmc_solve: Steady-state CTMC analysis  
- dtmc_solve: Steady-state DTMC analysis
- Transient analysis methods
- Simulation functions
- Time-reversal and other transformations

This module aggregates functions from both ctmc and dtmc modules.
"""

import jpype
import numpy as np
from .. import jlineMatrixToArray, jlineMatrixFromArray


def ctmc_makeinfgen(birth_rates, death_rates, max_states=None):
    """
    Construct infinitesimal generator for birth-death process.
    
    Creates a tridiagonal generator matrix for a birth-death CTMC
    with specified birth and death rates.
    
    Args:
        birth_rates: Array of birth rates (transitions i -> i+1).
        death_rates: Array of death rates (transitions i -> i-1).
        max_states: Maximum number of states (optional).
        
    Returns:
        numpy.ndarray: Infinitesimal generator matrix Q.
    """
    birth_rates = np.array(birth_rates)
    death_rates = np.array(death_rates)

    if max_states is None:
        max_states = max(len(birth_rates), len(death_rates))

    if len(birth_rates) < max_states:
        birth_rates = np.pad(birth_rates, (0, max_states - len(birth_rates)), 'constant')
    if len(death_rates) < max_states:
        death_rates = np.pad(death_rates, (0, max_states - len(death_rates)), 'constant')

    Q = np.zeros((max_states, max_states))

    for i in range(max_states):
        if i < max_states - 1 and birth_rates[i] > 0:
            Q[i, i+1] = birth_rates[i]

        if i > 0 and death_rates[i] > 0:
            Q[i, i-1] = death_rates[i]

        Q[i, i] = -np.sum(Q[i, :])

    return Q


def ctmc_solve(Q):
    """
    Solve for steady-state probabilities of a CTMC.
    
    Computes the stationary distribution π by solving πQ = 0
    with normalization constraint.
    
    Args:
        Q: Infinitesimal generator matrix.
        
    Returns:
        numpy.ndarray: Steady-state probability distribution.
    """
    Q = np.array(Q)
    n = Q.shape[0]

    java_Q = jlineMatrixFromArray(Q)
    java_result = jpype.JPackage('jline').api.mc.Ctmc_solveKt.ctmc_solve(java_Q)

    return jlineMatrixToArray(java_result).flatten()


def ctmc_transient(Q, initial_dist, time_points):
    """
    Compute transient probabilities of a CTMC.
    
    Calculates time-dependent state probabilities π(t) for
    specified time points using matrix exponential methods.
    
    Args:
        Q: Infinitesimal generator matrix.
        initial_dist: Initial probability distribution π(0).
        time_points: Array of time points to evaluate.
        
    Returns:
        numpy.ndarray: Transient probabilities at each time point.
    """
    Q = np.array(Q)
    initial_dist = np.array(initial_dist)
    time_points = np.array(time_points)

    java_Q = jlineMatrixFromArray(Q)
    java_initial = jlineMatrixFromArray(initial_dist.reshape(1, -1))
    java_times = jpype.JArray(jpype.JDouble)(time_points)

    java_result = jpype.JPackage('jline').api.mc.Ctmc_transientKt.ctmc_transient(
        java_Q, java_initial, java_times
    )

    return jlineMatrixToArray(java_result)


def ctmc_simulate(Q, initial_state, max_time, max_events=10000):
    """
    Simulate CTMC sample path.
    
    Generates a realization of the continuous-time Markov chain
    using the Gillespie algorithm or similar Monte Carlo method.
    
    Args:
        Q: Infinitesimal generator matrix.
        initial_state: Starting state (integer index).
        max_time: Maximum simulation time.
        max_events: Maximum number of transitions (default: 10000).
        
    Returns:
        dict: Simulation results with 'states', 'times', and 'sojourn_times'.
    """
    Q = np.array(Q)

    java_Q = jlineMatrixFromArray(Q)
    java_result = jpype.JPackage('jline').api.mc.Ctmc_simulateKt.ctmc_simulate(
        java_Q,
        jpype.JInt(initial_state),
        jpype.JDouble(max_time),
        jpype.JInt(max_events)
    )

    return {
        'states': np.array(java_result.getStates()),
        'times': np.array(java_result.getTimes()),
        'sojourn_times': np.array(java_result.getSojournTimes())
    }


def ctmc_randomization(Q, initial_dist, time_points, precision=1e-10):
    """
    Compute CTMC transient probabilities using randomization.
    
    Uses Jensen's randomization method to compute transient
    probabilities by converting the CTMC to a uniformized DTMC.
    
    Args:
        Q: Infinitesimal generator matrix.
        initial_dist: Initial probability distribution.
        time_points: Array of time points to evaluate.
        precision: Numerical precision for truncation (default: 1e-10).
        
    Returns:
        numpy.ndarray: Transient probabilities via randomization.
    """
    Q = np.array(Q)
    initial_dist = np.array(initial_dist)
    time_points = np.array(time_points)

    java_Q = jlineMatrixFromArray(Q)
    java_initial = jlineMatrixFromArray(initial_dist.reshape(1, -1))
    java_times = jpype.JArray(jpype.JDouble)(time_points)

    java_result = jpype.JPackage('jline').api.mc.Ctmc_randomizationKt.ctmc_randomization(
        java_Q, java_initial, java_times, jpype.JDouble(precision)
    )

    return jlineMatrixToArray(java_result)


def ctmc_uniformization(Q, lambda_rate=None):
    """
    Uniformize CTMC generator matrix.
    
    Converts CTMC to an equivalent uniformized discrete-time chain
    for numerical analysis and simulation purposes.
    
    Args:
        Q: Infinitesimal generator matrix.
        lambda_rate: Uniformization rate (optional, auto-computed if None).
        
    Returns:
        dict: Contains uniformized transition matrix 'P' and rate 'lambda'.
    """
    Q = np.array(Q)

    if lambda_rate is None:
        lambda_rate = -np.min(np.diag(Q))

    n = Q.shape[0]
    I = np.eye(n)
    P = I + Q / lambda_rate

    return {
        'P': P,
        'lambda': lambda_rate
    }


def ctmc_stochcomp(Q, keep_states, eliminate_states):
    """
    Compute stochastic complement of CTMC.

    Args:
        Q: Infinitesimal generator matrix.
        keep_states: States to retain in reduced model.
        eliminate_states: States to eliminate from model.

    Returns:
        numpy.ndarray: Reduced generator matrix.
    """
    Q = np.array(Q)
    keep_states = np.array(keep_states, dtype=int)
    eliminate_states = np.array(eliminate_states, dtype=int)

    java_Q = jlineMatrixFromArray(Q)
    java_keep = jpype.JArray(jpype.JInt)(keep_states)
    java_eliminate = jpype.JArray(jpype.JInt)(eliminate_states)

    java_result = jpype.JPackage('jline').api.mc.Ctmc_stochcompKt.ctmc_stochcomp(
        java_Q, java_keep, java_eliminate
    )

    return jlineMatrixToArray(java_result)


def ctmc_timereverse(Q, pi=None):
    """
    Compute time-reversed CTMC generator.

    Args:
        Q: Original infinitesimal generator matrix.
        pi: Steady-state distribution (optional, computed if None).

    Returns:
        numpy.ndarray: Time-reversed generator matrix.
    """
    Q = np.array(Q)

    if pi is None:
        pi = ctmc_solve(Q)
    else:
        pi = np.array(pi)

    java_Q = jlineMatrixFromArray(Q)
    java_pi = jlineMatrixFromArray(pi.reshape(1, -1))

    java_result = jpype.JPackage('jline').api.mc.Ctmc_timereverseKt.ctmc_timereverse(
        java_Q, java_pi
    )

    return jlineMatrixToArray(java_result)


def ctmc_rand(n, density=0.3, max_rate=10.0):
    """
    Generate random CTMC generator matrix.

    Args:
        n: Number of states.
        density: Sparsity density (default: 0.3).
        max_rate: Maximum transition rate (default: 10.0).

    Returns:
        numpy.ndarray: Random infinitesimal generator matrix.
    """
    java_result = jpype.JPackage('jline').api.mc.Ctmc_randKt.ctmc_rand(
        jpype.JInt(n),
        jpype.JDouble(density),
        jpype.JDouble(max_rate)
    )

    return jlineMatrixToArray(java_result)



def dtmc_solve(P):
    """
    Solve for steady-state probabilities of DTMC.

    Args:
        P: Transition probability matrix.

    Returns:
        numpy.ndarray: Steady-state probability distribution.
    """
    P = np.array(P)

    java_P = jlineMatrixFromArray(P)
    java_result = jpype.JPackage('jline').api.mc.Dtmc_solveKt.dtmc_solve(java_P)

    return jlineMatrixToArray(java_result).flatten()


def dtmc_makestochastic(A):
    """
    Convert matrix to stochastic transition matrix.

    Args:
        A: Input matrix to normalize.

    Returns:
        numpy.ndarray: Row-stochastic matrix.
    """
    A = np.array(A, dtype=float)

    row_sums = np.sum(A, axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1

    return A / row_sums


def dtmc_isfeasible(P, tolerance=1e-10):
    """
    Check if matrix is a valid DTMC transition matrix.

    Args:
        P: Candidate transition matrix.
        tolerance: Numerical tolerance (default: 1e-10).

    Returns:
        bool: True if matrix is valid DTMC.
    """
    P = np.array(P)

    if np.any(P < 0):
        return False

    row_sums = np.sum(P, axis=1)
    if np.any(np.abs(row_sums - 1.0) > tolerance):
        return False

    return True


def dtmc_simulate(P, initial_state, num_steps):
    """
    Simulate DTMC sample path.

    Args:
        P: Transition probability matrix.
        initial_state: Starting state index.
        num_steps: Number of simulation steps.

    Returns:
        numpy.ndarray: Sequence of visited states.
    """
    P = np.array(P)

    java_P = jlineMatrixFromArray(P)
    java_result = jpype.JPackage('jline').api.mc.Dtmc_simulateKt.dtmc_simulate(
        java_P,
        jpype.JInt(initial_state),
        jpype.JInt(num_steps)
    )

    return np.array(java_result)


def dtmc_rand(n, density=0.5):
    """
    Generate random DTMC transition matrix.

    Args:
        n: Number of states.
        density: Sparsity density (default: 0.5).

    Returns:
        numpy.ndarray: Random transition probability matrix.
    """
    java_result = jpype.JPackage('jline').api.mc.Dtmc_randKt.dtmc_rand(
        jpype.JInt(n),
        jpype.JDouble(density)
    )

    return jlineMatrixToArray(java_result)


def dtmc_stochcomp(P, keep_states, eliminate_states):
    """
    Compute stochastic complement of DTMC.

    Args:
        P: Transition probability matrix.
        keep_states: States to retain in reduced model.
        eliminate_states: States to eliminate from model.

    Returns:
        numpy.ndarray: Reduced transition matrix.
    """
    P = np.array(P)
    keep_states = np.array(keep_states, dtype=int)
    eliminate_states = np.array(eliminate_states, dtype=int)

    java_P = jlineMatrixFromArray(P)
    java_keep = jpype.JArray(jpype.JInt)(keep_states)
    java_eliminate = jpype.JArray(jpype.JInt)(eliminate_states)

    java_result = jpype.JPackage('jline').api.mc.Dtmc_stochcompKt.dtmc_stochcomp(
        java_P, java_keep, java_eliminate
    )

    return jlineMatrixToArray(java_result)


def dtmc_timereverse(P, pi=None):
    """
    Compute time-reversed DTMC transition matrix.

    Args:
        P: Original transition matrix.
        pi: Steady-state distribution (optional, computed if None).

    Returns:
        numpy.ndarray: Time-reversed transition matrix.
    """
    P = np.array(P)

    if pi is None:
        pi = dtmc_solve(P)
    else:
        pi = np.array(pi)

    java_P = jlineMatrixFromArray(P)
    java_pi = jlineMatrixFromArray(pi.reshape(1, -1))

    java_result = jpype.JPackage('jline').api.mc.Dtmc_timereverseKt.dtmc_timereverse(
        java_P, java_pi
    )

    return jlineMatrixToArray(java_result)

# Additional CTMC/MC functions for complete API coverage

def ctmc_courtois(Q, epsilon=1e-6):
    """
    Courtois decomposition for nearly completely decomposable CTMCs.
    
    Args:
        Q: Infinitesimal generator matrix
        epsilon: Decomposition threshold (default: 1e-6)
        
    Returns:
        dict: Decomposition results with aggregated chain and coupling matrices
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_courtoisKt.ctmc_courtois(
        jlineMatrixFromArray(Q),
        float(epsilon)
    )
    return {
        'aggregated': jlineMatrixToArray(result.getAggregated()),
        'coupling': jlineMatrixToArray(result.getCoupling()),
        'blocks': [jlineMatrixToArray(b) for b in result.getBlocks()]
    }


def ctmc_kms(Q, blocks):
    """
    Koury-McAllister-Stewart aggregation-disaggregation method for CTMCs.
    
    Args:
        Q: Infinitesimal generator matrix
        blocks: Block structure specification
        
    Returns:
        numpy.ndarray: Steady-state probability vector
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_kmsKt.ctmc_kms(
        jlineMatrixFromArray(Q),
        jlineMatrixFromArray(blocks)
    )
    return jlineMatrixToArray(result)


def ctmc_multi(Q, levels):
    """
    Multi-level aggregation method for CTMCs.
    
    Args:
        Q: Infinitesimal generator matrix
        levels: Number of aggregation levels
        
    Returns:
        dict: Multi-level aggregation results
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_multiKt.ctmc_multi(
        jlineMatrixFromArray(Q),
        int(levels)
    )
    return {
        'probabilities': jlineMatrixToArray(result.getProbabilities()),
        'levels': result.getLevels()
    }


def ctmc_pseudostochcomp(Q, partition):
    """
    Pseudo-stochastic complementation for CTMCs.
    
    Args:
        Q: Infinitesimal generator matrix
        partition: State partition specification
        
    Returns:
        numpy.ndarray: Reduced generator matrix
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_pseudostochcompKt.ctmc_pseudostochcomp(
        jlineMatrixFromArray(Q),
        jlineMatrixFromArray(partition)
    )
    return jlineMatrixToArray(result)


def ctmc_relsolve(Q, ref_states):
    """
    Equilibrium distribution of a continuous-time Markov chain re-normalized with respect to
    reference states.
    
    Args:
        Q: Infinitesimal generator matrix
        ref_states: Reference state indices
        
    Returns:
        numpy.ndarray: Re-normalized steady-state probabilities
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_relsolveKt.ctmc_relsolve(
        jlineMatrixFromArray(Q),
        jlineMatrixFromArray(ref_states)
    )
    return jlineMatrixToArray(result)


def ctmc_solve_reducible(Q):
    """
    Solve reducible CTMCs by converting to DTMC via randomization.
    
    Args:
        Q: Infinitesimal generator matrix (possibly reducible)
        
    Returns:
        numpy.ndarray: Steady-state probability vector
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_solve_reducibleKt.ctmc_solve_reducible(
        jlineMatrixFromArray(Q)
    )
    return jlineMatrixToArray(result)


def ctmc_stmonotone(Q):
    """
    Computes the stochastically monotone upper bound for a CTMC.
    
    Args:
        Q: Infinitesimal generator matrix
        
    Returns:
        numpy.ndarray: Stochastically monotone upper bound generator
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_stmonotoneKt.ctmc_stmonotone(
        jlineMatrixFromArray(Q)
    )
    return jlineMatrixToArray(result)


def ctmc_takahashi(Q, epsilon=1e-8):
    """
    Takahashi's aggregation-disaggregation method for CTMCs.
    
    Args:
        Q: Infinitesimal generator matrix
        epsilon: Convergence threshold (default: 1e-8)
        
    Returns:
        numpy.ndarray: Steady-state probability vector
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_takahashiKt.ctmc_takahashi(
        jlineMatrixFromArray(Q),
        float(epsilon)
    )
    return jlineMatrixToArray(result)


def ctmc_testpf_kolmogorov(Q):
    """
    Test if a CTMC has product form using Kolmogorov's criteria.
    
    Args:
        Q: Infinitesimal generator matrix
        
    Returns:
        bool: True if the CTMC has product form
    """
    result = jpype.JPackage('jline').api.mc.Ctmc_testpf_kolmogorovKt.ctmc_testpf_kolmogorov(
        jlineMatrixFromArray(Q)
    )
    return bool(result)


def dtmc_solve_reducible(P):
    """
    Solve reducible DTMCs.
    
    Args:
        P: Transition probability matrix (possibly reducible)
        
    Returns:
        numpy.ndarray: Steady-state probability vector
    """
    result = jpype.JPackage('jline').api.mc.Dtmc_solve_reducibleKt.dtmc_solve_reducible(
        jlineMatrixFromArray(P)
    )
    return jlineMatrixToArray(result)


def dtmc_uniformization(P, steps):
    """
    DTMC uniformization analysis.
    
    Args:
        P: Transition probability matrix
        steps: Number of uniformization steps
        
    Returns:
        dict: Uniformization results with probability vector and iterations
    """
    result = jpype.JPackage('jline').api.mc.Dtmc_uniformizationKt.dtmc_uniformization(
        jlineMatrixFromArray(P),
        int(steps)
    )
    return {
        'probabilities': jlineMatrixToArray(result.getProbabilities()),
        'iterations': result.getIterations()
    }
