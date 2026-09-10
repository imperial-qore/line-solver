"""
QMAM (Queueing Markov Analytical Methods) solvers.

Provides analytical solutions for specialized queueing systems:
- MAP/MAP/1 continuous-time queues
- MAP/M/c multiserver queues
- MMAP[K]/PH[K]/1 multiclass queues
- PH/PH/1 phase-type queues
- RAP/RAP/1 rational arrival process queues
- Sylvester equation solvers
"""

import numpy as np
from line_solver import native_to_array


# ========== MAP/MAP/1 QUEUE ==========

def lib_qmam_ct_map_map_1_steady_state(D0, D1):
    """Compute steady-state distribution for MAP/MAP/1 queue.

    Args:
        D0 (array-like): Arrival process generator matrix
        D1 (array-like): Arrival process transition matrix

    Returns:
        dict: Steady-state probability vector and related metrics
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = Queue.steadyState(native_to_array(D0), native_to_array(D1))
        return native_to_array(result) if hasattr(result, 'toArray') else result
    except Exception as e:
        raise RuntimeError(f"MAP/MAP/1 steady-state computation failed: {str(e)}")


def lib_qmam_ct_map_map_1_fundamental_matrix(D0, D1):
    """Compute fundamental matrix for MAP/MAP/1 queue.

    Args:
        D0 (array-like): Arrival process generator matrix
        D1 (array-like): Arrival process transition matrix

    Returns:
        ndarray: Fundamental matrix R
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = Queue.fundamentalMatrix(native_to_array(D0), native_to_array(D1))
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"MAP/MAP/1 fundamental matrix computation failed: {str(e)}")


# ========== MAP/M/c QUEUE ==========

def lib_qmam_ct_map_m_c_steady_state(D0, D1, mu, c):
    """Compute steady-state distribution for MAP/M/c queue.

    Args:
        D0 (array-like): Arrival process generator matrix
        D1 (array-like): Arrival process transition matrix
        mu (float): Service rate
        c (int): Number of servers

    Returns:
        dict: Steady-state probability vector
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    try:
        result = Queue.steadyState(
            native_to_array(D0), native_to_array(D1),
            float(mu), int(c)
        )
        return native_to_array(result) if hasattr(result, 'toArray') else result
    except Exception as e:
        raise RuntimeError(f"MAP/M/c steady-state computation failed: {str(e)}")


# ========== MMAP[K]/PH[K]/1 QUEUE ==========

def lib_qmam_ct_mmapk_phk_1_steady_state(D, D1, alpha, A):
    """Compute steady-state for MMAP[K]/PH[K]/1 queue.

    Args:
        D (list): List of K arrival process matrices
        D1 (array-like): Aggregate arrival transition matrix
        alpha (array-like): Service distribution initial vector
        A (array-like): Service distribution generator

    Returns:
        dict: Steady-state probability distribution
    """
    try:
        D_arr = [native_to_array(np.asarray(d)) for d in D]
        result = Queue.steadyState(
            D_arr, native_to_array(np.asarray(D1)),
            native_to_array(np.asarray(alpha)),
            native_to_array(np.asarray(A))
        )
        return native_to_array(result) if hasattr(result, 'toArray') else result
    except Exception as e:
        raise RuntimeError(f"MMAP[K]/PH[K]/1 steady-state computation failed: {str(e)}")


# ========== PH/PH/1 QUEUE ==========

def lib_qmam_ct_ph_ph_1_steady_state(alpha_a, A_a, alpha_s, A_s):
    """Compute steady-state distribution for PH/PH/1 queue.

    Args:
        alpha_a (array-like): Arrival phase-type initial vector
        A_a (array-like): Arrival phase-type generator
        alpha_s (array-like): Service phase-type initial vector
        A_s (array-like): Service phase-type generator

    Returns:
        ndarray: Steady-state probability vector
    """
    alpha_a = np.asarray(alpha_a)
    A_a = np.asarray(A_a)
    alpha_s = np.asarray(alpha_s)
    A_s = np.asarray(A_s)
    try:
        result = Queue.steadyState(
            native_to_array(alpha_a), native_to_array(A_a),
            native_to_array(alpha_s), native_to_array(A_s)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"PH/PH/1 steady-state computation failed: {str(e)}")


# ========== RAP/RAP/1 QUEUE ==========

def lib_qmam_ct_rap_rap_1_steady_state(D0, D1, D2):
    """Compute steady-state for RAP/RAP/1 queue.

    Args:
        D0 (array-like): Arrival RAP generator matrix
        D1 (array-like): Arrival RAP transition matrix
        D2 (array-like): Service RAP matrix (combined)

    Returns:
        ndarray: Steady-state probability distribution
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)
    D2 = np.asarray(D2)
    try:
        result = Queue.steadyState(
            native_to_array(D0), native_to_array(D1),
            native_to_array(D2)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"RAP/RAP/1 steady-state computation failed: {str(e)}")


# ========== SYLVESTER SOLVER ==========

def lib_qmam_sylvester_solve(A, B, C):
    """Solve Sylvester equation AX + XB = C.

    This solver is used internally for matrix equations in queueing analysis.

    Args:
        A (array-like): Left-hand matrix
        B (array-like): Right-hand matrix
        C (array-like): Right-hand side matrix

    Returns:
        ndarray: Solution matrix X
    """
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    try:
        result = Solver.solve(
            native_to_array(A), native_to_array(B),
            native_to_array(C)
        )
        return native_to_array(result)
    except Exception as e:
        raise RuntimeError(f"Sylvester equation solving failed: {str(e)}")
