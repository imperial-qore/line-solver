
"""
Discrete-Time Markov Chain (DTMC) analysis algorithms.

This module provides functions for analyzing discrete-time Markov chains,
including steady-state analysis, stochastic complement operations,
time-reversal, and simulation methods.

These functions are used by various solvers and analysis methods.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def dtmc_solve(matrix):
    """
    Solve for steady-state probabilities of a DTMC.
    
    Args:
        matrix: Transition probability matrix.
        
    Returns:
        numpy.ndarray: Steady-state probability distribution.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Dtmc_solveKt.dtmc_solve(
            jlineMatrixFromArray(matrix)
        )
    )


def dtmc_stochcomp(matrix, indexes):
    """
    Compute stochastic complement of a DTMC.
    
    Performs state space reduction by computing the stochastic complement,
    eliminating specified states while preserving transition probabilities.
    
    Args:
        matrix: Transition probability matrix.
        indexes: List of state indexes to eliminate.
        
    Returns:
        numpy.ndarray: Reduced transition matrix (stochastic complement).
    """
    ind = jpype.java.util.ArrayList()
    for i in range(len(indexes)):
        ind.add(jpype.JInt(indexes[i]))

    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Dtmc_stochcompKt.dtmc_stochcomp(
            jlineMatrixFromArray(matrix), ind
        )
    )


def dtmc_timereverse(matrix):
    """
    Compute time-reversed DTMC.
    
    Constructs the time-reversed discrete-time Markov chain using
    the detailed balance equations and steady-state probabilities.
    
    Args:
        matrix: Original transition probability matrix.
        
    Returns:
        numpy.ndarray: Time-reversed transition probability matrix.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Dtmc_timereverseKt.dtmc_timereverse(
            jlineMatrixFromArray(matrix)
        )
    )


def dtmc_makestochastic(matrix):
    """
    Normalize matrix to be row-stochastic.
    
    Converts a non-negative matrix to a proper stochastic matrix
    by normalizing rows to sum to 1.
    
    Args:
        matrix: Non-negative matrix to normalize.
        
    Returns:
        numpy.ndarray: Row-stochastic matrix.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Dtmc_makestochasticKt.dtmc_makestochastic(
            jlineMatrixFromArray(matrix)
        )
    )


def dtmc_rand(length):
    """
    Generate random DTMC transition matrix.
    
    Creates a random row-stochastic matrix suitable for use as
    a discrete-time Markov chain transition matrix.
    
    Args:
        length: Size of the square matrix (number of states).
        
    Returns:
        numpy.ndarray: Random stochastic transition matrix.
    """
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.mc.Dtmc_randKt.dtmc_rand(jpype.JInt(length))
    )


def dtmc_simulate(P, pi0, n):
    """
    Simulate DTMC state trajectory.
    
    Generates a sample path of the discrete-time Markov chain
    starting from an initial distribution.
    
    Args:
        P: Transition probability matrix.
        pi0: Initial state distribution.
        n: Number of time steps to simulate.
        
    Returns:
        numpy.ndarray: Sequence of visited states.
    """
    result = jpype.JPackage('jline').api.mc.Dtmc_simulateKt.dtmc_simulate(
        jlineMatrixFromArray(P),
        jlineMatrixFromArray(pi0),
        jpype.JInt(n)
    )
    return np.array([result[i] for i in range(len(result))])


def dtmc_isfeasible(P):
    """
    Check if matrix is a valid DTMC transition matrix.
    
    Verifies that the matrix satisfies the requirements for a
    discrete-time Markov chain: non-negative entries and row sums equal to 1.
    
    Args:
        P: Matrix to check.
        
    Returns:
        bool: True if matrix is a valid DTMC transition matrix.
    """
    return jpype.JPackage('jline').api.mc.Dtmc_isfeasibleKt.dtmc_isfeasible(
        jlineMatrixFromArray(P)
    )