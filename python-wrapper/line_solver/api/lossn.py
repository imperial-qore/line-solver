
"""
Loss network analysis functions.

This module provides functions for analyzing loss networks, where
customers are blocked and lost when all servers are busy. The primary
function implements the Erlang fixed-point algorithm for multi-service
loss networks.

Loss networks are used to model circuit-switched networks, call centers
with blocking, and other systems where customers are rejected when
resources are unavailable.
"""

import numpy as np
import jpype
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def lossn_erlangfp(nu_vec, amat, c_vec):
    """
    Erlang fixed-point algorithm for multi-service loss networks.
    
    Computes blocking probabilities and performance measures for
    loss networks where blocked customers are lost (not queued).
    
    Args:
        nu_vec: Vector of traffic intensities.
        amat: Service requirement matrix.
        c_vec: Vector of link capacities.
        
    Returns:
        tuple: (qlen, loss_prob, block_prob, niter) containing:
            - qlen: Mean queue lengths
            - loss_prob: Loss probabilities
            - block_prob: Blocking probabilities
            - niter: Number of iterations to convergence
    """
    nu_matrix = jlineMatrixFromArray(nu_vec)
    a_matrix = jlineMatrixFromArray(amat)
    c_matrix = jlineMatrixFromArray(c_vec)

    result = jpype.JPackage('jline').api.lossn.Lossn_erlangfpKt.lossn_erlangfp(
        nu_matrix, a_matrix, c_matrix
    )

    qlen = jlineMatrixToArray(result.qLen)
    loss_prob = jlineMatrixToArray(result.lossProb)
    block_prob = jlineMatrixToArray(result.blockProb)
    niter = result.niter

    return qlen, loss_prob, block_prob, niter
