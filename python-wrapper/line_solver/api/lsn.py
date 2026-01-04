
"""
Layered stochastic network (LSN) analysis functions.

This module provides utility functions for analyzing layered stochastic
networks, which model systems with hierarchical service dependencies
and multiplicities.

Layered stochastic networks extend basic queueing networks to capture
complex service dependencies and resource constraints in distributed
systems and software architectures.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def lsn_max_multiplicity(lsn):
    """
    Compute maximum multiplicity constraints for layered stochastic network.
    
    Determines the maximum allowable multiplicities for each layer
    in the stochastic network based on capacity and dependency constraints.
    
    Args:
        lsn: Layered stochastic network object.
        
    Returns:
        numpy.ndarray: Maximum multiplicity values for each layer.
    """
    result = jpype.JPackage('jline').api.lsn.LsnMaxMultiplicityKt.lsn_max_multiplicity(lsn)
    return jlineMatrixToArray(result)