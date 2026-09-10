"""
SOLVER_CUSTOM  The bare solution algorithm of the custom-solver template.

[Q,U,R,T,C,X] = solver_custom(sn, options)

This is the file to fill in: it receives the NetworkStruct and the options and
must return the six average-metric matrices. The template returns zeros and says
so, exactly as the MATLAB twin does.
"""

import numpy as np


def solver_custom(sn, options=None):
    M = sn.nstations                                   # number of stations
    K = sn.nclasses                                    # number of classes
    N = np.asarray(sn.njobs, dtype=float).ravel()      # job populations
    rates = sn.rates                                   # arrival and service rates
    V = sn.visits                                      # visits

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    print('The solution algorithm needs to be implemented in solver_custom.py: '
          'returning with no result.')
    return QN, UN, RN, TN, CN, XN
