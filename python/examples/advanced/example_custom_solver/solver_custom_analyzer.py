"""
SOLVER_CUSTOM_ANALYZER  Everything around the solution algorithm.

[Q,U,R,T,C,X,runtime] = solver_custom_analyzer(sn, options)

Any activity prior to or after launching the solution algorithm belongs here:
the template times the call and scrubs the NaNs the algorithm may leave behind.
"""

import time

import numpy as np

from solver_custom import solver_custom


def solver_custom_analyzer(sn, options=None):
    t0 = time.time()

    print('Any activity prior or after launching the solution algorithm needs '
          'to be implemented in solver_custom_analyzer.py.')
    QN, UN, RN, TN, CN, XN = solver_custom(sn, options)

    out = []
    for m in (QN, UN, RN, TN, CN, XN):
        m = np.asarray(m, dtype=float)
        m[np.isnan(m)] = 0.0
        out.append(m)

    return tuple(out) + (time.time() - t0,)
