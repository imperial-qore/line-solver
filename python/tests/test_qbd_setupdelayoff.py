"""
Regression test: qbd_setupdelayoff level summation.

The Python port replicated an old MATLAB defect that summed n+1 phase
probabilities per level while advancing the cursor by n, so one phase per level
was counted twice under two different weights and the terminating check dropped
the last level. The queue length came out high -- by 1.8% for a fast setup and
25% for a slow one, since the overlapped phases hold most mass when the setup is
slow. References are the MATLAB values, which agree with a 4e6-sample JMT
simulation of the same M/M/1 with setup (1.2290 and 3.4162).
"""
import numpy as np
import pytest

from line_solver.api.mam.qbd import qbd_setupdelayoff

LAMBDA, MU, BETA = 0.5, 1.0, 4.0


@pytest.mark.parametrize('alpharate,alphascv,expected', [
    (2.0, 1.0, 1.22727272),   # exponential setup, mean 0.5
    (0.2, 1.0, 3.41379310),   # exponential setup, mean 5.0 (slow: worst case)
    (2.0, 0.25, 1.21022727),  # Erlang-4 setup, mean 0.5
])
def test_matches_matlab(alpharate, alphascv, expected):
    qn = qbd_setupdelayoff(LAMBDA, MU, alpharate, alphascv, BETA, 1.0)
    assert np.abs(qn - expected) / expected < 1e-6


def test_monotone_in_setup_time():
    # A slower setup can only lengthen the queue.
    means = [0.1, 0.5, 1.0, 2.0, 5.0]
    qns = [qbd_setupdelayoff(LAMBDA, MU, 1.0 / m, 1.0, BETA, 1.0) for m in means]
    assert all(b > a for a, b in zip(qns, qns[1:]))


def test_setup_exceeds_setup_free_mm1():
    # With rho = 0.5 the setup-free M/M/1 holds rho/(1-rho) = 1 job; any setup
    # adds to that.
    assert qbd_setupdelayoff(LAMBDA, MU, 2.0, 1.0, BETA, 1.0) > 1.0
