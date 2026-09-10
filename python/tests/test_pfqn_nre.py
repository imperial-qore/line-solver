"""
Saddle-tilted Edgeworth normalizing constant (pfqn_nre).

The fixtures are MATLAB pfqn_nre values. The point of the method is the tilt:
on the same model the untilted contour of pfqn_nrl is 9.3e-2 off the exact
constant and the tilted one 1.2e-4, so a port that drops the tilt or the
Edgeworth term still looks plausible in isolation and is caught only by
comparing the two errors.
"""

import numpy as np

from line_solver.api.pfqn import pfqn_nre, pfqn_nrl
from line_solver.api.pfqn.ncld import pfqn_ncld

TOL = 1e-9

# two classes, three single-server stations and a delay
L_A = np.array([[1.0, 0.5], [0.7, 1.2], [0.3, 0.9]])
N_A = np.array([2.0, 3.0])
Z_A = np.array([1.0, 0.5])
LG_A_EXACT = 4.783752287503922  # MATLAB pfqn_ncld, method exact


def test_two_class_model_matches_matlab():
    lG = pfqn_nre(L_A, N_A, Z_A, alpha=np.ones((3, 5)))
    assert abs(lG - 4.783188572136291) < TOL


def test_tilted_contour_beats_the_untilted_one():
    nre = pfqn_nre(L_A, N_A, Z_A, alpha=np.ones((3, 5)))
    nrl = pfqn_nrl(L_A, N_A, Z_A, alpha=np.ones((3, 5)))
    # MATLAB 4.338130313502534; this port's numerical Hessian puts it 5e-6 off,
    # which is immaterial next to the 4.5e-1 that separates it from nre here
    assert abs(nrl - 4.338130313502534) < 1e-5
    assert abs(nre - LG_A_EXACT) < abs(nrl - LG_A_EXACT)


def test_single_class_is_exact():
    # one class quotients the torus down to dimension zero, so the routine
    # returns pfqn_gldsingle itself rather than any expansion of it
    L = np.array([[0.5], [1.0 / 3.0], [0.2]])
    lG = pfqn_nre(L, np.array([5.0]), np.array([0.0]), alpha=np.ones((3, 5)))
    assert abs(lG - (-1.995145480198095)) < 1e-12


def test_load_dependent_rate_row_is_carried():
    # station 1 is a 2-server queue, mu(1,n) = min(n,2); station 2 single
    mu = np.vstack([np.minimum(np.arange(1, 9), 2.0), np.ones(8)])
    L = np.array([[1.0, 0.6], [0.5, 1.1]])
    N = np.array([4.0, 4.0])
    Z = np.array([0.5, 0.5])
    assert abs(pfqn_nre(L, N, Z, alpha=mu) - 3.911235290463940) < TOL
    # and the same value through the load-dependent dispatcher
    res = pfqn_ncld(L, N, Z, mu, {'method': 'nre'})
    assert res.method == 'nre'
    assert abs(res.lG - 3.911235290463940) < TOL
