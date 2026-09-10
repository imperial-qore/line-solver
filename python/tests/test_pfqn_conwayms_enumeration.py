"""pfqn_conwayms must terminate, and its state sweep must visit each state once.

`pfqn_conwayms` sums over every way of distributing the busy servers of a station
among the R classes -- all length-R nonnegative vectors summing to `nservers` --
via the `sprod`/`sprod_next` iterator that the caller drives as
`while s >= 0: ...; s, nvec = sprod_next(s, SD, D)`. The returned index IS the
iterator, and `sprod_next` returned the literal 0 on every successful step, so `s`
oscillated 0 -> 1 -> 0 and the loop re-decoded state 1 forever. The routine never
returned, at any size: M=1, R=2, maxiter=5 was enough to hang it.

Fixing the index alone was not enough. `sprod` seeded the sweep at (n,0,...,0)
while the unranking puts (0,...,0,n) at rank 0, so the first state was visited
twice and one state never -- a wrong SUM rather than a hang, which is the harder
failure to notice. Both ends now go through `_sprod_unrank`.
"""
import itertools

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_conwayms
from line_solver.api.pfqn.linearizerms import sprod, sprod_next

GUARD = 200000   # any sweep longer than this is the non-termination bug, not a slow test


def _sweep(R, n):
    """Every state the iterator yields, in order, with a hard non-termination guard."""
    seen = []
    s, nvec, SD, D = sprod(R, n)
    while s >= 0:
        seen.append(tuple(int(x) for x in nvec))
        s, nvec = sprod_next(s, SD, D)
        assert len(seen) <= GUARD, 'sprod_next never returned -1 for R=%d, n=%d' % (R, n)
    return seen


@pytest.mark.parametrize('R', [1, 2, 3, 4])
@pytest.mark.parametrize('n', [0, 1, 2, 3, 5])
def test_sweep_visits_every_composition_exactly_once(R, n):
    seen = _sweep(R, n)
    expected = sorted(v for v in itertools.product(range(n + 1), repeat=R) if sum(v) == n)
    assert len(seen) == len(set(seen)), 'a state was visited twice: %r' % (seen,)
    assert sorted(seen) == expected


@pytest.mark.parametrize('R,n', [(2, 3), (3, 4)])
def test_sweep_length_is_the_stars_and_bars_count(R, n):
    from math import comb
    assert len(_sweep(R, n)) == comb(n + R - 1, R - 1)


def test_conwayms_terminates_and_returns_finite_metrics():
    D = np.array([[0.6, 0.4], [0.9, 0.3], [0.5, 0.8]])
    Q, U, R, C, X, totiter = pfqn_conwayms(
        D, np.array([3.0, 2.0]), np.array([1.0, 1.0]),
        np.array([1.0, 2.0, 3.0]), None, 1e-8, 1000)
    for name, val in (('Q', Q), ('U', U), ('R', R), ('X', X)):
        assert np.all(np.isfinite(np.asarray(val, dtype=float))), '%s is not finite' % name
    assert np.all(np.asarray(X, dtype=float) > 0)
    assert 0 < totiter <= 4 * 1000


def test_conwayms_single_server_reproduces_matlab():
    """With every station single-server the multiserver arm is bypassed, so this
    pins the shared Core/Estimate/ForwardMVA port against MATLAB pfqn_conwayms."""
    D = np.array([[0.6, 0.4], [0.9, 0.3], [0.5, 0.8]])
    _, _, _, _, X, _ = pfqn_conwayms(
        D, np.array([3.0, 2.0]), np.array([1.0, 1.0]),
        np.array([1.0, 1.0, 1.0]), None, 1e-8, 1000)
    matlab = np.array([0.616471377, 0.439592501])
    assert np.allclose(np.asarray(X, dtype=float).ravel(), matlab, atol=1e-9)


# (L, N, Z, nservers, MATLAB X). The multiserver arm reaches XR/XE, whose
# normalizing sums are EMPTY whenever the reduced population cannot fill the
# servers (row 2) or leaves a class with no job (row 5); MATLAB divided by that
# zero, and the resulting NaN made Core run to maxiter and starve the remaining
# population sweeps of iterations. See _kb/06-solver-catalog.md.
CONWAY_MS_CASES = [
    ([[0.5], [1.0]], [4.0], [1.0], [1.0, 3.0], [1.325424111]),
    ([[0.5], [1.0]], [2.0], [1.0], [1.0, 3.0], [0.763932023]),
    ([[0.6, 0.4], [0.9, 0.3]], [3.0, 2.0], [1.0, 1.0], [1.0, 2.0],
     [0.870071451, 0.711994388]),
    ([[0.6, 0.4], [0.9, 0.3]], [5.0, 4.0], [0.0, 0.0], [2.0, 3.0],
     [1.758200984, 1.931225373]),
    ([[0.6, 0.4, 0.2], [0.9, 0.3, 0.7], [0.2, 0.5, 0.4]], [3.0, 2.0, 2.0], [1.0, 1.0, 1.0],
     [1.0, 2.0, 3.0], [0.781242573, 0.568153763, 0.544295201]),
    ([[0.6, 0.4], [0.9, 0.3]], [1.0, 1.0], [1.0, 1.0], [1.0, 4.0],
     [0.381505806, 0.536169058]),
]


@pytest.mark.parametrize('L,N,Z,nservers,matlabX', CONWAY_MS_CASES)
def test_conwayms_multiserver_reproduces_matlab(L, N, Z, nservers, matlabX):
    L, N, Z = np.array(L), np.array(N), np.array(Z)
    Q, _, _, _, X, _ = pfqn_conwayms(L, N, Z, np.array(nservers), None, 1e-8, 1000)
    X = np.asarray(X, dtype=float).ravel()
    assert np.all(np.isfinite(X)), 'X is not finite: %r' % (X,)
    assert np.allclose(X, np.array(matlabX), atol=1e-9)
    # jobs are conserved: queued + thinking = N, per class
    assert np.allclose(np.asarray(Q, dtype=float).sum(axis=0) + X * Z, N, atol=1e-6)
