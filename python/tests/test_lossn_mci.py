"""Validate lossn_mci (Monte Carlo summation) against exact ground truth."""
import itertools
import numpy as np
from scipy import special

from line_solver.api.lossn import lossn_mci, lossn_erlangfp


def erlangB_ref(rho, C):
    inv = 1.0
    for k in range(1, C + 1):
        inv = 1.0 + inv * k / rho
    return 1.0 / inv


def lossn_exact(nu, A, C):
    """Exact normalization constant and class blocking by box enumeration."""
    nu = np.asarray(nu, float).ravel()
    A = np.asarray(A, float)
    C = np.asarray(C, float).ravel()
    R = len(nu)
    N = [int(np.floor(np.min(C[A[:, k] > 0] / A[A[:, k] > 0, k]))) for k in range(R)]
    states = np.array(list(itertools.product(*[range(n + 1) for n in N])), dtype=float)
    logq = states @ np.log(nu) - np.sum(special.gammaln(states + 1), axis=1)
    AV = states @ A.T
    feas = np.all(AV <= C[None, :], axis=1)
    mq = logq[feas].max()
    g = np.exp(mq) * np.sum(np.exp(logq[feas] - mq))
    beta = np.zeros(R)
    for r in range(R):
        Cr = (C - A[:, r])[None, :]
        feasR = np.all(AV <= Cr, axis=1)
        gr = np.exp(mq) * np.sum(np.exp(logq[feasR] - mq))
        beta[r] = 1.0 - gr / g
    return g, beta


def test_single_link_erlangB():
    rho, Ccap = 8.0, 10
    _, loss, lG, ci, _ = lossn_mci(rho, np.array([[1.0]]), np.array([Ccap]),
                                   samples=200000, seed=1)
    exactB = erlangB_ref(rho, Ccap)
    g, _ = lossn_exact([rho], [[1.0]], [Ccap])
    assert ci['loss'][0, 0] <= exactB <= ci['loss'][0, 1], (exactB, ci['loss'][0])
    assert abs(lG - np.log(g)) < 0.05


def test_two_link_multirate_exact():
    nu = np.array([3.0, 1.5])
    A = np.array([[1.0, 1.0], [1.0, 2.0]])
    C = np.array([4.0, 5.0])
    g, beta = lossn_exact(nu, A, C)
    qlen, loss, lG, ci, _ = lossn_mci(nu, A, C, samples=300000, seed=2)
    assert abs(lG - np.log(g)) < 0.05
    for r in range(2):
        assert ci['loss'][r, 0] <= beta[r] <= ci['loss'][r, 1], (r, beta[r], ci['loss'][r])
    assert np.max(np.abs(qlen - nu * (1 - loss))) < 1e-9


def test_reproducible_seed():
    nu = np.array([3.0, 1.5])
    A = np.array([[1.0, 1.0], [1.0, 2.0]])
    C = np.array([4.0, 5.0])
    _, l1, _, _, _ = lossn_mci(nu, A, C, samples=50000, seed=7)
    _, l2, _, _, _ = lossn_mci(nu, A, C, samples=50000, seed=7)
    assert np.array_equal(l1, l2)


def star_network(regime):
    C = np.array([90.0, 100.0, 110.0, 120.0])
    routes = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    base, hi = {'light': (9.0, 1.6), 'moderate': (10.0, 2.0),
                'heavy': (15.0, 3.0)}[regime]
    A = np.zeros((4, 12))
    nu = np.zeros(12)
    for p, (a, b) in enumerate(routes):
        A[a, p] = 1.0
        A[b, p] = 1.0
        nu[p] = base
        A[a, 6 + p] = 5.0
        A[b, 6 + p] = 5.0
        nu[6 + p] = hi
    return nu, A, C


def test_star_vs_erlangfp_acceptance():
    nu, A, C = star_network('light')
    _, lossFP, _, _ = lossn_erlangfp(nu, A, C)  # blocking probability
    accFP = 1.0 - lossFP
    _, loss, lG, ci, _ = lossn_mci(nu, A, C, samples=200000, seed=3)
    accMC = 1.0 - loss
    assert np.all((loss >= 0) & (loss <= 1))
    assert np.max(np.abs(accMC - accFP)) < 0.02


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_') and callable(fn):
            fn()
            print(f'PASS: {name}')
    print('all lossn_mci tests passed')
