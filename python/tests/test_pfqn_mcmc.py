"""Chen-O'Cinneide regularization (pfqn_mcmc).

W. Chen, C. A. O'Cinneide, "Towards a Polynomial-Time Randomized Algorithm for Closed
Product-Form Networks", ACM TOMACS 8(3):227-253, 1998.

Two kinds of assertion are made here, and only the first kind is statistical.

The throughput checks compare the estimator against the EXACT ratio G(N-e_r)/G(N) from
convolution -- pfqn_ca for the single-server models, pfqn_ncld with mu_i(k)=min(k,c_i) for
the multiserver one -- at a tolerance that the measured error clears with room to spare.
They are seeded, so they are deterministic runs of a random algorithm rather than flaky
tests, and the tolerance is set from the estimator's own two-sigma interval rather than
from the observed error, so a real regression still trips them.

The queue-length check is NOT statistical and holds to machine precision: every state of
the regularized chain satisfies sum_i Y(i,r) = N(r), and the reported Q is a weighted
average of those states, so with no delay station the columns of Q sum to N exactly
whatever the sample path was. That invariant is what catches an indexing or weighting
error, which a tolerance on a noisy mean cannot.
"""

import numpy as np
import pytest

from line_solver.api.pfqn.mcmc import pfqn_mcmc
from line_solver.api.pfqn.nc import pfqn_ca, pfqn_nc
from line_solver.api.pfqn.ncld import pfqn_ncld


def _exact_ratios(L, N, Z):
    """X(r) = G(N-e_r)/G(N) by convolution, the quantity pfqn_mcmc estimates."""
    N = np.asarray(N, dtype=float)
    lG = pfqn_ca(L, N, Z)[1]
    X = np.zeros(N.size)
    for r in range(N.size):
        Nr = N.copy()
        Nr[r] -= 1
        X[r] = np.exp(pfqn_ca(L, Nr, Z)[1] - lG)
    return X


def _example_5_1():
    """The 3 single-server stations plus IS station of Example 5.1 of the paper."""
    mu = np.array([0.2, 0.5, 0.8])
    sets = [[0, 1, 2], [0, 1], [0, 2], [1, 2]]
    L = np.zeros((3, 4))
    Z = np.zeros(4)
    for c, st in enumerate(sets):
        L[st, c] = 1.0 / mu[st]
        if c > 0:
            Z[c] = 1.0 / 0.5
    return L, 3 * np.ones(4), Z


def test_example_5_1_matches_the_exact_ratios():
    L, N, Z = _example_5_1()
    res = pfqn_mcmc(L, N, Z, None, {'samples': 200000, 'seed': 23000})
    exact = _exact_ratios(L, N, Z)
    assert np.allclose(res.X, exact, rtol=0.02)
    # Every interval is a real interval around the estimate, and every standard error is
    # a positive number rather than a NaN from a zero denominator.
    assert np.all(res.Xse > 0)
    assert np.all(res.Xlo < res.X) and np.all(res.X < res.Xhi)
    assert (res.batches, res.burnin) == (30, res.samples // 10)


def test_the_estimator_is_consistent():
    """Tenfold more work must buy roughly sqrt(10) less error, not a fixed floor.

    A biased estimator -- a mis-weighted holding time, a warm-up that never ends -- passes
    the tolerance check above at one sample size and then stops improving. Comparing two
    sizes separates noise from bias without pinning either run's value.
    """
    L, N, Z = _example_5_1()
    exact = _exact_ratios(L, N, Z)
    err = []
    for samples in (50000, 500000):
        res = pfqn_mcmc(L, N, Z, None, {'samples': samples, 'seed': 23000})
        err.append(np.max(np.abs(res.X - exact) / exact))
    assert err[1] < 0.5 * err[0]


def test_queue_lengths_conserve_the_population_exactly():
    """sum_i Y(i,r) = N(r) in every state, so the weighted average inherits it."""
    L = np.array([[0.6, 0.2], [0.3, 0.5], [0.1, 0.4]])
    N = np.array([4, 3])
    res = pfqn_mcmc(L, N, None, None, {'samples': 50000, 'seed': 23000})
    assert np.allclose(res.Q.sum(axis=0), N, atol=1e-12)


def test_multiserver_matches_the_exact_load_dependent_constant():
    """The paper's own selling point (its Tables IV and V): exact multiservers.

    The reference here is the load-dependent convolution with mu_i(k) = min(k, c_i), which
    is the same product form the regularized chain samples, so any disagreement beyond the
    Monte Carlo error is an error in the Psi_i(Y_i) = min(s_i, Y_i) rate.
    """
    L = np.array([[0.6, 0.2], [0.3, 0.5], [0.1, 0.4]])
    N = np.array([4, 3])
    c = np.array([2.0, 1.0, 3.0])
    Ntot = int(N.sum())
    mu = np.array([[min(k, c[i]) for k in range(1, Ntot + 1)] for i in range(3)])
    Z = np.zeros(2)
    lG = pfqn_ncld(L, N, Z, mu, {'method': 'exact'}).lG
    exact = np.zeros(2)
    for r in range(2):
        Nr = N.astype(float).copy()
        Nr[r] -= 1
        exact[r] = np.exp(pfqn_ncld(L, Nr, Z, mu, {'method': 'exact'}).lG - lG)
    res = pfqn_mcmc(L, N, None, c, {'samples': 400000, 'seed': 23000})
    assert np.allclose(res.X, exact, rtol=0.02)
    assert np.allclose(res.Q.sum(axis=0), N, atol=1e-12)


def test_a_fractional_population_is_refused_not_rounded():
    """The chain lives on the integer lattice, so a fractional N has no state space.

    This is what makes mcmc unusable as a SolverLN layer method, which is correct rather
    than unfortunate: rounding would answer a question nobody asked.
    """
    L = np.array([[0.6, 0.2], [0.3, 0.5]])
    with pytest.raises(ValueError, match='fractional'):
        pfqn_mcmc(L, np.array([2.5, 1.0]), None, None, {'samples': 1000, 'seed': 23000})


def test_an_infinite_population_is_refused():
    L = np.array([[0.6, 0.2], [0.3, 0.5]])
    with pytest.raises(ValueError, match='closed model'):
        pfqn_mcmc(L, np.array([np.inf, 1.0]), None, None, {'samples': 1000, 'seed': 23000})


def test_a_populated_class_with_no_demand_anywhere_is_refused():
    L = np.array([[0.6, 0.0], [0.3, 0.0]])
    with pytest.raises(ValueError, match='no demand'):
        pfqn_mcmc(L, np.array([2, 1]), None, None, {'samples': 1000, 'seed': 23000})


def test_an_empty_network_returns_zeros_without_simulating():
    L = np.array([[0.6, 0.2], [0.3, 0.5]])
    res = pfqn_mcmc(L, np.zeros(2), None, None, {'samples': 1000, 'seed': 23000})
    assert np.all(res.X == 0) and np.all(res.Q == 0)
    assert (res.batches, res.samples, res.burnin) == (0, 0, 0)


def test_pfqn_nc_mcmc_arm_returns_the_ble_constant_not_a_paper_result():
    """The method estimates ratios and never forms G, so pfqn_nc supplies lG from BLE.

    Pinning it here records that the number is deliberate: it cancels out of every mean
    value the analyzer reports, and only getProbNormConstAggr reads it.
    """
    L, N, Z = _example_5_1()
    from line_solver.api.pfqn.asymptotic import pfqn_ble
    assert pfqn_nc(L, N, Z, method='mcmc')[1] == pytest.approx(pfqn_ble(L, N, Z)[1])
