"""The Lam scaling of pfqn_ca must be read off a LOWER bound on log G that a
mixed state can reach, not off one station holding every class at once.

An LQN layer reaches SolverNC as a think-time class beside a zero-Z call class,
and the per-configuration estimate discarded the delay entirely as soon as one
class had no think time. On L=[1e-9,1], N=[99,1], Z=[1,0] -- the layer
lqn_twotasks hands the NC solver -- it returned the all-at-the-queue -2051.6
against a true log G of -359.134, so kscale was -30 and the scaling went the
WRONG WAY: Z/2^-30 = 1.07e9 made the delay column Z^n/n! peak at e^1699 and
overflow at n=[40,0], long before the full population. Python raised
OverflowError there; MATLAB, the JAR and C++ returned G = Inf, lG = Inf.
"""
import numpy as np
import pytest

from line_solver.api.pfqn.nc import pfqn_ca

# 120-digit mpmath convolution of the same recursion
LG_EXACT = -359.13420517157538927


def test_a_zero_think_time_class_does_not_discard_the_delay():
    L = np.array([[1e-9, 1.0]])
    N = np.array([99.0, 1.0])
    Z = np.array([1.0, 0.0])
    G, lG = pfqn_ca(L, N, Z)
    assert lG == pytest.approx(LG_EXACT, rel=1e-12)
    assert G == pytest.approx(np.exp(LG_EXACT), rel=1e-9)


@pytest.mark.parametrize("seed", range(40))
def test_the_estimate_stays_a_lower_bound_on_log_g(seed):
    # A scaling read off an estimate ABOVE log G would drive the recursion to
    # underflow, the other half of the range problem Reiser reports.
    rng = np.random.default_rng(seed)
    M, R = int(rng.integers(1, 4)), int(rng.integers(1, 3))
    L = np.exp(rng.uniform(-12.0, 2.0, size=(M, R)))
    N = rng.integers(1, 7, size=R).astype(float)
    Z = np.exp(rng.uniform(-4.0, 3.0, size=R))
    Z[rng.random(R) < 0.4] = 0.0
    _, lG = pfqn_ca(L, N, Z)
    best = 0.0
    for r in range(R):
        opts = [N[r] * np.log(L[i, r]) for i in range(M) if L[i, r] > 0]
        if Z[r] > 0:
            from math import lgamma
            opts.append(N[r] * np.log(Z[r]) - lgamma(N[r] + 1.0))
        best += max(opts)
    assert best <= lG + 1e-9


def test_the_delay_only_class_is_still_reached_through_a_queue():
    # Every class servable at the delay alone: the estimate must not fall back
    # to the queue column, which is what the old form did whenever one class
    # could not sit at the same single station as the rest.
    L = np.array([[1e-8, 1e-8]])
    N = np.array([30.0, 30.0])
    Z = np.array([5.0, 5.0])
    _, lG = pfqn_ca(L, N, Z)
    assert np.isfinite(lG)
    # 60 jobs at the delay alone already dominates the all-at-the-queue term
    from math import lgamma
    delay_only = sum(30.0 * np.log(5.0) - lgamma(31.0) for _ in range(2))
    assert lG >= delay_only - 1e-9
