"""
Analytic validation of the discrete-time qsys formulas (Geo/Geo/1, Geo^X/Geo/1).

Mirrors the JAR suites jline.api.qsys.QsysGeoGeo1Test and QsysGeoXGeo1Test.
Every assertion is an exact identity rather than a stored number:
  - the stationary pmf normalizes and reproduces the closed-form mean;
  - Little's law and the in-service decomposition hold in both conventions;
  - the conventions differ by exactly one slot of sojourn and agree on the
    waiting time;
  - the closed forms agree with a brute-force power iteration of the
    underlying chains;
  - the batch result collapses to the single-arrival result at beta = 1.
"""

import numpy as np
import pytest

from line_solver.api.dqsys import dqsys_geogeo1, dqsys_geoxgeo1, dqsys_geoxgeo1_moments

TOL = 1e-9
NUM_TOL = 1e-7

CONVENTIONS = ('LAS_DA', 'EAS')

# (a, s) with a < s so the queue is stable.
GRID = [(0.1, 0.9), (0.2, 0.5), (0.3, 0.4), (0.05, 0.95),
        (0.4, 0.8), (0.45, 0.5), (0.25, 1.0), (0.6, 0.75)]

# (a, beta, s) with a/beta < s.
BATCH_GRID = [(0.1, 0.5, 0.9), (0.2, 0.8, 0.6), (0.15, 0.4, 0.8),
              (0.3, 1.0, 0.7), (0.05, 0.2, 0.5), (0.4, 0.9, 0.9)]


# ----------------------------------------------------------------------
# Geo/Geo/1 distributional identities
# ----------------------------------------------------------------------

@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,s", GRID)
def test_pmf_normalizes(a, s, convention):
    r = dqsys_geogeo1(a, s, convention)
    mass = sum(r['pmf'](n) for n in range(0, 20001))
    assert abs(mass - 1.0) < NUM_TOL


@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,s", GRID)
def test_pmf_mean_matches_closed_form(a, s, convention):
    r = dqsys_geogeo1(a, s, convention)
    mean = sum(n * r['pmf'](n) for n in range(1, 20001))
    assert abs(mean - r['meanQueueLength']) < NUM_TOL


@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,s", GRID)
def test_empty_prob_matches_pmf_at_zero(a, s, convention):
    r = dqsys_geogeo1(a, s, convention)
    assert abs(r['emptyProb'] - r['pmf'](0)) < TOL


# ----------------------------------------------------------------------
# Operational laws
# ----------------------------------------------------------------------

@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,s", GRID)
def test_littles_law(a, s, convention):
    r = dqsys_geogeo1(a, s, convention)
    assert abs(r['meanQueueLength'] - r['throughput'] * r['meanSojournTime']) < TOL
    assert abs(r['meanWaitingQueue'] - r['throughput'] * r['meanWaitingTime']) < TOL


@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,s", GRID)
def test_sojourn_splits_into_waiting_plus_service(a, s, convention):
    r = dqsys_geogeo1(a, s, convention)
    assert abs(r['meanSojournTime']
               - (r['meanWaitingTime'] + r['meanServiceTime'])) < TOL
    assert abs((r['meanQueueLength'] - r['meanWaitingQueue'])
               - r['throughput'] * r['meanServiceTime']) < TOL


@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,s", GRID)
def test_utilization_is_load(a, s, convention):
    r = dqsys_geogeo1(a, s, convention)
    assert abs(r['utilization'] - a / s) < TOL


# ----------------------------------------------------------------------
# Relation between the two conventions
# ----------------------------------------------------------------------

@pytest.mark.parametrize("a,s", GRID)
def test_conventions_differ_by_exactly_one_slot(a, s):
    las = dqsys_geogeo1(a, s, 'LAS_DA')
    eas = dqsys_geogeo1(a, s, 'EAS')
    assert abs((las['meanSojournTime'] - eas['meanSojournTime']) - 1.0) < TOL
    assert abs((las['meanQueueLength'] - eas['meanQueueLength']) - a) < TOL
    assert abs(las['meanWaitingTime'] - eas['meanWaitingTime']) < TOL


def test_las_da_is_the_default():
    default = dqsys_geogeo1(0.2, 0.5)
    assert default['convention'] == 'LAS_DA'
    assert abs(default['meanSojournTime']
               - dqsys_geogeo1(0.2, 0.5, 'LAS_DA')['meanSojournTime']) < TOL


# ----------------------------------------------------------------------
# Degenerate limits
# ----------------------------------------------------------------------

def test_vanishing_load_leaves_only_the_service_time():
    s = 0.4
    las = dqsys_geogeo1(1e-9, s, 'LAS_DA')
    eas = dqsys_geogeo1(1e-9, s, 'EAS')
    assert abs(las['meanSojournTime'] - 1.0 / s) < 1e-6
    assert abs(eas['meanSojournTime'] - (1.0 - s) / s) < 1e-6
    assert abs(las['meanWaitingTime']) < 1e-6


@pytest.mark.parametrize("a", [0.05, 0.3, 0.7, 0.99])
def test_deterministic_unit_service_never_queues(a):
    r = dqsys_geogeo1(a, 1.0, 'LAS_DA')
    assert abs(r['meanSojournTime'] - 1.0) < TOL
    assert abs(r['meanWaitingTime']) < TOL
    assert abs(r['ratio']) < TOL


# ----------------------------------------------------------------------
# Brute-force cross-check of the underlying chains
# ----------------------------------------------------------------------

def _solve_birth_death(up0, down, up, n_max=4000):
    """Power-iterate a truncated discrete birth-death chain to stationarity."""
    pi = np.zeros(n_max + 1)
    pi[0] = 1.0
    for _ in range(400000):
        nxt = np.zeros(n_max + 1)
        nxt[0] += pi[0] * (1.0 - up0)
        nxt[1] += pi[0] * up0
        stay = 1.0 - down - up
        nxt[0:n_max] += pi[1:n_max + 1] * down
        nxt[1:n_max] += pi[1:n_max] * stay
        nxt[n_max] += pi[n_max] * (1.0 - down)
        nxt[2:n_max + 1] += pi[1:n_max] * up
        delta = np.abs(nxt - pi).sum()
        pi = nxt
        if delta < 1e-14:
            break
    return pi


def _assert_chain_agrees(pi, r, label):
    assert abs(pi.sum() - 1.0) < 1e-9, "chain mass, " + label
    assert abs(pi[0] - r['emptyProb']) < 1e-6, "pi_0, " + label
    mean = float((np.arange(pi.size) * pi).sum())
    assert abs(mean - r['meanQueueLength']) < 1e-5, "E[N], " + label
    for n in range(0, 11):
        assert abs(pi[n] - r['pmf'](n)) < 1e-6, "pi_%d, %s" % (n, label)


@pytest.mark.parametrize("a,s", GRID)
def test_las_da_matches_brute_force_chain(a, s):
    # Departure resolved before arrival: from n >= 1 the content falls with
    # probability s(1-a) and rises with a(1-s); the empty state rises with
    # probability a because no job can depart.
    pi = _solve_birth_death(a, s * (1.0 - a), a * (1.0 - s))
    _assert_chain_agrees(pi, dqsys_geogeo1(a, s, 'LAS_DA'),
                         "LAS-DA a=%s s=%s" % (a, s))


@pytest.mark.parametrize("a,s", GRID)
def test_eas_matches_brute_force_chain(a, s):
    # A job arriving into an empty system may depart in the same slot, so the
    # empty state rises with probability a(1-s) as well.
    pi = _solve_birth_death(a * (1.0 - s), s * (1.0 - a), a * (1.0 - s))
    _assert_chain_agrees(pi, dqsys_geogeo1(a, s, 'EAS'), "EAS a=%s s=%s" % (a, s))


# ----------------------------------------------------------------------
# Geo^X/Geo/1
# ----------------------------------------------------------------------

@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,s", [(0.2, 0.5), (0.1, 0.4), (0.45, 0.5), (0.25, 1.0)])
def test_degenerate_batch_reproduces_geogeo1(a, s, convention):
    single = dqsys_geogeo1(a, s, convention)
    batch = dqsys_geoxgeo1(a, 1.0, s, convention)
    for key in ('meanQueueLength', 'meanSojournTime', 'meanWaitingTime',
                'utilization', 'throughput'):
        assert abs(single[key] - batch[key]) < TOL, key


@pytest.mark.parametrize("convention", CONVENTIONS)
@pytest.mark.parametrize("a,beta,s", BATCH_GRID)
def test_batch_littles_law(a, beta, s, convention):
    r = dqsys_geoxgeo1(a, beta, s, convention)
    assert abs(r['meanQueueLength'] - r['throughput'] * r['meanSojournTime']) < 1e-9
    assert abs(r['meanWaitingQueue'] - r['throughput'] * r['meanWaitingTime']) < 1e-9
    assert abs(r['meanSojournTime']
               - (r['meanWaitingTime'] + r['meanServiceTime'])) < 1e-9


@pytest.mark.parametrize("a,beta,s", BATCH_GRID)
def test_batch_rates_and_empty_prob(a, beta, s):
    r = dqsys_geoxgeo1(a, beta, s)
    assert abs(r['arrivalRate'] - a / beta) < TOL
    assert abs(r['utilization'] - (a / beta) / s) < TOL
    assert abs(r['boundaryEmptyProb'] - (1.0 - r['utilization'])) < TOL


@pytest.mark.parametrize("a,beta,s", BATCH_GRID)
def test_batch_conventions_differ_by_one_slot(a, beta, s):
    las = dqsys_geoxgeo1(a, beta, s, 'LAS_DA')
    eas = dqsys_geoxgeo1(a, beta, s, 'EAS')
    assert abs((las['meanSojournTime'] - eas['meanSojournTime']) - 1.0) < 1e-9
    assert abs((las['meanQueueLength'] - eas['meanQueueLength'])
               - las['arrivalRate']) < 1e-9
    assert abs(las['meanWaitingTime'] - eas['meanWaitingTime']) < 1e-9


def test_batching_increases_congestion_at_equal_load():
    # Hold lambda fixed and enlarge the batches: the extra within-batch
    # queueing must show up as a longer sojourn.
    s, lambda_val, previous = 0.8, 0.4, 0.0
    for beta in (1.0, 0.8, 0.5, 0.25):
        r = dqsys_geoxgeo1(lambda_val * beta, beta, s)
        assert abs(r['arrivalRate'] - lambda_val) < 1e-12
        assert r['meanSojournTime'] > previous
        previous = r['meanSojournTime']


@pytest.mark.parametrize("a,beta,s", BATCH_GRID)
def test_pgf_normalizes_and_reproduces_the_mean(a, beta, s):
    r = dqsys_geoxgeo1(a, beta, s)
    assert abs(r['pgf'](1.0, 1.0) - 1.0) < TOL

    def slot_pgf(z):
        return 1.0 - a + a * (beta * z / (1.0 - (1.0 - beta) * z))

    # Both numerator and denominator vanish at z=1, so evaluating within
    # rounding distance of 1 loses every significant digit; stay a finite step
    # away and take a second-order backward stencil using P(1)=1 exactly.
    h = 1e-4
    p1 = r['pgf'](1.0 - h, slot_pgf(1.0 - h))
    p2 = r['pgf'](1.0 - 2.0 * h, slot_pgf(1.0 - 2.0 * h))
    slope = (3.0 - 4.0 * p1 + p2) / (2.0 * h)
    assert abs(slope - r['meanQueueLength']) < 1e-4 * max(1.0, r['meanQueueLength'])


@pytest.mark.parametrize("a,beta,s", BATCH_GRID)
def test_batch_closed_form_matches_brute_force_chain(a, beta, s):
    max_batch, n_max = 600, 3000
    batch = np.zeros(max_batch + 1)
    for k in range(1, max_batch + 1):
        batch[k] = beta * (1.0 - beta) ** (k - 1)
    # Fold the truncated tail onto the last supported size so the arrival law
    # stays a probability distribution.
    batch[max_batch] += 1.0 - batch.sum()

    pi = np.zeros(n_max + 1)
    pi[0] = 1.0
    for _ in range(200000):
        # Departure first.
        post = np.zeros(n_max + 1)
        post[0] += pi[0]
        post[0:n_max] += pi[1:n_max + 1] * s
        post[1:n_max + 1] += pi[1:n_max + 1] * (1.0 - s)
        # Then the batch arrival.
        nxt = post * (1.0 - a)
        for k in range(1, max_batch + 1):
            pk = batch[k]
            if pk == 0.0:
                continue
            shifted = np.zeros(n_max + 1)
            if k <= n_max:
                shifted[k:] = post[:n_max + 1 - k]
                shifted[n_max] += post[n_max + 1 - k:].sum()
            else:
                shifted[n_max] = post.sum()
            nxt += shifted * a * pk
        delta = np.abs(nxt - pi).sum()
        pi = nxt
        if delta < 1e-14:
            break

    r = dqsys_geoxgeo1(a, beta, s)
    mean = float((np.arange(pi.size) * pi).sum())
    label = "a=%s beta=%s s=%s" % (a, beta, s)
    assert abs(pi.sum() - 1.0) < 1e-8, "chain mass, " + label
    assert abs(pi[0] - r['boundaryEmptyProb']) < 1e-5, "p0, " + label
    assert abs(mean - r['meanQueueLength']) < 1e-4 * max(1.0, mean), "E[N], " + label


# ----------------------------------------------------------------------
# Argument validation
# ----------------------------------------------------------------------

def test_invalid_arguments_raise():
    with pytest.raises(ValueError):
        dqsys_geogeo1(0.5, 0.5)          # a = s is not stable
    with pytest.raises(ValueError):
        dqsys_geogeo1(0.7, 0.5)          # a > s is not stable
    with pytest.raises(ValueError):
        dqsys_geogeo1(0.0, 0.5)          # a must be positive
    with pytest.raises(ValueError):
        dqsys_geogeo1(0.2, 1.5)          # s must not exceed one
    with pytest.raises(ValueError):
        dqsys_geogeo1(0.2, 0.5, 'MIDSLOT')
    with pytest.raises(ValueError):
        dqsys_geogeo1(0.2, 0.5)['pmf'](-1)

    with pytest.raises(ValueError):
        dqsys_geoxgeo1(0.5, 0.5, 0.9)    # lambda = 1.0 > s
    with pytest.raises(ValueError):
        dqsys_geoxgeo1(0.1, 1.5, 0.9)    # beta must not exceed one
    with pytest.raises(ValueError):
        dqsys_geoxgeo1_moments(0.1, 0.5, 0.0, 0.9)   # E[X] < 1
    with pytest.raises(ValueError):
        dqsys_geoxgeo1_moments(0.1, 3.0, 0.0, 0.9)   # below E[X]^2-E[X]
    with pytest.raises(ValueError):
        dqsys_geoxgeo1(0.1, 0.5, 0.9)['pgf'](1.5, 1.0)
