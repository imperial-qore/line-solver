"""
Native-Python SSA tests for pass-and-swap (PAS) queues, validated against exact
CTMC-precision golden rows.

The PAS station has the order-independent product-form stationary distribution
(Dorsman & Gardner 2024). The native SSA solver drives the same
``State.afterEvent`` machinery as the CTMC solver, so PAS stations are simulated
through ``after_event_station_pas``. Golden QLen/Util/Tput values for the
PASQueue station rows match the MATLAB reference and agree across MATLAB,
Python-native, Java and Kotlin to CTMC precision.
"""

import numpy as np

from line_solver import Network, Source, Queue, Sink, OpenClass, Exp, SchedStrategy, SSA

# Absolute tolerance absorbing Monte-Carlo error at the SSA sample budget below.
SIM_TOL = 2e-2
SIM_SAMPLES = 200_000
SIM_SEED = 23000


def _build(name, lam, mu_fun, G, nservers, cap):
    model = Network(name)
    source = Source(model, 'Source')
    queue = Queue(model, 'PASQueue', SchedStrategy.PAS)
    sink = Sink(model, 'Sink')
    classes = [OpenClass(model, f'Class{r+1}') for r in range(len(lam))]
    for r in range(len(lam)):
        source.setArrival(classes[r], Exp(lam[r]))
    queue.setService(mu_fun)
    queue.setSwapGraph(G)
    queue.setNumberOfServers(nservers)
    queue.setCap(cap)
    P = model.initRoutingMatrix()
    for r in range(len(lam)):
        P[classes[r]] = Network.serialRouting(source, queue, sink)
    model.link(P)
    return model


def _ssa_queue_rows(model, cutoff, nclasses, samples=SIM_SAMPLES, seed=SIM_SEED):
    """PASQueue rows (QLen/Util/Tput) from the native SSA simulator.

    The native SSA solver drives the same ``State.afterEvent`` machinery as the
    CTMC solver, so PAS stations are simulated through ``after_event_station_pas``:
    the total rate mu(c) and the pass-and-swap departure class are evaluated on
    the actual ordered job list. Validated against the exact CTMC (JMT lacks PAS).
    """
    t = SSA(model, seed=seed, samples=samples, cutoff=cutoff).getAvgTable()
    q = np.asarray(t['QLen'], dtype=float)[-nclasses:]
    u = np.asarray(t['Util'], dtype=float)[-nclasses:]
    tp = np.asarray(t['Tput'], dtype=float)[-nclasses:]
    return q, u, tp


def _assert_sim(metric, got, ref, tol=SIM_TOL):
    assert np.allclose(np.asarray(got, dtype=float), np.asarray(ref, dtype=float), atol=tol), \
        f'SSA {metric} mismatch:\n got={np.round(np.asarray(got, float), 4)}\n ref={ref}'


def test_pas_mmk_ssa_matches_ctmc():
    """SSA reproduces the CTMC PASQueue rows for the M/M/K OI queue."""
    K = 2
    model = _build('PASmmk', [0.7, 0.5],
                   lambda c: float(min(len(np.atleast_1d(c)), K)),
                   np.zeros((2, 2)), K, 4)
    q, u, tp = _ssa_queue_rows(model, 4, 2)
    _assert_sim('QLen', q, [0.8032786885245902, 0.5737704918032788])
    _assert_sim('Util', u, [0.3248781568049634, 0.2320558262931842])
    _assert_sim('Tput', tp, [0.6497563136099268, 0.4641116525863684])


# Graph-aware SSA tests: the swapping graph (incl. self-loops) is exercised by
# the afterEvent pass-and-swap handler, validated against the CTMC golden rows.

def test_pas_selfloop_ssa_matches_ctmc():
    comp = np.eye(3)

    def mu_fun(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        return 0.0 if c.size == 0 else float(np.sum(np.any(comp[:, c], axis=1)))

    G = np.zeros((3, 3))
    G[0, 0] = 1
    G[1, 2] = 1
    G[2, 1] = 1
    model = _build('PASselfloop', [0.6, 0.4, 0.3], mu_fun, G, 3, 3)
    q, u, tp = _ssa_queue_rows(model, 3, 3, samples=120_000)
    _assert_sim('QLen', q, [0.7216685979, 0.4199304751, 0.2940903824], tol=3e-2)
    _assert_sim('Util', u, [0.1599073001, 0.1066048667, 0.07995365006], tol=3e-2)
    _assert_sim('Tput', tp, [0.4797219003, 0.3198146002, 0.2398609502], tol=3e-2)


def test_pas_compatibility_ssa_matches_ctmc():
    comp = np.array([[1, 0, 0, 1, 0],
                     [0, 1, 0, 1, 0],
                     [0, 0, 1, 0, 1]])

    def mu_fun(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        return 0.0 if c.size == 0 else float(np.sum(np.any(comp[:, c], axis=1)))

    G = np.zeros((5, 5))
    for i, j in [(0, 2), (0, 4), (1, 3), (2, 3), (3, 4)]:
        G[i, j] = 1
        G[j, i] = 1
    model = _build('PAScompatibility', [0.5, 0.4, 0.3, 0.2, 0.1], mu_fun, G, 3, 3)
    q, u, tp = _ssa_queue_rows(model, 3, 5, samples=90_000)
    _assert_sim('QLen', q, [0.561490291, 0.4214071462, 0.3053202126, 0.1123506815, 0.1017734042], tol=3e-2)
    _assert_sim('Tput', tp, [0.3909908962, 0.3127927169, 0.2345945377, 0.1563963585, 0.07819817923], tol=3e-2)
