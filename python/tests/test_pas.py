"""
Native-Python CTMC tests for pass-and-swap (PAS) / order-independent queues.

The PAS station has the order-independent product-form stationary distribution
(Dorsman & Gardner 2024). Golden QLen/Util/Tput values for the PASQueue station
rows match the MATLAB reference and agree across MATLAB, Python-native, Java and
Kotlin to CTMC precision. The three models mirror the examples in
``python/examples/advanced/passAndSwap/``.
"""

import numpy as np

from line_solver import Network, Source, Queue, Sink, OpenClass, Exp, SchedStrategy, CTMC, SSA

TOL = 1e-6
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


def _queue_rows(model, cutoff, nclasses):
    t = CTMC(model, cutoff=cutoff).getAvgTable()
    qlen = np.asarray(t['QLen'], dtype=float)[-nclasses:]
    util = np.asarray(t['Util'], dtype=float)[-nclasses:]
    tput = np.asarray(t['Tput'], dtype=float)[-nclasses:]
    return qlen, util, tput


def _assert(metric, got, ref):
    got = np.asarray(got, dtype=float)
    ref = np.asarray(ref, dtype=float)
    assert np.allclose(got, ref, atol=TOL), \
        f'PAS CTMC {metric} mismatch:\n got={got}\n ref={ref}'


def test_pas_mmk():
    K = 2
    model = _build('PASmmk', [0.7, 0.5],
                   lambda c: float(min(len(np.atleast_1d(c)), K)),
                   np.zeros((2, 2)), K, 4)
    q, u, tp = _queue_rows(model, 4, 2)
    _assert('QLen', q, [0.8032786885245902, 0.5737704918032788])
    _assert('Util', u, [0.3248781568049634, 0.2320558262931842])
    _assert('Tput', tp, [0.6497563136099268, 0.4641116525863684])


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


def test_pas_compatibility():
    comp = np.array([[1, 0, 0, 1, 0],
                     [0, 1, 0, 1, 0],
                     [0, 0, 1, 0, 1]])

    def mu_fun(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        if c.size == 0:
            return 0.0
        return float(np.sum(np.any(comp[:, c], axis=1)))

    G = np.zeros((5, 5))
    for i, j in [(0, 2), (0, 4), (1, 3), (2, 3), (3, 4)]:
        G[i, j] = 1
        G[j, i] = 1
    model = _build('PAScompatibility', [0.5, 0.4, 0.3, 0.2, 0.1], mu_fun, G, 3, 3)
    q, u, tp = _queue_rows(model, 3, 5)
    _assert('QLen', q, [0.561490291, 0.4214071462, 0.3053202126, 0.1123506815, 0.1017734042])
    _assert('Util', u, [0.1303302987, 0.104264239, 0.07819817923, 0.03366135172, 0.02606605974])
    _assert('Tput', tp, [0.3909908962, 0.3127927169, 0.2345945377, 0.1563963585, 0.07819817923])


def test_pas_selfloop():
    comp = np.eye(3)

    def mu_fun(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        if c.size == 0:
            return 0.0
        return float(np.sum(np.any(comp[:, c], axis=1)))

    G = np.zeros((3, 3))
    G[0, 0] = 1            # self-loop on class 1
    G[1, 2] = 1
    G[2, 1] = 1
    model = _build('PASselfloop', [0.6, 0.4, 0.3], mu_fun, G, 3, 3)
    q, u, tp = _queue_rows(model, 3, 3)
    _assert('QLen', q, [0.7216685979, 0.4199304751, 0.2940903824])
    _assert('Util', u, [0.1599073001, 0.1066048667, 0.07995365006])
    _assert('Tput', tp, [0.4797219003, 0.3198146002, 0.2398609502])


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
# ---------------------------------------------------------------------------

def _build_oi_mix(n_oi, use_delay, ps_servers):
    """Closed cyclic network: n_oi OI stations, an optional delay, and a PS queue
    with ps_servers servers (0 = no PS queue). Two classes, populations [2, 1]."""
    from line_solver import circul, Delay, ClosedClass
    model = Network('OImix')
    st = []
    oi = []
    for k in range(n_oi):
        q = Queue(model, 'OI%d' % (k + 1), SchedStrategy.OI)
        q.setCap(3)
        oi.append(q)
        st.append(q)
    d = None
    if use_delay:
        d = Delay(model, 'Delay')
        st.append(d)
    ps = None
    if ps_servers > 0:
        ps = Queue(model, 'PS', SchedStrategy.PS)
        ps.setNumberOfServers(ps_servers)
        st.append(ps)
    c1 = ClosedClass(model, 'C1', 2, st[0])
    c2 = ClosedClass(model, 'C2', 1, st[0])

    # Permutation-invariant (count-dependent) OI rates, distinct per station.
    def mu1(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        return 1.0 + 0.5 * np.sum(c == 0) + 0.2 * np.sum(c == 1)

    def mu2(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        return 1.5 * min(c.size, 2)

    def mu3(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        return 1.0 + 0.3 * c.size + 0.2 * np.sum(c == 1)

    oi[0].setService(mu1)
    if n_oi >= 2:
        oi[1].setService(mu2)
    if n_oi >= 3:
        oi[2].setService(mu3)
    if use_delay:
        d.setService(c1, Exp(2.0))
        d.setService(c2, Exp(1.0))
    if ps_servers > 0:
        ps.setService(c1, Exp(1.2))
        ps.setService(c2, Exp(0.8))
    P = model.initRoutingMatrix()
    P[c1] = circul(len(st))
    P[c2] = circul(len(st))
    model.link(P)
    return model


def _check_oi_mix(n_oi, use_delay, ps_servers):
    from line_solver import NC, MVA
    qc = np.array(CTMC(_build_oi_mix(n_oi, use_delay, ps_servers), cutoff=3).getAvgQLen())
    xc = np.array(CTMC(_build_oi_mix(n_oi, use_delay, ps_servers), cutoff=3).getAvgTput())
    qn = np.array(NC(_build_oi_mix(n_oi, use_delay, ps_servers)).getAvgQLen())
    xn = np.array(NC(_build_oi_mix(n_oi, use_delay, ps_servers)).getAvgTput())
    qm = np.array(MVA(_build_oi_mix(n_oi, use_delay, ps_servers)).getAvgQLen())
    xm = np.array(MVA(_build_oi_mix(n_oi, use_delay, ps_servers)).getAvgTput())
    np.testing.assert_allclose(qn, qc, atol=1e-9)
    np.testing.assert_allclose(xn, xc, atol=1e-9)
    np.testing.assert_allclose(qm, qc, atol=1e-9)
    np.testing.assert_allclose(xm, xc, atol=1e-9)


def test_oi_two_plus_is_ps():
    """Two OI stations + delay + single-server PS."""
    _check_oi_mix(2, True, 1)


def test_oi_three_plus_is():
    """Three OI stations + delay, no PS queue."""
    _check_oi_mix(3, True, 0)


def test_oi_two_plus_ps_no_delay():
    """Two OI stations + single-server PS, no delay node."""
    _check_oi_mix(2, False, 1)


def test_oi_multiserver_queue():
    """Two OI stations + delay + a MULTISERVER (c=2) PS queue. pfqn_mvaoi models
    an LI queue as a single server, so the multiserver station must be promoted to
    the OI representation mu(n) = (min(|n|,c)/|n|) sum_r n_r/D_r; otherwise MVA
    silently returns the c=1 answer."""
    _check_oi_mix(2, True, 2)


def test_oi_mva_rejects_non_product_form():
    """An OI station combined with a class-DEPENDENT FCFS queue is not product
    form, so it falls outside the exact OI path. AMVA cannot represent an OI rank
    rate mu(n) (it only sees sn.rates), so SolverMVA must refuse rather than
    silently return a zero queue-length at the OI station."""
    import pytest
    from line_solver import circul, Delay, ClosedClass, MVA
    from line_solver.solvers.solver_nc.solver_nc_oi_analyzer import nc_is_oi_model

    model = Network('OInonPF')
    oi1 = Queue(model, 'OI1', SchedStrategy.OI)
    d = Delay(model, 'Delay')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    c1 = ClosedClass(model, 'C1', 2, oi1)
    c2 = ClosedClass(model, 'C2', 1, oi1)

    def mu_fun(c):
        c = np.atleast_1d(np.asarray(c, dtype=int))
        return 1.0 + 0.5 * np.sum(c == 0) + 0.2 * np.sum(c == 1)

    oi1.setService(mu_fun)
    oi1.setCap(3)
    d.setService(c1, Exp(2.0))
    d.setService(c2, Exp(1.0))
    q.setService(c1, Exp(1.2))   # class-dependent FCFS
    q.setService(c2, Exp(0.6))
    P = model.initRoutingMatrix()
    P[c1] = circul(3)
    P[c2] = circul(3)
    model.link(P)

    assert not nc_is_oi_model(model.getStruct())
    with pytest.raises(ValueError):
        MVA(model).getAvgQLen()
