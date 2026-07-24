"""Regression tests for get_moment_table and get_moment_station_table, the
solver-level views over the pfqn_sens_* moment family.

Mirrors line-test.git/test_moment_table.m one for one.

The underlying algorithms are validated against brute-force enumeration, the
published tables of Strelen (1990) and simulation by the api-level harnesses
(test_pfqn_sens_mva.py, test_pfqn_sens_mom.py, test_pfqn_sens_respt.py,
test_pfqn_sens_linearizer.py). These tests therefore do NOT re-check the
mathematics. They check the solver-level plumbing, which is where a table method
can go wrong independently of correct algorithms:

  - that the demands, visit ratios and service times handed to the api layer are
    the ones the model actually expresses (asserted by requiring the table's
    means to reproduce get_avg_table, which is computed by a different code
    path);
  - that the per-class and per-station-total tables are mutually consistent,
    Var[Q_i] = sum over class pairs of the per-class covariances;
  - that RespTVar is populated at exactly the FCFS stations and NaN elsewhere,
    since the sojourn-time distribution is unknown at PS/LCFS centers and a
    number there would be wrong rather than missing;
  - that the dispatch by model type (closed / mixed / purely open, single-server
    / multiserver) picks a path that runs and agrees.
"""

import os

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, Source, Sink, ClosedClass,
                         OpenClass, Exp, SchedStrategy, SolverMVA)

# The numerical-derivative moment oracle (solveMeansForStruct) injects the
# perturbed demand column into model._sn and re-solves. Under lang=java the
# solver re-serializes the model's distributions to JSON for the jar, which does
# not carry the struct-level perturbation, so every central difference is exactly
# zero and the variance/skew collapse to 0. The oracle is a native-analytical
# path with no JSON-model counterpart; the hand-differentiated default-method
# tables (which operate on the struct natively) remain exercised under lang=java.
_skip_java_fd_oracle = pytest.mark.skipif(
    os.environ.get('LINE_SOLVER_LANG') == 'java',
    reason="finite-difference moment oracle perturbs the sn struct; the "
           "JSON-model java dispatch cannot serialize it (native-only path)")
from line_solver.api.pfqn.sens_mva import pfqn_sens_mva
from line_solver.api.pfqn.sens_mom import pfqn_sens_mom


# ---------- helpers ------------------------------------------------------

def closed_mixed_sched():
    """Closed, two classes. Q1 is FCFS with class-independent rates (as BCMP
    requires of FCFS); Q2 is PS with class-dependent rates, which is legal and
    which the FCFS gate must tolerate."""
    m = Network('mt_closed')
    d = Delay(m, 'Think')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', 3, d, 0)
    c2 = ClosedClass(m, 'C2', 2, d, 0)
    d.setService(c1, Exp(1 / 1.0))
    d.setService(c2, Exp(1 / 0.5))
    q1.setService(c1, Exp(1 / 0.4))
    q1.setService(c2, Exp(1 / 0.4))
    q2.setService(c1, Exp(1 / 0.3))
    q2.setService(c2, Exp(1 / 0.2))
    P = m.initRoutingMatrix()
    P[c1] = Network.serialRouting(d, q1, q2)
    P[c2] = Network.serialRouting(d, q1, q2)
    m.link(P)
    return m


def row_of(T, station, jobclass):
    hit = T[(T['Station'] == station) & (T['JobClass'] == jobclass)]
    assert len(hit) == 1, "row %s/%s not found" % (station, jobclass)
    return hit.iloc[0]


def avg_row(A, station, jobclass):
    # get_avg_table returns an IndexedTable (MATLAB-style formatting wrapper);
    # filter on the underlying DataFrame.
    A = getattr(A, 'data', A)
    hit = A[(A['Station'].astype(str) == station)
            & (A['JobClass'].astype(str) == jobclass)]
    assert len(hit) == 1, "avg row %s/%s not found" % (station, jobclass)
    return hit.iloc[0]


# ---------- tests --------------------------------------------------------

def test_means_reproduce_get_avg_table():
    """The table is only meaningful if the parameters handed to the api layer
    are the model's. get_avg_table reaches the same means by a different path,
    so agreement pins the translation."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    T, _ = s.get_moment_table()
    A = s.get_avg_table()
    for _, row in T.iterrows():
        a = avg_row(A, row['Station'], row['JobClass'])
        assert row['QLen'] == pytest.approx(float(a['QLen']), rel=1e-9)
        assert row['RespT'] == pytest.approx(float(a['RespT']), rel=1e-9)


def test_perclass_variance_matches_api():
    """The QLenVar column must be exactly what pfqn_sens_mva reports; nothing in
    the table method may rescale it."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    T, mom = s.get_moment_table()
    ref = pfqn_sens_mva(np.array([[0.4, 0.4], [0.3, 0.2]]),
                        np.array([3, 2]), np.array([1.0, 0.5]))
    np.testing.assert_allclose(mom['qlen'].QVar, ref.QVar, rtol=1e-9)
    r = row_of(T, 'Q1', 'C1')
    assert r['QLenVar'] == pytest.approx(ref.QVar[0, 0], rel=1e-9)
    assert r['QLenSCV'] == pytest.approx(ref.QVar[0, 0] / ref.Q[0, 0] ** 2,
                                         rel=1e-9)


def test_station_total_consistent_with_perclass():
    """Var[Q_i] of the station table must equal the sum of the per-class
    covariance block of the per-class table. The two come from different
    recursions (Strelen's column scaling vs de Souza e Silva-Muntz's per-class
    one), so this is a real cross-check, not a tautology."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    _, mom = s.get_moment_table()
    S, _ = s.get_moment_station_table()
    for ist in range(2):
        tot = float(np.sum(mom['qlen'].QCov[ist, :, :]))
        assert S['QLenVar'][ist] == pytest.approx(tot, rel=1e-8)
        assert S['QLen'][ist] == pytest.approx(
            float(np.sum(mom['qlen'].Q[ist, :])), rel=1e-9)


def test_respt_variance_only_at_fcfs():
    """RespTVar must be a number at the FCFS station and NaN at the PS station.
    The sojourn-time distribution at a PS center is not known in general, so a
    value there would be fabricated."""
    m = closed_mixed_sched()
    T, _ = SolverMVA(m).get_moment_table()
    for cl in ['C1', 'C2']:
        rf = row_of(T, 'Q1', cl)
        assert not np.isnan(rf['RespTVar']), "FCFS station must report RespTVar"
        assert rf['RespTVar'] > 0
        rp = row_of(T, 'Q2', cl)
        assert np.isnan(rp['RespTVar']), "PS station must NOT report RespTVar"
        assert np.isnan(rp['RespTSCV'])


def test_respt_mean_from_theorem41_matches_mva():
    """At the FCFS station the mean RespT is overwritten by the t=1 case of
    Strelen Theorem 4.1, which is a different expression from Little's law. They
    must agree, which also validates the service-time/visit-ratio factorization
    the table hands to pfqn_sens_respt."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    T, _ = s.get_moment_table()
    A = s.get_avg_table()
    for cl in ['C1', 'C2']:
        rf = row_of(T, 'Q1', cl)
        a = avg_row(A, 'Q1', cl)
        assert rf['RespT'] == pytest.approx(float(a['RespT']), rel=1e-8)


def test_open_mm1_closed_form():
    """A purely open M/M/1 has no population lattice, so the table takes the
    exact BCMP branch. Every moment is known in closed form: n ~ Geom(rho)."""
    lam, svc = 0.5, 1.0
    rho = lam * svc
    m = Network('mt_open')
    src = Source(m, 'S')
    snk = Sink(m, 'K')
    q = Queue(m, 'O1', SchedStrategy.FCFS)
    oc = OpenClass(m, 'O', 0)
    src.setArrival(oc, Exp(lam))
    q.setService(oc, Exp(1 / svc))
    m.link(Network.serialRouting(src, q, snk))
    T, _ = SolverMVA(m).get_moment_table()
    r = row_of(T, 'O1', 'O')
    assert r['QLen'] == pytest.approx(rho / (1 - rho), rel=1e-10)
    assert r['QLenVar'] == pytest.approx(rho / (1 - rho) ** 2, rel=1e-10)
    assert r['RespT'] == pytest.approx(svc / (1 - rho), rel=1e-10)


def test_open_multiclass_mm1_closed_form():
    """Two open classes at one station: the total is geometric in the aggregate
    utilization and the classes are multinomial given the total, so
    Var[n_r] = E[n](p_r - p_r^2) + p_r^2 Var[n]."""
    l1, l2, s1, s2 = 0.3, 0.2, 1.0, 1.5
    m = Network('mt_open2')
    src = Source(m, 'S')
    snk = Sink(m, 'K')
    q = Queue(m, 'O1', SchedStrategy.PS)
    o1 = OpenClass(m, 'O1c', 0)
    o2 = OpenClass(m, 'O2c', 0)
    src.setArrival(o1, Exp(l1))
    src.setArrival(o2, Exp(l2))
    q.setService(o1, Exp(1 / s1))
    q.setService(o2, Exp(1 / s2))
    P = m.initRoutingMatrix()
    P[o1] = Network.serialRouting(src, q, snk)
    P[o2] = Network.serialRouting(src, q, snk)
    m.link(P)
    T, _ = SolverMVA(m).get_moment_table()
    rho1, rho2 = l1 * s1, l2 * s2
    rho = rho1 + rho2
    En = rho / (1 - rho)
    Vn = rho / (1 - rho) ** 2
    p1 = rho1 / rho
    r = row_of(T, 'O1', 'O1c')
    assert r['QLen'] == pytest.approx(p1 * En, rel=1e-9)
    assert r['QLenVar'] == pytest.approx(En * (p1 - p1 ** 2) + p1 ** 2 * Vn,
                                         rel=1e-9)


def test_linearizer_tracks_exact():
    """The Linearizer is an approximation, so it is banded, not equated. The
    bands are the accuracy the reference claims."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    # The algorithm is a solver property, not a table argument: the Linearizer
    # is selected by constructing the solver with method 'lin'.
    E, _ = s.get_moment_station_table(3)
    L, _ = SolverMVA(m, method='lin').get_moment_station_table(3)
    for k in range(len(E)):
        assert L['QLen'][k] == pytest.approx(E['QLen'][k], rel=0.021)
        assert L['QLenM3'][k] == pytest.approx(E['QLenM3'][k], rel=0.062)


def test_multiserver_closed_dispatches():
    """A multiserver closed model must route through pfqn_sens_mvaldmx rather
    than pfqn_sens_mva, and still reproduce get_avg_table's means."""
    m = Network('mt_ms')
    d = Delay(m, 'Think')
    q1 = Queue(m, 'Q1', SchedStrategy.PS)
    q1.setNumberOfServers(2)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', 4, d, 0)
    d.setService(c1, Exp(1 / 1.0))
    q1.setService(c1, Exp(1 / 0.4))
    q2.setService(c1, Exp(1 / 0.3))
    m.link(Network.serialRouting(d, q1, q2))
    s = SolverMVA(m)
    T, _ = s.get_moment_table()
    A = s.get_avg_table()
    for _, row in T.iterrows():
        a = avg_row(A, row['Station'], row['JobClass'])
        assert row['QLen'] == pytest.approx(float(a['QLen']), rel=1e-8)
        assert row['QLenVar'] > 0


def test_station_table_rejects_open():
    """The higher-moment recursion of the reference is stated for closed
    networks; the method must say so rather than return something."""
    m = Network('mt_open3')
    src = Source(m, 'S')
    snk = Sink(m, 'K')
    q = Queue(m, 'O1', SchedStrategy.FCFS)
    oc = OpenClass(m, 'O', 0)
    src.setArrival(oc, Exp(0.5))
    q.setService(oc, Exp(1 / 1.0))
    m.link(Network.serialRouting(src, q, snk))
    with pytest.raises(ValueError):
        SolverMVA(m).get_moment_station_table()


def test_order_selects_columns():
    """order is a set of moment orders. A scalar k means 1:k; a vector is
    literal."""
    m = closed_mixed_sched()
    s = SolverMVA(m)

    T1, _ = s.get_moment_table(1)
    assert list(T1.columns) == ['Station', 'JobClass', 'QLen', 'RespT']

    T2, _ = s.get_moment_table()                    # default is 2
    assert list(T2.columns) == ['Station', 'JobClass', 'QLen', 'QLenVar',
                                'QLenSCV', 'RespT', 'RespTVar', 'RespTSCV']
    T12, _ = s.get_moment_table([1, 2])
    assert list(T12.columns) == list(T2.columns), \
        "scalar 2 must equal the set [1 2]"

    T3, _ = s.get_moment_table(3)
    assert 'RespTSkew' in T3.columns
    # order 3 DOES add a per-class queue-length skewness: scaling L(i,r) alone is
    # Akyildiz-Strelen Theorem 1 with T={r}, so the quantity exists
    assert 'QLenSkew' in T3.columns

    T23, _ = s.get_moment_table([2, 3])   # second moments and skew, no means
    assert 'QLen' not in T23.columns
    assert 'QLenVar' in T23.columns
    assert 'RespTSkew' in T23.columns

    S1, _ = s.get_moment_station_table(1)
    assert list(S1.columns) == ['Station', 'QLen']
    S2, _ = s.get_moment_station_table()
    assert list(S2.columns) == ['Station', 'QLen', 'QLenVar', 'QLenSCV']
    S3, _ = s.get_moment_station_table(3)
    assert 'QLenM3' in S3.columns
    assert 'QLenSkew' in S3.columns
    S13, _ = s.get_moment_station_table([1, 3])
    assert list(S13.columns) == ['Station', 'QLen', 'QLenM3', 'QLenSkew']

    # the values must not depend on which columns were asked for
    np.testing.assert_allclose(S3['QLen'], S2['QLen'], rtol=1e-12)
    np.testing.assert_allclose(T3['QLen'], T2['QLen'], rtol=1e-12)
    np.testing.assert_allclose(T3['RespTVar'], T2['RespTVar'], rtol=1e-12,
                               equal_nan=True)


def test_order_rejects_bad_values():
    m = closed_mixed_sched()
    s = SolverMVA(m)
    with pytest.raises(ValueError):
        s.get_moment_table(0)
    with pytest.raises(ValueError):
        s.get_moment_table(4)
    with pytest.raises(ValueError):
        s.get_moment_table(2.5)
    with pytest.raises(ValueError):
        s.get_moment_station_table([1, 5])
    # non-integers must be rejected on the VECTOR path too: rounding [1 2.5] to
    # [1 3] would silently answer a question that was not asked
    with pytest.raises(ValueError):
        s.get_moment_table([1, 2.5])


def test_perclass_skew_exists_and_matches_api():
    """Per-class queue-length skewness IS defined: scaling L(i,r) alone is the
    class-subset parameter T={r} of Akyildiz-Strelen Theorem 1. Order 3 must
    therefore expose QLenSkew, and it must equal pfqn_sens_mom with
    groups=1:R."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    T, _ = s.get_moment_table(3)
    assert 'QLenSkew' in T.columns
    ref = pfqn_sens_mom(np.array([[0.4, 0.4], [0.3, 0.2]]), np.array([3, 2]),
                        np.array([1.0, 0.5]), np.array([1, 1]),
                        np.array([1, 2]))
    r = row_of(T, 'Q1', 'C1')
    assert r['QLenSkew'] == pytest.approx(ref.Skew[0, 0], rel=1e-9)
    assert not np.isnan(r['QLenSkew'])


def test_chain_table_consistent_with_perclass():
    """Each chain's variance must equal the sum of the per-class covariance
    block over the classes of that chain, since Var[sum_{r in c} n_ir] =
    sum_{r,s in c} Cov. The two come from different groupings of the same
    recursion."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    _, mom = s.get_moment_table()
    sn = m.getStruct()
    chains = np.atleast_2d(np.asarray(sn.chains))
    Ch, _ = s.get_moment_chain_table(3)
    # both classes are closed and each forms its own chain in this model
    for _, row in Ch.iterrows():
        ist = ['Q1', 'Q2'].index(row['Station'])
        c = int(row['Chain'].replace('Chain', '')) - 1   # 1-based label
        cls = np.nonzero(chains[c, :])[0]
        tot = 0.0
        for r in cls:
            for s2 in cls:
                tot += mom['qlen'].QCov[ist, r, s2]
        assert row['QLenVar'] == pytest.approx(tot, rel=1e-8)
        assert row['QLen'] == pytest.approx(
            float(np.sum(mom['qlen'].Q[ist, cls])), rel=1e-9)


def test_single_chain_equals_station_total():
    """When every class sits in one chain, the chain grouping IS the station
    total, so get_moment_chain_table must reproduce get_moment_station_table
    exactly."""
    m = Network('mt_1chain')
    d = Delay(m, 'Think')
    q1 = Queue(m, 'Q1', SchedStrategy.PS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', 3, d, 0)
    c2 = ClosedClass(m, 'C2', 2, d, 0)
    d.setService(c1, Exp(1 / 1.0))
    d.setService(c2, Exp(1 / 0.5))
    q1.setService(c1, Exp(1 / 0.4))
    q1.setService(c2, Exp(1 / 0.6))
    q2.setService(c1, Exp(1 / 0.3))
    q2.setService(c2, Exp(1 / 0.2))
    # a class switch merges the two classes into a single chain. The python
    # RoutingMatrix takes a whole (nnodes x nnodes) block per class pair, so the
    # blocks are built here rather than assigned entry by entry as in MATLAB.
    # Network.getNodeIndex is 1-based (MATLAB-compatible); numpy is 0-based
    ti = m.getNodeIndex('Think') - 1
    q1i = m.getNodeIndex('Q1') - 1
    q2i = m.getNodeIndex('Q2') - 1
    nn = 3

    def block(*pairs):
        A = np.zeros((nn, nn))
        for (i, j) in pairs:
            A[i, j] = 1.0
        return A

    P = m.initRoutingMatrix()
    P[c1, c1] = block((ti, q1i), (q1i, q2i))
    P[c1, c2] = block((q2i, ti))     # C1 switches to C2 on the way back
    P[c2, c2] = block((ti, q1i), (q1i, q2i))
    P[c2, c1] = block((q2i, ti))     # and back again: one chain
    m.link(P)
    s = SolverMVA(m)
    sn = m.getStruct()
    assert np.atleast_2d(np.asarray(sn.chains)).shape[0] == 1, \
        "the model must have a single chain"
    Ch, _ = s.get_moment_chain_table(3)
    St, _ = s.get_moment_station_table(3)
    for col in ['QLen', 'QLenVar', 'QLenM3', 'QLenSkew']:
        np.testing.assert_allclose(Ch[col], St[col], rtol=1e-9)


# ---------- MATLAB cross-codebase parity --------------------------------

def test_matlab_parity_closed_mixed_sched():
    """Numerical agreement with MATLAB's getMomentTable /
    getMomentStationTable on the shared closedMixedSched model. Goldens are
    MATLAB's output (ground truth) for the same model."""
    m = closed_mixed_sched()
    s = SolverMVA(m)
    T, _ = s.get_moment_table()
    ml_T = [
        ('Q1', 'C1', 1.2590868104973123, 0.95915186471605907,
         0.60502876990150334, 1.1097277988535983, 0.67931527771890288,
         0.55161802791970638),
        ('Q1', 'C2', 1.0851704886408977, 0.58914886028854285,
         0.50029837557321333, 1.0257577857411551, 0.63730363144688584,
         0.605698849950644),
        ('Q2', 'C1', 0.60632251618801658, 0.63676521349464243,
         1.732095913838531, 0.5343975853570736, np.nan, np.nan),
        ('Q2', 'C2', 0.38586911557872094, 0.36619760863520778,
         2.4594356539929212, 0.36474291710388251, np.nan, np.nan),
    ]
    for (st, cl, q, qv, qs, rt, rv, rs) in ml_T:
        r = row_of(T, st, cl)
        assert r['QLen'] == pytest.approx(q, rel=1e-12)
        assert r['QLenVar'] == pytest.approx(qv, rel=1e-12)
        assert r['QLenSCV'] == pytest.approx(qs, rel=1e-12)
        assert r['RespT'] == pytest.approx(rt, rel=1e-12)
        if np.isnan(rv):
            assert np.isnan(r['RespTVar']) and np.isnan(r['RespTSCV'])
        else:
            assert r['RespTVar'] == pytest.approx(rv, rel=1e-12)
            assert r['RespTSCV'] == pytest.approx(rs, rel=1e-12)

    S, _ = s.get_moment_station_table(3)
    ml_S = [
        ('Q1', 2.34425729913821, 2.0485185486579778, 0.372760037605817,
         27.306873270057903, 0.0058474211658829474),
        ('Q2', 0.99219163176673753, 1.3177205563529943, 1.33854261180524,
         6.5608190827471295, 1.0985892990686912),
    ]
    for k, (st, q, qv, qs, m3, sk) in enumerate(ml_S):
        assert S['Station'][k] == st
        assert S['QLen'][k] == pytest.approx(q, rel=1e-12)
        assert S['QLenVar'][k] == pytest.approx(qv, rel=1e-12)
        assert S['QLenSCV'][k] == pytest.approx(qs, rel=1e-12)
        assert S['QLenM3'][k] == pytest.approx(m3, rel=1e-12)
        assert S['QLenSkew'][k] == pytest.approx(sk, rel=1e-10)

    # the Linearizer path is reached by the solver's method, set at construction
    L, _ = SolverMVA(m, method='lin').get_moment_station_table(3)
    ml_L = [
        ('Q1', 2.3439432294151681, 2.0434312059450681, 27.3802095587275),
        ('Q2', 0.99343471784425141, 1.309875354182517, 6.4821375562165207),
    ]
    for k, (st, q, qv, m3) in enumerate(ml_L):
        assert L['QLen'][k] == pytest.approx(q, rel=1e-10)
        assert L['QLenVar'][k] == pytest.approx(qv, rel=1e-10)
        assert L['QLenM3'][k] == pytest.approx(m3, rel=1e-10)


def test_matlab_parity_order3():
    """order=3 drives real work, not just display: it raises tmax to 3 in
    pfqn_sens_respt and calls pfqn_sens_mom with groups=1:R. RespTSkew and
    QLenSkew are the columns it adds, and both must match MATLAB. Column lists
    are MATLAB's too."""
    s = SolverMVA(closed_mixed_sched())
    T3, _ = s.get_moment_table(3)
    # (RespTSkew, QLenSkew) goldens straight from MATLAB's getMomentTable(3)
    ml = {('Q1', 'C1'): (1.162786468207976, 0.23040328361611218),
          ('Q1', 'C2'): (1.2423819292133318, -0.14591013674776285),
          ('Q2', 'C1'): (np.nan, 1.1468766640351271),
          ('Q2', 'C2'): (np.nan, 1.3184573516872069)}
    for (st, cl), (rsk, qsk) in ml.items():
        r = row_of(T3, st, cl)
        if np.isnan(rsk):
            assert np.isnan(r['RespTSkew'])
        else:
            assert r['RespTSkew'] == pytest.approx(rsk, rel=1e-12)
        # QLenSkew is reported at PS stations too: it needs no sojourn-time
        # distribution, only the queue-length law
        assert r['QLenSkew'] == pytest.approx(qsk, rel=1e-12)

    # MATLAB's exact column lists
    assert list(s.get_moment_table(1)[0].columns) == \
        ['Station', 'JobClass', 'QLen', 'RespT']
    assert list(T3.columns) == \
        ['Station', 'JobClass', 'QLen', 'QLenVar', 'QLenSCV', 'QLenSkew',
         'RespT', 'RespTVar', 'RespTSCV', 'RespTSkew']
    assert list(s.get_moment_table([2, 3])[0].columns) == \
        ['Station', 'JobClass', 'QLenVar', 'QLenSCV', 'QLenSkew', 'RespTVar',
         'RespTSCV', 'RespTSkew']
    assert list(s.get_moment_station_table([1, 3])[0].columns) == \
        ['Station', 'QLen', 'QLenM3', 'QLenSkew']


def test_matlab_parity_chain_table():
    """get_moment_chain_table against MATLAB's getMomentChainTable(3). Each
    class is its own chain in this model, so the chain moments coincide with the
    per-class ones -- which is itself the cross-check that the groups map is
    built and compacted correctly."""
    s = SolverMVA(closed_mixed_sched())
    Ch, _ = s.get_moment_chain_table(3)
    assert list(Ch.columns) == \
        ['Station', 'Chain', 'QLen', 'QLenVar', 'QLenSCV', 'QLenM3', 'QLenSkew']
    ml = [
        ('Q1', 'Chain1', 1.2590868104973123, 0.959151864716059,
         0.60502876990150323, 5.8354273290793657, 0.23040328361611218),
        ('Q1', 'Chain2', 1.0851704886408977, 0.58914886028854263,
         0.50029837557321311, 3.1298905718352059, -0.14591013674776285),
        ('Q2', 'Chain1', 0.60632251618801658, 0.63676521349464232,
         1.7320959138385306, 1.9639103838939835, 1.1468766640351271),
        ('Q2', 'Chain2', 0.38586911557872094, 0.36619760863520778,
         2.4594356539929212, 0.77353951782069441, 1.3184573516872069),
    ]
    assert len(Ch) == len(ml)
    for k, (st, ch, q, qv, qs, m3, sk) in enumerate(ml):
        assert Ch['Station'][k] == st
        assert Ch['Chain'][k] == ch
        assert Ch['QLen'][k] == pytest.approx(q, rel=1e-12)
        assert Ch['QLenVar'][k] == pytest.approx(qv, rel=1e-12)
        assert Ch['QLenSCV'][k] == pytest.approx(qs, rel=1e-12)
        assert Ch['QLenM3'][k] == pytest.approx(m3, rel=1e-12)
        assert Ch['QLenSkew'][k] == pytest.approx(sk, rel=1e-12)


def test_matlab_parity_sens_mom_groups():
    """The api-level groups=1:R call, against MATLAB's pfqn_sens_mom with the
    same arguments. This is the primitive the QLenSkew column rests on."""
    ref = pfqn_sens_mom(np.array([[0.4, 0.4], [0.3, 0.2]]), np.array([3, 2]),
                        np.array([1.0, 0.5]), np.array([1, 1]),
                        np.array([1, 2]))
    ml_m = np.array([[1.2590868104973123, 1.0851704886408977],
                     [0.60632251618801658, 0.38586911557872094]])
    ml_var = np.array([[0.959151864716059, 0.58914886028854263],
                       [0.63676521349464232, 0.36619760863520778]])
    ml_m3 = np.array([[5.8354273290793657, 3.1298905718352059],
                      [1.9639103838939835, 0.77353951782069441]])
    ml_skew = np.array([[0.23040328361611218, -0.14591013674776285],
                        [1.1468766640351271, 1.3184573516872069]])
    np.testing.assert_allclose(ref.m, ml_m, rtol=1e-12)
    np.testing.assert_allclose(ref.Var, ml_var, rtol=1e-12)
    np.testing.assert_allclose(ref.M3, ml_m3, rtol=1e-12)
    np.testing.assert_allclose(ref.Skew, ml_skew, rtol=1e-12)


def test_method_comes_from_the_solver_not_the_call():
    """The algorithm is a property of the solver object, set at construction,
    not an argument of the table call. Two solvers over the same model must
    therefore disagree exactly as exact-vs-approximate, and neither table method
    may accept a method argument."""
    m = closed_mixed_sched()
    ex, _ = SolverMVA(m).get_moment_station_table(3)
    lin, _ = SolverMVA(m, method='lin').get_moment_station_table(3)
    # same model, different solver method -> same quantity, different algorithm
    for k in range(len(ex)):
        assert lin['QLen'][k] == pytest.approx(ex['QLen'][k], rel=0.021)
        assert lin['QLenM3'][k] == pytest.approx(ex['QLenM3'][k], rel=0.062)
    # and they must not be bit-identical, else the method is being ignored
    assert not np.array_equal(np.asarray(lin['QLenM3']),
                              np.asarray(ex['QLenM3']))

    # passing a method positionally must NOT be silently accepted as an order
    with pytest.raises(TypeError):
        SolverMVA(m).get_moment_station_table(3, 'lin')
    with pytest.raises(TypeError):
        SolverMVA(m).get_moment_chain_table(3, 'lin')


def test_chain_table_rejects_linearizer_when_multichain():
    """A Linearizer-family solver approximates the per-station totals only, so
    it cannot express a per-chain grouping. closed_mixed_sched has two chains,
    so the chain table must refuse rather than return the station totals under a
    chain label."""
    m = closed_mixed_sched()
    with pytest.raises(ValueError, match='per-chain grouping'):
        SolverMVA(m, method='lin').get_moment_chain_table(3)


def test_nonproductform_is_refused():
    """The moment identity Cov = L dQ/dL is a theorem about the product form.
    Outside it, L dQ/dL is still computable and is simply NOT a covariance, so
    returning a number would be a confident wrong answer. Both tables must
    refuse. A heterogeneous-FCFS model (class-dependent rates at an FCFS
    station) is the canonical non-product-form case."""
    from line_solver.api.sn.predicates import sn_has_product_form
    m = Network('mt_npf')
    d = Delay(m, 'Think')
    f = Queue(m, 'F', SchedStrategy.FCFS)
    k1 = ClosedClass(m, 'C1', 2, d, 0)
    k2 = ClosedClass(m, 'C2', 1, d, 0)
    d.setService(k1, Exp(1 / 1.0))
    d.setService(k2, Exp(1 / 1.0))
    f.setService(k1, Exp(1 / 0.4))
    f.setService(k2, Exp(1 / 0.9))      # heterogeneous FCFS
    P = m.initRoutingMatrix()
    P[k1] = Network.serialRouting(d, f)
    P[k2] = Network.serialRouting(d, f)
    m.link(P)
    assert not sn_has_product_form(m.getStruct()), \
        'the fixture must be non-product-form'
    with pytest.raises(ValueError, match='requires a product-form model'):
        SolverMVA(m).get_moment_station_table(2)
    with pytest.raises(ValueError, match='requires a product-form model'):
        SolverMVA(m).get_moment_chain_table(2)


@_skip_java_fd_oracle
def test_fd_fallback_tracks_exact():
    """A product-form model solved by an approximate MVA method has no
    hand-differentiated moment implementation, so the solver's own means are
    differentiated numerically and the identity applied. The moments must track
    the exact ones, and must inherit the accuracy of THAT method's means: 'amva'
    has near-exact means here, so its moments must be near-exact too."""
    m = closed_mixed_sched()
    ex, _ = SolverMVA(m).get_moment_station_table(3)
    fd, _ = SolverMVA(m, method='amva').get_moment_station_table(3)
    for k in range(len(ex)):
        assert fd['QLen'][k] == pytest.approx(ex['QLen'][k], rel=0.02)
        assert fd['QLenVar'][k] == pytest.approx(ex['QLenVar'][k], rel=0.05)
        assert fd['QLenM3'][k] == pytest.approx(ex['QLenM3'][k], rel=0.05)
    # the numerical path must not be silently returning the analytic answer
    assert not np.array_equal(np.asarray(fd['QLenM3']),
                              np.asarray(ex['QLenM3']))
    # and the chain table must take the same path
    fdc, _ = SolverMVA(m, method='amva').get_moment_chain_table(3)
    exc, _ = SolverMVA(m).get_moment_chain_table(3)
    assert len(fdc) == len(exc)


@_skip_java_fd_oracle
def test_oracle_covers_every_solver_and_method():
    """The identity needs only the mean queue lengths as a function of the
    demands, so the numerical-derivative oracle re-runs THIS solver, whatever it
    is. That makes every method of every product-form solver usable, not just
    those with a hand-differentiated implementation. Restricting it to one
    analyzer would have restricted the moments to that analyzer's methods for no
    mathematical reason."""
    from line_solver import SolverNC
    m = closed_mixed_sched()
    ex, _ = SolverMVA(m).get_moment_station_table(3)

    # 'sum' is a SolverMVA method that no hand-differentiated implementation
    # reaches; it must now produce moments rather than an error.
    T, _ = SolverMVA(m, method='sum').get_moment_station_table(3)
    assert len(T) == len(ex)
    assert (np.asarray(T['QLenVar']) > 0).all()

    # SolverNC's convolution algorithm is EXACT, so differentiating its means
    # must reproduce the exact moments. This is the strongest evidence the
    # oracle is sound: a different solver, a different algorithm, the same
    # answer.
    NCca, _ = SolverNC(m, method='ca').get_moment_station_table(3)
    for k in range(len(ex)):
        assert NCca['QLen'][k] == pytest.approx(ex['QLen'][k], rel=1e-6)
        assert NCca['QLenVar'][k] == pytest.approx(ex['QLenVar'][k], rel=1e-4)
        assert NCca['QLenM3'][k] == pytest.approx(ex['QLenM3'][k], rel=1e-4)
