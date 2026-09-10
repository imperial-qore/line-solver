"""Regression tests for Network.set_global_dependence and the CTMC generator.

set_global_dependence declares a globally state-dependent rate scaling phi(n)
over the FULL (nstations, nclasses) population matrix, not the population local
to one station. Every existing hook -- lldscaling, cdscaling, jdscaling -- is
handed one station's slice, so none of them can express a rate that reads the
whole state. A Whittle network needs exactly that, and so does bandwidth
sharing, where one route holds several links at once.

Oracles, in increasing strength:
  - an identity phi must reproduce the model declared without it, bit for bit;
  - a phi reproducing a per-station alpha_i(n_i) must equal set_load_dependence,
    which pins the fold point;
  - two stations sharing one unit of capacity by phi_s(n) = n_s/|n| conserve the
    population, and one link shared by S open routes is multiclass M/M/1-PS with
    exact means rho_s/(1-rho);
  - the same balanced model with Erlang service of equal mean must give the same
    means: INSENSITIVITY, the defining property of a Whittle network and the only
    oracle here that fails if PHASE transitions go unscaled.

The MATLAB twin is line-test.git/test_ctmc_global_dependence.m, the JAR twin is
SolverCTMCGlobalDependenceTest.java and the C++ twin is
cpp/tests/test_ctmc_global_dependence.cpp.
"""
import numpy as np
import pytest

from line_solver import (ClosedClass, Disabled, Erlang, Exp, Network, OpenClass,
                         Queue, SchedStrategy, Sink, Source, SolverCTMC, SolverMVA)


def closed_pair(erlang=False):
    """Closed cyclic PS pair, one class, 3 jobs."""
    model = Network('gdclosed')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(model, 'C1', 3, q1, 0)
    if erlang:
        q1.setService(c1, Erlang.fitMeanAndOrder(1.0, 3))
        q2.setService(c1, Erlang.fitMeanAndOrder(0.5, 3))
    else:
        q1.setService(c1, Exp(1.0))
        q2.setService(c1, Exp(2.0))
    model.link(Network.serialRouting(q1, q2))
    return model


def qlen(model, **kw):
    return np.asarray(SolverCTMC(model, **kw).getAvgQLen(), dtype=float)


def test_identity_scaling_changes_nothing():
    plain = qlen(closed_pair())
    m = closed_pair()
    m.set_global_dependence(lambda n: 1.0, 1.0)
    assert np.max(np.abs(qlen(m) - plain)) == 0.0


def test_reproduces_per_station_load_dependence():
    ld = closed_pair()
    ld.getNodeByName('Q2').setLoadDependence(np.array([1.0, 2.0, 2.0, 2.0]))
    a = qlen(ld)

    gd = closed_pair()

    def phi(n):
        v = np.ones(n.shape)
        nj = n[1, :].sum()
        if nj > 0:
            v[1, :] = min(nj, 2.0)
        return v

    gd.set_global_dependence(phi, 2.0)
    assert np.allclose(qlen(gd), a, atol=1e-9)


def shared_link_pair(erlang=False):
    """phi_s(n) = n_s/|n|: the two stations share one unit of capacity (balanced)."""
    model = closed_pair(erlang)

    def phi(n):
        v = np.ones(n.shape)
        tot = n.sum()
        if tot > 0:
            for i in range(n.shape[0]):
                v[i, :] = n[i, :].sum() / tot
        return v

    model.set_global_dependence(phi, 1.0)
    return model


def test_shared_link_conserves_population_and_is_insensitive():
    a = qlen(shared_link_pair())
    b = qlen(shared_link_pair(erlang=True))
    assert a.sum() == pytest.approx(3.0, abs=1e-9)
    # Insensitivity: an Erlang-3 of the same mean must not move the means
    assert np.allclose(b, a, rtol=1e-6)


def test_single_link_open_is_processor_sharing():
    """S open routes over one link is multiclass M/M/1-PS: E[n_s] = rho_s/(1-rho)."""
    nu = [0.20, 0.15]
    mu = [1.00, 0.50]
    C = 1.0
    model = Network('gdopen')
    src = Source(model, 'Src')
    q1 = Queue(model, 'R1', SchedStrategy.PS)
    q2 = Queue(model, 'R2', SchedStrategy.PS)
    snk = Sink(model, 'Snk')
    c1 = OpenClass(model, 'C1')
    c2 = OpenClass(model, 'C2')
    src.setArrival(c1, Exp(nu[0]))
    src.setArrival(c2, Exp(nu[1]))
    q1.setService(c1, Exp(mu[0]))
    q1.setService(c2, Disabled())
    q2.setService(c1, Disabled())
    q2.setService(c2, Exp(mu[1]))
    P = model.initRoutingMatrix()
    P.set(c1, c1, src, q1, 1.0)
    P.set(c1, c1, q1, snk, 1.0)
    P.set(c2, c2, src, q2, 1.0)
    P.set(c2, c2, q2, snk, 1.0)
    model.link(P)

    sn = model.getStruct()
    # getNodeIndex is 1-BASED while nodeToStation is a 0-based array, so the
    # conversion needs the -1; indexing without it silently reads another
    # station's row and leaves the model effectively unscaled
    i1 = int(sn.nodeToStation[model.getNodeIndex('R1') - 1])
    i2 = int(sn.nodeToStation[model.getNodeIndex('R2') - 1])

    def phi(n):
        v = np.ones(n.shape)
        tot = n[i1, 0] + n[i2, 1]
        if tot > 0:
            v[i1, 0] = n[i1, 0] * C / tot
            v[i2, 1] = n[i2, 1] * C / tot
        return v

    model.set_global_dependence(phi, C)
    QN = qlen(model, cutoff=16)
    rho = np.array(nu) / np.array(mu)
    exact = rho / (1 - rho.sum())
    got = np.array([QN[i1, 0], QN[i2, 1]])
    # the residual is the geometric tail leaving the truncation box, not the rate
    assert np.allclose(got, exact, rtol=1e-4)


def test_refused_by_product_form_solvers():
    m = closed_pair()
    m.set_global_dependence(lambda n: 1.0, 1.0)
    with pytest.raises(Exception):
        SolverMVA(m).getAvgQLen()


def test_malformed_declaration_is_refused():
    m = closed_pair()
    with pytest.raises(Exception):
        m.set_global_dependence(1.0, 1.0)
    with pytest.raises(Exception):
        m.set_global_dependence(lambda n: 1.0, None)
    with pytest.raises(Exception):
        m.set_global_dependence(lambda n: 1.0, -1.0)
    # wrong output shape, refused at declaration rather than mid-generation
    with pytest.raises(Exception):
        m.set_global_dependence(lambda n: np.ones((7, 3)), 1.0)


# --------------------------------------------------------------------------
# SSA
# --------------------------------------------------------------------------

def test_ssa_serial_matches_ctmc():
    """SolverSSA carries the same factorization on the sample path.

    phi(n) is a constant within a state, so it is evaluated once per state and
    multiplies every station service rate there. An unscaled run would report
    [1.5, 1.5] against the exact [2, 1].
    """
    from line_solver import SolverSSA
    m = shared_link_pair()
    exact = qlen(m)
    sim = np.asarray(SolverSSA(m, seed=23000, samples=200000).getAvgQLen(), dtype=float)
    assert np.allclose(sim, exact, rtol=3e-2)


def test_ssa_nrm_falls_back_to_serial():
    """The NRM builds its propensities from the per-station population slice and
    never sees the whole population matrix phi reads, so an explicit
    method='nrm' must divert to the serial engine rather than run unscaled."""
    from line_solver import SolverSSA
    m = shared_link_pair()
    exact = qlen(m)
    sim = np.asarray(
        SolverSSA(m, method='nrm', seed=23000, samples=200000).getAvgQLen(), dtype=float)
    assert np.allclose(sim, exact, rtol=3e-2)


# --------------------------------------------------------------------------
# JSON wire
# --------------------------------------------------------------------------

def test_json_round_trip_preserves_the_answer():
    """phi(n) is a handle, so it reaches the wire only as a TABLE: the writer
    materializes it over the lattice of the whole network state, restricted to
    the (station,class) slots a class can occupy. The oracle is the ANSWER, not
    the document: a table read back at the wrong coordinate still parses."""
    from line_solver.io.linemodel_io import save_model, load_model
    import json as _json
    import tempfile, os
    m = shared_link_pair()
    fd, path = tempfile.mkstemp(suffix='.json')
    os.close(fd)
    try:
        save_model(m, path)
        blk = _json.load(open(path))['model']['globalDependence']
        assert blk['type'] == 'globalDependent'
        assert len(blk['slots']) == 2
        assert len(blk['scaling']) == 16   # the (3+1)^2 box over the two slots
        back = load_model(path)
        assert back.get_global_dependence() is not None
        assert np.allclose(qlen(back), qlen(m), rtol=1e-10)
    finally:
        os.unlink(path)
