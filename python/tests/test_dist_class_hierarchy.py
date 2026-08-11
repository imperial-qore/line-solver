"""Class-hierarchy divergence checks for the Markovian distribution family.

In MATLAB and the JAR, PH, APH, Coxian, ME and RAP are direct siblings under
Markovian. In native Python, APH subclasses PH, Coxian subclasses PH and Cox2
subclasses Coxian, while ME and RAP are flat. That divergence is deliberate and
is not restructured here; these tests pin the places where it could otherwise
leak into observable behaviour, so that an isinstance-ordered dispatch cannot
silently degrade a Coxian to a PH or a Cox2 to a Coxian again.

Ground truth for the process-type names is MATLAB (verified by running
model.getStruct() there) for PH, APH, Coxian, ME and RAP, and the JAR
(jline.lang.Network.getProcessType) for Cox2, which MATLAB has no instance of
because MATLAB's Cox2 is a sealed factory returning Coxian objects.
"""

import os
import tempfile

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp,
                         SchedStrategy, SolverFluid, SolverSSA)
from line_solver.constants import ProcessType
from line_solver.distributions.markovian import ME, RAP, PH, APH, Coxian, Cox2
from line_solver.io.linemodel_io import save_model, load_model


# Subgenerator shared by the PH and APH cases: two phases, rates 3 and 4.
_T = np.array([[-3.0, 3.0], [0.0, -4.0]])

# A genuine ME whose alpha is a proper probability vector, so the density stays
# nonnegative and the constructor issues no warning.
_ME_ALPHA = np.array([0.4, 0.6])
_ME_A = np.array([[-2.0, 1.0], [0.0, -3.0]])

# A genuine RAP: H1 carries a negative off-diagonal entry, so it is not a MAP.
_RAP_H0 = np.array([[-5.0, 1.0], [1.0, -6.0]])
_RAP_H1 = np.array([[4.5, -0.5], [4.5, 0.5]])


def _service_model(dist):
    """Open model whose single queue serves with the given distribution."""
    model = Network('svc')
    source = Source(model, 'S')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    sink = Sink(model, 'K')
    oclass = OpenClass(model, 'C1')
    source.setArrival(oclass, Exp(0.1))
    queue.setService(oclass, dist)
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _arrival_model(dist):
    """Open model whose source emits with the given distribution."""
    model = Network('arv')
    source = Source(model, 'S')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    sink = Sink(model, 'K')
    oclass = OpenClass(model, 'C1')
    source.setArrival(oclass, dist)
    queue.setService(oclass, Exp(6.0))
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _service_dists():
    return {
        'PH': PH(np.array([1.0, 0.0]), _T),
        'APH': APH(np.array([1.0, 0.0]), _T),
        'Coxian': Coxian(np.array([3.0, 4.0]), np.array([0.4, 1.0])),
        'Cox2': Cox2(3.0, 4.0, 0.4),
        'ME': ME(_ME_ALPHA, _ME_A),
    }


@pytest.mark.parametrize('name,expected', [
    ('PH', 'PH'),
    ('APH', 'APH'),
    ('Coxian', 'COXIAN'),
    ('Cox2', 'COX2'),
    ('ME', 'ME'),
])
def test_service_procid_matches_concrete_class(name, expected):
    """sn.procid names the concrete class, not the Python base class.

    Compared BY NAME: the ProcessType ordinals differ across codebases.
    Coxian and Cox2 subclass PH here, so an isinstance test against PH placed
    ahead of them reports ProcessType.PH for all three.
    """
    dist = _service_dists()[name]
    sn = _service_model(dist).getStruct()
    assert sn.procid[1, 0].name == expected


def test_rap_arrival_procid_is_rap():
    sn = _arrival_model(RAP(_RAP_H0, _RAP_H1)).getStruct()
    assert sn.procid[0, 0].name == 'RAP'


@pytest.mark.parametrize('name', ['PH', 'APH', 'Coxian', 'Cox2', 'ME'])
def test_service_proc_holds_two_matrices(name):
    """Every member of the family stores a two-entry process representation."""
    dist = _service_dists()[name]
    sn = _service_model(dist).getStruct()
    proc = sn.proc[1][0]
    assert proc is not None
    assert len(proc) == 2
    assert np.all(np.isfinite(np.asarray(proc[0], dtype=float)))
    assert np.all(np.isfinite(np.asarray(proc[1], dtype=float)))


@pytest.mark.parametrize('name', ['PH', 'APH', 'Coxian', 'ME'])
def test_json_round_trip_preserves_class_and_moments(name):
    """A save/load cycle must not downgrade the concrete distribution class.

    Cox2 is excluded deliberately: linemodel_save.m writes both Coxian and Cox2
    as type 'Coxian' and jline.io.LineModelIO accepts 'Cox2' only as an alias of
    it, so the wire format has no Cox2 identity to preserve. That case is
    covered by test_cox2_serializes_as_coxian.
    """
    dist = _service_dists()[name]
    model = _service_model(dist)
    path = tempfile.mktemp(suffix='.json')
    try:
        save_model(model, path)
        reloaded = load_model(path)
    finally:
        if os.path.exists(path):
            os.unlink(path)

    back = reloaded.getNodes()[1].getService(reloaded.getClasses()[0])
    assert type(back).__name__ == type(dist).__name__
    assert back.getMean() == pytest.approx(dist.getMean(), rel=1e-12)
    assert back.getSCV() == pytest.approx(dist.getSCV(), rel=1e-12)
    assert (reloaded.getStruct().procid[1, 0].name
            == model.getStruct().procid[1, 0].name)


def test_rap_arrival_json_round_trip_preserves_matrices():
    dist = RAP(_RAP_H0, _RAP_H1)
    model = _arrival_model(dist)
    path = tempfile.mktemp(suffix='.json')
    try:
        save_model(model, path)
        reloaded = load_model(path)
    finally:
        if os.path.exists(path):
            os.unlink(path)

    back = reloaded.getNodes()[0].getArrival(reloaded.getClasses()[0])
    assert type(back).__name__ == 'RAP'
    assert np.allclose(back.getD0(), dist.getD0())
    assert np.allclose(back.getD1(), dist.getD1())
    assert reloaded.getStruct().procid[0, 0].name == 'RAP'


def test_cox2_serializes_as_coxian():
    """Cox2 goes on the wire as a Coxian, matching linemodel_save.m.

    The moments must survive the alias, which is what makes the downgrade
    lossless rather than silent.
    """
    dist = Cox2(3.0, 4.0, 0.4)
    model = _service_model(dist)
    path = tempfile.mktemp(suffix='.json')
    try:
        save_model(model, path)
        reloaded = load_model(path)
    finally:
        if os.path.exists(path):
            os.unlink(path)

    back = reloaded.getNodes()[1].getService(reloaded.getClasses()[0])
    assert type(back).__name__ == 'Coxian'
    assert back.getMean() == pytest.approx(dist.getMean(), rel=1e-12)
    assert back.getSCV() == pytest.approx(dist.getSCV(), rel=1e-12)


@pytest.mark.parametrize('solver', [SolverFluid, SolverSSA])
def test_cox2_and_equivalent_coxian_agree(solver):
    """Relabelling Cox2 as COX2 must not change any solver answer.

    Cox2(mu1,mu2,phi1) is by construction the same law as
    Coxian([mu1,mu2],[phi1,1]), so every metric has to coincide. SSA is seeded
    so the two runs share a sample path.
    """
    cox2 = Cox2(3.0, 4.0, 0.4)
    coxian = Coxian(np.array([3.0, 4.0]), np.array([0.4, 1.0]))
    cols = ['QLen', 'Util', 'RespT']

    def solve(dist):
        model = _service_model(dist)
        if solver is SolverSSA:
            return solver(model, seed=23000, samples=20000).getAvgTable()[cols].to_numpy()
        return solver(model).getAvgTable()[cols].to_numpy()

    assert np.allclose(solve(cox2), solve(coxian))


def test_me_phase_decomposition_matches_matlab():
    """ME.getMu/getPhi/getInitProb reproduce the MATLAB Markovian formulas.

    Reference values were obtained by running MATLAB on the same
    representation: mu = -diag(A), phi = -D1*e ./ diag(A), alpha as given.
    getPhi leaves [0,1] here, which is the expected consequence of a
    matrix-exponential having no phase-type decomposition; it is recorded, not
    repaired.
    """
    me = ME(np.array([1.0, 0.0]), np.array([[-1.0, 2.0], [0.0, -3.0]]))
    assert np.allclose(me.getMu(), np.array([1.0, 3.0]))
    assert np.allclose(me.getPhi(), np.array([-1.0, 1.0]))
    assert np.allclose(me.getInitProb(), np.array([1.0, 0.0]))


def test_rap_init_prob_matches_matlab():
    """RAP.getInitProb is the arrival-embedded vector pie, as in MATLAB.

    Reference from MATLAB map_pie on the same (H0,H1). Its second entry is
    negative, which is exactly why a RAP has no phase-type initial
    distribution.
    """
    rap = RAP(_RAP_H0, _RAP_H1)
    assert np.allclose(rap.getInitProb(),
                       np.array([1.1020408163265307, -0.10204081632653073]))
