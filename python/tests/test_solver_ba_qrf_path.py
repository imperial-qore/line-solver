"""End-to-end coverage of the SolverBA -> solver_ctmc_qrf_analyzer -> qrf_noblo_*
adapter path on a native Network.

This is the path that kept the NLP QRF method names withheld from BA_METHODS. Exercising
it found three defects, all now fixed in both codebases:
  1. the adapter wrote v transposed relative to mu, which reversed the phase order
     of an Erlang and returned a degenerate UN = [0, 0] at K > 1;
  2. throughput was inverted as X = QN[refstat]/(V*stime_ref), valid only when the
     reference station is an infinite server, so a FCFS reference station produced
     X above the bottleneck capacity and UN > 1;
  3. an infinite-server station was silently modelled as a single-server queue.

WHAT THE ORACLE IS, AND WHAT IT IS NOT
--------------------------------------
qrf_noblo_mmi returns a POINT in the QRF polytope, chosen by a non-convex MMI
objective. It is not an optimized face, so its value carries NO guaranteed
direction relative to the exact solution. Measured on cyclic chains: at M=3 N=3 it
lands below the exact utilization, at M=3 N=4 and M=4 N=3 above it. An earlier
version of this file asserted UN <= exact, which held only on the M=2 fixtures and
was false in general.

The invariant that does hold is CONTAINMENT in the LP range of the same polytope.
REFERENCE below is that range, per station, computed with glpsol on the paper's
AMPL no-blocking model (noblo_skel.mod, the no-blocking specialisation of
qrboundsbas_skel.mod) for each fixture here.

Where the LP range is degenerate (min == max) the polytope is tight and the QRF
answer must equal the exact CTMC value. That is a property of the instance, not a
lucky fixture: at M == 2 the pairwise joint of a closed chain is fully determined
by the marginal, and it also happens at M = 3, N = 2. Elsewhere the polytope is a
strict relaxation and only containment is required.
"""
import numpy as np
import pytest

from line_solver import (Network, Queue, Delay, ClosedClass, SchedStrategy, Exp,
                         Erlang, SolverCTMC)
from line_solver.solvers.solver_ba.solver_ba import BA_METHODS
from line_solver.api.solvers.ctmc.solver_ctmc_qrf_analyzer import (
    solver_ctmc_qrf_analyzer)

# glpsol [min, max] of U[station] on noblo_skel.mod, per fixture.
REFERENCE = {
    'm2k1_N2': [(0.789473684, 0.789473684), (0.526315790, 0.526315790)],
    'm2k1_N3': [(0.876923077, 0.876923077), (0.584615385, 0.584615385)],
    'm2k2_N2': [(0.705882353, 0.857142857), (0.470588235, 0.571428571)],
    'm2k2_N3': [(0.735849057, 0.951219512), (0.490566038, 0.634146342)],
    'm3k1_N2': [(0.678260870, 0.678260870), (0.452173913, 0.452173913),
                (0.339130435, 0.339130435)],
    'm3k1_N3': [(0.775939850, 0.824385805), (0.517293233, 0.549590537),
                (0.387969925, 0.412192903)],
    'm3k1_N4': [(0.839042176, 0.894297964), (0.559361451, 0.596198643),
                (0.419521088, 0.447148982)],
    'm3mix_N3': [(0.677878517, 0.870629371), (0.451919012, 0.580419580),
                 (0.338939259, 0.435314685)],
}

# fixture -> (N, [(rate, phases), ...]). Service time of station i is 1/rate.
FIXTURES = {
    'm2k1_N2': (2, [(1.0, 1), (1.5, 1)]),
    'm2k1_N3': (3, [(1.0, 1), (1.5, 1)]),
    'm2k2_N2': (2, [(1.0, 1), (1.5, 2)]),
    'm2k2_N3': (3, [(1.0, 1), (1.5, 2)]),
    'm3k1_N2': (2, [(1.0, 1), (1.5, 1), (2.0, 1)]),
    'm3k1_N3': (3, [(1.0, 1), (1.5, 1), (2.0, 1)]),
    'm3k1_N4': (4, [(1.0, 1), (1.5, 1), (2.0, 1)]),
    'm3mix_N3': (3, [(1.0, 1), (1.5, 2), (2.0, 1)]),
}


def _chain(name):
    """Closed cyclic chain of single-server FCFS queues, Exp or Erlang service."""
    N, spec = FIXTURES[name]
    model = Network(name)
    queues = [Queue(model, 'Q%d' % (i + 1), SchedStrategy.FCFS)
              for i in range(len(spec))]
    cls = ClosedClass(model, 'C1', N, queues[0])
    for q, (rate, phases) in zip(queues, spec):
        q.setService(cls, Exp(rate) if phases == 1
                     else Erlang(rate * phases, phases))
    model.link(model.serialRouting(*queues))
    return model


def _run(name, method='qrf.mmi'):
    QN, UN, RN, TN, CN, XN, _rt = solver_ctmc_qrf_analyzer(
        _chain(name).getStruct(), _Opts(method))
    return (np.asarray(QN, dtype=float).flatten(),
            np.asarray(UN, dtype=float).flatten(),
            float(np.asarray(XN, dtype=float).flatten()[0]))


class _Opts(object):
    def __init__(self, method):
        self.method = method
        self.config = {}


def test_nlp_tokens_are_advertised():
    for tok in ('qr', 'qrf.mmi', 'qrf.mem', 'qrf.mmi.ld', 'qrf.mmi.linear'):
        assert tok in BA_METHODS


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_within_glpsol_bound_range(name):
    """The only invariant the MMI point guarantees: it lies in the LP range of
    the polytope it was drawn from."""
    _QN, UN, _XN = _run(name)
    for i, (lo, hi) in enumerate(REFERENCE[name]):
        assert lo - 1e-6 <= UN[i] <= hi + 1e-6, (
            '%s station %d: %.9f outside glpsol [%.9f, %.9f]'
            % (name, i + 1, UN[i], lo, hi))


# The stations whose glpsol range is degenerate, per fixture: there, and only
# there, the polytope is tight and the QRF point must equal the exact value.
TIGHT = dict((name, [i for i, (lo, hi) in enumerate(rng) if hi - lo < 1e-9])
             for name, rng in REFERENCE.items())


@pytest.mark.parametrize('name', sorted(n for n in FIXTURES if TIGHT[n]))
def test_exact_where_polytope_is_tight(name):
    """Where the glpsol range collapses to a point the relaxation is exact, so
    the adapter must reproduce the CTMC. This is what catches a marshalling
    error such as the transposed v.

    ONLY THE TIGHT FIXTURES ARE PARAMETRIZED. On the others every station is a
    strict relaxation, so there is no station this test could assert on -- it
    used to be generated for them and skip, which costs a line of the summary
    and says nothing. Containment is what those fixtures assert, in
    `test_within_glpsol_bound_range`, and every fixture is covered there.
    """
    tight = TIGHT[name]
    _QN, UN, _XN = _run(name)
    exact = np.asarray(SolverCTMC(_chain(name)).getAvgTable().Util,
                       dtype=float).flatten()
    for i in tight:
        assert UN[i] == pytest.approx(exact[i], rel=1e-7, abs=1e-9)


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_utilization_is_a_probability(name):
    """The regression that motivated the fix: UN must never exceed 1 at a
    finite-server station, and the population must be conserved."""
    QN, UN, XN = _run(name)
    N, _spec = FIXTURES[name]
    assert np.all(UN >= -1e-12)
    assert np.all(UN <= 1.0 + 1e-9), 'UN above 1 for %s' % name
    assert QN.sum() == pytest.approx(float(N), rel=1e-9)
    assert XN > 0.0


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_throughput_is_consistent_across_stations(name):
    """X is derived from one finite-server station via U_i = X*V_i*stime_i. That
    is only legitimate if every station agrees, which is the property the old
    QN[refstat] inversion lacked."""
    _QN, UN, XN = _run(name)
    _N, spec = FIXTURES[name]
    stimes = np.array([1.0 / rate for rate, _ph in spec])
    np.testing.assert_allclose(UN / stimes, np.full(len(spec), XN), rtol=1e-8)


@pytest.mark.parametrize('name', ['m2k2_N3', 'm3mix_N3'])
def test_phasetype_is_not_degenerate(name):
    """The transposed-v symptom was UN = [0, 0] at K > 1, which reads as an NLP
    failure rather than a marshalling bug."""
    _QN, UN, _XN = _run(name)
    assert np.all(UN > 1e-6), 'degenerate utilization for %s' % name


def _delay_model(N=3):
    model = Network('qrf_delay')
    think = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.FCFS)
    cls = ClosedClass(model, 'C1', N, think)
    think.setService(cls, Exp(1.0))
    q.setService(cls, Exp(1.5))
    model.link(model.serialRouting(think, q))
    return model


def test_infinite_server_is_rejected_by_the_alpha_free_arms():
    """The alpha-free arms build a population-free q, so they model every
    station as a single server and a Delay would silently change the model
    rather than approximate it. The refusal now names the arms that serve it."""
    with pytest.raises(ValueError, match="Use 'qrf.mmi.ld' or 'qrf.mmi.linear'"):
        solver_ctmc_qrf_analyzer(_delay_model().getStruct(), _Opts('qrf.mmi'))


def test_infinite_server_is_exact_on_the_ld_arm():
    """alpha(i,n) = n IS the rate law of an infinite server, so the ld arm runs
    the model's own chain. At M = 2 the pairwise joint of a closed chain is
    fully determined by the marginal, so the polytope is tight and the answer
    must be the exact CTMC one, not merely close to it."""
    model = _delay_model()
    QN, UN, _RN, TN, _CN, _XN, _rt = solver_ctmc_qrf_analyzer(
        model.getStruct(), _Opts('qrf.mmi.ld'))
    exact = SolverCTMC(model).getAvgTable()
    assert np.allclose(np.asarray(QN).ravel(), np.asarray(exact.QLen, dtype=float), atol=1e-9)
    assert np.allclose(np.asarray(UN).ravel(), np.asarray(exact.Util, dtype=float), atol=1e-9)
    assert np.allclose(np.asarray(TN).ravel(), np.asarray(exact.Tput, dtype=float), atol=1e-9)
