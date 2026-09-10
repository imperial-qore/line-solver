"""Regression tests for the shared fork-join driver and the SolverNC route on it.

The MMT transformation replaces every fork by a router and every join by a
zero-service delay, carrying the parallelism on auxiliary open classes. The
fixed point around it consumes only the metric matrices of an inner solve, so it
was moved out of SolverMVA into ``line_solver.solvers.fork_join_driver`` and is
now driven by both MVA and NC through ``_fj_inner_solver``.

The tests assert (i) that the move left SolverMVA numerically unchanged, (ii)
that SolverNC now solves fork-join models, and (iii) that the two routes agree.

Python mirror of line-test.git/test/testsFJ/test_fj_driver_nc.m and
jar/src/test/java/jline/solvers/fj/FJDriverNCTest.java.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Fork, Join, Network, OpenClass,
                         Queue, SchedStrategy, SolverCTMC, SolverMVA, SolverNC,
                         Sink, Source)


def _open_fj():
    model = Network('fjopen')
    source = Source(model, 'Source')
    fork = Fork(model, 'Fork')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    join = Join(model, 'Join', fork)
    sink = Sink(model, 'Sink')
    jobclass = OpenClass(model, 'C1')
    source.setArrival(jobclass, Exp(0.5))
    q1.setService(jobclass, Exp(2.0))
    q2.setService(jobclass, Exp(3.0))
    P = model.initRoutingMatrix()
    for a, b in ((source, fork), (fork, q1), (fork, q2),
                 (q1, join), (q2, join), (join, sink)):
        P.set(jobclass, jobclass, a, b, 1.0)
    model.link(P)
    return model


def _closed_fj(N):
    model = Network('fjclosed')
    delay = Delay(model, 'Think')
    fork = Fork(model, 'Fork')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    join = Join(model, 'Join', fork)
    jobclass = ClosedClass(model, 'C1', N, delay)
    delay.setService(jobclass, Exp(1.0))
    q1.setService(jobclass, Exp(2.0))
    q2.setService(jobclass, Exp(2.5))
    P = model.initRoutingMatrix()
    for a, b in ((delay, fork), (fork, q1), (fork, q2),
                 (q1, join), (q2, join), (join, delay)):
        P.set(jobclass, jobclass, a, b, 1.0)
    model.link(P)
    return model


def _plain_cqn():
    model = Network('cqn')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Q1', SchedStrategy.PS)
    jobclass = ClosedClass(model, 'C1', 3, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(2.0))
    P = model.initRoutingMatrix()
    P.set(jobclass, jobclass, Network.serialRouting(delay, queue))
    model.link(P)
    return model


def test_mva_forkjoin_values_unchanged():
    # Golden values measured against MATLAB, which they reproduce to 12 digits.
    #
    # Refreshed 2026-07-22 after the MVA method-dispatch fix (5c1a51fba): the
    # fixed point now converges further, so the open-model throughput obeys flow
    # balance to 1.2e-7 where the previous golden held it at 1.2e-4 (0.49987793
    # instead of the exact 0.5). MATLAB, native Python and lang="java" all agree
    # on the values below, so the old block was a stale golden, not a regression.
    Q = np.asarray(SolverMVA(_open_fj()).getAvgQLen()).ravel()
    T = np.asarray(SolverMVA(_open_fj()).getAvgTput()).ravel()
    assert Q == pytest.approx([0.0, 0.333332200208543, 0.199999847107973, 0.283332462218689], abs=1e-11)
    assert T == pytest.approx([0.5, 0.499999880790713, 0.499999880790713, 0.5], abs=1e-11)

    Q4 = np.asarray(SolverMVA(_closed_fj(4)).getAvgQLen()).ravel()
    T4 = np.asarray(SolverMVA(_closed_fj(4)).getAvgTput()).ravel()
    assert Q4 == pytest.approx([1.43153004164563, 2.18851387932703, 1.23958258137759, 1.84538508086124], abs=1e-10)
    assert T4 == pytest.approx([1.4315299485931, 1.43152994188856, 1.43152994188856, 1.4315299485931], abs=1e-10)


def test_driver_is_transparent_without_forks():
    # With no Fork in the model the driver must not engage at all
    Qm = np.asarray(SolverMVA(_plain_cqn()).getAvgQLen()).ravel()
    Qn = np.asarray(SolverNC(_plain_cqn()).getAvgQLen()).ravel()
    assert Qn == pytest.approx(Qm, rel=1e-6)


def test_nc_accepts_forkjoin():
    assert SolverNC(_open_fj()).supports(_open_fj())
    assert SolverNC(_closed_fj(3)).supports(_closed_fj(3))


def test_ctmc_rejects_class_switch_on_a_fork_output_edge():
    # A class switch declared ON a fork output edge is encoded by link() as an
    # auto-inserted ClassSwitch node, which lives in the sn index space but not
    # in the user-supplied routing matrix that fjtag writes sibling routing
    # into. It must be refused with the documented message, not with a raw
    # out-of-bounds error. Model of matlab/examples/basic/forkJoin/fj_cs_postfork.m.
    model = Network('cspostfork')
    delay = Delay(model, 'Delay')
    fork = Fork(model, 'Fork1')
    join = Join(model, 'Join1', fork)
    q1 = Queue(model, 'Queue1', SchedStrategy.PS)
    q2 = Queue(model, 'Queue2', SchedStrategy.PS)
    c1 = ClosedClass(model, 'class1', 1, delay)
    c2 = ClosedClass(model, 'class2', 1, delay)
    for c in (c1, c2):
        delay.setService(c, Exp(0.25))
        q1.setService(c, Exp(2.0))
        q2.setService(c, Exp(2.0))
    P = model.initRoutingMatrix()
    for a, b in ((delay, fork), (fork, q1), (fork, q2),
                 (q1, join), (q2, join), (join, delay)):
        P.set(c1, c1, a, b, 1.0)
    P.set(c2, c2, delay, fork, 1.0)
    P.set(c2, c1, fork, q1, 1.0)
    P.set(c2, c1, fork, q2, 1.0)
    model.link(P)
    with pytest.raises(RuntimeError, match='Class switching between fork and join'):
        SolverCTMC(model).getAvgQLen()


def test_nc_open_forkjoin_is_exact():
    # On an open transformed model NC applies exact open-network formulas, so it
    # reproduces the MATLAB reference.
    # Refreshed 2026-07-22 with the MVA golden above, same cause (5c1a51fba).
    Qn = np.asarray(SolverNC(_open_fj()).getAvgQLen()).ravel()
    Qm = np.asarray(SolverMVA(_open_fj()).getAvgQLen()).ravel()
    assert Qn == pytest.approx([0.0, 0.333319770103521, 0.1999926758261, 0.283319278568081], abs=1e-8)
    # NC runs its fork-join fixed point at iter_tol=1e-4 and MVA at 1e-6 (the
    # per-solver defaults, identical in MATLAB, the JAR and native Python), so
    # the two routes agree only to ~1.3e-5 on this model. Against the exact 1/3
    # at Q1, MVA errs 1.1e-6 and NC 1.35e-5; MATLAB shows the same spread. This
    # is a cross-solver agreement check, not a golden.
    assert Qn == pytest.approx(Qm, abs=1e-4)


def test_nc_closed_forkjoin_agrees_with_mva():
    # Both routes drive the same fixed point on the same transformation. The
    # auxiliary classes start at GlobalConstants.FineTol, so this is also the
    # test that the normalizing constant does not degenerate on them: a fixed
    # point stuck at its first iteration leaves the Join queue length at zero.
    Qn = np.asarray(SolverNC(_closed_fj(3)).getAvgQLen()).ravel()
    Qm = np.asarray(SolverMVA(_closed_fj(3)).getAvgQLen()).ravel()
    Tn = np.asarray(SolverNC(_closed_fj(3)).getAvgTput()).ravel()
    assert np.linalg.norm(Qn - Qm) / np.linalg.norm(Qm) < 0.05
    assert Qn[3] > 0.1
    assert Tn[0] == pytest.approx(Tn[3], rel=1e-6)


def test_ctmc_solves_closed_forkjoin_exactly():
    # SolverCTMC solves the same model natively on the tag-augmented copy, which
    # is the exact reference the two fixed-point routes approximate. The values
    # agree with MATLAB and the JAR to 12 significant digits.
    solver = SolverCTMC(_closed_fj(3))
    Qc = np.asarray(solver.getAvgQLen()).ravel()
    Tc = np.asarray(solver.getAvgTput()).ravel()
    assert Qc == pytest.approx([1.418899863876, 1.196465903699,
                                0.857919086961, 1.107815281588], abs=1e-9)
    # every station of this single-chain model sits on the same cycle, the Join
    # included: its DEP is the join firing, which fires only from the vanishing
    # marking in which the sibling set is complete and is therefore counted
    # through the rate complement rather than the tangible rows alone
    assert Tc == pytest.approx([Qc[0]] * 4, abs=1e-9)


def test_nc_closed_forkjoin_is_no_worse_than_mva_against_ctmc():
    # Mirrors the accuracy check of line-test.git/test/testsFJ/test_fj_driver_nc.m. Station 4
    # is the Join, whose queue length is the synchronisation delay of the
    # transformed formulation and is not comparable with the exact one.
    Qc = np.asarray(SolverCTMC(_closed_fj(3)).getAvgQLen()).ravel()
    Tc = np.asarray(SolverCTMC(_closed_fj(3)).getAvgTput()).ravel()
    Qm = np.asarray(SolverMVA(_closed_fj(3)).getAvgQLen()).ravel()
    Qn = np.asarray(SolverNC(_closed_fj(3)).getAvgQLen()).ravel()
    Tm = np.asarray(SolverMVA(_closed_fj(3)).getAvgTput()).ravel()
    Tn = np.asarray(SolverNC(_closed_fj(3)).getAvgTput()).ravel()

    err_mva = np.linalg.norm(Qm[:3] - Qc[:3]) / np.linalg.norm(Qc[:3])
    err_nc = np.linalg.norm(Qn[:3] - Qc[:3]) / np.linalg.norm(Qc[:3])
    assert err_nc < 0.20
    assert err_nc <= err_mva + 1e-6

    tput_err_mva = abs(Tm[0] - Tc[0]) / Tc[0]
    tput_err_nc = abs(Tn[0] - Tc[0]) / Tc[0]
    assert tput_err_nc < 0.15
    assert tput_err_nc <= tput_err_mva + 1e-6


def test_nc_forkjoin_metrics_are_finite():
    solver = SolverNC(_closed_fj(5))
    for metric in (solver.getAvgQLen(), solver.getAvgUtil(),
                   solver.getAvgRespT(), solver.getAvgTput()):
        arr = np.asarray(metric, dtype=float)
        assert np.all(np.isfinite(arr))
    assert np.all(np.asarray(solver.getAvgQLen(), dtype=float) >= 0)
