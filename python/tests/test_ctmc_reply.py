"""Synchronous calls (REPLY signals) under SolverCTMC.

A class whose jobs make a synchronous call (sn.syncreply >= 0) leaves the
caller for the callee but KEEPS its server there; the server is released only
when the matching REPLY signal class arrives back. LDES keys that hold by job
id (Solver_ssj.pendingReplyMap); a CTMC has no job identity, so the held
servers are counted per class in a local-state block (api.state.reply_block)
and the caller is selected structurally, as the station the REPLY class is
routed into.

Model (the canonical LQN shape RefTask -> ClientTask --synchCall--> ServerTask):

  ThinkDelay -> ClientQueue -> ServerQueue -> [switch to Reply]
       ^             ^                              |
       |             |______ Reply signal __________|
       |____________ [switch back to Job] ___________|

Conventions pinned here, all taken from LDES:
  - QLen and Util at the caller COUNT the held server and the job blocked out
    at the callee (simultaneous resource possession, the point of the feature);
    RespT does NOT, since it measures time spent at the station.
  - The REPLY does not queue at the caller: it takes the server it released and
    is routed on, so its queue length there is only the residual of its own
    Immediate service (Tput / GlobalConstants.Immediate).
  - Queue lengths therefore sum to N PLUS the mean number of outstanding calls
    (a blocked job is counted at the callee and at the caller alike); the token
    count is still N, since the reply IS the calling job.

The expected values are the MATLAB SolverCTMC goldens, which agree with LDES at
3e5 samples to within 0.41%. Companion of
line-test.git/test/testsAdvFeatures/des/test_ctmc_reply.m.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, ClosedSignal,
                         SignalType, SchedStrategy, Exp, Immediate, SolverCTMC,
                         VerboseLevel)
from line_solver.constants import GlobalConstants


def reply_model(N, muC, muS, Z, c):
    model = Network('ClosedReplySignal')
    delay = Delay(model, 'ThinkDelay')
    queue1 = Queue(model, 'ClientQueue', SchedStrategy.FCFS)
    queue1.setNumberOfServers(c)
    queue2 = Queue(model, 'ServerQueue', SchedStrategy.FCFS)

    jobClass = ClosedClass(model, 'Job', N, delay)
    delay.setService(jobClass, Exp(1.0 / Z))
    queue1.setService(jobClass, Exp(muC))
    queue2.setService(jobClass, Exp(muS))

    # The reply is the calling job, class-switched at the callee, so the
    # population is conserved; ClosedSignal shares the caller's reference station.
    replyClass = ClosedSignal(model, 'Reply', SignalType.REPLY, delay).forJobClass(jobClass)
    delay.setService(replyClass, Immediate())
    queue1.setService(replyClass, Immediate())
    queue2.setService(replyClass, Immediate())

    P = model.initRoutingMatrix()
    P[jobClass, jobClass] = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]], dtype=float)
    P[jobClass, replyClass] = np.array([[0, 0, 0], [0, 0, 0], [0, 1, 0]], dtype=float)
    P[replyClass, jobClass] = np.array([[0, 0, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
    model.link(P)
    return model


# (N, muC, muS, Z, c) -> MATLAB SolverCTMC goldens
#   X, client Util, client QLen, client RespT, server QLen
CONFIGS = [
    ((2, 4, 6, 1.0, 1), (1.316614409, 0.5485893369, 0.6833855782, 0.3523809555, 0.2194357348)),
    ((3, 4, 6, 1.0, 1), (1.795275570, 0.7480314876, 1.2047244120, 0.5043859739, 0.2992125950)),
    ((2, 2, 3, 0.5, 1), (1.100917421, 0.9174311841, 1.4495412790, 0.9833333403, 0.3669724736)),
    ((4, 5, 5, 1.0, 1), (2.179233428, 0.8716933713, 1.8207665500, 0.6355078104, 0.4358466856)),
    ((3, 4, 6, 1.0, 2), (2.038178615, 0.4594935554, 0.9618213644, 0.2710159469, 0.4094424571)),
]

# Stations and classes, in model declaration order.
DELAY, CLIENT, SERVER = 0, 1, 2
JOB, REPLY = 0, 1


def _solve(cfg):
    solver = SolverCTMC(reply_model(*cfg), verbose=VerboseLevel.SILENT)
    return (solver.getAvgQLen(), solver.getAvgUtil(),
            solver.getAvgRespT(), solver.getAvgTput())


@pytest.mark.parametrize('cfg,gold', CONFIGS)
def test_reply_matches_matlab_goldens(cfg, gold):
    """Throughput and the caller's QLen/Util/RespT against the MATLAB goldens."""
    Q, U, R, T = _solve(cfg)
    X, Uc, Qc, Rc, Qs = gold
    assert T[CLIENT, JOB] == pytest.approx(X, rel=1e-6)
    assert U[CLIENT, JOB] == pytest.approx(Uc, rel=1e-6)
    assert Q[CLIENT, JOB] == pytest.approx(Qc, rel=1e-6)
    assert R[CLIENT, JOB] == pytest.approx(Rc, rel=1e-6)
    assert Q[SERVER, JOB] == pytest.approx(Qs, rel=1e-6)


@pytest.mark.parametrize('cfg,gold', CONFIGS)
def test_reply_does_not_reside_at_the_caller(cfg, gold):
    """The reply takes the server it released and is routed on.

    Queueing it instead stole service capacity from the caller and cost about
    5% of the throughput. All that is left is the residual of the reply's own
    Immediate service, X / GlobalConstants.Immediate, which getAvg reports as
    zero: a response time below 10*FineTol masks queue length and utilization.
    Same assertion as the MATLAB and JAR twins.
    """
    Q, U, R, T = _solve(cfg)
    # RespT is QN/TN, two independently accumulated sums, so the bound carries a
    # relative slack: lang='java' lands one ulp above it, native exactly on it.
    assert R[CLIENT, REPLY] <= (1.0 + 1e-9) / GlobalConstants.Immediate
    assert Q[CLIENT, REPLY] < 1e-8
    assert U[CLIENT, REPLY] < 1e-8


@pytest.mark.parametrize('cfg,gold', CONFIGS)
def test_reply_throughput_is_uniform(cfg, gold):
    """One token per call, so throughput is uniform around the cycle."""
    _, _, _, T = _solve(cfg)
    tput = np.array([T[DELAY, JOB], T[CLIENT, JOB], T[CLIENT, REPLY], T[SERVER, JOB]])
    assert np.max(np.abs(tput - gold[0])) < 1e-6


@pytest.mark.parametrize('cfg,gold', CONFIGS)
def test_qlen_excess_equals_outstanding_calls(cfg, gold):
    """Simultaneous resource possession double-counts by construction.

    A job blocked out at the callee is counted BOTH there and at the caller
    that holds a server for it, so the station queue lengths sum to N plus the
    mean number of outstanding calls -- which here is exactly the callee's
    queue length (every call in progress sits at the ServerQueue). LDES sums
    the same way. The underlying token count is still N: the reply IS the
    calling job, class-switched.
    """
    Q, _, _, _ = _solve(cfg)
    N = cfg[0]
    excess = float(np.sum(Q)) - N
    assert excess == pytest.approx(Q[SERVER, JOB], abs=1e-6)


@pytest.mark.parametrize('cfg,gold', CONFIGS)
def test_util_counts_the_held_server(cfg, gold):
    """Utilization exceeds the carried load T*E[S] by the blocked holding time."""
    _, U, _, T = _solve(cfg)
    N, muC, muS, Z, c = cfg
    carried = T[CLIENT, JOB] / muC / c
    assert U[CLIENT, JOB] > carried + 0.05


def test_reply_block_is_declared_at_the_caller_only():
    """The holding station is the one the REPLY class is routed INTO.

    In this shape that is the client and NOT the server: the server's departure
    is the one that CREATES the reply, and LDES likewise does not block it
    (Solver_ssj: !classSwitchedToReply). The Delay is skipped because an INF
    station has a server for every job.
    """
    sn = reply_model(2, 4, 6, 1.0, 1).getStruct()
    assert list(np.asarray(sn.syncreply).ravel()) == [REPLY, -1]
    rb = np.asarray(sn.replyblock)
    client_node = int(sn.stationToNode[CLIENT])
    assert rb[client_node, JOB] == 1
    assert int(np.sum(rb)) == 1


def test_non_fcfs_caller_is_rejected():
    """Holding a server across a call has no state representation elsewhere."""
    model = Network('PSReply')
    delay = Delay(model, 'ThinkDelay')
    queue1 = Queue(model, 'ClientQueue', SchedStrategy.PS)
    queue2 = Queue(model, 'ServerQueue', SchedStrategy.FCFS)
    jobClass = ClosedClass(model, 'Job', 2, delay)
    delay.setService(jobClass, Exp(1.0))
    queue1.setService(jobClass, Exp(4.0))
    queue2.setService(jobClass, Exp(6.0))
    replyClass = ClosedSignal(model, 'Reply', SignalType.REPLY, delay).forJobClass(jobClass)
    delay.setService(replyClass, Immediate())
    queue1.setService(replyClass, Immediate())
    queue2.setService(replyClass, Immediate())
    P = model.initRoutingMatrix()
    P[jobClass, jobClass] = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]], dtype=float)
    P[jobClass, replyClass] = np.array([[0, 0, 0], [0, 0, 0], [0, 1, 0]], dtype=float)
    P[replyClass, jobClass] = np.array([[0, 0, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
    model.link(P)
    with pytest.raises(RuntimeError, match='REPLY signals'):
        model.getStruct()


def test_reply_block_stays_at_the_tail_under_a_map_service():
    """A persistent MAP server-phase variable must not displace the reply block.

    Native Python appends the MAP phase variable AFTER the state space is
    generated, i.e. behind the counter columns fromMarginal produced, so
    ctmc_ssg rotates the block back to the tail. MMPP2(4, 4, 1, 1) has equal
    phase rates, hence is Exp(4): the model must reproduce the Exp(4) goldens
    of the first configuration exactly.
    """
    from line_solver import MMPP2
    N, muC, muS, Z, c = CONFIGS[0][0]
    X, Uc, Qc, Rc, Qs = CONFIGS[0][1]

    model = Network('MapReply')
    delay = Delay(model, 'ThinkDelay')
    queue1 = Queue(model, 'ClientQueue', SchedStrategy.FCFS)
    queue2 = Queue(model, 'ServerQueue', SchedStrategy.FCFS)
    jobClass = ClosedClass(model, 'Job', N, delay)
    delay.setService(jobClass, Exp(1.0 / Z))
    queue1.setService(jobClass, MMPP2(muC, muC, 1.0, 1.0))
    queue2.setService(jobClass, Exp(muS))
    replyClass = ClosedSignal(model, 'Reply', SignalType.REPLY, delay).forJobClass(jobClass)
    delay.setService(replyClass, Immediate())
    queue1.setService(replyClass, Immediate())
    queue2.setService(replyClass, Immediate())
    P = model.initRoutingMatrix()
    P[jobClass, jobClass] = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]], dtype=float)
    P[jobClass, replyClass] = np.array([[0, 0, 0], [0, 0, 0], [0, 1, 0]], dtype=float)
    P[replyClass, jobClass] = np.array([[0, 0, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
    model.link(P)

    solver = SolverCTMC(model, verbose=VerboseLevel.SILENT)
    Q, U, R, T = (solver.getAvgQLen(), solver.getAvgUtil(),
                  solver.getAvgRespT(), solver.getAvgTput())
    assert T[CLIENT, JOB] == pytest.approx(X, rel=1e-6)
    assert U[CLIENT, JOB] == pytest.approx(Uc, rel=1e-6)
    assert Q[CLIENT, JOB] == pytest.approx(Qc, rel=1e-6)
    assert R[CLIENT, JOB] == pytest.approx(Rc, rel=1e-6)
    assert Q[SERVER, JOB] == pytest.approx(Qs, rel=1e-6)
