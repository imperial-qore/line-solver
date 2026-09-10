"""Tests for Fork.setTasksPerLink under the MMT fork-join transform.

A fork with tasksPerLink = w sends w IDENTICAL tasks down each of its B links, so
a firing creates w*B siblings. Two things follow, and both were wrong before:

  - the join synchronises on the order statistic of w*B branch times, each branch
    REPLICATED w times, not on w times the order statistic of B. The latter is
    w*H_B/mu where the answer is H_(w*B)/mu -- already 3.0/mu against 2.083/mu at
    B = w = 2 -- and it cost SolverMVA about ten points of accuracy against the
    simulators;
  - a branch station can hold w jobs per circulating parent, so the chain
    population is no longer its buffer bound.

The reference values below are SolverMVA on the model in `_tpl_model`; MATLAB, the
JAR, native Python and the C++ port agree on them to every printed digit, and
SolverJMT/SolverLDES (which simulate tasksPerLink directly and agree with each
other to 0.3%) put the true throughput about 10% above them, the same bias an
ordinary fork-join carries under this transform.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Fork, Join, Network, Queue,
                         SchedStrategy, SolverMVA)


def _tpl_model(N, w, nbranch):
    model = Network('FJ-TPL')
    delay = Delay(model, 'Delay')
    qs = [Queue(model, 'Queue%d' % (b + 1), SchedStrategy.PS) for b in range(nbranch)]
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    c = ClosedClass(model, 'class1', N, delay)
    delay.setService(c, Exp(1.0))
    for q in qs:
        q.setService(c, Exp(2.0))
    fork.setTasksPerLink(w)
    P = model.init_routing_matrix()
    P.set(c, c, delay, fork, 1.0)
    for q in qs:
        P.set(c, c, fork, q, 1.0)
        P.set(c, c, q, join, 1.0)
    P.set(c, c, join, delay, 1.0)
    model.link(P)
    return model


def _tput(model):
    return float(np.asarray(SolverMVA(model).getAvg()[3])[0, 0])


# (branches, w) -> SolverMVA throughput, agreed across MATLAB / JAR / python / C++
REFERENCE = {
    (2, 1): 1.167124,
    (2, 2): 0.702141,
    (2, 3): 0.507307,
    (3, 1): 1.058612,
    (3, 2): 0.662789,
    (3, 3): 0.487276,
}


@pytest.mark.parametrize('key', sorted(REFERENCE))
def test_tasks_per_link_throughput(key):
    nbranch, w = key
    assert _tput(_tpl_model(3, w, nbranch)) == pytest.approx(REFERENCE[key], rel=1e-5)


def test_one_task_per_link_is_the_ordinary_fork_join():
    """w = 1 replicates a branch to itself and scales the capacity by one, so a
    plain fork must answer exactly as it did before tasksPerLink was handled."""
    for nbranch in (2, 3):
        model = _tpl_model(3, 1, nbranch)
        assert _tput(model) == pytest.approx(REFERENCE[(nbranch, 1)], rel=1e-5)
        assert np.asarray(model.getStruct().classcap)[0, 0] == pytest.approx(3.0)


def test_more_tasks_per_link_lowers_throughput():
    """Each extra task per link is extra work per firing AND a later synchronisation."""
    for nbranch in (2, 3):
        xs = [_tput(_tpl_model(3, w, nbranch)) for w in (1, 2, 3)]
        assert xs[0] > xs[1] > xs[2]


def test_branch_capacity_scales_with_tasks_per_link():
    """A branch station holds w jobs per circulating parent, not one."""
    for nbranch in (2, 3):
        for w in (1, 2, 3):
            cap = np.asarray(_tpl_model(3, w, nbranch).getStruct().classcap)
            assert cap[0, 0] == pytest.approx(3 * w)


def test_sync_delay_is_not_the_scaled_order_statistic():
    """Pin the formula itself: w*E[X_(B)] would be strictly larger than
    E[X_(w*B)] on identical branches, so the two cannot both be right."""
    from line_solver.api.fjnative import fj_ordstat_exp

    ri = np.array([0.5, 0.5])
    scaled = 2 * fj_ordstat_exp(ri, 2)                     # the old rule
    replicated = fj_ordstat_exp(np.tile(ri, 2), 4)         # the rule in force
    assert scaled == pytest.approx(2 * 0.5 * (1 + 1.0 / 2), rel=1e-12)
    assert replicated == pytest.approx(0.5 * (1 + 1.0 / 2 + 1.0 / 3 + 1.0 / 4), rel=1e-12)
    assert replicated < scaled


def test_sibling_count_is_per_link_not_out_degree_times_a_scalar():
    """A VARIABLE FORKING LEVEL sets one link of one class, so the sibling count
    is the sum of the per-link counts and not the out-degree times a node-wide
    scalar. sn_join_siblings reads fanOutLink, which carries them."""
    from line_solver.api.fjnative import sn_join_siblings
    from line_solver.lang.base import NodeType

    def _nt(nt):
        return int(nt.value) if hasattr(nt, 'value') else int(nt)

    model = _tpl_model(3, 2, 2)                 # 2 links, 2 tasks each -> 4
    sn = model.getStruct()
    j = [i for i, nt in enumerate(sn.nodetype) if _nt(nt) == _nt(NodeType.JOIN)][0]
    assert sn_join_siblings(sn, j, 0) == 4

    # Now raise ONE link to 3 tasks: 2 + 3 = 5, not 2*2 = 4 and not 2*3 = 6.
    model2 = _tpl_model(3, 2, 2)
    nodes = model2.getNodes()
    fork = [n for n in nodes if isinstance(n, Fork)][0]
    q2 = [n for n in nodes if n.getName() == 'Queue2'][0]
    fork.setTasksPerLink(3, model2.getClasses()[0], q2)
    sn2 = model2.getStruct()
    j2 = [i for i, nt in enumerate(sn2.nodetype) if _nt(nt) == _nt(NodeType.JOIN)][0]
    assert sn_join_siblings(sn2, j2, 0) == 5
