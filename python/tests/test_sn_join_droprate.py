"""Tests for sn_join_siblings / sn_join_droprate and the Join row of getAvgLossTable.

A Join is the one station where LossRate = ArvR - Tput does not hold: ArvR counts
the SIBLINGS offered and Tput the PARENT jobs released. Reading the identity there
charges (N-1)/N of the offered traffic as lost at every join, standard joins
included; these tests pin the rule that replaces it.
"""

import numpy as np
import pytest

from line_solver import SolverMVA, SolverNC
from line_solver.api.fjnative import sn_join_droprate, sn_join_siblings
from line_solver.gallery import gallery_fj_closed, gallery_fj_quorum
from line_solver.lang.base import NodeType


def _nt(nt):
    return int(nt.value) if hasattr(nt, 'value') else int(nt)


def _join_node(sn):
    return [i for i, nt in enumerate(sn.nodetype) if _nt(nt) == _nt(NodeType.JOIN)][0]


def test_siblings_counted_at_the_fork():
    sn = gallery_fj_quorum().getStruct()
    assert sn_join_siblings(sn, _join_node(sn)) == 3
    sn2 = gallery_fj_closed().getStruct()
    assert sn_join_siblings(sn2, _join_node(sn2)) == 2


def test_standard_join_loses_nothing():
    """AN = N*TN at a standard join, and all N siblings are consumed."""
    model = gallery_fj_closed()
    sn = model.getStruct()
    ist = int(np.asarray(sn.nodeToStation).flatten()[_join_node(sn)])
    M, K = int(sn.nstations), int(sn.nclasses)
    TN = np.zeros((M, K))
    AN = np.zeros((M, K))
    TN[ist, 0] = 0.75
    AN[ist, 0] = 2 * 0.75
    d = sn_join_droprate(sn, TN, AN)
    assert d[ist, 0] == pytest.approx(0.0, abs=1e-12)


def test_quorum_discards_the_stragglers():
    """A 2-of-3 join consumes 2 siblings per firing and throws the third away."""
    model = gallery_fj_quorum()
    sn = model.getStruct()
    ist = int(np.asarray(sn.nodeToStation).flatten()[_join_node(sn)])
    M, K = int(sn.nstations), int(sn.nclasses)
    TN = np.zeros((M, K))
    AN = np.zeros((M, K))
    x = 1.6
    TN[ist, 0] = x
    AN[ist, 0] = 3 * x
    d = sn_join_droprate(sn, TN, AN)
    assert d[ist, 0] == pytest.approx((3 - 2) * x, rel=1e-12)


def test_droprate_is_zero_away_from_the_join():
    model = gallery_fj_quorum()
    sn = model.getStruct()
    ist = int(np.asarray(sn.nodeToStation).flatten()[_join_node(sn)])
    M, K = int(sn.nstations), int(sn.nclasses)
    TN = np.full((M, K), 1.6)
    AN = np.full((M, K), 4.8)
    d = sn_join_droprate(sn, TN, AN)
    for i in range(M):
        if i != ist:
            assert d[i, 0] == pytest.approx(0.0, abs=1e-12)


def test_droprate_never_negative():
    """A quorum met by more siblings than are forked is a full join, not a gain."""
    model = gallery_fj_quorum()
    sn = model.getStruct()
    ist = int(np.asarray(sn.nodeToStation).flatten()[_join_node(sn)])
    M, K = int(sn.nstations), int(sn.nclasses)
    TN = np.zeros((M, K))
    AN = np.zeros((M, K))
    TN[ist, 0] = 5.0   # inconsistent with AN on purpose
    AN[ist, 0] = 1.0
    d = sn_join_droprate(sn, TN, AN)
    assert d[ist, 0] >= 0.0


@pytest.mark.parametrize('solver_cls', [SolverMVA, SolverNC])
def test_loss_table_join_row_standard(solver_cls):
    """The old ArvR - Tput rule reported LossRatio = 1/2 here. Nothing is lost."""
    table = solver_cls(gallery_fj_closed()).getAvgLossTable()
    df = table.data if hasattr(table, 'data') else table
    row = df[df['Station'] == 'Join']
    assert len(row) == 1
    # The residual is the solver's own AN-vs-TN convergence gap, not a loss; the
    # rule under test is that it is not 1/2, which is what ArvR - Tput reported.
    assert row['LossRate'].iloc[0] == pytest.approx(0.0, abs=1e-6)
    assert row['LossRatio'].iloc[0] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize('solver_cls', [SolverMVA, SolverNC])
def test_loss_table_join_row_quorum(solver_cls):
    """One of the three siblings is discarded, so the loss ratio is exactly 1/3."""
    table = solver_cls(gallery_fj_quorum()).getAvgLossTable()
    df = table.data if hasattr(table, 'data') else table
    row = df[df['Station'] == 'Join']
    assert len(row) == 1
    assert row['LossRatio'].iloc[0] == pytest.approx(1.0 / 3.0, rel=1e-6)
    # and the rate is one sibling per synchronisation
    assert row['LossRate'].iloc[0] == pytest.approx(row['ArvR'].iloc[0] / 3.0, rel=1e-6)


def test_quorum_chain_capacity_is_unbounded():
    """A quorum join leaves stragglers in flight past their parent's next fork, so
    a branch station is NOT bounded by the class population."""
    snq = gallery_fj_quorum().getStruct()
    assert np.all(np.isinf(np.asarray(snq.classcap)))
    sns = gallery_fj_closed().getStruct()
    assert np.all(np.isfinite(np.asarray(sns.classcap)))
