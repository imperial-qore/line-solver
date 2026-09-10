"""
Admission constraints on an LN layer station.

A Task or Processor may declare ``A*n <= b`` on the station that represents it
in its layer (``addConstraint`` / ``setConstraint``), which is how a passive
resource (semaphore, buffer pool, admission method name) is expressed in an LQN.
``SolverLN._add_layer_admission_constraint`` emits the rows as a Region on the
layer's server station, expanding an entry column to the CALL classes targeting
that entry, and ``_region_capable_layer_solver`` escalates the layer solver to
one whose feature set covers Region.

There is no external oracle: JMT's FCR XML cannot express a general A*n <= b
and lqns has no admission-constraint concept, so these tests assert structural
and conservation properties rather than golden numbers. Twin of
``line-test.git/test/testsAdvFeatures/fcr/test_fcr_lincon_ln.m``.
"""
import numpy as np
import pytest

from line_solver import (LayeredNetwork, Processor, Task, Entry, Activity,
                         ActivityPrecedence, SchedStrategy, Exp, SolverLN)
from line_solver.solvers.solver_mva.solver_mva import SolverMVA

TOL_EQ = 1e-6
TOL_BAL = 1e-3


def mk_lqn():
    """Two-tier LQN: T1 (reference, N=2) calls entries E2 and E3 of T2 in sequence."""
    model = LayeredNetwork('fcrlqn')
    p1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    p2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    t1 = Task(model, 'T1', 2, SchedStrategy.REF).on(p1)
    t2 = Task(model, 'T2', 2, SchedStrategy.FCFS).on(p2)
    e1 = Entry(model, 'E1').on(t1)
    e2 = Entry(model, 'E2').on(t2)
    e3 = Entry(model, 'E3').on(t2)
    a1 = Activity(model, 'A1', Exp(1.0)).on(t1).bound_to(e1).synch_call(e2, 1.0)
    a1b = Activity(model, 'A1b', Exp(1.0)).on(t1).synch_call(e3, 1.0)
    t1.add_precedence(ActivityPrecedence.Serial(a1, a1b))
    Activity(model, 'A2', Exp(2.0)).on(t2).bound_to(e2).replies_to(e2)
    Activity(model, 'A3', Exp(2.0)).on(t2).bound_to(e3).replies_to(e3)
    return model, t2, e2, e3


def with_cap(cap):
    """Model whose T2 admits at most CAP calls across its two entries."""
    model, t2, e2, e3 = mk_lqn()
    t2.addConstraint([e2, e3], [1, 1], cap)
    return model


def avg_table(model):
    return SolverLN(model).getAvgTable()


def row_of(table, node_name, col):
    """Value of column COL on the row whose Node is NODE_NAME.

    Node is a column when the table has been printed and the index otherwise,
    so both are accepted.
    """
    names = list(table['Node']) if 'Node' in table.columns else list(table.index)
    assert node_name in names, f"row {node_name} not in {names}"
    # positional, since boolean masking misaligns on this table's index
    return float(table.iloc[names.index(node_name)][col])


def lincon_of(model, tidx_offset=1):
    lqn = model.getStruct()
    return lqn.lincon[lqn.tshift + tidx_offset]


def test_declaration_forms_agree():
    """Positional, named-by-handle and named-by-string produce identical rows."""
    m1, t2a, _, _ = mk_lqn()
    t2a.setConstraint([[1, 1], [0, 1]], [2, 1])
    a1, b1 = lincon_of(m1)

    m2, t2b, e2b, e3b = mk_lqn()
    t2b.addConstraint([e2b, e3b], [1, 1], 2)
    t2b.addConstraint(e3b, 1, 1)
    a2, b2 = lincon_of(m2)

    m3, t2c, _, _ = mk_lqn()
    t2c.addConstraint(['E2', 'E3'], None, 2)  # None coeffs default to all ones
    t2c.addConstraint('E3', 1, 1)
    a3, b3 = lincon_of(m3)

    assert a1.shape == (2, 2)
    assert np.array_equal(a1, a2) and np.array_equal(b1, b2)
    assert np.array_equal(a1, a3) and np.array_equal(b1, b3)


def test_foreign_operand_rejected():
    """An operand that is not an entry of this task is an error, not a silent mis-map."""
    model, t2, _, _ = mk_lqn()
    t2.addConstraint('E1', 1, 1)  # E1 belongs to T1
    with pytest.raises(ValueError, match='not one of the entries'):
        model.getStruct()


def test_duplicate_operand_rejected():
    _, t2, _, _ = mk_lqn()
    with pytest.raises(ValueError, match='same operand more than once'):
        t2.addConstraint(['E2', 'E2'], [1, 1], 2)


def test_wrong_column_count_rejected():
    """The positional form checks the column count, which is its only guard."""
    model, t2, _, _ = mk_lqn()
    t2.setConstraint([[1, 1, 1]], [2])
    with pytest.raises(ValueError, match='columns but there are'):
        model.getStruct()


def test_layer_carries_region_and_escalates_off_mva():
    solver = SolverLN(with_cap(1))
    solver.getAvgTable()
    constrained = [e for e, lay in enumerate(solver.ensemble)
                   if getattr(lay, 'regions', None)]
    assert len(constrained) == 1
    assert not isinstance(solver.solvers[constrained[0]], SolverMVA)


def test_constraint_is_enforced_and_satisfied():
    free = avg_table(mk_lqn()[0])
    con = avg_table(with_cap(1))
    q_free = row_of(free, 'E2', 'QLen')
    q_e2 = row_of(con, 'E2', 'QLen')
    q_e3 = row_of(con, 'E3', 'QLen')
    assert abs(q_e2 - q_free) > TOL_EQ, "constraint is inert"
    assert q_e2 + q_e3 <= 1 + TOL_BAL


def test_flow_balance_across_layer_boundary():
    """Caller and callee must agree on throughput.

    A job blocked at the region is counted at no station, so without the
    chain-level Little recovery in the call-metric update LN converges to a
    fixed point where the calling activity and the called entry disagree.
    """
    con = avg_table(with_cap(1))
    x_a1 = row_of(con, 'A1', 'Tput')
    x_e2 = row_of(con, 'E2', 'Tput')
    assert abs(x_e2 - x_a1) / abs(x_a1) <= TOL_BAL


def test_throughput_monotone_in_cap():
    """Relaxing the cap cannot reduce throughput.

    Compared only within the constrained family: a constrained layer is solved
    by SolverCTMC and an unconstrained one by SolverMVA, so the
    constrained-vs-free gap carries a solver discrepancy on top of the
    constraint effect.
    """
    x1 = row_of(avg_table(with_cap(1)), 'E2', 'Tput')
    x2 = row_of(avg_table(with_cap(2)), 'E2', 'Tput')
    x3 = row_of(avg_table(with_cap(3)), 'E2', 'Tput')
    assert x1 <= x2 + TOL_EQ
    assert x2 <= x3 + TOL_EQ


def test_slack_cap_matches_unconstrained():
    """A cap that cannot bind (3 method names for 2 jobs) reproduces the free solution
    up to the CTMC-vs-MVA layer-solver discrepancy."""
    x_slack = row_of(avg_table(with_cap(3)), 'E2', 'Tput')
    x_free = row_of(avg_table(mk_lqn()[0]), 'E2', 'Tput')
    assert abs(x_slack - x_free) / abs(x_free) <= 0.01


def test_json_round_trip_preserves_constraint():
    """The wire form is named, so it is order-independent by construction."""
    import json
    import os
    import tempfile
    from line_solver.io.linemodel_io import save_model, load_model

    model, t2, _, _ = mk_lqn()
    t2.setConstraint([[1, 1], [0, 1]], [2, 1])
    path = os.path.join(tempfile.gettempdir(), 'ln_lincon_rt.json')
    save_model(model, path)
    with open(path) as fh:
        wire = json.load(fh)
    rows = [t.get('admissionConstraints') for t in wire['model']['tasks']
            if t['name'] == 'T2'][0]
    assert [r['operands'] for r in rows] == [['E2', 'E3'], ['E3']]

    a0, b0 = lincon_of(model)
    a1, b1 = lincon_of(load_model(path))
    assert np.array_equal(a0, a1) and np.array_equal(b0, b1)
