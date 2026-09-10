"""QRF BAS-blocking NLP bounds: qrf.bas.mmi, qrf.bas.mem and qrf.bas.bethe.

Port of MATLAB qrf_bas_mmi.m / qrf_bas_mem.m / qrf_bas_bethe.m. The polytope is
NOT re-derived here: it is the one `qr_bounds_bas` already builds and that the
LP token `qrf.bas` is validated on against the AMPL model, extracted through
`MapqnLpModel._assemble()`. Only the objective differs -- the LP maximizes one
station's utilization, these minimize mutual information (MMI), negative
entropy (MEM) or the tree-reweighted free entropy (BETHE) over the same
feasible set.

ONE POLYTOPE FOR EVERY BAS OBJECTIVE, HERE AND IN MATLAB SINCE 2026-09-03. The
ports take the LP assembler's family set. MATLAB implements the same five
THM30-group families once each, phase-aggregated, in every BAS routine -- but
`qrf_bas_mmi.m` also carried a second, stale AMPL-literal transcription of each,
emitting every one of them twice under different index conventions, so its
phase 1 could not reach feasibility (maximum equality residual 2.469e-01 on
cqn_bas_blocking) and `qrf.bas.mmi` was unusable. `qrf_bas_mem.m` had commented
those duplicates out, which is why it worked. Deleting them put all three MATLAB
objectives on the same polytope as these. See _kb/06-solver-catalog.md.

Re-transcribing the fifteen families instead would duplicate roughly 500 lines
whose index conventions are exactly where the MATLAB twin twice went wrong: the
THM30/THM3 population-to-index shift, and a MARGINALS sum whose upper limit was
taken on the 1-based index rather than the population (both fixed 2026-07-29).

The decision vector is the LP's own, so the objective is written against the
`p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}_{m}` variable names rather than the flat
`sub_qrfvar` layout of the no-blocking family.
"""

import numpy as np

from .lpmodel import MapqnLpModel
from .parameters import QRBoundsBasParameters
from .qrf_noblo_common import LOGTOL, solve_qrf_nlp
from .qr_bounds_bas import (
    _add_balance_constraints_bas,
    _add_bound_constraints_bas,
    _add_definition_constraints_bas,
    _add_thm3_constraints_bas,
    _add_zero_constraints_bas,
    _register_variables_bas,
)


def _bas_polytope(params):
    """The BAS feasible set as (Aeq, beq, Aub, bub, model), dense.

    Built by the same calls `mapqn_qr_bounds_bas` makes, so the two tokens
    cannot drift apart: a family added to the LP reaches the NLP with it.
    """
    model = MapqnLpModel()
    _register_variables_bas(model, params)
    _add_zero_constraints_bas(model, params)
    _add_definition_constraints_bas(model, params)
    _add_balance_constraints_bas(model, params)
    _add_thm3_constraints_bas(model, params)
    _add_bound_constraints_bas(model, params)

    A_ub, b_ub, A_eq, b_eq = model._assemble()
    n = model.get_num_variables()
    Aeq = A_eq.toarray() if A_eq is not None else np.zeros((0, n))
    beq = np.asarray(b_eq, dtype=float) if A_eq is not None else np.zeros(0)
    Aub = A_ub.toarray() if A_ub is not None else np.zeros((0, n))
    bub = np.asarray(b_ub, dtype=float) if A_ub is not None else np.zeros(0)
    return Aeq, beq, Aub, bub, model


def _p2_columns(model, params):
    """Column of every p2 entry, as a dict keyed by (j,nj,kj,i,ni,hi,m).

    Entries the model never registered (a population above a station's
    capacity) are absent, and the objective skips them.
    """
    M, F, K, MR = params.M, params.F, params.K, params.MR
    cols = {}
    for j in range(M):
        for nj in range(F[j] + 1):
            for kj in range(K[j]):
                for i in range(M):
                    for ni in range(F[i] + 1):
                        for hi in range(K[i]):
                            for m in range(MR):
                                idx = model.get_variable_index(
                                    'p2_%d_%d_%d_%d_%d_%d_%d'
                                    % (j, nj, kj, i, ni, hi, m))
                                if idx is not None:
                                    cols[(j, nj, kj, i, ni, hi, m)] = idx
    return cols


def _mmi_terms(cols, params, n_from=0):
    """(ij, ii, jj) column triples of the MI objective, i != j, ni,nj >= n_from.

    n_from is 0 for both callers, which is the AMPL source's own range
    (`sum {nj in 0..F[j]} sum {ni in 0..F[i]}`) and MATLAB's. This port used 1
    until 2026-09-03, dropping the idle cells: since U_i = 1 - p_i(0) that
    excluded the strongest correlation in a closed chain from the very
    functional meant to measure coupling, and it made qrf.bas.mmi disagree with
    MATLAB. The identical repair was made to the NO-BLOCKING arm on 2026-09-02;
    see qrf_noblo_common.mmi_objective. The parameter is kept because the range
    is the one thing that differs between these objectives and MEM's.
    """
    M, F, K, MR = params.M, params.F, params.K, params.MR
    ij, ii, jj = [], [], []
    for m in range(MR):
        for i in range(M):
            for ki in range(K[i]):
                for j in range(M):
                    if i == j:
                        continue
                    for kj in range(K[j]):
                        for ni in range(n_from, F[i] + 1):
                            for nj in range(n_from, F[j] + 1):
                                a = cols.get((i, ni, ki, j, nj, kj, m))
                                b = cols.get((i, ni, ki, i, ni, ki, m))
                                c = cols.get((j, nj, kj, j, nj, kj, m))
                                if a is None or b is None or c is None:
                                    continue
                                ij.append(a); ii.append(b); jj.append(c)
    return np.array(ij, dtype=int), np.array(ii, dtype=int), np.array(jj, dtype=int)


def _mem_terms(cols, params, n_from=1):
    """Diagonal columns of the MEM objective, ni >= n_from.

    MEM itself starts at 1, the range the AMPL source states for it. BETHE
    passes n_from=0: the idle cell is the strongest correlation in a closed
    chain, its MI block already spans n >= 0, and an entropy taken over a
    different range than the mutual information it is combined with is not a
    free entropy of anything.
    """
    M, F, K, MR = params.M, params.F, params.K, params.MR
    out = []
    for m in range(MR):
        for i in range(M):
            for k in range(K[i]):
                for ni in range(n_from, F[i] + 1):
                    idx = cols.get((i, ni, k, i, ni, k, m))
                    if idx is not None:
                        out.append(idx)
    return np.array(out, dtype=int)


def _mmi_pair(ij, ii, jj):
    """MMI objective and gradient over the LP variable vector."""
    def objective(x):
        pij, pii, pjj = x[ij], x[ii], x[jj]
        return float(np.sum(pij * (np.log(LOGTOL + pij) - np.log(LOGTOL + pii)
                                   - np.log(LOGTOL + pjj))))

    def gradient(x):
        g = np.zeros_like(x)
        pij, pii, pjj = x[ij], x[ii], x[jj]
        np.add.at(g, ij, np.log(LOGTOL + pij) - np.log(LOGTOL + pii)
                  - np.log(LOGTOL + pjj) + pij / (LOGTOL + pij))
        np.add.at(g, ii, -pij / (LOGTOL + pii))
        np.add.at(g, jj, -pij / (LOGTOL + pjj))
        return g

    return objective, gradient


def _mem_pair(diag):
    """MEM objective and gradient over the LP variable vector.

    Returned as the NEGATIVE entropy, because _solve minimises: the AMPL model
    states this objective as `maximize H`. Returning +H (as every port did
    until 2026-08-29) selects the minimum-entropy face instead.
    """
    def objective(x):
        p = x[diag]
        return float(np.sum(p * np.log(LOGTOL + p)))

    def gradient(x):
        g = np.zeros_like(x)
        p = x[diag]
        np.add.at(g, diag, np.log(LOGTOL + p) + p / (LOGTOL + p))
        return g

    return objective, gradient


def _bethe_pair(ij, ii, jj, diag, M):
    """Tree-reweighted (Bethe) free entropy and its gradient, at lambda = 1/M.

        lam * sum_{i!=j} I(n_i;n_j) - sum_i H(n_i),      lam = 1/M

    which is lam times the MI body plus the MEM body, both over n >= 0. lam is
    the uniform point 2/M of the spanning-tree polytope of K_M halved, the
    largest uniform edge weight at which the tree-reweighted entropy is concave
    and so the program convex. See qrf_noblo_common.bethe_objective, whose
    functional this is over the BAS decision vector.
    """
    lam = 1.0 / M

    def objective(x):
        pij, pii, pjj = x[ij], x[ii], x[jj]
        mi = float(np.sum(pij * (np.log(LOGTOL + pij) - np.log(LOGTOL + pii)
                                 - np.log(LOGTOL + pjj))))
        p = x[diag]
        return lam * mi + float(np.sum(p * np.log(LOGTOL + p)))

    def gradient(x):
        g = np.zeros_like(x)
        pij, pii, pjj = x[ij], x[ii], x[jj]
        np.add.at(g, ij, lam * (np.log(LOGTOL + pij) - np.log(LOGTOL + pii)
                                - np.log(LOGTOL + pjj) + pij / (LOGTOL + pij)))
        np.add.at(g, ii, -lam * pij / (LOGTOL + pii))
        np.add.at(g, jj, -lam * pij / (LOGTOL + pjj))
        p = x[diag]
        np.add.at(g, diag, np.log(LOGTOL + p) + p / (LOGTOL + p))
        return g

    return objective, gradient


def _metrics(x, cols, model, params):
    """Utilization and queue length from the optimal point.

    UTILIZATION, not occupancy. UN used to sum the diagonal p2 over ALL blocking
    configurations, i.e. P(n_i >= 1) with the BLOCKED ones included. A blocked
    BAS server holds a job it has already finished and does no work, so that is
    occupancy: on the M=2, N=3, F=[2 3] cyclic model it reported U2 = 1 where the
    exact utilization is 7/15, which the LP over the SAME polytope already
    returns. The e variables carry the right quantity -- UEFF pins e(i,ki) to the
    mass with n_i >= 1 in the configurations where i is NOT blocked.

    NO 1/M, matching the LP builder this shares: `qr_bounds_bas` emits UEFF as
    one row per (j, i, ki), which leaves e unscaled. QN stays on the diagonal p2
    over every configuration, blocked included, because a blocked job is still
    held at the station and counts towards its population.
    """
    M, F, K, MR = params.M, params.F, params.K, params.MR
    UN = np.zeros(M)
    QN = np.zeros(M)
    for i in range(M):
        for ki in range(K[i]):
            idx = model.get_variable_index('e_%d_%d' % (i, ki))
            if idx is not None:
                UN[i] += x[idx]
        for m in range(MR):
            for ni in range(1, F[i] + 1):
                for ki in range(K[i]):
                    idx = cols.get((i, ni, ki, i, ni, ki, m))
                    if idx is None:
                        continue
                    QN[i] += ni * x[idx]
    return UN, QN


def qrf_bas_mmi(params: QRBoundsBasParameters):
    """Minimum-mutual-information bound on the BAS-blocking polytope.

    Args:
        params: the same QRBoundsBasParameters the LP token qrf.bas takes.

    Returns:
        (UN, QN): utilization and queue length per station.
    """
    Aeq, beq, Aub, bub, model = _bas_polytope(params)
    cols = _p2_columns(model, params)
    ij, ii, jj = _mmi_terms(cols, params)
    if ij.size == 0:
        raise ValueError("qrf_bas_mmi: the model has no off-diagonal joint "
                         "variables, so the MMI objective is empty.")
    objective, gradient = _mmi_pair(ij, ii, jj)
    x = _solve(objective, gradient, Aeq, beq, Aub, bub, 'qrf_bas_mmi')
    return _metrics(x, cols, model, params)


def qrf_bas_mem(params: QRBoundsBasParameters):
    """Maximum-entropy bound on the BAS-blocking polytope."""
    Aeq, beq, Aub, bub, model = _bas_polytope(params)
    cols = _p2_columns(model, params)
    diag = _mem_terms(cols, params)
    if diag.size == 0:
        raise ValueError("qrf_bas_mem: the model has no marginal variables, "
                         "so the MEM objective is empty.")
    objective, gradient = _mem_pair(diag)
    x = _solve(objective, gradient, Aeq, beq, Aub, bub, 'qrf_bas_mem')
    return _metrics(x, cols, model, params)


def qrf_bas_bethe(params: QRBoundsBasParameters):
    """Tree-reweighted (Bethe) free entropy bound on the BAS-blocking polytope.

    The BAS twin of qrf_noblo_bethe: same polytope as qrf.bas and qrf.bas.mmi,
    same phase-1 start, and the objective of qrf_noblo_bethe evaluated over the
    BAS decision vector, so the blocking configurations and the per-station
    capacities enter through the ranges alone.

    Args:
        params: the same QRBoundsBasParameters the LP token qrf.bas takes.

    Returns:
        (UN, QN): utilization and queue length per station.
    """
    Aeq, beq, Aub, bub, model = _bas_polytope(params)
    cols = _p2_columns(model, params)
    ij, ii, jj = _mmi_terms(cols, params, n_from=0)
    diag = _mem_terms(cols, params, n_from=0)
    if ij.size == 0:
        raise ValueError("qrf_bas_bethe: the model has no off-diagonal joint "
                         "variables, so the mutual-information term of the "
                         "Bethe objective is empty.")
    if diag.size == 0:
        raise ValueError("qrf_bas_bethe: the model has no marginal variables, "
                         "so the entropy term of the Bethe objective is empty.")
    objective, gradient = _bethe_pair(ij, ii, jj, diag, params.M)
    x = _solve(objective, gradient, Aeq, beq, Aub, bub, 'qrf_bas_bethe')
    return _metrics(x, cols, model, params)


def _solve(objective, gradient, Aeq, beq, Aub, bub, name):
    """Phase-1 LP onto the polytope, then the NLP in the null space of Aeq."""
    from scipy.optimize import linprog
    n = Aeq.shape[1] if Aeq.size else Aub.shape[1]
    res = linprog(np.zeros(n),
                  A_ub=Aub if Aub.shape[0] else None,
                  b_ub=bub if Aub.shape[0] else None,
                  A_eq=Aeq if Aeq.shape[0] else None,
                  b_eq=beq if Aeq.shape[0] else None,
                  bounds=[(0.0, 1.0)] * n, method='highs')
    if not res.success:
        raise ValueError("%s: the BAS polytope is infeasible (phase-1 LP "
                         "status %d: %s)" % (name, res.status, res.message))
    return solve_qrf_nlp(objective, gradient, res.x, Aeq, beq, Aub, bub, name)
