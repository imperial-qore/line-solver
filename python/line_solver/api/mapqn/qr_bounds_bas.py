"""
QR Bounds with Blocking-After-Service (BAS) Protocol.

Implements quadratic reduction bounds using the BAS blocking protocol
for MAP queueing networks with finite capacity queues.

This is a port of qrf_bas.m MATLAB model.
"""

import numpy as np
from typing import Optional, List

from .lpmodel import MapqnLpModel
from .parameters import QRBoundsBasParameters
from .solution import MapqnSolution


def mapqn_qr_bounds_bas(
    params: QRBoundsBasParameters,
    objective_queue: int,
    sense: str = "min"
) -> MapqnSolution:
    """
    Solve the QR bounds for BAS blocking networks.

    Computes bounds on utilization for closed queueing networks with
    finite capacity queues using Blocking-After-Service protocol.

    Args:
        params: BAS parameters including:
            - M: Number of queues
            - N: Population
            - MR: Number of blocking configurations
            - f: Finite capacity queue index (1-based)
            - K: Number of phases per queue [M]
            - F: Capacity per queue [M]
            - MM, MM1, ZZ, BB: Blocking configuration matrices
            - mu: Service rate matrices
            - v: Background transition matrices
            - r: Routing matrix
        objective_queue: Queue index to optimize (1-based).
        sense: Optimization sense - "min" or "max".

    Returns:
        MapqnSolution containing:
            - objective_value: Optimal utilization bound
            - variables: Dictionary with U, Ueff, pb per queue

    Raises:
        ValueError: If indices are out of range or sense is invalid.
        RuntimeError: If LP is infeasible or unbounded.

    Reference:
        Based on qrf_bas.m from qrf-revised repository
    """
    params.validate()

    M = params.M
    N = params.N
    F = params.F
    K = params.K
    MR = params.MR
    f = params.f - 1  # Convert to 0-based

    if not (1 <= objective_queue <= M):
        raise ValueError(f"Objective queue must be in range 1..{M}")
    if sense not in ("min", "max"):
        raise ValueError("Sense must be 'min' or 'max'")

    model = MapqnLpModel()

    # Register all variables
    _register_variables_bas(model, params)

    # Add constraints
    _add_zero_constraints_bas(model, params)
    _add_definition_constraints_bas(model, params)
    _add_balance_constraints_bas(model, params)
    _add_thm3_constraints_bas(model, params)
    _add_bound_constraints_bas(model, params)

    # Build objective: sum of diagonal p2 at target queue across all blocking configs
    target = objective_queue - 1  # 0-based
    objective_terms = {}
    for m in range(MR):
        for ki in range(K[target]):
            for ni in range(1, F[target] + 1):
                var_name = f'p2_{target}_{ni}_{ki}_{target}_{ni}_{ki}_{m}'
                objective_terms[var_name] = 1.0

    # see _kb/03-api-layer.md for rationale
    if objective_terms:
        c = np.zeros(model.get_num_variables())
        for var_name, coef in objective_terms.items():
            vidx = model.get_variable_index(var_name)
            if vidx is not None:
                c[vidx] += coef
        solution = model.solve_with_objective(c, minimize=(sense == "min"))

        # Compute actual objective from solution
        obj_value = 0.0
        for var_name, coef in objective_terms.items():
            obj_value += coef * solution.variables.get(var_name, 0.0)

        # Update solution with derived variables
        variables = dict(solution.variables)

        # Compute U, Ueff, pb for each queue
        for i in range(M):
            total_u = 0.0
            total_e = 0.0
            for m in range(MR):
                for ki in range(K[i]):
                    for ni in range(1, F[i] + 1):
                        p2_val = variables.get(f'p2_{i}_{ni}_{ki}_{i}_{ni}_{ki}_{m}', 0.0)
                        total_u += p2_val
                    e_val = variables.get(f'e_{i}_{ki}', 0.0)
                    total_e += e_val
            variables[f'U_{i + 1}'] = total_u
            variables[f'Ueff_{i + 1}'] = total_e
            variables[f'pb_{i + 1}'] = total_u - total_e

        return MapqnSolution(
            objective_value=obj_value,
            variables=variables
        )
    else:
        return MapqnSolution(objective_value=0.0, variables={})


def _register_variables_bas(model: MapqnLpModel, params: QRBoundsBasParameters) -> None:
    """Register all variables for BAS model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K
    MR = params.MR

    # p2 variables: p2[j,nj,kj,i,ni,hi,m]
    for j in range(M):
        for nj in range(N + 1):
            for kj in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for hi in range(K[i]):
                            for m in range(MR):
                                model.add_variable(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}_{m}', lb=0.0, ub=1.0)

    # e variables: e[i,ki] - effective utilization
    for i in range(M):
        for ki in range(K[i]):
            model.add_variable(f'e_{i}_{ki}', lb=0.0, ub=1.0)


def _add_zero_constraints_bas(model: MapqnLpModel, params: QRBoundsBasParameters) -> None:
    """Add ZERO constraints for infeasible states in BAS model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K
    MR = params.MR
    BB = params.BB
    f = params.f - 1  # 0-based

    for j in range(M):
        for nj in range(N + 1):
            for kj in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for hi in range(K[i]):
                            for m in range(MR):
                                is_zero = False

                                # ZERO1: i==j, nj==ni, h!=k
                                if i == j and nj == ni and hi != kj:
                                    is_zero = True

                                # ZERO2: i==j, nj!=ni
                                if i == j and nj != ni:
                                    is_zero = True

                                # ZERO3: i!=j, nj+ni > N
                                if i != j and nj + ni > N:
                                    is_zero = True

                                # see _kb/03-api-layer.md for rationale
                                if nj > F[j] or ni > F[i]:
                                    is_zero = True

                                # ZERO5: BB[m,j]==1 and nj==0
                                if m >= 1 and BB[m, j] == 1 and nj == 0:
                                    is_zero = True

                                # ZERO7: BB[m,j]==1 and i!=j and i!=f and ni+nj+F[f]>N
                                if m >= 1 and BB[m, j] == 1 and i != j and i != f and ni + nj + F[f] > N:
                                    is_zero = True

                                # ZERO8: finite queue not at capacity in blocking config
                                if j == f and 1 <= nj <= F[f] - 1 and m >= 1:
                                    is_zero = True

                                if is_zero:
                                    constraint = model.constraint_builder()
                                    constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}_{m}', 1.0)
                                    model.add_constraint(constraint.eq(0.0))

    # ZERO4: For m>=1 and j!=f, p2[j,nj,k,f,nf,h,m]=0 when nf < F[f]
    for j in range(M):
        if j == f:
            continue
        for nj in range(N + 1):
            for kj in range(K[j]):
                for m in range(1, MR):
                    for nf in range(F[f]):
                        for hf in range(K[f]):
                            constraint = model.constraint_builder()
                            constraint.add_term(f'p2_{j}_{nj}_{kj}_{f}_{nf}_{hf}_{m}', 1.0)
                            model.add_constraint(constraint.eq(0.0))


def _add_definition_constraints_bas(model: MapqnLpModel, params: QRBoundsBasParameters) -> None:
    """Add definition constraints for BAS model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K
    MR = params.MR
    BB = params.BB

    # ONE: Normalization
    for j in range(M):
        constraint = model.constraint_builder()
        for nj in range(N + 1):
            for kj in range(K[j]):
                for m in range(MR):
                    constraint.add_term(f'p2_{j}_{nj}_{kj}_{j}_{nj}_{kj}_{m}', 1.0)
        model.add_constraint(constraint.eq(1.0))

    # SYMMETRY
    for j in range(M):
        for nj in range(min(N, F[j]) + 1):
            for kj in range(K[j]):
                for i in range(j + 1, M):
                    for ni in range(min(N, F[i]) + 1):
                        if i != j and nj + ni > N:
                            continue
                        for hi in range(K[i]):
                            for m in range(MR):
                                constraint = model.constraint_builder()
                                constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}_{m}', 1.0)
                                constraint.add_term(f'p2_{i}_{ni}_{hi}_{j}_{nj}_{kj}_{m}', -1.0)
                                model.add_constraint(constraint.eq(0.0))

    # MARGINALS
    for j in range(M):
        for kj in range(K[j]):
            for nj in range(min(N, F[j]) + 1):
                for i in range(M):
                    if i == j:
                        continue
                    for m in range(MR):
                        constraint = model.constraint_builder()
                        constraint.add_term(f'p2_{j}_{nj}_{kj}_{j}_{nj}_{kj}_{m}', 1.0)
                        for ni in range(min(N - nj, F[i]) + 1):
                            for hi in range(K[i]):
                                constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}_{m}', -1.0)
                        model.add_constraint(constraint.eq(0.0))

    # see _kb/03-api-layer.md for rationale
    for j in range(M):
        for i in range(M):
            for ki in range(K[i]):
                constraint = model.constraint_builder()
                constraint.add_term(f'e_{i}_{ki}', -1.0)
                for nj in range(min(N, F[j]) + 1):
                    for kj in range(K[j]):
                        for m in range(MR):
                            if BB[m, i] == 0:
                                for ni in range(1, min(N, F[i]) + 1):
                                    constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{ki}_{m}', 1.0)
                model.add_constraint(constraint.eq(0.0))


def _add_balance_constraints_bas(model: MapqnLpModel, params: QRBoundsBasParameters) -> None:
    """Add balance constraints for BAS model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K
    MR = params.MR
    BB = params.BB
    MM = params.MM
    ZZ = params.ZZ
    ZM = params.ZM
    f = params.f - 1  # 0-based

    def q(i: int, j: int, k: int, h: int) -> float:
        return params.q(i, j, k, h)

    # THM1: Phase balance
    for i in range(M):
        for ki in range(K[i]):
            constraint = model.constraint_builder()
            # LHS
            for j in range(M):
                for hi in range(K[i]):
                    if j != i or hi != ki:
                        constraint.add_term(f'e_{i}_{ki}', q(i, j, ki, hi))
            # RHS (subtract)
            for j in range(M):
                for hi in range(K[i]):
                    if j != i or hi != ki:
                        constraint.add_term(f'e_{i}_{hi}', -q(i, j, hi, ki))
            model.add_constraint(constraint.eq(0.0))

    # THM2: Population constraint
    for j in range(M):
        for kj in range(K[j]):
            for nj in range(F[j] + 1):
                for m in range(MR):
                    constraint = model.constraint_builder()
                    constraint.add_term(f'p2_{j}_{nj}_{kj}_{j}_{nj}_{kj}_{m}', -float(N))
                    for i in range(M):
                        for ni in range(1, F[i] + 1):
                            for ki in range(K[i]):
                                constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{ki}_{m}', float(ni))
                    model.add_constraint(constraint.eq(0.0))

    # COR1: Second moment constraint
    constraint = model.constraint_builder()
    for m in range(MR):
        for i in range(M):
            for j in range(M):
                for nj in range(1, F[j] + 1):
                    for ni in range(1, F[i] + 1):
                        for ki in range(K[i]):
                            for kj in range(K[j]):
                                constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{ki}_{m}', float(ni * nj))
    model.add_constraint(constraint.eq(float(N * N)))


def _add_thm3_constraints_bas(model: MapqnLpModel, params: QRBoundsBasParameters) -> None:
    """Add the THM3 family of marginal-balance constraints (BAS protocol).

    Ports THM30, THM3, THM3f, THM3I and THM3L from the authoritative AMPL model
    qrboundsbas_skel.mod (lines 53-57). These were absent from this port; a
    missing constraint enlarges the feasible region, so their omission LOOSENED
    the bound (the same defect class as the QBAL gap in qr_bounds_rsrd).

    AMPL is 1-based in i/j/k/h/u/m/w/p and 0-based in the populations nj/ni;
    everything below is 0-based, so AMPL `MM[m,1]` is `MM[m,0]`, and a config
    index stored as a VALUE (MM[m,1], MM1[m,j]) is 1-based and must have 1
    subtracted before it indexes p2's m slot.
    """
    M = params.M
    F = params.F
    K = params.K
    MR = params.MR
    BB = params.BB
    MM = params.MM
    MM1 = params.MM1
    ZZ = params.ZZ
    ZM = params.ZM
    f = params.f - 1  # 0-based

    def q(i: int, j: int, k: int, h: int) -> float:
        return params.q(i, j, k, h)

    def p2(j, nj, kj, i, ni, hi, m) -> str:
        return f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}_{m}'

    _thm30(model, params, q, p2)
    _thm3(model, params, q, p2)
    _thm3f(model, params, q, p2)
    _thm3i(model, params, q, p2)
    _thm3l(model, params, q, p2)


def _unpack(params):
    return (params.M, params.F, params.K, params.MR, params.BB,
            params.MM, params.MM1, params.ZZ, params.ZM, params.f - 1)


def _thm30(model, params, q, p2):
    """THM30: ni = 0 boundary case of THM3 (AMPL qrboundsbas_skel.mod:53)."""
    M, F, K, MR, BB, MM, MM1, ZZ, ZM, f = _unpack(params)
    for i in range(M):
        if i == f:
            continue
        for u in range(K[i]):
            c = model.constraint_builder()
            # LHS term 1: j != i, j != f, BB[m,j] == 0
            # LHS term 2: j != i, j == f, MM[m,0] != i+1
            for j in range(M):
                if j == i:
                    continue
                for m in range(MR):
                    if j != f:
                        if BB[m, j] != 0:
                            continue
                    else:
                        if int(MM[m, 0]) == i + 1:
                            continue
                    for nj in range(1, F[j] + 1):
                        for k in range(K[j]):
                            coef = 0.0
                            for h in range(K[j]):
                                coef += q(j, i, k, h)
                            if coef != 0.0:
                                c.add_term(p2(j, nj, k, i, 0, u, m), coef)
            # RHS terms 1 and 2: q[i,j,k,u] * p2[j,nj,h, i,1,k, m]
            for j in range(M):
                if j == i:
                    continue
                for m in range(MR):
                    if BB[m, i] != 0:
                        continue
                    for nj in range(F[j] + 1):
                        if j == f and nj >= F[j]:
                            continue
                        for k in range(K[i]):
                            coef = q(i, j, k, u)
                            if coef == 0.0:
                                continue
                            for h in range(K[j]):
                                c.add_term(p2(j, nj, h, i, 1, k, m), -coef)
            # RHS term 3: blocked f at capacity handing over to i
            for m in range(MR):
                if BB[m, i] != 1 or int(MM[m, 0]) != i + 1:
                    continue
                j = f
                if j == i:
                    continue
                nj = F[j]
                # see _kb/03-api-layer.md for rationale
                for y in range(K[j]):
                    coef = 0.0
                    for w in range(M):
                        if w == f or w == i:
                            continue
                        for p in range(K[f]):
                            coef += q(f, w, y, p)
                    if coef != 0.0:
                        c.add_term(p2(f, nj, y, i, 1, u, m), -coef)
            model.add_constraint(c.eq(0.0))


def _thm3(model, params, q, p2):
    """THM3: marginal balance for ni in 0..F[i]-1, i != f (skel:54)."""
    M, F, K, MR, BB, MM, MM1, ZZ, ZM, f = _unpack(params)
    for i in range(M):
        if i == f:
            continue
        for ni in range(F[i]):
            c = model.constraint_builder()
            # LHS terms 1 and 2
            for j in range(M):
                if j == i:
                    continue
                for m in range(MR):
                    if j != f:
                        if BB[m, j] != 0:
                            continue
                    else:
                        if int(MM[m, 0]) == i + 1:
                            continue
                    for nj in range(1, F[j] + 1):
                        for k in range(K[j]):
                            coef = 0.0
                            for h in range(K[j]):
                                coef += q(j, i, k, h)
                            if coef == 0.0:
                                continue
                            for u in range(K[i]):
                                c.add_term(p2(j, nj, k, i, ni, u, m), coef)
            # RHS terms 1 and 2
            for j in range(M):
                if j == i:
                    continue
                for m in range(MR):
                    if BB[m, i] != 0:
                        continue
                    for nj in range(F[j] + 1):
                        if j == f and nj >= F[j]:
                            continue
                        for k in range(K[i]):
                            coef = 0.0
                            for h in range(K[i]):
                                coef += q(i, j, k, h)
                            if coef == 0.0:
                                continue
                            for u in range(K[j]):
                                c.add_term(p2(j, nj, u, i, ni + 1, k, m), -coef)
            # RHS term 3
            for m in range(MR):
                if BB[m, i] != 1 or int(MM[m, 0]) != i + 1:
                    continue
                j = f
                if j == i:
                    continue
                nj = F[j]
                for k in range(K[i]):
                    for u in range(K[j]):
                        coef = 0.0
                        for w in range(M):
                            if w == f or w == i:
                                continue
                            for p in range(K[f]):
                                coef += q(f, w, u, p)
                        if coef != 0.0:
                            c.add_term(p2(f, nj, u, i, ni + 1, k, m), -coef)
            model.add_constraint(c.eq(0.0))


def _thm3f(model, params, q, p2):
    """THM3f: marginal balance at the finite-capacity queue f (skel:55)."""
    M, F, K, MR, BB, MM, MM1, ZZ, ZM, f = _unpack(params)
    # see _kb/03-api-layer.md for rationale
    i = f
    for ni in range(F[i]):
        c = model.constraint_builder()
        for j in range(M):
            if j == i or j == f:
                continue
            for m in range(MR):
                if BB[m, j] != 0:
                    continue
                for nj in range(1, F[j] + 1):
                    for k in range(K[j]):
                        coef = 0.0
                        for h in range(K[j]):
                            coef += q(j, i, k, h)
                        if coef == 0.0:
                            continue
                        for u in range(K[i]):
                            c.add_term(p2(j, nj, k, i, ni, u, m), coef)
        # RHS is pinned to blocking configuration m = 1 (AMPL) -> 0 here.
        for j in range(M):
            if j == i:
                continue
            for nj in range(F[j] + 1):
                for k in range(K[i]):
                    coef = 0.0
                    for h in range(K[i]):
                        coef += q(i, j, k, h)
                    if coef == 0.0:
                        continue
                    for u in range(K[j]):
                        c.add_term(p2(j, nj, u, i, ni + 1, k, 0), -coef)
        model.add_constraint(c.eq(0.0))


def _thm3i(model, params, q, p2):
    """THM3I: couples blocking configs of successive depth z -> z+1 (skel:56)."""
    M, F, K, MR, BB, MM, MM1, ZZ, ZM, f = _unpack(params)
    for z in range(ZM):
        c = model.constraint_builder()
        for j in range(M):
            if j == f:
                continue
            for m in range(MR):
                if BB[m, j] != 0 or int(ZZ[m]) != z:
                    continue
                for nj in range(1, F[j] + 1):
                    for k in range(K[j]):
                        coef = 0.0
                        for h in range(K[j]):
                            coef += q(j, f, k, h)
                        if coef == 0.0:
                            continue
                        for u in range(K[f]):
                            c.add_term(p2(j, nj, k, f, F[f], u, m), coef)
        for j in range(M):
            if j == f:
                continue
            for m in range(MR):
                if int(ZZ[m]) != z + 1:
                    continue
                for nj in range(F[j] + 1):
                    for k in range(K[f]):
                        coef = 0.0
                        for h in range(K[f]):
                            coef += q(f, j, k, h)
                        if coef == 0.0:
                            continue
                        for u in range(K[j]):
                            c.add_term(p2(j, nj, u, f, F[f], k, m), -coef)
        model.add_constraint(c.eq(0.0))


def _thm3l(model, params, q, p2):
    """THM3L: deepest blocking configurations, ZZ[m] == ZM-1 (skel:57)."""
    M, F, K, MR, BB, MM, MM1, ZZ, ZM, f = _unpack(params)
    for m in range(MR):
        if int(ZZ[m]) != ZM - 1:
            continue
        c = model.constraint_builder()
        for j in range(M):
            if j == f or BB[m, j] != 0 or MM1[m, j] <= 0:
                continue
            for nj in range(1, F[j] + 1):
                for k in range(K[j]):
                    coef = 0.0
                    for h in range(K[j]):
                        coef += q(j, f, k, h)
                    if coef == 0.0:
                        continue
                    for u in range(K[f]):
                        c.add_term(p2(j, nj, k, f, F[f], u, m), coef)
        for j in range(M):
            if j == f or BB[m, j] != 0 or MM1[m, j] <= 0:
                continue
            mtarget = int(MM1[m, j]) - 1  # stored 1-based
            if not (0 <= mtarget < MR):
                continue
            for k in range(K[f]):
                coef = 0.0
                for w in range(M):
                    if w == f:
                        continue
                    for u in range(K[f]):
                        coef += q(f, w, k, u)
                if coef != 0.0:
                    c.add_term(p2(f, F[f], k, f, F[f], k, mtarget), -coef)
        model.add_constraint(c.eq(0.0))


def _add_bound_constraints_bas(model: MapqnLpModel, params: QRBoundsBasParameters) -> None:
    """Add bound constraints for BAS model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K
    MR = params.MR

    # THM4: Queue-length bound inequality
    for j in range(M):
        for kj in range(K[j]):
            for i in range(M):
                for m in range(MR):
                    constraint = model.constraint_builder()
                    # LHS: sum_t sum_ht sum_nj sum_nt nt * p2
                    for t in range(M):
                        for ht in range(K[t]):
                            for njt in range(F[j] + 1):
                                for nt in range(1, F[t] + 1):
                                    constraint.add_term(f'p2_{j}_{njt}_{kj}_{t}_{nt}_{ht}_{m}', float(nt))
                    # RHS: -N * sum
                    for hi in range(K[i]):
                        for njt in range(F[j] + 1):
                            for ni in range(1, F[i] + 1):
                                constraint.add_term(f'p2_{j}_{njt}_{kj}_{i}_{ni}_{hi}_{m}', -float(N))
                    model.add_constraint(constraint.geq(0.0))


__all__ = ['mapqn_qr_bounds_bas']
