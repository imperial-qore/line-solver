"""
QR Bounds with Repetitive-Service Random-Destination (RSRD) Blocking.

Implements quadratic reduction bounds using the RSRD blocking protocol
for MAP queueing networks with finite capacity queues.

This is a port of qrf_rsrd.m MATLAB model.
"""

import numpy as np
from typing import Optional

from .lpmodel import MapqnLpModel
from .parameters import QRBoundsRsrdParameters
from .solution import MapqnSolution


def mapqn_qr_bounds_rsrd(
    params: QRBoundsRsrdParameters,
    objective_queue: int,
    sense: str = "min"
) -> MapqnSolution:
    """
    Solve the QR bounds for RSRD blocking networks.

    Computes bounds on utilization for closed queueing networks with
    finite capacity queues using RSRD blocking protocol.

    Args:
        params: RSRD parameters including:
            - M: Number of queues
            - N: Population
            - F: Capacity per queue [M]
            - K: Number of phases per queue [M]
            - mu: Service rate matrices
            - v: Background transition matrices
            - alpha: Load-dependent rates
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
        Based on qrf_rsrd.m from qrf-revised repository
    """
    params.validate()

    M = params.M
    N = params.N
    F = params.F
    K = params.K

    if not (1 <= objective_queue <= M):
        raise ValueError(f"Objective queue must be in range 1..{M}")
    if sense not in ("min", "max"):
        raise ValueError("Sense must be 'min' or 'max'")

    model = MapqnLpModel()

    # Register all variables
    _register_variables_rsrd(model, params)

    # Add constraints with ZERO bounds applied
    _add_zero_constraints_rsrd(model, params)
    _add_definition_constraints_rsrd(model, params)
    _add_balance_constraints_rsrd(model, params)
    _add_qbal_constraints_rsrd(model, params)
    _add_bound_constraints_rsrd(model, params)

    # Build objective: sum of diagonal p2 at target queue
    target = objective_queue - 1  # 0-based
    objective_terms = {}
    for ki in range(K[target]):
        for ni in range(1, F[target] + 1):
            var_name = f'p2_{target}_{ni}_{ki}_{target}_{ni}_{ki}'
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

        # Update solution with correct objective and derived variables
        variables = dict(solution.variables)

        # Compute U, Ueff, pb for each queue
        for i in range(M):
            total_u = 0.0
            total_ueff = 0.0
            for ki in range(K[i]):
                for ni in range(1, F[i] + 1):
                    u_val = variables.get(f'U_{i}_{ki}_{ni}', 0.0)
                    ueff_val = variables.get(f'Ueff_{i}_{ki}_{ni}', 0.0)
                    total_u += u_val
                    total_ueff += ueff_val
            variables[f'U_{i + 1}'] = total_u
            variables[f'Ueff_{i + 1}'] = total_ueff
            variables[f'pb_{i + 1}'] = total_u - total_ueff

        return MapqnSolution(
            objective_value=obj_value,
            variables=variables
        )
    else:
        return MapqnSolution(objective_value=0.0, variables={})


def _register_variables_rsrd(model: MapqnLpModel, params: QRBoundsRsrdParameters) -> None:
    """Register all variables for RSRD model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K

    # p2 variables: p2[j,nj,kj,i,ni,hi]
    for j in range(M):
        for nj in range(F[j] + 1):
            for kj in range(K[j]):
                for i in range(M):
                    for ni in range(F[i] + 1):
                        for hi in range(K[i]):
                            model.add_variable(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}', lb=0.0, ub=1.0)

    # U and Ueff variables
    for i in range(M):
        for ki in range(K[i]):
            for ni in range(1, F[i] + 1):
                model.add_variable(f'U_{i}_{ki}_{ni}', lb=0.0, ub=1.0)
                model.add_variable(f'Ueff_{i}_{ki}_{ni}', lb=0.0, ub=1.0)

    # pb variables: blocking probability per queue
    # (qrboundsrsrd_skel.mod:51 `var pb{i in 1..M} >=0, <=1`)
    for i in range(M):
        model.add_variable(f'pbvar_{i}', lb=0.0, ub=1.0)


def _add_zero_constraints_rsrd(model: MapqnLpModel, params: QRBoundsRsrdParameters) -> None:
    """Add ZERO constraints for infeasible states."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K

    for j in range(M):
        for nj in range(F[j] + 1):
            for kj in range(K[j]):
                for i in range(M):
                    for ni in range(F[i] + 1):
                        for hi in range(K[i]):
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

                            # ZERO6: i!=j, N-nj-ni > sum of other capacities
                            if i != j:
                                sum_other_f = sum(F[y] for y in range(M) if y != i and y != j)
                                if N - nj - ni > sum_other_f:
                                    is_zero = True

                            if is_zero:
                                constraint = model.constraint_builder()
                                constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}', 1.0)
                                model.add_constraint(constraint.eq(0.0))

        # ZERO7: N-nj > sum of other capacities (for diagonal)
        for kj in range(K[j]):
            sum_other_f = sum(F[y] for y in range(M) if y != j)
            for nj in range(F[j] + 1):
                if N - nj > sum_other_f:
                    constraint = model.constraint_builder()
                    constraint.add_term(f'p2_{j}_{nj}_{kj}_{j}_{nj}_{kj}', 1.0)
                    model.add_constraint(constraint.eq(0.0))


def _add_definition_constraints_rsrd(model: MapqnLpModel, params: QRBoundsRsrdParameters) -> None:
    """Add definition constraints for RSRD model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K
    r = params.r

    # ONE: Normalization - sum of diagonal p2 = 1 for each queue
    for j in range(M):
        constraint = model.constraint_builder()
        for nj in range(F[j] + 1):
            for kj in range(K[j]):
                constraint.add_term(f'p2_{j}_{nj}_{kj}_{j}_{nj}_{kj}', 1.0)
        model.add_constraint(constraint.eq(1.0))

    # SYMMETRY: p2[j,nj,k,i,ni,h] = p2[i,ni,h,j,nj,k]
    for j in range(M):
        for nj in range(F[j] + 1):
            for kj in range(K[j]):
                for i in range(j + 1, M):
                    for ni in range(F[i] + 1):
                        for hi in range(K[i]):
                            constraint = model.constraint_builder()
                            constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}', 1.0)
                            constraint.add_term(f'p2_{i}_{ni}_{hi}_{j}_{nj}_{kj}', -1.0)
                            model.add_constraint(constraint.eq(0.0))

    # MARGINALS: p2[j,nj,k,j,nj,k] = sum{ni,h} p2[j,nj,k,i,ni,h] for each i!=j
    for j in range(M):
        for kj in range(K[j]):
            for nj in range(F[j] + 1):
                for i in range(M):
                    if i == j:
                        continue
                    constraint = model.constraint_builder()
                    constraint.add_term(f'p2_{j}_{nj}_{kj}_{j}_{nj}_{kj}', 1.0)
                    for ni in range(F[i] + 1):
                        for hi in range(K[i]):
                            constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}', -1.0)
                    model.add_constraint(constraint.eq(0.0))

    # UCLASSIC: U[i,k,ni] = p2[i,ni,k,i,ni,k]
    for i in range(M):
        for ki in range(K[i]):
            for ni in range(1, F[i] + 1):
                constraint = model.constraint_builder()
                constraint.add_term(f'U_{i}_{ki}_{ni}', 1.0)
                constraint.add_term(f'p2_{i}_{ni}_{ki}_{i}_{ni}_{ki}', -1.0)
                model.add_constraint(constraint.eq(0.0))

    # UEFFS: Ueff[i,k,ni] = U[i,k,ni] - sum{j: r[i,j]>0} r[i,j]*p2[i,ni,k,j,F[j],h]
    for i in range(M):
        for ki in range(K[i]):
            for ni in range(1, F[i] + 1):
                constraint = model.constraint_builder()
                constraint.add_term(f'Ueff_{i}_{ki}_{ni}', 1.0)
                constraint.add_term(f'p2_{i}_{ni}_{ki}_{i}_{ni}_{ki}', -1.0)
                for j in range(M):
                    if j != i and r[i, j] > 0:
                        for hj in range(K[j]):
                            constraint.add_term(f'p2_{i}_{ni}_{ki}_{j}_{F[j]}_{hj}', r[i, j])
                model.add_constraint(constraint.eq(0.0))

    # PBLOCK: pb[i] = sum_{k, ni>=1} (U[i,k,ni] - Ueff[i,k,ni])
    # (qrboundsrsrd_skel.mod:52)
    for i in range(M):
        constraint = model.constraint_builder()
        constraint.add_term(f'pbvar_{i}', 1.0)
        for ki in range(K[i]):
            for ni in range(1, F[i] + 1):
                constraint.add_term(f'U_{i}_{ki}_{ni}', -1.0)
                constraint.add_term(f'Ueff_{i}_{ki}_{ni}', 1.0)
        model.add_constraint(constraint.eq(0.0))

    # PBB: pb[i] <= sum_{j!=i, r[i,j]>0} sum_h p2[j,F[j],h,j,F[j],h]
    # (qrboundsrsrd_skel.mod:53)
    for i in range(M):
        constraint = model.constraint_builder()
        constraint.add_term(f'pbvar_{i}', 1.0)
        for j in range(M):
            if j != i and r[i, j] > 0:
                for hj in range(K[j]):
                    constraint.add_term(f'p2_{j}_{F[j]}_{hj}_{j}_{F[j]}_{hj}', -1.0)
        model.add_constraint(constraint.leq(0.0))


def _add_balance_constraints_rsrd(model: MapqnLpModel, params: QRBoundsRsrdParameters) -> None:
    """Add balance constraints for RSRD model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K

    # Compute q values
    def q(i: int, j: int, ki: int, hi: int, n: int) -> float:
        return params.q(i, j, ki, hi, n)

    # see _kb/03-api-layer.md for rationale
    for j in range(M):
        for kj in range(K[j]):
            constraint = model.constraint_builder()
            for nj in range(1, F[j] + 1):
                # RHS: -N * p2[j,nj,k,j,nj,k]
                constraint.add_term(f'p2_{j}_{nj}_{kj}_{j}_{nj}_{kj}', -float(N))
                # LHS: sum over all i (including i == j)
                for i in range(M):
                    for ni in range(1, F[i] + 1):
                        for hi in range(K[i]):
                            constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}', float(ni))
            model.add_constraint(constraint.eq(0.0))

    # see _kb/03-api-layer.md for rationale
    for j in range(M):
        for kj in range(K[j]):
            constraint = model.constraint_builder()
            constraint.add_term(f'p2_{j}_0_{kj}_{j}_0_{kj}', -float(N))
            for i in range(M):
                if i == j:
                    continue
                for ni in range(1, F[i] + 1):
                    for hi in range(K[i]):
                        constraint.add_term(f'p2_{j}_0_{kj}_{i}_{ni}_{hi}', float(ni))
            model.add_constraint(constraint.eq(0.0))

    # THM2: Phase balance
    for i in range(M):
        for ki in range(K[i]):
            constraint = model.constraint_builder()
            for ni in range(1, F[i] + 1):
                # Terms with Ueff (j!=i)
                for j in range(M):
                    if j == i:
                        continue
                    for hi in range(K[i]):
                        if hi == ki:
                            continue
                        coef = q(i, j, ki, hi, ni)
                        constraint.add_term(f'Ueff_{i}_{ki}_{ni}', coef)
                # Terms with U (j==i, h!=k)
                for hi in range(K[i]):
                    if hi == ki:
                        continue
                    coef = q(i, i, ki, hi, ni)
                    constraint.add_term(f'p2_{i}_{ni}_{ki}_{i}_{ni}_{ki}', coef)

            # RHS terms (subtract)
            for ni in range(1, F[i] + 1):
                # Terms with Ueff (j!=i)
                for j in range(M):
                    if j == i:
                        continue
                    for hi in range(K[i]):
                        if hi == ki:
                            continue
                        coef = q(i, j, hi, ki, ni)
                        constraint.add_term(f'Ueff_{i}_{hi}_{ni}', -coef)
                # Terms with U (j==i, h!=k)
                for hi in range(K[i]):
                    if hi == ki:
                        continue
                    coef = q(i, i, hi, ki, ni)
                    constraint.add_term(f'p2_{i}_{ni}_{hi}_{i}_{ni}_{hi}', -coef)
            model.add_constraint(constraint.eq(0.0))

    # THM3a: Marginal balance for ni in 1:F[i]-1
    for i in range(M):
        for ni in range(1, F[i]):
            constraint = model.constraint_builder()
            # LHS: arrivals to queue i at level ni
            for j in range(M):
                if j == i:
                    continue
                for kj in range(K[j]):
                    for hj in range(K[j]):
                        for ui in range(K[i]):
                            for nj in range(1, F[j] + 1):
                                constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{ui}', q(j, i, kj, hj, nj))
            # RHS: departures from queue i at level ni+1
            for j in range(M):
                if j == i:
                    continue
                for ki in range(K[i]):
                    for hi in range(K[i]):
                        for uj in range(K[j]):
                            for nj in range(F[j]):
                                if ni + 1 <= F[i]:
                                    constraint.add_term(f'p2_{i}_{ni + 1}_{ki}_{j}_{nj}_{uj}', -q(i, j, ki, hi, ni + 1))
            model.add_constraint(constraint.eq(0.0))

    # THM3b: Marginal balance for ni=0, per phase
    for i in range(M):
        for ui in range(K[i]):
            constraint = model.constraint_builder()
            # LHS: arrivals to queue i at ni=0
            for j in range(M):
                if j == i:
                    continue
                for kj in range(K[j]):
                    for hj in range(K[j]):
                        for nj in range(1, F[j] + 1):
                            constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_0_{ui}', q(j, i, kj, hj, nj))
            # RHS: departures from queue i at ni=1
            for j in range(M):
                if j == i:
                    continue
                for ki in range(K[i]):
                    for nj in range(F[j]):
                        for hj in range(K[j]):
                            constraint.add_term(f'p2_{i}_1_{ki}_{j}_{nj}_{hj}', -q(i, j, ki, ui, 1))
            model.add_constraint(constraint.eq(0.0))


def _add_qbal_constraints_rsrd(model: MapqnLpModel, params: QRBoundsRsrdParameters) -> None:
    """QBAL: queue balance constraint (qrboundsrsrd_skel.mod:64).

    Ported from matlab/lib/qrf/qrf_rsrd.m. This block was entirely absent from
    the native-Python port; a missing constraint enlarges the polytope and so
    loosens the bound (on the paper's M=5,N=20 RS-RD example its absence gave
    U1min=0.87001 against the published 0.87058).
    """
    M = params.M
    N = params.N
    F = params.F
    K = params.K

    def q(i, j, ki, hi, n):
        return params.q(i, j, ki, hi, n)

    for i in range(M):
        for ki in range(K[i]):
            constraint = model.constraint_builder()

            # LHS 1: sum_{h!=k, j!=i, ni>=1, u, nj<=F[j]-1} q[i,j,k,h,ni]*ni*p2[i,ni,k,j,nj,u]
            for hi in range(K[i]):
                if hi == ki:
                    continue
                for j in range(M):
                    if j == i:
                        continue
                    for ni in range(1, F[i] + 1):
                        coef = q(i, j, ki, hi, ni) * ni
                        if coef == 0.0:
                            continue
                        for uj in range(K[j]):
                            for nj in range(0, F[j]):
                                constraint.add_term(
                                    f'p2_{i}_{ni}_{ki}_{j}_{nj}_{uj}', coef)

            # LHS 2: sum_{h!=k, ni>=1} q[i,i,k,h,ni]*ni*p2[i,ni,k,i,ni,k]
            for hi in range(K[i]):
                if hi == ki:
                    continue
                for ni in range(1, F[i] + 1):
                    coef = q(i, i, ki, hi, ni) * ni
                    if coef == 0.0:
                        continue
                    constraint.add_term(f'p2_{i}_{ni}_{ki}_{i}_{ni}_{ki}', coef)

            # LHS 3: sum_{j!=i, h, ni>=1, u, nj<=min(F[j]-1, N-ni)} q[i,j,h,k,ni]*p2[i,ni,h,j,nj,u]
            for j in range(M):
                if j == i:
                    continue
                for hi in range(K[i]):
                    for ni in range(1, F[i] + 1):
                        coef = q(i, j, hi, ki, ni)
                        if coef == 0.0:
                            continue
                        max_nj = min(F[j] - 1, N - ni)
                        for uj in range(K[j]):
                            for nj in range(0, max_nj + 1):
                                constraint.add_term(
                                    f'p2_{i}_{ni}_{hi}_{j}_{nj}_{uj}', coef)

            # RHS 1: sum_{j!=i, h, ni<=F[i]-1, u, nj>=1} q[j,i,h,u,nj]*p2[i,ni,k,j,nj,h]
            for j in range(M):
                if j == i:
                    continue
                for hj in range(K[j]):
                    for uj in range(K[j]):
                        for nj in range(1, F[j] + 1):
                            coef = q(j, i, hj, uj, nj)
                            if coef == 0.0:
                                continue
                            for ni in range(0, F[i]):
                                constraint.add_term(
                                    f'p2_{i}_{ni}_{ki}_{j}_{nj}_{hj}', -coef)

            # RHS 2: sum_{h!=k, ni>=1} q[i,i,h,k,ni]*ni*p2[i,ni,h,i,ni,h]
            for hi in range(K[i]):
                if hi == ki:
                    continue
                for ni in range(1, F[i] + 1):
                    coef = q(i, i, hi, ki, ni) * ni
                    if coef == 0.0:
                        continue
                    constraint.add_term(f'p2_{i}_{ni}_{hi}_{i}_{ni}_{hi}', -coef)

            # RHS 3: sum_{h!=k, j!=i, ni>=1, u, nj<=F[j]-1} q[i,j,h,k,ni]*ni*p2[i,ni,h,j,nj,u]
            for hi in range(K[i]):
                if hi == ki:
                    continue
                for j in range(M):
                    if j == i:
                        continue
                    for ni in range(1, F[i] + 1):
                        coef = q(i, j, hi, ki, ni) * ni
                        if coef == 0.0:
                            continue
                        for uj in range(K[j]):
                            for nj in range(0, F[j]):
                                constraint.add_term(
                                    f'p2_{i}_{ni}_{hi}_{j}_{nj}_{uj}', -coef)

            model.add_constraint(constraint.eq(0.0))


def _add_bound_constraints_rsrd(model: MapqnLpModel, params: QRBoundsRsrdParameters) -> None:
    """Add bound constraints for RSRD model."""
    M = params.M
    N = params.N
    F = params.F
    K = params.K

    # THM4: Bound constraint
    for j in range(M):
        for kj in range(K[j]):
            for i in range(M):
                constraint = model.constraint_builder()
                # LHS: sum{t,h,nj,nt} nt*p2[j,nj,k,t,nt,h]
                for t in range(M):
                    for ht in range(K[t]):
                        for nj in range(F[j] + 1):
                            for nt in range(1, F[t] + 1):
                                constraint.add_term(f'p2_{j}_{nj}_{kj}_{t}_{nt}_{ht}', float(nt))
                # RHS (negated): -N*sum{h,nj,ni} p2[j,nj,k,i,ni,h]
                for hi in range(K[i]):
                    for nj in range(F[j] + 1):
                        for ni in range(1, F[i] + 1):
                            constraint.add_term(f'p2_{j}_{nj}_{kj}_{i}_{ni}_{hi}', -float(N))
                model.add_constraint(constraint.geq(0.0))


__all__ = ['mapqn_qr_bounds_rsrd']
