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
    _add_bound_constraints_bas(model, params)

    # Build objective: sum of diagonal p2 at target queue across all blocking configs
    target = objective_queue - 1  # 0-based
    objective_terms = {}
    for m in range(MR):
        for ki in range(K[target]):
            for ni in range(1, F[target] + 1):
                var_name = f'p2_{target}_{ni}_{ki}_{target}_{ni}_{ki}_{m}'
                objective_terms[var_name] = 1.0

    if objective_terms:
        first_var = list(objective_terms.keys())[0]
        solution = model.solve(first_var, minimize=(sense == "min"))

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

                                # ZERO6: nj > F[j]
                                if nj > F[j]:
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

    # UEFF: e[i,ki] = sum of p2 where queue i is not blocked
    for i in range(M):
        for ki in range(K[i]):
            constraint = model.constraint_builder()
            constraint.add_term(f'e_{i}_{ki}', -1.0)
            for j in range(M):
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
