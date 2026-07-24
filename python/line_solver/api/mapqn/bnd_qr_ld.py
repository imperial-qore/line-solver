"""
Quadratic Reduction Bounds for Load-Dependent Systems.

Implements quadratic reduction bounds for load-dependent MAP queueing
networks. Provides tighter performance bounds through quadratic
approximation methods for systems with load-dependent service rates.

This is a port of bnd_quadraticreduction_ld.mod AMPL model.
"""

import numpy as np
from typing import Optional

from .lpmodel import MapqnLpModel
from .parameters import QuadraticLDParameters
from .solution import MapqnSolution


def mapqn_bnd_qr_ld(
    params: QuadraticLDParameters,
    objective_queue: int,
    objective_phase: int,
    objective_n: int
) -> MapqnSolution:
    """
    Solve the quadratic reduction bound for load-dependent systems.

    Computes bounds on marginal probabilities for closed queueing networks
    with load-dependent service rates using p2 LP formulation.

    Args:
        params: Quadratic LD parameters including:
            - M: Number of queues
            - N: Population
            - K: Number of phases per queue [M]
            - mu: Service rate matrices
            - v: Background transition matrices
            - alpha: Load-dependent rates [M x N]
            - r: Routing matrix
        objective_queue: Queue index for objective (1-based).
        objective_phase: Phase index for objective (1-based).
        objective_n: Population level for objective (0 to N).

    Returns:
        MapqnSolution containing:
            - objective_value: Maximum marginal probability
            - variables: All LP variable values

    Raises:
        ValueError: If indices are out of range.
        RuntimeError: If LP is infeasible or unbounded.

    Reference:
        Based on AMPL model bnd_quadraticreduction_ld.mod
    """
    params.validate()

    M = params.M
    N = params.N
    K = params.K

    if not (1 <= objective_queue <= M):
        raise ValueError(f"Objective queue must be in range 1..{M}")
    if not (1 <= objective_phase <= K[objective_queue - 1]):
        raise ValueError(f"Objective phase must be in range 1..{K[objective_queue - 1]}")
    if not (0 <= objective_n <= N):
        raise ValueError(f"Objective N must be in range 0..{N}")

    model = MapqnLpModel()

    # Register all variables (p2 only)
    _register_variables_ld(model, M, N, K)

    # Add all constraints
    _add_definition_constraints_ld(model, params)
    _add_littles_law_constraints_ld(model, params)
    _add_balance_constraints_ld(model, params)
    _add_bound_constraints_ld(model, params)

    # Solve: maximize p2[objective_queue, objective_n, objective_phase, ...]
    objective_var = f'p2_{objective_queue}_{objective_n}_{objective_phase}_{objective_queue}_{objective_n}_{objective_phase}'
    solution = model.solve(objective_var, minimize=False)

    return solution


def _register_variables_ld(model: MapqnLpModel, M: int, N: int, K: np.ndarray) -> None:
    """Register p2 variables for load-dependent model."""
    for j in range(1, M + 1):
        for nj in range(N + 1):
            for k in range(1, K[j - 1] + 1):
                for i in range(1, M + 1):
                    for ni in range(N + 1):
                        for h in range(1, K[i - 1] + 1):
                            model.add_variable(f'p2_{j}_{nj}_{k}_{i}_{ni}_{h}', lb=0.0, ub=1.0)


def _add_definition_constraints_ld(model: MapqnLpModel, params: QuadraticLDParameters) -> None:
    """Add definition constraints for load-dependent model."""
    M = params.M
    N = params.N
    K = params.K

    # ONE: sum{nj,k} p2[j,nj,k,j,nj,k] = 1 for each j
    for j in range(1, M + 1):
        constraint = model.constraint_builder()
        for nj in range(N + 1):
            for k in range(1, K[j - 1] + 1):
                constraint.add_term(f'p2_{j}_{nj}_{k}_{j}_{nj}_{k}', 1.0)
        model.add_constraint(constraint.eq(1.0))

    # ZERO1: p2[j,nj,k,i,ni,h] = 0 when i==j, nj==ni, h!=k
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for nj in range(N + 1):
                for h in range(1, K[j - 1] + 1):
                    if h != k:
                        constraint = model.constraint_builder()
                        constraint.add_term(f'p2_{j}_{nj}_{k}_{j}_{nj}_{h}', 1.0)
                        model.add_constraint(constraint.eq(0.0))

    # ZERO2: p2[j,nj,k,i,ni,h] = 0 when i==j, nj!=ni
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for nj in range(N + 1):
                for ni in range(N + 1):
                    if nj != ni:
                        for h in range(1, K[j - 1] + 1):
                            constraint = model.constraint_builder()
                            constraint.add_term(f'p2_{j}_{nj}_{k}_{j}_{ni}_{h}', 1.0)
                            model.add_constraint(constraint.eq(0.0))

    # ZERO3: p2[j,nj,k,i,ni,h] = 0 when i!=j, nj+ni > N
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for nj in range(N + 1):
                for i in range(1, M + 1):
                    if i != j:
                        for ni in range(N + 1):
                            if nj + ni > N:
                                for h in range(1, K[i - 1] + 1):
                                    constraint = model.constraint_builder()
                                    constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_{ni}_{h}', 1.0)
                                    model.add_constraint(constraint.eq(0.0))

    # SYMMETRY: p2[i,ni,h,j,nj,k] = p2[j,nj,k,i,ni,h]
    for j in range(1, M + 1):
        for nj in range(N + 1):
            for k in range(1, K[j - 1] + 1):
                for i in range(1, M + 1):
                    for ni in range(N + 1):
                        for h in range(1, K[i - 1] + 1):
                            constraint = model.constraint_builder()
                            constraint.add_term(f'p2_{i}_{ni}_{h}_{j}_{nj}_{k}', 1.0)
                            constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_{ni}_{h}', -1.0)
                            model.add_constraint(constraint.eq(0.0))

    # MARGINALS: p2[j,nj,k,j,nj,k] = sum{ni,h} p2[j,nj,k,i,ni,h] for i!=j
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for nj in range(N + 1):
                for i in range(1, M + 1):
                    if i != j:
                        constraint = model.constraint_builder()
                        constraint.add_term(f'p2_{j}_{nj}_{k}_{j}_{nj}_{k}', 1.0)
                        for ni in range(max(0, N - nj) + 1):
                            for h in range(1, K[i - 1] + 1):
                                constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_{ni}_{h}', -1.0)
                        model.add_constraint(constraint.eq(0.0))


def _add_littles_law_constraints_ld(model: MapqnLpModel, params: QuadraticLDParameters) -> None:
    """Add Little's law constraints for load-dependent model."""
    M = params.M
    N = params.N
    K = params.K

    # THM1: Queue length theorem
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            constraint = model.constraint_builder()
            # Left side: sum ni*p2[j,nj,k,i,ni,h]
            for i in range(1, M + 1):
                for nj in range(1, N + 1):
                    for ni in range(1, N + 1):
                        for h in range(1, K[i - 1] + 1):
                            constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_{ni}_{h}', float(ni))
            # Right side: N*sum p2[j,nj,k,j,nj,k]
            for nj in range(1, N + 1):
                constraint.add_term(f'p2_{j}_{nj}_{k}_{j}_{nj}_{k}', -float(N))
            model.add_constraint(constraint.eq(0.0))

    # THM1c: Empty queue theorem
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            constraint = model.constraint_builder()
            # Left side
            for i in range(1, M + 1):
                for ni in range(1, N + 1):
                    for h in range(1, K[i - 1] + 1):
                        constraint.add_term(f'p2_{j}_0_{k}_{i}_{ni}_{h}', float(ni))
            # Right side
            constraint.add_term(f'p2_{j}_0_{k}_{j}_0_{k}', -float(N))
            model.add_constraint(constraint.eq(0.0))

    # PC2: sum nj*ni*p2[j,nj,k,i,ni,h] = N*N
    constraint = model.constraint_builder()
    for i in range(1, M + 1):
        for j in range(1, M + 1):
            for ni in range(1, N + 1):
                for nj in range(1, N + 1):
                    for h in range(1, K[i - 1] + 1):
                        for k in range(1, K[j - 1] + 1):
                            constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_{ni}_{h}', float(nj * ni))
    model.add_constraint(constraint.eq(float(N * N)))


def _add_balance_constraints_ld(model: MapqnLpModel, params: QuadraticLDParameters) -> None:
    """Add balance constraints for load-dependent model."""
    M = params.M
    N = params.N
    K = params.K

    # THM2: Phase balance
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            constraint = model.constraint_builder()
            for j in range(1, M + 1):
                for h in range(1, K[i - 1] + 1):
                    if not (h == k and i == j):
                        for ni in range(1, N + 1):
                            q_out = params.q(i - 1, j - 1, k - 1, h - 1, ni)
                            q_in = params.q(i - 1, j - 1, h - 1, k - 1, ni)
                            constraint.add_term(f'p2_{i}_{ni}_{k}_{i}_{ni}_{k}', q_out)
                            constraint.add_term(f'p2_{i}_{ni}_{h}_{i}_{ni}_{h}', -q_in)
            model.add_constraint(constraint.eq(0.0))

    # THM3a: Flow balance for ni in 1..N-1
    for i in range(1, M + 1):
        for ni in range(1, N):
            constraint = model.constraint_builder()
            # Incoming flow
            for j in range(1, M + 1):
                if j != i:
                    for k in range(1, K[j - 1] + 1):
                        for h in range(1, K[j - 1] + 1):
                            for u in range(1, K[i - 1] + 1):
                                for nj in range(1, N - ni + 1):
                                    q = params.q(j - 1, i - 1, k - 1, h - 1, nj)
                                    constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_{ni}_{u}', q)
            # Outgoing flow
            for j in range(1, M + 1):
                if j != i:
                    for k in range(1, K[i - 1] + 1):
                        for h in range(1, K[i - 1] + 1):
                            q = params.q(i - 1, j - 1, k - 1, h - 1, ni + 1)
                            constraint.add_term(f'p2_{i}_{ni + 1}_{k}_{i}_{ni + 1}_{k}', -q)
            model.add_constraint(constraint.eq(0.0))

    # THM3b: Flow balance for ni=0
    for i in range(1, M + 1):
        for u in range(1, K[i - 1] + 1):
            constraint = model.constraint_builder()
            # Incoming to empty queue
            for j in range(1, M + 1):
                if j != i:
                    for k in range(1, K[j - 1] + 1):
                        for h in range(1, K[j - 1] + 1):
                            for nj in range(1, N + 1):
                                q = params.q(j - 1, i - 1, k - 1, h - 1, nj)
                                constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_0_{u}', q)
            # Outgoing from single customer
            for j in range(1, M + 1):
                if j != i:
                    for k in range(1, K[i - 1] + 1):
                        q = params.q(i - 1, j - 1, k - 1, u - 1, 1)
                        constraint.add_term(f'p2_{i}_1_{k}_{i}_1_{k}', -q)
            model.add_constraint(constraint.eq(0.0))

    # QBAL: Queue balance
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            constraint = model.constraint_builder()
            # Outgoing from phase k to other phases
            for h in range(1, K[i - 1] + 1):
                if h != k:
                    for j in range(1, M + 1):
                        for ni in range(1, N + 1):
                            q = params.q(i - 1, j - 1, k - 1, h - 1, ni)
                            constraint.add_term(f'p2_{i}_{ni}_{k}_{i}_{ni}_{k}', q * ni)
            # Incoming from other phases
            for j in range(1, M + 1):
                if j != i:
                    for h in range(1, K[i - 1] + 1):
                        for ni in range(1, N + 1):
                            q = params.q(i - 1, j - 1, h - 1, k - 1, ni)
                            constraint.add_term(f'p2_{i}_{ni}_{h}_{i}_{ni}_{h}', q)
            model.add_constraint(constraint.eq(0.0))


def _add_bound_constraints_ld(model: MapqnLpModel, params: QuadraticLDParameters) -> None:
    """Add bound constraints for load-dependent model."""
    M = params.M
    N = params.N
    K = params.K

    # THM4: QMIN-like constraint
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                constraint = model.constraint_builder()
                # Left side: total queue length
                for t in range(1, M + 1):
                    for h in range(1, K[t - 1] + 1):
                        for nj in range(N + 1):
                            for nt in range(N + 1):
                                constraint.add_term(f'p2_{j}_{nj}_{k}_{t}_{nt}_{h}', float(nt))
                # Right side: N times probability mass at queue i
                for h in range(1, K[i - 1] + 1):
                    for nj in range(N + 1):
                        for ni in range(N + 1):
                            constraint.add_term(f'p2_{j}_{nj}_{k}_{i}_{ni}_{h}', -float(N))
                model.add_constraint(constraint.geq(0.0))


__all__ = ['mapqn_bnd_qr_ld']
