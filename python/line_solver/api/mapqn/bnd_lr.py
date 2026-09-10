"""
General Linear Reduction Bounds for MAP Queueing Networks.

Implements general linear reduction methods for computing performance
bounds in MAP queueing networks with phase-type service processes.

This is a port of bnd_linearreduction_new.mod AMPL model.
"""

import numpy as np
from typing import Optional

from .lpmodel import MapqnLpModel
from .parameters import LinearReductionParameters
from .solution import MapqnSolution


def mapqn_bnd_lr(
    params: LinearReductionParameters,
    objective_queue: int,
    objective_phase: int,
    sense: str = 'max',
    objective_kind: str = 'U'
) -> MapqnSolution:
    """
    Solve the general linear reduction bound for MAP queueing networks.

    Computes bounds on utilization for closed queueing networks with
    phase-type service distributions using a linear programming formulation.

    Args:
        params: Linear reduction parameters including:
            - M: Number of queues
            - N: Population
            - K: Number of phases per queue [M]
            - mu: Service rate matrices [M] list of K[i] x K[i]
            - r: Routing matrix [M x M]
            - v: Background transition matrices [M] list of K[i] x K[i]
        objective_queue: Queue index to maximize utilization (1-based).
        objective_phase: Phase index to maximize utilization (1-based).

    Returns:
        MapqnSolution containing:
            - objective_value: Maximum utilization bound
            - variables: All LP variable values

    Raises:
        ValueError: If indices are out of range.
        RuntimeError: If LP is infeasible or unbounded.

    Reference:
        Based on AMPL model bnd_linearreduction_new.mod
    """
    params.validate()

    M = params.M
    N = params.N
    K = params.K

    if not (1 <= objective_queue <= M):
        raise ValueError(f"Objective queue must be in range 1..{M}")
    if not (1 <= objective_phase <= K[objective_queue - 1]):
        raise ValueError(f"Objective phase must be in range 1..{K[objective_queue - 1]}")

    model = MapqnLpModel()

    # Register all variables
    _register_variables_lr(model, M, N, K)

    # Add all constraints
    _add_definition_constraints_lr(model, params)
    _add_mean_indices_constraints_lr(model, params)
    _add_balance_constraints_lr(model, params)
    _add_bound_constraints_lr(model, params)

    # see _kb/03-api-layer.md for rationale
    if sense not in ('min', 'max'):
        raise ValueError("sense must be 'min' or 'max', got %r" % (sense,))
    if objective_kind not in ('U', 'Q'):
        raise ValueError("objective_kind must be 'U' or 'Q', got %r"
                         % (objective_kind,))
    objective_var = f'{objective_kind}_{objective_queue}_{objective_phase}'
    solution = model.solve(objective_var, minimize=(sense == 'min'))

    return solution


def _register_variables_lr(model: MapqnLpModel, M: int, N: int, K: np.ndarray) -> None:
    """Register all variables for linear reduction model."""

    # U variables: utilization at each queue-phase
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            model.add_variable(f'U_{i}_{k}', lb=0.0, ub=1.0)

    # IT variables: idle time at each queue-phase
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            model.add_variable(f'IT_{i}_{k}', lb=0.0, ub=1.0)

    # Q variables: mean queue length
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            model.add_variable(f'Q_{i}_{k}', lb=0.0, ub=float(N))

    # UP variables: utilization products
    for j in range(1, M + 1):
        for kj in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                for hi in range(1, K[i - 1] + 1):
                    model.add_variable(f'UP_{j}_{kj}_{i}_{hi}', lb=0.0, ub=1.0)

    # QP variables: queue-length products
    for j in range(1, M + 1):
        for kj in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                for hi in range(1, K[i - 1] + 1):
                    model.add_variable(f'QP_{j}_{kj}_{i}_{hi}', lb=0.0, ub=float(N))

    # C variables: conditional queue lengths
    for j in range(1, M + 1):
        for kj in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                model.add_variable(f'C_{j}_{kj}_{i}', lb=0.0, ub=float(N))

    # I variables: conditional idle lengths
    for j in range(1, M + 1):
        for kj in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                model.add_variable(f'I_{j}_{kj}_{i}', lb=0.0, ub=float(N))

    # p1 variables: marginal probabilities
    for j in range(1, M + 1):
        for kj in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                for ni in range(N + 1):
                    for hi in range(1, K[i - 1] + 1):
                        model.add_variable(f'p1_{j}_{kj}_{i}_{ni}_{hi}', lb=0.0, ub=1.0)

    # p1c variables: complementary marginal probabilities
    for j in range(1, M + 1):
        for kj in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                for ni in range(N + 1):
                    for hi in range(1, K[i - 1] + 1):
                        model.add_variable(f'p1c_{j}_{kj}_{i}_{ni}_{hi}', lb=0.0, ub=1.0)


def _add_definition_constraints_lr(model: MapqnLpModel, params: LinearReductionParameters) -> None:
    """Add definition constraints for linear reduction model."""
    M = params.M
    N = params.N
    K = params.K

    # ZER1: p1[j,k,j,0,k] = 0
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            constraint = model.constraint_builder()
            constraint.add_term(f'p1_{j}_{k}_{j}_0_{k}', 1.0)
            model.add_constraint(constraint.eq(0.0))

    # ZER2: p1[j,k,j,nj,h] = 0 for h != k
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for nj in range(N + 1):
                for h in range(1, K[j - 1] + 1):
                    if h != k:
                        constraint = model.constraint_builder()
                        constraint.add_term(f'p1_{j}_{k}_{j}_{nj}_{h}', 1.0)
                        model.add_constraint(constraint.eq(0.0))

    # ZER3: p1[j,k,i,N,h] = 0 for j != i
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                if j != i:
                    for h in range(1, K[i - 1] + 1):
                        constraint = model.constraint_builder()
                        constraint.add_term(f'p1_{j}_{k}_{i}_{N}_{h}', 1.0)
                        model.add_constraint(constraint.eq(0.0))

    # ZER4: p1c[j,k,j,nj,h] = 0 for nj >= 1
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for nj in range(1, N + 1):
                for h in range(1, K[j - 1] + 1):
                    constraint = model.constraint_builder()
                    constraint.add_term(f'p1c_{j}_{k}_{j}_{nj}_{h}', 1.0)
                    model.add_constraint(constraint.eq(0.0))

    # CEQU: C[j,k,j] = Q[j,k]
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            constraint = model.constraint_builder()
            constraint.add_term(f'C_{j}_{k}_{j}', 1.0)
            constraint.add_term(f'Q_{j}_{k}', -1.0)
            model.add_constraint(constraint.eq(0.0))

    # ONE1: sum of all p1 and p1c = 1 for each (j, i) pair
    for j in range(1, M + 1):
        for i in range(1, M + 1):
            constraint = model.constraint_builder()
            for kj in range(1, K[j - 1] + 1):
                for hi in range(1, K[i - 1] + 1):
                    for ni in range(N + 1):
                        constraint.add_term(f'p1_{j}_{kj}_{i}_{ni}_{hi}', 1.0)
                        constraint.add_term(f'p1c_{j}_{kj}_{i}_{ni}_{hi}', 1.0)
            model.add_constraint(constraint.eq(1.0))


def _add_mean_indices_constraints_lr(model: MapqnLpModel, params: LinearReductionParameters) -> None:
    """Add mean indices constraints for linear reduction model."""
    M = params.M
    N = params.N
    K = params.K

    # UTLB: U[i,k] = sum{t,nt,h} p1[i,k,t,nt,h]
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            for t in range(1, M + 1):
                constraint = model.constraint_builder()
                constraint.add_term(f'U_{i}_{k}', 1.0)
                for nt in range(N + 1):
                    for h in range(1, K[t - 1] + 1):
                        constraint.add_term(f'p1_{i}_{k}_{t}_{nt}_{h}', -1.0)
                model.add_constraint(constraint.eq(0.0))

    # UTLC: IT[i,k] = sum{t,nt,h} p1c[i,k,t,nt,h]
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            for t in range(1, M + 1):
                constraint = model.constraint_builder()
                constraint.add_term(f'IT_{i}_{k}', 1.0)
                for nt in range(N + 1):
                    for h in range(1, K[t - 1] + 1):
                        constraint.add_term(f'p1c_{i}_{k}_{t}_{nt}_{h}', -1.0)
                model.add_constraint(constraint.eq(0.0))

    # QLEN: Q[i,k] = sum{ni} ni * p1[i,k,i,ni,k]
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            constraint = model.constraint_builder()
            constraint.add_term(f'Q_{i}_{k}', 1.0)
            for ni in range(N + 1):
                constraint.add_term(f'p1_{i}_{k}_{i}_{ni}_{k}', -float(ni))
            model.add_constraint(constraint.eq(0.0))

    # CLEN: C[j,k,i] = sum{ni,h} ni * p1[j,k,i,ni,h]
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                constraint = model.constraint_builder()
                constraint.add_term(f'C_{j}_{k}_{i}', 1.0)
                for ni in range(N + 1):
                    for h in range(1, K[i - 1] + 1):
                        constraint.add_term(f'p1_{j}_{k}_{i}_{ni}_{h}', -float(ni))
                model.add_constraint(constraint.eq(0.0))


def _add_balance_constraints_lr(model: MapqnLpModel, params: LinearReductionParameters) -> None:
    """Add balance constraints for linear reduction model."""
    M = params.M
    N = params.N
    K = params.K

    # ONE: sum{k} (U[j,k] + IT[j,k]) = 1
    for j in range(1, M + 1):
        constraint = model.constraint_builder()
        for k in range(1, K[j - 1] + 1):
            constraint.add_term(f'U_{j}_{k}', 1.0)
            constraint.add_term(f'IT_{j}_{k}', 1.0)
        model.add_constraint(constraint.eq(1.0))

    # POPC: sum{i,k} Q[i,k] = N
    constraint = model.constraint_builder()
    for i in range(1, M + 1):
        for k in range(1, K[i - 1] + 1):
            constraint.add_term(f'Q_{i}_{k}', 1.0)
    model.add_constraint(constraint.eq(float(N)))

    # MPCB: sum{i} C[j,k,i] = N * U[j,k]
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            constraint = model.constraint_builder()
            for i in range(1, M + 1):
                constraint.add_term(f'C_{j}_{k}_{i}', 1.0)
            constraint.add_term(f'U_{j}_{k}', -float(N))
            model.add_constraint(constraint.eq(0.0))

    # see _kb/03-api-layer.md for rationale
    for i in range(1, M + 1):
        if K[i - 1] < 2:
            continue
        for k in range(1, K[i - 1] + 1):
            constraint = model.constraint_builder()
            for j in range(1, M + 1):
                for h in range(1, K[i - 1] + 1):
                    constraint.add_term(
                        f'U_{i}_{k}', params.q(i - 1, j - 1, k - 1, h - 1))
                    constraint.add_term(
                        f'U_{i}_{h}', -params.q(i - 1, j - 1, h - 1, k - 1))
            model.add_constraint(constraint.eq(0.0))

    # see _kb/03-api-layer.md for rationale
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            for i in range(1, M + 1):
                for h in range(1, K[i - 1] + 1):
                    if j < i or (j == i and k < h):
                        constraint = model.constraint_builder()
                        for ni in range(1, N + 1):
                            constraint.add_term(f'p1_{j}_{k}_{i}_{ni}_{h}', 1.0)
                        for nj in range(1, N + 1):
                            constraint.add_term(f'p1_{i}_{h}_{j}_{nj}_{k}', -1.0)
                        model.add_constraint(constraint.eq(0.0))

    # GFFL0: level-crossing balance at an empty station, per phase (AMPL THM30)
    for i in range(1, M + 1):
        for u in range(1, K[i - 1] + 1):
            constraint = model.constraint_builder()
            for j in range(1, M + 1):
                if j == i:
                    continue
                for k in range(1, K[j - 1] + 1):
                    for h in range(1, K[j - 1] + 1):
                        constraint.add_term(
                            f'p1_{j}_{k}_{i}_0_{u}',
                            params.q(j - 1, i - 1, k - 1, h - 1))
                for k in range(1, K[i - 1] + 1):
                    constraint.add_term(
                        f'p1_{i}_{k}_{i}_1_{k}',
                        -params.q(i - 1, j - 1, k - 1, u - 1))
            model.add_constraint(constraint.eq(0.0))

    # see _kb/03-api-layer.md for rationale
    for i in range(1, M + 1):
        for ni in range(N):
            constraint = model.constraint_builder()
            for j in range(1, M + 1):
                if j == i:
                    continue
                for k in range(1, K[j - 1] + 1):
                    for h in range(1, K[j - 1] + 1):
                        w = params.q(j - 1, i - 1, k - 1, h - 1)
                        for u in range(1, K[i - 1] + 1):
                            constraint.add_term(f'p1_{j}_{k}_{i}_{ni}_{u}', w)
                for k in range(1, K[i - 1] + 1):
                    for h in range(1, K[i - 1] + 1):
                        constraint.add_term(
                            f'p1_{i}_{k}_{i}_{ni + 1}_{k}',
                            -params.q(i - 1, j - 1, k - 1, h - 1))
            model.add_constraint(constraint.eq(0.0))

    # see _kb/03-api-layer.md for rationale
    for i in range(1, M + 1):
        constraint = model.constraint_builder()
        for j in range(1, M + 1):
            if j == i:
                continue
            for k in range(1, K[i - 1] + 1):
                for h in range(1, K[i - 1] + 1):
                    constraint.add_term(
                        f'U_{i}_{k}', params.q(i - 1, j - 1, k - 1, h - 1))
            for k in range(1, K[j - 1] + 1):
                for h in range(1, K[j - 1] + 1):
                    w = params.q(j - 1, i - 1, k - 1, h - 1)
                    for u in range(1, K[i - 1] + 1):
                        for nj in range(1, N + 1):
                            constraint.add_term(f'p1_{i}_{u}_{j}_{nj}_{k}', -w)
                    for u in range(1, K[i - 1] + 1):
                        constraint.add_term(f'p1_{j}_{k}_{i}_0_{u}', -w)
        model.add_constraint(constraint.eq(0.0))


def _add_bound_constraints_lr(model: MapqnLpModel, params: LinearReductionParameters) -> None:
    """Add bound constraints for linear reduction model."""
    M = params.M
    N = params.N
    K = params.K

    # UUB1: sum{k} U[i,k] <= 1
    for i in range(1, M + 1):
        constraint = model.constraint_builder()
        for k in range(1, K[i - 1] + 1):
            constraint.add_term(f'U_{i}_{k}', 1.0)
        model.add_constraint(constraint.leq(1.0))

    # QUB1: Q[j,k] <= N * U[j,k]
    for j in range(1, M + 1):
        for k in range(1, K[j - 1] + 1):
            constraint = model.constraint_builder()
            constraint.add_term(f'Q_{j}_{k}', 1.0)
            constraint.add_term(f'U_{j}_{k}', -float(N))
            model.add_constraint(constraint.leq(0.0))


__all__ = ['mapqn_bnd_lr']
