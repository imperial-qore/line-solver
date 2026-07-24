"""
MVA-based Linear Reduction Bounds for MAP Queueing Networks.

Implements linear reduction bounds using Mean Value Analysis (MVA)
techniques for MAP queueing networks.

This is a port of bnd_mvaversion.mod AMPL model.
"""

import numpy as np
from typing import Optional

from .lpmodel import MapqnLpModel
from .parameters import MVAVersionParameters
from .solution import MapqnSolution


def mapqn_bnd_lr_mva(
    params: MVAVersionParameters,
    objective_queue: int,
    objective_level: int
) -> MapqnSolution:
    """
    Solve the MVA-based linear reduction bound for MAP queueing networks.

    Computes bounds on utilization for closed queueing networks with
    a MAP queue using MVA-based linear programming formulation.

    Args:
        params: MVA version parameters including:
            - M: Number of queues
            - N: Population
            - K: Number of levels (scalar)
            - muM: Service rates for queues 1..M-1
            - muMAP: Service rate matrix for MAP queue [K x K]
            - r: Routing matrix [M x M]
            - v: Level change rate matrix [K x K]
        objective_queue: Queue index to maximize utilization (1-based).
        objective_level: Level index to maximize utilization (1-based).

    Returns:
        MapqnSolution containing:
            - objective_value: Maximum utilization bound
            - variables: All LP variable values

    Raises:
        ValueError: If indices are out of range.
        RuntimeError: If LP is infeasible or unbounded.

    Reference:
        Based on AMPL model bnd_mvaversion.mod
    """
    params.validate()

    M = params.M
    N = params.N
    K = params.K

    if not (1 <= objective_queue <= M):
        raise ValueError(f"Objective queue must be in range 1..{M}")
    if not (1 <= objective_level <= K):
        raise ValueError(f"Objective level must be in range 1..{K}")

    model = MapqnLpModel()

    # Register all variables
    _register_variables_mva(model, M, N, K)

    # Add all constraints
    _add_constraints_mva(model, params)

    # Solve: maximize UN[objective_queue, objective_level]
    objective_var = f'UN_{objective_queue}_{objective_level}'
    solution = model.solve(objective_var, minimize=False)

    return solution


def _register_variables_mva(model: MapqnLpModel, M: int, N: int, K: int) -> None:
    """Register all variables for MVA version model."""

    # UN variables: utilization at each queue-level
    for i in range(1, M + 1):
        for k in range(1, K + 1):
            model.add_variable(f'UN_{i}_{k}', lb=0.0, ub=1.0)

    # QN variables: queue length at each queue-level
    for i in range(1, M + 1):
        for k in range(1, K + 1):
            model.add_variable(f'QN_{i}_{k}', lb=0.0, ub=float(N))

    # B variables: auxiliary variables B[j,k,i]
    for j in range(1, M + 1):
        for k in range(1, K + 1):
            for i in range(1, M + 1):
                model.add_variable(f'B_{j}_{k}_{i}', lb=0.0, ub=float(N))


def _add_constraints_mva(model: MapqnLpModel, params: MVAVersionParameters) -> None:
    """Add all constraints for MVA version model."""
    M = params.M
    N = params.N
    K = params.K

    # QNB: QN[i,k] >= B[j,k,i]
    for i in range(1, M + 1):
        for k in range(1, K + 1):
            for j in range(1, M + 1):
                constraint = model.constraint_builder()
                constraint.add_term(f'QN_{i}_{k}', 1.0)
                constraint.add_term(f'B_{j}_{k}_{i}', -1.0)
                model.add_constraint(constraint.geq(0.0))

    # UMAX: sum{k} UN[i,k] <= 1
    for i in range(1, M + 1):
        constraint = model.constraint_builder()
        for k in range(1, K + 1):
            constraint.add_term(f'UN_{i}_{k}', 1.0)
        model.add_constraint(constraint.leq(1.0))

    # POPCONSTR: sum{i,k} QN[i,k] = N
    constraint = model.constraint_builder()
    for i in range(1, M + 1):
        for k in range(1, K + 1):
            constraint.add_term(f'QN_{i}_{k}', 1.0)
    model.add_constraint(constraint.eq(float(N)))

    # FLOW: Flow balance equations
    for i in range(1, M + 1):
        constraint = model.constraint_builder()
        for k in range(1, K + 1):
            for m in range(1, K + 1):
                for w in range(1, M + 1):
                    q_in = params.q(w - 1, i - 1, k - 1, m - 1)
                    q_out = params.q(i - 1, w - 1, m - 1, k - 1)
                    constraint.add_term(f'UN_{w}_{k}', q_in)
                    constraint.add_term(f'UN_{i}_{m}', -q_out)
        model.add_constraint(constraint.eq(0.0))

    # UBAL: Balance for MAP queue levels
    for k in range(1, K + 1):
        constraint = model.constraint_builder()
        for h in range(1, K + 1):
            if h != k:
                for w in range(1, M + 1):
                    q_out = params.q(M - 1, w - 1, k - 1, h - 1)
                    q_in = params.q(M - 1, w - 1, h - 1, k - 1)
                    constraint.add_term(f'UN_{M}_{k}', q_out)
                    constraint.add_term(f'UN_{M}_{h}', -q_in)
        model.add_constraint(constraint.eq(0.0))

    # QBAL: Queue balance for MAP queue
    for k in range(1, K + 1):
        constraint = model.constraint_builder()

        # Outgoing from level k
        for h in range(1, K + 1):
            if h != k:
                for w in range(1, M + 1):
                    q = params.q(M - 1, w - 1, k - 1, h - 1)
                    constraint.add_term(f'QN_{M}_{k}', q)

        # Incoming to level k from other queues
        for m in range(1, K + 1):
            for j in range(1, M):  # j < M (non-MAP queues)
                q = params.q(M - 1, j - 1, m - 1, k - 1)
                constraint.add_term(f'UN_{M}_{m}', q)

        # Incoming from other queues j to MAP
        for j in range(1, M):  # j < M
            q = params.q(j - 1, M - 1, k - 1, k - 1)
            constraint.add_term(f'UN_{j}_{k}', -q)

        # Incoming from other levels
        for h in range(1, K + 1):
            if h != k:
                for w in range(1, M + 1):
                    q = params.q(M - 1, w - 1, h - 1, k - 1)
                    constraint.add_term(f'QN_{M}_{h}', -q)

        model.add_constraint(constraint.eq(0.0))

    # MCC: Mean customer count constraint
    for i in range(1, M + 1):
        constraint = model.constraint_builder()

        # First sum: outgoing from queue i
        for k in range(1, K + 1):
            for m in range(1, K + 1):
                for w in range(1, M + 1):
                    if w != i:
                        q = params.q(i - 1, w - 1, k - 1, m - 1)
                        constraint.add_term(f'QN_{i}_{k}', q)

        # Second sum: incoming to queue i
        for k in range(1, K + 1):
            for m in range(1, K + 1):
                for j in range(1, M + 1):
                    if j != i:
                        q = params.q(j - 1, i - 1, k - 1, m - 1)
                        constraint.add_term(f'QN_{j}_{k}', q)

        # Third sum: B variables
        for k in range(1, K + 1):
            for m in range(1, K + 1):
                for j in range(1, M + 1):
                    if j != i:
                        q = params.q(j - 1, i - 1, k - 1, m - 1)
                        for wp in range(1, M + 1):
                            if wp != i and wp != j:
                                constraint.add_term(f'B_{j}_{k}_{wp}', q)

        # RHS: (N+1) * incoming utilization
        for k in range(1, K + 1):
            for m in range(1, K + 1):
                for j in range(1, M + 1):
                    if j != i:
                        q = params.q(j - 1, i - 1, k - 1, m - 1)
                        constraint.add_term(f'UN_{j}_{k}', -(N + 1) * q)

        model.add_constraint(constraint.eq(0.0))

    # MCC2: Second mean customer count constraint
    for i in range(1, M + 1):
        constraint = model.constraint_builder()

        # LHS: outgoing
        for k in range(1, K + 1):
            for m in range(1, K + 1):
                for w in range(1, M + 1):
                    if w != i:
                        q = params.q(i - 1, w - 1, k - 1, m - 1)
                        constraint.add_term(f'QN_{i}_{k}', q)

        # RHS: incoming
        for k in range(1, K + 1):
            for m in range(1, K + 1):
                for j in range(1, M + 1):
                    if j != i:
                        q = params.q(j - 1, i - 1, k - 1, m - 1)
                        constraint.add_term(f'B_{j}_{k}_{i}', -q)
                        constraint.add_term(f'UN_{j}_{k}', -q)

        model.add_constraint(constraint.eq(0.0))

    # QMAX: QN[w,k] <= N * UN[w,k]
    for w in range(1, M + 1):
        for k in range(1, K + 1):
            constraint = model.constraint_builder()
            constraint.add_term(f'QN_{w}_{k}', 1.0)
            constraint.add_term(f'UN_{w}_{k}', -float(N))
            model.add_constraint(constraint.leq(0.0))

    # QMIN: sum{w} QN[w,k] >= N * UN[j,k]
    for k in range(1, K + 1):
        for j in range(1, M + 1):
            constraint = model.constraint_builder()
            for w in range(1, M + 1):
                constraint.add_term(f'QN_{w}_{k}', 1.0)
            constraint.add_term(f'UN_{j}_{k}', -float(N))
            model.add_constraint(constraint.geq(0.0))


__all__ = ['mapqn_bnd_lr_mva']
