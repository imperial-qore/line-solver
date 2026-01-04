"""
Linear Reduction Bounds for Product Form Networks.

Implements linear reduction bounds specialized for product-form MAP
queueing networks. Leverages product-form properties to provide
efficient bounds computation for networks with separable solutions.

This is a port of bnd_linearreduction_pf.mod AMPL model.
"""

import numpy as np
from typing import Optional

from .lpmodel import MapqnLpModel
from .parameters import PFParameters
from .solution import MapqnSolution


def mapqn_bnd_lr_pf(params: PFParameters, objective_queue: int = 1) -> MapqnSolution:
    """
    Solve the linear reduction bound for product-form networks.

    Computes lower bounds on utilization for closed queueing networks
    using a linear programming formulation that exploits product-form
    structure (no phases).

    Args:
        params: Product-form parameters including:
            - M: Number of queues
            - N: Population
            - mu: Service rates [M]
            - r: Routing matrix [M x M]
        objective_queue: Queue index to minimize utilization (1-based).
                        Default is 1.

    Returns:
        MapqnSolution containing:
            - objective_value: Minimum utilization bound
            - variables: All LP variable values

    Raises:
        ValueError: If objective_queue is out of range.
        RuntimeError: If LP is infeasible or unbounded.

    Example:
        >>> params = PFParameters(
        ...     _M=2, _N=5,
        ...     mu=np.array([1.0, 2.0]),
        ...     r=np.array([[0.0, 1.0], [1.0, 0.0]])
        ... )
        >>> sol = mapqn_bnd_lr_pf(params, objective_queue=1)
        >>> print(f"Min utilization at queue 1: {sol.objective_value:.4f}")

    Reference:
        Based on AMPL model bnd_linearreduction_pf.mod
    """
    params.validate()

    M = params.M
    N = params.N

    if not (1 <= objective_queue <= M):
        raise ValueError(f"Objective queue must be in range 1..{M}")

    model = MapqnLpModel()

    # Register all variables
    _register_variables_pf(model, M, N)

    # Add all constraints
    _add_definition_constraints_pf(model, params)
    _add_mean_indices_constraints_pf(model, params)
    _add_balance_constraints_pf(model, params)

    # Solve: minimize U[objective_queue]
    objective_var = f'U_{objective_queue}'
    solution = model.solve(objective_var, minimize=True)

    return solution


def _register_variables_pf(model: MapqnLpModel, M: int, N: int) -> None:
    """Register all variables for the product-form model."""

    # U variables: utilization at each queue
    for i in range(1, M + 1):
        model.add_variable(f'U_{i}', lb=0.0, ub=1.0)

    # Q variables: mean queue length at each queue
    for i in range(1, M + 1):
        model.add_variable(f'Q_{i}', lb=0.0, ub=float(N))

    # C variables: conditional queue lengths C[j,i]
    for j in range(1, M + 1):
        for i in range(1, M + 1):
            model.add_variable(f'C_{j}_{i}', lb=0.0, ub=float(N))

    # p1 variables: marginal probabilities p1[j,i,ni]
    for j in range(1, M + 1):
        for i in range(1, M + 1):
            for ni in range(N + 1):
                model.add_variable(f'p1_{j}_{i}_{ni}', lb=0.0, ub=1.0)

    # p1c variables: complementary marginal probabilities
    for j in range(1, M + 1):
        for i in range(1, M + 1):
            for ni in range(N + 1):
                model.add_variable(f'p1c_{j}_{i}_{ni}', lb=0.0, ub=1.0)


def _add_definition_constraints_pf(model: MapqnLpModel, params: PFParameters) -> None:
    """Add definition constraints for product-form model."""
    M = params.M
    N = params.N

    # ZER1: p1[j,j,0] = 0 for all j
    for j in range(1, M + 1):
        constraint = model.constraint_builder()
        constraint.add_term(f'p1_{j}_{j}_0', 1.0)
        model.add_constraint(constraint.eq(0.0))

    # ZER3: p1[j,i,N] = 0 for j != i
    for j in range(1, M + 1):
        for i in range(1, M + 1):
            if j != i:
                constraint = model.constraint_builder()
                constraint.add_term(f'p1_{j}_{i}_{N}', 1.0)
                model.add_constraint(constraint.eq(0.0))

    # CEQU: C[j,j] = Q[j]
    for j in range(1, M + 1):
        constraint = model.constraint_builder()
        constraint.add_term(f'C_{j}_{j}', 1.0)
        constraint.add_term(f'Q_{j}', -1.0)
        model.add_constraint(constraint.eq(0.0))

    # ONE1: sum{ni} (p1[j,i,ni] + p1c[j,i,ni]) = 1
    for j in range(1, M + 1):
        for i in range(1, M + 1):
            constraint = model.constraint_builder()
            for ni in range(N + 1):
                constraint.add_term(f'p1_{j}_{i}_{ni}', 1.0)
                constraint.add_term(f'p1c_{j}_{i}_{ni}', 1.0)
            model.add_constraint(constraint.eq(1.0))


def _add_mean_indices_constraints_pf(model: MapqnLpModel, params: PFParameters) -> None:
    """Add mean indices constraints for product-form model."""
    M = params.M
    N = params.N

    # UTIL: U[i] = sum{t,nt} p1[i,t,nt]
    for i in range(1, M + 1):
        for t in range(1, M + 1):
            constraint = model.constraint_builder()
            constraint.add_term(f'U_{i}', 1.0)
            for nt in range(N + 1):
                constraint.add_term(f'p1_{i}_{t}_{nt}', -1.0)
            model.add_constraint(constraint.eq(0.0))

    # QLEN: Q[i] = sum{ni} ni * p1[i,i,ni]
    for i in range(1, M + 1):
        constraint = model.constraint_builder()
        constraint.add_term(f'Q_{i}', 1.0)
        for ni in range(N + 1):
            constraint.add_term(f'p1_{i}_{i}_{ni}', -float(ni))
        model.add_constraint(constraint.eq(0.0))

    # CLEN: C[j,i] = sum{ni} ni * p1[j,i,ni]
    for j in range(1, M + 1):
        for i in range(1, M + 1):
            constraint = model.constraint_builder()
            constraint.add_term(f'C_{j}_{i}', 1.0)
            for ni in range(N + 1):
                constraint.add_term(f'p1_{j}_{i}_{ni}', -float(ni))
            model.add_constraint(constraint.eq(0.0))


def _add_balance_constraints_pf(model: MapqnLpModel, params: PFParameters) -> None:
    """Add balance constraints for product-form model."""
    M = params.M
    N = params.N

    # MPCB: sum{i} C[j,i] = N * U[j]
    for j in range(1, M + 1):
        constraint = model.constraint_builder()
        for i in range(1, M + 1):
            constraint.add_term(f'C_{j}_{i}', 1.0)
        constraint.add_term(f'U_{j}', -float(N))
        model.add_constraint(constraint.eq(0.0))

    # POPC: sum{i} Q[i] = N
    constraint = model.constraint_builder()
    for i in range(1, M + 1):
        constraint.add_term(f'Q_{i}', 1.0)
    model.add_constraint(constraint.eq(float(N)))

    # GFFL0: Global flow for ni = 0
    for i in range(1, M + 1):
        constraint = model.constraint_builder()
        for j in range(1, M + 1):
            if j != i:
                qji = params.q(j - 1, i - 1)
                qij = params.q(i - 1, j - 1)
                constraint.add_term(f'p1_{j}_{i}_0', qji)
                constraint.add_term(f'p1_{i}_{i}_1', -qij)
        model.add_constraint(constraint.eq(0.0))

    # GFFL: Global flow for ni in 1..N-1
    for i in range(1, M + 1):
        for ni in range(1, N):
            constraint = model.constraint_builder()
            for j in range(1, M + 1):
                if j != i:
                    qji = params.q(j - 1, i - 1)
                    qij = params.q(i - 1, j - 1)
                    constraint.add_term(f'p1_{j}_{i}_{ni}', qji)
                    constraint.add_term(f'p1_{i}_{i}_{ni + 1}', -qij)
            model.add_constraint(constraint.eq(0.0))

    # UJNT: Joint probability symmetry
    for i in range(1, M + 1):
        for j in range(1, M + 1):
            constraint = model.constraint_builder()
            for ni in range(1, N + 1):
                constraint.add_term(f'p1_{j}_{i}_{ni}', 1.0)
            for nj in range(1, N + 1):
                constraint.add_term(f'p1_{i}_{j}_{nj}', -1.0)
            model.add_constraint(constraint.eq(0.0))

    # QBAL: Queue balance
    for i in range(1, M + 1):
        constraint = model.constraint_builder()
        for j in range(1, M + 1):
            if j != i:
                qij = params.q(i - 1, j - 1)
                qji = params.q(j - 1, i - 1)

                constraint.add_term(f'U_{i}', qij)

                # sum{nj} p1[i,j,nj]
                for nj in range(1, N + 1):
                    constraint.add_term(f'p1_{i}_{j}_{nj}', -qji)

                # p1[j,i,0]
                constraint.add_term(f'p1_{j}_{i}_0', -qji)
        model.add_constraint(constraint.eq(0.0))


__all__ = ['mapqn_bnd_lr_pf']
