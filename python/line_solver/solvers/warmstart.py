"""
Warm-start placement computation shared by the solvers that can start from an
auxiliary solver's steady-state solution (LDES, SSA, JMT, FLD).

The placement is an integer (nstations x nclasses) job assignment decided by
the steady-state distribution of the auxiliary solver: with SolverCTMC it is
the mode of the exact stationary distribution over the aggregate state space,
with any other solver the rounded steady-state mean queue lengths, conserving
each closed-class population. Only service stations (Queue/Delay) receive an
initial population.
"""

import numpy as np

from ..api.sn.network_struct import NodeType


def warm_start_placement(init_solver, sn) -> np.ndarray:
    """Integer job placement (nstations x nclasses) decided by the
    steady-state solution of the auxiliary solver.

    Args:
        init_solver: auxiliary solver used to compute the steady-state
            distribution (e.g. SolverCTMC or SolverMVA on the same model)
        sn: struct of the model to be warm-started

    Returns:
        placement matrix (nstations x nclasses)
    """
    from .solver_ctmc.solver_ctmc import SolverCTMC

    if isinstance(init_solver, SolverCTMC):
        return _placement_from_ctmc_steady_state(init_solver, sn)
    return _placement_from_mean_qlen(init_solver, sn)


def _placement_from_ctmc_steady_state(ctmc_solver, sn) -> np.ndarray:
    """Placement from the exact CTMC stationary distribution: aggregate the
    stationary probabilities over the aggregate (per-station, per-class job
    count) state space and return the aggregate state of maximum stationary
    probability."""
    pi = np.asarray(ctmc_solver.getSteadyState()).flatten()
    ssq = np.asarray(ctmc_solver.getStateSpaceAggr())

    # Aggregate the stationary probability over identical aggregate states
    # and locate the mode of the aggregate distribution.
    u_rows, inverse = np.unique(ssq, axis=0, return_inverse=True)
    aggr_prob = np.bincount(inverse, weights=pi[:len(inverse)])
    mode_state = u_rows[int(np.argmax(aggr_prob))]

    # space_aggr is already station-major (nstations*nclasses columns); only
    # service stations (Queue/Delay) receive an initial population.
    M = int(sn.nstations)
    K = int(sn.nclasses)
    placement = np.zeros((M, K))
    station_to_node = np.asarray(sn.stationToNode).flatten()
    for i in range(M):
        ind = int(station_to_node[i])
        if sn.nodetype[ind] not in (NodeType.QUEUE, NodeType.DELAY):
            continue
        for r in range(K):
            col = i * K + r
            if col < len(mode_state):
                placement[i, r] = mode_state[col]
    return placement


def _placement_from_mean_qlen(init_solver, sn) -> np.ndarray:
    """Placement from the steady-state mean queue lengths of a generic
    network solver: floor the per-station means and distribute the residual
    closed-class jobs by largest remainder so each closed population is
    conserved."""
    qlen = np.asarray(init_solver.getAvgQLen())
    M = int(sn.nstations)
    K = int(sn.nclasses)
    if qlen.ndim == 1:
        qlen = qlen.reshape(M, K)
    njobs = np.asarray(sn.njobs).flatten()
    refstat = np.asarray(sn.refstat).flatten()
    station_to_node = np.asarray(sn.stationToNode).flatten()

    placement = np.zeros((M, K))
    for r in range(K):
        is_closed = np.isfinite(njobs[r])
        floors = np.zeros(M, dtype=int)
        fracs = -np.ones(M)
        eligible = np.zeros(M, dtype=bool)
        for i in range(M):
            ind = int(station_to_node[i])
            eligible[i] = sn.nodetype[ind] in (NodeType.QUEUE, NodeType.DELAY)
            if not eligible[i]:
                continue
            mean = max(0.0, float(qlen[i, r]))
            if is_closed:
                floors[i] = int(np.floor(mean))
                fracs[i] = mean - floors[i]
            else:
                floors[i] = int(round(mean))
            placement[i, r] = floors[i]
        if is_closed:
            # Largest-remainder apportionment of the residual jobs.
            residual = int(round(njobs[r])) - int(floors.sum())
            while residual > 0:
                bi = int(np.argmax(fracs))
                if fracs[bi] < 0:
                    # All remainders consumed: place the leftover jobs at the
                    # reference station to conserve the population.
                    ref_station = int(refstat[r])
                    if not eligible[ref_station]:
                        ref_station = int(np.argmax(eligible))
                    placement[ref_station, r] += residual
                    residual = 0
                    break
                placement[bi, r] += 1
                fracs[bi] = -1.0
                residual -= 1
    return placement
