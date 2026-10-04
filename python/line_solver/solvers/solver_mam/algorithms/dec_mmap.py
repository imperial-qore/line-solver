"""
Decomposition with MMAP departures (dec.mmap) algorithm.

Extends dec.source with service-scaled departure processes.
The key difference is that departures from each queue are modeled as
MMAPprocesses and routed through the network, providing more accurate
inter-departure time distributions.

Algorithm:
1. Initialize departure processes from service distributions
2. Iterate: traffic() → arrivals → solve queues → build departures
3. Extract response time PH from queue solutions
4. Compress MMAP states when exceeding space_max
5. Route departures through network

References:
    MATLAB: matlab/src/solvers/MAM/solver_mam.m (187 lines)
"""

import numpy as np
import time
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass

from . import MAMAlgorithm, MAMResult
from ..utils.network_adapter import (
    extract_mam_params,
    extract_visit_counts,
    build_routing_matrix,
    check_closed_network,
)
from ....api.sn import SchedStrategy
from ....api.mam import (
    map_lambda,
    mmap_super_safe,
    mmap_compress,
)
from ....api.solvers.mam.handler import (
    solver_mam as handler_solver_mam,
    SolverMAMOptions as HandlerOptions,
)


@dataclass
class DecMMAPOptions:
    """Options for dec.mmap solver."""
    tol: float = 1e-6
    max_iter: int = 100
    space_max: int = 500  # Max MMAP state space
    verbose: bool = False


class DecMMAPAlgorithm(MAMAlgorithm):
    """Decomposition with MMAP departures."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Structural applicability of dec.mmap; there is none beyond the
        feature set.

        THE TWO RESTRICTIONS OF THIS ALGORITHM ARE DECLARED, NOT STRUCTURAL:
        it solves OPEN models only and serves EXT/FCFS/HOL/FCFSPRPRIO/PS
        stations only, and both are things the model HAS, so both are expressed
        as SolverMAM.getMethodFeatureSet deltas for 'dec.mmap' (ClosedClass,
        SelfLoopingClass and SchedStrategy_INF dropped) rather than here. This
        predicate answers about topology and class mix, of which dec.mmap asks
        nothing further, so it stays unconditionally true -- and it is not the
        whole gate: supportsModelMethod falls through to the feature set after
        asking it.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        # see _kb/06-solver-catalog.md (MAM: "dec.mmap/dec.source on Fork-Join")
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve using decomposition with ETAQA MMAP departures.

        Delegates to the handler port of `solver_mam.m`, exactly as
        DecSourceAlgorithm delegates to the port of `solver_mam_basic.m`. The
        implementation that used to live here was an acknowledged approximation
        -- it never built a departure process ("For now, approximate as
        exponential"), so no ETAQA step ran under the dec.mmap name and the
        reported utilizations were inflated by the visit ratio.

        Args:
            sn: NetworkStruct
            options: DecMMAPOptions or SolverMAMOptions

        Returns:
            MAMResult
        """
        start_time = time.time()

        handler_opts = HandlerOptions()
        handler_opts.method = 'dec.mmap'
        if options is not None:
            for src, dst in (('tol', 'tol'), ('max_iter', 'iter_max'),
                             ('iter_max', 'iter_max'), ('iter_tol', 'iter_tol'),
                             ('verbose', 'verbose'), ('space_max', 'space_max')):
                if hasattr(options, src):
                    value = getattr(options, src)
                    if value is not None:
                        setattr(handler_opts, dst, value)
            # config carries etaqa_trunc and space_max, which the handler reads
            # off the options object itself (MATLAB reads options.config.*)
            cfg = getattr(options, 'config', None)
            if isinstance(cfg, dict):
                for key in ('etaqa_trunc', 'space_max'):
                    if cfg.get(key) is not None:
                        setattr(handler_opts, key, cfg[key])

        # The handler RAISES on a closed model and on a discipline its station
        # ladder does not serve, as solver_mam.m now does. It used to return an
        # empty result here, which this method turned into a MAMResult with no
        # metrics -- and the caller then died unpacking it, which is how a
        # refusal reached the user as "MAMResult.__init__() missing 4 required
        # positional arguments".
        res = handler_solver_mam(sn, handler_opts)
        runtime = time.time() - start_time

        TN = res.T if res.T is not None else np.zeros_like(res.Q)
        return MAMResult(
            QN=res.Q,
            UN=res.U,
            RN=res.R,
            TN=TN,
            CN=res.C,
            XN=res.X,
            totiter=res.it,
            method="dec.mmap",
            runtime=runtime
        )




def _init_departures(M: int, K: int, rates: np.ndarray, scv: np.ndarray) -> List:
    """Initialize departure processes from service distributions.

    Args:
        M: Number of stations
        K: Number of classes
        rates: Service rates
        scv: Service SCVs

    Returns:
        List of departure MMAP representations
    """
    departures = []
    for m in range(M):
        # Simple approximation: exponential departures with rate = service rate
        rate_m = rates[m, 0] if rates.shape[1] > 0 else 1.0
        D0 = np.array([[-rate_m]])
        D1 = [np.array([[rate_m]])]
        departures.append((D0, D1))

    return departures


def _solve_station_mmap(station_idx: int,
                       M: int,
                       K: int,
                       arrival_rates: np.ndarray,
                       service_rates: np.ndarray,
                       service_scv: np.ndarray,
                       nservers: float,
                       options: DecMMAPOptions) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solve a station with MMAP arrivals.

    Args:
        station_idx: Station index
        M: Total stations
        K: Classes
        arrival_rates: Arrival rates per class
        service_rates: Service rates
        service_scv: Service SCVs
        nservers: Number of servers
        options: Solver options

    Returns:
        (QN, UN, RN, XN) metrics for this station
    """
    S = 1.0 / np.maximum(service_rates, 1e-10)

    QN = np.zeros(K)
    UN = np.zeros(K)
    RN = np.zeros(K)

    lambda_total = np.sum(arrival_rates)
    ca2 = 1.0  # MMAP -> assume exponential inter-arrivals

    for k in range(K):
        if arrival_rates[k] <= 1e-10:
            UN[k] = 0.0
            QN[k] = 0.0
            RN[k] = S[k]
        else:
            rho = arrival_rates[k] * S[k] / nservers
            UN[k] = rho

            if rho >= 1.0:
                QN[k] = np.inf
                RN[k] = np.inf
            else:
                cs2 = service_scv[k]
                if nservers <= 1:
                    mean_wait = (ca2 + cs2) * rho / (2.0 * (1.0 - rho)) * S[k]
                else:
                    mean_wait = rho / (1.0 - rho) * S[k]

                RN[k] = S[k] + mean_wait
                QN[k] = arrival_rates[k] * mean_wait

    XN = np.array([lambda_total])

    return QN, UN, RN, XN
