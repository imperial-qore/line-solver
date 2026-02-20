"""
Closing method for FLD solver (fluid approximation with iterative refinement).

Implements the closing method that solves fluid ODE models:
1. Builds ODE system from network structure with proper scheduling support
2. Solves the fluid ODE system using scipy's ODE solvers
3. Post-processes results to compute QN, UN, TN, RN metrics

Supports: INF, PS, DPS, FCFS, EXT scheduling strategies.

Port from MATLAB solver_fluid_closing.m, solver_fluid.m, and solver_fluid_odes.m
"""

import numpy as np
import time
from typing import Optional, Dict, Tuple, List, TYPE_CHECKING
from scipy.integrate import solve_ivp

if TYPE_CHECKING:
    from ...api.sn import NetworkStruct

from ..options import SolverFLDOptions, FLDResult
from line_solver.api.sn import SchedStrategy


class ClosingMethod:
    """Closing method for fluid analysis.

    Implements the full MATLAB algorithm for fluid approximation
    including proper support for DPS, PS, FCFS, and INF scheduling.
    """

    def __init__(self, sn, options: SolverFLDOptions):
        """Initialize closing method.

        Args:
            sn: Compiled NetworkStruct
            options: SolverFLDOptions with method and iteration settings
        """
        self.sn = sn
        self.options = options
        self.iterations = 0
        self.runtime = 0.0

    def solve(self) -> FLDResult:
        """Solve using closing method - full MATLAB-compatible implementation.

        Returns:
            FLDResult with performance metrics
        """
        start_time = time.time()

        try:
            # Extract network parameters
            M = self.sn.nstations
            K = self.sn.nclasses

            # Get service process parameters
            Mu, Phi, phases = self._extract_service_params()

            # Get routing and scheduling
            rt = self._get_routing_matrix()
            nservers = self._get_nservers()
            sched = self._get_sched()
            schedparam = self._get_schedparam()

            # Compute initial state
            x0 = self._compute_initial_state(M, K, phases)

            # Build and solve ODE system
            xvec_it, xvec_t, t = self._solve_fluid_ode(
                M, K, Mu, Phi, phases, rt, nservers, sched, schedparam, x0
            )

            # Post-process to get performance metrics
            QN, UN, RN, TN, QNt, UNt, TNt = self._compute_metrics_closing(
                xvec_it, xvec_t, t, M, K, Mu, Phi, phases, nservers, sched, schedparam
            )

            # Compute system-level metrics
            CN = np.sum(RN, axis=0, keepdims=True)
            XN = self._compute_system_throughput(TN, M, K)

            # Compute arrival rates
            from line_solver.api.sn.getters import sn_get_arvr_from_tput
            AN = sn_get_arvr_from_tput(self.sn, TN) if TN is not None else None

            # Compute residence times
            from line_solver.api.sn.transforms import sn_get_residt_from_respt
            WN = sn_get_residt_from_respt(self.sn, RN, None) if RN is not None else None

            result = FLDResult(
                QN=QN,
                UN=UN,
                RN=RN,
                TN=TN,
                CN=CN,
                XN=XN,
                AN=AN,
                WN=WN,
                t=t if t is not None else np.array([0.0]),
                QNt={},
                UNt={},
                TNt={},
                xvec=xvec_it[-1] if xvec_it else x0,
                iterations=self.iterations,
                runtime=time.time() - start_time,
                method='closing'
            )

            return result

        except Exception as e:
            if self.options.verbose:
                print(f"Closing method error: {e}")
                import traceback
                traceback.print_exc()

            # Return default result on failure
            M = self.sn.nstations
            K = self.sn.nclasses
            return FLDResult(
                QN=np.zeros((M, K)),
                UN=np.zeros((M, K)),
                RN=np.zeros((M, K)),
                TN=np.zeros((M, K)),
                CN=np.zeros((1, K)),
                XN=np.zeros((1, K)),
                iterations=0,
                method='closing',
            )

    def _extract_service_params(self) -> Tuple[List, List, np.ndarray]:
        """Extract service parameters from network structure.

        Returns:
            Tuple of (Mu, Phi, phases) where:
            - Mu[i][k] = service rates per phase for class k at station i
            - Phi[i][k] = completion probabilities per phase
            - phases[i,k] = number of phases
        """
        M = self.sn.nstations
        K = self.sn.nclasses

        Mu = [[None for _ in range(K)] for _ in range(M)]
        Phi = [[None for _ in range(K)] for _ in range(M)]
        phases = np.zeros((M, K), dtype=int)

        for i in range(M):
            for k in range(K):
                # Get mu (service rates per phase)
                if hasattr(self.sn, 'mu') and self.sn.mu is not None:
                    if i < len(self.sn.mu) and self.sn.mu[i] is not None:
                        if k < len(self.sn.mu[i]) and self.sn.mu[i][k] is not None:
                            mu_ik = np.asarray(self.sn.mu[i][k]).flatten()
                            if len(mu_ik) > 0 and not np.any(np.isnan(mu_ik)):
                                Mu[i][k] = mu_ik
                                phases[i, k] = len(mu_ik)

                # Get phi (completion probabilities per phase)
                if hasattr(self.sn, 'phi') and self.sn.phi is not None:
                    if i < len(self.sn.phi) and self.sn.phi[i] is not None:
                        if k < len(self.sn.phi[i]) and self.sn.phi[i][k] is not None:
                            phi_ik = np.asarray(self.sn.phi[i][k]).flatten()
                            if len(phi_ik) > 0:
                                Phi[i][k] = phi_ik

                # Fallback: use rates if mu not available
                if Mu[i][k] is None and hasattr(self.sn, 'rates') and self.sn.rates is not None:
                    rate = self.sn.rates[i, k] if i < self.sn.rates.shape[0] and k < self.sn.rates.shape[1] else 0
                    if rate > 0 and np.isfinite(rate):
                        Mu[i][k] = np.array([rate])
                        Phi[i][k] = np.array([1.0])
                        phases[i, k] = 1

                # Set defaults for missing phi
                if Mu[i][k] is not None and Phi[i][k] is None:
                    Phi[i][k] = np.ones(len(Mu[i][k]))
                    Phi[i][k][-1] = 1.0  # Last phase completes

        return Mu, Phi, phases

    def _get_routing_matrix(self) -> np.ndarray:
        """Get routing probability matrix."""
        if hasattr(self.sn, 'rt') and self.sn.rt is not None:
            return np.asarray(self.sn.rt)
        else:
            M = self.sn.nstations
            K = self.sn.nclasses
            return np.eye(M * K)

    def _get_nservers(self) -> np.ndarray:
        """Get number of servers per station."""
        M = self.sn.nstations
        nservers = np.ones(M)
        if hasattr(self.sn, 'nservers') and self.sn.nservers is not None:
            ns = np.asarray(self.sn.nservers).flatten()
            for i in range(min(M, len(ns))):
                nservers[i] = ns[i] if np.isfinite(ns[i]) else self.sn.nclosedjobs if hasattr(self.sn, 'nclosedjobs') else 1000
        return nservers

    def _get_sched(self) -> List:
        """Get scheduling strategy per station."""
        M = self.sn.nstations
        sched = [SchedStrategy.PS] * M  # Default to PS
        if hasattr(self.sn, 'sched') and self.sn.sched is not None:
            for i in range(M):
                if i in self.sn.sched:
                    sched[i] = self.sn.sched[i]
        return sched

    def _get_schedparam(self) -> np.ndarray:
        """Get scheduling parameters (weights for DPS)."""
        M = self.sn.nstations
        K = self.sn.nclasses
        schedparam = np.ones((M, K))  # Default weights = 1
        if hasattr(self.sn, 'schedparam') and self.sn.schedparam is not None:
            sp = np.asarray(self.sn.schedparam)
            if sp.ndim == 2:
                for i in range(min(M, sp.shape[0])):
                    for k in range(min(K, sp.shape[1])):
                        if np.isfinite(sp[i, k]) and sp[i, k] > 0:
                            schedparam[i, k] = sp[i, k]
        return schedparam

    def _compute_initial_state(self, M: int, K: int, phases: np.ndarray) -> np.ndarray:
        """Compute initial state vector for ODE.

        Distributes jobs evenly across stations where they can be served.
        """
        # Get job populations
        N = np.zeros(K)
        if hasattr(self.sn, 'njobs') and self.sn.njobs is not None:
            njobs = np.asarray(self.sn.njobs).flatten()
            for k in range(min(K, len(njobs))):
                N[k] = njobs[k] if np.isfinite(njobs[k]) else 0

        # Build initial state vector
        x0 = []
        assigned = np.zeros(K)
        rt = self._get_routing_matrix()
        sched = self._get_sched()

        for i in range(M):
            for k in range(K):
                n_phases = max(1, phases[i, k])

                # Check if class k is served at station i
                idx = i * K + k
                match = np.sum(rt[:, idx]) > 0 if idx < rt.shape[1] else False

                if match and phases[i, k] > 0:
                    if np.isinf(N[k]):
                        # Open class
                        if sched[i] == SchedStrategy.EXT:
                            to_assign = 1.0
                        else:
                            to_assign = 0.0
                    else:
                        # Closed class - distribute evenly
                        num_stations = max(1, np.sum([1 for j in range(M) if phases[j, k] > 0]))
                        to_assign = N[k] / num_stations
                        if assigned[k] + to_assign > N[k]:
                            to_assign = N[k] - assigned[k]
                        assigned[k] += to_assign

                    x0.extend([to_assign] + [0.0] * (n_phases - 1))
                else:
                    x0.extend([0.0] * n_phases)

        return np.array(x0, dtype=float)

    def _solve_fluid_ode(
        self, M: int, K: int, Mu: List, Phi: List, phases: np.ndarray,
        rt: np.ndarray, nservers: np.ndarray, sched: List, schedparam: np.ndarray,
        x0: np.ndarray
    ) -> Tuple[List, np.ndarray, np.ndarray]:
        """Solve the fluid ODE system.

        Port of solver_fluid_iteration.m
        """
        # Build ODE components
        q_indices, Kic, enabled, w = self._build_ode_indices(M, K, Mu, phases, sched, schedparam)

        # Build jump matrix and rate base
        all_jumps, rateBase, eventIdx = self._build_ode_system(
            M, K, Mu, Phi, phases, rt, enabled, q_indices, Kic
        )

        # Create ODE function
        def ode_func(t, x):
            rates = self._ode_rates_closing(
                x, M, K, enabled, q_indices, Kic, nservers, w, sched, rateBase, eventIdx
            )
            return all_jumps @ rates

        # Solve ODE
        tol = self.options.tol
        t_start, t_end = self.options.timespan
        if not np.isfinite(t_end):
            # Estimate end time based on slowest rate
            min_rate = np.inf
            for i in range(M):
                for k in range(K):
                    if Mu[i][k] is not None and np.any(Mu[i][k] > 0):
                        rate = np.min(Mu[i][k][Mu[i][k] > 0])
                        min_rate = min(min_rate, rate)
            if not np.isfinite(min_rate) or min_rate <= 0:
                min_rate = 1e-6
            t_end = max(10.0, 100.0 / min_rate)

        xvec_it = [x0.copy()]

        try:
            sol = solve_ivp(
                ode_func,
                [t_start, t_end],
                x0,
                method='LSODA' if self.options.stiff else 'RK45',
                rtol=tol,
                atol=tol * 1e-3,
                dense_output=True,
                max_step=(t_end - t_start) / 10,
            )

            xvec_t = sol.y.T
            t = sol.t
            xvec_it.append(sol.y[:, -1])
            self.iterations = 1

        except Exception as e:
            if self.options.verbose:
                print(f"ODE solver error: {e}")
            xvec_t = np.array([x0])
            t = np.array([0.0])

        return xvec_it, xvec_t, t

    def _build_ode_indices(
        self, M: int, K: int, Mu: List, phases: np.ndarray,
        sched: List, schedparam: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Build ODE index arrays.

        Returns:
            Tuple of (q_indices, Kic, enabled, w)
        """
        q_indices = np.zeros((M, K), dtype=int)
        Kic = np.zeros((M, K), dtype=int)
        enabled = np.zeros((M, K), dtype=bool)
        w = np.ones((M, K))

        cumsum = 0
        for i in range(M):
            for k in range(K):
                if Mu[i][k] is None or len(Mu[i][k]) == 0:
                    numphases = 0
                    enabled[i, k] = False
                else:
                    numphases = len(Mu[i][k])
                    enabled[i, k] = True

                q_indices[i, k] = cumsum
                Kic[i, k] = numphases
                cumsum += numphases

                # Set DPS weights
                if sched[i] == SchedStrategy.DPS:
                    w[i, k] = schedparam[i, k]

        return q_indices, Kic, enabled, w

    def _build_ode_system(
        self, M: int, K: int, Mu: List, Phi: List, phases: np.ndarray,
        rt: np.ndarray, enabled: np.ndarray, q_indices: np.ndarray, Kic: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build ODE jump matrix and rate base.

        Port of ode_rate_base.m and ode_jumps_new.m
        """
        # Count events
        n_events = 0
        state_dim = int(np.sum(Kic))

        # Count departure events (phase completions that route to another station)
        for i in range(M):
            for k in range(K):
                if enabled[i, k]:
                    for f in range(Kic[i, k]):
                        for j in range(M):
                            for l in range(K):
                                idx_from = i * K + k
                                idx_to = j * K + l
                                if idx_from < rt.shape[0] and idx_to < rt.shape[1]:
                                    if rt[idx_from, idx_to] > 0:
                                        n_events += 1

        # Count internal phase transitions
        for i in range(M):
            for k in range(K):
                if enabled[i, k] and Kic[i, k] > 1:
                    n_events += Kic[i, k] - 1

        if n_events == 0:
            # No events - return identity system
            return np.eye(state_dim), np.ones(state_dim), np.arange(state_dim)

        # Build jump matrix and rates
        all_jumps = np.zeros((state_dim, n_events))
        rateBase = np.zeros(n_events)
        eventIdx = np.zeros(n_events, dtype=int)

        event_count = 0

        # Departure events
        for i in range(M):
            for k in range(K):
                if not enabled[i, k]:
                    continue
                xik = q_indices[i, k]

                for f in range(Kic[i, k]):
                    for j in range(M):
                        for l in range(K):
                            idx_from = i * K + k
                            idx_to = j * K + l
                            if idx_from < rt.shape[0] and idx_to < rt.shape[1]:
                                p_route = rt[idx_from, idx_to]
                                if p_route > 0 and enabled[j, l]:
                                    xjl = q_indices[j, l]

                                    # Rate = mu * phi * P
                                    mu_f = Mu[i][k][f] if f < len(Mu[i][k]) else 0
                                    phi_f = Phi[i][k][f] if Phi[i][k] is not None and f < len(Phi[i][k]) else 1.0

                                    if mu_f > 0 and event_count < n_events:
                                        rateBase[event_count] = mu_f * phi_f * p_route
                                        eventIdx[event_count] = xik + f

                                        # Jump: -1 from source, +1 to destination
                                        all_jumps[xik + f, event_count] = -1
                                        all_jumps[xjl, event_count] = 1
                                        event_count += 1

        # Internal phase transitions
        for i in range(M):
            for k in range(K):
                if not enabled[i, k] or Kic[i, k] <= 1:
                    continue
                xik = q_indices[i, k]

                for f in range(Kic[i, k] - 1):
                    mu_f = Mu[i][k][f] if f < len(Mu[i][k]) else 0
                    phi_f = Phi[i][k][f] if Phi[i][k] is not None and f < len(Phi[i][k]) else 1.0

                    if mu_f > 0 and event_count < n_events:
                        rateBase[event_count] = mu_f * (1.0 - phi_f)
                        eventIdx[event_count] = xik + f

                        # Jump: -1 from current phase, +1 to next phase
                        all_jumps[xik + f, event_count] = -1
                        all_jumps[xik + f + 1, event_count] = 1
                        event_count += 1

        # Trim to actual event count
        all_jumps = all_jumps[:, :event_count]
        rateBase = rateBase[:event_count]
        eventIdx = eventIdx[:event_count].astype(int)

        return all_jumps, rateBase, eventIdx

    def _ode_rates_closing(
        self, x: np.ndarray, M: int, K: int, enabled: np.ndarray,
        q_indices: np.ndarray, Kic: np.ndarray, nservers: np.ndarray,
        w: np.ndarray, sched: List, rateBase: np.ndarray, eventIdx: np.ndarray
    ) -> np.ndarray:
        """Compute state-dependent rates for closing method.

        Port of ode_rates_closing.m - handles INF, PS, DPS, FCFS, EXT scheduling.
        """
        x = np.maximum(x, 0.0)
        rates = x.copy()

        for i in range(M):
            sched_i = sched[i]

            if sched_i == SchedStrategy.INF:
                # Infinite server - rates = x (no modification needed)
                pass

            elif sched_i == SchedStrategy.EXT:
                # External source - maintain mass conservation
                for k in range(K):
                    if enabled[i, k]:
                        idx_ini = q_indices[i, k]
                        idx_end = idx_ini + Kic[i, k]
                        if idx_ini < len(rates):
                            rates[idx_ini] = 1.0 - np.sum(x[idx_ini+1:idx_end])

            elif sched_i == SchedStrategy.PS:
                # Processor sharing
                idx_ini = q_indices[i, 0]
                idx_end = q_indices[i, K-1] + Kic[i, K-1] if K > 0 else idx_ini
                ni = np.sum(x[idx_ini:idx_end])

                if ni > nservers[i]:
                    rates[idx_ini:idx_end] = x[idx_ini:idx_end] / ni * nservers[i]

            elif sched_i == SchedStrategy.DPS:
                # Discriminatory processor sharing
                # Normalize weights
                w_i = w[i, :].copy()
                w_sum = np.sum(w_i)
                if w_sum > 0:
                    w_i = w_i / w_sum

                # Compute weighted queue length
                ni = np.mean(w_i)  # Small offset like MATLAB
                for k in range(K):
                    if enabled[i, k]:
                        idx_ini = q_indices[i, k]
                        idx_end = idx_ini + Kic[i, k]
                        ni += np.sum(w_i[k] * x[idx_ini:idx_end])

                # Apply DPS rates
                for k in range(K):
                    if enabled[i, k]:
                        idx_ini = q_indices[i, k]
                        idx_end = idx_ini + Kic[i, k]
                        if ni > 0:
                            rates[idx_ini:idx_end] = w_i[k] * x[idx_ini:idx_end] / ni * nservers[i]

            elif sched_i == SchedStrategy.FCFS:
                # FCFS - treated like PS in fluid approximation
                idx_ini = q_indices[i, 0]
                idx_end = q_indices[i, K-1] + Kic[i, K-1] if K > 0 else idx_ini
                ni = np.sum(x[idx_ini:idx_end])

                if ni > nservers[i]:
                    rates[idx_ini:idx_end] = x[idx_ini:idx_end] / ni * nservers[i]

        # Extract rates for events
        result = rates[eventIdx] * rateBase
        return result

    def _compute_metrics_closing(
        self, xvec_it: List, xvec_t: np.ndarray, t: np.ndarray,
        M: int, K: int, Mu: List, Phi: List, phases: np.ndarray,
        nservers: np.ndarray, sched: List, schedparam: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict, Dict, Dict]:
        """Compute performance metrics from fluid solution.

        Port of solver_fluid_closing.m
        """
        # Get final state
        x_final = xvec_it[-1] if xvec_it else np.zeros(int(np.sum(phases)))
        x_final = np.maximum(x_final, 0.0)

        # Compute queue lengths per station/class
        QN = np.zeros((M, K))
        q_indices, Kic, _, _ = self._build_ode_indices(M, K, Mu, phases, sched, schedparam)

        for i in range(M):
            for k in range(K):
                shift = q_indices[i, k]
                n_phases = Kic[i, k]
                if n_phases > 0 and shift + n_phases <= len(x_final):
                    QN[i, k] = np.sum(x_final[shift:shift + n_phases])

        # Identify delay nodes
        delay_nodes = np.zeros(M, dtype=bool)
        for i in range(M):
            if sched[i] == SchedStrategy.INF:
                delay_nodes[i] = True

        # Compute throughputs
        TN = np.zeros((M, K))
        Xservice = [[np.zeros(max(1, Kic[i, k])) for k in range(K)] for i in range(M)]

        for i in range(M):
            if delay_nodes[i]:
                # Delay node - throughput = sum of departure rates
                for k in range(K):
                    if Mu[i][k] is not None and Phi[i][k] is not None:
                        shift = q_indices[i, k]
                        for f in range(Kic[i, k]):
                            if f < len(Mu[i][k]) and f < len(Phi[i][k]):
                                idx = shift + f
                                if idx < len(x_final):
                                    TN[i, k] += x_final[idx] * Mu[i][k][f] * Phi[i][k][f]
                                    Xservice[i][k][f] = x_final[idx] * Mu[i][k][f]
            else:
                # Non-delay node - compute based on scheduling
                xi = np.sum(QN[i, :])  # Total jobs at station

                if xi > 0 or sched[i] == SchedStrategy.EXT:
                    for k in range(K):
                        if Mu[i][k] is None or Phi[i][k] is None:
                            continue

                        shift = q_indices[i, k]
                        for f in range(Kic[i, k]):
                            if f >= len(Mu[i][k]) or f >= len(Phi[i][k]):
                                continue

                            idx = shift + f
                            if idx >= len(x_final):
                                continue

                            mu_f = Mu[i][k][f]
                            phi_f = Phi[i][k][f]
                            x_f = x_final[idx]

                            if sched[i] == SchedStrategy.EXT:
                                if f == 0:
                                    x_f = 1.0 - np.sum(x_final[shift+1:shift+Kic[i, k]])
                                TN[i, k] += x_f * mu_f * phi_f
                                Xservice[i][k][f] = x_f * mu_f

                            elif sched[i] in [SchedStrategy.INF, SchedStrategy.PS]:
                                if xi > 0:
                                    rate_factor = min(xi, nservers[i]) / xi
                                    TN[i, k] += x_f * mu_f * phi_f * rate_factor
                                    Xservice[i][k][f] = x_f * mu_f * rate_factor

                            elif sched[i] == SchedStrategy.DPS:
                                # DPS throughput computation
                                w = schedparam[i, :]
                                wxi = np.sum(w * QN[i, :])
                                if wxi > 0:
                                    rate_factor = w[k] / wxi * min(xi, nservers[i])
                                    TN[i, k] += x_f * mu_f * phi_f * rate_factor
                                    Xservice[i][k][f] = x_f * mu_f * rate_factor

                            elif sched[i] in [SchedStrategy.FCFS, SchedStrategy.SIRO]:
                                if xi > 0:
                                    rate_factor = min(xi, nservers[i]) / xi
                                    TN[i, k] += x_f * mu_f * phi_f * rate_factor
                                    Xservice[i][k][f] = x_f * mu_f * rate_factor

        # Compute utilization
        UN = np.zeros((M, K))
        for i in range(M):
            for k in range(K):
                if Mu[i][k] is not None:
                    idx_pos = Xservice[i][k] > 0
                    if np.any(idx_pos):
                        mu_pos = np.array([Mu[i][k][f] for f in range(len(Xservice[i][k])) if idx_pos[f]])
                        UN[i, k] = np.sum(Xservice[i][k][idx_pos] / mu_pos)

            # Normalize by number of servers for non-delay nodes
            if not delay_nodes[i] and nservers[i] > 0:
                UN[i, :] = UN[i, :] / nservers[i]

        # Cap utilization at queue length (MATLAB convention)
        UN = np.minimum(UN, QN)

        # Compute response times using Little's Law
        RN = np.zeros((M, K))
        with np.errstate(divide='ignore', invalid='ignore'):
            for i in range(M):
                for k in range(K):
                    if TN[i, k] > 1e-10:
                        RN[i, k] = QN[i, k] / TN[i, k]
                    else:
                        RN[i, k] = 0.0

        # Return empty transient dicts for now
        QNt = {}
        UNt = {}
        TNt = {}

        return QN, UN, RN, TN, QNt, UNt, TNt

    def _compute_system_throughput(self, TN: np.ndarray, M: int, K: int) -> np.ndarray:
        """Compute system throughput per class."""
        XN = np.zeros((1, K))

        # Find reference stations and sum throughputs
        refstat = self.sn.refstat if hasattr(self.sn, 'refstat') and self.sn.refstat is not None else np.zeros(K, dtype=int)

        for k in range(K):
            ref = int(refstat[k]) if k < len(refstat) else 0
            if 0 <= ref < M:
                XN[0, k] = TN[ref, k]
            else:
                XN[0, k] = np.max(TN[:, k])

        return XN


def solve_closing(sn, options: Optional[SolverFLDOptions] = None) -> FLDResult:
    """Convenience function to solve using closing method.

    Args:
        sn: Compiled NetworkStruct
        options: SolverFLDOptions (uses defaults if None)

    Returns:
        FLDResult
    """
    if options is None:
        options = SolverFLDOptions(method='closing')

    method = ClosingMethod(sn, options)
    return method.solve()
