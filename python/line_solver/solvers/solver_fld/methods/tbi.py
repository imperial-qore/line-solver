"""
Trajectory-based iteration (TBI) method for the FLD solver.

Decomposes the transient fluid ODE by partitioning the station set into cells
and applying Jacobi waveform relaxation: on each growing-horizon time segment,
the IVP of every cell is solved with the state of the other cells frozen at the
trajectory computed in the previous sweep. Cross-cell inflows are therefore
evaluated on frozen trajectories, interior flows on the live cell state,
matching the decomposed ODEs of the TBI method. Sweeps repeat until the
trajectory sup-norm gap falls below tbi_tol.

Reference: Sheldon, Tuncer, Casale, "TBI: Transient Hierarchical Modeling of
Large-Scale Vehicle Sharing Systems", IEEE T-ITS.

Port from MATLAB solver_fluid_tbi_iteration.m and tbi_partition.m.
"""

import numpy as np
import time
from typing import List

from scipy.integrate import solve_ivp

from ..options import SolverFLDOptions, FLDResult
from ..utils.ratemult import fluid_interpcols
from .closing import ClosingMethod
from line_solver.api.sn import SchedStrategy, NodeType


def tbi_partition(sn, options: SolverFLDOptions) -> List[List[int]]:
    """Partition the station set into cells for trajectory-based iteration.

    Honors options.config['tbi_cells'], a list of disjoint station-index lists
    covering range(nstations). Otherwise stations are agglomerated greedily on
    the symmetrized station-level routing weights, targeting
    options.config['tbi_cellsize'] stations per cell (default 5), so that
    strongly coupled stations share a cell and connecting flows stay weak.

    Port of MATLAB tbi_partition.m (station indices are 0-based here).
    """
    M = sn.nstations
    K = sn.nclasses

    config = getattr(options, 'config', None) or {}

    tbi_cells = config.get('tbi_cells', None)
    if tbi_cells is not None and len(tbi_cells) > 0:
        cells = [list(np.asarray(c, dtype=int).flatten()) for c in tbi_cells]
        covered = sorted(int(i) for c in cells for i in c)
        if covered != list(range(M)):
            raise ValueError(
                "options.config['tbi_cells'] must be a partition of the "
                "station set range(%d)." % M)
        return cells

    cellsize = config.get('tbi_cellsize', None)
    if cellsize is None:
        cellsize = 5

    # station-level coupling weights, aggregated over classes and symmetrized
    rt = np.asarray(sn.rt)
    W = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            W[i, j] = np.sum(rt[i * K:(i + 1) * K, j * K:(j + 1) * K])
    A = W + W.T
    np.fill_diagonal(A, 0.0)

    ncells_target = max(1, int(np.ceil(M / cellsize)))
    cells = [[i] for i in range(M)]
    C = A.copy()  # inter-cell coupling weights
    while len(cells) > ncells_target:
        n = len(cells)
        szs = np.array([len(c) for c in cells])
        # most coupled cell pair, respecting the cell size cap
        maxc = -1.0
        besta = -1
        bestb = -1
        for a in range(n):
            for b in range(a + 1, n):
                if szs[a] + szs[b] <= 2 * cellsize and C[a, b] > maxc:
                    maxc = C[a, b]
                    besta = a
                    bestb = b
        if besta < 0:  # every merge exceeds the size cap: merge the two smallest
            order = np.argsort(szs, kind='stable')
            besta = min(order[0], order[1])
            bestb = max(order[0], order[1])

        cells[besta] = cells[besta] + cells[bestb]
        C[besta, :] = C[besta, :] + C[bestb, :]
        C[:, besta] = C[:, besta] + C[:, bestb]
        C[besta, besta] = 0.0
        del cells[bestb]
        C = np.delete(C, bestb, axis=0)
        C = np.delete(C, bestb, axis=1)

    return cells


class TBIMethod(ClosingMethod):
    """Trajectory-based iteration method for closed fluid networks.

    Reuses the closing-method ODE machinery (index construction, jump matrix,
    state-dependent rates, and metric post-processing) and overrides only the
    ODE integration with the TBI segment/sweep waveform-relaxation loop.
    """

    def solve(self) -> FLDResult:
        start_time = time.time()

        M = self.sn.nstations
        K = self.sn.nclasses

        # Gates: closed models only, and no cache nodes (mirrors other native
        # FLD gates in this package).
        self._assert_closed()
        self._assert_no_cache()

        # Extract network parameters (reused from ClosingMethod)
        Mu, Phi, phases = self._extract_service_params()
        rt = self._get_routing_matrix()
        nservers = self._get_nservers()
        sched = self._get_sched()
        schedparam = self._get_schedparam()

        # Initial state: honor an explicitly supplied init_sol (as MATLAB
        # solver_fluid.m does), else distribute jobs over reachable points.
        x0 = self._resolve_initial_state(M, K, phases)

        xvec_it, xvec_t, t = self._solve_fluid_ode_tbi(
            M, K, Mu, Phi, phases, rt, nservers, sched, schedparam, x0, start_time
        )

        QN, UN, RN, TN, _, _, _ = self._compute_metrics_closing(
            xvec_it, xvec_t, t, M, K, Mu, Phi, phases, nservers, sched, schedparam
        )

        CN = np.sum(RN, axis=0, keepdims=True)
        XN = self._compute_system_throughput(TN, M, K)

        from line_solver.api.sn.getters import sn_get_arvr_from_tput
        from line_solver.api.sn.transforms import sn_get_residt_from_respt
        AN = sn_get_arvr_from_tput(self.sn, TN) if TN is not None else None
        WN = sn_get_residt_from_respt(self.sn, RN, None) if RN is not None else None

        # Transient queue-length trajectories from the accumulated state history.
        QNt = self._build_qnt(xvec_t, M, K, Mu, phases, sched, schedparam)

        return FLDResult(
            QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN, AN=AN, WN=WN,
            t=t if t is not None else np.array([0.0]),
            QNt=QNt, UNt={}, TNt={},
            xvec=xvec_it[-1] if xvec_it else x0,
            iterations=self.iterations,
            runtime=time.time() - start_time,
            method='tbi'
        )

    def _assert_closed(self):
        if not hasattr(self.sn, 'njobs') or self.sn.njobs is None:
            raise ValueError("TBI method requires a closed network (njobs unset)")
        njobs = np.asarray(self.sn.njobs).flatten()
        if len(njobs) == 0 or not np.all(np.isfinite(njobs)):
            raise ValueError("TBI method requires a closed network (open classes present)")

    def _assert_no_cache(self):
        nodetype = getattr(self.sn, 'nodetype', None)
        if nodetype is not None:
            for nt in nodetype:
                if nt == NodeType.CACHE:
                    raise ValueError("TBI method does not support Cache nodes")

    def _resolve_initial_state(self, M: int, K: int, phases: np.ndarray) -> np.ndarray:
        """Initial ODE state, using options.init_sol when supplied."""
        init_sol = getattr(self.options, 'init_sol', None)
        if init_sol is not None:
            x0 = np.asarray(init_sol, dtype=float).flatten()
            if len(x0) == int(np.sum(phases)):
                return x0
        return self._compute_initial_state(M, K, phases)

    def _solve_fluid_ode_tbi(
        self, M, K, Mu, Phi, phases, rt, nservers, sched, schedparam, x0, start_time
    ):
        """TBI segment/sweep waveform-relaxation loop.

        Port of MATLAB solver_fluid_tbi_iteration.m.
        """
        # ODE machinery on the full model (reused from ClosingMethod)
        q_indices, Kic, enabled, w = self._build_ode_indices(
            M, K, Mu, phases, sched, schedparam)
        all_jumps, rateBase, eventIdx = self._build_ode_system(
            M, K, Mu, Phi, phases, rt, enabled, q_indices, Kic)
        ndim = int(np.sum(Kic))

        def rates_h(x):
            return self._ode_rates_closing(
                x, M, K, enabled, q_indices, Kic, nservers, w, sched, rateBase, eventIdx)

        # see _kb/06-solver-catalog.md (Fluid: method='tbi') for the cell
        # partitioning / restricted-jump-matrix rationale
        cells = tbi_partition(self.sn, self.options)
        ncells = len(cells)
        cellmask = []
        Jcell = []
        for kc in range(ncells):
            mask = np.zeros(ndim, dtype=bool)
            for i in cells[kc]:
                for c in range(K):
                    if Kic[i, c] > 0:
                        mask[q_indices[i, c]:q_indices[i, c] + Kic[i, c]] = True
            idx = np.where(mask)[0]
            cellmask.append(idx)
            Jcell.append(all_jumps[idx, :])

        # TBI configuration
        config = getattr(self.options, 'config', None) or {}
        tbi_tol = config.get('tbi_tol', None)
        if tbi_tol is None:
            tbi_tol = 1e-3
        tbi_iter_max = config.get('tbi_iter_max', None)
        if tbi_iter_max is None:
            tbi_iter_max = 50

        tol = self.options.tol
        stiff = self.options.stiff
        method = 'LSODA' if stiff else 'RK45'
        t_start, t_end = self.options.timespan
        iter_max = self.options.iter_max
        max_time = getattr(self.options, 'timeout', float('inf'))

        # horizon growth heuristic, mirrors solver_fluid_iteration: min service
        # completion (exit) rate over enabled phases
        slowrate = []
        for i in range(M):
            for c in range(K):
                if Mu[i][c] is not None and len(Mu[i][c]) > 0:
                    slowrate.append(np.min(Mu[i][c]))
        slowrate = np.asarray(slowrate, dtype=float)
        nonZeroRates = slowrate[(slowrate > tol) & np.isfinite(slowrate)]
        if nonZeroRates.size == 0:
            nonZeroRates = np.array([1.0])  # fallback when all rates zero/infinite
        min_rate = np.min(nonZeroRates)

        xvec_it = [x0.copy()]
        t = None
        xvec_t = None
        goon = True
        it = 0
        T0 = t_start
        T = 0.0

        while (np.isfinite(t_end) and T < t_end) or (goon and it < iter_max):
            it += 1
            if time.time() - start_time > max_time:
                goon = False
                break

            y0 = xvec_it[it - 1].copy()
            if it == 1:
                T = min(t_end, abs(10.0 / min_rate))
            else:
                T = min(t_end, abs(10.0 * it / min_rate))
            trange = [T0, T]
            max_step = (T - T0) / 10.0 if T > T0 else np.inf

            # frozen trajectory on this segment, initialized constant at the
            # segment entry state (warm start of the waveform relaxation)
            tprev = np.array([T0, T])
            Yprev = np.array([y0, y0])

            delta = np.inf
            for sweep in range(tbi_iter_max):
                # see _kb/06-solver-catalog.md (Fluid: method='tbi') -- shared
                # clamped piecewise-linear interpolator, columns are time samples
                frozen_cols = np.ascontiguousarray(Yprev.T)

                cell_t = [None] * ncells
                cell_y = [None] * ncells
                tgrid = tprev.copy()
                for kc in range(ncells):
                    idx = cellmask[kc]

                    def ode_c(tt, xc, idx=idx, kc=kc):
                        xfull = fluid_interpcols(tprev, frozen_cols, tt)
                        xfull[idx] = xc
                        return Jcell[kc] @ rates_h(xfull)

                    sol = solve_ivp(
                        ode_c, trange, y0[idx], method=method,
                        rtol=tol, atol=tol * 1e-3, max_step=max_step)
                    tc = sol.t
                    yc = sol.y.T  # (ntimes, len(idx))
                    cell_t[kc] = tc
                    cell_y[kc] = yc
                    tgrid = np.union1d(tgrid, tc)

                # assemble the new full trajectory on the union time grid
                Ynew = np.zeros((len(tgrid), ndim))
                for kc in range(ncells):
                    tc = cell_t[kc]
                    tq = np.clip(tgrid, tc[0], tc[-1])
                    for jc, col in enumerate(cellmask[kc]):
                        Ynew[:, col] = np.interp(tq, tc, cell_y[kc][:, jc])

                # sup-norm gap against the previous sweep trajectory
                Yold = np.zeros((len(tgrid), ndim))
                tq = np.clip(tgrid, tprev[0], tprev[-1])
                for col in range(ndim):
                    Yold[:, col] = np.interp(tq, tprev, Yprev[:, col])
                delta = np.max(np.abs(Ynew - Yold))
                tprev = tgrid
                Yprev = Ynew
                if delta < tbi_tol or time.time() - start_time > max_time:
                    break

            if delta >= tbi_tol and self.options.verbose:
                print("TBI sweeps did not converge within tbi_iter_max=%d on "
                      "segment [%g,%g], residual gap %g." % (tbi_iter_max, T0, T, delta))

            if xvec_t is None:
                xvec_t = Yprev.copy()
                t = tprev.copy()
            else:
                xvec_t = np.vstack([xvec_t, Yprev])
                t = np.concatenate([t, tprev])
            xvec_it.append(xvec_t[-1].copy())
            self.iterations = it
            T0 = T  # for next segment

            if T >= t_end:
                goon = False

        if xvec_t is None:
            xvec_t = np.array([x0])
            t = np.array([t_start])

        return xvec_it, xvec_t, t

    def _build_qnt(self, xvec_t, M, K, Mu, phases, sched, schedparam):
        """Transient queue-length trajectories QNt[(i,r)] from state history."""
        QNt = {}
        if xvec_t is None or xvec_t.ndim != 2:
            return QNt
        q_indices, Kic, _, _ = self._build_ode_indices(
            M, K, Mu, phases, sched, schedparam)
        for i in range(M):
            for r in range(K):
                shift = q_indices[i, r]
                nph = Kic[i, r]
                if nph > 0 and shift + nph <= xvec_t.shape[1]:
                    QNt[(i, r)] = np.sum(xvec_t[:, shift:shift + nph], axis=1)
        return QNt


def solve_tbi(sn, options=None) -> FLDResult:
    """Convenience function to solve using the TBI method."""
    if options is None:
        options = SolverFLDOptions(method='tbi')
    return TBIMethod(sn, options).solve()
