"""
Ko-Pender fluid and diffusion limits for the (MAP_t/Ph_t/inf)^N network.

Implements Y. M. Ko and J. Pender, "Diffusion limits for the
(MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45 (2017) 248-253.

The paper writes the network as a superposition of unit-rate time-changed
Poisson processes (its (3.1)-(3.2)) and derives a fluid limit (Thm 3.1) and a
diffusion limit (Thm 3.3). This module integrates the mean and the covariance of
that limit jointly:

    dq/dt     = F(t, q) = A f(t, q)
    dSigma/dt = J Sigma + Sigma J' + G,   J = A df/dq,  G = A diag(f) A'

with A the jump matrix whose column e is the jump vector l_e of event e and f
the event rate vector. G is exactly dH dH' of Thm 3.3, each independent Poisson
term contributing l_e l_e' f_e. Where f is affine in q -- infinite-server
stations and the arrival phase process -- J does not depend on q and both
equations close exactly, so for the (MAP_t/Ph_t/inf)^N case the mean and the
covariance are exact rather than asymptotic. Finite-server stations are admitted
through the usual fluid min(x, c) term, where the covariance degrades to a
linear-noise approximation and a warning is emitted.

Why this method does not reuse the closing ODE
----------------------------------------------
The closing method routes a departure from the source to the destination station
and lets mass return from the network through the STATIONARY arrival-instant
vector pie, which replaces the D1' operator by the rank-one map
pie (x) (D1 e). That is the PH renewal process with representation (pie, D0):
its stationary arrival rate is exact but its autocorrelation is gone. Since a
non-renewal arrival stream is the entire point of a MAP, this method keeps the
paper's own events, in which a phase change of the arrival MAP that generates an
arrival (rate d1_{kj}) moves the phase from k to j AND starts a job, so D1 acts
as itself.

State layout
------------
Station-major, arrival phases before service phases, matching the paper's
q_m = (u_{1m}, ..., u_{hA,m}, x_{1m}, ..., x_{hS,m}):

  u-block  one per (EXT station, class) with an arrival process; the phase
           occupancy of the arrival MAP, summing to 1 at all times
  x-block  one per (queueing station, class); the fluid count in each service
           phase

An open network has no mass returning to the source, so departures out of the
network simply leave the state; routing to an EXT station is what "leaves the
network" means in sn.rt, whose rows are closed through the Source.
"""

import time
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from ..options import FLDResult, SolverFLDOptions
from ..utils.phase_type import (extract_mu_phi_from_phase_type, is_mapt, is_pht,
                                prepare_phase_type_structures, schedule_segments)
from ..utils.ratemult import ratemult_max_step
from line_solver.api.sn import SchedStrategy
from line_solver.constants import GlobalConstants


class KoPenderMethod:
    """Fluid + diffusion limits of a (MAP_t/Ph_t/inf)^N network."""

    def __init__(self, sn, options: SolverFLDOptions):
        self.sn = sn
        self.options = options
        self.runtime = 0.0
        self.linear = True  # cleared when a finite-server station is present

    # ------------------------------------------------------------------
    # structure
    # ------------------------------------------------------------------

    def _blocks(self):
        """Index the u- and x-blocks of the state vector.

        Returns:
            (ublocks, xblocks, dim) where each block list holds
            (station, class, offset, nphases).
        """
        sn = self.sn
        M, K = sn.nstations, sn.nclasses
        proc_matrix, pie_dict, phases = self._nominal_structs()
        ublocks, xblocks = [], []
        off = 0
        for i in range(M):
            is_ext = int(sn.sched[i]) == int(SchedStrategy.EXT)
            for c in range(K):
                h = int(phases[i, c])
                if h <= 0:
                    continue
                rate = sn.rates[i, c] if sn.rates is not None else np.nan
                if not np.isfinite(rate) or rate <= 0:
                    continue  # disabled class at this station
                if is_ext:
                    ublocks.append((i, c, off, h))
                else:
                    xblocks.append((i, c, off, h))
                off += h
        return ublocks, xblocks, off

    def _nominal_structs(self):
        cached = getattr(self, '_nom_cache', None)
        if cached is None:
            cached = prepare_phase_type_structures(self.sn)
            self._nom_cache = cached
        return cached

    def _segment_at(self, t: float) -> Optional[Dict[Tuple[int, int], int]]:
        """Which schedule segment each MAPt/PHt station-class is in at time t.

        Returns None when no station carries a schedule, so the caller can reuse
        the nominal structures unchanged.
        """
        sn = self.sn
        entries = getattr(self, '_sched_entries', None)
        if entries is None:
            entries = []
            for i in range(sn.nstations):
                for c in range(sn.nclasses):
                    kind = 'MAPt' if is_mapt(sn, i, c) else ('PHt' if is_pht(sn, i, c) else None)
                    if kind is None:
                        continue
                    bp, pairs, cyclic = schedule_segments(sn.proc[i][c], kind)
                    entries.append((i, c, bp, len(pairs), cyclic))
            self._sched_entries = entries
        if not entries:
            return None
        seg = {}
        for (i, c, bp, nseg, cyclic) in entries:
            period = float(bp[-1] - bp[0])
            off = float(t) - bp[0]
            if cyclic and period > 0:
                off = off % period
            elif off < 0.0 or off >= period:
                seg[(i, c)] = None  # silent past a non-cyclic horizon
                continue
            idx = int(np.searchsorted(bp[1:], bp[0] + off, side='right'))
            seg[(i, c)] = min(idx, nseg - 1)
        return seg

    def _structs_at(self, t: float):
        """(proc, pie) in force at t, cached per segment tuple."""
        seg = self._segment_at(t)
        if seg is None:
            # _nominal_structs caches the FULL (proc, pie, ph) triple that
            # prepare_phase_type_structures returns; this accessor's contract is
            # the (proc, pie) pair, as the cached branch below returns. Returning
            # the triple here made every caller's `proc, pie = ...` raise
            # "too many values to unpack" on any model WITHOUT a MAPt/PHt
            # schedule -- which is every ordinary open network, i.e. the whole
            # non-time-varying half of the method's domain.
            nominal = self._nominal_structs()
            return nominal[0], nominal[1]
        key = tuple(sorted((k, v) for k, v in seg.items()))
        cache = getattr(self, '_seg_cache', None)
        if cache is None:
            cache = {}
            self._seg_cache = cache
        if key not in cache:
            if any(v is None for v in seg.values()):
                # A silent process: zero every rate of the affected block.
                proc, pie, ph = prepare_phase_type_structures(self.sn)
                for (i, c), v in seg.items():
                    if v is None:
                        D0, D1 = proc[i][c]
                        proc[i][c] = [np.zeros_like(D0), np.zeros_like(D1)]
                cache[key] = (proc, pie, ph)
            else:
                cache[key] = prepare_phase_type_structures(self.sn, segment=seg)
        return cache[key][0], cache[key][1]

    def _breakpoint_grid(self, t0: float, tend: float) -> np.ndarray:
        """Union of every schedule's segment boundaries inside (t0, tend)."""
        self._segment_at(t0)  # populate _sched_entries
        bounds = []
        for (_, _, bp, _, cyclic) in (self._sched_entries or []):
            period = float(bp[-1] - bp[0])
            if cyclic and period > 0:
                kmax = int(np.ceil((tend - t0) / period)) + 2
                grid = np.concatenate([bp + k * period for k in range(-1, kmax + 1)])
            else:
                grid = bp
            bounds.extend(grid[(grid > t0) & (grid < tend)].tolist())
        return np.unique(np.asarray([t0] + sorted(bounds) + [tend], dtype=float))

    # ------------------------------------------------------------------
    # events
    # ------------------------------------------------------------------

    def _build_events(self, ublocks, xblocks, dim):
        """Enumerate the paper's five event families.

        Each event carries its jump vector and a descriptor from which the rate
        at time t is assembled: (kind, i, c, k, j, dest_i, dest_c, dest_k).
        The rates themselves are time- and state-dependent and are evaluated by
        _rates.
        """
        sn = self.sn
        K = sn.nclasses
        rt = np.asarray(sn.rt)
        ublk = {(i, c): (off, h) for (i, c, off, h) in ublocks}
        xblk = {(i, c): (off, h) for (i, c, off, h) in xblocks}

        jumps = []
        descs = []

        def add(jump_pairs, desc):
            col = np.zeros(dim)
            for idx, delta in jump_pairs:
                col[idx] += delta
            jumps.append(col)
            descs.append(desc)

        # (A0) arrival-MAP phase change without an arrival
        for (i, c, off, h) in ublocks:
            for k in range(h):
                for j in range(h):
                    if k != j:
                        add([(off + k, -1.0), (off + j, 1.0)], ('A0', i, c, k, j, -1, -1, -1))

        # (A1) arrival-MAP phase change WITH an arrival, routed into a service phase
        for (i, c, off, h) in ublocks:
            for (n, l, noff, nh) in xblocks:
                p = rt[i * K + c, n * K + l]
                if p <= 0:
                    continue
                for k in range(h):
                    for j in range(h):
                        for ip in range(nh):
                            add([(off + k, -1.0), (off + j, 1.0), (noff + ip, 1.0)],
                                ('A1', i, c, k, j, n, l, ip))

        # (S) service phase change inside a station
        for (i, c, off, h) in xblocks:
            for p in range(h):
                for q in range(h):
                    if p != q:
                        add([(off + p, -1.0), (off + q, 1.0)], ('S', i, c, p, q, -1, -1, -1))

        # (D) service completion that leaves the network, and
        # (R) service completion routed to another station
        for (i, c, off, h) in xblocks:
            p_out = 0.0
            for j in range(sn.nstations):
                for l in range(K):
                    if (j, l) not in xblk:
                        p_out += rt[i * K + c, j * K + l]
            for p in range(h):
                if p_out > 0:
                    add([(off + p, -1.0)], ('D', i, c, p, -1, -1, -1, -1))
            for (n, l, noff, nh) in xblocks:
                p = rt[i * K + c, n * K + l]
                if p <= 0:
                    continue
                for pp in range(h):
                    for ip in range(nh):
                        add([(off + pp, -1.0), (noff + ip, 1.0)],
                            ('R', i, c, pp, -1, n, l, ip))

        A = np.array(jumps).T if jumps else np.zeros((dim, 0))
        return A, descs

    def _rates(self, t, q, descs, ublocks, xblocks):
        """Rate vector f(t, q) of the enumerated events."""
        sn = self.sn
        K = sn.nclasses
        rt = np.asarray(sn.rt)
        proc, pie = self._structs_at(t)
        ublk = {(i, c): (off, h) for (i, c, off, h) in ublocks}
        xblk = {(i, c): (off, h) for (i, c, off, h) in xblocks}

        # server-capacity factor per x-block, theta_p = x_p * min(n,c)/n
        theta_scale = {}
        for (i, c, off, h) in xblocks:
            sched = int(sn.sched[i])
            if sched == int(SchedStrategy.INF):
                theta_scale[(i, c)] = 1.0
                continue
            nsrv = float(sn.nservers[i]) if sn.nservers is not None else 1.0
            if not np.isfinite(nsrv):
                theta_scale[(i, c)] = 1.0
                continue
            # station total across classes sharing the server
            ni = 0.0
            for (j, l, joff, jh) in xblocks:
                if j == i:
                    ni += float(np.sum(q[joff:joff + jh]))
            theta_scale[(i, c)] = 1.0 if ni <= nsrv else nsrv / ni

        f = np.zeros(len(descs))
        for e, (kind, i, c, k, j, n, l, ip) in enumerate(descs):
            if kind == 'A0':
                off, h = ublk[(i, c)]
                D0 = proc[i][c][0]
                f[e] = D0[k, j] * max(q[off + k], 0.0)
            elif kind == 'A1':
                off, h = ublk[(i, c)]
                D1 = proc[i][c][1]
                beta = pie[n][l]
                p = rt[i * K + c, n * K + l]
                f[e] = D1[k, j] * p * beta[ip] * max(q[off + k], 0.0)
            elif kind == 'S':
                off, h = xblk[(i, c)]
                D0 = proc[i][c][0]
                f[e] = D0[k, j] * max(q[off + k], 0.0) * theta_scale[(i, c)]
            elif kind == 'D':
                off, h = xblk[(i, c)]
                D1 = proc[i][c][1]
                p_out = 0.0
                for jj in range(sn.nstations):
                    for ll in range(K):
                        if (jj, ll) not in xblk:
                            p_out += rt[i * K + c, jj * K + ll]
                f[e] = float(np.sum(D1[k, :])) * p_out * max(q[off + k], 0.0) * theta_scale[(i, c)]
            else:  # 'R'
                off, h = xblk[(i, c)]
                D1 = proc[i][c][1]
                beta = pie[n][l]
                p = rt[i * K + c, n * K + l]
                f[e] = float(np.sum(D1[k, :])) * p * beta[ip] * max(q[off + k], 0.0) \
                    * theta_scale[(i, c)]
        return f

    # ------------------------------------------------------------------
    # ODE
    # ------------------------------------------------------------------

    def _drift(self, t, q, A, descs, ublocks, xblocks):
        return A @ self._rates(t, q, descs, ublocks, xblocks)

    def _jacobian(self, t, q, A, descs, ublocks, xblocks):
        """dF/dq by central differences on the rate function.

        Differentiating the assembled rates rather than hand-coding a Jacobian
        keeps every capacity term (min(x,c) and the class-sharing denominator)
        consistent with the drift actually integrated.
        """
        n = q.size
        J = np.zeros((n, n))
        scale = max(1.0, float(np.max(np.abs(q))) if q.size else 1.0)
        hstep = 1e-6 * scale
        for m in range(n):
            qp = q.copy(); qp[m] += hstep
            qm = q.copy(); qm[m] -= hstep
            J[:, m] = (self._drift(t, qp, A, descs, ublocks, xblocks)
                       - self._drift(t, qm, A, descs, ublocks, xblocks)) / (2.0 * hstep)
        return J

    def _augmented_rhs(self, t, z, A, descs, ublocks, xblocks, n):
        q = z[:n]
        Sigma = z[n:].reshape(n, n)
        f = self._rates(t, q, descs, ublocks, xblocks)
        dq = A @ f
        J = self._jacobian(t, q, A, descs, ublocks, xblocks)
        G = A @ np.diag(f) @ A.T
        dS = J @ Sigma + Sigma @ J.T + G
        return np.concatenate([dq, dS.ravel()])

    def _initial_state(self, ublocks, xblocks, dim, t0):
        """u(0) = stationary phase vector of the arrival process at t0; x(0) = 0.

        The paper starts the network empty with the modulating chain in some
        phase; taking its stationary distribution is the choice that makes a
        single-segment MAP_t reproduce the stationary MAP exactly at every t.
        """
        proc, _ = self._structs_at(t0)
        q0 = np.zeros(dim)
        for (i, c, off, h) in ublocks:
            D0, D1 = proc[i][c]
            Q = np.asarray(D0) + np.asarray(D1)
            try:
                Amat = np.vstack([Q.T, np.ones(h)])
                b = np.zeros(h + 1); b[-1] = 1.0
                theta, *_ = np.linalg.lstsq(Amat, b, rcond=None)
                theta = np.maximum(theta, 0.0)
                total = float(np.sum(theta))
                theta = theta / total if total > 0 else np.ones(h) / h
            except np.linalg.LinAlgError:
                theta = np.ones(h) / h
            q0[off:off + h] = theta
        # NOT options.init_sol: that is laid out for the CLOSING state vector,
        # which differs from this method's (u-blocks before x-blocks, no mass
        # returning to the source). Consuming it would silently zero the source
        # phase mass and with it the whole network.
        #
        # A WRONG-SIZED SEED IS REFUSED, not ignored: dropping it would integrate
        # from the default initial condition under the caller's name and return a
        # plausible trajectory for a model the caller did not ask about.
        cfg = getattr(self.options, 'config', None)
        init = cfg.get('kp_init_sol') if isinstance(cfg, dict) else None
        if init is not None and np.size(init) > 0:
            if np.size(init) != dim:
                raise ValueError(
                    "config['kp_init_sol'] has %d entries but the 'kp' state vector of this "
                    "model has %d, laid out station-major over the (station, class) blocks. "
                    "It is NOT laid out like options.init_sol." % (int(np.size(init)), dim))
            q0 = np.asarray(init, dtype=float).ravel().copy()
        return q0

    # ------------------------------------------------------------------
    # driver
    # ------------------------------------------------------------------

    def solve(self) -> FLDResult:
        start = time.time()
        sn = self.sn
        M, K = sn.nstations, sn.nclasses

        if sn.njobs is not None and np.any(np.isfinite(np.asarray(sn.njobs, dtype=float))):
            raise ValueError(
                "the 'kp' method analyses the open (MAP_t/Ph_t/inf)^N network of "
                "Ko and Pender (2017); a closed class has no arrival process to "
                "modulate. Use 'closing' or 'matrix' for closed models.")

        ublocks, xblocks, dim = self._blocks()
        if not ublocks:
            raise ValueError(
                "the 'kp' method needs at least one Source with an arrival process")
        A, descs = self._build_events(ublocks, xblocks, dim)

        for (i, c, off, h) in xblocks:
            if int(sn.sched[i]) != int(SchedStrategy.INF):
                nsrv = float(sn.nservers[i]) if sn.nservers is not None else 1.0
                if np.isfinite(nsrv):
                    self.linear = False
        if not self.linear:
            from line_solver.api.io.logging import line_warning
            line_warning('solver_fld_kp',
                         'a finite-server station makes the rate functions nonlinear, so '
                         'the covariance is a linear-noise approximation rather than the '
                         'exact second moment; it is exact for infinite-server stations.')

        t0, tend = self.options.timespan
        unbounded = not np.isfinite(tend)
        if not np.isfinite(t0):
            t0 = 0.0
        if unbounded:
            rates = np.asarray(sn.rates, dtype=float)
            finite = rates[np.isfinite(rates) & (rates > 0)]
            slow = float(np.min(finite)) if finite.size else 1.0
            tend = t0 + max(10.0, 30.0 / slow)
            grid = self._breakpoint_grid(t0, tend)
            period = 0.0
            for (_, _, bp, _, cyclic) in (self._sched_entries or []):
                if cyclic:
                    period = max(period, float(bp[-1] - bp[0]))
            if period > 0:
                tend = max(tend, t0 + 10.0 * period)

        q0 = self._initial_state(ublocks, xblocks, dim, t0)
        # The arrival phase is not known at t0, it is drawn from the stationary
        # vector, so its indicator carries covariance diag(u0) - u0 u0'. Starting
        # the covariance at exactly zero would assert a known initial phase and
        # understate the variance until the initial condition has washed out; the
        # error is visible early and decays, and it does not affect steady state.
        Sigma0 = np.zeros((dim, dim))
        for (i, c, off, h) in ublocks:
            u0 = q0[off:off + h]
            Sigma0[off:off + h, off:off + h] = np.diag(u0) - np.outer(u0, u0)
        # Companion seed for the covariance. A caller that carries a DISTRIBUTION
        # across a handoff -- the ENV blend is the case this exists for -- supplies
        # the second moment beside the mean, so the next stage does not restart from
        # a point mass it never had. Same layout as kp_init_sol, and refused rather
        # than ignored when it does not fit.
        init_cov = None
        cfg = getattr(self.options, 'config', None)
        if isinstance(cfg, dict):
            init_cov = cfg.get('init_cov')
        if init_cov is not None and np.size(init_cov) > 0:
            init_cov = np.asarray(init_cov, dtype=float)
            if init_cov.shape != (dim, dim):
                raise ValueError(
                    "config['init_cov'] is %s but the 'kp' state vector of this model has %d "
                    "entries, so the covariance must be (%d, %d)."
                    % (init_cov.shape, dim, dim, dim))
            asym = np.linalg.norm(init_cov - init_cov.T)
            if asym > 1e-6 * max(1.0, float(np.linalg.norm(init_cov))):
                raise ValueError("config['init_cov'] must be symmetric.")
            Sigma0 = init_cov
        z0 = np.concatenate([q0, Sigma0.ravel()])

        tol = self.options.tol
        grid = self._breakpoint_grid(t0, tend)
        max_step = ratemult_max_step(grid) if grid.size > 2 else np.inf
        max_step = min(max_step, (tend - t0) / 10.0)
        odemaxstep = getattr(self.options, 'odemaxstep', None)
        if odemaxstep is not None and np.isfinite(odemaxstep) and odemaxstep > 0:
            max_step = min(max_step, float(odemaxstep))
        sol = solve_ivp(
            lambda t, z: self._augmented_rhs(t, z, A, descs, ublocks, xblocks, dim),
            [t0, tend], z0,
            method='LSODA' if self.options.stiff else 'RK45',
            rtol=tol, atol=tol * 1e-3, dense_output=True, max_step=max_step)
        if not sol.success:
            raise RuntimeError('kp: the augmented ODE failed to integrate (%s)' % sol.message)

        # The period average below is a quadrature, and the integrator's adaptive grid
        # is chosen for ACCURACY OF THE SOLUTION, not for quadrature: at a tight tol it
        # takes large steps through the smooth stretches and the trapezoid over them
        # loses more than the integration gained. Refine the reported grid with a
        # uniform mesh over the averaging window, taken from the dense output.
        t = sol.t
        if unbounded and period > 0:
            lo = max(t0, tend - period)
            refine = np.linspace(lo, tend, 2001)
            # Straddle no discontinuity: the source rate series jumps where the
            # schedule switches, so a trapezoid interval spanning a breakpoint would
            # smear the jump. Bracket each boundary instead.
            bounds = self._breakpoint_grid(lo, tend)
            eps = max(1e-9, 1e-7 * (tend - lo))
            brackets = np.concatenate([bounds - eps, bounds, bounds + eps])
            brackets = brackets[(brackets > lo) & (brackets < tend)]
            t = np.union1d(t, np.union1d(refine, brackets))
        Z = sol.sol(t) if sol.sol is not None else sol.y
        Qtraj = Z[:dim, :]
        Straj = Z[dim:, :].reshape(dim, dim, -1)

        # A cyclic schedule has no fixed point, so a steady-state request is
        # answered by the time average over the last full period of the periodic
        # regime; taking the value at tend would report an arbitrary point of the
        # cycle and make the source and station throughputs disagree.
        period = 0.0
        for (_, _, bp, _, cyclic) in (self._sched_entries or []):
            if cyclic:
                period = max(period, float(bp[-1] - bp[0]))
        avg_window = (max(t0, tend - period), tend) if (unbounded and period > 0) else None

        def summarise(series):
            if avg_window is None:
                return float(series[-1])
            lo, hi = avg_window
            mask = (t >= lo) & (t <= hi)
            if np.count_nonzero(mask) < 2:
                return float(series[-1])
            return float(np.trapezoid(series[mask], t[mask]) / (t[mask][-1] - t[mask][0]))

        QN = np.zeros((M, K)); UN = np.zeros((M, K))
        RN = np.zeros((M, K)); TN = np.zeros((M, K))
        QNt, UNt, TNt, QVart = {}, {}, {}, {}

        for (i, c, off, h) in xblocks:
            qser = Qtraj[off:off + h, :].sum(axis=0)
            var = np.array([float(np.sum(Straj[off:off + h, off:off + h, n]))
                            for n in range(t.size)])
            QNt[(i, c)] = qser
            QVart[(i, c)] = var
            QN[i, c] = summarise(qser)
            nsrv = float(sn.nservers[i]) if sn.nservers is not None else 1.0
            if int(sn.sched[i]) == int(SchedStrategy.INF) or not np.isfinite(nsrv):
                UNt[(i, c)] = qser.copy()
            else:
                UNt[(i, c)] = np.minimum(qser, nsrv) / nsrv
            UN[i, c] = summarise(UNt[(i, c)])

        # throughput: completion rate out of each service block
        for (i, c, off, h) in xblocks:
            ser = np.zeros(t.size)
            for n in range(t.size):
                proc, _ = self._structs_at(t[n])
                D1 = np.asarray(proc[i][c][1])
                x = Qtraj[off:off + h, n]
                sched = int(sn.sched[i])
                nsrv = float(sn.nservers[i]) if sn.nservers is not None else 1.0
                if sched == int(SchedStrategy.INF) or not np.isfinite(nsrv):
                    scale = 1.0
                else:
                    ni = 0.0
                    for (j, l, joff, jh) in xblocks:
                        if j == i:
                            ni += float(np.sum(Qtraj[joff:joff + jh, n]))
                    scale = 1.0 if ni <= nsrv else nsrv / ni
                ser[n] = float(np.sum(D1, axis=1) @ np.maximum(x, 0.0)) * scale
            TNt[(i, c)] = ser
            TN[i, c] = summarise(ser)
            # TN is zero only to the integrator's accuracy; see minnormal.py.
            RN[i, c] = QN[i, c] / TN[i, c] if TN[i, c] > GlobalConstants.Zero else 0.0

        # the source reports its own arrival rate
        for (i, c, off, h) in ublocks:
            arr = np.zeros(t.size)
            for n in range(t.size):
                proc, _ = self._structs_at(t[n])
                D1 = np.asarray(proc[i][c][1])
                arr[n] = float(np.sum(D1, axis=1) @ np.maximum(Qtraj[off:off + h, n], 0.0))
            TNt[(i, c)] = arr
            TN[i, c] = summarise(arr)
            QNt[(i, c)] = np.zeros(t.size)
            UNt[(i, c)] = np.zeros(t.size)
            QVart[(i, c)] = np.zeros(t.size)

        CN = np.sum(RN, axis=0, keepdims=True)
        XN = np.zeros((1, K))
        for c in range(K):
            for (i, cc, off, h) in ublocks:
                if cc == c:
                    XN[0, c] = TN[i, c]

        # AN/WN come from the shared getters, as in every other FLD method. Leaving
        # AN unset made the reported Source ArvR fall back to its throughput, which
        # is the OFFERED-vs-CARRIED convention MATLAB does not use, and showed up as
        # a 100% parity disagreement on that one cell.
        from line_solver.api.sn.getters import sn_get_arvr_from_tput
        from line_solver.api.sn.transforms import sn_get_residt_from_respt
        AN = sn_get_arvr_from_tput(self.sn, TN) if TN is not None else None
        WN = sn_get_residt_from_respt(self.sn, RN, None) if RN is not None else None

        result = FLDResult(
            QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN, AN=AN, WN=WN,
            t=t, QNt=QNt, UNt=UNt, TNt=TNt,
            xvec=Qtraj[:, -1], iterations=1,
            runtime=time.time() - start, method='kp')
        # variance trajectory and full covariance, which no other FLD method has
        result.QVart = QVart
        result.Sigmat = Straj
        result.kpBlocks = {'u': ublocks, 'x': xblocks, 'dim': dim}
        self.runtime = result.runtime
        return result


def solve_kp(sn, options: SolverFLDOptions) -> FLDResult:
    """Module-level entry point matching the other FLD methods."""
    return KoPenderMethod(sn, options).solve()


__all__ = ['KoPenderMethod', 'solve_kp']
