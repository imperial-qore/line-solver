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
from scipy.integrate import solve_ivp, BDF, LSODA, RK45, Radau

if TYPE_CHECKING:
    from ...api.sn import NetworkStruct

from ..options import SolverFLDOptions, FLDResult
from ..utils.phase_type import is_nhpp
from ..utils.ratemult import (
    fluid_interpcols,
    merge_multipliers,
    ratemult_max_step,
    solver_fluid_ratemult,
)
from line_solver.api.sn import SchedStrategy
from ....constants import GlobalConstants
from ..utils.closures import (capacity_closure, share_closure, gps_share, lld_scaling,
                              project_rate)
from ..utils.metrics import fluid_visited_pairs


from ..ode.native_lsoda import NativeLSODA, NativeLSODAStiff

# scipy's compiled integrators, plus the in-tree LSODA of line_solver.lib.lsoda
# under its own names. NOTHING here changes what is selected by default; the two
# native entries are reached only through options.odesolver.
_ODE_CLASSES = {'LSODA': LSODA, 'BDF': BDF, 'RK45': RK45, 'Radau': Radau,
                'LSODA_NATIVE': NativeLSODA, 'LSODA_NATIVE_STIFF': NativeLSODAStiff}


def _ode_class(method):
    """The integrator class for a method name, or the class itself if given one."""
    if isinstance(method, type):
        return method
    return _ODE_CLASSES.get(method, LSODA)


class FluidIntegrationFailure(RuntimeError):
    """The fluid ODE did not integrate over the requested window."""


def _nonnegative_rhs(fun):
    """Port of MATLAB `odenonnegative.m` (lines 29-33).

    ``yp(i) = max(yp(i), 0)`` wherever ``y(i) <= 0``. The STATE handed to the
    drift is deliberately untouched: `_ode_rate_factors` starts from ``rates =
    x`` and clamping its argument would make the drift discontinuous at x = 0
    and collapse the step size. This projects only the returned derivative, so
    the flow cannot push a coordinate that is already at or below zero further
    down -- the tangent-cone projection that defines the fluid limit of a
    nonnegative process.
    """
    def projected(t, y):
        yp = np.asarray(fun(t, y), dtype=float).copy()
        at_floor = np.asarray(y) <= 0.0
        if at_floor.any():
            yp[at_floor] = np.maximum(yp[at_floor], 0.0)
        return yp
    return projected


def _integrate_nonnegative(fun, t0, y0, t_end, method, rtol, atol, max_step,
                           tranpoints=None, guard=None, fixed_point=None):
    """Integrate [t0, t_end] under MATLAB's ``odeset('NonNegative')`` semantics.

    scipy's `solve_ivp` has no NonNegative option, so the stepping loop is
    driven here to reproduce what `ode15s.m` does with one:

    - the derivative is projected on the floor (`_nonnegative_rhs`);
    - a step landing negative by more than the tolerated band has the excursion
      charged as error, as `ode15s.m:678-687` charges ``max(0,-ynew)/AbsTol``
      against `rtol`, and the step cap is halved -- a rejection in effect;
    - such a step is then clipped to zero and the integrator restarted from it,
      which is `ode15s.m:765-771` clipping ``ynew`` and resetting the
      divided-difference table. A fresh solver instance IS that reset, and it
      is what keeps the multistep history from ever straddling the kink at
      x = 0. Clipping without the reset leaves the history inconsistent, which
      is why the projection alone does not fix a stiff LN layer.
    - an excursion INSIDE the band is round-off on a coordinate resting at zero,
      which `errNN <= rtol` accepts: it is clipped in the recorded state and the
      integrator is left running.

    The restart CARRIES THE LAST ACCEPTED STEP, as `ode15s.m` does: its reset
    drops the order to 1 but keeps `h`, it does not re-enter the initial-step
    heuristic. A fresh scipy solver defaults to `select_initial_step`, which on
    a settled trajectory is orders of magnitude below the step just accepted;
    since a coordinate sitting at zero dips negative by round-off on nearly
    every step, that would re-select the tiny initial step every step and the
    window would never advance.

    `tranpoints`, when given, are instants the CALLER wants the trajectory at.
    They are added to the recorded grid from the integrator's own continuous
    extension as each step is accepted, which is what MATLAB's ode solvers do
    when handed an output vector (`options.tranpoints` in
    solver_fluid_iteration.m). It costs no extra steps and no accuracy, and it
    is what lets SolverENV sum its exit average on the sojourn grid without
    LINEARLY interpolating a step grid that has no resolution there.

    `guard`, when given, is called with each accepted state and returns the
    index of a closed class whose conserved population has left the model, or
    -1. It is what stops this loop from running forever: the clipping above is
    exactly MATLAB's `NonNegative`, so a moment-closure drift that leaves the
    simplex is CLAMPED rather than reported here too -- which injects mass,
    halves the cap on every step, and leaves the window never returning. See
    FLUID_CONSERVATION_GUARD (MATLAB) and FluidConservationGuard (JAR).

    `fixed_point`, when given, is the pair ``(residual_tol, slowest_rate)`` that
    decides when a state IS a fixed point rather than merely near one. A state
    below it ENDS THE WINDOW IN CLOSED FORM: see below. It is tested at the
    window's entry AND after every accepted step, because a window that settles
    mid-flight stalls exactly as one entered at its fixed point does.

    Returns (t, y, status, message) with `status` 0 on a completed window, -1
    on a failure and -2 when the guard tripped, never raising, so the caller
    decides what each means.
    """
    y = np.asarray(y0, dtype=float).ravel().copy()
    t0 = float(t0)
    t_end = float(t_end)
    if not (t_end > t0):
        return np.array([t0]), y.reshape(-1, 1), 0, 'empty window'

    cls = _ode_class(method)
    atol_floor = float(np.max(atol)) if np.ndim(atol) else float(atol)
    # The band a negative excursion may sit in without being charged: MATLAB's
    # errNN = max(0,-ynew)/AbsTol tested against rtol, where its AbsTol IS the
    # solver tolerance (solver_fluid_iteration.m sets AbsTol = RelTol = tol).
    # The `atol` handed to the integrator here is deliberately TIGHTER than
    # that, for accuracy, and reading the band off it charges round-off as
    # error: on mqn_singleserver_ps the settled trajectory holds ~99 jobs, dips
    # 1.7e-11 on four consecutive steps, and each charge halves the step cap
    # until the window dies below the spacing of doubles -- 386 windows into a
    # trajectory that had already converged.
    tolerated = rtol * max(atol_floor, rtol)
    f = _nonnegative_rhs(fun)

    cap = max_step
    ts = [t0]
    ys = [y.copy()]
    req = np.empty(0)
    if tranpoints is not None and len(tranpoints) > 0:
        req = np.asarray(tranpoints, dtype=float).ravel()
        req = np.unique(req[(req > t0) & (req < t_end)])
    jreq = 0

    # A FIXED POINT ENDS THE WINDOW IN CLOSED FORM, and this is what keeps a
    # window that has already converged from becoming a window that never
    # returns. The drift here is AUTONOMOUS -- the caller passes `fixed_point`
    # only when it built `ode_func` without a rate schedule -- so f(y*) = 0
    # means y(t) = y* for EVERY later t. The rest of the span is then known
    # exactly, and stepping through it computes a constant the caller already
    # has.
    #
    # THE THRESHOLD IS ROUND-OFF, NOT THE SOLVER TOLERANCE, and the difference
    # is the whole safety of this. The window loop's `drift_displ < drift_tol`
    # (1e-4 by default) says "converged to what the caller asked for", and a
    # state that merely satisfies THAT is still moving: cutting the window there
    # was measured to shift results by 1.8e-5 and broke the three
    # test_lsoda_fluid agreement tests, which hold two integrators to 1e-8. A
    # normalized residual below GlobalConstants.Zero = 1e-14 is a different
    # claim -- the drift is zero to double precision, so no horizon can move the
    # state -- and that claim is what makes skipping the rest of the span exact.
    #
    # WHY THAT WINDOW DOES NOT END ON ITS OWN. scipy's BDF stalls when started
    # AT an equilibrium: on the LN layer of test_LQN_13 it advanced t by 0.011
    # in 20000 steps from a state with |f| = 1.5e-16, and 60 steps covered the
    # whole 1000-unit span once the same state was nudged 1e-6 off it. The
    # layer carries an Immediate() coordinate, i.e. an eigenvalue of exactly
    # -GlobalConstants.Immediate = -1e8 that ODE_ELIMINATE_IMMEDIATE did not
    # fold out, so the step controller is pinned near 1/1e8 while the window
    # runs to 10*iter/min_rate. 272 windows of that layer took 3.0 s between
    # them and the 273rd had not returned after 143 s.
    def _settled(t_at, y_at):
        if fixed_point is None:
            return False
        residual_tol, slowest = fixed_point
        total = float(np.sum(y_at))
        if not (total > 0.0) or not (slowest > 0.0):
            return False
        resid = float(np.abs(np.asarray(f(t_at, y_at), dtype=float)).sum())
        return resid / 2.0 / total / slowest < residual_tol

    def _constant_to_end(y_at):
        """Emit the requested instants and the endpoint at the settled state."""
        j = jreq
        while j < req.size:
            tq = float(req[j])
            if tq > ts[-1]:
                ts.append(tq)
                ys.append(y_at.copy())
            j += 1
        if t_end > ts[-1]:
            ts.append(t_end)
            ys.append(y_at.copy())
        return (np.array(ts), np.array(ys).T, 0,
                'window completed at a fixed point')

    if _settled(t0, y):
        return _constant_to_end(np.maximum(y, 0.0))

    solver = cls(f, t0, y, t_bound=t_end, rtol=rtol, atol=atol, max_step=cap)
    while solver.status == 'running':
        message = solver.step()
        if solver.status == 'failed':
            return np.array(ts), np.array(ys).T, -1, (message or 'ODE step failed')
        tprev = ts[-1]
        tn = float(solver.t)
        yn = np.asarray(solver.y, dtype=float)
        if jreq < req.size and req[jreq] < tn:
            dense = solver.dense_output()
            while jreq < req.size and req[jreq] < tn:
                tq = float(req[jreq])
                if tq > tprev:
                    ts.append(tq)
                    # clipped like every recorded state below: the continuous
                    # extension is of the unprojected step
                    ys.append(np.maximum(np.asarray(dense(tq), dtype=float), 0.0))
                jreq += 1
        if guard is not None:
            bad = guard(yn)
            if bad >= 0:
                ts.append(tn)
                ys.append(np.maximum(yn, 0.0))
                return (np.array(ts), np.array(ys).T, -2,
                        'closed chain %d left the model at t = %g' % (bad, tn))
        # THE ENTRY TEST ABOVE IS TAKEN ONCE, at t0, and a window that reaches
        # the fixed point AFTER its first step is left to grind out the rest of
        # its span at the reciprocal of the FASTEST rate while the span runs to
        # 10*iter/min_rate, set by the SLOWEST. That is the same stall, entered
        # one step later, so it is tested for at every accepted step and ended
        # the same way -- in closed form, at the same round-off threshold, which
        # is what keeps this exact rather than a second convergence criterion.
        # Measured on test_LQN_13 arbitraryMultiplicity, whose T1 multiplicity
        # of 20 drives the single-server FCFS T4 layer into saturation: the
        # entry test alone does not fire there and the window never returns.
        settled_state = np.maximum(yn, 0.0)
        if _settled(tn, settled_state):
            ts.append(tn)
            ys.append(settled_state)
            return _constant_to_end(settled_state)
        below = yn < 0.0
        if not below.any():
            ts.append(tn)
            ys.append(yn.copy())
            cap = max_step
            continue
        clipped = np.where(below, 0.0, yn)
        if float(np.max(-yn[below])) <= tolerated:
            # inside the band `errNN = max(0,-ynew)/AbsTol <= rtol` accepts: this
            # is round-off on a coordinate resting at zero, not a crossing of the
            # boundary. Record the clipped state but leave the integrator alone --
            # a coordinate pinned at zero dips this way on nearly every step, and
            # resetting each time pins the step size at whatever it held when the
            # first dip occurred, so the window never advances.
            ts.append(tn)
            ys.append(clipped)
            cap = max_step
            continue
        # charged as error: shrink the cap, which is what rejecting the step
        # and retrying it smaller accomplishes
        span = cap if np.isfinite(cap) else abs(t_end - t0)
        cap = max(span / 2.0, np.finfo(float).eps * max(1.0, abs(tn)))
        ts.append(tn)
        ys.append(clipped)
        if not (t_end > tn):
            break
        h = getattr(solver, 'h_abs', tn - tprev)
        h = min(float(h), cap, t_end - tn)
        kwargs = {'first_step': h} if h > 0.0 else {}
        solver = cls(f, tn, clipped, t_bound=t_end, rtol=rtol, atol=atol,
                     max_step=cap, **kwargs)
    return np.array(ts), np.array(ys).T, 0, 'window completed'


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
        # verdict of the last _solve_fluid_ode call; 0 is a completed horizon
        self.ode_status = 0
        self.ode_message = ''

    # scheduling policies with a branch in _ode_rate_factors; anything else
    # would keep rates = x, i.e. be integrated as an INFINITE SERVER
    _DRIFT_SCHEDS = (SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS,
                     SchedStrategy.FCFS, SchedStrategy.DPS, SchedStrategy.GPS)

    def _check_scheds(self):
        """Refuse a discipline the drift cannot represent.

        Checked BEFORE the broad handler in solve(), which would otherwise
        swallow the error and return a zero result -- worse than the silent
        infinite-server answer it replaces. On Delay(Z=1) -> Queue(c=1), N=4 the
        fall-through returns 2.0000 against an exact 3.0154.
        """
        sched = self._get_sched()
        for i in range(self.sn.nstations):
            if sched[i] not in self._DRIFT_SCHEDS:
                raise ValueError(
                    "Station %d uses a scheduling policy with no fluid drift branch (%s). "
                    "The closing family would integrate it as an infinite server, which is "
                    "silently wrong. Use options.method='matrix'." % (i, sched[i]))

    def solve(self) -> FLDResult:
        """Solve using closing method - full MATLAB-compatible implementation.

        Returns:
            FLDResult with performance metrics
        """
        start_time = time.time()
        self._check_scheds()

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
                QNt=QNt,
                UNt=UNt,
                TNt=TNt,
                xvec=xvec_it[-1] if xvec_it else x0,
                iterations=self.iterations,
                runtime=time.time() - start_time,
                method='closing'
            )

            if self.ode_status < 0:
                # The horizon was never reached, so `result` reports metrics read
                # off a trajectory that stops short -- at worst off x0 itself.
                # Refuse it. An all-zeros table would be just as false and reads
                # as an answer; see _kb/11-conventions-and-gotchas.md.
                raise FluidIntegrationFailure(
                    "closing: the fluid ODE failed to integrate (%s). The trajectory "
                    "does not reach the requested horizon, so no steady state can be "
                    "read off it." % self.ode_message)

            return result

        except Exception as e:
            if self.options.verbose:
                print(f"Closing method error: {e}")
                import traceback
                traceback.print_exc()
            # A failure must not come back wearing a result. Returning an
            # all-zeros FLDResult here used to turn every exception -- a bad
            # model, an unsupported discipline, a dead integrator -- into a
            # plausible-looking table of zeros, which is the same
            # silent-wrong-answer shape the status propagation above exists to
            # remove, one layer up.
            raise

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
        mu_phi_from_proc = None

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

                # Fallback: derive the per-phase rates from sn.proc. MATLAB fills
                # sn.mu/sn.phi in Station.getServiceRates but the native struct
                # leaves them None, so without this every process collapsed to a
                # single exponential phase at sn.rates and the closing ODE solved
                # the wrong model for Erlang, Coxian, PH, MAP and MMPP2 alike.
                if Mu[i][k] is None:
                    if mu_phi_from_proc is None:
                        mu_phi_from_proc = self._mu_phi_from_proc()
                    mu_ik, phi_ik = mu_phi_from_proc
                    cand = mu_ik.get(i, {}).get(k)
                    if cand is not None and cand.size > 0 and np.all(np.isfinite(cand)) \
                            and np.all(cand > 0):
                        Mu[i][k] = cand
                        Phi[i][k] = phi_ik[i][k]
                        phases[i, k] = cand.size

                # Last resort: a single exponential phase at the mean rate.
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

    def _phase_type_structures(self):
        """Cached (proc_matrix, pie, phases) derived from sn.proc.

        Reuses the extraction the matrix method already relies on, so the two
        FLD paths read one representation of the process.
        """
        cached = getattr(self, '_ph_cache', None)
        if cached is None:
            from ..utils.phase_type import prepare_phase_type_structures
            cached = prepare_phase_type_structures(self.sn)
            self._ph_cache = cached
        return cached

    def _mu_phi_from_proc(self):
        """Per-phase (Mu, Phi) derived from sn.proc, keyed [station][class]."""
        from ..utils.phase_type import extract_mu_phi_from_phase_type
        proc_matrix, _, ph = self._phase_type_structures()
        return extract_mu_phi_from_phase_type(proc_matrix, ph)

    def _get_routing_matrix(self) -> np.ndarray:
        """Get routing probability matrix.

        sn.rt is the pseudo-closed routing matrix: an open class leaving at the
        Sink re-enters at the Source (see Network._refresh_routing). The closing
        ODE has no Sink state, so it relies on that closure -- without it the
        open-class rows are substochastic and the fluid accumulates without
        bound.

        It is indexed by STATEFUL NODE, while the drift below indexes it by
        STATION. The two coincide only when every stateful node is a station: a
        model with a Cache, a Router or a Logger has stateful non-stations, and
        reading the matrix station-major then returns another node pair's
        routing. sn_rt_stations absorbs those nodes exactly (they hold no jobs),
        and returns sn.rt unchanged when there are none.
        """
        if hasattr(self.sn, 'rt') and self.sn.rt is not None:
            from ....api.sn import sn_rt_stations
            rt, _ = sn_rt_stations(self.sn)
            return np.asarray(rt)
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

        An explicitly supplied `options.init_sol` wins, as it does in
        `tbi.py:_resolve_initial_state`: without that, a caller that seeds a
        starting marginal (SolverENV couples its stages by entering each one at
        the previous stage's exit state) is ignored, the transient starts at
        this stage's OWN steady state and stays flat, and the stage-exit average
        collapses to the quasi-stationary blend. A vector of full ODE length is
        taken as is; one of length M*K is a per-(station, class) marginal and is
        laid onto the first phase of each entry.

        Otherwise jobs are distributed evenly across stations where they can be
        served.
        """
        init_sol = getattr(self.options, 'init_sol', None)
        if init_sol is not None:
            x0_in = np.asarray(init_sol, dtype=float).flatten()
            if len(x0_in) == int(np.sum(phases)):
                return x0_in
            if len(x0_in) == M * K:
                x0 = []
                for i in range(M):
                    for k in range(K):
                        n_phases = int(phases[i, k])
                        if n_phases == 0:
                            continue
                        x0.append(x0_in[i * K + k])
                        x0.extend([0.0] * (n_phases - 1))
                return np.array(x0)

        # Get job populations
        N = np.zeros(K)
        if hasattr(self.sn, 'njobs') and self.sn.njobs is not None:
            njobs = np.asarray(self.sn.njobs).flatten()
            for k in range(min(K, len(njobs))):
                # The INFINITY IS THE TEST below, so it must survive the read:
                # folding it to zero here made the EXT branch unreachable and
                # started every Source coordinate at 0 rather than at the 1
                # `solver_fluid.m` assigns an EXT station for an open class
                # ("open job pool"). At 0 that coordinate opens on the kink the
                # NonNegative projection puts at the floor, and on
                # oqn_multichain_cs the BDF step collapsed there once the
                # horizon reached iter_max = 1000 windows.
                N[k] = njobs[k]

        # Build initial state vector
        x0 = []
        assigned = np.zeros(K)
        rt = self._get_routing_matrix()
        sched = self._get_sched()

        # see _kb/06-solver-catalog.md (Fluid: "Python closing method: initial
        # state must match Kic") for why reachability must match Kic exactly
        reachable = np.zeros((M, K), dtype=bool)
        for i in range(M):
            for k in range(K):
                if phases[i, k] > 0:
                    idx = i * K + k
                    if idx < rt.shape[1] and np.sum(rt[:, idx]) > 0:
                        reachable[i, k] = True
        num_reach = reachable.sum(axis=0)

        for i in range(M):
            for k in range(K):
                n_phases = int(phases[i, k])
                if n_phases == 0:
                    continue  # disabled class: no ODE state (Kic==0)

                if reachable[i, k]:
                    if np.isinf(N[k]):
                        # Open class
                        to_assign = 1.0 if sched[i] == SchedStrategy.EXT else 0.0
                    else:
                        # Closed class - distribute evenly over reachable points
                        nst = max(1, int(num_reach[k]))
                        to_assign = N[k] / nst
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

        # Fold out the immediate transitions when asked. They carry
        # GlobalConstants.Immediate, a rate large enough to make the drift stiff
        # and short-lived enough to add no dynamics, so stochastic
        # complementation removes the coordinates without changing the
        # stationary law of the ones that remain. Port of the
        # solver_fluid_odes.m block that calls ode_eliminate_immediate; the flag
        # used to be accepted here and silently ignored.
        from ..immediate import fluid_hide_immediate, ode_eliminate_immediate
        self._immediate_absorb = None
        self._immediate_tput = None
        eventIdx0 = eventIdx
        if fluid_hide_immediate(self.sn, self.options):
            (all_jumps, rateBase, eventIdx, state_map,
             emap, self._immediate_absorb) = ode_eliminate_immediate(
                all_jumps, rateBase, eventIdx, self.sn, self.options)
            # THE COMPLETIONS AN ELIMINATED COORDINATE MAKES ARE NOT LOST. The
            # metrics below read throughputs off the STATE, as x_f * mu_f * phi_f
            # summed over phases, and an eliminated coordinate holds no mass there
            # -- so its completions, which are finite because mu_f is InfRate,
            # would silently vanish and the station would stop balancing against
            # its neighbours. Their total rate is exactly what Emap carries: the
            # composed event that replaced the inflow stands for the original
            # completion too. Built once here, added to TN once the rates are
            # known at the solved state.
            self._immediate_tput = self._build_immediate_tput(
                M, K, q_indices, Kic, eventIdx0, state_map, emap)

            # THE INITIAL POINT HAS TO BE PROJECTED TOO. Once the instantaneous
            # coordinates are complemented away no event moves them any more, so
            # whatever mass the initial condition parked there -- solver_fluid_initsol
            # starts everything in phase 1, and a warm start from an earlier LN
            # iterate parks its own -- would be frozen for the whole integration and
            # lost from its chain. absorb sends it where the eliminated coordinate
            # would have sent it instantaneously. Port of the solver_fluid_iteration.m
            # block; leaving it out is what made an LN reference task with an
            # Immediate() think time report exactly half its throughput.
            if self._immediate_absorb is not None:
                x0 = (np.asarray(self._immediate_absorb, dtype=float).T
                      @ np.asarray(x0, dtype=float).ravel())

        # see _kb/06-solver-catalog.md (Fluid: "Time-varying rates + NHPP
        # transient") -- port of solver_fluid_odes.m lines 104-112
        nominal = np.zeros((M, K))
        for i in range(M):
            for k in range(K):
                if Mu[i][k] is not None and len(Mu[i][k]) > 0:
                    nominal[i, k] = Mu[i][k][0]
        rt_tgrid, rt_Mmat = solver_fluid_ratemult(
            len(rateBase), enabled, q_indices, Kic, nominal, eventIdx, self.options
        )

        # MAPt/PHt: the schedule modulates individual matrix entries, so the
        # multiplier cannot be one scalar per (station,class) as it is for an
        # NHPP. Rebuild the base rates under each segment's matrices and take the
        # ratio to the nominal, which gives an exact per-event factor and needs
        # no event-to-entry bookkeeping. The constant-support rule enforced by
        # the MAPt/PHt constructors is what makes the ratio well defined.
        sch_tgrid, sch_Mmat = self._schedule_event_multiplier(
            M, K, phases, rt, enabled, q_indices, Kic, rateBase, t_end_hint=None)
        if sch_Mmat is not None:
            rt_tgrid, rt_Mmat = merge_multipliers(
                rt_tgrid, rt_Mmat, sch_tgrid, sch_Mmat, len(rateBase))

        # Moment closure, carried in options.config exactly as MATLAB carries it
        # into solver_fluid_odes.m:84-94. Absent means the first-order closure and
        # leaves the plain 'closing' drift bit-identical.
        cfg = getattr(self.options, 'config', None) or {}
        moment_sigma2 = cfg.get('moment_sigma2', None)
        moment_cov = cfg.get('moment_cov', None)
        # sn.lldscaling is read UNCONDITIONALLY, as solver_fluid_odes.m passes it
        # to ode_rates_closing. Gating it on a closure being present left plain
        # method='closing' integrating a load-dependent station at its nominal
        # rate: measured QLen 2/3/4/5 on the alpha = [1 1.7 2.2 2.5 2.6 2.65]
        # sweep of example_fluid_momentclosure, against MATLAB's
        # 1.5882/2.2/2.8667/3.6154. See _kb/06-solver-catalog.md.
        moment_lld = getattr(self.sn, 'lldscaling', None)
        if moment_lld is not None and np.size(moment_lld) == 0:
            moment_lld = None

        # Rate vector at an arbitrary state, kept so the metrics can add back the
        # completions the eliminated coordinates make (see _build_immediate_tput).
        self._immediate_rates_fn = (lambda x: self._ode_rates_closing(
            x, M, K, enabled, q_indices, Kic, nservers, w, sched, rateBase, eventIdx,
            sigma2=moment_sigma2, lld=moment_lld, covblk=moment_cov))

        # Create ODE function
        if rt_Mmat is None:
            def ode_func(t, x):
                rates = self._ode_rates_closing(
                    x, M, K, enabled, q_indices, Kic, nservers, w, sched, rateBase, eventIdx,
                    sigma2=moment_sigma2, lld=moment_lld, covblk=moment_cov
                )
                return all_jumps @ rates
        else:
            def ode_func(t, x):
                rates = self._ode_rates_closing(
                    x, M, K, enabled, q_indices, Kic, nservers, w, sched, rateBase, eventIdx,
                    sigma2=moment_sigma2, lld=moment_lld, covblk=moment_cov
                )
                return all_jumps @ (fluid_interpcols(rt_tgrid, rt_Mmat, t) * rates)

        # Solve ODE. Port of solver_fluid_iteration.m: the horizon is reached by
        # a sequence of windows, each a FRESH integrator run started from the
        # previous window's end state, NOT by one call over the whole span. The
        # restart is load-bearing on a stiff model: it resets the multistep
        # order and history at every window boundary.
        tol = self.options.tol
        t_start, t_end = self.options.timespan

        # slowrate: the slowest per-(station,class) service rate, as
        # solver_fluid.m builds it and ClosingAndStateDepMethodsAnalyzer.java
        # mirrors. Entries at or below CoarseTol, and non-finite ones, are
        # excluded before the minimum, exactly as both do.
        min_rate = np.inf
        for i in range(M):
            for k in range(K):
                if Mu[i][k] is None or not np.size(Mu[i][k]):
                    continue
                cand = np.min(np.asarray(Mu[i][k], dtype=float))
                if np.isfinite(cand) and cand > GlobalConstants.CoarseTol and cand < min_rate:
                    min_rate = cand
        if not np.isfinite(min_rate) or min_rate <= 0:
            min_rate = 1.0

        xvec_it = [x0.copy()]

        # Wall-clock budget (options.timeout, seconds; inf = none), checked at the
        # top of each window as MATLAB checks toc(Tstart) > max_time.
        max_time = getattr(self.options, 'timeout', float('inf'))
        deadline = None
        if np.isfinite(max_time) and max_time > 0:
            deadline = time.monotonic() + max_time
            if time.monotonic() >= deadline:
                # budget gone before the first step, as MATLAB breaks the
                # extension loop before integrating anything
                xvec_it.append(x0.copy())
                self.iterations = 0
                return xvec_it, np.array([x0]), np.array([t_start])

        odemaxstep = getattr(self.options, 'odemaxstep', None)
        if odemaxstep is not None and np.isfinite(odemaxstep):
            max_step = min(odemaxstep, ratemult_max_step(rt_tgrid))
        else:
            max_step = ratemult_max_step(rt_tgrid)
        # MATLAB's stiff slot is @ode15s (SolverOptions.m accurateStiffOdeSolver),
        # a variable-order BDF/NDF method; scipy's analogue of it is BDF. LSODA is
        # a DIFFERENT algorithm -- Adams/BDF auto-switching -- that this path never
        # used in the reference, and whose corrector is the one that fails on an LN
        # layer carrying Immediate rates.
        method = getattr(self.options, 'odesolver', None) \
            or ('BDF' if self.options.stiff else 'RK45')
        iter_max = max(1, int(getattr(self.options, 'iter_max', 200) or 200))

        y = np.asarray(x0, dtype=float).ravel().copy()
        T0 = float(t_start)
        T = T0
        tchunks = [np.array([T0])]
        xchunks = [y.reshape(-1, 1)]
        self.ode_status = 0
        self.ode_message = 'window completed'
        it = 0
        goon = True
        # ARM THE CONSERVATION GUARD ONLY WHERE THE MOMENT CLOSURE IS ACTIVE.
        # minnormal runs its first pass at sigma2 = 0 -- that pass IS the
        # first-order solve -- and only the LATER passes integrate a drift that
        # can leave the simplex. Testing sigma2 keeps the guard off every
        # first-order call, so 'closing' and 'matrix' (which are also the
        # ladder's own fallback) keep identical behaviour and cannot be sent
        # down a fallback by their own watchdog.
        # ONLY WHEN THE STATE IS THE ONE THE PHASES DESCRIBE. hide_immediate
        # folds immediate coordinates out through ode_eliminate_immediate, which
        # SHRINKS the vector, and the block offsets read off phases would then
        # address the wrong coordinates and trip on a sum that was never that
        # class's population. The lengths agreeing is the exact test for that,
        # so the guard stays off rather than guessing.
        _guard = None
        if _closure_active(self.options) and int(np.asarray(phases).sum()) == y.size:
            _guard = _conservation_guard(self.sn, phases)
        _guard_trip = None
        # Early stop on the GEOMETRIC TAIL of the window iteration. The mass
        # moved over ONE window underestimates the distance still left to the
        # fixed point by exactly the tail it drops, r*rho/(1-rho) for a mode
        # contracting by rho per window, so that tail is what is tested; rho is
        # read off the iteration itself, since on a slowly mixing model no rate
        # in the model stands in for the slowest system mode. The drift, zero AT
        # a fixed point, is an independent second bound. Both must hold on two
        # consecutive windows, past the slowest relaxation time.
        earlystop = bool(getattr(self.options, 'config', {}).get('fluid_earlystop', True))
        iter_tol = float(getattr(self.options, 'iter_tol', 0.0) or 0.0)
        drift_tol = max(iter_tol, tol)
        drift_safety = 0.01     # headroom, since rho is estimated, not known
        min_horizon = 10.0 / min_rate
        moved_prev = np.inf
        rho_hist = np.full(3, np.nan)
        drift_below = 0
        from line_solver.api.io import console as _console
        _console.loop('integrating the fluid ODEs over successive time windows')
        try:
            while (np.isfinite(t_end) and T < t_end) or (goon and it < iter_max):
                it += 1
                if deadline is not None and time.monotonic() >= deadline:
                    goon = False
                    break
                # T_i = 10*i/min(slowrate), capped by the requested horizon
                T = min(t_end, abs(10.0 * it / min_rate)) if np.isfinite(t_end) \
                    else abs(10.0 * it / min_rate)
                if not (T > T0):
                    if np.isfinite(t_end) and T >= t_end:
                        goon = False
                        break
                    continue
                # No further cap: MATLAB's odeopt carries AbsTol/RelTol and
                # NonNegative and NO MaxStep, so the integrator is free to stride
                # over a window once the fast modes have settled. Only a rate
                # SCHEDULE (ratemult) and an explicit options.odemaxstep bound it.
                _console.iter_line(it, 'window %d up to t = %g', it, T)
                step_cap = max_step
                y_prev = y.copy()
                # The fixed-point completion is offered only for an AUTONOMOUS
                # drift: with a rate schedule (rt_Mmat) a zero residual now says
                # nothing about the next segment, so the window must be stepped.
                # It rides on `earlystop` because it is the same residual the
                # early stop reads, so switching that off still steps every
                # window in full.
                tw, xw, st, msg = _integrate_nonnegative(
                    ode_func, T0, y, T, method, tol, tol * 1e-3, step_cap,
                    tranpoints=getattr(self.options, 'tranpoints', None),
                    guard=_guard,
                    fixed_point=(GlobalConstants.Zero, min_rate)
                    if (earlystop and rt_Mmat is None
                        and not np.isfinite(t_end)) else None)
                if tw.size > 1:
                    tchunks.append(tw[1:])
                    xchunks.append(xw[:, 1:])
                if st == -2:
                    # The guard tripped. Recorded rather than raised HERE because
                    # this loop sits inside a blanket `except Exception` that
                    # turns any raise into ode_status = -1 and a returned state;
                    # the ladder needs the typed error, so it is raised past it.
                    _guard_trip = msg
                    break
                if st < 0:
                    self.ode_status = -1
                    self.ode_message = msg
                    break
                y = xw[:, -1].copy()
                T0 = T
                total_prev = float(y_prev.sum())
                ratio = float(np.abs(y - y_prev).sum()) / 2.0 / total_prev \
                    if total_prev > 0 else 0.0
                if earlystop and it > 1 and not np.isfinite(t_end) and T >= min_horizon:
                    rho_hist[it % rho_hist.size] = ratio / max(moved_prev, GlobalConstants.Zero)
                    rho = np.nanmax(rho_hist)
                    total = float(y.sum())
                    dy = np.asarray(ode_func(T, y), dtype=float).ravel()
                    drift_displ = float(np.abs(dy).sum()) / 2.0 / total / min_rate \
                        if total > 0 else 0.0
                    # a non-contracting iteration has no tail to sum
                    #
                    # THE 1e-6 GATE IS NOT AN OVERSIGHT, even though it sits below
                    # the integrator's own tol. Relaxing it to "the moved mass
                    # reached the integrator floor, so trust the drift residual
                    # alone" was TRIED and REVERTED: it stops the M/M/1 rho = 0.9
                    # minnormal solve at 7.018088 against the 7.021524680 all four
                    # codebases agree on, and it truncates the statedep response-time
                    # trajectory to t = 60 instead of its 2000. The accuracy of this
                    # loop comes from running the windows, so the stop has to stay
                    # conservative. It is also not what makes a solve hang: see
                    # _kb/06-solver-catalog.md, where the minnormal closure diverges
                    # outright on a bounded multiserver station.
                    if rho < 1.0 and ratio * rho / (1.0 - rho) < drift_safety * drift_tol \
                            and drift_displ < drift_tol:
                        drift_below += 1
                        if drift_below >= 2:
                            goon = False
                    else:
                        drift_below = 0
                moved_prev = ratio
                if np.isfinite(t_end) and T >= t_end:
                    goon = False
        except Exception as e:
            if self.options.verbose:
                print(f"ODE solver error: {e}")
            self.ode_status = -1
            self.ode_message = str(e)

        if _guard_trip is not None:
            from .minnormal import FluidNonHyperbolicError
            raise FluidNonHyperbolicError(
                'The moment-closure drift left the model: %s, and the population '
                'the drift conserves exactly is what moved, so the excursion is a '
                'divergence rather than a solution. Integrating on would not '
                'return. Falling back to a first-order closure.' % _guard_trip)

        t = np.concatenate(tchunks)
        xvec_t = np.concatenate(xchunks, axis=1).T
        xvec_it.append(xvec_t[-1, :].copy())
        self.iterations = it

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

                # DPS and GPS both carry their weights in schedparam, and both
                # branches of ode_rates_closing_factors read w. Setting only DPS
                # here left the GPS mean solve integrating at UNIT weights while
                # the closure that read the fixed point back used the real ones,
                # so flow did not balance: measured QLen [1.494 1.510] against
                # MATLAB's [1.610 1.394] at weights [1 2]. MATLAB assigns both
                # in solver_fluid_odes.m and fluid_moment_terms.m.
                if sched[i] in (SchedStrategy.DPS, SchedStrategy.GPS):
                    w[i, k] = schedparam[i, k]

        return q_indices, Kic, enabled, w

    def _build_ode_system(
        self, M: int, K: int, Mu: List, Phi: List, phases: np.ndarray,
        rt: np.ndarray, enabled: np.ndarray, q_indices: np.ndarray, Kic: np.ndarray,
        ph_structs=None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build ODE jump matrix and rate base.

        Port of MATLAB ode_jumps_new.m and ode_rate_base.m.
        Matches JAR PassageTimeODE.calculateJumps/calculateRateBaseAndEventIdxs.
        """
        state_dim = int(np.sum(Kic))
        if ph_structs is None:
            ph_structs = self._phase_type_structures()
        derived_proc, derived_pie = ph_structs[0], ph_structs[1]

        # Get pie (initial phase probabilities) and proc (PH matrices) per station-class
        pie_all = [[None for _ in range(K)] for _ in range(M)]
        proc_all = [[None for _ in range(K)] for _ in range(M)]
        for j in range(M):
            for l in range(K):
                # Extract pie from sn.pie
                if hasattr(self.sn, 'pie') and self.sn.pie is not None:
                    if j < len(self.sn.pie) and self.sn.pie[j] is not None:
                        if l < len(self.sn.pie[j]) and self.sn.pie[j][l] is not None:
                            pie_all[j][l] = np.asarray(self.sn.pie[j][l]).flatten()
                # The derived structures already normalise every process to a
                # [D0, D1] pair, including the schedule-bearing ones whose raw
                # sn.proc slot holds breakpoints and per-segment matrices and
                # must never be read as phase-type matrices.
                cand = derived_proc.get(j, {}).get(l)
                if cand is not None and len(cand) >= 2 \
                        and np.asarray(cand[0]).shape[0] == Kic[j, l]:
                    proc_all[j][l] = cand

                # Fall back to the pie derived from sn.proc, which is map_pie
                # of the process. Defaulting to phase 0 is right only for an
                # acyclic PH; for a MAP or MMPP2 the arrival-instant phase
                # distribution is not concentrated on phase 1, and using e_1
                # there put the source in the wrong phase mix.
                if pie_all[j][l] is None:
                    cand = derived_pie.get(j, {}).get(l)
                    if cand is not None and np.size(cand) == Kic[j, l] \
                            and np.all(np.isfinite(cand)):
                        pie_all[j][l] = np.asarray(cand, dtype=float).ravel()
                if pie_all[j][l] is None:
                    n_ph = Kic[j, l]
                    if n_ph > 0:
                        pie_all[j][l] = np.zeros(n_ph)
                        pie_all[j][l][0] = 1.0
                    else:
                        pie_all[j][l] = np.array([1.0])

        # Count events
        # Departure events: for each (i,k,ki) -> (j,l,kj) where rt > 0
        n_events = 0
        for i in range(M):
            for k in range(K):
                if enabled[i, k]:
                    for j in range(M):
                        for l in range(K):
                            idx_from = i * K + k
                            idx_to = j * K + l
                            if idx_from < rt.shape[0] and idx_to < rt.shape[1]:
                                if rt[idx_from, idx_to] > 0:
                                    # One event per source phase * destination phase
                                    n_events += Kic[i, k] * Kic[j, l]

        # Internal phase transitions: for each (i,k,ki) -> (i,k,kip) where ki != kip.
        # Every source phase is enumerated, the last one included: bounding ki at
        # Kic-2 is valid only for an acyclic PH, whose last row of D0 has no
        # off-diagonal, and silently drops those transitions for a general MAP or
        # an MMPP2, whose D0 is cyclic. Dropping them cost 14% of the arrival rate
        # on a 2-phase MAP source. Acyclic representations gain only zero-rate
        # candidates here, which the rate>0 guard below discards.
        for i in range(M):
            for k in range(K):
                if enabled[i, k] and Kic[i, k] > 1:
                    for ki in range(Kic[i, k]):
                        for kip in range(Kic[i, k]):
                            if ki != kip:
                                n_events += 1

        if n_events == 0:
            return np.eye(state_dim), np.ones(state_dim), np.arange(state_dim)

        # Build jump matrix and rates
        all_jumps = np.zeros((state_dim, n_events))
        rateBase = np.zeros(n_events)
        eventIdx = np.zeros(n_events, dtype=int)

        event_count = 0

        # Departure events (matching MATLAB ode_jumps_new + ode_rate_base)
        sched = self._get_sched()
        for i in range(M):
            for k in range(K):
                if not enabled[i, k]:
                    continue
                xik = q_indices[i, k]

                for j in range(M):
                    for l in range(K):
                        idx_from = i * K + k
                        idx_to = j * K + l
                        if idx_from >= rt.shape[0] or idx_to >= rt.shape[1]:
                            continue
                        p_route = rt[idx_from, idx_to]
                        if p_route <= 0 or not enabled[j, l]:
                            continue

                        xjl = q_indices[j, l]
                        pie_jl = pie_all[j][l]

                        for ki in range(Kic[i, k]):
                            mu_f = Mu[i][k][ki] if ki < len(Mu[i][k]) else 0
                            phi_f = Phi[i][k][ki] if Phi[i][k] is not None and ki < len(Phi[i][k]) else 1.0

                            for kj in range(Kic[j, l]):
                                pie_kj = pie_jl[kj] if kj < len(pie_jl) else 0.0

                                # Rate = mu * phi * P * pie
                                rateBase[event_count] = mu_f * phi_f * p_route * pie_kj
                                eventIdx[event_count] = xik + ki

                                # Jump: -1 from source phase, +1 to destination phase.
                                # ACCUMULATED, not assigned: a self-transition
                                # (same station, class and phase) has both ends on
                                # one coordinate and must net to zero. Assigning
                                # turned it into a pure +1 birth, so a station with
                                # a routing self-loop created population out of
                                # nothing (p_self=0.7 on cqn_repairmen diverged to
                                # 1e26). Mirrors MATLAB ode_jumps_new.
                                #
                                # EXCEPT AT A SOURCE, WHICH GENERATES RATHER THAN
                                # HOLDS. Its coordinates are the arrival process's
                                # phase indicator, held at total mass one by the EXT
                                # rate rule (the rate factor is 1 identically, so the
                                # arrival is a CONSTANT-rate event and the -1 feeds
                                # back into nothing). Debiting it makes that
                                # indicator fall to zero at t = 1/lambda and then sit
                                # ON the floor, where `_nonnegative_rhs` projects a
                                # derivative of -lambda to 0 and the field has a
                                # kink: `ode15s` crosses it, scipy's BDF collapses
                                # its step there ("Required step size is less than
                                # spacing between numbers") on a trajectory that has
                                # long since settled. Nothing reads the coordinate --
                                # `_compute_metrics_closing` and `minnormal` both
                                # force `QN[EXT, :] = 0` -- so holding it constant
                                # changes no reported quantity, only the smoothness
                                # of the field it lives in.
                                #
                                # THAT ARGUMENT IS ABOUT A ONE-PHASE SOURCE ONLY.
                                # With several phases the rate factor is NOT 1: the
                                # EXT rule pins phase 1 at 1 - sum(others) and reads
                                # the rest off the state, so the phase MIX decides
                                # the arrival rate. Dropping the debit then lets the
                                # returning flow credit pie into phases that never
                                # pay for their arrivals, and the mix drifts: on the
                                # 2-phase MAP D0=[-5 1; 2 -4], D1=[3 1; 1 1] the
                                # source delivered 2.7586 against the exact 3.2,
                                # while MATLAB, which always debits, delivers 3.2.
                                if sched[i] != SchedStrategy.EXT or Kic[i, k] > 1:
                                    all_jumps[xik + ki, event_count] -= 1
                                all_jumps[xjl + kj, event_count] += 1
                                event_count += 1

        # Everything emitted so far is a service COMPLETION; what follows is an
        # intra-PH phase change. The boundary is what lets a caller say which
        # (station,class) an event is a completion of.
        self._n_departure_events = event_count

        # Internal phase transitions (matching MATLAB ode_rate_base using D0 matrix)
        for i in range(M):
            for k in range(K):
                if not enabled[i, k] or Kic[i, k] <= 1:
                    continue
                xik = q_indices[i, k]

                for ki in range(Kic[i, k]):
                    for kip in range(Kic[i, k]):
                        if ki == kip:
                            continue

                        # Use D0 matrix entry if available, else mu*(1-phi) for sequential
                        if proc_all[i][k] is not None:
                            D0 = np.asarray(proc_all[i][k][0])
                            rate = D0[ki, kip] if ki < D0.shape[0] and kip < D0.shape[1] else 0.0
                        else:
                            # Fallback: sequential transitions only
                            if kip == ki + 1:
                                mu_f = Mu[i][k][ki] if ki < len(Mu[i][k]) else 0
                                phi_f = Phi[i][k][ki] if Phi[i][k] is not None and ki < len(Phi[i][k]) else 1.0
                                rate = mu_f * (1.0 - phi_f)
                            else:
                                rate = 0.0

                        if rate > 0:
                            rateBase[event_count] = rate
                            eventIdx[event_count] = xik + ki
                            all_jumps[xik + ki, event_count] = -1
                            all_jumps[xik + kip, event_count] = 1
                            event_count += 1

        # Trim to actual event count
        all_jumps = all_jumps[:, :event_count]
        rateBase = rateBase[:event_count]
        eventIdx = eventIdx[:event_count].astype(int)

        return all_jumps, rateBase, eventIdx

    def _schedule_event_multiplier(self, M, K, phases, rt, enabled, q_indices,
                                   Kic, rate_base_nom, t_end_hint=None):
        """Per-event multiplier trajectory for the MAPt/PHt stations, if any.

        Returns (tgrid, Mmat) with Mmat of shape (nEvents x len(tgrid)), or
        (None, None) when the model has no schedule-bearing matrix process.
        """
        from ..utils.phase_type import (is_mapt, is_pht, schedule_segments,
                                        extract_mu_phi_from_phase_type,
                                        prepare_phase_type_structures)

        # Steady state uses the time-averaged nominal, matching the NHPP
        # convention: a cyclic schedule's stationary regime IS its time average,
        # and a multiplier grid clamped past its last sample would instead freeze
        # the final segment. The schedule is honoured whenever the caller asks
        # for a bounded horizon, which is what getTranAvg sets up.
        if not np.isfinite(self.options.timespan[1]):
            return None, None

        entries = []
        for i in range(M):
            for k in range(K):
                if not enabled[i, k]:
                    continue
                kind = 'MAPt' if is_mapt(self.sn, i, k) else ('PHt' if is_pht(self.sn, i, k) else None)
                if kind is None:
                    continue
                slot = self.sn.proc[i][k]
                bp, pairs, cyclic = schedule_segments(slot, kind)
                entries.append((i, k, bp, len(pairs), cyclic))
        if not entries:
            return None, None

        t0, tend = self.options.timespan
        if not np.isfinite(t0):
            t0 = 0.0
        if not np.isfinite(tend):
            tend = t_end_hint
        if tend is None or not np.isfinite(tend):
            # Cover a few periods so a cyclic schedule is represented rather
            # than clamped after its first segment.
            period = max(float(e[2][-1] - e[2][0]) for e in entries)
            tend = t0 + 3.0 * period if period > 0 else t0 + 1.0

        # Union of every schedule's segment boundaries over the horizon.
        bounds = [t0, tend]
        for (_, _, bp, _, cyclic) in entries:
            period = float(bp[-1] - bp[0])
            if cyclic and period > 0:
                kmax = int(np.ceil((tend - t0) / period)) + 2
                grid = np.concatenate([bp + kk * period for kk in range(-1, kmax + 1)])
            else:
                grid = bp
            bounds.extend(grid[(grid > t0) & (grid < tend)].tolist())
        bounds = np.unique(np.asarray(bounds, dtype=float))

        rate_base_nom = np.asarray(rate_base_nom, dtype=float)
        neps = max(1e-9, 1e-6 * (tend - t0))
        nb = bounds.size - 1
        seg_t = np.zeros(2 * nb)
        Mmat = np.ones((rate_base_nom.size, 2 * nb))
        for b in range(nb):
            a, z = bounds[b], bounds[b + 1]
            mid = 0.5 * (a + z)
            segment = {}
            for (i, k, bp, nseg, cyclic) in entries:
                period = float(bp[-1] - bp[0])
                off = mid - bp[0]
                if cyclic and period > 0:
                    off = off % period
                elif off < 0.0 or off >= period:
                    segment[(i, k)] = None
                    continue
                idx = int(np.searchsorted(bp[1:], bp[0] + off, side='right'))
                segment[(i, k)] = min(idx, nseg - 1)
            if any(v is None for v in segment.values()):
                # Past a non-cyclic horizon the process is silent.
                col = np.zeros(rate_base_nom.size)
            else:
                ph_seg = prepare_phase_type_structures(self.sn, segment=segment)
                Mu_s, Phi_s = extract_mu_phi_from_phase_type(ph_seg[0], ph_seg[2])
                Mu_l = [[Mu_s.get(i, {}).get(k) for k in range(K)] for i in range(M)]
                Phi_l = [[Phi_s.get(i, {}).get(k) for k in range(K)] for i in range(M)]
                rb = self._build_ode_system(M, K, Mu_l, Phi_l, phases, rt,
                                            enabled, q_indices, Kic,
                                            ph_structs=ph_seg)[1]
                if rb.size != rate_base_nom.size:
                    raise ValueError(
                        'the MAPt/PHt schedule changed the event count between '
                        'segments (%d vs %d); the matrix sparsity pattern must be '
                        'constant across segments' % (rb.size, rate_base_nom.size))
                col = np.ones(rate_base_nom.size)
                nz = rate_base_nom != 0.0
                col[nz] = rb[nz] / rate_base_nom[nz]
            seg_t[2 * b] = a
            seg_t[2 * b + 1] = max(a + neps, z - neps)
            Mmat[:, 2 * b] = col
            Mmat[:, 2 * b + 1] = col
        return seg_t, Mmat

    def _ode_rates_closing(
        self, x: np.ndarray, M: int, K: int, enabled: np.ndarray,
        q_indices: np.ndarray, Kic: np.ndarray, nservers: np.ndarray,
        w: np.ndarray, sched: List, rateBase: np.ndarray, eventIdx: np.ndarray,
        sigma2=None, lld=None, covblk=None
    ) -> np.ndarray:
        """Event-rate vector of the closing method.

        Port of ode_rates_closing.m: the per-coordinate service share gathered
        over the event index set and scaled by the rate base.
        """
        rates = self._ode_rate_factors(
            x, M, K, enabled, q_indices, Kic, nservers, w, sched,
            sigma2=sigma2, lld=lld, covblk=covblk)
        return rates[eventIdx] * rateBase

    def _ode_rate_factors(
        self, x: np.ndarray, M: int, K: int, enabled: np.ndarray,
        q_indices: np.ndarray, Kic: np.ndarray, nservers: np.ndarray,
        w: np.ndarray, sched: List, sigma2=None, lld=None, covblk=None
    ) -> np.ndarray:
        """Per-coordinate service share, before event indexing and rate base.

        Port of ode_rates_closing_factors.m. Kept separate from
        ``_ode_rates_closing`` so the moment-closure methods read the same
        service shares the ODE integrated.

        ``sigma2`` is the per-station closure variance (None or zero: the
        first-order closure). ``lld`` is ``sn.lldscaling``. ``covblk`` is a list
        of per-station coordinate covariance blocks, which closes the PS/DPS
        capacity-share ratio. All three leave the legacy code path
        bit-identical when absent, so the untouched methods are unaffected.
        """
        # No projection of x here: MATLAB ode_rates_closing_factors.m starts from
        # `rates = x` and leaves nonnegativity to odeset('NonNegative'), which acts
        # on the ACCEPTED step. Clamping the drift's argument instead makes it
        # discontinuous at x=0 and collapses the integrator step size.
        rates = np.array(x, dtype=float, copy=True)
        if sigma2 is None:
            sigma2 = np.zeros(M)
        gaussian = bool(np.any(np.asarray(sigma2) > 0))
        has_lld = lld is not None and len(lld) > 0

        for i in range(M):
            sched_i = sched[i]

            lldrow = None
            if has_lld and i < len(lld):
                row = np.asarray(lld[i], dtype=float).ravel()
                if row.size > 0 and np.any(np.abs(row - 1.0) > 1e-14):
                    lldrow = row

            Ci = None
            if covblk is not None and i < len(covblk):
                Ci = covblk[i]

            if sched_i == SchedStrategy.INF:
                # without load dependence each job is served at its own rate and
                # the share is the identity; alpha(n_i) scales the whole station
                if lldrow is not None:
                    lo = q_indices[i, 0]
                    hi = q_indices[i, K - 1] + Kic[i, K - 1] if K > 0 else lo
                    ni = float(x[lo:hi].sum())
                    if ni > 0:
                        h, _ = capacity_closure(ni, nservers[i], sigma2[i], lldrow, True)
                        rates[lo:hi] = x[lo:hi] / ni * h

            elif sched_i == SchedStrategy.EXT:
                # keep total mass 1 into the source for all classes at all times
                for k in range(K):
                    if enabled[i, k]:
                        idx_ini = q_indices[i, k]
                        idx_end = idx_ini + Kic[i, k]
                        if idx_ini < len(rates):
                            rates[idx_ini] = 1.0 - x[idx_ini + 1:idx_end].sum()

            elif sched_i in (SchedStrategy.PS, SchedStrategy.FCFS):
                lo = q_indices[i, 0]
                hi = q_indices[i, K - 1] + Kic[i, K - 1] if K > 0 else lo
                ni = float(x[lo:hi].sum())
                if (gaussian or lldrow is not None) and ni > 0:
                    h, dh = capacity_closure(ni, nservers[i], sigma2[i], lldrow, False)
                    if Ci is None:
                        rates[lo:hi] = x[lo:hi] / ni * h
                    else:
                        # THE SHARE AND THE CAPACITY ARE CLOSED JOINTLY. What the
                        # station clears is S_j*psi(N), and both factors move with
                        # N, so the product needs Cov(S_j,N)*psi'(n) on top of the
                        # two separate closures; see share_closure. With unit
                        # weights this is the DPS branch below.
                        s, cn = share_closure(x[lo:hi], np.ones(hi - lo), Ci,
                                              want_cov=True)
                        rates[lo:hi] = project_rate(s * h + dh * cn, x[lo:hi],
                                                    lldrow is None, h)
                elif ni > nservers[i]:  # case min = ni handled by rates = x
                    rates[lo:hi] = x[lo:hi] / ni * nservers[i]

            elif sched_i == SchedStrategy.DPS:
                # DPS is PS with a weighted share: the class-k coordinates get
                # w_k*x/nbar of the station capacity psi(xi) instead of x/xi.
                #
                # The denominator used to carry an ADDITIVE mean(w) term as a
                # divide-by-zero guard. It never cancels, so the class shares
                # summed to 1 - mean(w)/nbar instead of 1 and utilization was
                # depressed by exactly that factor (24.9% measured with weights
                # [1 4]). The capacity was also the full c rather than
                # psi(xi) = min(xi,c)*alpha(xi). Both now match PS/FCFS.
                #
                # CONSISTENCY CHECK: with equal weights this collapses to
                # x/xi*psi(xi), the PS branch verbatim, so equal-weight DPS must
                # be bit-identical to PS.
                #
                # The share itself is a RATIO, so evaluating it at the mean is a
                # separate closure from the min(): share_closure corrects it at
                # second order when a covariance block is supplied.
                w_i = w[i, :].copy()
                w_sum = np.sum(w_i)
                if w_sum > 0:
                    w_i = w_i / w_sum

                lo = q_indices[i, 0]
                hi = q_indices[i, K - 1] + Kic[i, K - 1] if K > 0 else lo
                wv = np.zeros(hi - lo)
                for k in range(K):
                    if enabled[i, k]:
                        a0 = q_indices[i, k] - lo
                        wv[a0:a0 + Kic[i, k]] = w_i[k]

                xi = float(x[lo:hi].sum())
                if xi > 0 and float(wv @ x[lo:hi]) > 0:
                    psi, dpsi = capacity_closure(xi, nservers[i], sigma2[i], lldrow, False)
                    s, cn = share_closure(x[lo:hi], wv, Ci, want_cov=True)
                    rates[lo:hi] = project_rate(s * psi + dpsi * cn, x[lo:hi],
                                                lldrow is None, psi)

            elif sched_i == SchedStrategy.GPS:
                # GPS splits the server by WEIGHT among the BACKLOGGED classes,
                # then equally among that class's own jobs. The share is a
                # function of the backlog indicator, so gps_share closes it over
                # the 2^K patterns using P(N_k >= 1). No capacity term
                # multiplies it: GPS is single-server and the indicator already
                # carries the idle server, so the shares sum to
                # 1 - P(station empty) by design.
                if nservers[i] > 1:
                    raise ValueError('Multi-server GPS stations are not supported yet.')
                lo = q_indices[i, 0]
                hi = q_indices[i, K - 1] + Kic[i, K - 1] if K > 0 else lo
                xk = np.zeros(K)
                vk = np.zeros(K)
                blk = [None] * K
                for k in range(K):
                    if enabled[i, k]:
                        a0 = q_indices[i, k]
                        blk[k] = slice(a0, a0 + Kic[i, k])
                        xk[k] = float(x[blk[k]].sum())
                        if Ci is not None:
                            loc = slice(a0 - lo, a0 - lo + Kic[i, k])
                            vk[k] = max(0.0, float(np.sum(np.asarray(Ci)[loc, loc])))
                sk = gps_share(xk, w[i, :], vk)
                a = 1.0
                if lldrow is not None:
                    aa, _ = lld_scaling(lldrow, float(x[lo:hi].sum()))
                    a = float(aa[0])
                for k in range(K):
                    if enabled[i, k] and xk[k] > 0:
                        rates[blk[k]] = x[blk[k]] / xk[k] * sk[k] * a

            else:
                # A station without a case above would keep rates = x, i.e. it
                # would be integrated as an INFINITE SERVER, and the answer
                # would be wrong with no warning: on Delay(Z=1) -> Queue(c=1),
                # N=4 the fall-through returns 2.0000 against an exact 3.0154.
                # Refuse instead. SIRO/LCFS/LCFSPR are declared in the solver
                # feature set because the matrix method handles them as PS,
                # which is the right aggregate for any work-conserving
                # discipline; this drift is the closing family only.
                raise ValueError(
                    "Station %d uses a scheduling policy with no fluid drift branch "
                    "(%s). The closing family integrates such a station as an infinite "
                    "server, which is silently wrong. Use options.method='matrix'."
                    % (i, sched_i))

        return rates

    def _add_immediate_tput(self, TN, x, M, K):
        """Add the completions the eliminated coordinates make to TN, in place."""
        C = getattr(self, '_immediate_tput', None)
        fn = getattr(self, '_immediate_rates_fn', None)
        if C is None or fn is None:
            return TN
        r = np.asarray(fn(np.asarray(x, dtype=float).ravel()), dtype=float).ravel()
        if r.size != C.shape[1]:
            return TN
        extra = (C @ r).reshape(M, K)
        return TN + extra

    def _build_immediate_tput(self, M, K, q_indices, Kic, eventIdx0, state_map, emap):
        """Per-(station,class) completion rate carried by the ELIMINATED coordinates.

        Returns a matrix ``C`` of shape ``(M*K, n_reduced_events)`` such that
        ``C @ r`` is the throughput those coordinates contribute, given the reduced
        rate vector ``r`` at the solved state. ``Emap[e, o]`` is the expected number
        of times the original event ``o`` fires per firing of the reduced event
        ``e``, so summing it over the original COMPLETIONS sourced at an eliminated
        coordinate of ``(i,k)`` is exactly that station's missing flow.
        """
        n_events0 = int(np.size(eventIdx0))
        n_dep = int(getattr(self, '_n_departure_events', n_events0))
        nstate = int(np.sum(Kic))
        coord_station = np.full(nstate, -1, dtype=int)
        coord_class = np.full(nstate, -1, dtype=int)
        for i in range(M):
            for k in range(K):
                lo = int(q_indices[i, k])
                for f in range(int(Kic[i, k])):
                    coord_station[lo + f] = i
                    coord_class[lo + f] = k
        is_kept = np.zeros(nstate, dtype=bool)
        is_kept[np.asarray(state_map, dtype=int)] = True

        C = np.zeros((M * K, emap.shape[0]))
        for o in range(min(n_dep, n_events0)):
            c = int(eventIdx0[o])
            if c < 0 or c >= nstate or is_kept[c]:
                continue
            i, k = coord_station[c], coord_class[c]
            if i < 0 or k < 0:
                continue
            C[i * K + k, :] += emap[:, o]
        return C

    def _compute_metrics_closing(
        self, xvec_it: List, xvec_t: np.ndarray, t: np.ndarray,
        M: int, K: int, Mu: List, Phi: List, phases: np.ndarray,
        nservers: np.ndarray, sched: List, schedparam: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict, Dict, Dict]:
        """Compute performance metrics from fluid solution.

        Port of solver_fluid_closing.m. Returns the steady-state (final-state)
        metrics together with the transient trajectories QNt/UNt/TNt, keyed by
        (station, class), obtained by evaluating the same state-to-metrics map
        at every integration time point.
        """
        # Get final state
        x_final = xvec_it[-1] if xvec_it else np.zeros(int(np.sum(phases)))
        x_final = np.maximum(x_final, 0.0)
        q_indices, Kic, _, _ = self._build_ode_indices(M, K, Mu, phases, sched, schedparam)

        QN, UN, RN, TN = self._metrics_from_state(
            x_final, M, K, Mu, Phi, nservers, sched, schedparam, q_indices, Kic)
        TN = self._add_immediate_tput(TN, x_final, M, K)
        # RN follows TN, so it is recomputed rather than left over the pre-correction
        # throughput; a station whose only completions were the instantaneous ones
        # would otherwise report an infinite response time.
        with np.errstate(divide='ignore', invalid='ignore'):
            RN = np.where(TN > GlobalConstants.Zero, QN / np.where(TN > 0, TN, 1.0), RN)

        # Transient trajectories: the same map applied along the integrated path.
        QNt = {}
        UNt = {}
        TNt = {}
        if xvec_t is not None and np.ndim(xvec_t) == 2 and xvec_t.shape[0] > 1:
            nt = xvec_t.shape[0]
            Qtr = np.zeros((nt, M, K))
            Utr = np.zeros((nt, M, K))
            Ttr = np.zeros((nt, M, K))
            for it in range(nt):
                xq, xu, _, xt = self._metrics_from_state(
                    np.maximum(xvec_t[it, :], 0.0), M, K, Mu, Phi, nservers,
                    sched, schedparam, q_indices, Kic)
                Qtr[it] = xq
                Utr[it] = xu
                Ttr[it] = xt
            for i in range(M):
                for k in range(K):
                    QNt[(i, k)] = Qtr[:, i, k]
                    UNt[(i, k)] = Utr[:, i, k]
                    TNt[(i, k)] = Ttr[:, i, k]

        return QN, UN, RN, TN, QNt, UNt, TNt

    def _metrics_from_state(
        self, x_final: np.ndarray, M: int, K: int, Mu: List, Phi: List,
        nservers: np.ndarray, sched: List, schedparam: np.ndarray,
        q_indices: np.ndarray, Kic: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Map one fluid state vector to (QN, UN, RN, TN)."""
        # Compute queue lengths per station/class
        QN = np.zeros((M, K))

        for i in range(M):
            for k in range(K):
                shift = q_indices[i, k]
                n_phases = Kic[i, k]
                if n_phases > 0 and shift + n_phases <= len(x_final):
                    QN[i, k] = np.sum(x_final[shift:shift + n_phases])
            if sched[i] == SchedStrategy.EXT:
                # A Source holds no jobs. Its ODE coordinates are the arrival
                # process's phase indicator, held at total mass one by the EXT
                # rate rule, and the station-level routing folds every departure
                # of an open class back onto it -- so from an initial state
                # above the fixed point the coordinate accumulates the mass that
                # LEFT the system. Reporting that as a queue length is wrong on
                # its own and, when a caller carries the marginal back in
                # (SolverENV), feeds it to the next solve as jobs at the Source.
                QN[i, :] = 0.0

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

        # Compute response times using Little's Law. A CLASS THE MODEL NEVER
        # ROUTES HERE HAS NO RESPONSE TIME, and neither QN nor TN says so by its
        # size: both hold a remnant of the initial state the integrator was still
        # draining. See fluid_visited_pairs.
        RN = np.zeros((M, K))
        visited = fluid_visited_pairs(self.sn, M, K)
        with np.errstate(divide='ignore', invalid='ignore'):
            for i in range(M):
                for k in range(K):
                    if visited[i, k] and TN[i, k] > 1e-10:
                        RN[i, k] = QN[i, k] / TN[i, k]
                    else:
                        RN[i, k] = 0.0

        return QN, UN, RN, TN

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


def _closure_active(options):
    """True when the moment closure is what is being integrated."""
    cfg = getattr(options, 'config', None) or {}
    if isinstance(cfg, dict):
        sigma2 = cfg.get('moment_sigma2')
    else:
        sigma2 = getattr(cfg, 'moment_sigma2', None)
    if sigma2 is None:
        return False
    return bool(np.any(np.asarray(sigma2, dtype=float) != 0.0))


def _chain_partition(sn, K):
    """The classes of each chain, as a list of index arrays over 0..K-1.

    ``sn.chains`` comes in either the index format (1D, ``chains[class]`` is the
    chain) or the membership format (2D, ``chains[chain, class]``); see
    _kb/04-networkstruct.md. A struct declaring neither degrades to one chain
    per class, which is the safe reading rather than a guess: with no chain map
    there is no class switching to merge classes, so each class IS its own
    conserved unit.
    """
    raw = getattr(sn, 'chains', None)
    if raw is None:
        return [np.array([k]) for k in range(K)]
    raw = np.asarray(raw)
    if raw.size == 0:
        return [np.array([k]) for k in range(K)]
    if raw.ndim == 1:
        return [np.where(raw == c)[0] for c in np.unique(raw)]
    return [np.where(raw[c, :] != 0)[0] for c in range(raw.shape[0])]


def _conservation_guard(sn, phases, tol=0.1):
    """Detect a moment-closure trajectory that has left the model.

    THE TEST IS AN EXACT INVARIANT, not a heuristic bound on time or magnitude.
    The drift conserves the population of every CLOSED CHAIN exactly, so any
    deviation is the `NonNegative` clipping injecting mass and nothing else.
    ``tol`` is therefore a generous fraction of that population rather than a
    numerical tolerance: the integrator's own error is ~1e-4 relative, while
    the documented excursion reaches 5.2e4 against a true population of 0.05. A
    closed model whose population has moved by ``tol`` is no longer solving the
    model, whatever it is converging to.

    THE CHAIN IS THE CONSERVED UNIT, NOT THE CLASS, and the difference is the
    whole correctness of this check. ``sn.njobs[k]`` is the population class k
    STARTS with; class switching then moves jobs between the classes of one
    chain, so only the chain total is invariant. Watching classes instead
    condemns every class-switching model out of hand -- measured on
    ``cqn_twoclass_hyperl`` (313 of 447 accepted states), on ``init_state_ps``
    (286 of 310) and on every one of the 162 fluid layers an LQN builds under
    the ``srvn.cs`` encoding, where the chain sum never moved at all. A cache
    model is the same story with the hit/miss classes.

    A wall-clock budget would have caught the same thing and was rejected: it
    makes the answer depend on how busy the host is, so the same model would
    fall back on one machine and not on another. This invariant is
    deterministic.

    Returns a callable mapping a state vector to the index of the offending
    closed chain, or -1.
    """
    phases = np.asarray(phases)
    M, K = phases.shape
    njobs = np.asarray(getattr(sn, 'njobs', np.array([]))).ravel()
    watched = []
    for c, members in enumerate(_chain_partition(sn, K)):
        total = 0.0
        for k in members:
            total += njobs[k] if k < njobs.size else np.inf
        if not np.isfinite(total) or total <= 0:
            continue  # open, or absent: no conserved population to check
        cols = []
        for k in members:
            for i in range(M):
                ph = int(phases[i, k])
                if ph <= 0:
                    continue
                shift = int(phases[:i, :].sum() + phases[i, :k].sum())
                cols.extend(range(shift, shift + ph))
        if cols:
            watched.append((c, np.sort(np.asarray(cols, dtype=int)), float(total)))

    def guard(state):
        if not watched:
            return -1
        x = np.asarray(state, dtype=float).ravel()
        for c, cols, target in watched:
            if cols[-1] >= x.size:
                continue  # a caller integrating a different state vector
            if abs(float(x[cols].sum()) - target) > tol * max(1.0, target):
                return c
        return -1

    return guard
