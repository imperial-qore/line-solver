"""
Second-order moment-closure fluid analysis, backing ``options.method='minnormal'``
and ``options.method='refined'``.

The default fluid methods close the moment hierarchy at first order: the drift
of the mean depends on ``E[min(X_i,c_i)]``, which they replace by
``min(E[X_i],c_i)``. No second moment ever enters, so no variance is produced
and the mean itself is biased wherever ``min()`` is not locally linear. This
method reinstates the second moment.

``minnormal``
    Min-normal closure (Guenther, Stefanek, Bradley). The drift uses
    ``E[min(X_i,c_i)]`` under a normal marginal whose variance is produced by
    the covariance equation, so mean and covariance are solved self-consistently
    by fixed-point iteration. This corrects the mean, most visibly near rho = 1
    where the first-order closure is worst.

``refined``
    Refined mean field (Gast, POMACS 2017). The O(1/N) term of the expansion of
    the true stationary mean about the MEAN-FIELD fixed point, obtained by
    carrying the stationary covariance through the Hessian of the drift. It is
    an expansion about the mean-field point, NOT about the Gaussian one: adding
    it to the ``minnormal`` fixed point would count the same O(1/N) term twice,
    since the Gaussian closure already resums it. So the base point is recomputed
    with the first-order closure, while the Hessian and the Jacobian are taken
    from the smooth Gaussian drift -- the hard ``min`` being only piecewise
    linear and, at saturation, kinked exactly at the fixed point.

All performance measures are read back from the same event representation that
defines the drift, so throughputs balance flow at the fixed point under
whichever closure was used.

Port of ``solver_fluid_moments.m``, ``fluid_moment_terms.m``,
``fluid_drift_jacobian.m``, ``fluid_lyapunov.m`` and
``fluid_refine_meanfield.m``.
"""

import time
from typing import List, Optional

import numpy as np
from scipy.linalg import orth, solve_sylvester

from line_solver.api.sn import SchedStrategy
from line_solver.constants import GlobalConstants
from ..options import SolverFLDOptions, FLDResult
from .closing import ClosingMethod, FluidIntegrationFailure
from ..utils.closures import capacity_closure, share_closure, gps_share, lld_scaling
from ..utils.metrics import fluid_visited_pairs


class FluidNonHyperbolicError(ValueError):
    """The fluid fixed point has no stationary linear noise approximation.

    Typed so the caller can tell it apart from any other failure:
    fluid_minnormal_applicable cannot see a non-hyperbolic fixed point in
    advance (it exists only once the mean is solved), so SolverFLD switches a
    RESOLVED 'minnormal' to the first-order method on this class alone. An
    explicit method='minnormal' still propagates. MATLAB twin: the
    'LINE:FluidNonHyperbolic' identifier raised by fluid_lyapunov.m.
    """


class FluidIntegrationError(FluidIntegrationFailure):
    """The mean-field ODE did not integrate, so there is no fixed point to close at.

    Typed apart from FluidNonHyperbolicError, which is a property of a fixed
    point that WAS reached. This one says no fixed point was reached at all, and
    it must not be demoted to the first-order method behind the caller's back:
    the trajectory the integrator returns in this case is the initial guess.
    """


def fluid_lyapunov(A, Qdiff, D, tol=None):
    """Stationary covariance of the linear noise approximation.

    Around a fixed point of the fluid drift the fluctuation obeys
    ``dZ = A*Z*dt + sqrt(Qdiff)*dW``, whose stationary covariance solves

        A*Sigma + Sigma*A' + Qdiff = 0,   Qdiff = D*diag(r)*D'

    ``A`` is singular whenever the model conserves population: every closed
    class contributes a left null vector. It does have a unique solution on the
    reachable subspace ``range(D)``: the state can only move along jump
    directions, so the fluctuation lives there and nowhere else. Both ``A`` and
    ``Qdiff`` map into ``range(D)`` as well, so restricting to an orthonormal
    basis of it is an EXACT reduction, not an approximation.
    """
    if tol is None:
        tol = np.sqrt(np.finfo(float).eps)
    n = A.shape[0]
    V = orth(D)
    if V.size == 0:
        return np.zeros((n, n))

    Ar = V.T @ A @ V
    Qr = V.T @ Qdiff @ V
    Qr = (Qr + Qr.T) / 2

    ev = np.linalg.eigvals(Ar)
    maxRe = float(np.max(ev.real))
    if not (maxRe < -tol):
        raise FluidNonHyperbolicError(
            "The fluid fixed point is not exponentially stable on the reachable subspace "
            "(largest Jacobian eigenvalue has real part %g), so the linear noise "
            "approximation has no stationary covariance. This happens at an unstable model "
            "or at a drift kink; use options.method='closing' for the mean only." % maxRe)

    W = solve_sylvester(Ar, Ar.T, -Qr)
    W = (W + W.T) / 2
    Sigma = V @ W @ V.T
    return (Sigma + Sigma.T) / 2


def fluid_refine_meanfield(x, sigma2, Sigma, drift_fcn, jac_fcn, D, covblk=None,
                           epsrel=1e-4):
    """Refined mean field correction of a fluid fixed point (Gast, POMACS 2017).

    The mean-field fixed point x* is the leading term of an expansion of the true
    stationary mean in powers of the system size. The next term is obtained by
    carrying the second moment through the drift: writing A for the Jacobian at
    x* and B for its Hessian tensor, the correction V solves

        A*V + (1/2) * sum_{j,k} Sigma_{jk} * d2F/dx_j dx_k = 0

    with Sigma the stationary covariance from :func:`fluid_lyapunov`. Because
    Sigma scales with the population, V is the O(1/N) term written directly in
    job counts, so no explicit density rescaling is needed. The Hessian
    contraction is evaluated without ever forming the tensor: with
    ``Sigma = sum_m lam_m v_m v_m'`` by eigendecomposition,

        sum_{jk} Sigma_{jk} d2F/dx_j dx_k = sum_m lam_m d2F/dv_m^2

    and each directional second derivative is one central second difference, so
    the cost is O(rank(Sigma)) drift evaluations rather than O(n^2).

    The drift must be twice differentiable for this to mean anything. The
    first-order closure is only piecewise linear -- its second derivative is zero
    away from the kink and a delta at it -- so this must be called on the
    Gaussian-closed drift. A zero ``sigma2`` is REJECTED rather than silently
    returning zero.

    Port of ``fluid_refine_meanfield.m``.

    Args:
        x: the mean-field fixed point
        sigma2: (M,) station population variances defining the smooth drift
        Sigma: (n,n) stationary covariance
        drift_fcn: F(z, sigma2, covblk) -> (n,) drift
        jac_fcn: A(z, sigma2, covblk) -> (n,n) Jacobian
        D: the event incidence matrix, whose range is the reachable subspace
        covblk: per-station covariance blocks closing the DPS share ratio
        epsrel: relative step of the second difference

    Returns:
        (V, info) with V the correction to add to x.
    """
    sigma2 = np.asarray(sigma2, dtype=float).ravel()
    if not np.any(sigma2 > 0):
        raise ValueError(
            "The refined mean field expansion needs a twice-differentiable drift, but the "
            "first-order closure is only piecewise linear. Reach this through "
            "options.method='refined', which converges the Gaussian closure first.")

    x = np.asarray(x, dtype=float).ravel()
    n = x.size

    def F(z):
        return np.asarray(drift_fcn(z, sigma2, covblk), dtype=float).ravel()

    # eigendecomposition of the covariance, dropping numerically null directions
    Sigma = (np.asarray(Sigma, dtype=float) + np.asarray(Sigma, dtype=float).T) / 2
    lam, Vec = np.linalg.eigh(Sigma)
    if lam.size == 0 or float(np.max(lam)) <= 0:
        keep = np.zeros(lam.shape, dtype=bool)
    else:
        keep = (lam > float(np.max(lam)) * np.sqrt(np.finfo(float).eps)) & (lam > 0)
    Vec = Vec[:, keep]
    lam = lam[keep]

    scale = max(1.0, float(np.linalg.norm(x)))
    step = epsrel * scale
    b = np.zeros(n)
    F0 = F(x)
    for m in range(lam.size):
        d = Vec[:, m]
        b += lam[m] * (F(x + step * d) - 2.0 * F0 + F(x - step * d)) / step ** 2
    b *= 0.5

    # solve A*V = -b on the reachable subspace, where A is invertible
    A = np.asarray(jac_fcn(x, sigma2, covblk), dtype=float)
    Vbasis = orth(np.asarray(D, dtype=float))
    Ar = Vbasis.T @ A @ Vbasis
    cond_ar = float(np.linalg.cond(Ar)) if Ar.size else np.inf
    if not np.isfinite(cond_ar) or cond_ar > 1.0 / np.sqrt(np.finfo(float).eps):
        raise FluidNonHyperbolicError(
            "The fluid Jacobian is numerically singular on the reachable subspace (condition "
            "number %.3g), so the refinement equation A*V = -b has no meaningful solution. The "
            "fixed point sits at a drift kink or the model is marginally stable; use "
            "options.method='minnormal', which resums the same correction without inverting A."
            % cond_ar)
    V = Vbasis @ np.linalg.solve(Ar, -(Vbasis.T @ b))

    # the refinement is the next term of an asymptotic expansion, so it is only
    # meaningful while it stays small against the leading term; a correction of
    # the same size as the fixed point means the expansion has not kicked in at
    # this population, and returning it would be worse than refusing
    nv = float(np.linalg.norm(V))
    nx = float(np.linalg.norm(x))
    if nv > 0.5 * max(nx, np.sqrt(np.finfo(float).eps)):
        raise ValueError(
            "The 1/N refinement (norm %.3g) is not small against the mean-field fixed point "
            "(norm %.3g), so the asymptotic expansion is outside its range of validity at this "
            "population. Use options.method='minnormal'." % (nv, nx))

    info = {'rank': int(lam.size), 'stepsize': step,
            'residual': float(np.linalg.norm(A @ V + b)), 'condition': cond_ar}
    return V, info



def _blend_cov(cov_a, cov_b, step):
    """``cov_a + step*(cov_b - cov_a)`` entry by entry.

    A ``None`` entry is read as the zero matrix, and an entry stays ``None`` when
    both sides are. The matrix half of the damped variance step; the twin of the
    scalar ``sigma2ok + step*(sigma2 - sigma2ok)``, so the two stay consistent.
    """
    out = [None] * len(cov_b)
    for i in range(len(cov_b)):
        a = cov_a[i] if cov_a is not None and i < len(cov_a) else None
        b = cov_b[i]
        if a is None and b is None:
            out[i] = None
        elif a is None:
            out[i] = step * b
        elif b is None:
            out[i] = (1.0 - step) * a
        else:
            out[i] = a + step * (b - a)
    return out


def _kink_stations(sigma2, sched, nservers, lldscaling, stationBlock, x, M):
    """Every station whose population sits ON the saturation kink n_i = c_i of
    the first-order rate factor, in increasing order.

    Only the branches that take the indicator derivative can sit on one: a
    positive sigma2 or a load-dependent row makes the closure smooth, and an
    infinite server never saturates. Twin of
    ``FluidRateFactors.driftKinkStations`` in the JAR."""
    if sigma2 is not None and np.any(np.asarray(sigma2) > 0):
        return []  # the Gaussian closure is smooth, it has no kink
    tol = np.sqrt(np.finfo(float).eps)
    out = []
    for i in range(M):
        if sched[i] in (SchedStrategy.INF, SchedStrategy.EXT):
            continue
        if lldscaling is not None and i < np.shape(lldscaling)[0] \
                and np.any(np.asarray(lldscaling)[i, :] != 1.0):
            continue  # psi is piecewise quadratic and closed smoothly
        c = float(nservers[i])
        if not np.isfinite(c) or c <= 0:
            continue
        blk = stationBlock[i]
        if not np.size(blk):
            continue
        ni = float(np.sum(x[blk]))
        if ni <= 0:
            continue  # g = x on an empty station, no saturation term
        if abs(ni - c) <= tol * max(1.0, c):
            out.append(i)
    return out


def _nudge_off_kink(x, kink, nservers, stationBlock, rel):
    """A copy of ``x`` with every station in ``kink`` moved to ``c_i*(1+rel)``,
    i.e. strictly onto one side of it. The station's coordinates are scaled
    together, so the phase mix and every other station are untouched."""
    y = np.array(x, dtype=float, copy=True)
    for i in kink:
        blk = stationBlock[i]
        ni = float(np.sum(y[blk]))
        if not (ni > 0):
            continue
        y[blk] *= float(nservers[i]) * (1.0 + rel) / ni
    return y


class MomentTerms:
    """Everything the Gaussian closure is defined by, for one model.

    A plain attribute bag rather than a dataclass: it is built once, read
    everywhere, and mirrors the ``terms`` struct the MATLAB reference passes
    around field for field, so the two can be compared by name when a parity
    failure has to be traced.
    """

    __slots__ = (
        "M", "K", "Mu", "Phi", "phases", "rt",
        "nservers", "sched", "schedparam", "q_indices", "Kic", "enabled",
        "w", "D", "rateBase", "eventIdx", "nstate", "lld",
        "evIsDeparture", "evStation", "evClass", "Emap", "immediateAbsorb",
        "covIdx", "stationBlock", "classBlock",
        "shareSched", "minExact", "factors", "rates", "jac",
    )

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)

    def drift(self, x, s2, cb):
        """D r(x, sigma2): the closure's vector field."""
        return self.D @ self.rates(x, s2, cb)


class MinNormalSolver(ClosingMethod):
    """Gaussian moment-closure fluid solver, backing 'minnormal' and 'refined'."""

    def _is_refined(self) -> bool:
        """True when the caller asked for the O(1/N) refined mean field."""
        return str(getattr(self.options, 'method', '')).lower() in (
            'refined', 'fluid.refined')

    def build_moment_terms(self):
        """
        The closure's event representation, gathered once.

        The twin of MATLAB's ``fluid_moment_terms.m`` and the JAR's
        ``FluidMomentTerms.java``, which are standalone there and were inline here.
        It is factored out because ``dae.py`` needs exactly the same drift, rate
        factors and Jacobian as ``minnormal`` -- the DAE method IS the min-normal
        closure, differing only in how the coupled equations are discharged -- and
        a second copy of this would be a second closure, free to drift from the
        first the next time either is touched.

        THE STATE-SIZE CAP IS NOT APPLIED HERE. ``minnormal`` reads
        ``moment_maxstate`` (200) and the DAE route its own ``dae_maxstate`` (100),
        which is lower because its finite-difference Jacobian is over the unknowns
        rather than over the states, so the cost is quartic and not cubic.
        """
        M = self.sn.nstations
        K = self.sn.nclasses
        Mu, Phi, phases = self._extract_service_params()
        rt = self._get_routing_matrix()
        nservers = self._get_nservers()
        sched = self._get_sched()
        schedparam = self._get_schedparam()
        q_indices, Kic, enabled, w = self._build_ode_indices(M, K, Mu, phases, sched, schedparam)

        # GPS carries its weights in schedparam exactly as DPS does, and
        # _build_ode_indices now sets both; patching only here left the mean
        # solve, which builds its own w, integrating at unit GPS weights.
        D, rateBase, eventIdx = self._build_ode_system(
            M, K, Mu, Phi, phases, rt, enabled, q_indices, Kic)
        nstate = int(np.sum(Kic))

        # THE MOMENT CLOSURE READS THE SAME REDUCED EVENT SET AS EVERY OTHER ROUTE.
        # It used to refuse the reduction, on the grounds that it needs the
        # untransformed event set; what it actually needs is to be able to say which
        # (station,class) each event is a completion of, and Emap carries exactly
        # that across the composition -- an event folded through an immediate
        # coordinate keeps a row with weight on every original event it stands for,
        # including the two completions a pass-through realises at once. The
        # diffusion D*diag(r)*D' is then the diffusion of the reduced process, which
        # is the right one: the eliminated coordinate holds O(1/InfRate) mass and
        # contributes noise of the same order.
        from ..immediate import fluid_hide_immediate, ode_eliminate_immediate
        eventIdx0 = eventIdx
        D0 = D
        rateBase0 = rateBase
        Emap = np.eye(len(rateBase))
        immediate_absorb = None
        if fluid_hide_immediate(self.sn, self.options):
            D, rateBase, eventIdx, _, Emap, immediate_absorb = ode_eliminate_immediate(
                D, rateBase, eventIdx, self.sn, self.options)

        lld = getattr(self.sn, 'lldscaling', None)
        if lld is not None and np.size(lld) == 0:
            lld = None

        # the covariance is a dense nstate-by-nstate object and the Lyapunov
        # -- event classification, so throughputs can be read off the rates ----
        # Classified on the ORIGINAL events, which is what Emap maps onto. Without a
        # reduction Emap is the identity and the two indexings coincide.
        evIsDeparture, evStation, evClass = self._classify_events(
            M, K, enabled, q_indices, Kic, rt, D0, rateBase0, eventIdx0)

        # -- open and mixed models: the covariance lives on the QUEUE
        # coordinates only. The closing form models a Source as an EXT
        # pseudo-station holding unit mass, so its coordinate is a
        # normalisation constant, not a job count, and D*diag(r)*D' over it
        # would invent noise for a direction with no population. Dropping those
        # rows leaves exactly the right open event set, because with a
        # SINGLE-PHASE source the EXT rate factor is 1 identically: an arrival
        # is a CONSTANT-rate event whose jump, once the source row is dropped,
        # is a lone +1 into the destination queue, the canonical exogenous
        # Poisson arrival. A MULTI-PHASE source is refused: those coordinates
        # track the phase of one arrival process, not a population.
        covMask = np.ones(nstate, dtype=bool)
        for i in range(M):
            if sched[i] != SchedStrategy.EXT:
                continue
            for k in range(K):
                if Kic[i, k] > 0:
                    if Kic[i, k] > 1:
                        raise ValueError(
                            "The moment-closure method needs a Poisson arrival stream, but the "
                            "source of class %d is a %d-phase process. Those coordinates track "
                            "the phase of a single arrival process rather than a population, so "
                            "they carry no linear noise approximation. Use an exponential "
                            "inter-arrival time, or options.method='matrix'." % (k, Kic[i, k]))
                    covMask[q_indices[i, k]:q_indices[i, k] + Kic[i, k]] = False
        covIdx = np.where(covMask)[0]

        stationBlock = [None] * M
        classBlock = [[None] * K for _ in range(M)]
        for i in range(M):
            lo = q_indices[i, 0]
            hi = q_indices[i, K - 1] + Kic[i, K - 1] if K > 0 else lo
            stationBlock[i] = np.arange(lo, hi)
            for k in range(K):
                classBlock[i][k] = np.arange(q_indices[i, k], q_indices[i, k] + Kic[i, k])

        shareSched = np.array([sched[i] in (SchedStrategy.PS, SchedStrategy.FCFS,
                                            SchedStrategy.DPS, SchedStrategy.GPS)
                               for i in range(M)])

        # A STATION THAT CANNOT FILL ITS SERVERS HAS NOTHING TO CLOSE.
        # min(n_i,c_i) is the identity on the whole support whenever the
        # occupancy of station i is bounded above by its server count, and there
        # the Gaussian closure is not an improvement on the first-order one, it
        # is an ERROR: it spreads a normal marginal over n_i > c_i, mass the
        # station can never hold, and returns E[min(n_i,c_i)] < n_i. On a closed
        # model with one job per chain the exact answer is R = D at every queue
        # (a job cannot queue behind itself), which the first-order closure
        # reproduces to machine precision while the closure reads 0.4758 against
        # 0.5 on the queue length. The bound is the total population of every
        # chain that VISITS the station -- a station may declare a service time
        # for every class while the routing never sends most of them there --
        # and an open chain contributes an infinite population and never
        # qualifies. Their drift variance is held at zero below, exactly as at
        # the delay stations, whose min() is likewise absent.
        minExact = np.zeros(M, dtype=bool)
        chains = getattr(self.sn, 'chains', None)
        njobs = getattr(self.sn, 'njobs', None)
        visits = getattr(self.sn, 'visits', None)
        # sn.visits is indexed by STATEFUL node, not by station
        st2sf = getattr(self.sn, 'stationToStateful', None)
        st2sf = np.asarray(st2sf).ravel().astype(int) if st2sf is not None else None
        if chains is not None and np.size(chains) > 0 and njobs is not None:
            chains = np.asarray(chains)
            njobs = np.asarray(njobs, dtype=float).ravel()
            for i in range(M):
                if sched[i] in (SchedStrategy.EXT, SchedStrategy.INF) \
                        or not np.isfinite(nservers[i]):
                    continue
                bound = 0.0
                covered = np.zeros(K, dtype=bool)
                for ch in range(chains.shape[0]):
                    inch = np.where(chains[ch, :] > 0)[0]
                    if inch.size == 0:
                        continue
                    covered[inch] = True
                    vis = None
                    if visits is not None and len(visits) > ch:
                        vis = np.asarray(visits[ch])
                    if vis is not None and vis.size > 0 and st2sf is not None:
                        here = bool(np.any(vis[st2sf[i], inch] > 0))
                    else:
                        here = bool(np.any(enabled[i, inch]))
                    if here:
                        bound += float(np.sum(njobs[inch]))
                if np.any(np.asarray(enabled[i, :], dtype=bool) & ~covered):
                    continue  # a class outside every chain carries no bound
                minExact[i] = np.isfinite(bound) and bound <= nservers[i] + 1e-8

        def factors(x, s2, cb):
            return self._ode_rate_factors(x, M, K, enabled, q_indices, Kic, nservers,
                                          w, sched, sigma2=s2, lld=lld, covblk=cb)

        def rates(x, s2, cb):
            return factors(x, s2, cb)[eventIdx] * rateBase

        def jac(x, s2, cb):
            return self._drift_jacobian(x, M, K, enabled, q_indices, Kic, nservers,
                                        w, sched, rateBase, eventIdx, D, s2, lld, cb)

        return MomentTerms(
            M=M,
            K=K,
            Mu=Mu,
            Phi=Phi,
            phases=phases,
            rt=rt,
            nservers=nservers,
            sched=sched,
            schedparam=schedparam,
            q_indices=q_indices,
            Kic=Kic,
            enabled=enabled,
            w=w,
            D=D,
            rateBase=rateBase,
            eventIdx=eventIdx,
            nstate=nstate,
            lld=lld,
            evIsDeparture=evIsDeparture,
            evStation=evStation,
            evClass=evClass,
            Emap=Emap,
            immediateAbsorb=immediate_absorb,
            covIdx=covIdx,
            stationBlock=stationBlock,
            classBlock=classBlock,
            shareSched=shareSched,
            minExact=minExact,
            factors=factors,
            rates=rates,
            jac=jac,
        )

    def solve(self) -> FLDResult:
        start_time = time.time()
        # fail before integration rather than from inside the ODE right-hand
        # side, where the LSODA callback reports it as a solver failure
        self._check_scheds()

        terms = self.build_moment_terms()
        M = terms.M
        K = terms.K
        Mu = terms.Mu
        Phi = terms.Phi
        phases = terms.phases
        rt = terms.rt
        nservers = terms.nservers
        sched = terms.sched
        schedparam = terms.schedparam
        q_indices = terms.q_indices
        Kic = terms.Kic
        enabled = terms.enabled
        w = terms.w
        D = terms.D
        rateBase = terms.rateBase
        eventIdx = terms.eventIdx
        nstate = terms.nstate
        lld = terms.lld
        evIsDeparture = terms.evIsDeparture
        evStation = terms.evStation
        evClass = terms.evClass
        covIdx = terms.covIdx
        stationBlock = terms.stationBlock
        classBlock = terms.classBlock
        shareSched = terms.shareSched
        minExact = terms.minExact
        factors = terms.factors
        rates = terms.rates
        jac = terms.jac

        # the covariance is a dense nstate-by-nstate object and the Lyapunov
        # solve is cubic in it, so refuse rather than silently crawl
        maxstate = getattr(self.options, 'moment_maxstate', None) or 200
        if nstate > maxstate:
            raise ValueError(
                "The moment-closure method solves a %dx%d Lyapunov equation, above the limit "
                "of %d. Raise that limit or use options.method='closing'." % (nstate, nstate, maxstate))

        # -- outer fixed point ------------------------------------------------
        x0 = self._compute_initial_state(M, K, phases)
        sigma2 = np.zeros(M)
        covblk = [None] * M
        outer_max = 20
        # THE CLOSURE IS JUDGED FAR TIGHTER THAN 1e-3, so it must not stop there.
        # Converged only to 1e-3 this alternation is not a fixed point to two
        # machines: on mqn_singleserver_ps the MATLAB twin answered 42.3962 on two
        # hosts and 42.8207 on a third, 1e-2 relative apart, because the transient
        # iterate below fell on opposite sides. 1e-6 is the loosest that
        # reproduces; min(), not assignment, so a caller may still ask tighter.
        #
        # It governs the INNER mean solve too, through self.options.iter_tol below.
        # That is not incidental: the transient non-hyperbolic iterate this method
        # used to abort on was an artefact of a loosely converged inner solve.
        mom_tol = 1e-6
        _opt_iter_tol = getattr(self.options, 'iter_tol', None)
        if _opt_iter_tol:
            mom_tol = min(mom_tol, float(_opt_iter_tol))
        outer_tol = mom_tol
        # The last variance whose Lyapunov solve SUCCEEDED, and the floor on the
        # step taken toward the next one. See the damping in the loop below.
        sigma2ok = np.zeros(M)
        covok = [None] * M
        damp_min = 1.0 / 64.0
        Sigma = np.zeros((nstate, nstate))
        sigma2solve = sigma2
        covsolve = covblk
        xvec_t = None
        t = None
        x = x0

        from line_solver.api.io import console as _console
        _console.loop('iterating the moment closure (at most %d passes)', outer_max)
        # A DELAY STATION HAS NO min() TO CLOSE, so its variance must never reach
        # the drift -- only the report. This mask used to be applied to
        # sigma2drift after the loop and nowhere inside it, so every mean solve
        # of the fixed point ran with the delay variance switched on. The rate
        # factor there is mu*n, which the Gaussian correction turns into
        # something that does not vanish with n: the coordinate is driven
        # NEGATIVE, the drift is conservative so another coordinate grows to
        # match, and the trajectory leaves the simplex for good. On CQN_Cox_CS_9
        # (Delay + PS + PS(c=5), N=6) the first window past sigma2 = 0 moved
        # 8.7e3 of mass and the drift norm reached 5.4e9.
        no_drift_var = [sched[i] in (SchedStrategy.INF, SchedStrategy.EXT)
                        for i in range(M)]
        for outer in range(outer_max):
            _console.iter_line(outer + 1, 'closure pass %d', outer + 1)
            # A TRANSIENT ITERATE MUST NOT VETO THE METHOD. The Lyapunov gate asks
            # whether the linear noise approximation has a stationary covariance at
            # the point THIS iterate landed on; a fixed point that fails it is a
            # model the closure cannot answer, but an intermediate iterate that
            # fails it is only a variance step that overshot. On mqn_singleserver_ps
            # iterate 1 was stable at -4.93e-03, iterate 2 declined at +1.07e+01,
            # and the fixed point the fallback then found was stable at -4.92e-03.
            # So a failing iterate RETREATS toward the last variance that succeeded,
            # halving until the LNA is defined again; only a step below damp_min, or
            # a failure at the seed where there is nothing to retreat toward, is the
            # model's own non-hyperbolicity and still raises.
            step = 1.0
            while True:
                sigma2try = sigma2ok + step * (sigma2 - sigma2ok)
                covtry = _blend_cov(covok, covblk, step)
                sigma2solve = sigma2try.copy()
                covsolve = list(covtry)
                # the closure stays first order where min(n,c) is the identity; the
                # share closure follows, since mu_r*(n_r/n)*min(n,c) collapses to
                # mu_r*n_r there. The covariance is still solved and still reported,
                # it just does not enter the drift.
                for i in range(M):
                    if minExact[i] or no_drift_var[i]:
                        sigma2solve[i] = 0.0
                        covsolve[i] = None
                x, xvec_t, t = self._integrate(M, K, Mu, Phi, phases, rt, nservers,
                                               sched, schedparam, x0,
                                               sigma2solve, covsolve,
                                               iter_tol=mom_tol)

                r = rates(x, sigma2solve, covsolve)

                # A POINT ON A SATURATION KINK HAS NO JACOBIAN. The rate factor
                # min(n_i,c_i) has slope 1 below c_i and 0 above, and ``jac`` resolves
                # the tie onto the saturated side, so a verdict read off it would
                # depend on which side the integrator stopped. The VERDICT, not the
                # point, has to be side-independent: both one-sided Jacobians are
                # ordinary matrices, so ASK BOTH and decline only when a side fails.
                # Refusing at every kink instead throws away models the reference
                # solves -- the first outer iterate runs at sigma2 = 0 and a saturated
                # model's first-order fixed point lands on the kink by construction.
                # Later iterates carry a positive sigma2 and are smooth, so this costs
                # two Jacobians on the seed and nothing after it. Twin of
                # FluidRateFactors.driftKinkStation/nudgedOffKink in the JAR.
                kink = _kink_stations(sigma2solve, sched, nservers, lld,
                                      stationBlock, x, M)
                if kink:
                    for rel in (-1e-6, 1e-6):
                        xs = _nudge_off_kink(x, kink, nservers, stationBlock, rel)
                        try:
                            Ds = D[covIdx, :]
                            rs = rates(xs, sigma2solve, covsolve)
                            fluid_lyapunov(jac(xs, sigma2solve, covsolve)[np.ix_(covIdx, covIdx)],
                                           Ds @ np.diag(rs) @ Ds.T, Ds)
                        except FluidNonHyperbolicError as err:
                            raise FluidNonHyperbolicError(
                                'The fluid fixed point sits on the saturation kink of station %d '
                                '(population equals its %g servers) and the two one-sided drift '
                                'Jacobians there disagree on hyperbolicity, so which of them the '
                                'linear noise approximation would use is decided by the '
                                "integrator's rounding residue rather than by the model. This is "
                                'the saturated boundary of a continuum of equilibria; use '
                                "method='closing' for the mean only. Underlying: %s"
                                % (kink[0] + 1, float(nservers[kink[0]]), err))

                A = jac(x, sigma2solve, covsolve)
                Dc = D[covIdx, :]
                Qc = Dc @ np.diag(r) @ Dc.T
                try:
                    Sc = fluid_lyapunov(A[np.ix_(covIdx, covIdx)], Qc, Dc)
                    break
                except FluidNonHyperbolicError:
                    at_seed = (float(np.linalg.norm(sigma2 - sigma2ok, 1)) == 0.0
                               and all(c is None for c in covblk))
                    if step <= damp_min or at_seed:
                        raise
                    step = step / 2.0
            sigma2ok = sigma2try
            covok = covtry
            Sigma = np.zeros((nstate, nstate))
            Sigma[np.ix_(covIdx, covIdx)] = Sc

            sigma2new = np.zeros(M)
            covnew = [None] * M
            for i in range(M):
                blk = stationBlock[i]
                if blk.size:
                    sigma2new[i] = max(0.0, float(np.sum(Sigma[np.ix_(blk, blk)])))
                    if shareSched[i]:
                        covnew[i] = Sigma[np.ix_(blk, blk)]

            delta = np.linalg.norm(sigma2new - sigma2try, 1) / max(1.0, np.linalg.norm(sigma2new, 1))
            for i in range(M):
                # sigma2 is the SUM of a block, so it can converge while the
                # off-diagonals the share closure reads are still moving
                if covnew[i] is not None:
                    dc = covnew[i] if covtry[i] is None else covnew[i] - covtry[i]
                    delta = max(delta, np.linalg.norm(dc, 1) / max(1.0, np.linalg.norm(covnew[i], 1)))
            sigma2 = sigma2new
            covblk = covnew
            if delta < outer_tol:
                break

        # Metrics must be read at the SAME variance the mean solve used, not at
        # the variance that solve produced: the latter evaluates the rate
        # functions at a point that is not their fixed point and throughput then
        # fails to balance. The delay stations have no min() to close, so their
        # variance must not enter the drift; keep it for reporting only.
        # no_drift_var already masked sigma2solve on the way in, so this IS the
        # variance the mean solve used, which is what the paragraph above asks for.
        sigma2drift = sigma2solve.copy()
        covdrift = list(covsolve)

        refinement = None
        if self._is_refined():
            # The refinement expands about the MEAN-FIELD fixed point, not the
            # Gaussian one: adding it to the minnormal point would count the same
            # O(1/N) term twice, since the Gaussian closure already resums it.
            # So the base point is recomputed with the first-order closure, while
            # the Hessian and the Jacobian read the SMOOTH Gaussian drift.
            sigma2mf = np.zeros(M)
            covmf = [None] * M
            xmf, xvec_t, t = self._integrate(M, K, Mu, Phi, phases, rt, nservers,
                                             sched, schedparam, x0, sigma2mf, covmf)
            A = jac(xmf, sigma2drift, covdrift)
            rmf = rates(xmf, sigma2drift, covdrift)
            Dc = D[covIdx, :]
            Sc = fluid_lyapunov(A[np.ix_(covIdx, covIdx)],
                                Dc @ np.diag(rmf) @ Dc.T, Dc)
            Sigma = np.zeros((nstate, nstate))
            Sigma[np.ix_(covIdx, covIdx)] = Sc

            def drift(z, s2, cb):
                return D @ rates(z, s2, cb)

            # A LINEAR DRIFT NEEDS NO REFINEMENT, and that is not the same as
            # the degenerate call fluid_refine_meanfield refuses. When EVERY
            # station is either an infinite server or minExact -- min(n,c) is
            # the identity on the reachable set, because the population bound
            # never reaches c -- the drift is exactly affine there, its Hessian
            # vanishes and the O(1/N) correction is identically zero. The mask
            # above then zeroes all of sigma2drift, which the refinement reads
            # as "the caller handed me the first-order closure" and rejects. So
            # settle it here, where the reason for the zero is known: a null
            # correction, not an error. (Delay + PS(c=2) at N=2 is the smallest
            # case; without this, method='refined' raised on a model whose
            # mean field is already exact.)
            if all(bool(minExact[i]) or bool(no_drift_var[i]) for i in range(M)):
                refinement = np.zeros_like(xmf)
                _info = None
            else:
                refinement, _info = fluid_refine_meanfield(
                    xmf, sigma2drift, Sigma, drift, jac, D, covdrift)
            x = xmf + refinement
            x[x < 0] = 0.0
            # the corrected point is a correction OF the mean-field fixed point,
            # so its rates are read with the mean-field (zero) variance
            sigma2drift = sigma2mf
            covdrift = [None] * M
            sigma2 = np.zeros(M)
            for i in range(M):
                blk = stationBlock[i]
                if blk.size:
                    sigma2[i] = max(0.0, float(np.sum(Sigma[np.ix_(blk, blk)])))

        r = rates(x, sigma2drift, covdrift)
        gfac = factors(x, sigma2drift, covdrift)

        # a load-dependent station clears alpha(n) times the nominal work, so
        # its utilization normalises by the peak scaling (T*S/peak, as in CTMC)
        Seff = np.array(nservers, dtype=float).copy()
        N = float(np.sum(getattr(self.sn, 'njobs', np.zeros(1))))
        for i in range(M):
            if not np.isfinite(Seff[i]):
                Seff[i] = max(N, 1.0)
        if lld is not None:
            for i in range(min(M, len(lld))):
                row = np.asarray(lld[i], dtype=float).ravel()
                if row.size:
                    Seff[i] = max(Seff[i], float(np.max(row)))

        QN = np.zeros((M, K)); UN = np.zeros((M, K))
        TN = np.zeros((M, K)); RN = np.zeros((M, K))
        for i in range(M):
            for k in range(K):
                blk = classBlock[i][k]
                if blk.size == 0:
                    continue
                QN[i, k] = float(np.sum(x[blk]))
                if sched[i] == SchedStrategy.INF:
                    UN[i, k] = QN[i, k]
                else:
                    UN[i, k] = float(np.sum(gfac[blk])) / Seff[i]
                # Summed over ORIGINAL events through Emap: a reduced event folded
                # through an immediate coordinate is a completion at more than one
                # (station,class), and its rate has to reach every one of them.
                sel = evIsDeparture & (evStation == i) & (evClass == k)
                TN[i, k] = float(np.dot(r, terms.Emap @ sel.astype(float)))
            if sched[i] == SchedStrategy.EXT:
                # A Source holds no jobs, exactly as `_compute_metrics_closing`
                # rules: its coordinates are the arrival process's phase
                # indicator, and the station-level routing folds every departure
                # of an open class back onto them, so the coordinate accumulates
                # the mass that LEFT the system. Reading that as a queue length
                # also poisons RN = QN/TN below.
                QN[i, :] = 0.0
        # A CLASS THE MODEL NEVER ROUTES HERE HAS NO RESPONSE TIME, and neither
        # QN nor TN says so by its size: both hold a remnant of the initial state
        # the integrator was still draining. See fluid_visited_pairs.
        nz = fluid_visited_pairs(self.sn, M, K) & (TN > GlobalConstants.Zero)
        RN[nz] = QN[nz] / TN[nz]

        QVar = np.zeros((M, K))
        for i in range(M):
            for k in range(K):
                blk = classBlock[i][k]
                if blk.size:
                    QVar[i, k] = max(0.0, float(np.sum(Sigma[np.ix_(blk, blk)])))

        CN = np.sum(RN, axis=0, keepdims=True)
        XN = self._compute_system_throughput(TN, M, K)
        from line_solver.api.sn.getters import sn_get_arvr_from_tput
        from line_solver.api.sn.transforms import sn_get_residt_from_respt
        AN = sn_get_arvr_from_tput(self.sn, TN)
        WN = sn_get_residt_from_respt(self.sn, RN, None)

        result = FLDResult(
            QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN, AN=AN, WN=WN,
            t=t if t is not None else np.array([0.0]),
            xvec=x, iterations=outer + 1,
            runtime=time.time() - start_time,
            method='refined' if self._is_refined() else 'minnormal')
        result.moments = {
            'Sigma': Sigma, 'QVar': QVar, 'QStd': np.sqrt(QVar),
            'sigma2': sigma2, 'outerIters': outer + 1,
            # the same variance as it entered the DRIFT (zero at the delay
            # stations): a later solve on this fixed point, i.e. the passage-time
            # ODE, must close its capacity term at this one or it evaluates a
            # first-order drift at a second-order fixed point
            'sigma2Drift': sigma2drift,
            'stationBlock': stationBlock, 'classBlock': classBlock,
            # None for 'minnormal', which resums the same term instead of
            # adding it; the (n,) vector for 'refined'
            'refinement': refinement}
        return result

    # ------------------------------------------------------------------ utils
    def _integrate(self, M, K, Mu, Phi, phases, rt, nservers, sched, schedparam,
                   x0, sigma2, covblk, iter_tol=None):
        """Solve the mean under the given closure, through the closing method.

        MATLAB `solver_fluid_moments.m` (meanopt.method = 'closing' with the
        closure in config.moment_sigma2/moment_cov), the JAR `MinNormalAnalyzer`
        and the C++ `solver_fluid_moments` all run the mean solve through the
        ordinary closing path. Doing the same here inherits its horizon, its
        `options.stiff`, `timespan` and `odemaxstep` handling and its max-step
        bound, none of which a private LSODA call at a fixed tmax honoured: on an
        LN layer carrying Immediate (1e8) rates that call failed its very first
        corrector step and returned x0 as the fixed point.
        """
        saved = getattr(self.options, 'config', None)
        saved_iter_tol = getattr(self.options, 'iter_tol', None)
        cfg = dict(saved) if saved else {}
        cfg['moment_sigma2'] = np.asarray(sigma2, dtype=float)
        cfg['moment_cov'] = list(covblk)
        self.options.config = cfg
        # The closure's own convergence tolerance governs THIS solve too; see
        # mom_tol at the outer fixed point. Restored with config, because
        # self.options is shared with everything else the solver does.
        if iter_tol is not None:
            self.options.iter_tol = iter_tol
        try:
            _, xvec_t, t = self._solve_fluid_ode(
                M, K, Mu, Phi, phases, rt, nservers, sched, schedparam, x0)
        finally:
            self.options.config = saved
            if iter_tol is not None:
                self.options.iter_tol = saved_iter_tol

        if self.ode_status < 0:
            raise FluidIntegrationError(
                "minnormal: the mean-field ODE failed to integrate (%s). The state it "
                "returned is the initial guess, not a fixed point, so no closure can be "
                "read off it." % self.ode_message)

        xt = np.asarray(xvec_t, dtype=float)
        x = np.maximum(xt[-1, :], 0.0)
        return x, xt, t

    def _classify_events(self, M, K, enabled, q_indices, Kic, rt, D, rateBase, eventIdx):
        """Which events are service completions, and where they are sourced.

        ``_build_ode_system`` emits every departure first and then every
        intra-PH phase change, so counting the departures classifies the vector.
        Summing the completion rates sourced at (i,c) gives the class-c
        throughput at station i exactly, because the routing probabilities and
        the PH entry vector each sum to one over the destinations enumerated.
        """
        nDeparture = 0
        for i in range(M):
            for k in range(K):
                if not enabled[i, k]:
                    continue
                for j in range(M):
                    for l in range(K):
                        if i * K + k >= rt.shape[0] or j * K + l >= rt.shape[1]:
                            continue
                        if rt[i * K + k, j * K + l] > 0 and enabled[j, l]:
                            nDeparture += Kic[i, k] * Kic[j, l]
        nevents = len(rateBase)
        evIsDeparture = np.zeros(nevents, dtype=bool)
        evIsDeparture[:min(nDeparture, nevents)] = True

        nstate = int(np.sum(Kic))
        coordStation = np.zeros(nstate, dtype=int)
        coordClass = np.zeros(nstate, dtype=int)
        for i in range(M):
            for k in range(K):
                idx = np.arange(q_indices[i, k], q_indices[i, k] + Kic[i, k])
                coordStation[idx] = i
                coordClass[idx] = k
        return evIsDeparture, coordStation[eventIdx], coordClass[eventIdx]

    def _drift_jacobian(self, x, M, K, enabled, q_indices, Kic, nservers, w, sched,
                        rateBase, eventIdx, D, sigma2=None, lld=None, covblk=None):
        """Analytic Jacobian of the fluid drift.

        Must mirror ``_ode_rate_factors`` branch by branch: any policy without a
        case there keeps g = x and contributes the identity here.
        """
        x = np.maximum(np.asarray(x, dtype=float).ravel(), 0.0)
        n = len(x)
        if sigma2 is None:
            sigma2 = np.zeros(M)
        gaussian = bool(np.any(np.asarray(sigma2) > 0))
        has_lld = lld is not None and len(lld) > 0
        G = np.eye(n)  # INF, EXT phases 2.., and every policy without a case

        for i in range(M):
            lldrow = None
            if has_lld and i < len(lld):
                row = np.asarray(lld[i], dtype=float).ravel()
                if row.size > 0 and np.any(np.abs(row - 1.0) > 1e-14):
                    lldrow = row
            Ci = covblk[i] if (covblk is not None and i < len(covblk)) else None

            lo = q_indices[i, 0]
            hi = q_indices[i, K - 1] + Kic[i, K - 1] if K > 0 else lo
            blk = np.arange(lo, hi)

            if sched[i] == SchedStrategy.INF:
                if lldrow is not None and blk.size:
                    ni = float(np.sum(x[blk]))
                    if ni > 0:
                        h, dh = capacity_closure(ni, nservers[i], sigma2[i], lldrow, True)
                        f = h / ni
                        fp = (ni * dh - h) / ni ** 2
                        G[np.ix_(blk, np.arange(n))] = 0
                        G[np.ix_(blk, blk)] = f * np.eye(blk.size) + np.outer(x[blk] * fp, np.ones(blk.size))

            elif sched[i] == SchedStrategy.EXT:
                for k in range(K):
                    if enabled[i, k]:
                        a0 = q_indices[i, k]
                        a1 = a0 + Kic[i, k]
                        G[a0, :] = 0
                        if a1 > a0 + 1:
                            G[a0, a0 + 1:a1] = -1

            elif sched[i] in (SchedStrategy.PS, SchedStrategy.FCFS):
                if blk.size == 0:
                    continue
                ni = float(np.sum(x[blk]))
                if ni <= 0:
                    continue
                if gaussian or lldrow is not None:
                    h, dh, d2h = capacity_closure(ni, nservers[i], sigma2[i], lldrow,
                                                  False, want_d2=True)
                    if Ci is not None:
                        # g = s(x_blk)*h(ni) + h'(ni)*cn(x_blk), the joint closure of
                        # the share and the capacity; differentiating it with the
                        # covariance held fixed adds h'*dcn and h''*cn to the
                        # product rule.
                        s, ds, cn, dcn = share_closure(x[blk], np.ones(blk.size), Ci,
                                                       want_jac=True, want_cov=True)
                        G[np.ix_(blk, np.arange(n))] = 0
                        G[np.ix_(blk, blk)] = (ds * h + np.outer(s, dh * np.ones(blk.size))
                                               + dh * dcn
                                               + np.outer(cn, d2h * np.ones(blk.size)))
                        continue
                elif ni > nservers[i] - GlobalConstants.FineTol * max(1.0, ni):
                    # THE SATURATION TEST CARRIES A BAND, and it is a cross-codebase
                    # requirement: a saturated fixed point sits exactly at ni = c, and
                    # each engine's ODE stops on its own residual (MATLAB 1.0004, the
                    # C++ port 1 - 1.8e-13 on the same model). A strict ni > c reads
                    # saturated in one and unsaturated in the other, which flips this
                    # whole station block between a zero row and the identity, and
                    # with it the hyperbolicity verdict and the method SolverFLD
                    # answers with. See ``min_closure``, which carries the same band.
                    h, dh = float(nservers[i]), 0.0
                else:
                    continue
                f = h / ni
                fp = (ni * dh - h) / ni ** 2
                G[np.ix_(blk, np.arange(n))] = 0
                G[np.ix_(blk, blk)] = f * np.eye(blk.size) + np.outer(x[blk] * fp, np.ones(blk.size))

            elif sched[i] == SchedStrategy.DPS:
                wi = w[i, :] / np.sum(w[i, :])
                wv = np.zeros(blk.size)
                for k in range(K):
                    if enabled[i, k]:
                        a0 = q_indices[i, k] - lo
                        wv[a0:a0 + Kic[i, k]] = wi[k]
                xi = float(np.sum(x[blk]))
                if xi <= 0 or float(wv @ x[blk]) <= 0:
                    continue
                psi, dpsi, d2psi = capacity_closure(xi, nservers[i], sigma2[i], lldrow,
                                                    False, want_d2=True)
                s, ds, cn, dcn = share_closure(x[blk], wv, Ci, want_jac=True,
                                               want_cov=True)
                G[np.ix_(blk, np.arange(n))] = 0
                G[np.ix_(blk, blk)] = (ds * psi + np.outer(s, dpsi * np.ones(blk.size))
                                       + dpsi * dcn
                                       + np.outer(cn, d2psi * np.ones(blk.size)))

            elif sched[i] == SchedStrategy.GPS:
                # g_j = (x_j/x_k)*s_k*a, so
                #   dg_j/dx_l = [delta_jl/x_k - x_j/x_k^2]*s_k*a   (l in class k)
                #             + (x_j/x_k)*ds_k/dx_m*a              (l in class m)
                #             + (x_j/x_k)*s_k*da/dxi               (l in station)
                # The weighted share of a GPS server has no rate law in LINE past
                # one server, so refuse here as the drift does in closing.py and
                # as fluid_drift_jacobian.m, FluidRateFactors and fluid_moments.h
                # do: the Jacobian is the second entry to this closure.
                if nservers[i] > 1:
                    raise ValueError('Multi-server GPS stations are not supported yet.')
                xk = np.zeros(K); vk = np.zeros(K); bk = [None] * K
                for k in range(K):
                    if enabled[i, k]:
                        a0 = q_indices[i, k]
                        bk[k] = np.arange(a0, a0 + Kic[i, k])
                        xk[k] = float(np.sum(x[bk[k]]))
                        if Ci is not None:
                            loc = np.arange(a0 - lo, a0 - lo + Kic[i, k])
                            vk[k] = max(0.0, float(np.sum(np.asarray(Ci)[np.ix_(loc, loc)])))
                sk, dsk = gps_share(xk, w[i, :], vk, want_jac=True)
                xi = float(np.sum(x[blk]))
                a, da = 1.0, 0.0
                if lldrow is not None:
                    aa, dda = lld_scaling(lldrow, xi)
                    a, da = float(aa[0]), float(dda[0])
                G[np.ix_(blk, np.arange(n))] = 0
                for k in range(K):
                    if not enabled[i, k] or xk[k] <= 0:
                        continue
                    b = bk[k]
                    G[np.ix_(b, b)] += (sk[k] * a / xk[k]) * np.eye(b.size) \
                        - (sk[k] * a / xk[k] ** 2) * np.outer(x[b], np.ones(b.size))
                    for m in range(K):
                        if enabled[i, m]:
                            G[np.ix_(b, bk[m])] += (a * dsk[k, m] / xk[k]) * np.outer(x[b], np.ones(bk[m].size))
                    if da != 0:
                        G[np.ix_(b, blk)] += (sk[k] * da / xk[k]) * np.outer(x[b], np.ones(blk.size))

        return D @ (rateBase[:, None] * G[eventIdx, :])


def solve_minnormal(sn, options: Optional[SolverFLDOptions] = None) -> FLDResult:
    """Entry point for ``options.method='minnormal'`` and ``'refined'``."""
    if options is None:
        options = SolverFLDOptions()
    return MinNormalSolver(sn, options).solve()
