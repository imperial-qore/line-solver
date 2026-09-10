"""
Differential-algebraic formulation of the min-normal closure: ``options.method='dae'``.

The port of MATLAB's ``solver_fluid_dae.m``, and the twin of the JAR's
``DaeAnalyzer`` and of C++'s ``fluid_dae.h``.

``minnormal`` already solves a differential system (the mean) coupled to an
algebraic one (the covariance). It solves them by SUCCESSIVE SUBSTITUTION:
integrate the mean to its fixed point at a held variance, solve the Lyapunov
equation there, extract sigma2, repeat, up to 20 times and only to a coarse
tolerance. This method states the same closure as one system and solves it as
one system. NOTHING ABOUT THE CLOSURE CHANGES -- the drift, the rate factors
and the Lyapunov equation are the same functions, taken unmodified from
``MomentTerms`` and ``fluid_lyapunov`` -- only the way the coupled equations are
discharged.

Two modes, chosen by the horizon, because the DAE is genuinely different in
each:

  STEADY STATE (``timespan[1]`` infinite, the usual case)
      Solve the algebraic system

          0 = D r(x, sigma2)              drift residual, nstate rows
          0 = C x - Nchain                population conservation, one row per
                                          closed chain
          0 = sigma2 - sigmaOf(x,sigma2)  closure consistency, one row per
                                          closable station

      simultaneously by a damped projected Newton. Convergence is quadratic
      near the root instead of the linear rate of substitution, and the answer
      is converged to ``options.tol`` rather than to the outer coarse tolerance.

      THE FIXED POINT IS NOT FOUND BY INTEGRATING TO IT. Integrating a stable
      ODE until it stops moving is a poor way to solve f(x)=0: the cost is set
      by the slowest mode of the model rather than by the accuracy wanted. One
      seed trajectory is still integrated, cheaply and at the first-order
      closure, because Newton needs a point inside the basin; everything after
      that is algebraic.

  TRANSIENT (finite horizon)
      Integrate the index-1 DAE

          d/dt x     = D r(x, sigma2(t))     differential rows
          0          = C x - Nchain          algebraic, one per closed chain
          d/dt Sigma = A Sigma + Sigma A' + Q  differential

      with a SINGULAR mass matrix. This and ``kp`` are the only fluid methods
      that produce a time-varying second moment, and unlike ``minnormal`` --
      which evaluates its whole transient at the single STATIONARY variance --
      the variance here is the one the trajectory actually had at each instant.

WHY RODAS AND NOT scipy. scipy has no mass matrix at all: ``solve_ivp`` solves
y' = f, and BDF/Radau take no M. An index-1 DAE with a singular M therefore has
no scipy route, so this integrates with the vendored Hairer-Wanner RODAS in
``..ode.rodas`` -- the same integrator the C++ and JAR ports use and the one
MATLAB reaches through ``options.odesolvers.daeSolver = @rodas``.

WHAT THE ALGEBRAIC CONSTRAINT BUYS. Population conservation otherwise holds
only to integrator tolerance: it is a CONSEQUENCE of the drift (the rows of D
sum to zero on a closed chain), never an equation. Writing it as a constraint
enforces it to solver tolerance, and it is also what makes the Newton system
solvable at all -- the drift Jacobian is singular along exactly the conserved
directions, the same singularity ``fluid_lyapunov`` works around by projecting
onto range(D), so the constraint rows supply the missing rank instead of a
pseudo-inverse hiding it.

WHY THE COVARIANCE IS NOT A NEWTON UNKNOWN. Sigma is nstate^2 entries, so a
Jacobian over it is quartic work -- strictly worse than the cubic Lyapunov
solves it would replace. Sigma is LINEAR in itself for a held x, so it is
eliminated by one Lyapunov solve per residual evaluation and only sigma2, M
numbers, joins x in the unknown vector.

NO sigma2=0 SEED, SO NO KINK WORKAROUND. ``minnormal`` must start its
alternation at sigma2 = 0, where min(n,c) has no derivative, and a saturated
model's first-order fixed point lands on that kink by construction; it carries a
two-sided-Jacobian probe to decide hyperbolicity side-independently there. The
simultaneous solve never adopts sigma2 = 0 as an iterate, so that probe is not
needed. It does NOT rescue a fixed point that genuinely sits on the kink at the
converged variance -- that is a real continuum of equilibria and no formulation
removes it.
"""

import time
import warnings
from typing import Optional

import numpy as np

from line_solver.api.sn import SchedStrategy

from ..options import FLDResult, SolverFLDOptions
from ..ode.rodas import rodas
from .minnormal import (FluidNonHyperbolicError, MinNormalSolver,
                        fluid_lyapunov)
from line_solver.constants import GlobalConstants
from ..utils.metrics import fluid_visited_pairs

__all__ = ["DaeSolver", "solve_dae", "capacity_constraints", "capacity_staging",
           "capacity_gates", "CapacityConstraints", "CapacityStaging",
           "CapacityGates"]

_ZERO = 1e-14
_FINE = 1e-8

# HOW A BLOCKED JOB IS HELD, which is not one answer but three, and the three
# need three different equations: held at the upstream station whose departure is
# disabled, lost (the arrival fires and is not admitted), or staged in an entry
# queue. A cap carries `staged`; the other two are per EVENT, as gates.held and
# gates.loss. See CAPACITY_GATES for which model each cap gets and why.
# FiniteCapacityRegion.UNBOUNDED, and DropStrategy.WaitingQueue, as sn stores
# them. Kept as literals rather than imported so that this module does not
# depend on the lang layer for two integers the interchange format pins anyway.
_UNBOUNDED = -1.0
_WAITQ = -1.0
_DROP = 1.0
# DropStrategy, for the refusal messages. Kept as literals for the same reason
# the two above are: the interchange format pins the integers anyway.
_DROP_NAMES = {-1.0: "a waiting queue", 1.0: "a drop", 2.0: "BAS blocking",
               3.0: "BBS blocking", 4.0: "RSRD blocking", 5.0: "retrial",
               6.0: "retrial with a limit"}


# ---------------------------------------------------------------------------
# capacity limits, as linear admission constraints
# ---------------------------------------------------------------------------

class CapacityConstraints(object):
    """``A x <= b``, one row per cap, with the limit each row came from."""

    __slots__ = ("A", "As", "b", "region", "station", "cls", "staged", "label",
                 "member", "klass", "nregions")

    def __init__(self, nstate):
        self.A = np.zeros((0, nstate))
        self.As = np.zeros((0, 0))   # the same rows over the staging coordinates
        self.b = np.zeros(0)
        self.region = np.zeros(0, dtype=int)     # -1 for a station cap
        self.station = np.zeros(0, dtype=int)    # -1 for a region cap
        self.cls = np.zeros(0, dtype=int)        # -1 when it limits several
        self.staged = np.zeros(0, dtype=bool)    # True: the entry-queue model
        self.label = []
        self.member = np.zeros((0, nstate), dtype=bool)
        self.klass = np.zeros(nstate, dtype=int)
        self.nregions = 0

    @property
    def empty(self):
        return self.b.size == 0


def _eliminated_coords(terms, nstate):
    """
    Mark the coordinates the immediate reduction folded away, where the reduced
    drift holds no mass and no event lands.

    ``terms.immediateAbsorb`` is the projector the reduction returns: the
    identity on a surviving coordinate and the absorption distribution on an
    eliminated one, so a zero diagonal is exactly the eliminated case. It is None
    when nothing was eliminated, where every coordinate survives.
    """
    gone = np.zeros(int(nstate), dtype=bool)
    absorb = getattr(terms, 'immediateAbsorb', None)
    if absorb is None:
        return gone
    absorb = np.asarray(absorb, dtype=float)
    if absorb.ndim != 2 or absorb.shape[0] == 0:
        return gone
    n = min(int(nstate), absorb.shape[0], absorb.shape[1])
    gone[:n] = np.diag(absorb)[:n] == 0.0
    return gone


def _reach(row, coord_class, coord_station, njobs, K):
    """
    The largest ``row @ x`` the population can produce, ignoring the coupling.

    Used to drop a cap that can never bind. The bound is per CLASS rather than
    the whole population times the heaviest weight, which is what the region-only
    version did: that bound never prunes a PER-CLASS cap set to its own class
    population -- exactly the row ``refreshCapacity`` derives for every closed
    model -- so every such model carried a row that could not bind and an
    unknown with no equation to pin it.
    """
    total = 0.0
    for r in range(K):
        sel = (coord_class == r) & (coord_station >= 0)
        if not sel.any():
            continue
        w = float(np.max(row[sel]))
        if w <= 0:
            continue
        nr = float(njobs[r]) if r < njobs.size else np.inf
        if not np.isfinite(nr):
            return np.inf
        total += w * nr
    return total


def capacity_constraints(sn, terms):
    """
    Every capacity limit in the model as a linear constraint on the fluid state.

    TWO DECLARATIONS, ONE FAMILY OF ROWS. A finite capacity region caps a SET of
    stations jointly (``addRegion``); a station cap limits one station's own
    buffer (``setCapacity``/``setClassCapacity``). LINE stores the region form
    four ways -- a region-global job cap, a per-class job cap, a memory budget
    with per-class sizes, and an arbitrary linear pair (A,b) -- and the station
    form two, a total and a per-class buffer. All six are the same object once
    written against the state:

        Arow x <= b,   Arow(s) = the weight of the class of coordinate s,
                                 zero outside the stations the limit covers

    so the solver carries one mechanism rather than six.

    WHAT DIFFERS IS NOT THE ROW BUT WHERE THE BLOCKED JOB GOES, which is decided
    per admission event in CAPACITY_GATES, not here.

    WHAT IS REFUSED, AND WHY IT IS NOT A CONSTRAINT. Only a waiting queue is a
    constraint on this drift. Under WAITQ the job waits and is admitted later, so
    the population is conserved and only the admission FLOW is throttled -- which
    is what an algebraic equation on the same state can express. DROP destroys
    the job for a class that may not lose one, BAS/BBS/RSRD add a blocked-server
    state to the upstream station, and the retrial rules move the job to an
    orbit. Each changes the event set itself, so each needs a different drift
    rather than a constraint on this one, and is refused by name.
    """
    M, K, nstate = terms.M, terms.K, terms.nstate
    con = CapacityConstraints(nstate)
    nregions = int(getattr(sn, 'nregions', 0) or 0)

    coord_station = np.full(nstate, -1, dtype=int)
    coord_class = np.zeros(nstate, dtype=int)
    for i in range(M):
        for c in range(K):
            if terms.Kic[i, c] > 0:
                lo = terms.q_indices[i, c]
                hi = lo + terms.Kic[i, c]
                coord_station[lo:hi] = i
                coord_class[lo:hi] = c
    con.klass = coord_class
    con.nregions = nregions
    con.member = np.zeros((max(nregions, 1), nstate), dtype=bool)

    rows, bs, rgs, sts, cls, stgd, labels = [], [], [], [], [], [], []
    regionrule = np.asarray(getattr(sn, 'regionrule', np.zeros((0, 0))), dtype=float)
    regionweight = np.asarray(getattr(sn, 'regionweight', np.zeros((0, 0))), dtype=float)
    regionsz = np.asarray(getattr(sn, 'regionsz', np.zeros((0, 0))), dtype=float)
    regionmembers = getattr(sn, 'regionmembers', []) or []
    regionmaxmem = getattr(sn, 'regionmaxmem', []) or []
    regionlincon = getattr(sn, 'regionlincon', []) or []
    region = getattr(sn, 'region', []) or []

    for f in range(nregions):
        members = np.zeros(M, dtype=bool)
        if f < len(regionmembers) and regionmembers[f] is not None:
            members = np.asarray(regionmembers[f], dtype=bool).ravel()
        if not members.any():
            continue
        in_region = (coord_station >= 0) & members[np.maximum(coord_station, 0)]
        con.member[f, :] = in_region

        # -- refusals, named before anything is built ------------------------
        for r in range(K):
            if regionrule.ndim == 2 and regionrule.shape[0] > f and regionrule.shape[1] > r:
                if float(regionrule[f, r]) != _WAITQ:
                    raise ValueError(
                        "Region %d applies a drop rule other than a waiting queue to class %d. "
                        "Only a waiting queue is a constraint on the fluid drift: it conserves "
                        "the population and throttles the admission flow. The other rules change "
                        "the event set instead, so they need a different drift rather than an "
                        "algebraic equation on this one." % (f + 1, r + 1))
        if regionweight.ndim == 2 and regionweight.shape[0] > f:
            wf = regionweight[f, :]
            if np.any(np.abs(wf - 1.0) > _ZERO):
                raise ValueError(
                    "Region %d sets per-class admission weights, which decide WHICH blocked "
                    "class enters when capacity frees up. The constraint form throttles the "
                    "admission flow in proportion to its own rate and carries no such priority, "
                    "so the weights would be ignored silently. Use SolverCTMC, SolverJMT or "
                    "SolverSSA." % (f + 1))

        Rmat = np.asarray(region[f], dtype=float) if f < len(region) and region[f] is not None \
            else np.zeros((0, 0))
        member_row = int(np.where(members)[0][0])

        # -- 1. region-global job cap ----------------------------------------
        if Rmat.ndim == 2 and Rmat.shape[1] >= K + 1:
            gcap = float(Rmat[member_row, K])
            if np.isfinite(gcap) and gcap != _UNBOUNDED and gcap >= 0:
                row = np.zeros(nstate)
                row[in_region] = 1.0
                rows.append(row); bs.append(gcap); rgs.append(f); sts.append(-1)
                cls.append(-1); stgd.append(True)
                labels.append("region %d global job cap" % (f + 1))

        # -- 2. per-class job caps -------------------------------------------
        if Rmat.ndim == 2 and Rmat.size:
            for r in range(min(K, Rmat.shape[1])):
                ccap = float(Rmat[member_row, r])
                if not np.isfinite(ccap) or ccap == _UNBOUNDED or ccap < 0:
                    continue
                sel = in_region & (coord_class == r)
                if not sel.any():
                    continue
                row = np.zeros(nstate)
                row[sel] = 1.0
                rows.append(row); bs.append(ccap); rgs.append(f); sts.append(-1)
                cls.append(r); stgd.append(True)
                labels.append("region %d class %d job cap" % (f + 1, r + 1))

        # -- 3. region-global memory budget ----------------------------------
        # The per-class memory limit has already been folded into the per-class
        # job cap above by the struct refresh (memjobs = floor(maxMem/size)), so
        # only the region-global budget is left: the job-count row with each
        # class weighted by its size.
        if f < len(regionmaxmem) and regionmaxmem[f] is not None:
            mem = np.asarray(regionmaxmem[f], dtype=float).ravel()
            if member_row < mem.size:
                gmem = float(mem[member_row])
                if np.isfinite(gmem) and gmem != _UNBOUNDED and gmem >= 0:
                    sz = np.ones(K)
                    if regionsz.ndim == 2 and regionsz.shape[0] > f:
                        sz = regionsz[f, :]
                    row = np.zeros(nstate)
                    for r in range(K):
                        sel = in_region & (coord_class == r)
                        row[sel] = sz[r] if r < sz.size else 1.0
                    if np.any(row != 0):
                        rows.append(row); bs.append(gmem); rgs.append(f); sts.append(-1)
                        cls.append(-1); stgd.append(True)
                        labels.append("region %d memory budget" % (f + 1))

        # -- 4. explicit linear constraints ----------------------------------
        if f < len(regionlincon) and regionlincon[f] is not None:
            Alin = np.atleast_2d(np.asarray(regionlincon[f][0], dtype=float))
            blin = np.asarray(regionlincon[f][1], dtype=float).ravel()
            for c in range(Alin.shape[0]):
                row = np.zeros(nstate)
                for r in range(min(K, Alin.shape[1])):
                    if Alin[c, r] != 0:
                        sel = in_region & (coord_class == r)
                        row[sel] = Alin[c, r]
                if np.any(row != 0):
                    rows.append(row); bs.append(float(blin[c])); rgs.append(f); sts.append(-1)
                    cls.append(-1); stgd.append(True)
                    labels.append("region %d linear constraint %d" % (f + 1, c + 1))

    # -- 5. per-station buffers ----------------------------------------------
    # A STATION CAP IS THE ONE-STATION CASE OF THE SAME ROW, and every fluid
    # method other than this one ignores it outright: nothing in the FLD tree
    # reads sn.cap or sn.classcap, so a capped station was integrated as an
    # unbounded one and reported more jobs in the buffer than the buffer holds.
    cap = np.asarray(getattr(sn, 'cap', np.zeros(0)), dtype=float).ravel()
    classcap = np.atleast_2d(np.asarray(getattr(sn, 'classcap', np.zeros((0, 0))), dtype=float))
    for i in range(M):
        if terms.sched[i] == SchedStrategy.EXT:
            continue   # a source holds no jobs, so its cap caps nothing
        at_station = coord_station == i
        if not at_station.any():
            continue
        if i < cap.size:
            gcap = float(cap[i])
            if np.isfinite(gcap) and gcap >= 0:
                row = np.zeros(nstate)
                row[at_station] = 1.0
                rows.append(row); bs.append(gcap); rgs.append(-1); sts.append(i)
                cls.append(-1); stgd.append(False)
                labels.append("station %d buffer" % (i + 1))
        if i < classcap.shape[0]:
            for r in range(min(K, classcap.shape[1])):
                ccap = float(classcap[i, r])
                if not np.isfinite(ccap) or ccap < 0:
                    continue
                sel = at_station & (coord_class == r)
                if not sel.any():
                    continue
                row = np.zeros(nstate)
                row[sel] = 1.0
                rows.append(row); bs.append(ccap); rgs.append(-1); sts.append(i)
                cls.append(r); stgd.append(False)
                labels.append("station %d class %d buffer" % (i + 1, r + 1))

    if not rows:
        return con

    A = np.vstack(rows)
    b = np.asarray(bs, dtype=float)
    rg = np.asarray(rgs, dtype=int)
    st = np.asarray(sts, dtype=int)
    kl = np.asarray(cls, dtype=int)
    sd = np.asarray(stgd, dtype=bool)

    # A CAP THE POPULATION CANNOT REACH IS NOT A CAP, and dropping it here keeps
    # it out of the active-set loop, where it would be tested on every pass and
    # never bind while its multiplier stayed an unknown with no equation to pin
    # it. An open class makes the reachable total unbounded, so nothing on a row
    # an open class weighs is pruned.
    njobs = np.asarray(getattr(sn, 'njobs', np.zeros(0)), dtype=float).ravel()
    keep = np.ones(b.size, dtype=bool)
    for c in range(b.size):
        if _reach(A[c, :], coord_class, coord_station, njobs, K) <= b[c] + _ZERO:
            keep[c] = False

    # A CAP ON COORDINATES THE IMMEDIATE REDUCTION FOLDED AWAY IS VACUOUS, NOT
    # MALFORMED. hide_immediate stochastic-complements an Immediate-rate
    # coordinate out of the event set -- the MMT transform's zero-service Join is
    # one, until the fork-join fixed point gives it a synchronisation delay -- and
    # the reduced drift then holds no mass there and has no event landing on it.
    # Left in, such a row reaches CapacityGates with no gating event and is
    # reported as a limit the model cannot approach, which is the right message
    # for a station that really does hold jobs and the wrong one here: this
    # station holds none, so its buffer is satisfied identically.
    gone = _eliminated_coords(terms, A.shape[1])
    for c in range(b.size):
        if keep[c] and not np.any((A[c, :] != 0.0) & ~gone):
            keep[c] = False

    # THE SAME ROW TWICE IS A SINGULAR NEWTON SYSTEM, not a redundancy the least
    # squares absorbs: two identical rows both bind, each takes a multiplier, and
    # nothing distinguishes them. refreshCapacity derives classcap from cap, so a
    # single-class model declares the station total and the class buffer as the
    # same row -- the common case, not a corner one. The TIGHTER bound survives;
    # on a tie the region row does, because a region is an explicit construct
    # with a waiting room of its own while a station cap is a buffer length.
    for c in range(b.size):
        if not keep[c]:
            continue
        for d in range(c + 1, b.size):
            if not keep[d] or not np.allclose(A[c, :], A[d, :], atol=_ZERO, rtol=0.0):
                continue
            takeover = b[d] < b[c] - _ZERO or (abs(b[d] - b[c]) <= _ZERO and sd[d] and not sd[c])
            if takeover:
                b[c] = b[d]; rg[c] = rg[d]; st[c] = st[d]; kl[c] = kl[d]; sd[c] = sd[d]
                labels[c] = labels[d]
            keep[d] = False

    # A ROW ITS OWN PER-CLASS ROWS ALREADY IMPLY IS RANK, NOT INFORMATION. LINE
    # derives the station total from the per-class buffers (cap = sum classcap), so
    # a two-class station capped 3 and 3 declares a total of 6 as well -- and the
    # total row is EXACTLY the sum of the two class rows, with a bound exactly
    # their sum. All three then bind together, the constraint block has rank 2
    # with 3 multipliers, and nothing distinguishes them: the least-squares step
    # still lands on the right STATE, but the multipliers it reports are one
    # arbitrary point of a line of solutions. Dropping the implied row makes them
    # unique again. Only an exact implication is dropped -- the per-class rows must
    # cover every class the candidate weighs, at no smaller a weight, and their
    # bounds must sum to no more than its own -- so a total TIGHTER than the sum of
    # its parts (which is a genuine extra constraint) survives.
    for c in range(b.size):
        if not keep[c] or kl[c] != -1:
            continue
        classes = [r for r in range(K)
                   if np.any((coord_class == r) & (coord_station >= 0) & (A[c, :] > 0))]
        if not classes:
            continue
        implied, budget = True, 0.0
        for r in classes:
            sel = (coord_class == r) & (coord_station >= 0) & (A[c, :] > 0)
            part = -1
            for d in range(b.size):
                if not keep[d] or d == c or kl[d] != r:
                    continue
                if rg[d] != rg[c] or st[d] != st[c]:
                    continue
                if np.all(A[d, sel] >= A[c, sel] - _ZERO):
                    part = d
                    break
            if part < 0:
                implied = False
                break
            budget += b[part]
        if implied and budget <= b[c] + _ZERO:
            keep[c] = False

    A, b, rg, st, kl, sd = A[keep, :], b[keep], rg[keep], st[keep], kl[keep], sd[keep]
    labels = [labels[c] for c in range(len(labels)) if keep[c]]

    # WHICH STATION RULES ARE A CONSTRAINT ON THIS DRIFT, and which are a
    # different event set. THE RULE IS NOT WHAT DECIDES THE SEMANTICS -- the class
    # type is, exactly as State.arrivalIsLost decides it for every other solver:
    # a closed job is never lost (the upstream departure is disabled instead) and
    # an open one is never held (it simply never enters). So a waiting queue and a
    # drop declaration are BOTH constraints here and differ only in which of those
    # two the class already implies; refreshCapacity in fact declares DROP by
    # default at a capped station reached by an open class, so refusing DROP would
    # refuse every open loss model. What is refused is the rules that add STATE:
    # BAS/BBS/RSRD give the upstream station a blocked-server state and the
    # retrial rules add an orbit, neither of which this drift carries.
    #
    # A rule is only a contradiction where the cap can BIND, so this runs on the
    # surviving rows: a declaration on a buffer the population can never fill
    # describes nothing, and refusing it would reject models that behave
    # identically with and without it.
    droprule = np.atleast_2d(np.asarray(getattr(sn, 'droprule', np.zeros((0, 0))), dtype=float))
    for c in range(b.size):
        if st[c] < 0 or st[c] >= droprule.shape[0]:
            continue
        for r in range(min(K, droprule.shape[1])):
            sel = (coord_station == st[c]) & (coord_class == r)
            if not sel.any() or float(np.max(A[c, sel])) <= 0:
                continue
            rule = float(droprule[st[c], r])
            if rule not in (_WAITQ, _DROP):
                raise ValueError(
                    "Station %d applies %s to class %d, and its buffer binds. Only a waiting "
                    "queue or a drop is a constraint on this drift: the first conserves the "
                    "population and throttles the admission flow, the second discards the flow "
                    "the cap will not take. BAS/BBS/RSRD add a blocked-server state to the "
                    "upstream station and the retrial rules add an orbit, so each needs a "
                    "different drift rather than an algebraic equation on this one. Use "
                    "SolverCTMC, SolverJMT, SolverSSA or SolverLDES."
                    % (int(st[c]) + 1, _DROP_NAMES.get(rule, "drop rule %g" % rule), r + 1))

    con.A, con.b, con.region, con.station, con.cls, con.staged = A, b, rg, st, kl, sd
    con.label = labels
    return con


class CapacityGates(object):
    """Which events each cap throttles, and what happens to the mass it stops."""

    __slots__ = ("gate", "held", "loss", "Dn", "DnExt", "Dp")

    def __init__(self, ncon, nevents):
        self.gate = np.zeros((ncon, nevents), dtype=bool)
        # the two ways a cap that does not stage can stop an event: hold the job
        # at the upstream station, or lose it. Masks rather than a mode code so
        # that the four ports read the same way.
        self.held = np.zeros((ncon, nevents), dtype=bool)
        self.loss = np.zeros((ncon, nevents), dtype=bool)
        # the jump matrix split in THREE, so that an admission can remove mass
        # upstream at one rate and deliver it at another. Dn carries the removal
        # at real stations, DnExt the removal from the EXT source pool -- which a
        # LOST arrival must be returned to, since that coordinate is a
        # normalisation and not a population -- and Dp the arrival. None where no
        # cap exists, which keeps the uncapped drift a single matrix product.
        self.Dn = None
        self.DnExt = None
        self.Dp = None


def capacity_gates(sn, terms, con):
    """
    Per cap and per event: is this event an admission the cap throttles, and
    WHERE DOES THE STOPPED MASS GO.

    An admission is an event that pushes the constrained quantity UP, read off
    ``A D`` rather than off the topology so that a per-class or memory-weighted
    row picks out its own admissions with no extra code.

    THE THREE ANSWERS ARE THE MODEL, AND THEY ARE NOT INTERCHANGEABLE. LINE's own
    semantics decides which one a cap gets, and the choice is visible in the
    answer -- a job held upstream is still counted at that station, a lost job is
    counted nowhere and breaks flow balance across the cap on purpose, a staged
    job is counted in neither and shows up as `blocked`:

      STAGED   a finite capacity region under a waiting queue. The job COMPLETES
               upstream service and waits outside the region, which is what JMT
               and LDES simulate; the upstream station empties exactly as it
               would with no region. See CAPACITY_STAGING.
      HELD     a station buffer reached by a CLOSED class. State.arrivalIsLost
               refuses to lose a closed job -- population conservation is a
               defining invariant -- and returns an empty successor instead,
               which DISABLES the upstream departure until room frees. The job is
               therefore still at the upstream station, in service as far as that
               station's own metrics are concerned, so the fluid analogue is to
               scale the WHOLE event: both the removal upstream and the arrival.
      LOSS     a station buffer reached by an OPEN class. The same predicate
               loses it: the external stream is memoryless, so a job that finds
               the buffer full simply never enters, the arrival event still fires
               (ArvR counts the offered job) and only the carried flow is
               admitted. The fluid analogue scales the ARRIVAL leg alone and
               destroys the difference.
    """
    nevents = int(np.size(terms.rateBase))
    ncon = int(con.b.size)
    gates = CapacityGates(ncon, nevents)
    if ncon == 0:
        return gates
    tol = np.sqrt(_ZERO)
    gates.Dn = np.minimum(terms.D, 0.0)
    gates.Dp = np.maximum(terms.D, 0.0)
    gates.DnExt = np.zeros_like(gates.Dn)
    for i in range(terms.M):
        if terms.sched[i] == SchedStrategy.EXT:
            blk = terms.stationBlock[i]
            if blk.size:
                gates.DnExt[blk, :] = gates.Dn[blk, :]
                gates.Dn[blk, :] = 0.0
    Dp = gates.Dp
    njobs = np.asarray(getattr(sn, 'njobs', np.zeros(0)), dtype=float).ravel()
    for c in range(ncon):
        delta = con.A[c, :] @ terms.D
        for e in range(nevents):
            if delta[e] <= tol:
                continue
            gates.gate[c, e] = True
            if con.staged[c]:
                continue
            # the class is read off the coordinate the mass LANDS on, inside the
            # capped station: a class switch on entry would otherwise ask the
            # class the job is leaving behind whether it may be lost
            landing = np.where((Dp[:, e] > tol) & (con.A[c, :] > 0))[0]
            k = int(con.klass[landing[0]]) if landing.size else -1
            open_class = k >= 0 and k < njobs.size and not np.isfinite(njobs[k])
            gates.loss[c, e] = open_class
            gates.held[c, e] = not open_class
        if not gates.gate[c, :].any():
            raise ValueError(
                "No event increases %s, so the cap can never be approached and there is no "
                "admission flow for the constraint to throttle. This is a malformed limit "
                "rather than a solvable one." % con.label[c])
    return gates


class CapacityStaging(object):
    """The waiting room outside each capped region, as fluid coordinates."""

    __slots__ = ("n", "region", "klass", "idx", "adm", "adm_region",
                 "adm_stage", "gated_by", "Dn", "Dp")

    def __init__(self, nevents, K, nregions, ncon=0):
        self.n = 0
        self.region = np.zeros(0, dtype=int)
        self.klass = np.zeros(0, dtype=int)
        self.idx = np.zeros((nregions, K), dtype=int)
        self.adm = np.zeros(nevents, dtype=bool)
        self.adm_region = np.zeros(nevents, dtype=int)
        self.adm_stage = np.zeros(nevents, dtype=int)
        self.gated_by = np.zeros((ncon, 0), dtype=bool)
        self.Dn = None
        self.Dp = None


def capacity_staging(terms, con):
    """
    The waiting room outside a capped region.

    WHY THE CONSTRAINT ALONE IS NOT ENOUGH, FOR A REGION. Throttling a region's
    admission events does hold its population at the cap, but it holds it by
    slowing the UPSTREAM STATION'S COMPLETIONS -- an admission event IS that
    station finishing a job -- so blocked mass piles up at a station it has
    already finished being served by. Where that station is a delay the error
    shows as a broken Little's law. Under a waiting queue the job COMPLETES
    upstream service and then waits; it is somewhere else, and the model needs
    somewhere else to put it.

    A STATION BUFFER GETS NO ROOM, and that is not an omission. There the job
    genuinely does stay where it was: LINE disables the upstream departure
    (State.arrivalIsLost) rather than moving the job out, so the blocked mass is
    still at the upstream station and still counted there. Giving it a room would
    move mass the reference keeps in place. See CAPACITY_GATES.

    Each capped region gains one coordinate per class, and every admission into
    it splits in two: upstream -> staging at the nominal rate, so the upstream
    station empties exactly as it would with no region, and staging -> region at
    ``theta * s``, the throttled leg. The two jumps sum to the original, so only
    where the mass rests changes.
    """
    nevents = int(np.size(terms.rateBase))
    K = terms.K
    nregions = con.nregions
    ncon = int(con.b.size)
    stg = CapacityStaging(nevents, K, max(nregions, 1), ncon)
    if nregions == 0 or con.empty or not np.any(con.staged):
        return stg

    D = terms.D
    stg.Dn = np.minimum(D, 0.0)
    stg.Dp = np.maximum(D, 0.0)
    stg.idx = np.zeros((nregions, K), dtype=int)
    tol = np.sqrt(_ZERO)

    staged_regions = set(int(f) for f in con.region[con.staged] if f >= 0)
    region, klass = [], []
    for f in range(nregions):
        if f not in staged_regions:
            continue
        member_row = con.member[f, :]
        if not member_row.any():
            continue
        # net change of this region's population per event: positive means the
        # event brings mass in from outside, which is what the queue feeds
        delta = member_row.astype(float) @ D
        for e in range(nevents):
            if delta[e] <= tol:
                continue
            # the class is read off the coordinate the mass LANDS on, inside the
            # region: a class switch on entry would otherwise stage the job under
            # the class it is leaving behind
            landing = np.where((stg.Dp[:, e] > tol) & member_row)[0]
            if landing.size == 0:
                continue
            c = int(con.klass[landing[0]])
            if c < 0 or c >= K:
                continue
            if stg.idx[f, c] == 0:
                region.append(f)
                klass.append(c)
                stg.idx[f, c] = len(region)   # 1-based; 0 means none
            stg.adm[e] = True
            stg.adm_region[e] = f
            stg.adm_stage[e] = stg.idx[f, c] - 1

    stg.n = len(region)
    stg.region = np.asarray(region, dtype=int)
    stg.klass = np.asarray(klass, dtype=int)

    # WHICH ROWS GATE WHICH ROOM. One drain rate per room and one equality per
    # row, and the two are not in bijection: a region-global cap gates every room
    # of its region, a per-class cap only the room of its class. A room gated by
    # several ACTIVE rows drains at the harmonic composition of their rates --
    # see _RESIDUAL -- which is what lets two caps of one region bind at once,
    # the case a single per-region throttle had to refuse.
    stg.gated_by = np.zeros((ncon, stg.n), dtype=bool)
    for c in range(ncon):
        if not con.staged[c] or con.region[c] < 0:
            continue
        for j in range(stg.n):
            if stg.region[j] != con.region[c]:
                continue
            sel = con.member[con.region[c], :] & (con.klass == stg.klass[j])
            if sel.any() and float(np.max(con.A[c, sel])) > 0:
                stg.gated_by[c, j] = True
    return stg


def capacity_extend(con, stg, terms):
    """
    Extend every cap to the staging coordinates that hold mass INSIDE it.

    A waiting room is outside the region it feeds, which is the whole point of
    it -- but it is not outside every OTHER limit. Where two regions overlap, an
    admission into the inner one is an INTERNAL move of the outer one: the job
    leaves a station of the outer region, waits, and re-enters a station of the
    same outer region, never having left it. Counting only the state coordinates
    would take that mass out of the outer cap for as long as it waits, so the
    outer cap would be met on paper while the region actually held more; and the
    Newton system that results is inconsistent rather than merely inexact -- two
    overlapping regions stalled at residual 5e-1, with the inner cap exceeded.

    A room counts toward a row when the row weighs the room's DESTINATION and
    also weighs every station that FEEDS it -- that is exactly "the job was
    inside and stays inside". A room fed from outside is a queue at the door and
    counts nowhere, as before.
    """
    ncon = int(con.b.size)
    con.As = np.zeros((ncon, stg.n))
    if ncon == 0 or stg.n == 0:
        return con
    nstate = con.A.shape[1]
    tol = np.sqrt(_ZERO)
    Dn = np.minimum(terms.D, 0.0)
    # which stations feed each room, and which coordinate its mass lands on
    feeds = [set() for _ in range(stg.n)]
    dest = [-1] * stg.n
    for e in np.where(stg.adm)[0]:
        j = int(stg.adm_stage[e])
        for s in np.where(Dn[:, e] < -tol)[0]:
            feeds[j].add(int(s))
        landing = np.where(np.maximum(terms.D[:, e], 0.0) > tol)[0]
        for s in landing:
            if con.member[stg.region[j], s]:
                dest[j] = int(s)
                break
    for c in range(ncon):
        for j in range(stg.n):
            if dest[j] < 0:
                continue
            w = float(con.A[c, dest[j]])
            if w <= 0:
                continue
            if all(float(np.max(con.A[c, s])) > 0 for s in feeds[j]) and feeds[j]:
                con.As[c, j] = w
    return con


# ---------------------------------------------------------------------------
# the solver
# ---------------------------------------------------------------------------

class DaeSolver(MinNormalSolver):
    """The min-normal closure solved as one differential-algebraic system."""

    def solve(self) -> FLDResult:
        start_time = time.time()
        self._check_scheds()
        terms = self.build_moment_terms()
        M, K, nstate = terms.M, terms.K, terms.nstate
        sn = self.sn
        cfg = getattr(self.options, 'config', None) or {}

        # The simultaneous solve carries one Lyapunov solve per residual
        # evaluation and takes a finite-difference Jacobian over
        # nstate + nclosable unknowns, so its cost is cubic per evaluation and
        # quartic overall. Affordable at the scale the moment methods run at,
        # but the crossover is lower than the 200 `minnormal` permits, so it
        # gets its own limit.
        maxstate = int(cfg.get('dae_maxstate') or 100)
        if nstate > maxstate:
            raise ValueError(
                "The dae method solves a %d-unknown algebraic system with a finite-difference "
                "Jacobian, above the limit of %d set by options.config['dae_maxstate']. Raise "
                "that limit, or use options.method='minnormal' for the same closure by "
                "successive substitution." % (nstate, maxstate))

        # DPS and GPS close on the covariance BETWEEN a station's class
        # coordinates, not on the station total, so their closure state is a
        # matrix block rather than the scalar this solves for. Carrying those
        # blocks as Newton unknowns puts the quartic term back; carrying them as
        # a chord reintroduces the alternation this method exists to remove.
        for i in range(M):
            if terms.sched[i] in (SchedStrategy.DPS, SchedStrategy.GPS):
                raise ValueError(
                    "The dae method closes on the per-station variance only, but this model has "
                    "a DPS or GPS station whose share closes on the covariance BETWEEN its class "
                    "coordinates. Use options.method='minnormal', which carries those blocks "
                    "through its outer iteration.")

        # Caps first: staging depends on them, and conservation must count the
        # blocked mass staging holds. GATES decides, per cap and per event, where
        # the mass a cap stops actually goes -- held upstream, lost, or staged.
        con = capacity_constraints(sn, terms)
        ncon = int(con.b.size)
        gates = capacity_gates(sn, terms, con)
        stg = capacity_staging(terms, con)
        con = capacity_extend(con, stg, terms)
        C, Nvec = self._conservation(terms, stg)
        if C.size:
            leak = float(np.max(np.abs(C[:, :nstate] @ terms.D))) if C.shape[0] else 0.0
            if leak > np.sqrt(_ZERO):
                raise ValueError(
                    "The event set does not conserve a closed chain: the largest population leak "
                    "per unit rate is %g, where it must be zero. The conservation constraint "
                    "would contradict the drift rather than complete it." % leak)

        # -- which stations have a variance that enters the drift ------------
        # A delay or a source has no min() to close. A station that cannot fill
        # its servers has min(n,c) = n on its whole support, so closing it is not
        # an improvement but an error. Those are held at zero and are not
        # unknowns, which keeps the Newton system as small as the closure is.
        closable = []
        for i in range(M):
            if terms.sched[i] in (SchedStrategy.EXT, SchedStrategy.INF):
                continue
            if not np.isfinite(terms.nservers[i]):
                continue
            if terms.minExact[i] or terms.stationBlock[i].size == 0:
                continue
            closable.append(i)
        cidx = np.asarray(closable, dtype=int)
        covblk = [None] * M   # always empty here: DPS/GPS are refused above

        # A SOURCE COORDINATE HAS NO FIXED POINT, so its drift row is not an
        # equation and must not be asked to vanish. In the closing
        # representation an open class's departures are routed back onto the EXT
        # pseudo-station, so that coordinate ACCUMULATES the mass that left the
        # system: its drift is the throughput, permanently, and `0 = D r` is
        # unsatisfiable there by construction. Dropping those rows is exact
        # rather than a concession, because the EXT rate factor of a
        # single-phase source is 1 identically -- verified above by
        # `build_moment_terms`, which refuses a multi-phase source outright --
        # so no other row depends on the value of that coordinate.
        ext_coord = np.zeros(nstate, dtype=bool)
        for i in range(M):
            if terms.sched[i] == SchedStrategy.EXT:
                blk = terms.stationBlock[i]
                if blk.size:
                    ext_coord[blk] = True

        # -- seed ------------------------------------------------------------
        # Newton needs a point in the basin, not an answer. One first-order
        # solve supplies it, at the cost of the single integration this method
        # exists to avoid repeating twenty times.
        x0 = self._compute_initial_state(M, K, terms.phases)
        xseed, xvec_t_seed, t_seed = self._integrate(
            M, K, terms.Mu, terms.Phi, terms.phases, terms.rt, terms.nservers,
            terms.sched, terms.schedparam, x0, np.zeros(M), [None] * M)
        x = np.asarray(xseed, dtype=float).ravel().copy()
        # and reset them to the unit normalisation mass the layout intends: the
        # seed integration ran to a horizon, so what it left there is however
        # much work happened to have crossed the system, not a state.
        x[ext_coord] = np.asarray(x0, dtype=float).ravel()[ext_coord]

        # THE VARIANCE IS SEEDED POSITIVE, which is why this route has no kink
        # probe. sigma2 = 0 is where min(n,c) has no derivative and a saturated
        # model's first-order fixed point sits exactly there. The station mean is
        # an O(N) starting value in the right units -- the variance a Poisson
        # population of that mean would have -- and costs no Lyapunov solve.
        sigma2 = np.zeros(M)
        for i in cidx:
            blk = terms.stationBlock[i]
            if blk.size:
                sigma2[i] = max(_FINE, float(np.sum(x[blk])))

        tol = float(getattr(self.options, 'tol', 0) or 0)
        if not np.isfinite(tol) or tol <= 0:
            tol = _FINE
        newton_max = 50
        it_max = getattr(self.options, 'iter_max', None)
        if it_max:
            newton_max = max(newton_max, int(it_max))

        # -- steady state: one simultaneous algebraic solve, per active set ---
        # A capacity constraint is an inequality, and an inequality has no
        # residual to hand a Newton solver -- only the caps that actually bind
        # become equations. Each pass is a complete simultaneous solve, so the
        # loop iterates over WHICH caps bind, not over the closure.
        # THE ACTIVE SET STARTS EMPTY, and the first pass therefore asks for the
        # UNCONSTRAINED fixed point -- which is what makes the answer of a model
        # whose caps bind independent of how far outside the seed happened to land.
        # It is also not always solvable: an overloaded M/M/1/K has NO equilibrium
        # without its cap, so that pass fails rather than converging, and the caps
        # the seed violates are then seeded into the set instead (see below). Doing
        # that up front for every model looked cheaper and changed answers: a region
        # of 8 over two identical queues settled 7.43/0.57 from the clipped seed and
        # 4/4 from the unconstrained one, and only the second is the reference's.
        active = np.zeros(0, dtype=int)
        seeded_from_failure = False
        active_solved = np.zeros(0, dtype=int)
        iters = 0
        mult = np.zeros(0)
        sg = np.zeros(stg.n)
        aset_max = max(4, 2 * ncon + 2)
        converged, resnorm = False, np.inf
        clamped = False
        # the last iterate that WAS a fixed point, and the whole answer read off
        # it. A pass that fails to converge leaves an iterate that is not a fixed
        # point of anything, and reading the next active set off it is how a single
        # bad pass turned into a walk through unrelated capacity combinations.
        x_ok, sigma2_ok = x.copy(), sigma2.copy()
        best = None
        u = None
        aset = 0
        for aset in range(aset_max):
            # ONE MULTIPLIER PER ACTIVE ROW, not per region. A region's waiting
            # room used to carry a single drain rate, so two caps of one region
            # were two equalities against one control and the case was refused.
            # A room gated by several active rows now drains at the harmonic
            # composition of their rates and a held cap composes as a product of
            # fractions, so each row keeps a control of its own; see _RESIDUAL.
            has_clamp = any(not con.staged[int(c)] for c in active)
            # Seed each waiting room with the mass that does not fit, and every
            # multiplier at unity -- an unthrottled fraction for a held or lost
            # cap, and the drain rate the region route has always started from.
            sg0 = np.zeros(stg.n)
            for k in range(active.size):
                c = int(active[k])
                if not con.staged[c]:
                    continue
                excess = max(0.0, float(con.A[c, :] @ x - con.b[c]))
                sel = np.where(stg.gated_by[c, :])[0]
                if sel.size:
                    sg0[sel] = excess / sel.size
            u0 = np.concatenate([x, sg0, sigma2[cidx] if cidx.size else np.zeros(0),
                                 np.ones(active.size)])
            # THE CLAMPED COVARIANCE IS A FALLBACK, NOT THE DEFAULT. A cap that
            # holds or loses fixes its own combination of the state, so the honest
            # LNA puts no fluctuation there (_CLAMP_TANGENT) -- but the truth is
            # neither that nor the unprojected variance: the population under a cap
            # follows a TRUNCATED distribution, whose variance is smaller than the
            # unconstrained one and larger than zero. The unprojected solve is what
            # every other fluid method computes, so it is what runs first and what
            # keeps the answers of the models that already worked; the projection
            # is tried only when the unprojected system has no stationary
            # covariance at all, which is exactly the neutral case an overloaded
            # loss station produces. Deciding once per PASS rather than per
            # residual matters: a projection that switched on and off between
            # iterates would give Newton a discontinuous system to converge on.
            nohyp = None
            for use_clamp in ((False, True) if has_clamp else (False,)):
                ctx = dict(terms=terms, C=C, Nvec=Nvec, cidx=cidx, covblk=covblk,
                           nstate=nstate, conA=con.A, conAs=con.As, conb=con.b,
                           conStaged=con.staged, gates=gates, stg=stg, active=active,
                           M=M, driftRows=~ext_coord,
                           clampT=self._clamp_tangent(con, active, terms) if use_clamp else None)
                try:
                    u, nit, converged, resnorm = self._newton(ctx, u0, tol, newton_max, nstate)
                except FluidNonHyperbolicError as exc:
                    if use_clamp or not has_clamp:
                        nohyp = exc
                        break
                    continue
                iters += nit
                clamped = use_clamp
                if converged:
                    break
            if nohyp is not None:
                # THE UNCONSTRAINED FIXED POINT NEED NOT EXIST. An overloaded open
                # station has no equilibrium until its buffer bounds it, so the pass
                # that asks for one fails and the caps the iterate violates -- or
                # left infinite, which no comparison catches -- are seeded into the
                # active set instead. Once. A second failure with caps already
                # bound is the model's answer and not a starting point to improve.
                cand = [c for c in range(ncon)
                        if c not in set(active.tolist())
                        and (not np.all(np.isfinite(x[con.A[c, :] > 0]))
                             or float(con.A[c, :] @ x) > con.b[c] + max(1e-9, tol))]
                if seeded_from_failure or not cand:
                    raise nohyp
                for c in cand:
                    val = float(con.A[c, :] @ x)
                    sel = con.A[c, :] > 0
                    npos = int(np.count_nonzero(sel))
                    if not npos or not con.b[c] > 0:
                        continue
                    if np.isfinite(val) and val > con.b[c]:
                        x[sel] *= con.b[c] / val
                    elif not np.isfinite(val):
                        x[sel] = con.b[c] / npos
                x[~np.isfinite(x)] = 0.0
                active = np.asarray(sorted(set(active.tolist()) | set(cand)), dtype=int)
                seeded_from_failure = True
                continue
            x = u[:nstate].copy()
            sg = u[nstate:nstate + stg.n].copy()
            sigma2 = np.zeros(M)
            if cidx.size:
                sigma2[cidx] = np.maximum(0.0, u[nstate + stg.n:nstate + stg.n + cidx.size])
            mult = np.maximum(0.0, u[nstate + stg.n + cidx.size:])
            # the set that produced THIS x and these multipliers, which is what
            # the metrics are read at; `active` below is the set to try NEXT
            active_solved = active
            if ncon == 0:
                break

            slack = con.b - con.A @ x - con.As @ sg
            if converged:
                x_ok, sigma2_ok = x.copy(), sigma2.copy()
                best = dict(x=x.copy(), sg=sg.copy(), sigma2=sigma2.copy(),
                            mult=mult.copy(), active=active, resnorm=resnorm,
                            clamped=clamped,
                            feasible=bool(np.all(slack >= -max(1e-9, tol))))
            violated = np.array([c for c in range(ncon)
                                 if slack[c] < -max(1e-9, tol) and c not in set(active.tolist())],
                                dtype=int)
            # THE RELEASE SIGNAL IS THE MULTIPLIER'S OWN UNITS, and the two kinds
            # do not share them. A held or lost cap throttles by a FRACTION, so a
            # fraction that came back above one was holding the flow down for no
            # reason and the cap never bound. A staged cap throttles by a RATE,
            # which has no such scale -- there the signal is a waiting room with
            # no blocked mass in it at all.
            released = []
            for k in range(active.size):
                c = int(active[k])
                if con.staged[c]:
                    rooms = np.where(stg.gated_by[c, :])[0]
                    if rooms.size == 0 or float(np.sum(sg[rooms])) < 1e-9:
                        released.append(c)
                elif mult[k] > 1.0 + max(1e-9, tol):
                    released.append(c)
            if not converged:
                # A FAILED PASS SAYS NOTHING ABOUT WHICH CAPS BIND. Its release
                # signal is still information -- a multiplier that ran above one
                # was on its way out of the active set -- but its state is not, so
                # no row is ADDED from it and the next pass restarts from the last
                # point that was a fixed point.
                violated = np.zeros(0, dtype=int)
                x, sigma2 = x_ok.copy(), sigma2_ok.copy()
                if not released:
                    break
            if violated.size == 0 and not released:
                break
            nxt = sorted((set(active.tolist()) | set(violated.tolist())) - set(released))
            active = np.asarray(nxt, dtype=int)

        # A CONVERGED POINT BEATS THE LAST ITERATE. The loop can end on a pass that
        # did not converge -- a cap the closure cannot hold at any multiplier will
        # be added, fail and be released for as long as the loop runs -- and the
        # last iterate of such a pass is not a fixed point of anything. Report the
        # last converged one instead, and say plainly which caps it does not meet
        # rather than presenting a cap-violating point as the answer.
        if best is not None and not converged:
            x, sg = best['x'], best['sg']
            sigma2, mult = best['sigma2'], best['mult']
            active_solved, resnorm, clamped = best['active'], best['resnorm'], best['clamped']
            converged = True
            if not best['feasible']:
                over = np.where((con.A @ x + con.As @ sg) > con.b + max(1e-9, tol))[0]
                warnings.warn(
                    "No fixed point of the closure satisfies %s. The closure wants more jobs there "
                    "than the cap allows and no admission multiplier holds it: the reported point "
                    "is the converged UNCONSTRAINED-in-that-cap fixed point, and it exceeds the "
                    "cap. Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES for this model."
                    % ", ".join(con.label[int(c)] for c in over), RuntimeWarning)

        if ncon > 0 and aset == aset_max - 1:
            warnings.warn("The active set did not settle: the same capacity constraints kept "
                          "binding and releasing. The reported point satisfies the last set tried.",
                          RuntimeWarning)
        if not converged:
            # Report rather than return a point that is not a fixed point. The
            # substitution route would simply have stopped at outer_max with no
            # indication, which is the failure mode this replaces.
            warnings.warn("The simultaneous closure solve stopped at residual %.3e after %d "
                          "Newton steps without reaching %.3e. The reported point is the last "
                          "iterate." % (resnorm, iters, tol), RuntimeWarning)

        # the variance as it entered the DRIFT, which is what any later solve on
        # this fixed point must close at
        sigma2drift = sigma2.copy()
        covdrift = list(covblk)

        ctx = dict(terms=terms, C=C, Nvec=Nvec, cidx=cidx, covblk=covblk,
                   nstate=nstate, conA=con.A, conAs=con.As, conb=con.b,
                   conStaged=con.staged,
                   gates=gates, stg=stg, active=active_solved, M=M,
                   driftRows=~ext_coord,
                   clampT=(self._clamp_tangent(con, active_solved, terms)
                           if clamped else None))
        # RMET is the flow that actually CROSSES each event, which is what the
        # throughput table must report, and RFIRE the rate the event fires at,
        # which is what the diffusion counts: a lost arrival fires and lands
        # nowhere, so the two differ by exactly the loss.
        _, rmet, rfire = self._residual(
            np.concatenate([x, sg, sigma2[cidx] if cidx.size else np.zeros(0), mult]),
            ctx, quiet=False)
        _, Sigma = self._sigma(x, sigma2, terms, cidx, covblk, rfire,
                               self._clamp_tangent(con, active_solved, terms)
                               if clamped else None)

        # -- transient -------------------------------------------------------
        Sigmat, QVart, tvar = None, None, None
        timespan = getattr(self.options, 'timespan', (0.0, np.inf))
        t = t_seed
        xvec_t = np.asarray(xvec_t_seed, dtype=float)
        switches = []
        if timespan is not None and np.isfinite(timespan[1]):
            # UNDER A CAP THE TRANSIENT IS A HYBRID DAE, integrated segment by
            # segment with the binding set updated at each located crossing. See
            # _TRANSIENT: what the steady state settles once with an active-set
            # loop, the trajectory settles again at every fill and every drain.
            t, xvec_t, Sigmat, tvar, switches = self._transient(
                terms, C, Nvec, cidx, covblk, sigma2drift, timespan, tol,
                np.asarray(xvec_t_seed, dtype=float), t_seed, con, gates, stg)
            if Sigmat is not None:
                QVart = self._qvar_t(Sigmat, terms)

        return self._assemble(terms, x, sg, sigma2, sigma2drift, covdrift, Sigma,
                              rmet, t, xvec_t, Sigmat, QVart, tvar, con, stg,
                              active_solved, mult, C, Nvec, iters, resnorm,
                              converged, start_time, switches)

    # ------------------------------------------------------------------ parts

    def _conservation(self, terms, stg):
        """
        Population conservation, one row per CLOSED chain.

        An open chain has no conserved population and contributes nothing. The
        EXT coordinates are excluded because the closing representation holds
        unit mass there as a normalisation constant, not as a job count.
        """
        M, K, nstate = terms.M, terms.K, terms.nstate
        coord_class = np.zeros(nstate, dtype=int)
        coord_station = np.full(nstate, -1, dtype=int)
        for i in range(M):
            for c in range(K):
                if terms.Kic[i, c] > 0:
                    lo = terms.q_indices[i, c]
                    hi = lo + terms.Kic[i, c]
                    coord_class[lo:hi] = c
                    coord_station[lo:hi] = i
        is_ext = np.array([terms.sched[i] == SchedStrategy.EXT for i in range(M)])

        nstg = stg.n
        chains = getattr(self.sn, 'chains', None)
        njobs = np.asarray(getattr(self.sn, 'njobs', []), dtype=float).ravel()
        rows, nvec = [], []
        if chains is not None and np.size(chains) > 0:
            chains = np.atleast_2d(np.asarray(chains))
            for ch in range(chains.shape[0]):
                inch = np.where(chains[ch, :] > 0)[0]
                if inch.size == 0:
                    continue
                Nch = float(np.sum(njobs[inch])) if njobs.size else 0.0
                if not np.isfinite(Nch) or Nch <= 0:
                    continue      # open chain: nothing is conserved
                sel = np.isin(coord_class, inch) & (coord_station >= 0) \
                    & ~is_ext[np.maximum(coord_station, 0)]
                if not sel.any():
                    continue
                row = np.zeros(nstate + nstg)
                row[:nstate][sel] = 1.0
                # A JOB IN THE WAITING QUEUE IS STILL IN THE CHAIN. Leaving the
                # staging coordinates out of this row would let the constraint
                # balance while the blocked mass quietly left the model.
                if nstg:
                    row[nstate + np.where(np.isin(stg.klass, inch))[0]] = 1.0
                rows.append(row)
                nvec.append(Nch)
        C = np.vstack(rows) if rows else np.zeros((0, nstate + nstg))
        return C, np.asarray(nvec, dtype=float)

    def _sigma(self, x, sigma2in, terms, cidx, covblk, r=None, clampT=None):
        """
        One Lyapunov solve. Sigma is LINEAR in itself for a held x, so this is a
        SOLVE and not an iteration -- which is why Sigma stays out of the Newton
        unknowns and only the M numbers it projects onto go in.

        R IS THE FIRING RATE VECTOR when a cap is active, so the diffusion matrix
        D diag(r) D' counts the events that actually happen. The drift Jacobian is
        NOT re-derived for the multiplier: it is constant within a mode, so it
        rescales the admission columns, and treating it as constant here is a
        first-order approximation of A rather than the exact one. It is flagged
        rather than hidden because the Gaussian closure over a capacity-truncated
        support is the larger approximation of the two.

        CLAMPT IS THE TANGENT SPACE OF THE CAPS THAT CLAMP. A cap that holds the
        job upstream or loses it fixes the constrained combination A x at b for as
        long as it binds -- blocking answers the state instantaneously -- so that
        combination does not fluctuate and the noise lives on the subspace
        orthogonal to A. WITHOUT THIS THE LYAPUNOV SOLVE HAS NO SOLUTION AT ALL,
        and not by accident: with the multiplier held constant the drift along the
        constrained direction is neutral (an overloaded M/M/1/K sits on a whole
        line of equilibria), so the Jacobian carries a zero eigenvalue there and
        FluidNonHyperbolicError is the correct verdict for the UNPROJECTED system.
        Projecting the jump directions is enough to state the reduced problem,
        because fluid_lyapunov restricts everything to range(D) already.

        A STAGED cap is NOT clamped and must not be projected: there the drain
        rate is constant and the room population is a state, so the region
        population genuinely fluctuates and its restoring force is the room.
        """
        A = terms.jac(x, sigma2in, covblk)
        if r is None:
            r = terms.rates(x, sigma2in, covblk)
        idx = terms.covIdx
        Dc = terms.D[idx, :]
        if clampT is not None:
            Dc = clampT @ Dc
        Qc = Dc @ np.diag(r) @ Dc.T
        Sc = fluid_lyapunov(A[np.ix_(idx, idx)], Qc, Dc)
        Sigma = np.zeros((terms.nstate, terms.nstate))
        Sigma[np.ix_(idx, idx)] = Sc
        return self._sigma2_from(Sigma, terms, cidx), Sigma

    @staticmethod
    def _clamp_tangent(con, active, terms):
        """
        Orthogonal projector onto the subspace the CLAMPING caps leave free.

        One row per active cap that holds or loses, restricted to the covariance
        coordinates, and the projector that annihilates them: ``I - R'(RR')^-1 R``.
        None when no active cap clamps, which is every model without a station
        buffer and every region model -- so the uncapped and the staged paths form
        the same Lyapunov system they always did.
        """
        rows = [int(c) for c in active if not con.staged[int(c)]]
        if not rows:
            return None
        idx = terms.covIdx
        R = con.A[np.ix_(rows, idx)]
        if R.size == 0 or not np.any(np.abs(R) > _ZERO):
            return None
        return np.eye(len(idx)) - R.T @ np.linalg.pinv(R @ R.T) @ R

    @staticmethod
    def _sigma2_from(Sigma, terms, cidx):
        sigma2 = np.zeros(terms.M)
        for i in cidx:
            blk = terms.stationBlock[i]
            if blk.size:
                sigma2[i] = max(0.0, float(np.sum(Sigma[np.ix_(blk, blk)])))
        return sigma2

    def _legs(self, terms, gates, stg, staged, active, x, sg, r, mult, staged_flow):
        """
        The two legs of every event under the active caps, and the waiting rooms.

        Shared by the steady-state residual and the transient right-hand side so
        that the two solve the SAME model and not two spellings of it. What differs
        between them is only what a STAGED cap's multiplier means:

          staged_flow=False (steady state) -- a drain RATE theta. The room mass is
              then pinned by the cap (s = inflow/theta at the fixed point), which
              is what makes the algebraic system square, and the proportional split
              across a region's rooms falls out of the common rate.
          staged_flow=True (transient) -- the admitted FLOW itself. A rate cannot
              start the constrained phase at all: at the instant the region fills
              the room is EMPTY, so theta*s is zero however large theta is, and
              holding the cap needs a finite admitted flow immediately. The flow is
              split across the room masses, or across their inflows while the rooms
              are still empty -- and the two rules AGREE at a fixed point, where
              s_j is proportional to inflow_j, so the transient and the steady
              state describe one model.

        A held or lost cap composes as a product of fractions either way.

        Returns (rup, rin, ds, drain): the rate each event FIRES at, the rate mass
        LANDS at, the derivative of each room, and each room's total outflow.
        """
        nevents = int(np.size(terms.rateBase))
        nstg = stg.n
        whole = np.ones(nevents)
        entry = np.ones(nevents)
        inv_theta = np.zeros(nstg)
        flow_of = np.zeros(nstg)
        throttled = np.zeros(nstg, dtype=bool)
        rooms_of = []
        for k in range(len(active)):
            c = int(active[k])
            m = float(mult[k])
            if staged[c]:
                rooms = np.where(stg.gated_by[c, :])[0]
                rooms_of.append((c, m, rooms))
                throttled[rooms] = True
                if not staged_flow:
                    inv_theta[rooms] += (1.0 / m) if m > _ZERO else np.inf
            else:
                whole[gates.held[c, :]] *= m
                entry[gates.loss[c, :]] *= m

        adm_split = np.zeros(nevents, dtype=bool)
        if nstg:
            for e in np.where(stg.adm)[0]:
                if throttled[stg.adm_stage[e]]:
                    adm_split[e] = True
        idx_split = np.where(adm_split)[0]

        # A staged event is NOT suppressed upstream even when a held cap gates it:
        # the job completes upstream service into the room, so the held fraction
        # moves to the room's exit leg below.
        rup = r * np.where(adm_split, 1.0, whole)
        rin = rup * entry

        inflow = np.zeros(nstg)
        R = np.zeros(nstg)
        for e in idx_split:
            inflow[stg.adm_stage[e]] += rup[e]
            R[stg.adm_stage[e]] += rup[e]

        if staged_flow:
            # each row's admitted flow, split across its rooms by mass while any
            # is held and by inflow while they are all empty
            for c, m, rooms in rooms_of:
                if rooms.size == 0:
                    continue
                mass = float(np.sum(sg[rooms]))
                if mass > _FINE:
                    w = sg[rooms] / mass
                else:
                    tot = float(np.sum(inflow[rooms]))
                    w = inflow[rooms] / tot if tot > _FINE else np.full(rooms.size, 1.0 / rooms.size)
                flow_of[rooms] += m * w
            drain = flow_of.copy()
        else:
            pos = inv_theta > 0
            theta_of = np.zeros(nstg)
            theta_of[pos] = 1.0 / inv_theta[pos]
            flow_of = theta_of * sg
            drain = flow_of.copy()

        left = np.zeros(nstg)
        for e in idx_split:
            j = stg.adm_stage[e]
            q = flow_of[j] * rup[e] / R[j] if R[j] > _FINE else 0.0
            # the held fraction gates the room's EXIT: what it stops stays in the
            # room, which is what a waiting queue does with it. The lost fraction
            # gates the ARRIVAL: that mass leaves the room and is destroyed, so it
            # is drained but never delivered.
            rin[e] = q * whole[e] * entry[e]
            left[j] += q * whole[e]
        hot = R > _FINE
        drain[hot] = left[hot]

        ds = inflow - drain
        # a waiting room with no cap above it must be EMPTY, not merely balanced
        ds[~throttled] = sg[~throttled]
        return rup, rin, ds, drain

    def _residual(self, u, ctx, quiet):
        """
        The coupled algebraic system, stacked. Returns ``(None, None, None)`` when
        the closure cannot be evaluated at this iterate, so that the line search
        can back off; the caller evaluates the seed with ``quiet=False``, where a
        genuine failure must surface.

        THE UNKNOWNS ARE ``[x; staging; sigma2; mult]``. Staging carries the
        blocked mass a capped region will not admit yet; MULT is one multiplier
        per ACTIVE cap, and which multiplier it is depends on where that cap sends
        the mass it stops:

          held or lost cap -> a FRACTION in [0,1] of the admissions it allows.
                              Several caps gating one event compose as a PRODUCT,
                              which is what independent blocking gives and what
                              keeps every active cap present in the Jacobian: a
                              product has a live derivative in each factor.
          staged cap       -> the RATE its waiting room drains at. Several caps
                              gating one room compose HARMONICALLY, 1/theta =
                              sum 1/theta_f, because the waits a job serves in
                              turn ADD. One cap is then exactly theta_f, which is
                              what the single-throttle region route computed.

        The two legs of an admission are tracked separately: RUP is the rate the
        event actually FIRES at, RIN the rate at which mass reaches the
        destination. They differ only under a cap -- lost mass has no
        destination, staged mass rests in the room first -- so where nothing is
        active this is exactly ``D r``.
        """
        terms = ctx['terms']
        nstate, cidx, covblk, stg = ctx['nstate'], ctx['cidx'], ctx['covblk'], ctx['stg']
        active, M = ctx['active'], ctx['M']
        gates, staged = ctx['gates'], ctx['conStaged']
        clampT = ctx.get('clampT')
        nstg, nact, ncl = stg.n, active.size, cidx.size
        nevents = int(np.size(terms.rateBase))

        x = u[:nstate]
        sg = u[nstate:nstate + nstg]
        s2 = np.zeros(M)
        # The Newton unknowns are unconstrained but a variance, a population and
        # a multiplier are not. `_newton` projects the iterate into the feasible
        # box, so these clamps are guards rather than the mechanism; clamping
        # alone would flatten the Jacobian at the boundary and pin the unknown.
        if ncl:
            s2[cidx] = np.maximum(0.0, u[nstate + nstg:nstate + nstg + ncl])
        mult = np.maximum(0.0, u[nstate + nstg + ncl:])

        if quiet:
            try:
                r = terms.rates(x, s2, covblk)
            except Exception:
                return None, None, None
        else:
            r = terms.rates(x, s2, covblk)

        # THE FIRING RATE, not the nominal one: a cap that holds the job upstream
        # suppresses the event itself, so the upstream station keeps the mass and
        # reports it -- which is what LINE's own semantics does with a closed job
        # that finds no room (State.arrivalIsLost returns an empty successor and
        # the departure does not fire).
        rup, rin, ds, _ = self._legs(terms, gates, stg, staged, active, x, sg, r,
                                     mult, staged_flow=False)

        if quiet:
            try:
                s2new, _ = self._sigma(x, s2, terms, cidx, covblk, rup, clampT)
            except Exception:
                return None, None, None
        else:
            s2new, _ = self._sigma(x, s2, terms, cidx, covblk, rup, clampT)

        if gates.Dn is None:
            drift = terms.D @ rup
        else:
            # THE THREE LEGS. The removal leaves at the rate the event FIRES and
            # the arrival lands at the rate mass actually ARRIVES; Dn + Dp = D, so
            # this collapses to D @ rup wherever the two agree. A LOST arrival is
            # returned to the source pool rather than destroyed -- the EXT
            # coordinate is a normalisation, and its row is a real equation in the
            # MATLAB representation -- while a job lost leaving a REAL station is
            # destroyed, so that station's row keeps the firing rate.
            drift = gates.Dn @ rup + gates.DnExt @ rin + gates.Dp @ rin

        parts = [drift[ctx['driftRows']], ds]
        C = ctx['C']
        if C.size:
            parts.append(C @ np.concatenate([x, sg]) - ctx['Nvec'])
        if ncl:
            parts.append(s2[cidx] - s2new[cidx])
        if nact:
            # THE UNDIFFERENTIATED CONSTRAINT, on purpose. At a fixed point every
            # flow already balances, so the differentiated form A D r = 0 is
            # satisfied by anything and pins no multiplier. A x = b does pin it,
            # through the dependence of the fixed point on the multiplier.
            parts.append(ctx['conA'][active, :] @ x + ctx['conAs'][active, :] @ sg
                         - ctx['conb'][active])
        return np.concatenate(parts), rin, rup

    def _newton(self, ctx, u, tol, maxit, nfree):
        """
        Damped PROJECTED Newton with a finite-difference Jacobian and an Armijo
        backtrack on the residual norm. The step is solved in LEAST SQUARES: the
        drift block is rank deficient by exactly the number of conserved chains,
        and the constraint rows restore that rank, so the stacked system is
        consistent and overdetermined rather than square.

        WHY PROJECTED, AND NOT MERELY CLAMPED. The unknowns past NFREE are a
        variance and an admission throttle, neither of which may go negative, and
        the residual reads them through max(0,.). Clamping inside the residual
        ALONE is a trap: once an iterate goes negative the residual stops
        depending on it, so the finite-difference column is exactly zero, there
        is no derivative to climb back on, and the unknown is pinned at the
        boundary for good. Projecting the ITERATE keeps every evaluation inside
        the box, where the forward difference across max(0,.) is live even at
        zero.
        """
        u = self._project(u, nfree)
        G, _, _ = self._residual(u, ctx, quiet=False)
        resnorm = float(np.max(np.abs(G))) if G.size else 0.0
        converged = resnorm < tol
        it = 0
        while not converged and it < maxit:
            it += 1
            J = self._fdjac(ctx, u, G)
            try:
                du = -np.linalg.lstsq(J, G, rcond=None)[0]
            except np.linalg.LinAlgError:
                break
            if not np.all(np.isfinite(du)):
                du = -np.linalg.pinv(J) @ G
            lam = 1.0
            stepped = False
            for _ls in range(25):
                # project BEFORE evaluating, so the accepted point and the
                # residual that measured it are the same feasible point
                un = self._project(u + lam * du, nfree)
                Gn, _, _ = self._residual(un, ctx, quiet=True)
                if Gn is not None and np.all(np.isfinite(Gn)) \
                        and float(np.max(np.abs(Gn))) < resnorm * (1 - 1e-4 * lam):
                    u, G = un, Gn
                    resnorm = float(np.max(np.abs(Gn)))
                    stepped = True
                    break
                lam *= 0.5
            if not stepped:
                break     # no descent along this direction; report the iterate
            converged = resnorm < tol
        return u, it, resnorm < tol, resnorm

    @staticmethod
    def _project(u, nfree):
        """Onto the feasible box: the state is free, the variances are not."""
        u = np.asarray(u, dtype=float).copy()
        if u.size > nfree:
            u[nfree:] = np.maximum(0.0, u[nfree:])
        return u

    def _fdjac(self, ctx, u, G):
        n, m = u.size, G.size
        J = np.zeros((m, n))
        for j in range(n):
            h = max(1e-7 * abs(u[j]), 1e-9)
            up = u.copy()
            up[j] += h
            Gp, _, _ = self._residual(up, ctx, quiet=True)
            if Gp is None or not np.all(np.isfinite(Gp)):
                um = u.copy()
                um[j] -= h
                Gm, _, _ = self._residual(um, ctx, quiet=True)
                if Gm is None or not np.all(np.isfinite(Gm)):
                    continue   # column left at zero; the least-squares step absorbs it
                J[:, j] = (G - Gm) / h
            else:
                J[:, j] = (Gp - G) / h
        return J

    # -------------------------------------------------------------- transient

    def _hold_multipliers(self, terms, gates, stg, con, active, x, sg, s2, m0):
        """
        The multipliers that hold the active caps at this state, by small Newton.

        Each active cap contributes one equation, ``d/dt (A x + As s) = 0``: the
        constrained quantity is AT the cap, so holding it there means the flow
        across the cap balances. Each contributes one unknown too -- a fraction for
        a cap that holds or loses, an admitted flow for one that stages -- so the
        system is square and small (one entry per binding cap, never more than a
        handful), and a finite-difference Newton on it costs a few drift
        evaluations.

        This is what makes an event RESTART consistent: RODAS and ode15s both need
        the algebraic unknowns to satisfy their equations at the initial point of
        an index-1 DAE, and the integrator's own Newton cannot be asked for them
        before the first step. It is also the FEASIBILITY test at an activation --
        a fraction above one, or an admitted flow above what arrives, means the cap
        is not binding after all and must not be activated.
        """
        active = np.asarray(active, dtype=int)
        if active.size == 0:
            return np.zeros(0), True

        def resid(m):
            r = terms.rates(x, s2, [None] * terms.M)
            rup, rin, ds, _ = self._legs(terms, gates, stg, con.staged, active, x, sg,
                                         r, m, staged_flow=True)
            dx = (gates.Dn @ rup + gates.DnExt @ rin + gates.Dp @ rin) \
                if gates.Dn is not None else terms.D @ rup
            return np.array([float(con.A[c, :] @ dx + con.As[c, :] @ ds) for c in active])

        m = np.asarray(m0, dtype=float).copy()
        F = resid(m)
        for _ in range(40):
            if float(np.max(np.abs(F))) < max(1e-12, 1e-10 * float(np.max(np.abs(m)) + 1.0)):
                break
            J = np.zeros((F.size, m.size))
            for j in range(m.size):
                h = max(1e-7 * abs(m[j]), 1e-9)
                mp = m.copy(); mp[j] += h
                J[:, j] = (resid(mp) - F) / h
            step = np.linalg.lstsq(J, -F, rcond=None)[0]
            lam = 1.0
            for _ls in range(20):
                mn = np.maximum(0.0, m + lam * step)
                Fn = resid(mn)
                if float(np.max(np.abs(Fn))) < float(np.max(np.abs(F))):
                    m, F = mn, Fn
                    break
                lam *= 0.5
            else:
                break
        return m, bool(float(np.max(np.abs(F))) < 1e-6)

    def _transient(self, terms, C, Nvec, cidx, covblk, sigma2ss, timespan, tol,
                   xfall, tfall, con=None, gates=None, stg=None):
        """
        The transient closure as an index-1 DAE, integrated with a SINGULAR mass
        matrix by the vendored RODAS.

        RODAS IS NAMED RATHER THAN DISPATCHED, and it has to be: the ordinary fluid
        path integrates y' = f with LSODA, which cannot carry a mass matrix at all,
        and scipy's BDF/Radau take no M either. The MATLAB reference reaches ode15s
        here and can reach RODAS through ``options.odesolvers.daeSolver``; this has
        only the one route, which is also what makes the two agree in the last
        digits.

        UNDER A CAP THIS IS A HYBRID DAE, and it is integrated as one rather than
        refused. The system SWITCHES every time a cap starts or stops binding, so
        the horizon is covered by segments: each segment is one index-1 DAE with a
        fixed set of binding caps, the segment ends at a located crossing, and the
        next one starts from that state with the set updated. Integrating with the
        set frozen instead would silently report the unconstrained trajectory
        through a cap the model declares, which is why this was refused before it
        was implemented.

        THE MULTIPLIER IS AN ALGEBRAIC UNKNOWN HERE, one per cap, carried in the
        state vector with a zero mass row, and its equation is the DIFFERENTIATED
        constraint ``d/dt (A x + As s) = 0``. The undifferentiated form ``A x = b``
        would be an index-2 DAE -- the multiplier appears only after one
        differentiation -- and neither RODAS nor ode15s solves index 2. A cap that
        is not binding keeps its unknown pinned at the inert value (one for a
        fraction, zero for a flow) by an equation of its own, so the state layout
        is the same on every segment and a restart costs no re-indexing.

        A STAGED cap's unknown is the admitted FLOW and not the drain rate the
        steady state solves for; see _LEGS for why a rate cannot start a
        constrained phase.
        """
        nstate = terms.nstate
        idx = terms.covIdx
        nc = int(np.size(idx))
        cfg = getattr(self.options, 'config', None) or {}
        maxcov = int(cfg.get('dae_maxcov') or 25)
        ncon = 0 if con is None else int(con.b.size)
        nstg = 0 if stg is None else int(stg.n)

        # The covariance adds nc^2 differential states and the Jacobian is
        # formed by finite differences over all of them, so this grows as nc^4.
        # Above the cap the mean is still integrated as a DAE -- conservation
        # stays an equation -- but the variance is held at its stationary value,
        # which is what `minnormal` does for the whole transient anyway.
        withcov = 0 < nc <= maxcov
        if nc > maxcov:
            warnings.warn(
                "The transient covariance would add %d differential states, above the limit of "
                "%d set by options.config['dae_maxcov']. Integrating the mean as a DAE with the "
                "variance held at its stationary value; raise that limit for a time-varying "
                "second moment." % (nc * nc, maxcov * maxcov), RuntimeWarning)

        # The initial condition is taken from the seed trajectory rather than
        # from options.init_sol so that it is certainly in the closing state
        # layout that TERMS indexes; the two agree whenever every class carries a
        # rate, and this does not depend on that.
        x0 = np.asarray(xfall, dtype=float)[0, :].copy()
        if x0.size != nstate:
            return tfall, xfall, None, None, []

        # One differential equation per closed chain is redundant -- the rows of
        # D sum to zero there -- so one is REPLACED by the constraint rather than
        # added to it. The row dropped is the one carrying the most mass at t=0,
        # which keeps the algebraic equation away from a coordinate that is
        # identically zero. Chains partition the classes, so no index is claimed
        # twice.
        algrow = []
        for rr in range(C.shape[0]):
            sup = np.where(C[rr, :nstate] != 0)[0]
            if sup.size == 0:
                continue
            algrow.append(int(sup[int(np.argmax(x0[sup]))]))
        algrow = np.asarray(algrow, dtype=int)

        ncov = nc * nc if withcov else 0
        nz = nstate + nstg + ncon + ncov
        o_sg, o_m, o_cov = nstate, nstate + nstg, nstate + nstg + ncon
        Dc = terms.D[idx, :]
        inert = np.ones(ncon)
        if ncon:
            inert[con.staged[:ncon]] = 0.0   # a flow is inert at zero, a fraction at one

        # A STATE ABOVE A CAP IS NOT A STATE THE MODEL CAN BE IN. The initial point
        # comes from the caller -- the default for a closed model spreads the
        # population over the stations -- so it may sit outside a cap; holding the
        # cap from there would freeze the violation for the whole horizon, since
        # the equation says the constrained quantity does not MOVE. Start on the
        # cap instead, and say so.
        #
        # THE EXCESS IS MOVED, NOT DROPPED. Conservation is an algebraic row here,
        # so an initial state that does not satisfy it is an INCONSISTENT
        # initialisation and no index-1 solver may be handed one. The mass goes
        # where the model would have put it: into the waiting room of a staged cap,
        # and otherwise onto the coordinates that feed the capped stations, which
        # is where a held job waits.
        sg0 = np.zeros(nstg)
        excess = np.zeros(max(ncon, 1))
        over = np.zeros(0, dtype=int)
        if ncon:
            over = np.where((con.A @ x0) > con.b + max(1e-9, tol))[0]
            for c in over:
                val = float(con.A[c, :] @ x0)
                if not (val > con.b[c] > 0):
                    continue
                sel = con.A[c, :] > 0
                excess[c] = float(np.sum(x0[sel])) * (1.0 - con.b[c] / val)
                x0[sel] *= con.b[c] / val

        def unpack(z):
            x = z[:nstate]
            sg = z[o_sg:o_m] if nstg else np.zeros(0)
            m = z[o_m:o_cov] if ncon else np.zeros(0)
            if withcov:
                Sc = z[o_cov:].reshape(nc, nc)
                Sc = 0.5 * (Sc + Sc.T)
                Sigma = np.zeros((nstate, nstate))
                Sigma[np.ix_(idx, idx)] = Sc
                s2 = self._sigma2_from(Sigma, terms, cidx)
            else:
                Sc, s2 = None, sigma2ss
            return x, sg, m, Sc, s2

        def derivs(x, sg, m, s2, active):
            r = terms.rates(x, s2, covblk)
            if ncon:
                rup, rin, ds, _ = self._legs(terms, gates, stg, con.staged, active,
                                             x, sg, r, m[active], staged_flow=True)
                dx = (gates.Dn @ rup + gates.DnExt @ rin + gates.Dp @ rin) \
                    if gates.Dn is not None else terms.D @ rup
            else:
                rup, rin, ds = r, r, np.zeros(0)
                dx = terms.D @ r
            return dx, ds, rup

        def make_rhs(active, gated_rooms):
            def rhs(_t, z, out):
                x, sg, m, Sc, s2 = unpack(z)
                dx, ds, rup = derivs(x, sg, m, s2, active)
                out[:nstate] = dx
                for rr in range(algrow.size):
                    # the algebraic rows: the mass matrix has zeroed these, so what
                    # is written here is the CONSTRAINT RESIDUAL, not a derivative
                    out[algrow[rr]] = float(C[rr, :nstate] @ x
                                            + (C[rr, nstate:] @ sg if nstg else 0.0)) - Nvec[rr]
                if nstg:
                    out[o_sg:o_m] = np.where(gated_rooms, ds, sg)
                if ncon:
                    mm = np.empty(ncon)
                    for c in range(ncon):
                        if c in active:
                            mm[c] = float(con.A[c, :] @ dx + con.As[c, :] @ ds)
                        else:
                            mm[c] = m[c] - inert[c]
                    out[o_m:o_cov] = mm
                if withcov:
                    A = terms.jac(x, s2, covblk)
                    Ac = A[np.ix_(idx, idx)]
                    dS = Ac @ Sc + Sc @ Ac.T + Dc @ np.diag(rup) @ Dc.T
                    out[o_cov:] = dS.ravel()
            return rhs

        def cap_value(x, sg):
            if not ncon:
                return np.zeros(0)
            return con.A @ x + (con.As @ sg if nstg else 0.0)

        rtol = max(tol, 1e-10)
        atol = max(tol, 1e-12)
        t0 = float(timespan[0]) if timespan[0] is not None else 0.0
        t1 = float(timespan[1])

        z0 = np.zeros(nz)
        # Sigma(0) = 0 is the consistent initialisation, and the physically right
        # one: the population at t=0 is a known deterministic state, so it has no
        # variance. C x0 = N holds by construction, so the algebraic rows are
        # satisfied at t=0 and no separate consistency solve is needed.

        # THE INITIAL MODE. A cap the initial state sits ON is already binding, so
        # the first segment must carry it: starting inactive would integrate one
        # step of the unconstrained system straight through the bound. Feasibility
        # decides -- a multiplier above one means the drift is pulling the
        # constrained quantity DOWN and the cap is not binding after all.
        m_init = inert.copy()
        active = set()
        if over.size:
            trial = sorted(int(c) for c in over)
            mm, ok = self._hold_multipliers(terms, gates, stg, con,
                                           np.asarray(trial, dtype=int), x0, sg0,
                                           sigma2ss, np.array([0.0 if con.staged[c] else 1.0
                                                               for c in trial]))
            keep = [c for pos, c in enumerate(trial)
                    if ok and (con.staged[c] or mm[pos] <= 1.0 + 1e-9)]
            if keep and len(keep) < len(trial):
                mm, ok = self._hold_multipliers(terms, gates, stg, con,
                                               np.asarray(keep, dtype=int), x0, sg0, sigma2ss,
                                               np.array([0.0 if con.staged[c] else 1.0
                                                         for c in keep]))
            if ok:
                active = set(keep)
                for pos, c in enumerate(keep):
                    m_init[c] = mm[pos]
        # and the mass the clip removed goes where the model would hold it: in the
        # waiting room of a staged cap that is binding, and otherwise at the
        # stations that feed the capped one, which is where a held job waits
        for c in (int(v) for v in over):
            if excess[c] <= 0:
                continue
            rooms = np.where(stg.gated_by[c, :])[0] if (nstg and con.staged[c] and c in active) \
                else np.zeros(0, dtype=int)
            if rooms.size:
                sg0[rooms] += excess[c] / rooms.size
                continue
            sel = con.A[c, :] > 0
            pool = ~sel
            if gates is not None and gates.Dn is not None:
                feeders = np.zeros(nstate, dtype=bool)
                for e in np.where(gates.gate[c, :])[0]:
                    feeders |= gates.Dn[:, e] < 0
                if np.any(feeders & pool):
                    pool = feeders & pool
            w = x0[pool]
            if float(np.sum(w)) > _FINE:
                x0[pool] += excess[c] * w / float(np.sum(w))
            elif np.any(pool):
                x0[pool] += excess[c] / float(np.count_nonzero(pool))
        if over.size:
            warnings.warn(
                "The initial state holds more jobs than %s allows; the transient starts on the "
                "cap, with the excess where the model would hold it -- in the waiting room, or at "
                "the stations feeding the cap."
                % ", ".join(con.label[int(c)] for c in over), RuntimeWarning)

        armed = set()          # a staged cap whose room has actually filled
        z0[:nstate] = x0
        if nstg:
            z0[o_sg:o_m] = sg0
            for c in active:
                if con.staged[c] and float(np.sum(sg0[stg.gated_by[c, :]])) > 1e-8:
                    armed.add(c)
        if ncon:
            z0[o_m:o_cov] = m_init
        switches = []
        ts, zs = [], []
        t_cur, z = t0, z0.copy()
        seg_max = 4 * ncon + 8
        frozen = False
        for seg in range(seg_max):
            gated_rooms = np.zeros(nstg, dtype=bool)
            for c in active:
                if con.staged[c]:
                    gated_rooms |= stg.gated_by[c, :]
            diag = np.ones(nz)
            diag[algrow] = 0.0
            if nstg:
                diag[o_sg:o_m] = gated_rooms.astype(float)
            if ncon:
                diag[o_m:o_cov] = 0.0
            rhs = make_rhs(np.asarray(sorted(active), dtype=int), gated_rooms)

            def mas(am, _d=diag):
                am[0, :] = _d

            # the event functions, all of them, as one vector: a cap that is not
            # active is watched for REACHING its bound, one that is active for
            # stopping to bind -- a fraction above one, or a room that has emptied
            def events(zz):
                x, sg, m, _Sc, _s2 = unpack(zz)
                val = cap_value(x, sg)
                g = np.ones(ncon)
                for c in range(ncon):
                    if c in active:
                        if con.staged[c]:
                            rooms = np.where(stg.gated_by[c, :])[0]
                            mass = float(np.sum(sg[rooms])) if rooms.size else 0.0
                            g[c] = mass if c in armed else 1.0
                        else:
                            g[c] = 1.0 - m[c]
                    else:
                        g[c] = con.b[c] - val[c]
                return g

            hit = {'c': -1, 't': None, 'z': None}
            gprev = events(z)
            last = {'t': t_cur, 'z': z.copy()}

            def solout(_nr, xold, xcur, y, dense, _gprev=gprev):
                yy = np.asarray(y, dtype=float)
                if not ncon:
                    ts.append(float(xcur))
                    zs.append(yy.copy())
                    return 0
                # arm a staged release only once its room has really filled, or the
                # release would fire at the activation instant, where the room is
                # empty by construction
                _x, _sg, _m, _, _ = unpack(yy)
                for c in list(active):
                    if con.staged[c] and c not in armed:
                        rooms = np.where(stg.gated_by[c, :])[0]
                        if rooms.size and float(np.sum(_sg[rooms])) > 1e-8:
                            armed.add(c)
                gcur = events(yy)
                cross = np.where((_gprev > 0) & (gcur <= 0))[0]
                _gprev[:] = gcur
                if cross.size == 0:
                    ts.append(float(xcur))
                    zs.append(yy.copy())
                    last['t'], last['z'] = float(xcur), yy.copy()
                    return 0
                # THE OVERSHOOTING STEP IS NOT REPORTED. The crossing is found at
                # the END of the step that passed it, so recording that step first
                # would put points beyond the cap in the trajectory -- measured at
                # 0.13 jobs above a cap of 5, present in the table and impossible
                # in the model. The located point is recorded instead, and the next
                # segment starts from it.
                # BISECT ON THE INTERPOLANT rather than on the integration: RODAS
                # carries a third-order interpolant over the step it has just
                # accepted, so locating the crossing costs no extra step and the
                # state at the crossing is the integrator's own, not a re-solve.
                lo, hi = float(xold), float(xcur)
                cbest = int(cross[0])
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    ymid = np.array([dense.value(i, mid) for i in range(nz)])
                    if events(ymid)[cbest] > 0:
                        lo = mid
                    else:
                        hi = mid
                    if hi - lo <= 1e-12 * max(1.0, abs(hi)):
                        break
                hit['c'] = cbest
                hit['t'] = hi
                hit['z'] = np.array([dense.value(i, hi) for i in range(nz)])
                ts.append(float(hi))
                zs.append(hit['z'].copy())
                return -1

            y = z.copy()
            try:
                res = rodas(nz, rhs, t_cur, y, t1, h=1e-6, rtol=rtol, atol=atol, itol=0,
                            ijac=0, mljac=nz, mujac=nz, ifcn=1,
                            mas=mas, imas=1, mlmas=0, mumas=0,
                            solout=solout, iout=1)
                if res.idid != 1 and hit['c'] < 0:
                    raise RuntimeError("RODAS returned idid=%d" % res.idid)
            except Exception as exc:
                warnings.warn("The transient DAE failed (%s); reporting the seed trajectory and "
                              "the stationary variance instead." % exc, RuntimeWarning)
                return tfall, xfall, None, None, switches

            if hit['c'] < 0:
                break   # the horizon was reached with this set of caps binding

            c = hit['c']
            t_cur, z = hit['t'], hit['z'].copy()
            x, sg, m, _Sc, s2 = unpack(z)
            if c in active:
                active.discard(c)
                armed.discard(c)
                z[o_m + c] = inert[c]
                switches.append((t_cur, int(c), 'release'))
            else:
                # ACTIVATING NEEDS THE MULTIPLIER THAT HOLDS THE CAP, both because
                # the restart must satisfy the algebraic rows and because that
                # multiplier IS the feasibility test: above one, the cap would have
                # to admit more than arrives, so it is not binding after all.
                trial = sorted(active | {c})
                m0 = np.array([z[o_m + k] if k in active else (0.0 if con.staged[k] else 1.0)
                               for k in trial])
                mm, ok = self._hold_multipliers(terms, gates, stg, con,
                                                np.asarray(trial, dtype=int), x, sg, s2, m0)
                feasible = ok
                for pos, k in enumerate(trial):
                    if not con.staged[k] and mm[pos] > 1.0 + 1e-9:
                        feasible = False
                if feasible:
                    active = set(trial)
                    for pos, k in enumerate(trial):
                        z[o_m + k] = mm[pos]
                    switches.append((t_cur, int(c), 'activate'))
                else:
                    # a crossing whose cap cannot be held: keep integrating and
                    # let the trajectory carry the state past it rather than
                    # freezing a violation into the algebraic rows
                    switches.append((t_cur, int(c), 'unheld'))
                    frozen = True
                    break
            if t_cur >= t1 - 1e-12:
                break
        else:
            frozen = True

        if frozen:
            warnings.warn(
                "The transient stopped switching after %d segments at t=%g of %g: the reported "
                "trajectory ends there rather than continuing with a capacity constraint that "
                "does not hold. Shorten options.timespan, or use SolverCTMC/SolverJMT/SolverSSA/"
                "SolverLDES." % (seg + 1, t_cur, t1), RuntimeWarning)

        t = np.asarray(ts, dtype=float)
        Z = np.vstack(zs) if zs else np.zeros((0, nz))
        xvec_t = Z[:, :nstate]
        if withcov:
            Sigmat = np.zeros((nstate, nstate, t.size))
            for s in range(t.size):
                Sc = Z[s, o_cov:].reshape(nc, nc)
                # the equation preserves symmetry; rounding does not
                Sc = 0.5 * (Sc + Sc.T)
                Sg = np.zeros((nstate, nstate))
                Sg[np.ix_(idx, idx)] = Sc
                Sigmat[:, :, s] = Sg
            return t, xvec_t, Sigmat, t, switches
        return t, xvec_t, None, None, switches

    @staticmethod
    def _qvar_t(Sigmat, terms):
        M, K = terms.M, terms.K
        nt = Sigmat.shape[2]
        QVart = {}
        for i in range(M):
            for k in range(K):
                blk = terms.classBlock[i][k]
                if blk.size == 0:
                    QVart[(i, k)] = np.zeros(nt)
                    continue
                v = np.zeros(nt)
                for s in range(nt):
                    v[s] = max(0.0, float(np.sum(Sigmat[np.ix_(blk, blk)][:, :, s])))
                QVart[(i, k)] = v
        return QVart

    # --------------------------------------------------------------- assembly

    def _assemble(self, terms, x, sg, sigma2, sigma2drift, covdrift, Sigma, r,
                  t, xvec_t, Sigmat, QVart, tvar, con, stg, active_solved, mult,
                  C, Nvec, iters, resnorm, converged, start_time, switches=()):
        """
        Performance measures, read back from the SAME event representation that
        defines the drift, so throughput balances flow at the fixed point. The
        THROTTLED rates, so the throughput reported at a blocked region is the
        flow that actually crosses it rather than the nominal one.
        """
        M, K = terms.M, terms.K
        gfac = terms.factors(x, sigma2drift, covdrift)

        Seff = np.asarray(terms.nservers, dtype=float).copy()
        if terms.lld is not None and np.size(terms.lld):
            lld = np.atleast_2d(terms.lld)
            for i in range(min(M, lld.shape[0])):
                Seff[i] = max(Seff[i], float(np.max(lld[i, :])))

        QN = np.zeros((M, K)); UN = np.zeros((M, K))
        TN = np.zeros((M, K)); RN = np.zeros((M, K))
        for i in range(M):
            for k in range(K):
                blk = terms.classBlock[i][k]
                if blk.size == 0:
                    continue
                QN[i, k] = float(np.sum(x[blk]))
                if terms.sched[i] == SchedStrategy.INF:
                    UN[i, k] = QN[i, k]
                else:
                    UN[i, k] = float(np.sum(gfac[blk])) / Seff[i]
                sel = terms.evIsDeparture & (terms.evStation == i) & (terms.evClass == k)
                TN[i, k] = float(np.sum(r[sel]))
            if terms.sched[i] == SchedStrategy.EXT:
                # A Source holds no jobs: its coordinates are the arrival
                # process's phase indicator, and the station-level routing folds
                # every departure of an open class back onto them.
                QN[i, :] = 0.0
        # UTILIZATION MUST BE READ FROM THE FLOW THAT ACTUALLY CROSSES once a cap
        # binds. Away from a constraint the in-service fluid sum(gfac)/s and the
        # carried utilization T/(mu s) are the same number, because the drift
        # balances; under an ACTIVE cap they are not -- the multiplier throttles
        # the departures (TN) and leaves the in-service fluid (gfac) alone, so the
        # two columns disagreed: a closed tandem capped at 1 reported Util 0.688
        # at a station whose own Tput/mu was 0.550, and the exact answer is
        # neither. Every other solver reports the CARRIED utilization (Util = X*D,
        # the utilization law), and SolverCTMC gives 0.444 for the same station,
        # so that is the convention the throttled point has to keep. Unconstrained
        # runs are bit-identical: the branch is entered only when the active set is
        # non-empty. see _kb/06-solver-catalog.md
        if np.size(active_solved) > 0:
            rates = np.asarray(getattr(self.sn, 'rates', np.zeros((M, K))), dtype=float)
            for i in range(M):
                if terms.sched[i] in (SchedStrategy.INF, SchedStrategy.EXT):
                    continue
                for k in range(K):
                    if terms.classBlock[i][k].size == 0:
                        continue
                    mu = rates[i, k] if (i < rates.shape[0] and k < rates.shape[1]) else 0.0
                    if np.isfinite(mu) and mu > 0.0:
                        UN[i, k] = TN[i, k] / (mu * Seff[i])

        # A CLASS THE MODEL NEVER ROUTES HERE HAS NO RESPONSE TIME, and neither
        # QN nor TN says so by its size: both hold a remnant of the initial state
        # the integrator was still draining. See fluid_visited_pairs.
        nz = fluid_visited_pairs(self.sn, M, K) & (TN > GlobalConstants.Zero)
        RN[nz] = QN[nz] / TN[nz]

        QVar = np.zeros((M, K))
        for i in range(M):
            for k in range(K):
                blk = terms.classBlock[i][k]
                if blk.size:
                    QVar[i, k] = max(0.0, float(np.sum(Sigma[np.ix_(blk, blk)])))

        # -- transient measures ----------------------------------------------
        # On the DAE route sigma2 varies along the trajectory, so the rate
        # factors are read at the variance the trajectory HAD at each instant.
        # `minnormal` reads its whole transient at the single stationary
        # variance, which is the quasi-steady-state approximation of this.
        QNt, UNt, TNt = {}, {}, {}
        xt = np.atleast_2d(np.asarray(xvec_t, dtype=float))
        nt = xt.shape[0]
        if nt:
            Gt = np.zeros((nt, terms.nstate))
            Rt = np.zeros((nt, int(np.size(terms.rateBase))))
            for s in range(nt):
                xs = xt[s, :]
                if Sigmat is not None and s < Sigmat.shape[2]:
                    s2s = self._sigma2_from(Sigmat[:, :, s], terms,
                                            np.arange(M))
                else:
                    s2s = sigma2drift
                Gt[s, :] = terms.factors(xs, s2s, covdrift)
                Rt[s, :] = terms.rates(xs, s2s, covdrift)
            for i in range(M):
                for k in range(K):
                    blk = terms.classBlock[i][k]
                    if blk.size == 0:
                        QNt[(i, k)] = np.zeros(nt)
                        UNt[(i, k)] = np.zeros(nt)
                        TNt[(i, k)] = np.zeros(nt)
                        continue
                    QNt[(i, k)] = np.sum(xt[:, blk], axis=1)
                    if terms.sched[i] == SchedStrategy.INF:
                        UNt[(i, k)] = QNt[(i, k)]
                    else:
                        UNt[(i, k)] = np.sum(Gt[:, blk], axis=1) / Seff[i]
                    sel = terms.evIsDeparture & (terms.evStation == i) & (terms.evClass == k)
                    TNt[(i, k)] = np.sum(Rt[:, sel], axis=1)

        CN = np.sum(RN, axis=0, keepdims=True)
        XN = self._compute_system_throughput(TN, M, K)
        from line_solver.api.sn.getters import sn_get_arvr_from_tput
        from line_solver.api.sn.transforms import sn_get_residt_from_respt
        AN = sn_get_arvr_from_tput(self.sn, TN)
        WN = sn_get_residt_from_respt(self.sn, RN, None)

        result = FLDResult(
            QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN, AN=AN, WN=WN,
            t=np.asarray(t) if t is not None else np.array([0.0]),
            xvec=x, iterations=iters,
            runtime=time.time() - start_time, method='dae')
        result.QNt, result.UNt, result.TNt = QNt, UNt, TNt
        conserr = 0.0
        if C.size:
            conserr = float(np.max(np.abs(C @ np.concatenate([x, sg]) - Nvec)))
        result.moments = {
            'Sigma': Sigma, 'QVar': QVar, 'QStd': np.sqrt(QVar),
            'sigma2': sigma2, 'sigma2Drift': sigma2drift,
            'refinement': None, 'outerIters': iters,
            'stationBlock': terms.stationBlock, 'classBlock': terms.classBlock,
            # the DAE-only outputs: a second moment along the trajectory rather
            # than at the fixed point alone
            'Sigmat': Sigmat, 'QVart': QVart, 'tvar': tvar,
            'residual': resnorm, 'converged': converged,
            'conservation': conserr,
            # What the capacity limits did: which bound, how much mass each
            # holds outside its region, and the multiplier each one settled at --
            # a fraction of the admissions allowed for a cap that holds the job
            # upstream or loses it, a drain rate for one that stages it. STAGING
            # IS REPORTED HERE AND NOT FOLDED INTO QN, because a STAGED job is at
            # no station -- the same choice LDES makes, whose station queues
            # likewise sum to less than N. A HELD job is not staged and IS folded
            # in, at the upstream station where the reference counts it.
            'capacity': {
                'label': list(con.label), 'b': con.b,
                'value': (con.A @ x + con.As @ sg) if con.b.size else np.zeros(0),
                'active': active_solved, 'staged': con.staged,
                'region': con.region, 'station': con.station,
                'staging': sg,
                'stagingRegion': stg.region, 'stagingClass': stg.klass,
                'blocked': float(np.sum(sg)) if sg.size else 0.0,
                'drain': mult, 'multiplier': mult,
                # every time the trajectory made a cap start or stop binding, as
                # (t, row, 'activate'|'release'|'unheld'); empty for a steady-state
                # solve, and the record of a hybrid transient's segments
                'switches': list(switches)}}
        return result


def solve_dae(sn, options: Optional[SolverFLDOptions] = None) -> FLDResult:
    """Entry point for ``options.method='dae'``."""
    if options is None:
        options = SolverFLDOptions()
    return DaeSolver(sn, options).solve()
