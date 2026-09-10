"""
State-vector analyzer for SolverENV (options.method='statevec').

Carries the full per-stage state distribution across environment switches
instead of collapsing it to marginal mean queue lengths. Mirrors
matlab/src/solvers/ENV/solver_env_statevec_analyzer.m and the Java
ports. Supports a CTMC backend (explicit enumerated generator, reusing
solver_ctmc_basic + _compute_metrics_sync) and a MAM/LDQBD backend
(single-class Delay+Queue closed or Source+Queue open, exact M/M/c boundary,
PH service at any number of servers via the busy-server phase multiset of
ldqbd_mphc, load-dependent scaling, open handled by level truncation).
"""

import warnings

import numpy as np

from ...api.mc.ctmc import ctmc_timeaverage, ctmc_solve_reducible
from ...api.solvers.ctmc.handler import solver_ctmc_basic, SolverCTMCOptions
from ...api.state.marginal import toMarginal
from ...api.mam.map_analysis import map_pie, map_mean, map_cdf
from ...api.mam.ldqbd_mphc import ldqbd_mphc
from ...api.sn.predicates import sn_has_load_dependence
from ...constants import SchedStrategy


# ---------------------------------------------------------------------------
# MAM/LDQBD backend: block construction, flattening and metric mapping.
# ---------------------------------------------------------------------------
def _sched_name(sched, i):
    s = sched.get(i, None) if isinstance(sched, dict) else (sched[i] if sched is not None else None)
    if s is None:
        return None
    return s.name if hasattr(s, 'name') else str(s)


def _proc_d0d1(pr):
    """(D0, D1) of a service/arrival process given as either a (D0,D1) pair or a
    scalar-rate dict {'rate': mu} (the un-phase-expanded get_struct() form)."""
    if isinstance(pr, dict):
        mu = float(pr.get('rate'))
        return np.array([[-mu]]), np.array([[mu]])
    return np.asarray(pr[0]), (np.asarray(pr[1]) if len(pr) > 1 else None)


def ldqbd_ld(sn, cutoff=None):
    """Build the LD-QBD blocks and parameters (the ``ld`` dict) for a single-class
    Delay/Queue (closed) or Source/Queue (open) model. Mirrors solver_mam_ldqbd.m."""
    M = sn.nstations
    K = sn.nclasses
    if K != 1:
        raise RuntimeError("LDQBD method requires a single-class model.")
    njobs = np.asarray(sn.njobs).flatten()
    Npop = njobs[0]
    is_open = not np.isfinite(Npop)

    n_delay = n_queue = n_source = 0
    delay_idx = queue_idx = src_idx = -1
    for i in range(M):
        name = _sched_name(sn.sched, i)
        if name == 'INF':
            n_delay += 1; delay_idx = i
        elif name == 'FCFS':
            n_queue += 1; queue_idx = i
        elif name == 'EXT':
            n_source += 1; src_idx = i
    if is_open:
        if n_source != 1 or n_queue != 1 or M != 2:
            raise RuntimeError("Open LDQBD method requires exactly one Source and one Queue station.")
    else:
        if n_delay != 1 or n_queue != 1 or M != 2:
            raise RuntimeError("Closed LDQBD method requires exactly one Delay and one Queue station.")

    rates = np.asarray(sn.rates)
    nservers = np.asarray(sn.nservers).flatten()
    n_servers = int(nservers[queue_idx])

    D0, D1 = _proc_d0d1(sn.proc[queue_idx][0])
    is_ph = not (D0.shape[0] == 1 and D0.shape[1] == 1)
    mu = np.nan
    if not is_ph:
        mu = -D0[0, 0]
        n_phases = 1
        mean_service = 1.0 / mu
        alpha = None
    else:
        n_phases = D0.shape[0]
        alpha = np.asarray(map_pie(D0, D1)).reshape(1, -1)
        mean_service = map_mean(D0, D1)

    lld_row = None
    lldscaling = getattr(sn, 'lldscaling', None)
    has_lld = (sn_has_load_dependence(sn) and lldscaling is not None
               and np.asarray(lldscaling).size > 0
               and np.asarray(lldscaling).shape[0] > queue_idx
               and np.any(np.asarray(lldscaling)[queue_idx, :] != 1))
    if has_lld:
        lld_row = np.asarray(lldscaling)[queue_idx, :]
        lldlimit = lld_row.size
        sf_max = lld_row[lldlimit - 1]
        # Peak capacity normalizes the utilization: the LARGEST factor declared,
        # not the saturated one, matching CTMC's
        # ceff = max(nservers, max(lldscaling[ist, :])).
        util_peak = max(float(n_servers), float(np.max(lld_row)))
    else:
        lldlimit = 0
        sf_max = n_servers
        util_peak = float(n_servers)

    rt = np.asarray(sn.rt)
    if is_open:
        arr_D0, _ = _proc_d0d1(sn.proc[src_idx][0])
        if np.asarray(arr_D0).shape[0] > 1:
            raise RuntimeError("Open LDQBD method currently supports Poisson (exponential) arrivals only.")
        lam = rates[src_idx, 0]
        lambda_eff = lam * rt[src_idx, queue_idx]
        rho = lambda_eff * mean_service / sf_max
        if rho >= 1:
            raise RuntimeError("Open LDQBD method requires a stable queue (rho = %.4f >= 1)." % rho)
        if cutoff is not None and np.isscalar(cutoff) and np.isfinite(cutoff):
            Nlev = max(n_servers + 1, int(round(cutoff)))
        elif cutoff is not None and np.size(cutoff) >= 1 and np.isfinite(np.asarray(cutoff).flatten()[0]):
            Nlev = max(n_servers + 1, int(round(np.asarray(cutoff).flatten()[0])))
        else:
            tail_tol = 1e-10
            Nlev = n_servers + int(np.ceil(np.log(tail_tol) / np.log(rho)))
            Nlev = min(max(Nlev, n_servers + 10), 100000)
        arr_rate = np.full(Nlev + 1, lambda_eff)
        arr_rate[Nlev] = 0.0
    else:
        N = int(Npop)
        lambda_d = rates[delay_idx, 0]
        lambda_eff = lambda_d * rt[delay_idx, queue_idx]
        Nlev = N
        arr_rate = np.array([(N - n) * lambda_eff for n in range(Nlev + 1)])

    sf = np.zeros(Nlev)
    for n in range(1, Nlev + 1):
        sf[n - 1] = lld_row[min(n, lldlimit) - 1] if has_lld else min(n, n_servers)

    Q0, Q1, Q2 = [], [], []
    if not is_ph:
        for n in range(Nlev):
            Q0.append(np.array([[arr_rate[n]]]))
        for n in range(Nlev + 1):
            dep = sf[n - 1] * mu if n > 0 else 0.0
            Q1.append(np.array([[-(arr_rate[n] + dep)]]))
        for n in range(1, Nlev + 1):
            Q2.append(np.array([[sf[n - 1] * mu]]))
    else:
        # PH service: the level carries the multiset of the phases the min(n,c)
        # busy servers sit in, exact at any c and identical to the plain phase
        # indexing at c == 1. Same builder the steady-state LDQBD solver uses.
        Q0, Q1, Q2 = ldqbd_mphc(D0, D1, alpha, n_servers, arr_rate, sf)

    ld = {
        'Q0': Q0, 'Q1': Q1, 'Q2': Q2, 'Nlev': Nlev, 'nPhases': n_phases,
        'isPH': is_ph, 'isOpen': is_open, 'queueIdx': queue_idx, 'M': M,
        'nServers': n_servers, 'meanService': mean_service, 'hasLLD': has_lld,
        'sf': sf, 'utilPeak': util_peak,
        'lambdaEff': lambda_eff,
    }
    if is_open:
        ld['refIdx'] = src_idx; ld['delayRate'] = np.nan; ld['N'] = np.inf
    else:
        ld['refIdx'] = delay_idx; ld['delayRate'] = rates[delay_idx, 0]; ld['N'] = Npop
    return ld


def ldqbd_flatten(ld):
    """Assemble a dense generator from the block-tridiagonal representation."""
    Nlev = ld['Nlev']
    Q1 = ld['Q1']; Q0 = ld['Q0']; Q2 = ld['Q2']
    level_size = [Q1[n].shape[0] for n in range(Nlev + 1)]
    level_start = [0]
    for n in range(Nlev + 1):
        level_start.append(level_start[n] + level_size[n])
    dim = level_start[Nlev + 1]

    Q = np.zeros((dim, dim))
    level_of = np.zeros(dim, dtype=int)
    for n in range(Nlev + 1):
        r0 = level_start[n]
        level_of[r0:r0 + level_size[n]] = n
        Q[r0:r0 + level_size[n], r0:r0 + level_size[n]] = Q1[n]
        if n < Nlev:
            c0 = level_start[n + 1]
            Q[r0:r0 + level_size[n], c0:c0 + level_size[n + 1]] = Q0[n]
        if n >= 1:
            c0 = level_start[n - 1]
            Q[r0:r0 + level_size[n], c0:c0 + level_size[n - 1]] = Q2[n - 1]
    return Q, level_of


def ldqbd_avg(ld, piflat, level_of):
    """Map a flat LD-QBD distribution to per-(station,class) mean metrics."""
    Nlev = ld['Nlev']; M = ld['M']; qi = ld['queueIdx']; ri = ld['refIdx']
    pf = np.asarray(piflat).ravel().astype(float)
    pf[pf < 0] = 0
    psum = pf.sum()
    if psum > 0:
        pf = pf / psum
    pLevel = np.zeros(Nlev + 1)
    for i in range(len(pf)):
        pLevel[level_of[i]] += pf[i]
    mean_queue = float(np.dot(np.arange(Nlev + 1), pLevel))

    # Utilization is the fraction of PEAK capacity in use,
    # sum_n p(n)*sf(n)/utilPeak, the work-based convention CTMC/MVA/NC report.
    # Without load dependence sf(n) = min(n,c) and utilPeak = c, so this is the
    # average fraction of c servers in use; at c = 1 it collapses to 1 - p(0).
    sf = ld['sf']
    util_peak = ld['utilPeak']
    util = 0.0
    for n in range(1, Nlev + 1):
        util += (sf[n - 1] / util_peak) * pLevel[n]

    QN = np.zeros((M, 1)); UN = np.zeros((M, 1)); RN = np.zeros((M, 1)); TN = np.zeros((M, 1))
    if ld['isOpen']:
        X = ld['lambdaEff'] * (1 - pLevel[Nlev])
        R_queue = mean_queue / X if X > 0 else 0.0
        QN[ri, 0] = 0; UN[ri, 0] = 0; RN[ri, 0] = 0; TN[ri, 0] = X
        QN[qi, 0] = mean_queue; UN[qi, 0] = util; RN[qi, 0] = R_queue; TN[qi, 0] = X
    else:
        mean_delay = ld['N'] - mean_queue
        X = mean_delay * ld['lambdaEff']
        R_queue = mean_queue / X if X > 0 else 0.0
        R_delay = 1.0 / ld['delayRate']
        QN[ri, 0] = mean_delay; UN[ri, 0] = mean_delay; RN[ri, 0] = R_delay; TN[ri, 0] = X
        # util is already per-server: the /utilPeak is inside the sum above
        QN[qi, 0] = mean_queue; UN[qi, 0] = util; RN[qi, 0] = R_queue; TN[qi, 0] = X
    return QN, UN, RN, TN


def _station_has_lldcd(sn, ist):
    lld = getattr(sn, 'lldscaling', None)
    if lld is not None:
        lld = np.asarray(lld)
        if lld.ndim >= 2 and ist < lld.shape[0] and not np.allclose(lld[ist, :], 1.0):
            return True
    cd = getattr(sn, 'cdscaling', None)
    if cd is not None:
        cdf = (cd.get(ist) if isinstance(cd, dict) else (cd[ist] if ist < len(cd) else None))
        if callable(cdf):
            return True
    return False


def _avg_from_pi(sn, pi, SSaggr, hashed, arvRates, depRates):
    """Discipline-aware mapping of a distribution to mean metrics. Mirrors
    matlab solver_ctmc_avg_from_pi.m."""
    M = sn.nstations
    K = sn.nclasses
    pi = np.asarray(pi, dtype=float).ravel()
    # Skipped under a matrix exponential: the vector is then a genuinely signed
    # measure and the clamp would delete real mass. See ctmc handler.
    from ...api.solvers.ctmc.handler import _sn_all_phasetype
    if _sn_all_phasetype(sn):
        pi[pi < 1e-14] = 0
    s = pi.sum()
    if s > 0:
        pi = pi / s
    nservers = np.asarray(sn.nservers).flatten()

    QN = np.zeros((M, K)); UN = np.zeros((M, K)); RN = np.zeros((M, K)); TN = np.zeros((M, K))
    for ist in range(M):
        isf = int(sn.stationToStateful[ist])
        ind = int(sn.stationToNode[ist])
        S = nservers[ist]
        for k in range(K):
            TN[ist, k] = float(pi @ depRates[:, isf, k])
            QN[ist, k] = float(pi @ SSaggr[:, ist * K + k])
        name = _sched_name(sn.sched, ist)
        if name == 'EXT':
            continue
        if name == 'INF':
            for k in range(K):
                UN[ist, k] = QN[ist, k]
            continue
        is_ps = name in ('PS', 'DPS', 'GPS')
        has_sd = _station_has_lldcd(sn, ist)
        if not has_sd:
            for k in range(K):
                pr = sn.proc[ist][k]
                if pr is None or len(pr) < 2:
                    continue
                # Departure-rate estimator only; see _kb/06-solver-catalog.md
                mean = map_mean(np.asarray(pr[0]), np.asarray(pr[1]))
                UN[ist, k] = TN[ist, k] * mean / S
        else:
            # Load/class-dependent (or PAS): per-state marginal in-service share.
            # Weight the per-class capacity share by the current scaling and
            # normalize by the peak (effective capacity), as the CTMC analyzer
            # does: the unweighted share is a P(busy)-style value, not a
            # busy-server fraction, and would overstate a load-dependent station.
            lld = getattr(sn, 'lldscaling', None)
            lld_row = None
            ceff = S
            if lld is not None and np.asarray(lld).size > 0 and ist < np.asarray(lld).shape[0]:
                lld_row = np.asarray(lld)[ist, :]
                ceff = max(ceff, float(np.max(lld_row)))
            space_isf = (np.atleast_2d(sn.space[isf])
                         if (sn.space is not None and isf in sn.space and sn.space[isf] is not None)
                         else None)
            if space_isf is None:
                continue
            for s_idx in range(len(pi)):
                h = int(hashed[s_idx, isf])
                tm = toMarginal(sn, ind, space_isf[h:h + 1])
                ni = np.atleast_1d(tm[0]).flatten()
                nir = np.atleast_1d(tm[1]).flatten()
                sir = np.atleast_1d(tm[2]).flatten()
                if is_ps:
                    if np.any(ni <= 0):
                        continue
                    lldnow = 1.0
                    if lld_row is not None:
                        col = int(min(max(float(np.sum(ni)), 1), lld_row.size))
                        lldnow = float(lld_row[col - 1])
                    sp = np.array([float(sn.schedparam[ist, r]) for r in range(K)])
                    denom = float(np.dot(nir, sp))
                    if denom > 0:
                        for k in range(K):
                            UN[ist, k] += pi[s_idx] * nir[k] * sp[k] / denom * lldnow / ceff
                else:
                    if np.any(ni <= 0):
                        continue
                    for k in range(K):
                        if k < len(sir):
                            UN[ist, k] += pi[s_idx] * sir[k] / S
    for k in range(K):
        for ist in range(M):
            RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > 0 else 0.0
    return QN, UN, RN, TN


# ---------------------------------------------------------------------------
# Main state-vector analyzer.
# ---------------------------------------------------------------------------
def _row_normalize_nonneg(v):
    r = np.asarray(v, dtype=float).ravel().copy()
    r[r < 0] = 0
    s = r.sum()
    if s > 0:
        r = r / s
    return r


def _is_exp(dist):
    return type(dist).__name__ == 'Exp'


def _map_rep(dist):
    """(D0, D1) representation of an environment transition distribution."""
    from ...environment import _get_map_representation
    return _get_map_representation(dist)


def _opt(options, key, default):
    if isinstance(options, dict):
        return options.get(key, default)
    return getattr(options, key, default)


def _inner_opt(solver, key, default=None):
    """Read an inner-solver option from whichever holder carries it
    (.options, .solveropt, or the native solver's options)."""
    holders = [getattr(solver, 'options', None), getattr(solver, 'solveropt', None)]
    native = getattr(solver, '_native_solver', None)
    if native is not None:
        holders.append(getattr(native, 'options', None))
    for h in holders:
        if h is not None:
            v = getattr(h, key, None) if not isinstance(h, dict) else h.get(key, None)
            if v is not None:
                return v
    return default


def solver_env_statevec(env, solvers, options):
    """Run the state-vector ENV analyzer; returns (QN, UN, RN, TN) as M x K arrays."""
    env.init()
    E = env.num_stages
    ensemble = env.getEnsemble()

    # Raw env rate matrix (zero diagonal): E0rate[e,h] = rate e->h.
    E0rate = np.zeros((E, E))
    for e in range(E):
        for h in range(E):
            d = env.env[e][h]
            if d is not None and hasattr(d, 'getRate'):
                E0rate[e, h] = d.getRate()

    reset_fun = getattr(env, 'resetStateFun', None)

    def reset_state(h, e, pex):
        if reset_fun is not None and reset_fun[h][e] is not None:
            return np.asarray(reset_fun[h][e](pex)).ravel()
        return np.asarray(pex).ravel()

    sojourn_det = (_opt(options, 'sojourn', None) is not None
                   and str(_opt(options, 'sojourn', None)).lower() == 'deterministic')

    # ---- pre: build per-stage generator + metric-mapping data ----
    stages = []
    pi_enter = [None] * E
    for e in range(E):
        solver_e = solvers[e]
        timespan = _inner_opt(solver_e, 'timespan')
        cutoff = _inner_opt(solver_e, 'cutoff')
        if timespan is None or len(timespan) < 2 or not np.isfinite(timespan[1]):
            raise RuntimeError("The statevec analyzer requires a finite inner-solver timespan "
                               "for stage %d, e.g. CTMC(model,'timespan',[0,T])." % e)
        sn_e = ensemble[e].get_struct()
        backend = type(solver_e).__name__
        st = {'timespan': list(timespan)}
        if 'MAM' in backend:
            ld = ldqbd_ld(sn_e, cutoff)
            Qf, level_of = ldqbd_flatten(ld)
            st.update(backend='mam', Q=Qf, ld=ld, levelOf=level_of)
        else:
            opts = SolverCTMCOptions()
            opts.force = True
            opts.cutoff = cutoff
            res = solver_ctmc_basic(sn_e, opts)
            st.update(backend='ctmc', Q=res.infgen, SS=res.space, SSaggr=res.space_aggr,
                      arvRates=res.arvRates, depRates=res.depRates, hashed=res.space_hashed,
                      sn=res.sn if res.sn is not None else sn_e)
        stages.append(st)

    # Warm start. A stage's OWN stationary distribution is not a usable seed
    # here: a stage that is individually unstable or critical (arrival rate >=
    # its own service rate) has no stationary law at all, and
    # ctmc_solve_reducible then returns the stationary law of the TRUNCATED
    # generator, which piles mass against the truncation wall and whose mean
    # grows linearly with the cutoff (for a critical M/M/1 truncated at N it is
    # uniform, with mean N/2).
    #
    # The fixed point chained in post() is exact -- it is the stationary
    # equation of the joint (queue,stage) chain,
    # phi_e = (sum_h phi_h q_he) (s_e I - Q_e)^-1 -- and it does contract to the
    # right answer from that seed, but the number of sweeps needed grows with
    # the cutoff. At a finite iter_max the reported result therefore drifts
    # further from the truth as the cutoff is RAISED, i.e. the natural user
    # response to a suspect number makes it worse.
    #
    # Seed instead from the environment-averaged generator sum_e probEnv[e]*Q_e,
    # which is positive recurrent exactly when the model is stable on average --
    # the regime in which the answer exists -- so its stationary law is
    # cutoff-independent. Averaging needs one common state space; when the
    # stages differ in size (reset_state is what bridges them) fall back to the
    # per-stage law, which is no worse than before.
    dims = [stages[e]['Q'].shape[0] for e in range(E)]
    pi_shared = None
    if E > 1 and all(d == dims[0] for d in dims):
        w = np.asarray(env.probEnv, dtype=float).ravel()
        if w.size != E or not np.all(np.isfinite(w)) or w.sum() <= 0:
            w = np.ones(E) / E  # stage probabilities unavailable: weight equally
        else:
            w = w / w.sum()
        Qbar = sum(w[e] * stages[e]['Q'] for e in range(E))
        pi_shared = _row_normalize_nonneg(ctmc_solve_reducible(Qbar))
    for e in range(E):
        if pi_shared is None:
            pi_enter[e] = _row_normalize_nonneg(ctmc_solve_reducible(stages[e]['Q']))
        else:
            pi_shared = pi_shared.copy()
            pi_enter[e] = pi_shared

    pi_exit_dest = [[None] * E for _ in range(E)]
    pi_time_avg = [None] * E
    pi_enter_prev = list(pi_enter)

    # ---- analyze stage e ----
    def analyze(e):
        st = stages[e]
        Q = st['Q']
        pi0 = pi_enter[e]
        t0, t1 = st['timespan'][0], st['timespan'][1]

        if sojourn_det:
            d_e = max(map_mean(env.holdTime[e][0], env.holdTime[e][1]), np.finfo(float).eps)
            pi_avg, pi_ex = ctmc_timeaverage(pi0, Q, d_e)
            for h in range(E):
                pi_exit_dest[e][h] = pi_ex if E0rate[e, h] > 0 else None
            pi_time_avg[e] = pi_avg
            return

        exp_sojourn = all(not (E0rate[e, h] > 0 and not _is_exp(env.env[e][h])) for h in range(E))
        if exp_sojourn:
            s_e = float(np.sum(E0rate[e, :]))
            d = Q.shape[0]
            A = s_e * np.eye(d) - Q
            pi_res = np.linalg.solve(A.T, s_e * pi0)
            for h in range(E):
                pi_exit_dest[e][h] = pi_res if E0rate[e, h] > 0 else None
            pi_time_avg[e] = pi_res
            return

        # General Markovian sojourn: adaptive transient pi(t)=pi0*exp(Q t).
        from scipy.integrate import solve_ivp
        sol = solve_ivp(lambda t, y: y @ Q, [t0, t1], pi0, method='RK23',
                        rtol=1e-3, atol=1e-6, dense_output=False)
        tgrid = sol.t
        pit = sol.y.T
        for h in range(E):
            if E0rate[e, h] <= 0:
                pi_exit_dest[e][h] = None
                continue
            D0, D1 = _map_rep(env.env[e][h])
            w = _cdf_increment_weights(D0, D1, tgrid)
            pi_exit_dest[e][h] = _weighted_avg_rows(pit, w)
        D0h, D1h = env.holdTime[e][0], env.holdTime[e][1]
        wh = _cdf_increment_weights(D0h, D1h, tgrid)
        avg = _weighted_avg_rows(pit, wh)
        pi_time_avg[e] = avg if avg is not None else pit[-1, :]

    # ---- post: chain entry distributions ----
    def post():
        nonlocal pi_enter, pi_enter_prev
        pi_enter_prev = [None if p is None else p.copy() for p in pi_enter]
        new = [None] * E
        for e in range(E):
            ns = stages[e]['Q'].shape[0]
            acc = np.zeros(ns)
            wsum = 0.0
            for h in range(E):
                po = env.probOrig[h, e]
                if po > 0 and pi_exit_dest[h][e] is not None:
                    pex = reset_state(h, e, pi_exit_dest[h][e])
                    if pex.size != ns:
                        raise RuntimeError("resetStateFun[%d][%d] returned a %d-element vector but stage %d "
                                           "has %d states." % (h, e, pex.size, e, ns))
                    acc = acc + po * pex
                    wsum += po
            if wsum > 0:
                acc = acc / wsum
            else:
                acc = pi_enter[e].copy()
            new[e] = _row_normalize_nonneg(acc)
        pi_enter = new

    def converged():
        if pi_enter_prev is None:
            return False
        l1 = 0.0
        for e in range(E):
            a, b = pi_enter[e], pi_enter_prev[e]
            if a is None or b is None or a.size != b.size:
                return False
            l1 = max(l1, float(np.sum(np.abs(a - b))))
        if not np.isfinite(l1):
            return False
        return l1 < _opt(options, 'iter_tol', 1e-4)

    # ---- iterate ----
    max_iter = _opt(options, 'iter_max', 100)
    has_converged = False
    for _ in range(int(max_iter)):
        for e in range(E):
            analyze(e)
        post()
        if converged():
            has_converged = True
            break
    if not has_converged:
        # Exiting on iter_max is a NON-convergence: the last iterate can be far
        # from the fixed point (and, before the warm start above was fixed,
        # was reliably so). Returning it silently is what makes a wrong number
        # indistinguishable from a right one.
        warnings.warn(
            "The SolverENV statevec fixed point did not converge in "
            "options.iter_max=%d iterations; the returned solution is the last "
            "iterate and may be far from the fixed point. Raise options.iter_max, "
            "or loosen options.iter_tol only if the residual is already small."
            % int(max_iter), RuntimeWarning, stacklevel=2)

    # ---- finish: environment-averaged blend ----
    M = ensemble[0].get_struct().nstations
    K = ensemble[0].get_struct().nclasses
    Qval = np.zeros((M, K)); Uval = np.zeros((M, K)); Tval = np.zeros((M, K))
    for e in range(E):
        piF = pi_time_avg[e]
        if piF is None:
            continue
        st = stages[e]
        if st['backend'] == 'mam':
            QN, UN, RN, TN = ldqbd_avg(st['ld'], piF, st['levelOf'])
        else:
            QN, UN, RN, TN = _avg_from_pi(st['sn'], piF, st['SSaggr'], st['hashed'],
                                          st['arvRates'], st['depRates'])
        pe = env.probEnv[e]
        Qval += pe * QN
        Uval += pe * UN
        Tval += pe * TN

    Rval = np.zeros((M, K))
    for i in range(M):
        for k in range(K):
            Rval[i, k] = Qval[i, k] / Tval[i, k] if Tval[i, k] > 0 else 0.0
    return Qval, Uval, Rval, Tval


def _cdf_increment_weights(D0, D1, t):
    nt = len(t)
    w = np.zeros(nt)
    if nt < 2:
        return w
    later = map_cdf(D0, D1, np.asarray(t[1:]))
    earlier = map_cdf(D0, D1, np.asarray(t[:-1]))
    w[1:] = np.asarray(later).ravel() - np.asarray(earlier).ravel()
    return w


def _weighted_avg_rows(pit, w):
    sw = float(np.sum(w))
    if np.any(np.isnan(w)) or not (sw > 0):
        return None
    return (w @ pit) / sw
