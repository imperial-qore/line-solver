"""
Stationary analysis by decision-diagram level aggregation.

Port of MATLAB solver_ctmc_mdd_analyzer.m. Analyses a closed single-class
network whose CTMC state space is held in a decision diagram and solved by
level aggregation, after A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact
state space of a Markov model to compute approximate stationary measures",
SIGMETRICS 2000.

This is the 'mdd' method of SolverCTMC. It never forms the |S|-state generator:
the reachable set is stored in an MDD and K coupled level-CTMCs are iterated to
a fixed point, so the memory cost is O(sum_k |M_k|) rather than O(|S|). The
saving grows with the number of stations, and is negative at K=3, where the
diagram compresses nothing.

Exactness
---------
The single approximation is Pr{i_k | alpha} = Pr{i_k | p}. It is EXACT on
product-form networks (paper Sec. 5), which covers exponential service under any
work-conserving discipline and general service at PS or IS stations (BCMP types
2 and 3). It is an approximation otherwise, notably phase-type service at FCFS
or LCFS, where errors of a fraction of a percent on the mean queue lengths have
been observed.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

import numpy as np

from ....api.io.logging import line_error, line_printf
from ....api.mdd import mdd_descriptor, mdd_mcd, mdd_ps, mdd_reachset
from .solver_ctmc_qrf_analyzer import _proc_entry_to_map

# Disciplines compared by NAME: the SchedStrategy numeric values are not shared
# across the codebases.
_SHARED = ('PS', 'DPS', 'GPS', 'INF')
_NONPREEMPTIVE = ('FCFS', 'LCFS', 'SIRO', 'HOL')


def _sched_name(s):
    if s is None:
        return ''
    if hasattr(s, 'name'):
        return str(s.name).upper()
    return str(s).upper().split('.')[-1]


def solver_ctmc_mdd_supports(sn):
    """Can the mdd decision-diagram method be asked for this model?

    The model-shape gate asked as a predicate rather than raised. The analyzer
    below refuses with it before it builds anything, and
    SolverCTMC.supportsModelMethod asks the very same call so that a caller
    (model.help, findSolver, SolverAUTO) sees the verdict without paying for a
    run. One predicate with two callers is what stops the report and the
    analyzer from disagreeing about which models the method serves.

    A STOCHASTIC PETRI NET IS EXEMPT: a Place model is read through spn_mdd,
    which builds the reachable set and the Kronecker descriptor from the
    marking rather than from the (station,class) encoding, so neither the
    single-class rule nor the closed-population rule applies to it.

    The deeper refusals the analyzer still raises -- a station-to-station chain
    that is not stochastic, and a phase-type law at a discipline neither local
    encoding represents -- are not restated here: they are decided from
    quantities the analyzer computes on its way through, not from the model
    shape, so a caller cannot be told about them without doing the work.

    Args:
        sn: NetworkStruct of the model.

    Returns:
        (ok, reason); reason is '' when ok is True.
    """
    from ...sn.network_struct import NodeType
    nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
    if np.any(nodetype == int(NodeType.PLACE)):
        return True, ''
    R = int(sn.nclasses)
    if R != 1:
        return False, ('the mdd method analyses single-class networks; this model has %d '
                       'classes. The Kronecker descriptor would need one level per '
                       '(station,class).' % R)
    if np.any(nodetype == int(NodeType.SOURCE)) or np.any(nodetype == int(NodeType.SINK)):
        return False, ('the mdd method analyses CLOSED networks; an open stream makes the '
                       'marking unbounded, so the reachable set has no finite decision '
                       'diagram')
    N = float(np.ravel(sn.njobs)[0])
    if not np.isfinite(N) or N <= 0:
        return False, 'the mdd method needs a finite positive closed population'
    return True, ''


def solver_ctmc_mdd_analyzer(sn, options=None, model=None):
    """Analyse a closed single-class network, or a stochastic Petri net, by MDD
    level aggregation.

    Passing ``model`` routes a net holding Places and Transitions through
    spn_mdd instead of mdd_descriptor: the levels are then (place, class) pairs
    plus one phase level per phase-type mode, and the measures come back per
    place. The approximation is the same Eq. 5 as for a queueing network, and it
    is exact on a product-form net, which SolverNC's 'rec' method
    (solver_nc_spn_analyzer) solves exactly and far more cheaply -- the
    aggregation earns its place on the nets that have NO product form.

    A net carries no per-station service rate, so mdd_mcd returns only the level
    marginals. The token throughput is then assembled here from the mode rates
    and those marginals, under the SAME independence across levels that the
    aggregation already assumes: it is the method's own approximation applied
    once more, not a second one layered on top.

    Returns (QN, UN, RN, TN, CN, XN, info).
    """
    verbose = bool(getattr(options, 'verbose', 0)) if options is not None else False

    # The model-shape rules live in solver_ctmc_mdd_supports, which
    # SolverCTMC.supportsModelMethod also asks: the analyzer must refuse
    # exactly what the report refuses, and one predicate with two callers is
    # what keeps the two from drifting apart.
    shape_ok, shape_reason = solver_ctmc_mdd_supports(sn)
    if not shape_ok:
        line_error('solver_ctmc_mdd_analyzer', shape_reason)

    from ...sn.network_struct import NodeType as _NT
    if np.any(np.ravel(np.asarray(sn.nodetype, dtype=int)) == int(_NT.PLACE)):
        if model is None:
            line_error('solver_ctmc_mdd_analyzer',
                       'a stochastic Petri net is read from the model object, not from the '
                       'network structure; call solver_ctmc_mdd_analyzer(sn, options, model)')
        return _spn(model, sn, options)

    M = int(sn.nstations)
    R = int(sn.nclasses)
    N = int(round(float(np.ravel(sn.njobs)[0])))

    # ---- service laws and the station-to-station routing chain
    mu = np.zeros(M)
    servers = np.zeros(M)
    proc = [None] * M
    for i in range(M):
        mu[i] = float(sn.rates[i, 0])           # 1/E[S] by LINE convention
        servers[i] = float(np.ravel(sn.nservers)[i])
        nph = 1
        phases = getattr(sn, 'phases', None)
        if phases is not None and np.size(phases) > 0:
            nph = int(np.asarray(phases).reshape(M, -1)[i, 0])
        if nph > 1:
            entry = sn.proc[i][0] if (sn.proc is not None and i < len(sn.proc)) else None
            pair = _proc_entry_to_map(entry)
            if pair is None:
                line_error('solver_ctmc_mdd_analyzer',
                           'station %d declares %d phases but carries no representable service '
                           'law' % (i + 1, nph))
            proc[i] = (np.asarray(pair[0], dtype=float), np.asarray(pair[1], dtype=float))
    P = _routing(sn, M)

    # ---- pick the encoding from the disciplines present
    sched = [sn.sched[i] if sn.sched is not None else None for i in range(M)]
    anyPH = any(p is not None for p in proc)
    kind = _encoding(sched, servers, proc, anyPH)

    if kind == 'ps':
        desc = mdd_ps(mu, P, servers, N, {'proc': proc})
    else:
        desc = mdd_descriptor(mu, P, servers, N, {'proc': proc, 'sched': sched})

    mdd = mdd_reachset(desc['domain'], desc['init'], desc['nextfun'])
    # The level iteration is an INNER numerical solve and needs a far tighter
    # tolerance than the reported means: mdd_mcd verifies the population
    # invariant at 1e-6, so a loose tolerance converges short of the fixed point
    # and trips that guard. options.iter_tol is the solver-level fixed-point
    # tolerance (sized for AMVA outer loops) and must NOT be reused here; the
    # level knobs are taken from options.config instead.
    mcdopt = {}
    cfg = getattr(options, 'config', None) if options is not None else None
    if isinstance(cfg, dict):
        if cfg.get('mdd_maxiter'):
            mcdopt['maxiter'] = int(cfg['mdd_maxiter'])
        if cfg.get('mdd_tol'):
            mcdopt['tol'] = float(cfg['mdd_tol'])
    out = mdd_mcd(mdd.to_struct(), desc, mcdopt)

    # ---- pack the analyzer contract
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    RN = np.zeros((M, R))
    TN = np.zeros((M, R))
    for i in range(M):
        QN[i, 0] = out['QLen'][i]
        UN[i, 0] = out['U'][i]
        TN[i, 0] = out['X'][i]
        if TN[i, 0] > 0:
            RN[i, 0] = QN[i, 0] / TN[i, 0]      # Little's law at the station

    # system throughput at the reference station, per unit visit
    ref = int(np.ravel(sn.refstat)[0])
    vis = 1.0
    visits = getattr(sn, 'visits', None)
    if visits is not None and len(visits) >= 1 and visits[0] is not None:
        v = np.atleast_2d(np.asarray(visits[0], dtype=float))
        # visits is indexed by STATEFUL node, refstat by station
        isf = int(sn.stationToStateful[ref]) if hasattr(sn, 'stationToStateful') else ref
        if v.shape[0] > isf and v[isf, 0] > 0:
            vis = float(v[isf, 0])
    XN = np.zeros((1, R))
    CN = np.zeros((1, R))
    XN[0, 0] = TN[ref, 0] / vis
    if XN[0, 0] > 0:
        CN[0, 0] = N / XN[0, 0]

    info = {
        'mdd': mdd,
        'desc': desc,
        'levelSizes': out['levelSizes'],
        'iters': out['iters'],
        'numStates': mdd.cardinality(),
        'encoding': kind,
        # Structural certificate: when no diagram node is shared, conditioning
        # on the node equals conditioning on the whole path and the single
        # approximation is an identity, so the answer is EXACT without a
        # reference solve. False means "not certified" rather than
        # "approximate": a product-form model is exact however much it shares.
        'noAggregation': out['noAggregation'],
        'pathsPerLevel': out['pathsPerLevel'],
    }

    if verbose:
        line_printf('\nCTMC-mdd: %d levels, |S| = %d held as %d level states (%.1fx), '
                    '%d fixed-point sweeps\n'
                    % (desc['K'], info['numStates'], int(sum(out['levelSizes'])),
                       info['numStates'] / max(int(sum(out['levelSizes'])), 1), out['iters']))
    return QN, UN, RN, TN, CN, XN, info


# ---------------------------------------------------------------------------
def _routing(sn, M):
    """Station-to-station routing probabilities of the single class.

    Read from sn.rt, which is indexed over stateful nodes, class-major.
    """
    P = np.zeros((M, M))
    R = int(sn.nclasses)
    sts = np.ravel(np.asarray(sn.stationToStateful)).astype(int)
    rt = np.asarray(sn.rt, dtype=float)
    for i in range(M):
        ni = int(sts[i])
        for j in range(M):
            nj = int(sts[j])
            P[i, j] = rt[ni * R, nj * R]
    rs = P.sum(axis=1)
    if np.any(np.abs(rs - 1.0) > 1e-8):
        line_error('solver_ctmc_mdd_analyzer',
                   'the station-to-station routing chain is not stochastic; the mdd method needs '
                   'every completion to move the job to another station')
    return P


# ---------------------------------------------------------------------------
def _encoding(sched, servers, proc, anyPH):
    """Which local-state encoding represents these disciplines exactly.

    Exponential service is discipline-insensitive for the queue-length law, so
    the compact count encoding serves any work-conserving station. Phase-type
    service is not: the count-plus-one-phase encoding is non-preemptive, while
    processor sharing needs the per-phase counts of every job present.
    """
    M = len(sched)
    if not anyPH:
        return 'np'
    isShared = np.zeros(M, dtype=bool)
    isNP = np.zeros(M, dtype=bool)
    for i in range(M):
        nm = _sched_name(sched[i])
        isShared[i] = nm in _SHARED
        isNP[i] = (nm in _NONPREEMPTIVE) and servers[i] == 1
    phAt = [i for i in range(M) if proc[i] is not None]
    if np.all(isShared):
        return 'ps'
    if all(isNP[i] for i in phAt):
        return 'np'
    bad = [i for i in phAt if not isNP[i] and not isShared[i]]
    if not bad:
        bad = [phAt[0]]
    line_error('solver_ctmc_mdd_analyzer',
               'station %d combines a phase-type service law with a discipline that neither local '
               'encoding represents: the count-plus-phase encoding is non-preemptive, and the '
               'per-phase-count encoding covers only shared servers (PS/DPS/GPS/INF). Mixing a '
               'shared and a non-preemptive phase-type station in one model is likewise '
               'unsupported.' % (bad[0] + 1))


def _spn(model, sn, options):
    """The Petri-net route: spn_mdd supplies the reachable set and the Kronecker
    descriptor, mdd_mcd aggregates, and the measures are read back per place."""
    from ...spn import spn_mdd
    from ...mdd import mdd_mcd

    verbose = int(getattr(options, 'verbose', 0) or 0)
    mdds, desc, spninfo = spn_mdd(model, {'verbose': verbose > 1})

    mcdopt = {}
    cfg = getattr(options, 'config', None)
    if cfg is not None:
        get = cfg.get if isinstance(cfg, dict) else (lambda k, d=None: getattr(cfg, k, d))
        if get('mdd_maxiter') is not None:
            mcdopt['maxiter'] = get('mdd_maxiter')
        if get('mdd_tol') is not None:
            mcdopt['tol'] = get('mdd_tol')
    out = mdd_mcd(mdds, desc, mcdopt)

    M = int(sn.nstations); R = int(sn.nclasses); L = int(spninfo['nplacelevels'])
    QN = np.zeros((M, R)); UN = np.zeros((M, R))
    RN = np.zeros((M, R)); TN = np.zeros((M, R))

    node_to_station = np.ravel(np.asarray(sn.nodeToStation)).astype(int)
    places = spninfo['places']
    for pp, node in enumerate(places):
        ist = int(node_to_station[node])
        if ist < 0:
            continue
        for k in range(R):
            QN[ist, k] = out['QLen'][pp * R + k]
            UN[ist, k] = QN[ist, k]              # a Place is an INF station: U = Q

    # Mode throughputs from the level marginals. P(m_l = v) is read off the level
    # chain; the enabling degree of a mode is then treated as independent across
    # its input levels, which is Eq. 5 of the paper applied once more rather than
    # a fresh approximation.
    pl = _level_marginals(out, mdds, L)
    md = spninfo['modes']
    for e in range(len(md)):
        if md[e]['nph'] > 1:
            continue                             # no single rate; read the phase level
        x = float(np.asarray(md[e]['D1']).ravel()[0]) * _mean_servers(pl, md[e], L)
        # FIRING EVENTS, not tokens: sn_pn_avg_rates converts the Place rows to a
        # token rate afterwards, exactly as it does for the explicit CTMC path,
        # and weighting here as well would count a weighted arc twice.
        for l in range(L):
            if md[e]['enab'][l] > 0:
                pp, k = divmod(l, R)
                ist = int(node_to_station[places[pp]])
                if ist >= 0:
                    TN[ist, k] += x
    with np.errstate(divide='ignore', invalid='ignore'):
        RN = np.where(TN > 0, QN / np.where(TN > 0, TN, 1.0), 0.0)

    XN = np.zeros(R); CN = np.zeros(R)
    refstat = np.ravel(np.asarray(sn.refstat)).astype(int)
    for k in range(R):
        ref = int(refstat[k]) if k < refstat.size else -1
        if 0 <= ref < M:
            XN[k] = TN[ref, k]
        Nk = QN[:, k].sum()
        if XN[k] > 0 and Nk > 0:
            CN[k] = Nk / XN[k]

    info = {'mdd': spninfo['mdd'], 'desc': desc, 'spn': spninfo,
            'levelSizes': out['levelSizes'], 'iters': out['iters'], 'marginal': pl}
    return QN, UN, RN, TN, CN, XN, info


def _level_marginals(out, mdds, L):
    """P(level l = v) for each place level, from the converged level chains.
    mdd_mcd works in the paper's orientation, paper level k <-> level K+1-k."""
    K = int(mdds.K if hasattr(mdds, 'K') else mdds['K'])
    domain = np.asarray(mdds.domain if hasattr(mdds, 'domain') else mdds['domain'], dtype=int)
    pl = []
    for l in range(L):
        k = K - 1 - l
        p = np.zeros(int(domain[l]))
        rows = np.asarray(out['Mrows'][k]); pk = np.asarray(out['pik'][k]).ravel()
        for r in range(rows.shape[0]):
            p[int(rows[r, 1])] += pk[r]
        tot = p.sum()
        if tot > 0:
            p = p / tot
        pl.append(p)
    return pl


def _mean_servers(pl, mde, L):
    """E[min(enabling degree, servers)] under independence across input levels."""
    lv = [l for l in range(L) if mde['enab'][l] > 0]
    if not lv:
        return 1.0
    kmax = np.inf
    for l in lv:
        kmax = min(kmax, np.floor((len(pl[l]) - 1) / mde['enab'][l]))
    if np.isfinite(mde['srv']):
        kmax = min(kmax, mde['srv'])
    n = 0.0
    for k in range(1, int(kmax) + 1):
        ge = 1.0                                 # P(deg >= k) = prod_l P(m_l >= k*I_l)
        for l in lv:
            thr = int(k * mde['enab'][l])
            if thr >= len(pl[l]):
                ge = 0.0
                break
            ge *= pl[l][thr:].sum()
        n += ge                                  # E[min(deg,srv)] = sum_k P(min >= k)
    return n
