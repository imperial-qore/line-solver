"""
Perfect-sampling steady-state analysis of closed single-class product-form
networks, the sampling alternative to CTMC state-space enumeration.

States are drawn iid from the exact stationary distribution by monotone
Coupling From The Past (pfqn_cftp), so the estimator carries Monte Carlo error
O(samples^(-1/2)) but never enumerates the state space.

Reference: S. Kijima and T. Matsui, "Approximate/Perfect Samplers for Closed
Jackson Networks", Proc. Winter Simulation Conference, 2005.
"""

import numpy as np

from ....api.pfqn import pfqn_cftp
from ....api.sn import sn_get_demands_chain
from ....lang.base import NodeType, SchedStrategy, RoutingStrategy


def _cftp_sampler(method):
    """Map the solver method string onto the pfqn_cftp sampler name."""
    m = (method or '').lower()
    if m in ('cftp', 'cftp.exact'):
        return 'cftp'
    if m == 'cftp.approx':
        return 'approx'
    raise ValueError("Unknown cftp variant '%s'. Use 'cftp' or 'cftp.approx'." % method)


def solver_ctmc_cftp_supports(sn, options):
    """Can the cftp perfect sampler be asked for this model?

    The model-class gate asked as a predicate rather than raised. The analyzer
    below refuses with it, and SolverCTMC.supportsModelMethod asks the very
    same call so that a caller (model.help, findSolver, SolverAUTO) sees the
    verdict before paying for a run. A second copy of the rules in the analyzer
    is how the report and the run drift into two different answers.

    The sampler is exact only on the closed single-class product form its
    balance function encodes; anything else must be refused, not approximated.
    What the feature registry CAN name is also declared in
    SolverCTMC.getMethodFeatureSet, which usually reports the offending FEATURE
    first; this predicate is reached for the structural rules the registry has
    no name for -- the class count, the station count and the steady-state
    restriction.

    Args:
        sn: NetworkStruct of the model.
        options: solver options, read for the timespan.

    Returns:
        (ok, reason); reason is '' when ok is True.
    """
    timespan = getattr(options, 'timespan', None)
    if timespan is not None and len(timespan) > 0 and np.isfinite(timespan[0]):
        return False, ('The cftp method supports steady-state analysis only, '
                       'not transient analysis.')
    if int(sn.nclasses) != 1:
        return False, ('The cftp method supports single-class models only, '
                       'this model has %d classes.' % int(sn.nclasses))
    njobs = np.ravel(sn.njobs)
    if np.any(np.isinf(njobs)) or njobs[0] < 1:
        return False, ('The cftp method supports closed models only, with a '
                       'finite positive population.')
    if int(sn.nstations) < 2:
        return False, ('The cftp method requires at least two stations.')
    for ind in range(int(sn.nnodes)):
        if int(sn.nodetype[ind]) not in (int(NodeType.QUEUE), int(NodeType.DELAY),
                                         int(NodeType.ROUTER)):
            return False, ('The cftp method supports Queue, Delay and Router '
                           'nodes only, node %d is of a different type.' % ind)
    pf_sched = (int(SchedStrategy.INF), int(SchedStrategy.PS), int(SchedStrategy.FCFS),
                int(SchedStrategy.SIRO), int(SchedStrategy.LCFSPR))
    for i in range(int(sn.nstations)):
        if int(sn.sched[i]) not in pf_sched:
            return False, ('The cftp method requires a product-form scheduling '
                           'strategy (INF, PS, FCFS, SIRO, LCFSPR) at station %d.' % i)
        if sn.phases is not None and sn.phases[i, 0] > 1:
            return False, ('The cftp method requires exponential service times, '
                           'station %d has %d phases.' % (i, int(sn.phases[i, 0])))
        if np.isfinite(sn.cap[i]) and sn.cap[i] < njobs[0]:
            return False, ('The cftp method requires infinite buffers, station %d '
                           'has capacity %d.' % (i, int(sn.cap[i])))
        if not np.isfinite(sn.rates[i, 0]) or sn.rates[i, 0] <= 0:
            return False, ('The cftp method requires a finite positive service '
                           'rate at station %d.' % i)
    for name in ('lldscaling', 'cdscaling', 'jdscaling'):
        val = getattr(sn, name, None)
        if val is not None and np.size(val) > 0:
            return False, ('The cftp method does not support load-dependent, '
                           'class-dependent or joint-dependent service rates.')
    if getattr(sn, 'nregions', 0):
        return False, ('The cftp method does not support finite capacity regions.')
    pf_routing = (int(RoutingStrategy.PROB), int(RoutingStrategy.RAND),
                  int(RoutingStrategy.DISABLED))
    for ind in range(int(sn.nnodes)):
        for r in range(int(sn.nclasses)):
            if int(sn.routing[ind, r]) not in pf_routing:
                return False, ('The cftp method requires Markovian routing (PROB, '
                               'RAND), node %d uses a state-dependent strategy.' % ind)
    return True, ''


def solver_ctmc_cftp_analyzer(sn, options):
    """
    Steady-state metrics of a closed single-class product-form network from iid
    perfect samples of its stationary distribution.

    Args:
        sn: NetworkStruct of the model
        options: solver options carrying method, samples and timespan

    Returns:
        (QN, UN, RN, TN, CN, XN, Xs, Ts, pAggr, SSq) where Xs holds the sampled
        states, Ts the per-sample coalescence horizon, and (SSq, pAggr) the
        distinct sampled states with their empirical probabilities.
    """
    sampler = _cftp_sampler(getattr(options, 'method', 'cftp'))
    ok, reason = solver_ctmc_cftp_supports(sn, options)
    if not ok:
        raise ValueError(reason)

    M = int(sn.nstations)
    K = int(sn.nclasses)

    dem = sn_get_demands_chain(sn)
    L = np.ravel(dem.Lchain[:, 0]).astype(float)
    STchain = dem.STchain
    Vchain = dem.Vchain
    N = int(np.ravel(dem.Nchain)[0])
    iref = int(np.ravel(dem.refstatchain)[0])

    S = np.array([float(s) for s in np.ravel(sn.nservers)[:M]])
    for i in range(M):
        if int(sn.sched[i]) == int(SchedStrategy.INF):
            S[i] = np.inf

    nsamples = getattr(options, 'samples', 10000)
    if nsamples is None or not np.isfinite(nsamples) or nsamples < 1:
        raise ValueError('The cftp method requires a finite positive options.samples.')
    nsamples = int(round(nsamples))

    Qs, Xs, Ts = pfqn_cftp(L, N, S, nsamples, sampler)
    Qs = np.ravel(Qs)

    busy = np.array([np.mean(np.minimum(Xs[:, i], S[i])) for i in range(M)])

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros(K)
    XN = np.zeros(K)

    # The utilization law X = mu_i*E[min(n_i,c_i)]/V_i holds at every station, but
    # each station estimates it with its own Monte Carlo error. Taking the estimate
    # at the reference station and propagating it through the visit ratios matches
    # the CTMC convention (XN is the arrival rate at the reference station) and
    # keeps flow balance, Little's law and C = N/X exact in the reported table.
    if STchain[iref, 0] > 0 and Vchain[iref, 0] > 0:
        XN[0] = busy[iref] / STchain[iref, 0] / Vchain[iref, 0]
    if XN[0] > 0:
        CN[0] = N / XN[0]

    for i in range(M):
        QN[i, 0] = Qs[i]
        TN[i, 0] = Vchain[i, 0] * XN[0]
        # Utilization keeps its own estimator E[min(n_i,c_i)]/c_i: it is unbiased
        # and confined to [0,1] by construction, whereas deriving it from the
        # reference-station throughput lets Monte Carlo error push a saturated
        # station above 1.
        if int(sn.sched[i]) == int(SchedStrategy.INF):
            UN[i, 0] = QN[i, 0]
        else:
            UN[i, 0] = busy[i] / S[i]
        if TN[i, 0] > 0:
            RN[i, 0] = QN[i, 0] / TN[i, 0]

    SSq, counts = np.unique(Xs, axis=0, return_counts=True)
    pAggr = counts.astype(float) / nsamples

    return QN, UN, RN, TN, CN, XN, Xs, Ts, pAggr, SSq
