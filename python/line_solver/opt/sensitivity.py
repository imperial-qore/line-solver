"""
Analytic performance sensitivities for line-opt.

Computes exact derivatives of per-(station, class) performance measures with
respect to continuous model parameters (service rates), for product-form
networks, so the optimizer can assemble analytic gradients without extra
solver evaluations.

Supports open product-form networks (single-server queueing and infinite-server
delay stations, which decouple) and closed single-server product-form networks
with unit visit ratios (via the exact differentiated MVA primitive pfqn_sens,
full cross-station Jacobian). Mixed, multiserver, and non-unit-visit topologies
return ``None`` (the optimizer then falls back to finite differences).

References:
    Z. Liu and P. Nain, INRIA RR-1144, 1989 (mixed BCMP, Thm 3.2 for open).
    X.-R. Cao and D.-J. Ma, Performance Evaluation 26:181-199, 1996.
"""

import logging
from typing import Any, Dict, Optional

import numpy as np

logger = logging.getLogger(__name__)


def compute_model_sensitivities(model: Any, use_ctmc: bool = False) -> Optional[Dict]:
    """Analytic d(metric)/d(rate) for a product-form model, keyed by name.

    Parameters
    ----------
    model : Network
        A solved-or-solvable LINE network (the per-evaluation copy).
    use_ctmc : bool, optional
        Enable the generator-derivative CTMC fallback (default ``False``).
        The open and closed product-form branches are fast but narrow: they
        return ``None`` for open, multiserver, or non-unit-visit networks.
        When they do and ``use_ctmc`` is true, the exact generator-derivative
        fallback ``_ctmc_sensitivities`` is tried instead, which is limited by
        the state-space size rather than by product-form assumptions. Mirrors
        MATLAB ``opt.sens.computeModelSensitivities``.

    Returns
    -------
    dict or None
        ``{metric_kind: {metric_key: {param_key: value}}}`` where metric_kind
        is 'RespT', 'QLen', 'Tput' or 'Util'; metric_key is (station_name,
        class_name) except for 'Util' which is keyed by station_name; and
        param_key is ('rate', station_name, class_name). ``None`` when the
        network is not a supported open product-form model. The CTMC fallback
        produces only 'QLen'.
    """
    try:
        s = _open_sensitivities(model)
        if s is not None:
            return s
        s = _closed_sensitivities(model)
    except Exception as e:  # never let sensitivity failure break evaluation
        logger.debug("analytic sensitivity unavailable: %s", e)
        s = None

    if s is None and use_ctmc:
        try:
            s = _ctmc_sensitivities(model)
        except Exception as e:
            logger.debug("CTMC sensitivity unavailable: %s", e)
            s = None
    return s


# SchedStrategy ids (see line_solver.constants). EXT marks a Source; INF an
# infinite-server delay station. Others are treated as single-server FCFS/PS.
_SCHED_EXT = 16
_SCHED_INF = 14


def _open_sensitivities(model: Any) -> Optional[Dict]:
    sn = model.getStruct()
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    # supported only when every chain is open (infinite population)
    if njobs.size == 0 or not np.all(np.isinf(njobs)):
        return None

    nstations = int(sn.nstations)
    R = int(sn.nclasses)
    rates = np.asarray(sn.rates, dtype=float)            # (nstations x R)
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    sched = getattr(sn, 'sched', {})
    nodenames = list(getattr(sn, 'nodenames', []))
    classnames = list(getattr(sn, 'classnames', []))

    # station index -> node index (for names). Fall back to identity.
    def station_node_name(ist):
        try:
            nidx = int(sn.stationToNode(ist))
        except Exception:
            nidx = ist
        if 0 <= nidx < len(nodenames):
            return str(nodenames[nidx])
        return str(ist)

    def sched_of(ist):
        s = sched.get(ist) if isinstance(sched, dict) else None
        return int(getattr(s, 'value', s)) if s is not None else -1

    # visits per (station, class): sn.visits is a dict chain -> (nstations x R)
    visits = np.zeros((nstations, R))
    v = getattr(sn, 'visits', None)
    if isinstance(v, dict):
        for _, vm in v.items():
            visits += np.asarray(vm, dtype=float)
    else:
        visits[:] = 1.0

    # arrival rate per class from the Source/EXT station(s)
    lam = np.zeros(R)
    source_stations = [i for i in range(nstations) if sched_of(i) == _SCHED_EXT]
    for i in source_stations:
        for r in range(R):
            if np.isfinite(rates[i, r]):
                lam[r] += rates[i, r]

    # only single-server queueing (or delay) stations are supported here
    queue_stations = []
    for i in range(nstations):
        sc = sched_of(i)
        if sc == _SCHED_EXT:
            continue
        if sc != _SCHED_INF and nservers[i] > 1:
            return None            # multiserver M/M/c: Erlang-C, not yet here
        queue_stations.append(i)

    sens = {'RespT': {}, 'QLen': {}, 'Tput': {}, 'Util': {}}

    def classname(r):
        return str(classnames[r]) if r < len(classnames) else str(r)

    for i in queue_stations:
        st_name = station_node_name(i)
        is_delay = (sched_of(i) == _SCHED_INF)
        # per-class demand and load at station i
        D = np.zeros(R)
        rho = np.zeros(R)
        for r in range(R):
            mu = rates[i, r]
            if np.isfinite(mu) and mu > 0 and visits[i, r] > 0:
                D[r] = visits[i, r] / mu
                rho[r] = lam[r] * D[r]
        U = 0.0 if is_delay else float(np.sum(rho))
        denom = 1.0 - U                          # M/M/1 idle probability
        if not is_delay and denom <= 0:
            return None                          # unstable: sensitivity blows up

        # station utilization U_i = sum_s rho_s; d U_i/d mu_s = -rho_s/mu_s
        util_entry = {}
        for s in range(R):
            mu_s = rates[i, s]
            if np.isfinite(mu_s) and mu_s > 0 and visits[i, s] > 0:
                util_entry[('rate', st_name, classname(s))] = -rho[s] / mu_s
        if util_entry:
            sens['Util'][st_name] = util_entry

        for r in range(R):
            mu = rates[i, r]
            if visits[i, r] <= 0 or not (np.isfinite(mu) and mu > 0):
                continue
            cl = classname(r)
            dResp = {}
            dQ = {}
            for s in range(R):
                mu_s = rates[i, s]
                if not (np.isfinite(mu_s) and mu_s > 0 and visits[i, s] > 0):
                    continue
                pkey_s = ('rate', st_name, classname(s))
                dU = 0.0 if is_delay else (-rho[s] / mu_s)     # d U_i/d mu_s
                dD_r = (-D[r] / mu_s) if s == r else 0.0        # d D_r/d mu_s
                drho_r = (-rho[r] / mu_s) if s == r else 0.0    # d rho_r/d mu_s
                if is_delay:
                    dResp[pkey_s] = dD_r
                    dQ[pkey_s] = drho_r
                else:
                    # RespT_r = D_r/(1-U), QLen_r = rho_r/(1-U): quotient rule
                    dResp[pkey_s] = (dD_r * denom + D[r] * dU) / (denom * denom)
                    dQ[pkey_s] = (drho_r * denom + rho[r] * dU) / (denom * denom)
            sens['RespT'][(st_name, cl)] = dResp
            sens['QLen'][(st_name, cl)] = dQ
            sens['Tput'][(st_name, cl)] = {('rate', st_name, cl): 0.0}

    return sens


def _closed_sensitivities(model: Any) -> Optional[Dict]:
    """Analytic d(metric)/d(rate) for a closed single-server product-form model.

    Full cross-station Jacobian (closed networks couple stations), computed via
    the exact differentiated MVA primitive pfqn_sens. Assumes unit visit ratios
    per station (serial/tandem topologies); returns None otherwise-unsupported.
    Util is aggregated per station to match EvaluationResult.utilizations.
    """
    from ..api.sn.transforms import sn_get_product_form_params
    from ..api.pfqn.sens import pfqn_sens
    from ..api.sn.network_struct import NodeType

    sn = model.getStruct()
    N = np.asarray(sn.njobs, dtype=float).flatten()
    if N.size == 0 or not np.all(np.isfinite(N)):
        return None                                   # closed only
    pf = sn_get_product_form_params(sn)
    D = np.atleast_2d(np.asarray(pf.D, dtype=float))
    S = np.asarray(pf.S, dtype=float).flatten()
    if np.any(S[np.isfinite(S)] > 1):
        return None                                   # single-server only
    R = int(sn.nclasses)
    Z = np.atleast_2d(np.asarray(pf.Z, dtype=float)).sum(axis=0)
    sens = pfqn_sens(D, N, Z)

    node_to_station = np.asarray(sn.nodeToStation).flatten()
    nodenames = list(sn.nodenames)
    classnames = list(sn.classnames)
    queue_nodes = [i for i, nt in enumerate(sn.nodetype) if nt == NodeType.QUEUE]
    rates = np.asarray(sn.rates, dtype=float)

    # precompute rate parameter (queue j, class s) -> (param index p, chain, key)
    params = []
    for j, node_j in enumerate(queue_nodes):
        sj = int(node_to_station[node_j])
        for s in range(R):
            ratej = rates[sj, s]
            if not np.isfinite(ratej) or ratej <= 0 or D[j, s] <= 0:
                continue
            params.append((j, s, j * R + s, -D[j, s] / ratej,
                           ('rate', str(nodenames[node_j]), str(classnames[s]))))

    out = {'RespT': {}, 'QLen': {}, 'Tput': {}, 'Util': {}}
    for i, node_i in enumerate(queue_nodes):
        name_i = str(nodenames[node_i])
        util_entry = out['Util'].setdefault(name_i, {})
        for r in range(R):
            if D[i, r] <= 0:
                continue
            key = (name_i, str(classnames[r]))
            dR = {}; dQ = {}; dT = {}
            for (j, s, p, chain, pkey) in params:
                dR[pkey] = float(sens.dR[i, r, p] * chain)
                dQ[pkey] = float(sens.dQ[i, r, p] * chain)
                dT[pkey] = float(sens.dX[r, p] * chain)        # unit visits
                util_entry[pkey] = util_entry.get(pkey, 0.0) \
                    + float(sens.dU[i, r, p] * chain)
            out['RespT'][key] = dR
            out['QLen'][key] = dQ
            out['Tput'][key] = dT
    return out


def _ctmc_sensitivities(model: Any) -> Optional[Dict]:
    """Analytic d(QLen)/d(rate) for an arbitrary Markovian model, via the
    generator-derivative equation of Trivedi and Bobbio (2017), Eq. (9.81).

    This is the fallback for the cases the differentiated-MVA primitive cannot
    reach: open, multiserver, and non-unit-visit networks, for which
    ``_open_sensitivities`` and ``_closed_sensitivities`` return ``None``. It is
    exact wherever SolverCTMC is exact, and correspondingly it is limited by the
    state-space size rather than by the product-form assumptions. Mirrors MATLAB
    ``opt.sens.ctmcSensitivities``.

    Only queue-length sensitivities are produced. Utilization, response time,
    and throughput are reward rates whose definition involves the solved metrics
    themselves, so they need the reward-derivative term of Eq. (9.83) and are not
    covered here. Returns ``None`` when the model has no CTMC representation of
    tractable size or no perturbable exponential rate.
    """
    from ..constants import ProcessType
    from ..distributions.continuous import Exp
    from ..solvers.solver_ctmc.solver_ctmc import SolverCTMC

    sn = model.getStruct()
    M = int(sn.nstations)
    K = int(sn.nclasses)
    rates = np.asarray(sn.rates, dtype=float)
    nodenames = list(sn.nodenames)
    classnames = list(sn.classnames)
    station_to_node = np.asarray(sn.stationToNode).ravel()
    procid = getattr(sn, 'procid', None)

    def station_node_idx(ist):
        return int(station_to_node[ist])

    solver = SolverCTMC(model)
    try:
        space_aggr = np.asarray(solver.getStateSpaceAggr(), dtype=float)
    except Exception:
        return None                          # state space unavailable / over the gate
    if space_aggr.size == 0:
        return None
    space_aggr = np.atleast_2d(space_aggr)

    # One parameter per finite positive exponential service rate. The rate setter
    # substitutes an Exp of the perturbed rate, which is a perturbation of theta
    # only where the nominal process is itself exponential; for an Erlang or
    # Coxian service the substitution would change the distribution family and the
    # difference quotient would not be dQ/dtheta, so those stations are skipped.
    params = []
    for ist in range(M):
        for k in range(K):
            rate = rates[ist, k]
            if not np.isfinite(rate) or rate <= 0:
                continue
            if procid is None or procid[ist, k] != ProcessType.EXP:
                continue
            node_idx = station_node_idx(ist)
            params.append({
                'value': float(rate),
                'node': node_idx,
                'class': k,
                'pkey': ('rate', str(nodenames[node_idx]), str(classnames[k])),
            })

    out: Dict = {'QLen': {}}
    for pr in params:
        param = {'value': pr['value'],
                 'set': _make_rate_setter(pr['node'], pr['class'], Exp)}
        try:
            _, _, dpi, _ = solver.getSensitivity(param)
        except Exception:
            continue                         # perturbation changed the state space; skip
        dpi = np.asarray(dpi, dtype=float).flatten()

        for ist in range(M):
            name_i = str(nodenames[station_node_idx(ist)])
            for r in range(K):
                col = ist * K + r
                if col >= space_aggr.shape[1]:
                    continue
                rvec = space_aggr[:, col]
                if not np.all(np.isfinite(rvec)):
                    continue                 # Source stations carry an infinite population
                mkey = (name_i, str(classnames[r]))
                out['QLen'].setdefault(mkey, {})[pr['pkey']] = float(dpi @ rvec)

    return out if out['QLen'] else None


def _make_rate_setter(node_idx: int, k: int, Exp):
    """Handle setting the service rate of node ``node_idx`` class ``k`` on a
    model copy. The distribution is replaced by an exponential of the requested
    rate, so this is only meaningful where the nominal service is itself
    exponential. Node and class are resolved inside the copy, since the objects
    of the original model do not belong to it.
    """
    def _set(m, value):
        nodes = m.get_nodes()
        classes = m.get_classes()
        nodes[node_idx].set_service(classes[k], Exp(value))
    return _set
