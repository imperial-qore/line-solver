"""
RCAT-based INAP solver for SolverMAM.

Implements the RCAT (Reversed Compound Agent Theorem) iterative fixed-point
solver (INAP and INAP+) matching the MATLAB solver_mam_ag.m implementation.

Each (station, class) pair becomes an isolated component, and the components
are coupled only through the reversed rates of the synchronizing actions. A
component is a QBD whose LEVEL is the queue length and whose PHASE is the pair
(arrival phase, service phase), laid out in the Kronecker order of
qbd_mapmap1: an arrival moves the level up carrying kron(D1^a, I), a service
completion moves it down carrying kron(I, D1^s), the busy levels evolve under
krons(D0^a, D0^s) and level zero under kron(D0^a, I), because no server is
running there. With exponential processes every block is 1x1 and the QBD
collapses to the scalar birth-death chain this analyzer used before.

Algorithm:
1. Build RCAT model from network structure
2. Initialize action rates randomly
3. Iterate:
   - Compute equilibrium distribution for each process using CTMC solve
   - Update action rates based on equilibrium terms:
     - INAP: x(a) = mean over the support of (pi Aa)_j / pi_j
     - INAP+: x(a) = sum_ij Aa(i,j) pi(i)
   - Check convergence
4. Compute performance metrics from equilibrium distributions

References:
    MATLAB: matlab/src/solvers/MAM/solver_mam_ag.m
"""

import numpy as np
import time
from typing import Tuple, Optional, List, Dict, Any
from dataclasses import dataclass

from . import AGAlgorithm, AGResult
from ...solver_mam.utils.network_adapter import extract_mam_params
from ....api.mc import ctmc_solve, ctmc_makeinfgen
from ....constants import ProcessType


# The process types the RCAT construction can give a phase dimension to. After
# sn_nonmarkov_toph (which the RCAT methods run with phfit='ph' and
# preserveDet=False) each of these holds a genuine (D0,D1) pair with
# non-negative off-diagonal rates and a single arrival per epoch.
#
# The list is an ALLOW-list on purpose: a process type nobody has checked
# against this construction must be refused, not answered. Refused here are the
# laws whose matrices are not a generator (ME, RAP), those that are not
# time-homogeneous (NHPP, MAPt, PHt), those that are not continuous-time
# (DMAP), and those that arrive in batches (BMAP, MMAP), since a batch moves
# the level by more than one.
RCAT_PROCESS_TYPES = frozenset([
    'EXP', 'ERLANG', 'HYPEREXP', 'PH', 'APH', 'COXIAN', 'COX2', 'MAP',
    'MMPP2', 'DET', 'UNIFORM', 'GAMMA', 'PARETO', 'WEIBULL', 'LOGNORMAL',
    'REPLAYER', 'TRACE', 'IMMEDIATE', 'DISABLED',
])


def rcat_supports_processes(sn, method='inap') -> Tuple[bool, Optional[str]]:
    """Reject processes the RCAT phase construction cannot represent.

    build_rcat gives every component a service-phase and an arrival-phase
    dimension, so any law with a genuine (D0,D1) Markovian representation is
    admissible. What is not is a law whose matrices are not a generator, one
    that is not time-homogeneous, one that is not continuous-time, and one that
    arrives in batches. See RCAT_PROCESS_TYPES above.

    Compare the ProcessType BY NAME: the enum's numeric values differ across
    codebases (MATLAB EXP=0, Python EXP=1).
    """
    procid = getattr(sn, 'procid', None)
    if procid is None:
        return True, None
    rates = np.asarray(getattr(sn, 'rates', None), dtype=float)
    # A signal is a trigger with no service, so its service entry is never
    # read; only its Source arrival rate is, and that one must stay exponential
    # because the removal is folded into the component as a scalar rate.
    issignal = np.asarray(getattr(sn, 'issignal', []), dtype=float).ravel()
    nodetype = getattr(sn, 'nodetype', None)
    station_to_node = getattr(sn, 'stationToNode', None)

    def _is_source(ist):
        # NodeType.SOURCE is value 0, the same convention _build_rcat uses.
        if nodetype is None or station_to_node is None:
            return False
        try:
            ntype = np.asarray(nodetype).ravel()[
                int(np.asarray(station_to_node).ravel()[ist])]
        except (IndexError, TypeError, ValueError):
            return False
        if hasattr(ntype, 'value'):
            ntype = ntype.value
        elif hasattr(ntype, 'ID'):
            ntype = ntype.ID
        try:
            return int(ntype) == 0
        except (TypeError, ValueError):
            return False

    for ist in range(sn.nstations):
        is_src = _is_source(ist)
        for r in range(sn.nclasses):
            is_sig = r < len(issignal) and bool(issignal[r])
            if is_sig and not is_src:
                continue
            try:
                p = procid[ist][r] if not hasattr(procid, 'shape') else procid[ist, r]
            except (IndexError, TypeError, KeyError):
                continue
            if p is None:
                continue
            pname = p.name if hasattr(p, 'name') else str(p)
            # Only a process that is actually in use can mis-answer.
            if rates.size and (not np.isfinite(rates[ist, r]) or rates[ist, r] <= 0):
                continue
            if is_sig:
                if pname == ProcessType.EXP.name:
                    continue
                return False, (
                    "The %s method needs an exponential signal arrival process "
                    "(a removal signal is folded into the component as a scalar "
                    "rate), but station %d class %d is %s. Use the dec.source "
                    "method for such models." % (method, ist + 1, r + 1, pname)
                )
            if pname in RCAT_PROCESS_TYPES:
                continue
            return False, (
                "The %s method supports processes with a Markovian (D0,D1) "
                "representation only (RCAT builds a phase dimension per "
                "component out of it), but station %d class %d is %s. Use the "
                "dec.source method for such models." % (method, ist + 1, r + 1, pname)
            )

    # see _kb/06-solver-catalog.md (MAM: "AG / RCAT methods") -- RCAT is
    # single-server-only, gated since 2026-07-17
    nservers = np.asarray(getattr(sn, 'nservers', []), dtype=float).ravel()
    for ist in range(min(sn.nstations, nservers.size)):
        c = nservers[ist]
        if np.isfinite(c) and c > 1:
            return False, (
                "The %s method supports single-server stations only (RCAT does "
                "not model sn.nservers, so a multiserver station is driven at "
                "rho = lambda/mu instead of lambda/(c*mu)), but station %d has "
                "%d servers. Use the dec.source method for multiserver models."
                % (method, ist + 1, int(c))
            )
    return True, None


def _krons(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Kronecker sum, MATLAB's krons: kron(A, I_nb) + kron(I_na, B)."""
    A = np.atleast_2d(np.asarray(A, dtype=float))
    B = np.atleast_2d(np.asarray(B, dtype=float))
    return np.kron(A, np.eye(B.shape[0])) + np.kron(np.eye(A.shape[0]), B)


def _blk(n: int, m: int) -> slice:
    """Row/column range of level N (0-based) in a component with M phases."""
    return slice(n * m, n * m + m)


def _is_markovian_map(D0, D1) -> bool:
    """True when (D0,D1) is a genuine MAP: D0 has non-negative off-diagonal
    rates, D1 is non-negative, and (D0+D1) has zero row sums. A RAP or ME
    violates the sign conditions while still defining a valid point process, so
    a CTMC assembled from it is a rational generator whose stationary solution
    is a signed vector."""
    if D0 is None or D1 is None:
        return False
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    if D0.ndim != 2 or D0.shape[0] != D0.shape[1] or D1.shape != D0.shape:
        return False
    if not (np.all(np.isfinite(D0)) and np.all(np.isfinite(D1))):
        return False
    tol = 1e-9 * max(1.0, float(np.max(np.abs(np.concatenate([D0.ravel(), D1.ravel()])))))
    off = D0 - np.diag(np.diag(D0))
    return (bool(np.all(off >= -tol)) and bool(np.all(D1 >= -tol))
            and bool(np.all(np.abs((D0 + D1).sum(axis=1)) <= tol)))


def _proc_map(sn, ist: int, r: int, rates: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """(D0,D1) of the process at (IST,R), or the exponential pair built from
    sn.rates when the struct carries no usable matrix representation.

    A non-Markovian pair (a RAP or an ME) would assemble a rational generator
    rather than a CTMC, so it is refused here and answered as its mean rate;
    the solver-level gate rejects those models before they reach this point.
    """
    D0 = D1 = None
    proc = getattr(sn, 'proc', None)
    if proc is not None:
        try:
            entry = proc[ist][r]
        except (IndexError, TypeError, KeyError):
            entry = None
        if entry is not None:
            from ....api.sn.proc_form import proc_to_map
            try:
                c0, c1 = proc_to_map(entry)
            except Exception:
                c0, c1 = None, None
            if _is_markovian_map(c0, c1):
                D0, D1 = np.asarray(c0, dtype=float), np.asarray(c1, dtype=float)
    if D0 is None:
        rate = rates[ist, r] if (ist < rates.shape[0] and r < rates.shape[1]) else 0.0
        if not np.isfinite(rate) or rate <= 0:
            rate = 0.0
        D0 = np.array([[-rate]])
        D1 = np.array([[rate]])
    return D0, D1


@dataclass
class RCATModelData:
    """RCAT model data structure.

    Attributes:
        R: Rate matrices dictionary
        AP: Action-process mapping (num_actions x 2)
        process_map: Mapping from (station, class) to process index
        action_map: List of action details
        N: State space sizes for each process (nlev * mph)
        num_processes: Number of processes
        num_actions: Number of actions
        nlev: Number of QBD levels per process
        mph: Number of phases per level per process
        level: Level index of every state, per process
        svcrate: Service completion rate out of every state, per process
                 (zero on level 0, where no server runs)
        svcdown: Service completion rate per phase at a busy level, per process
    """
    R: Dict[Tuple[int, int], np.ndarray]
    AP: np.ndarray
    process_map: np.ndarray
    action_map: List[Dict]
    N: np.ndarray
    num_processes: int
    num_actions: int
    nlev: np.ndarray = None
    mph: np.ndarray = None
    level: List[np.ndarray] = None
    svcrate: List[np.ndarray] = None
    svcdown: List[np.ndarray] = None


def _build_rcat(sn, max_states: int = 100) -> RCATModelData:
    """Build RCAT model from LINE network structure.

    Args:
        sn: NetworkStruct
        max_states: Maximum number of states for truncation

    Returns:
        RCATModelData instance
    """
    # maxStates is a STATE COUNT, and the callers below read it straight off
    # options.config, which does not have to hand one over as an int. The
    # lang='python' bridge is the case that bites: MATLAB has no integer scalar,
    # so options.config.maxStates crosses py.* as a float, an open class then
    # gets nlev = 100.0, and np.tile(svcrow, nlev - 1) raises "'float' object
    # cannot be interpreted as an integer" from inside the truncation. The JAR
    # (int maxStates) and C++ (std::size_t max_states) are typed and never see
    # it; Python is the one codebase where the declared type has to be enforced.
    max_states = int(max_states)
    M = sn.nstations
    K = sn.nclasses

    # Get routing table
    if hasattr(sn, 'rt') and sn.rt is not None:
        rt = np.asarray(sn.rt, dtype=np.float64)
    else:
        rt = np.zeros((M * K, M * K))

    # Get rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates, dtype=np.float64)
    else:
        rates = np.ones((M, K))

    # Get node types
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        nodetype = list(sn.nodetype) if not isinstance(sn.nodetype, np.ndarray) else sn.nodetype
    else:
        nodetype = [1] * M  # Default to Queue

    # Get station to node mapping
    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
        stationToNode = np.asarray(sn.stationToNode).flatten()
    else:
        stationToNode = np.arange(M)

    # Identify station types (iterate over stations)
    source_stations = []
    queue_stations = []

    for ist in range(M):
        node_idx = int(stationToNode[ist]) if ist < len(stationToNode) else ist
        if node_idx < len(nodetype):
            ntype = nodetype[node_idx]
            # Check for NodeType
            if hasattr(ntype, 'value'):
                ntype_val = ntype.value
            elif hasattr(ntype, 'ID'):
                ntype_val = ntype.ID
            else:
                ntype_val = int(ntype) if not np.isnan(ntype) else 1

            # Source = 0, Sink = 1, Queue = 2, Delay = 3
            if ntype_val == 0:
                source_stations.append(ist)
            elif ntype_val == 1:
                pass  # Sink is handled separately at node level
            else:
                queue_stations.append(ist)
        else:
            queue_stations.append(ist)

    # Identify Sink nodes (iterate over ALL nodes, not just stations)
    sink_nodes = []
    for node_idx in range(len(nodetype)):
        ntype = nodetype[node_idx]
        if hasattr(ntype, 'value'):
            ntype_val = ntype.value
        elif hasattr(ntype, 'ID'):
            ntype_val = ntype.ID
        else:
            ntype_val = int(ntype) if not np.isnan(ntype) else 1
        if ntype_val == 1:  # Sink = 1
            sink_nodes.append(node_idx)

    # Check for signal classes
    issignal = np.zeros(K, dtype=bool)
    if hasattr(sn, 'issignal') and sn.issignal is not None:
        issignal_raw = sn.issignal
        if callable(issignal_raw):
            for r in range(K):
                try:
                    issignal[r] = issignal_raw(r)
                except:
                    issignal[r] = False
        else:
            issignal = np.asarray(issignal_raw, dtype=bool).flatten()
            if len(issignal) < K:
                issignal = np.pad(issignal, (0, K - len(issignal)), constant_values=False)

    # Check for signal types (for G-networks with negative customers)
    # signaltype[r] will be True if class r is a negative signal
    signaltype_is_negative = [False] * K
    if hasattr(sn, 'signaltype') and sn.signaltype is not None:
        for r in range(K):
            if r < len(sn.signaltype):
                st = sn.signaltype[r]
                if st is not None:
                    # Check if this is a negative signal type
                    # SignalType.NEGATIVE has value 'negative' or could be enum member
                    if hasattr(st, 'value'):
                        # Enum case: check if value is 'negative' or equals NEGATIVE
                        signaltype_is_negative[r] = (st.value == 'negative' or
                                                     str(st.name).upper() == 'NEGATIVE')
                    elif hasattr(st, 'name'):
                        signaltype_is_negative[r] = str(st.name).upper() == 'NEGATIVE'
                    elif isinstance(st, str):
                        signaltype_is_negative[r] = st.lower() == 'negative'
                    else:
                        # MATLAB uses numeric SignalType - 1 = NEGATIVE
                        try:
                            signaltype_is_negative[r] = int(st) == 1
                        except:
                            pass

    # CATASTROPHE is a distinct SignalType value from NEGATIVE, so it has to be
    # tested on its own; the two arrays are read together everywhere below.
    signaltype_is_catastrophe = [False] * K
    if hasattr(sn, 'signaltype') and sn.signaltype is not None:
        for r in range(K):
            if r < len(sn.signaltype):
                st = sn.signaltype[r]
                if st is not None:
                    if hasattr(st, 'value'):
                        signaltype_is_catastrophe[r] = (st.value == 'catastrophe' or
                                                        str(st.name).upper() == 'CATASTROPHE')
                    elif isinstance(st, str):
                        signaltype_is_catastrophe[r] = st.lower() == 'catastrophe'

    is_catastrophe_src = [False] * K
    if hasattr(sn, 'iscatastrophe') and sn.iscatastrophe is not None:
        isCat = np.asarray(sn.iscatastrophe).flatten()
        for r_idx in range(min(K, len(isCat))):
            is_catastrophe_src[r_idx] = bool(isCat[r_idx])

    # Create process mapping: each (station, class) pair at queue stations
    process_idx = 0
    process_map = np.zeros((M, K), dtype=int)

    for ist in queue_stations:
        for r in range(K):
            if r < len(issignal) and issignal[r]:
                continue
            if ist < rates.shape[0] and r < rates.shape[1]:
                rate = rates[ist, r]
                if not np.isnan(rate) and rate > 0:
                    process_idx += 1
                    process_map[ist, r] = process_idx

    num_processes = process_idx

    if num_processes == 0:
        return RCATModelData(
            R={},
            AP=np.zeros((0, 2), dtype=int),
            process_map=process_map,
            action_map=[],
            N=np.array([]),
            num_processes=0,
            num_actions=0
        )

    # Identify Sink nodes was done above; the QBD shape of each component needs
    # both the service MAP of its station and the arrival MAP of the external
    # streams reaching it, so they are resolved before anything is assembled.
    N = np.zeros(num_processes, dtype=int)
    njobs = sn.njobs if hasattr(sn, 'njobs') and sn.njobs is not None else np.full(K, np.inf)
    njobs = np.asarray(njobs).flatten()

    nlev = np.zeros(num_processes, dtype=int)
    mph = np.zeros(num_processes, dtype=int)
    level = [None] * num_processes
    svcrate = [None] * num_processes
    svcdown = [None] * num_processes
    pinfo = [None] * num_processes

    for p in range(1, num_processes + 1):
        positions = np.where(process_map == p)
        ist, r = int(positions[0][0]), int(positions[1][0])
        info = {'ist': ist, 'r': r}
        info['Ds0'], info['Ds1'] = _proc_map(sn, ist, r, rates)
        (info['Da0'], info['Da1'], info['lamNeg'], info['lamCat'],
         info['batch']) = _arrival_map(sn, ist, r, rt, source_stations, K, rates,
                                       issignal, signaltype_is_negative,
                                       signaltype_is_catastrophe, is_catastrophe_src)
        info['ns'] = info['Ds0'].shape[0]
        info['na'] = info['Da0'].shape[0]
        info['mph'] = info['na'] * info['ns']
        if r < len(njobs) and njobs[r] < np.inf:
            info['nlev'] = int(njobs[r]) + 1
        else:
            info['nlev'] = max_states
        # Service completion: level down, arrival phase untouched (qbd_mapmap1).
        info['Dsvc'] = np.kron(np.eye(info['na']), info['Ds1'])
        info['N'] = info['nlev'] * info['mph']
        pinfo[p - 1] = info

        N[p - 1] = info['N']
        nlev[p - 1] = info['nlev']
        mph[p - 1] = info['mph']
        level[p - 1] = np.repeat(np.arange(info['nlev']), info['mph'])
        svcrow = info['Dsvc'].sum(axis=1)
        svcrate[p - 1] = np.concatenate(
            [np.zeros(info['mph']), np.tile(svcrow, info['nlev'] - 1)])
        svcdown[p - 1] = svcrow

    # Count actions
    action_map = []

    # Check for catastrophe flags
    is_catastrophe = [False] * K
    if hasattr(sn, 'iscatastrophe') and sn.iscatastrophe is not None:
        isCat = np.asarray(sn.iscatastrophe).flatten()
        for r_idx in range(min(K, len(isCat))):
            is_catastrophe[r_idx] = bool(isCat[r_idx])

    for ist in queue_stations:
        for r in range(K):
            if process_map[ist, r] > 0:
                # NEGATIVE and CATASTROPHE are distinct SignalType values; both tested
                is_negative_class = False
                is_catastrophe_class = False
                removal_dist = None
                if r < len(issignal) and issignal[r]:
                    if r < len(is_catastrophe):
                        is_catastrophe_class = is_catastrophe[r]
                    if r < len(signaltype_is_negative):
                        is_negative_class = signaltype_is_negative[r] or is_catastrophe_class
                    # Get removal distribution
                    if (hasattr(sn, 'signalremdist') and sn.signalremdist is not None
                            and r < len(sn.signalremdist)):
                        removal_dist = sn.signalremdist[r]

                for jst in queue_stations:
                    for s in range(K):
                        if process_map[jst, s] > 0:
                            from_idx = ist * K + r
                            to_idx = jst * K + s
                            if from_idx < rt.shape[0] and to_idx < rt.shape[1]:
                                prob = rt[from_idx, to_idx]
                                if prob > 0 and (ist != jst or r != s):
                                    action_map.append({
                                        'from_station': ist,
                                        'from_class': r,
                                        'to_station': jst,
                                        'to_class': s,
                                        'prob': prob,
                                        'isNegative': is_negative_class,
                                        'isCatastrophe': is_catastrophe_class,
                                        'removalDistribution': removal_dist,
                                    })

    num_actions = len(action_map)

    # Initialize R and AP
    R = {}
    AP = np.zeros((num_actions, 2), dtype=int) if num_actions > 0 else np.zeros((0, 2), dtype=int)

    # Build local rate matrices L for each process
    for p in range(1, num_processes + 1):
        R[(num_actions, p - 1)] = _build_local_rates(
            sn, pinfo[p - 1], rt, sink_nodes, K, rates)

    # Build active and passive matrices for each action
    for a, am in enumerate(action_map):
        ist = am['from_station']
        r = am['from_class']
        p_active = process_map[ist, r] - 1
        AP[a, 0] = p_active
        pa = pinfo[p_active]
        prob = am['prob']

        # Active matrix: level n -> n-1 carrying the service completion block
        # kron(I, D1^s), scaled by the routing probability of this action.
        Aa = np.zeros((pa['N'], pa['N']))
        for n in range(1, pa['nlev']):
            Aa[_blk(n, pa['mph']), _blk(n - 1, pa['mph'])] = pa['Dsvc'] * prob
        # see _kb/06-solver-catalog.md (MAM: "AG / RCAT methods",
        # "Boundary self-loop is closed-class only"). It is written on the
        # diagonal so it stays inert in the generator while still contributing
        # the pi(i)/pi(i) = 1 ratio the INAP estimators read off the top level.
        if r < len(njobs) and njobs[r] < np.inf:
            top = _blk(pa['nlev'] - 1, pa['mph'])
            Aa[top, top] += np.diag(pa['Dsvc'].sum(axis=1) * prob)
        R[(a, 0)] = Aa

        # Passive process (arrival or signal effect)
        jst = am['to_station']
        sc = am['to_class']
        p_passive = process_map[jst, sc] - 1
        AP[a, 1] = p_passive
        pp = pinfo[p_passive]
        Im = np.eye(pp['mph'])
        Pb = np.zeros((pp['N'], pp['N']))

        if am.get('isNegative', False):
            # NEGATIVE: Job removal at destination (G-network negative customer)
            if am.get('isCatastrophe', False):
                # CATASTROPHE: every level drops to level 0
                for n in range(pp['nlev']):
                    Pb[_blk(n, pp['mph']), _blk(0, pp['mph'])] = Im
            elif am.get('removalDistribution') is not None:
                # BATCH REMOVAL: remove a random number of jobs
                _add_batch_removal(Pb, am['removalDistribution'], 1.0,
                                   pp['nlev'], pp['mph'])
                Pb[_blk(0, pp['mph']), _blk(0, pp['mph'])] = Im
            else:
                # SINGLE REMOVAL: level n -> n-1, empty queue absorbs the signal
                Pb[_blk(0, pp['mph']), _blk(0, pp['mph'])] = Im
                for n in range(1, pp['nlev']):
                    Pb[_blk(n, pp['mph']), _blk(n - 1, pp['mph'])] = Im
        else:
            # POSITIVE: Normal job arrival at destination. The phase is
            # untouched: a job joining does not restart the server, and the
            # service phase frozen at level 0 is the one the last completion
            # left behind, which for a phase-type is already its entry
            # distribution.
            for n in range(pp['nlev'] - 1):
                Pb[_blk(n, pp['mph']), _blk(n + 1, pp['mph'])] = Im
            top = _blk(pp['nlev'] - 1, pp['mph'])
            Pb[top, top] = Im
        R[(a, 1)] = Pb

    return RCATModelData(
        R=R,
        AP=AP,
        process_map=process_map,
        action_map=action_map,
        N=N,
        num_processes=num_processes,
        num_actions=num_actions,
        nlev=nlev,
        mph=mph,
        level=level,
        svcrate=svcrate,
        svcdown=svcdown,
    )


def _arrival_map(sn, ist: int, r: int, rt: np.ndarray, source_stations: List[int],
                 K: int, rates: np.ndarray, issignal, signaltype_is_negative,
                 signaltype_is_catastrophe, is_catastrophe_src):
    """External (Source) streams reaching (IST,R), as one arrival MAP for the
    positive customers plus the scalar rates of the removal signals.

    Each stream is thinned by its routing probability -- a MAP thinned with
    probability p is (D0 + (1-p) D1, p D1) -- and the streams are superposed by
    the Kronecker sum, so several Poisson sources still collapse to the single
    rate sum this analyzer used before. Removal signals stay scalar: a signal is
    a trigger with no service, and its arrival process is required exponential.
    """
    Da0 = np.zeros((1, 1))
    Da1 = np.zeros((1, 1))
    have_arrival = False
    lam_neg = 0.0     # negative signal arrivals (single removal)
    lam_cat = 0.0     # catastrophe arrivals (remove all)
    batch = []        # batch removal arrivals: (rate, distribution)

    for isrc in source_stations:
        for s_src in range(K):
            is_signal = s_src < len(issignal) and bool(issignal[s_src])

            if is_signal:
                # A signal routes to itself, so its effect on this component is
                # its total probability of reaching this STATION in any class.
                prob_src = 0.0
                for s_dst in range(K):
                    from_idx = isrc * K + s_src
                    to_idx = ist * K + s_dst
                    if from_idx < rt.shape[0] and to_idx < rt.shape[1]:
                        prob_src += rt[from_idx, to_idx]
            else:
                from_idx = isrc * K + s_src
                to_idx = ist * K + r
                prob_src = (rt[from_idx, to_idx]
                            if from_idx < rt.shape[0] and to_idx < rt.shape[1] else 0.0)

            if prob_src <= 0:
                continue
            if not (isrc < rates.shape[0] and s_src < rates.shape[1]):
                continue
            src_rate = rates[isrc, s_src]
            if np.isnan(src_rate):
                continue

            is_removal = is_signal and (
                (s_src < len(signaltype_is_negative) and signaltype_is_negative[s_src])
                or (s_src < len(signaltype_is_catastrophe) and signaltype_is_catastrophe[s_src]))
            if is_removal:
                is_cat = ((s_src < len(is_catastrophe_src) and is_catastrophe_src[s_src])
                          or (s_src < len(signaltype_is_catastrophe)
                              and signaltype_is_catastrophe[s_src]))
                if is_cat:
                    lam_cat += src_rate * prob_src
                    continue
                removal_dist = None
                if (hasattr(sn, 'signalremdist') and sn.signalremdist is not None
                        and s_src < len(sn.signalremdist)):
                    removal_dist = sn.signalremdist[s_src]
                if removal_dist is not None:
                    batch.append((src_rate * prob_src, removal_dist))
                else:
                    lam_neg += src_rate * prob_src
                continue

            if src_rate <= 0:
                continue
            S0, S1 = _proc_map(sn, isrc, s_src, rates)
            if prob_src < 1:
                S0 = S0 + (1.0 - prob_src) * S1
                S1 = prob_src * S1
            if have_arrival:
                Da0 = _krons(Da0, S0)
                Da1 = _krons(Da1, S1)
            else:
                Da0, Da1 = S0, S1
                have_arrival = True

    return Da0, Da1, lam_neg, lam_cat, batch


def _add_batch_removal(B: np.ndarray, dist, rate: float, nlev: int, mph: int) -> np.ndarray:
    """Accumulate the level n -> level m block of a batch removal into B,
    scaled by RATE. Landing on the empty level absorbs the whole upper tail of
    the pmf, which is what keeps the block stochastic once the batch exceeds
    the queue length. The phase is untouched: a removal takes a waiting job,
    not the one in service."""
    Im = np.eye(mph)
    for n in range(1, nlev):
        for m in range(1, n + 1):
            pk = dist.evalPMF(n - m) if hasattr(dist, 'evalPMF') else 0.0
            if pk > 0:
                B[_blk(n, mph), _blk(m, mph)] += rate * pk * Im
        cdf = sum(dist.evalPMF(j) if hasattr(dist, 'evalPMF') else 0.0
                  for j in range(n))
        tail = 1.0 - cdf
        if tail > 0:
            B[_blk(n, mph), _blk(0, mph)] += rate * tail * Im
    return B


def _build_local_rates(sn, p: Dict, rt: np.ndarray, sink_nodes: List[int],
                       K: int, rates: np.ndarray) -> np.ndarray:
    """Build local/hidden transition matrix for the component P.

    Note: sink_nodes contains node indices (not station indices) for Sink nodes.
    """
    ist = p['ist']
    r = p['r']
    mph = p['mph']
    nlev = p['nlev']
    Im = np.eye(mph)
    L = np.zeros((p['N'], p['N']))

    # Level-local blocks: the arrival phase always runs, the service phase only
    # while the server is busy (qbd_mapmap1's Lbar = kron(D0^a, I) at level 0
    # and L = krons(D0^a, D0^s) above it). With one phase each these are pure
    # diagonals, which ctmc_makeinfgen discards and rebuilds from the row sums.
    L[_blk(0, mph), _blk(0, mph)] = np.kron(p['Da0'], np.eye(p['ns']))
    Lbusy = _krons(p['Da0'], p['Ds0'])
    for n in range(1, nlev):
        L[_blk(n, mph), _blk(n, mph)] = Lbusy

    # Positive arrival transitions: level n -> n+1, carrying kron(D1^a, I).
    Aup = np.kron(p['Da1'], np.eye(p['ns']))
    for n in range(nlev - 1):
        L[_blk(n, mph), _blk(n + 1, mph)] += Aup
    # At the truncation the job is lost but the arrival process still moves on,
    # so the block stays on the top level. With a single arrival phase this is a
    # pure diagonal and is discarded, exactly as before.
    top = _blk(nlev - 1, mph)
    L[top, top] += Aup

    # Catastrophe arrival transitions: every busy level drops to level 0
    if p['lamCat'] > 0:
        for n in range(1, nlev):
            L[_blk(n, mph), _blk(0, mph)] += p['lamCat'] * Im

    # Batch removal arrival transitions
    for batch_rate, dist in p['batch']:
        _add_batch_removal(L, dist, batch_rate, nlev, mph)

    # Single removal negative arrivals: level n -> n-1 (busy levels only)
    if p['lamNeg'] > 0:
        for n in range(1, nlev):
            L[_blk(n, mph), _blk(n - 1, mph)] += p['lamNeg'] * Im

    # Service completions that are not synchronizing actions: departures to a
    # Sink (level down) and self-routing (level unchanged, service restarted).
    mu_ir = rates[ist, r] if (ist < rates.shape[0] and r < rates.shape[1]) else 0.0
    if not np.isnan(mu_ir) and mu_ir > 0:
        stationToNode = (sn.stationToNode if hasattr(sn, 'stationToNode')
                         and sn.stationToNode is not None else np.arange(sn.nstations))
        stationToNode = np.asarray(stationToNode).flatten()
        node_idx = int(stationToNode[ist]) if ist < len(stationToNode) else ist

        prob_sink = 0.0
        if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
            rtnodes = np.asarray(sn.rtnodes)
            for jsnk in sink_nodes:
                for s in range(K):
                    from_idx = node_idx * K + r
                    to_idx = jsnk * K + s
                    if from_idx < rtnodes.shape[0] and to_idx < rtnodes.shape[1]:
                        prob_sink += rtnodes[from_idx, to_idx]

        self_idx = ist * K + r
        prob_self = (rt[self_idx, self_idx]
                     if self_idx < rt.shape[0] and self_idx < rt.shape[1] else 0.0)

        if prob_sink > 0:
            for n in range(1, nlev):
                L[_blk(n, mph), _blk(n - 1, mph)] += p['Dsvc'] * prob_sink
        if prob_self > 0:
            for n in range(1, nlev):
                L[_blk(n, mph), _blk(n, mph)] += p['Dsvc'] * prob_self

    return L


def _birth_death_solve(Q: np.ndarray) -> np.ndarray:
    """Solve equilibrium of a birth-death (tridiagonal) CTMC using recursion.

    For a birth-death chain with birth rate lambda_n = Q[n, n+1] and
    death rate mu_n = Q[n, n-1], the equilibrium is:
        pi(n) = pi(0) * prod(lambda_k / mu_{k+1}, k=0..n-1)

    This is numerically stable and avoids the ill-conditioned linear system
    that plagues eigenvalue/null-space methods for large state spaces.
    """
    n = Q.shape[0]
    if n <= 1:
        return np.ones(1)

    pi = np.zeros(n)
    pi[0] = 1.0

    for i in range(1, n):
        birth_rate = Q[i - 1, i]
        death_rate = Q[i, i - 1]
        if death_rate > 0:
            pi[i] = pi[i - 1] * birth_rate / death_rate
        else:
            pi[i] = 0.0

    total = np.sum(pi)
    if total > 0:
        pi /= total
    else:
        pi[:] = 1.0 / n

    return pi


def _is_tridiagonal(Q: np.ndarray) -> bool:
    """Check if a matrix is tridiagonal (only main diagonal and ±1 diagonals)."""
    n = Q.shape[0]
    for i in range(n):
        for j in range(n):
            if abs(i - j) > 1 and abs(Q[i, j]) > 1e-14:
                return False
    return True


def _is_block_tridiagonal(Q: np.ndarray, lvl: np.ndarray) -> bool:
    """True when every transition of Q stays within the neighbouring level,
    LVL being the level index of each state."""
    n = Q.shape[0]
    far = np.abs(lvl[:, None] - lvl[None, :]) > 1
    return not bool(np.any(np.abs(Q[far]) > 1e-14)) if n else True


def _stat_vector(C: np.ndarray) -> np.ndarray:
    """Stationary vector of a generator C, allowing a reducible one.

    Replaces the first balance equation by the normalization, which is the
    equation it is redundant with (the columns of a generator sum to zero), and
    solves the resulting square system. Unlike a null-space solve this stays
    well posed when the chain is reducible with ONE closed class, which the
    level-0 chain of a phase-expanded component routinely is: a phase-type
    restarts in the support of alpha, so every service phase outside that
    support is unreachable once the queue has emptied at least once.
    """
    m = C.shape[0]
    A = C.copy()
    A[:, 0] = 1.0
    rhs = np.zeros(m)
    rhs[0] = 1.0
    try:
        return np.linalg.solve(A.T, rhs)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(A.T, rhs, rcond=None)[0]


def _qbd_finite_solve(Q: np.ndarray, m: int, nlev: int) -> np.ndarray:
    """Stationary vector of a finite block-tridiagonal generator.

    Linear level reduction: censor the chain level by level from the top,
        C(nlev-1) = B(nlev-1),  C(n) = B(n) + F(n) (-C(n+1))^-1 D(n+1),
    with B, F and D the diagonal, up and down blocks. C(0) is the generator of
    the chain censored on level 0, so pi_0 is its stationary vector and the rest
    follows from pi_(n+1) = pi_n F(n) (-C(n+1))^-1. This is the block form of
    _birth_death_solve and reduces to it entry for entry when m == 1.
    """
    if nlev <= 1:
        return ctmc_solve(Q)

    C = [None] * nlev
    C[nlev - 1] = Q[_blk(nlev - 1, m), _blk(nlev - 1, m)]
    for n in range(nlev - 2, -1, -1):
        F = Q[_blk(n, m), _blk(n + 1, m)]
        D = Q[_blk(n + 1, m), _blk(n, m)]
        C[n] = Q[_blk(n, m), _blk(n, m)] + F @ np.linalg.solve(-C[n + 1], D)

    pi = np.zeros(nlev * m)
    pi[_blk(0, m)] = _stat_vector(C[0])
    for n in range(nlev - 1):
        F = Q[_blk(n, m), _blk(n + 1, m)]
        pi[_blk(n + 1, m)] = np.linalg.solve((-C[n + 1]).T, (pi[_blk(n, m)] @ F).T).T

    total = pi.sum()
    if total > 0:
        pi = pi / total
    else:
        pi = np.ones(nlev * m) / (nlev * m)
    return pi


def _solve_component(Qk: np.ndarray, mph: int, nlev: int, lvl: np.ndarray) -> np.ndarray:
    """Stationary vector of one isolated component.

    A component with a single phase per level is the birth-death chain the
    analyzer has always built, and the ratio recursion is both exact and stable
    there; a phase-expanded component is block tridiagonal instead, and the
    matrix analogue of that recursion (linear level reduction) keeps the same
    stability at the 100-level truncation, where a null-space solve is already
    ill-conditioned. Anything that reaches beyond the neighbouring level -- a
    catastrophe, a batch removal -- is neither, and falls back to ctmc_solve.
    """
    if mph == 1:
        if _is_tridiagonal(Qk):
            return _birth_death_solve(Qk)
    elif _is_block_tridiagonal(Qk, lvl):
        return _qbd_finite_solve(Qk, mph, nlev)
    return ctmc_solve(Qk)


def _agent_generator(k: int, x: np.ndarray, model: RCATModelData) -> np.ndarray:
    """Agent k's generator at the current reversed rates.

    Split from the stationary solve because the halves cost different orders:
    assembling the generator is O(N^2) and solving it is O(N^3). The cluster
    backend therefore ships only the stationary vector back and rebuilds the
    generator on the coordinator, rather than putting an N-by-N payload on the
    wire to save the cheaper half.
    """
    num_actions = model.num_actions
    N = model.N
    AP = model.AP
    R = model.R

    # Start with local rates
    L = R.get((num_actions, k), np.zeros((N[k], N[k])))
    Qk = L - np.diag(L @ np.ones(N[k]))

    # Add contributions from each action
    for c in range(num_actions):
        if AP[c, 1] == k:
            # Process k is passive for action c
            Pb = R.get((c, 1), np.zeros((N[k], N[k])))
            # MATLAB: Qk = Qk + x(c) * Pb{c} - diag(Pb{c} * ones(N(k), 1));
            # Note: diagonal adjustment is NOT scaled by x[c]
            Qk = Qk + x[c] * Pb - np.diag(Pb @ np.ones(N[k]))
        elif AP[c, 0] == k:
            # Process k is active for action c
            Aa = R.get((c, 0), np.zeros((N[k], N[k])))
            Qk = Qk + Aa - np.diag(Aa @ np.ones(N[k]))

    # Convert to valid generator
    return ctmc_makeinfgen(Qk)


def _agent_solve(k: int, x: np.ndarray, model: RCATModelData) -> Tuple[np.ndarray, np.ndarray]:
    """Agent k's generator and stationary vector, given the reversed rates.

    The single point where an agent is evaluated. Every execution backend routes
    through it -- the serial loop, the thread pool and the remote worker alike --
    so there is one definition of what an agent's answer is.
    """
    Qk = _agent_generator(k, x, model)
    pik = _solve_component(Qk, int(model.mph[k]), int(model.nlev[k]), model.level[k])
    return Qk, pik


def _compute_equilibrium(x: np.ndarray, model: RCATModelData,
                         exec_backend=None) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """One sweep: every agent's generator and stationary vector at the rates x.

    Agent k is solved in ISOLATION -- its generator reads the rest of the model
    only through the scalar reversed rates x, and it writes only its own slot --
    so this is a fan-out and not a recurrence. That is what lets a thread pool or
    a set of remote workers evaluate the agents and still walk the same iterates
    as the serial loop. Any cross-agent read added here would silently make those
    backends race.
    """
    num_processes = model.num_processes

    if exec_backend is not None:
        return exec_backend.sweep(x, model, _agent_solve, _agent_generator,
                                  _solve_component)

    Q_list = []
    pi_list = []
    for k in range(num_processes):
        Qk, pik = _agent_solve(k, x, model)
        Q_list.append(Qk)
        pi_list.append(pik)

    return pi_list, Q_list


def _exec_backend_for(options, method):
    """Resolve the execution backend named by options.config, or None (serial).

    Kept next to the algorithms rather than in the solver so that a caller
    driving an algorithm class directly gets the same backend selection, and the
    'cluster' + 'inapinf' refusal, as one going through SolverAG.
    """
    if options is None:
        return None
    config = getattr(options, 'config', None)
    if config is None:
        return None
    if not isinstance(config, dict):
        config = {k: getattr(config, k) for k in dir(config)
                  if not k.startswith('_') and not callable(getattr(config, k))}
    from ..exec_backend import create
    return create(config, method)

def _inap_solve(model: RCATModelData, tol: float = 1e-6, max_iter: int = 1000,
                method: str = 'inap', verbose: bool = False,
                exec_backend=None) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray], int]:
    """INAP iterative solver for RCAT model.

    ``exec_backend`` decides who evaluates the agents of each sweep; None is the
    serial loop and is the reference the other backends are asserted against.
    """
    num_actions = model.num_actions
    num_processes = model.num_processes
    N = model.N
    AP = model.AP
    R = model.R

    if num_actions == 0:
        # No actions - solve using local rates only, through the same dispatcher
        # the fixed point uses: this branch carries a whole M/PH/1 on its own,
        # whose marginal spans tens of orders of magnitude over the truncation,
        # and the level recursions are stable there where a null-space solve is
        # not.
        x = np.array([])
        pi = []
        Q = []
        for p in range(num_processes):
            L = R.get((num_actions, p), np.zeros((N[p], N[p])))
            Qp = L - np.diag(L @ np.ones(N[p]))
            Qp = ctmc_makeinfgen(Qp)
            Q.append(Qp)
            pi.append(_solve_component(Qp, int(model.mph[p]), int(model.nlev[p]),
                                       model.level[p]))
        return x, pi, Q, 0

    # see _kb/06-solver-catalog.md (MAM: "AG / RCAT methods") -- non-birth-death
    # processes (catastrophe/batch-removal) take the INAP+ rate-conservation
    # estimator. The test is on the LEVEL distance, not the state distance: with
    # a phase block per level the within-level phase transitions of a PH sit far
    # off the diagonal and are not a departure from birth-death structure.
    not_birth_death = [False] * num_processes
    for k in range(num_processes):
        Lk = R.get((num_actions, k), np.zeros((N[k], N[k])))
        lvl = model.level[k]
        far = np.abs(lvl[:, None] - lvl[None, :]) > 1
        # A PHASE-EXPANDED component takes the same estimator, for the same
        # reason. On a birth-death chain every state-wise reversed rate equals
        # lambda, so their mean is exact; with a phase block per level they do
        # not, the deep truncation levels dominate the unweighted mean, and the
        # mean-of-ratios overestimates the departure rate exactly as it does on
        # a catastrophe (measured on a tandem with Erlang(2) service at Q1: the
        # reversed rate came out 1.27 against the exact 0.5, so flow was not
        # conserved). Rate conservation has no such failure mode.
        not_birth_death[k] = bool(np.any(Lk[far] > 0)) or int(model.mph[k]) > 1

    # Columns of each active matrix that carry any rate. Aa does not depend on
    # x, so this is fixed for the whole fixed point.
    active_cols = []
    for a in range(num_actions):
        k = AP[a, 0]
        Aa = R.get((a, 0), np.zeros((N[k], N[k])))
        active_cols.append(np.flatnonzero(Aa.sum(axis=0) > 0))

    # see _kb/06-solver-catalog.md (MAM: "AG / RCAT methods") -- deterministic
    # start keeps results reproducible across back-ends
    x = np.arange(1, num_actions + 1) / (num_actions + 1)

    # Compute initial equilibrium
    pi, Q = _compute_equilibrium(x, model, exec_backend)

    from ....api.da import da_fpi

    # Reversed-rate fixed point on the isolated-component equilibria, driven
    # by the generic DA driver
    def inap_sweep(pi_prev, itnum):
        nonlocal pi, Q

        # Update each action rate
        for a in range(num_actions):
            k = AP[a, 0]
            Aa = R.get((a, 0), np.zeros((N[k], N[k])))

            if method == 'inapplus' or not_birth_death[k]:
                # INAP+: x(a) = sum_ij Aa(i,j) pi(i), the departure rate of the
                # active component.
                lambda_sum = float(np.sum(pi[k] @ Aa))
                if lambda_sum > 0:
                    x[a] = lambda_sum
            else:
                # INAP: x(a) = mean over the support of the STATE-WISE reversed
                # rate (pi Aa)_j / pi_j, which RCAT requires to be independent
                # of j. On a birth-death component every column of Aa holds one
                # entry, so this is the reference's entrywise mean of
                # Aa(i,j) pi(i) / pi(j) term for term; with a phase block per
                # level a column holds one entry per phase, and only the column
                # form is the reversed rate.
                cols = active_cols[a]
                if cols.size:
                    num = (pi[k] @ Aa)[cols]
                    den = pi[k][cols]
                    ok = (den > 0) & (num > 0)
                    if np.any(ok):
                        ratio = num[ok] / den[ok]
                        ratio = ratio[np.isfinite(ratio)]
                        if ratio.size:
                            x[a] = float(np.mean(ratio))

        # Recompute equilibrium
        pi, Q = _compute_equilibrium(x, model, exec_backend)
        return [p.copy() for p in pi], pi_prev

    def pi_blocknorm(xn, xr):
        e = 0.0
        for k in range(num_processes):
            e = max(e, float(np.sum(np.abs(xn[k] - xr[k]))))
        return e

    _, iteration, _ = da_fpi(inap_sweep, [p.copy() for p in pi], max_iter, tol,
                             norm=pi_blocknorm)

    return x, pi, Q, iteration


def _qbd_scalar_rho(f: float, b: float, g: float) -> float:
    """Sub-unit root rho of the scalar QBD characteristic equation

        b*rho^2 - (f+b+g)*rho + f = 0,

    where f is the up-1 rate, b the down-1 rate and g the extra local outflow
    (catastrophe drain to the empty state). This is the block-size-1 instance
    of Neuts' rate matrix R. For b == 0 it degenerates to rho = f/(f+g).
    """
    if b <= 1e-14:
        if f + g <= 0:
            return np.inf
        return f / (f + g)
    c1 = -(f + b + g)
    disc = c1 * c1 - 4.0 * b * f
    if disc < 0:
        return np.inf
    sq = np.sqrt(disc)
    r1 = (-c1 - sq) / (2.0 * b)
    r2 = (-c1 + sq) / (2.0 * b)
    cands = sorted([r1, r2])
    return cands[0] if cands[0] > 0 else cands[1]


def _qbd_matrix_tail(Qk: np.ndarray, A0: np.ndarray, A1: np.ndarray,
                     A2: np.ndarray, mph: int):
    """Neuts' matrix-geometric solution of one open component with MPH phases
    per level: R from qbd_R_logred, then the boundary equations of levels 0
    and 1,
        pi_0 B00 + pi_1 A2 = 0,   pi_0 A0 + pi_1 (A1 + R A2) = 0,
    normalized by pi_0 e + pi_1 (I - R)^-1 e = 1. Returns None when R has no
    sub-unit spectral radius, i.e. when the isolated component is unstable and
    has no stationary tail to report.

    Logarithmic reduction rather than successive substitutions: this runs once
    per component per fixed-point sweep, and the quadratic convergence is what
    keeps that affordable.
    """
    from ....api.mam.qbd import qbd_R_logred
    try:
        R = qbd_R_logred(A2, A1, A0, iter_max=500)
    except Exception:
        return None
    # The minimal solution of a QBD is NON-NEGATIVE; anything else is the
    # iteration having failed rather than a rate matrix.
    if not np.all(np.isfinite(R)) or np.any(R < -1e-12):
        return None

    B00 = Qk[_blk(0, mph), _blk(0, mph)]
    Sys = np.block([[B00, A0], [A2, A1 + R @ A2]])
    IR = np.eye(mph) - R
    try:
        tail_mass = np.linalg.solve(IR, np.ones(mph))
    except np.linalg.LinAlgError:
        return None
    # STABILITY WITHOUT AN EIGENSOLVER. (I-R)^-1 = I + R + R^2 + ... converges
    # exactly when the spectral radius is below one, and every row of that
    # series is e_i plus non-negative terms, so (I-R)^-1 e >= 1 entrywise. When
    # the isolated component is unstable the series diverges and the inverse
    # picks up negative entries, so this is the spectral condition without an
    # eigensolve (the C++ twin has no LAPACK to call).
    if not np.all(np.isfinite(tail_mass)) or np.any(tail_mass < 1.0 - 1e-9):
        return None
    nrm = np.concatenate([np.ones(mph), tail_mass])
    Sys[:, 0] = nrm
    rhs = np.zeros(2 * mph)
    rhs[0] = 1.0
    try:
        v = np.linalg.solve(Sys.T, rhs)
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(v)):
        return None

    pi0 = v[:mph]
    pi1 = v[mph:]
    busy = np.linalg.solve(IR.T, pi1.T).T        # sum_{n>=1} pi_n
    qlen = float(np.linalg.solve(IR.T, busy.T).T @ np.ones(mph))
    return {'R': R, 'pi0': pi0, 'pi1': pi1, 'busy': busy, 'qlen': qlen}


def _qbd_tail_expand(g: Dict, nlev: int, mph: int) -> np.ndarray:
    """Materialize the matrix-geometric tail over NLEV levels, so the block norm
    of the fixed point and the RCAT residual read one vector shape for every
    component. The metrics use the closed forms in G instead."""
    pi = np.zeros(nlev * mph)
    pi[_blk(0, mph)] = g['pi0']
    v = g['pi1']
    for n in range(1, nlev):
        pi[_blk(n, mph)] = v
        v = v @ g['R']
    return pi


def _compute_equilibrium_qbd(x: np.ndarray, model: RCATModelData,
                             is_open_proc: np.ndarray, exec_backend=None):
    """Solve each isolated component given the current reversed rates x.

    Open components are solved on their infinite state space: with one phase
    per level by the scalar matrix-geometric (QBD / catastrophe) decomposition,
    yielding the exact geometric marginal pi_n = (1-rho) rho^n, and with a phase
    block per level by Neuts' rate matrix R closed by the boundary equations.
    Closed (finite) components fall back to the explicit finite solve. Mirrors
    solver_mam_ag.m (compute_equilibrium_qbd).
    """
    num_processes = model.num_processes
    num_actions = model.num_actions
    N = model.N
    AP = model.AP
    R = model.R

    Q_list = []
    pi_list = []
    rho_proc = np.zeros(num_processes)
    is_geom_proc = np.zeros(num_processes, dtype=bool)
    geom_data = [None] * num_processes

    if exec_backend is not None and hasattr(exec_backend, "sweep_qbd"):
        return exec_backend.sweep_qbd(x, model, is_open_proc, _agent_solve_qbd)

    for k in range(num_processes):
        Qk, pik, rho, isg, gd = _agent_solve_qbd(k, x, model, is_open_proc)
        Q_list.append(Qk)
        pi_list.append(pik)
        rho_proc[k] = rho
        is_geom_proc[k] = isg
        geom_data[k] = gd
    return pi_list, Q_list, rho_proc, is_geom_proc, geom_data



def _agent_solve_qbd(k: int, x: np.ndarray, model: RCATModelData,
                     is_open_proc: np.ndarray):
    """One agent under the matrix-geometric ('inapinf') treatment.

    The 'inapinf' twin of _agent_solve, and the same rule holds: an agent
    reads the rest of the model only through x, so every execution backend
    routes through this one definition.
    """
    num_actions = model.num_actions
    N = model.N
    AP = model.AP
    R = model.R

    Qk_out = None
    pik = None
    rho_out = 0.0
    is_geom = False
    geom = None
    Nk = N[k]
    mph = int(model.mph[k])
    nlev = int(model.nlev[k])

    # Assemble strictly off-diagonal rate matrix for component k
    L = R.get((num_actions, k), np.zeros((Nk, Nk)))
    Off = L.copy()
    np.fill_diagonal(Off, 0.0)
    for c in range(num_actions):
        if AP[c, 1] == k:
            Pb = R.get((c, 1), np.zeros((Nk, Nk)))
            Off = Off + x[c] * Pb
        elif AP[c, 0] == k:
            Aa = R.get((c, 0), np.zeros((Nk, Nk)))
            Off = Off + Aa
    np.fill_diagonal(Off, 0.0)

    Qk = ctmc_makeinfgen(Off - np.diag(Off @ np.ones(Nk)))
    Qk_out = Qk

    solved_geom = False
    if is_open_proc[k] and nlev >= 5:
        if mph == 1:
            s0 = Nk - 2                    # interior state (level s0), 0-indexed
            row = Off[s0, :]
            f = row[s0 + 1]                # up-1 rate (arrival)
            b = row[s0 - 1]                # down-1 rate (service + single removal)
            g0 = row[0]                    # drain to empty state (catastrophe)
            inter_down = float(np.sum(row[1:s0 - 1]))  # batch to interior levels
            if inter_down <= 1e-11 and f > 0:
                rho = _qbd_scalar_rho(f, b, g0)
                if np.isfinite(rho) and 0.0 < rho < 1.0 - 1e-12:
                    rho_out = rho
                    is_geom = True
                    pik = ((1.0 - rho) * rho ** np.arange(Nk))
                    solved_geom = True
        elif _is_block_tridiagonal(Qk, model.level[k]):
            # Read the homogeneous interior blocks one level below the
            # truncation boundary, for the same reason the scalar branch
            # reads the interior row there.
            s0 = nlev - 2
            A0 = Qk[_blk(s0, mph), _blk(s0 + 1, mph)]   # up: arrival
            A1 = Qk[_blk(s0, mph), _blk(s0, mph)]       # local, with diagonal
            A2 = Qk[_blk(s0, mph), _blk(s0 - 1, mph)]   # down: departure
            if np.any(A0.sum(axis=1) > 0):
                g = _qbd_matrix_tail(Qk, A0, A1, A2, mph)
                if g is not None:
                    is_geom = True
                    geom = g
                    pik = (_qbd_tail_expand(g, nlev, mph))
                    solved_geom = True

    if not solved_geom:
        pik = (_solve_component(Qk, mph, nlev, model.level[k]))


    return Qk_out, pik, rho_out, is_geom, geom

def _inap_inf_solve(model: RCATModelData, is_open_proc: np.ndarray,
                    tol: float = 1e-6, max_iter: int = 1000, verbose: bool = False,
                    exec_backend=None):
    """Matrix-geometric INAP fixed point (no truncation) for RCAT product forms.

    Mirrors solver_mam_ag.m (inap_inf): the reversed rate is updated by the
    weighted-mean formula Eq. (4) evaluated in closed form on the geometric
    tail, and the RCAT product-form residual (Remark 2) is returned.

    Reference: A. Marin, S. Rota Bulo, S. Balsamo, "A Numerical Algorithm for
    the Decomposition of Cooperating Structured Markov Processes", MASCOTS 2012.
    """
    num_actions = model.num_actions
    num_processes = model.num_processes
    N = model.N
    AP = model.AP
    R = model.R

    if num_actions == 0:
        x = np.array([])
        pi = []
        Q = []
        for p in range(num_processes):
            L = R.get((num_actions, p), np.zeros((N[p], N[p])))
            Qp = ctmc_makeinfgen(L - np.diag(L @ np.ones(N[p])))
            Q.append(Qp)
            pi.append(_solve_component(Qp, int(model.mph[p]), int(model.nlev[p]),
                                       model.level[p]))
        return (x, pi, Q, 0, np.zeros(num_processes),
                np.zeros(num_processes, dtype=bool), [None] * num_processes, 0.0)

    # Active-transition row sums (rate of the active label out of each state)
    a_row_sum = []
    for a in range(num_actions):
        k = AP[a, 0]
        Aa = R.get((a, 0), np.zeros((N[k], N[k])))
        a_row_sum.append(np.sum(Aa, axis=1))

    # Deterministic initial guess (reproducible across back-ends).
    x = np.arange(1, num_actions + 1) / (num_actions + 1)

    pi, Q, rho_proc, is_geom_proc, geom_data = _compute_equilibrium_qbd(
        x, model, is_open_proc, exec_backend)

    from ....api.da import da_fpi

    # Reversed-rate fixed point (matrix-geometric variant), driven by the
    # generic DA driver
    def inapinf_sweep(pi_prev, itnum):
        nonlocal pi, Q, rho_proc, is_geom_proc, geom_data

        # Reversed-rate update, Eq. (4): x_l = pi^(alpha_l) T^(l) e.
        for a in range(num_actions):
            k = AP[a, 0]
            if is_geom_proc[k]:
                mph = int(model.mph[k])
                if mph == 1:
                    # Geometric tail: the active label fires only in occupied
                    # states, so x = (per-occupied-state rate) * P(occupied).
                    occ = a_row_sum[a][min(1, N[k] - 1)]
                    x[a] = occ * rho_proc[k]
                else:
                    # Matrix-geometric tail: sum_{n>=1} pi_n = pi_1 (I - R)^-1,
                    # and the active label has the same row sums at every busy
                    # level.
                    occ = a_row_sum[a][_blk(1, mph)]
                    x[a] = float(np.dot(geom_data[k]['busy'], occ))
            else:
                x[a] = float(np.dot(pi[k], a_row_sum[a]))

        pi, Q, rho_proc, is_geom_proc, geom_data = _compute_equilibrium_qbd(
            x, model, is_open_proc, exec_backend)
        return [p.copy() for p in pi], pi_prev

    def pi_blocknorm_trunc(xn, xr):
        e = 0.0
        for k in range(num_processes):
            m = min(len(xn[k]), len(xr[k]))
            e = max(e, float(np.sum(np.abs(xn[k][:m] - xr[k][:m]))))
        return e

    _, iteration, _ = da_fpi(inapinf_sweep, [p.copy() for p in pi], max_iter, tol,
                             norm=pi_blocknorm_trunc)

    # RCAT product-form residual (Remark 2): max_l || pi^(alpha_l)(x_l I - T^(l)) ||
    rcat_res = 0.0
    for a in range(num_actions):
        k = AP[a, 0]
        Aa = R.get((a, 0), np.zeros((N[k], N[k])))
        v = pi[k]
        res_vec = x[a] * v - v @ Aa
        rcat_res = max(rcat_res, float(np.linalg.norm(res_vec, 2)))

    return x, pi, Q, iteration, rho_proc, is_geom_proc, geom_data, rcat_res


def _rcat_metrics(sn, x: np.ndarray, pi: List[np.ndarray], Q: List[np.ndarray],
                  model: RCATModelData, rho_proc: Optional[np.ndarray] = None,
                  is_geom_proc: Optional[np.ndarray] = None,
                  geom_data: Optional[List] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Convert RCAT solution to LINE performance metrics.

    When rho_proc/is_geom_proc/geom_data are supplied ('inapinf'), processes
    flagged geometric use the exact closed-form moments of the infinite
    marginal instead of the truncated explicit vector pi.
    """
    M = sn.nstations
    K = sn.nclasses
    process_map = model.process_map
    N = model.N

    rates = sn.rates if hasattr(sn, 'rates') and sn.rates is not None else np.ones((M, K))
    rates = np.asarray(rates, dtype=np.float64)

    njobs = sn.njobs if hasattr(sn, 'njobs') and sn.njobs is not None else np.full(K, np.inf)
    njobs = np.asarray(njobs).flatten()

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))

    # Compute metrics for each (station, class) pair
    for ist in range(M):
        for r in range(K):
            p = int(process_map[ist, r]) - 1
            if p >= 0 and p < len(pi) and len(pi[p]) > 0:
                mph = int(model.mph[p])

                if is_geom_proc is not None and p < len(is_geom_proc) and is_geom_proc[p]:
                    if mph == 1:
                        # Infinite geometric marginal pi_n = (1-rho) rho^n:
                        #   E[N] = rho/(1-rho),  P(N>0) = rho.
                        rho = rho_proc[p]
                        QN[ist, r] = rho / (1.0 - rho)
                        UN[ist, r] = rho
                        mu_ir = rates[ist, r] if ist < rates.shape[0] and r < rates.shape[1] else 0.0
                        if not np.isnan(mu_ir) and mu_ir > 0:
                            TN[ist, r] = mu_ir * rho
                    else:
                        # Matrix-geometric tail pi_(n+1) = pi_n R.
                        g = geom_data[p]
                        QN[ist, r] = g['qlen']
                        UN[ist, r] = float(np.sum(g['busy']))
                        TN[ist, r] = float(np.dot(g['busy'], model.svcdown[p]))
                else:
                    v = np.asarray(pi[p], dtype=float).ravel()

                    # Queue length: E[N] = sum over states of level*pi
                    QN[ist, r] = float(np.dot(model.level[p], v))

                    # Utilization: P(N > 0) = 1 - P(level 0)
                    UN[ist, r] = 1.0 - float(np.sum(v[:mph]))

                    # Throughput: the rate of service completions, i.e. the
                    # phase-dependent departure rate averaged over the marginal.
                    # With one phase this is the mean rate times P(N>0).
                    TN[ist, r] = float(np.dot(model.svcrate[p], v))

    # Response times from Little's law
    for ist in range(M):
        for r in range(K):
            if TN[ist, r] > 0:
                RN[ist, r] = QN[ist, r] / TN[ist, r]

    # System metrics
    CN = np.zeros(K)
    XN = np.zeros(K)

    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        nodetype = list(sn.nodetype) if not isinstance(sn.nodetype, np.ndarray) else sn.nodetype
    else:
        nodetype = [1] * M

    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
        stationToNode = np.asarray(sn.stationToNode).flatten()
    else:
        stationToNode = np.arange(M)

    refstat = sn.refstat if hasattr(sn, 'refstat') and sn.refstat is not None else np.zeros(K)
    refstat = np.asarray(refstat).flatten()

    for r in range(K):
        if r < len(njobs) and njobs[r] >= np.inf:
            # Open class
            for ist in range(M):
                node_idx = int(stationToNode[ist]) if ist < len(stationToNode) else ist
                if node_idx < len(nodetype):
                    ntype = nodetype[node_idx]
                    if hasattr(ntype, 'value'):
                        ntype_val = ntype.value
                    elif hasattr(ntype, 'ID'):
                        ntype_val = ntype.ID
                    else:
                        ntype_val = int(ntype) if not np.isnan(ntype) else 1
                    if ntype_val == 0:  # Source
                        if ist < rates.shape[0] and r < rates.shape[1]:
                            XN[r] = rates[ist, r]
                        break
            CN[r] = np.sum(RN[:, r])
        else:
            # Closed class
            refst = int(refstat[r]) if r < len(refstat) else 0
            if 0 <= refst < M:
                XN[r] = TN[refst, r]
                if XN[r] > 0 and r < len(njobs):
                    CN[r] = njobs[r] / XN[r]

    return QN, UN, RN, TN, CN, XN


class INAPAlgorithm(AGAlgorithm):
    """RCAT iterative fixed-point (INAP) solver.

    Designed for open queueing networks with independent queue states.
    For closed networks, the algorithm may produce approximate results
    since it doesn't model the state space constraint (sum of jobs = N).
    """

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if INAP can solve this network.

        INAP is designed for open networks. For closed networks,
        results may be approximate.
        """
        ok, reason = rcat_supports_processes(sn, 'inap')
        if not ok:
            return False, reason
        # Check if network has closed classes
        njobs = getattr(sn, 'njobs', None)
        if njobs is not None:
            njobs = np.asarray(njobs).flatten()
            has_closed = any(n < np.inf for n in njobs if not np.isnan(n))
            if has_closed:
                return True, "INAP is designed for open networks; closed network results are approximate"
        return True, None

    def solve(self, sn, options=None) -> AGResult:
        """Solve using RCAT INAP method.

        Note: For closed networks, this produces approximate results
        since the RCAT model doesn't enforce the state constraint.
        """
        import warnings
        start_time = time.time()

        # Check for closed classes and warn
        njobs = getattr(sn, 'njobs', None)
        if njobs is not None:
            njobs_arr = np.asarray(njobs).flatten()
            has_closed = any(n < np.inf for n in njobs_arr if not np.isnan(n))
            if has_closed:
                warnings.warn(
                    "INAP is designed for open networks. "
                    "Closed network results are approximate; consider using MVA instead.",
                    UserWarning
                )

        tol = getattr(options, 'tol', 1e-4) if options else 1e-4
        max_iter = getattr(options, 'max_iter', 100) if options else 100
        verbose = getattr(options, 'verbose', False) if options else False

        max_states = 100
        if options and hasattr(options, 'config') and options.config is not None:
            if hasattr(options.config, 'maxStates'):
                max_states = options.config.maxStates
            elif isinstance(options.config, dict) and 'maxStates' in options.config:
                max_states = options.config['maxStates']

        M = sn.nstations
        K = sn.nclasses

        # Build RCAT model
        model = _build_rcat(sn, max_states)

        if model.num_processes == 0:
            # Fall back to simple M/M/c analysis
            return self._solve_simple(sn, options)

        # Solve using INAP
        backend = _exec_backend_for(options, 'inap')
        try:
            x, pi, Q, iterations = _inap_solve(model, tol, max_iter, 'inap', verbose,
                                               backend)
        finally:
            if backend is not None:
                backend.close()

        # Convert to metrics
        QN, UN, RN, TN, CN, XN = _rcat_metrics(sn, x, pi, Q, model)

        runtime = time.time() - start_time

        return AGResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=XN.reshape(1, -1) if XN.ndim == 1 else XN,
            totiter=iterations,
            method="inap",
            runtime=runtime
        )

    def _solve_simple(self, sn, options) -> AGResult:
        """Simple M/M/c fallback when RCAT model can't be built."""
        params = extract_mam_params(sn)
        M = params['nstations']
        K = params['nclasses']
        rates = params['rates']
        nservers = params['nservers']

        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))

        S = 1.0 / np.maximum(rates, 1e-10)

        for m in range(M):
            for k in range(K):
                rho = 0.5  # Default utilization
                UN[m, k] = rho
                RN[m, k] = S[m, k] / (1.0 - rho)
                QN[m, k] = rho * RN[m, k]
                TN[m, k] = rates[m, k] * rho

        return AGResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=TN.copy(),
            totiter=0,
            method="inap",
            runtime=0.0
        )


class INAPPlusAlgorithm(AGAlgorithm):
    """RCAT iterative weighted variant (INAP+)."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """INAP+ shares build_rcat, so it is exponential-only too."""
        return rcat_supports_processes(sn, 'inapplus')

    def solve(self, sn, options=None) -> AGResult:
        """Solve using RCAT INAP+ method."""
        start_time = time.time()

        tol = getattr(options, 'tol', 1e-4) if options else 1e-4
        max_iter = getattr(options, 'max_iter', 100) if options else 100
        verbose = getattr(options, 'verbose', False) if options else False

        max_states = 100
        if options and hasattr(options, 'config') and options.config is not None:
            if hasattr(options.config, 'maxStates'):
                max_states = options.config.maxStates
            elif isinstance(options.config, dict) and 'maxStates' in options.config:
                max_states = options.config['maxStates']

        M = sn.nstations
        K = sn.nclasses

        # Build RCAT model
        model = _build_rcat(sn, max_states)

        if model.num_processes == 0:
            # Fall back to INAP
            algo = INAPAlgorithm()
            return algo.solve(sn, options)

        # Solve using INAP+
        backend = _exec_backend_for(options, 'inapplus')
        try:
            x, pi, Q, iterations = _inap_solve(model, tol, max_iter, 'inapplus',
                                               verbose, backend)
        finally:
            if backend is not None:
                backend.close()

        # Convert to metrics
        QN, UN, RN, TN, CN, XN = _rcat_metrics(sn, x, pi, Q, model)

        runtime = time.time() - start_time

        return AGResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=XN.reshape(1, -1) if XN.ndim == 1 else XN,
            totiter=iterations,
            method="inapplus",
            runtime=runtime
        )


class INAPInfAlgorithm(AGAlgorithm):
    """Matrix-geometric INAP (infinite state space, no truncation).

    Each isolated open component is solved directly on its infinite state
    space by a scalar matrix-geometric (QBD / catastrophe) decomposition,
    giving the exact geometric marginal pi_n = (1-rho) rho^n. Mirrors the
    MATLAB 'inapinf' method (solver_mam_ag.m). See Marin, Rota Bulo, Balsamo,
    MASCOTS 2012.
    """

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """INAPinf shares build_rcat, so it is exponential-only too."""
        return rcat_supports_processes(sn, 'inapinf')

    def solve(self, sn, options=None) -> AGResult:
        start_time = time.time()

        tol = getattr(options, 'tol', 1e-4) if options else 1e-4
        max_iter = getattr(options, 'max_iter', 100) if options else 100
        verbose = getattr(options, 'verbose', False) if options else False

        max_states = 100
        if options and hasattr(options, 'config') and options.config is not None:
            if hasattr(options.config, 'maxStates'):
                max_states = options.config.maxStates
            elif isinstance(options.config, dict) and 'maxStates' in options.config:
                max_states = options.config['maxStates']

        # Build RCAT model
        model = _build_rcat(sn, max_states)

        if model.num_processes == 0:
            # Fall back to INAP
            return INAPAlgorithm().solve(sn, options)

        # Open/closed flag per process (open classes have infinite population).
        njobs = getattr(sn, 'njobs', None)
        njobs = np.full(sn.nclasses, np.inf) if njobs is None else np.asarray(njobs).flatten()
        is_open_proc = np.zeros(model.num_processes, dtype=bool)
        for p in range(1, model.num_processes + 1):
            idx = np.argwhere(model.process_map == p)
            if idx.size > 0:
                r = int(idx[0, 1])
                is_open_proc[p - 1] = r < len(njobs) and np.isinf(njobs[r])

        backend = _exec_backend_for(options, 'inapinf')
        try:
            (x, pi, Q, iterations, rho_proc, is_geom_proc, geom_data,
             rcat_res) = _inap_inf_solve(model, is_open_proc, tol, max_iter, verbose,
                                         backend)
        finally:
            if backend is not None:
                backend.close()

        QN, UN, RN, TN, CN, XN = _rcat_metrics(sn, x, pi, Q, model, rho_proc,
                                               is_geom_proc, geom_data)

        runtime = time.time() - start_time

        return AGResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=XN.reshape(1, -1) if XN.ndim == 1 else XN,
            totiter=iterations,
            method="inapinf",
            runtime=runtime
        )


__all__ = [
    'INAPAlgorithm',
    'INAPPlusAlgorithm',
    'INAPInfAlgorithm',
    'RCATModelData',
]
