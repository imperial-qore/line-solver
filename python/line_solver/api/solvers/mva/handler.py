"""
MVA Solver handler.

Native Python implementation of MVA solver handler that orchestrates
chain aggregation, core MVA computation, and result disaggregation.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    sn_get_demands_chain,
    sn_deaggregate_chain_results,
    sn_has_product_form,
    sn_has_product_form_not_het_fcfs,
    sn_has_open_classes,
)
from ...pfqn.mva import pfqn_mva
from ...pfqn.mvac import pfqn_mvac
from ...pfqn.mvald import pfqn_mvams
from ...pfqn.lcfs import pfqn_lcfsqn_mva


def _sched_matches(sched_value, *strategies) -> bool:
    """Check if a scheduling strategy matches any of the given strategies.

    Handles comparison across different enum class implementations by comparing
    enum names as strings.
    """
    if sched_value is None:
        return False

    # Get the name of the scheduling strategy
    if hasattr(sched_value, 'name'):
        sched_name = sched_value.name
    else:
        sched_name = str(sched_value)

    # Compare against each target strategy
    for strategy in strategies:
        if hasattr(strategy, 'name'):
            target_name = strategy.name
        else:
            target_name = str(strategy)
        if sched_name == target_name:
            return True

    return False


# ---------------------------------------------------------------------------
# Per-method structural gates, shared by SolverMVA.supportsModelMethod and by
# the analyzers below. ONE predicate with TWO callers: a rule kept in two
# places is how the report comes to offer a (solver, method) pair that the
# analyzer then refuses -- or, worse, answers with a table of zeros.
# ---------------------------------------------------------------------------

#: The AMVA algorithms whose recursion is over a CLOSED population vector. Each
#: approximates the arrival-instant queue length E[Q(N-1_r)] from E[Q(N)] and is
#: handed (L, N, Z) alone, with no arrival rate and no rate function, so an open
#: chain gives it nothing to recur on. Canonical spellings only; an 'amva.'
#: prefix is stripped before the list is consulted, exactly as the dispatcher
#: strips it before it selects an algorithm.
MVA_CLOSED_POPULATION_METHODS = (
    'bs', 'aql', 'qsa', 'sqni', 'tay', 'scat', 'lcp', 'chow',
    'pamb', 'pami', 'pamt', 'clust', 'dmlin', 'ab', 'schmidt', 'schmidt-ext',
)

#: Scheduling feature names OUTSIDE the BCMP set {INF, PS, FCFS, SIRO, LCFS-PR}
#: that the base MVA envelope declares. The chain algorithms that walk the
#: stations one by one (the summation method, MVAC, QNA) accept the BCMP set and
#: refuse the rest, so each drops these from its own envelope.
MVA_NON_BCMP_SCHED_FEATURES = (
    'SchedStrategy_HOL', 'SchedStrategy_DPS', 'SchedStrategy_FCFSPRPRIO',
    'SchedStrategy_LCFS', 'SchedStrategy_POLLING', 'SchedStrategy_SJF',
    'SchedStrategy_SRPT', 'SchedStrategy_PSJF', 'SchedStrategy_FB',
    'SchedStrategy_LRPT', 'SchedStrategy_SETF',
    'SchedStrategy_OI', 'SchedStrategy_PAS',
)


def mva_base_method(method) -> str:
    """The method name with any leading 'amva.' alias stripped."""
    m = (method or '').lower()
    return m[5:] if m.startswith('amva.') else m


def mva_is_closed_population_method(method) -> bool:
    """True when `method` names one of the closed-population AMVA algorithms."""
    return mva_base_method(method) in MVA_CLOSED_POPULATION_METHODS


def mva_supports_closed_population(sn, method):
    """(ok, reason) for the closed-population AMVA family.

    Open chains are also expressed in the registry (SolverMVA.getMethodFeatureSet
    drops OpenClass for these methods), which is what keeps them off the report;
    they are repeated here because the analyzer must refuse by name with a
    sentence rather than fall through and answer under a method nobody asked
    for. Strict product form has no registry name at all, so this is its only
    home. Mirrors MATLAB SolverMVA.supportsClosedPopulation.
    """
    if not mva_is_closed_population_method(method):
        return True, ''
    base = mva_base_method(method)
    if sn is None:
        return True, ''
    if sn_has_open_classes(sn):
        return False, (
            "the '%s' method approximates the arrival-instant queue length as a "
            "function of the closed population vector N, so it is defined for "
            "closed models only; use 'default', 'lin', 'qd' or 'qna' for a model "
            "with open classes" % base)
    # ab, schmidt and schmidt-ext ARE the class-dependent FCFS algorithms, so
    # heterogeneous FCFS service means are their subject matter rather than a
    # disqualification.
    check_means = base not in ('ab', 'schmidt', 'schmidt-ext')
    if sn_has_product_form_not_het_fcfs(sn, check_means):
        return True, ''
    return False, (
        "the '%s' method is defined for strict product-form, load-independent "
        "models; use 'default', 'lin' or 'qd' for this model" % base)


def mva_supports_single_class_open(sn, method):
    """(ok, reason) for RQNA and RQT.

    Both decompose an open network into GI/G/1 queues and build one uncertainty
    set per flow out of the first two moments of a SINGLE stream, so a
    multiclass model has no counterpart in their equations. No registry feature
    names a class count, so that half of the rule is structural.

    A FORK-JOIN model is refused too. A Join is a synchronisation node, not a
    queue: it carries no service process, so the index-of-dispersion curve these
    analyzers read off every station does not exist for it, and neither has a
    synchronisation term to put in its place. That half IS nameable, so
    getMethodFeatureSet drops Fork/Join for these two methods as well and this
    is the analyzer's half of it -- without it RQNA dereferenced the absent
    service process and RQT reported an infinite queue length at the Join.

    Mirrors MATLAB SolverMVA.supportsSingleClassOpen.
    """
    base = mva_base_method(method)
    if base not in ('rqna', 'rqt'):
        return True, ''
    if sn is None:
        return True, ''
    label = base.upper()
    nodetype = getattr(sn, 'nodetype', None) or []
    for nt in nodetype:
        if nt in (NodeType.FORK, NodeType.JOIN):
            return False, ("%s decomposes an open network into GI/G/1 queues and has no "
                           "synchronisation term; a Join carries no service process for its "
                           "index of dispersion to be read from. Use SolverMVA's 'default' "
                           "method for a fork-join model." % label)
    if int(getattr(sn, 'nclasses', 1)) != 1:
        return False, ("%s supports single-class open networks only. Use the 'qna' "
                       "method for multiclass models." % label)
    return True, ''


def mva_supports_mvac(sn, method):
    """(ok, reason) for MVAC.

    MVAC (Conway-de Souza e Silva-Lavenberg) is the exact chain recursion over
    single-server fixed-rate (SSFR) queues and infinite-server centres of a
    product-form network. Neither the server count nor product form has a
    registry feature name, so both are structural; the scheduling restriction IS
    nameable and lives in getMethodFeatureSet. Mirrors MATLAB
    SolverMVA.supportsMvac.
    """
    if mva_base_method(method) != 'mvac' or sn is None:
        return True, ''
    if not sn_has_product_form(sn):
        return False, 'MVAC requires a product-form model.'
    sched = getattr(sn, 'sched', None) or {}
    nservers = np.asarray(getattr(sn, 'nservers', []), dtype=float).ravel()
    for ist in range(int(getattr(sn, 'nstations', 0))):
        st = sched.get(ist) if isinstance(sched, dict) else None
        if _sched_matches(st, SchedStrategy.INF, SchedStrategy.EXT):
            continue
        # An infinite count is refused here too, as the reference does: a station
        # scheduled FCFS with infinitely many servers is not an IS centre to MVAC,
        # and infSET below is built from the DISCIPLINE, not from the count.
        if ist < nservers.size and float(nservers[ist]) != 1.0:
            return False, ("MVAC supports single-server (SSFR) queues only; use "
                           "method 'exact' for multiserver stations.")
    return True, ''


def mva_supports_schmidt_ext(njobs, fcfs_rows, method):
    """(ok, reason) for the extended Schmidt method.

    Schmidt's EXTENSION over plain Schmidt is an alpha correction applied at an
    FCFS station, and the correction is computed from the network with ONE
    class-r customer TAGGED, that is at population N - 1_r. A class holding no
    customer has none to tag: the sub-problem is formed at a negative
    population, whose state lattice prod(N+1) collapses to zero and the
    recursion indexes an empty array. Plain `schmidt` forms no such
    sub-problem, which is why the requirement is the -ext arm's alone.

    THE TEST IS STATED AT THE FCFS STATION AND NOT AT A CLASS-DEPENDENT ONE,
    because the four kernels differ on when they form the correction: MATLAB,
    C++ and native python form it only where the station's demands differ by
    class, the JAR forms it at every FCFS station. Stating the union is what
    keeps one rule safe for all four; the case it costs -- an FCFS station whose
    demands are identical across classes, one of them empty -- is one where the
    extension reduces to plain `schmidt`, which stays offered.

    `njobs` and `fcfs_rows` are the population vector and the per-row discipline
    the CALLER'S OWN arm hands the kernel: CHAIN-indexed in MATLAB and the C++
    port, CLASS-indexed in the JAR and native python. That difference belongs to
    those arms and not to this rule, which is the same statement everywhere.
    """
    if mva_base_method(method) != 'schmidt-ext':
        return True, ''
    if njobs is None or fcfs_rows is None or not any(fcfs_rows):
        return True, ''
    N = np.asarray(njobs, dtype=float).ravel()
    for r in range(N.size):
        if np.isfinite(N[r]) and N[r] < 1.0:
            return False, (
                "the 'schmidt-ext' method corrects an FCFS station from the network with "
                "one customer of that class tagged, so it needs every class to hold at "
                "least one customer; class %d holds none. Use 'schmidt' for the "
                "uncorrected recursion." % (r + 1))
    return True, ''


def mva_supports_qna_scheduling(sn):
    """(ok, reason) for QNA's station update.

    The update has an arm for INF, PS and FCFS and none for any other
    discipline, so a SIRO, LCFS, LCFS-PR, HOL or priority station used to leave
    its whole row of Q, U, R and T at zero and the table was returned as a
    solution. The registry expresses this as well (getMethodFeatureSet drops the
    disciplines from QNA's envelope); this is the analyzer's half of it.
    """
    if sn is None:
        return True, ''
    sched = getattr(sn, 'sched', None) or {}
    nodetype = getattr(sn, 'nodetype', None)
    station_to_node = getattr(sn, 'stationToNode', None)
    for ist in range(int(getattr(sn, 'nstations', 0))):
        if nodetype is not None and station_to_node is not None:
            ind = int(station_to_node[ist])
            if 0 <= ind < len(nodetype) and nodetype[ind] == NodeType.JOIN:
                continue  # a Join station carries no service and is skipped
        st = sched.get(ist) if isinstance(sched, dict) else None
        if _sched_matches(st, SchedStrategy.EXT, SchedStrategy.INF,
                          SchedStrategy.PS, SchedStrategy.FCFS):
            continue
        name = getattr(st, 'name', str(st))
        return False, ("QNA decomposes every station as a GI/G/m centre and has "
                       "no arm for %s scheduling. Use the 'default' or 'lin' "
                       "methods." % name)
    return True, ''


@dataclass
class SolverMVAOptions:
    """Options for MVA solver."""
    method: str = 'exact'
    tol: float = 1e-4      # Match MATLAB lineDefaults tol=1e-4
    iter_tol: float = 1e-4 # Match MATLAB lineDefaults iter_tol=1e-4
    iter_max: int = 100    # Match MATLAB lineDefaults iter_max=100
    verbose: bool = False
    # Interlock matrix of Franks (1999), Eq. (4.7), CLASS-indexed: interlock[r,s] is the
    # share of the class-s queue that a class-r arrival must not see, because that work was
    # itself caused by the class-r request. None for every model but the layers of SolverLN.
    interlock: object = None


@dataclass
class SolverMVAReturn:
    """
    Result of MVA solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        lG: Log normalizing constant
        runtime: Runtime in seconds
        method: Method used
        it: Number of iterations
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    lG: float = 0.0
    runtime: float = 0.0
    method: str = "exact"
    it: int = 0


def _solver_mva_lcfsqn(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions],
    lcfs_stat: int,
    lcfspr_stat: int
) -> SolverMVAReturn:
    """
    Specialized MVA solver for LCFS + LCFS-PR 2-station networks.

    Wraps the pfqn_lcfsqn_mva algorithm and maps LINE's data structures
    to/from the algorithm's expected format.

    Args:
        sn: Network structure
        options: Solver options
        lcfs_stat: Index of the LCFS station
        lcfspr_stat: Index of the LCFS-PR station

    Returns:
        SolverMVAReturn with performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    M = sn.nstations
    nclasses = sn.nclasses
    njobs = sn.njobs if sn.njobs is not None else np.zeros(nclasses)

    # Extract service times for each class at each station
    # alpha(r) = mean service time at LCFS station for class r
    # beta(r) = mean service time at LCFS-PR station for class r
    alpha = np.zeros(nclasses)
    beta = np.zeros(nclasses)

    rates = sn.rates
    for r in range(nclasses):
        if njobs[r] > 0:
            mu_lcfs = rates[lcfs_stat, r]
            mu_lcfspr = rates[lcfspr_stat, r]

            if mu_lcfs <= 0 or not np.isfinite(mu_lcfs):
                raise RuntimeError(f"Invalid service rate at LCFS station for class {r+1}.")
            if mu_lcfspr <= 0 or not np.isfinite(mu_lcfspr):
                raise RuntimeError(f"Invalid service rate at LCFS-PR station for class {r+1}.")

            alpha[r] = 1.0 / mu_lcfs
            beta[r] = 1.0 / mu_lcfspr

    # Get population vector
    N = njobs.copy()

    # Call the LCFS MVA algorithm
    result = pfqn_lcfsqn_mva(alpha, beta, N)

    # Map results back to LINE format
    Q = np.zeros((M, nclasses))
    U = np.zeros((M, nclasses))
    T = np.zeros((M, nclasses))
    R = np.zeros((M, nclasses))
    X = np.zeros((1, nclasses))
    C = np.zeros((1, nclasses))

    # Map queue lengths
    Q[lcfs_stat, :] = result.Q[0, :]
    Q[lcfspr_stat, :] = result.Q[1, :]

    # Map utilizations
    U[lcfs_stat, :] = result.U[0, :]
    U[lcfspr_stat, :] = result.U[1, :]

    # Throughput is the same at all stations in a closed network
    for r in range(nclasses):
        if njobs[r] > 0:
            X[0, r] = result.T[0, r]
            T[lcfs_stat, r] = result.T[0, r]
            T[lcfspr_stat, r] = result.T[0, r]

    # Compute response times: R = Q / T (Little's Law)
    for k in [lcfs_stat, lcfspr_stat]:
        for r in range(nclasses):
            if T[k, r] > 0:
                R[k, r] = Q[k, r] / T[k, r]

    # Compute cycle times: C = sum of response times at all stations
    for r in range(nclasses):
        if njobs[r] > 0:
            C[0, r] = R[lcfs_stat, r] + R[lcfspr_stat, r]

    return SolverMVAReturn(
        Q=Q,
        U=U,
        R=R,
        T=T,
        C=C,
        X=X,
        lG=np.nan,  # Log of normalizing constant not computed by this method
        runtime=time.time() - start_time,
        method=options.method,
        it=0
    )


def solver_mva(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    MVA solver handler.

    Performs Mean Value Analysis by:
    1. Aggregating class-level parameters into chains
    2. Separating delay stations from queueing stations
    3. Computing chain-level performance using pfqn_mvams
    4. Disaggregating results back to class level

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMVAReturn with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations (non-product-form,
                     LCFS without LCFS-PR, etc.)
    """
    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    M = sn.nstations
    K = sn.nclasses

    # Check for LCFS scheduling
    sched_dict = sn.sched if sn.sched else {}
    lcfs_stats = []
    lcfspr_stats = []
    for i in range(M):
        station_sched = sched_dict.get(i)
        if _sched_matches(station_sched, SchedStrategy.LCFS):
            lcfs_stats.append(i)
        elif _sched_matches(station_sched, SchedStrategy.LCFSPR):
            lcfspr_stats.append(i)

    # Handle LCFS + LCFS-PR 2-station network
    if lcfs_stats and lcfspr_stats:
        if len(lcfs_stats) != 1 or len(lcfspr_stats) != 1:
            raise RuntimeError("LCFS MVA requires exactly one LCFS and one LCFS-PR station.")

        # Check for closed network
        Nchain = sn.njobs if sn.njobs is not None else np.zeros(K)
        if np.any(np.isinf(Nchain)):
            raise RuntimeError("LCFS MVA requires a closed queueing network.")

        # Check for self-loops
        rt = sn.rt
        nclasses = sn.nclasses
        for ist in lcfs_stats + lcfspr_stats:
            for r in range(nclasses):
                if rt is not None and rt[(ist) * nclasses + r, (ist) * nclasses + r] > 0:
                    raise RuntimeError("LCFS MVA does not support self-loops at stations.")

        # Call specialized LCFS MVA solver
        return _solver_mva_lcfsqn(sn, options, lcfs_stats[0], lcfspr_stats[0])
    elif lcfs_stats:
        raise RuntimeError("LCFS scheduling requires a paired LCFS-PR station.")

    # For non-LCFS models, check product-form requirement. METHOD 'mva' IS THE
    # DELIBERATE APPROXIMATION: the dispatch warns that the exact recursion is
    # being run outside its hypotheses and promises an answer, so raising here
    # would contradict its own message. Only an implicit or 'exact' request is
    # refused.
    if not sn_has_product_form(sn) and getattr(options, 'method', None) != 'mva':
        raise RuntimeError(
            "Unsupported exact MVA analysis, the model does not have a product form"
        )

    # Separate infinite server (delay) stations from queueing stations
    infSET = []  # Delay stations
    qSET = []    # Queueing stations

    for i in range(M):
        station_sched = sched_dict.get(i)

        if _sched_matches(station_sched, SchedStrategy.EXT):
            continue
        elif _sched_matches(station_sched, SchedStrategy.INF):
            infSET.append(i)
        elif _sched_matches(station_sched, SchedStrategy.PS, SchedStrategy.LCFSPR,
                           SchedStrategy.FCFS, SchedStrategy.SIRO):
            qSET.append(i)
        else:
            sched_name = str(station_sched) if station_sched else "unknown"
            raise RuntimeError(f"Unsupported exact MVA analysis for {sched_name} scheduling")

    # Get chain-level demands
    chain_result = sn_get_demands_chain(sn)
    Lchain = chain_result.Lchain
    STchain = chain_result.STchain
    Vchain = chain_result.Vchain
    alpha = chain_result.alpha
    Nchain = chain_result.Nchain.flatten()
    refstatchain = chain_result.refstatchain

    nservers = sn.nservers
    if nservers is None:
        nservers = np.ones(M)
    else:
        nservers = nservers.flatten()

    C_chains = sn.nchains

    # Initialize chain-level result matrices
    Uchain = np.zeros((M, C_chains))
    Tchain = np.zeros((M, C_chains))
    Wchain = np.zeros((M, C_chains))
    Qchain = np.zeros((M, C_chains))
    Xchain = np.zeros((1, C_chains))
    C_result = np.zeros((1, C_chains))
    lG = 0.0

    # Handle open classes
    ocl = [i for i in range(len(Nchain)) if np.isinf(Nchain[i])]
    lambda_chain = np.zeros((1, C_chains))

    for r in ocl:
        ref_stat = int(refstatchain[r, 0]) if refstatchain.ndim > 1 else int(refstatchain[r])
        if STchain[ref_stat, r] > 0:
            lambda_chain[0, r] = 1.0 / STchain[ref_stat, r]
        Qchain[ref_stat, r] = np.inf

    # Classes with non-zero population
    rset = [i for i in range(C_chains) if Nchain[i] != 0.0]

    # Build Lp and Zp for queueing and delay stations
    Lp = np.zeros((len(qSET), C_chains))
    for idx, q_idx in enumerate(qSET):
        for j in range(C_chains):
            Lp[idx, j] = STchain[q_idx, j] * Vchain[q_idx, j]

    Zp = np.zeros((len(infSET), C_chains))
    for idx, inf_idx in enumerate(infSET):
        for j in range(C_chains):
            Zp[idx, j] = STchain[inf_idx, j] * Vchain[inf_idx, j]

    # Get nservers for queueing stations
    nserversp = np.array([nservers[q] for q in qSET]).reshape(-1, 1) if qSET else np.ones((1, 1))

    # Compute Z (think times) as sum of delay demands
    Z_total = Zp.sum(axis=0) if Zp.size > 0 else np.zeros(C_chains)

    # Call MVA core algorithm
    _use_mvac = getattr(options, 'method', 'exact') == 'mvac'
    if _use_mvac:
        # MVAC (Conway-de Souza e Silva-Lavenberg 1989): closed product-form
        # networks of single-server fixed-rate queues plus IS centers only.
        # One predicate for the gate and the run: SolverMVA.supportsModelMethod
        # asks mva_supports_mvac before the report offers 'mvac', so a listed
        # row is a row that runs and the refusal reads the same either way.
        if any(np.isinf(Nchain)):
            raise RuntimeError("MVAC supports closed models only; use method 'exact' for open/mixed networks.")
        _mvac_ok, _mvac_reason = mva_supports_mvac(sn, 'mvac')
        if not _mvac_ok:
            raise RuntimeError(_mvac_reason)
    if len(qSET) > 0:
        if _use_mvac:
            Xchain_out, Qpf, _Upf_mvac, _Cpf_mvac = pfqn_mvac(Lp, Nchain, Z_total)
            lG = np.nan
        else:
            # Interlocked flow (Franks 1999, Eq. 4.7): a request cannot queue behind work
            # that its own submission caused, so the arrival-instant queue drops the
            # interlocked share of the other chains. SolverLN supplies the matrix.
            from ...sn import sn_interlock_chain
            _IL = sn_interlock_chain(sn, getattr(options, 'interlock', None))
            # the interlocked recursion is a separate entry point: pfqn_mvams and the
            # pfqn_mva family it dispatches to carry the standard arrival theorem only
            if _IL is None or np.size(_IL) == 0:
                Xchain_out, Qpf, Upf, Cpf, lG = pfqn_mvams(
                    lambda_chain.flatten(), Lp, Nchain, Z_total,
                    mi=np.ones(len(qSET)),
                    S=nserversp.flatten().astype(int)
                )
            else:
                from ...pfqn import pfqn_mvams_ilock
                Xchain_out, Qpf, Upf, Cpf, lG = pfqn_mvams_ilock(
                    lambda_chain.flatten(), Lp, Nchain, Z_total,
                    mi=np.ones(len(qSET)),
                    S=nserversp.flatten().astype(int),
                    IL=_IL
                )

        # Map results back to full station indices
        for idx, q_idx in enumerate(qSET):
            for j in range(C_chains):
                Qchain[q_idx, j] = Qpf[idx, j] if Qpf.ndim > 1 else Qpf[idx]

        Xchain = Xchain_out.reshape(1, -1) if Xchain_out.ndim == 1 else Xchain_out
    else:
        # No queueing stations - compute throughput from delays only
        for r in range(C_chains):
            if np.isfinite(Nchain[r]) and Z_total[r] > 0:
                Xchain[0, r] = Nchain[r] / Z_total[r]

    # Compute queue lengths at delay stations
    for idx, inf_idx in enumerate(infSET):
        for j in range(C_chains):
            Qchain[inf_idx, j] = Xchain[0, j] * STchain[inf_idx, j] * Vchain[inf_idx, j]

    # Compute waiting times
    ccl = [i for i in range(len(Nchain)) if np.isfinite(Nchain[i])]
    for r in rset:
        for k in infSET:
            Wchain[k, r] = STchain[k, r]
        for k in qSET:
            if np.isinf(nservers[k]):
                Wchain[k, r] = STchain[k, r]
            else:
                if Vchain[k, r] == 0.0 or Xchain[0, r] == 0.0:
                    Wchain[k, r] = 0.0
                else:
                    Wchain[k, r] = Qchain[k, r] / (Xchain[0, r] * Vchain[k, r])

    # Compute response times and throughputs
    for r in rset:
        W_sum = np.sum(Vchain[:, r] * Wchain[:, r])
        if W_sum == 0.0:
            Xchain[0, r] = 0.0
        else:
            if np.isinf(Nchain[r]):
                C_result[0, r] = W_sum
            elif Nchain[r] == 0.0:
                Xchain[0, r] = 0.0
                C_result[0, r] = 0.0
            else:
                C_result[0, r] = W_sum
                Xchain[0, r] = Nchain[r] / C_result[0, r]

        for k in range(M):
            Qchain[k, r] = Xchain[0, r] * Vchain[k, r] * Wchain[k, r]
            Tchain[k, r] = Xchain[0, r] * Vchain[k, r]

    # Compute utilizations
    for k in range(M):
        for r in rset:
            if np.isinf(nservers[k]):
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[0, r]
            else:
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[0, r] / nservers[k]

    # Utilization capping for FCFS/PS stations
    for k in range(M):
        station_sched = sched_dict.get(k)
        for r in range(C_chains):
            if Vchain[k, r] * STchain[k, r] > options.tol:
                if _sched_matches(station_sched, SchedStrategy.FCFS, SchedStrategy.PS):
                    Urow_sum = np.sum(Uchain[k, :])
                    if Urow_sum > 1 + options.tol:
                        denom = np.sum(Vchain[k, :] * STchain[k, :] * Xchain[0, :])
                        if denom > 0:
                            Uchain[k, r] = (min(1.0, Urow_sum) * Vchain[k, r] *
                                           STchain[k, r] * Xchain[0, r] / denom)

    # Compute Rchain from Qchain and Tchain
    Rchain = np.zeros((M, C_chains))
    for i in range(M):
        for j in range(C_chains):
            if Tchain[i, j] > 0:
                Rchain[i, j] = Qchain[i, j] / Tchain[i, j]

    # Clean up non-finite values
    Xchain = np.where(~np.isfinite(Xchain), 0.0, Xchain)
    Uchain = np.where(~np.isfinite(Uchain), 0.0, Uchain)
    Qchain = np.where(~np.isfinite(Qchain), 0.0, Qchain)
    Rchain = np.where(~np.isfinite(Rchain), 0.0, Rchain)

    # Zero out results for zero population chains
    Nzero = [i for i in range(C_chains) if Nchain[i] == 0.0]
    for j in Nzero:
        Xchain[0, j] = 0.0
        Uchain[:, j] = 0.0
        Qchain[:, j] = 0.0
        Rchain[:, j] = 0.0
        Tchain[:, j] = 0.0
        Wchain[:, j] = 0.0

    # Disaggregate chain results to class level
    deagg_result = sn_deaggregate_chain_results(
        sn, Lchain, None, STchain, Vchain, alpha,
        None, None, Rchain, Tchain, None, Xchain
    )

    result = SolverMVAReturn(
        Q=deagg_result.Q,
        U=deagg_result.U,
        R=deagg_result.R,
        T=deagg_result.T,
        C=deagg_result.C,
        X=deagg_result.X,
        lG=lG,
        runtime=time.time() - start_time,
        method=options.method,
        it=0
    )

    return result


__all__ = [
    'solver_mva',
    'SolverMVAReturn',
    'SolverMVAOptions',
]
