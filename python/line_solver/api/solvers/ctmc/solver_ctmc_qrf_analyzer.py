"""
Adapter bridging NetworkStruct to QRF library functions within SolverCTMC.

Port of MATLAB solver_ctmc_qrf_analyzer.m.
"""

import time
import numpy as np


def _proc_entry_to_map(entry):
    """Convert a native sn.proc[i] entry to a MAP [D0, D1] pair.

    sn.proc entries are heterogeneous compact PH parameter forms:
      {'rate': r}                 -> exponential
      {'k': k, 'mu': mu}          -> Erlang-k (k phases, rate mu)
      {'probs': p, 'rates': r}    -> hyperexponential
      [alpha, T]                  -> explicit PH (initial row, subgenerator)
    Each is first reduced to a phase-type (alpha, T), then mapped to the
    renewal MAP D0 = T, D1 = t0 @ alpha with t0 = -T @ 1 (the exit-rate
    column re-initialized by alpha). Returns None for empty/disabled service.
    """
    if entry is None:
        return None
    e = entry
    # Unwrap a single-element list wrapper: proc[i] is [<form>].
    if isinstance(e, (list, tuple)) and len(e) == 1 and \
            isinstance(e[0], (list, tuple, dict)):
        e = e[0]
    alpha = None
    T = None
    if isinstance(e, dict):
        if 'rate' in e:
            r = float(e['rate'])
            if r <= 0:
                return None
            alpha = np.array([1.0])
            T = np.array([[-r]])
        elif 'k' in e and 'mu' in e:
            k = int(e['k'])
            mu = float(e['mu'])
            alpha = np.zeros(k); alpha[0] = 1.0
            T = np.zeros((k, k))
            for j in range(k):
                T[j, j] = -mu
                if j + 1 < k:
                    T[j, j + 1] = mu
        elif 'probs' in e and 'rates' in e:
            p = np.asarray(e['probs'], dtype=float).ravel()
            r = np.asarray(e['rates'], dtype=float).ravel()
            alpha = p / p.sum() if p.sum() > 0 else p
            T = np.diag(-r)
        else:
            return None
    elif isinstance(e, (list, tuple)) and len(e) >= 2 and \
            not isinstance(e[0], dict):
        a0 = np.asarray(e[0], dtype=float).ravel()
        M0 = np.atleast_2d(np.asarray(e[1], dtype=float))
        # PH may be stored as [alpha, T] or [D0, D1]; disambiguate by whether
        # the first array is a probability row (nonnegative, sums ~1).
        if a0.size == M0.shape[0] and a0.min() >= -1e-12 and \
                abs(a0.sum() - 1.0) < 1e-6:
            alpha = a0; T = M0
        else:
            # Treat as [D0, D1] already.
            D0 = a0.reshape(M0.shape) if a0.ndim == 1 and \
                a0.size == M0.size else np.atleast_2d(np.asarray(e[0], float))
            D1 = M0
            return [np.atleast_2d(D0), np.atleast_2d(D1)]
    else:
        return None
    t0 = -T @ np.ones(T.shape[0])
    D1 = np.outer(t0, alpha)
    return [np.atleast_2d(T), np.atleast_2d(D1)]


def solver_ctmc_qrf_analyzer(sn, options):
    """Adapter for QRF library functions within SolverCTMC.

    Bridges the LINE sn struct to QRF (Quadratic Reduction Framework) library
    functions for approximating performance metrics of single-class closed
    queueing networks with PH service.

    Args:
        sn: NetworkStruct with network parameters
        options: SolverCTMCOptions with method and config

    Returns:
        QN: Queue lengths [M, K]
        UN: Utilizations [M, K]
        RN: Response times [M, K]
        TN: Throughputs [M, K]
        CN: Cycle times [1, K]
        XN: System throughputs [1, K]
        runtime: Elapsed time in seconds
    """
    t_start = time.time()

    M = sn.nstations
    K = sn.nclasses
    N = int(np.sum(sn.njobs))
    S = np.asarray(sn.nservers).flatten()
    PH = sn.proc

    # QRF only supports single-class closed networks
    if K != 1:
        raise ValueError(
            f"QRF methods only support single-class networks (found {K} classes)."
        )
    njobs = np.asarray(sn.njobs).flatten()
    if np.any(np.isinf(njobs)):
        raise ValueError("QRF methods only support closed networks.")

    # see _kb/06-solver-catalog.md (BA: QRF) -- qrf_noblo_* has no infinite-server notion
    if np.any(np.isinf(S)):
        raise ValueError(
            "QRF methods do not support infinite-server (Delay) stations: the "
            "qrf_noblo_* formulation models every station as a single server. "
            "Use a finite-server model, or SolverMVA/SolverNC for models with "
            "think time.")

    # see _kb/06-solver-catalog.md (BA: QRF "Correction (2026-07-24): native
    # Python HAS ported the no-blocking QRF family") for the fix history
    method_norm = options.method if hasattr(options, 'method') else 'qrf.mmi'
    if method_norm == 'qr':
        method_norm = 'qrf.mmi'

    # Extract MAPs as list of [D0, D1] pairs
    MAPs = []
    K_phases = np.zeros(M, dtype=int)
    for i in range(M):
        mp = _proc_entry_to_map(PH[i])
        if mp is not None:
            D0, D1 = mp[0], mp[1]
            MAPs.append([D0, D1])
            K_phases[i] = D0.shape[0]
        else:
            K_phases[i] = 1
            MAPs.append([np.array([[-1.0]]), np.array([[1.0]])])

    # Build routing matrix (M x M) from sn.rt
    rt_raw = np.asarray(sn.rt)
    rt = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            rt[i, j] = rt_raw[i, j]

    # Extract mu and v arrays from MAPs
    Kmax = int(max(K_phases))
    mu = np.zeros((M, Kmax, Kmax))
    v = np.zeros((M, Kmax, Kmax))
    for i in range(M):
        D0 = MAPs[i][0]
        D1 = MAPs[i][1]
        for h in range(K_phases[i]):
            for k in range(K_phases[i]):
                # see _kb/06-solver-catalog.md (BA: QRF numerical robustness
                # notes) for the mu/v (from-phase,to-phase) indexing convention
                mu[i, h, k] = D1[h, k]
                if h == k:
                    v[i, h, k] = 0.0
                else:
                    v[i, h, k] = D0[h, k]

    method = options.method if hasattr(options, 'method') else 'qrf.mmi'
    config = options.config if hasattr(options, 'config') else {}

    # Dispatch based on method
    if method == 'qrf.mmi':
        from ...mapqn.qrf_noblo_mmi import qrf_noblo_mmi
        MR = 1
        UN_qrf, QN_qrf = qrf_noblo_mmi(M, MR, K_phases, N, mu, v, rt)

    elif method == 'qrf.mem':
        from ...mapqn.qrf_noblo_mem import qrf_noblo_mem
        UN_qrf, QN_qrf = qrf_noblo_mem(MAPs, N, rt)

    elif method == 'qrf.mmi.ld':
        from ...mapqn.qrf_noblo_mmi_ld import qrf_noblo_mmi_ld
        alpha = config.get('qrf_alpha', np.ones((M, N)))
        UN_qrf, QN_qrf = qrf_noblo_mmi_ld(MAPs, N, rt, alpha)

    elif method == 'qrf.mmi.linear':
        from ...mapqn.qrf_noblo_mmi_linear import qrf_noblo_mmi_linear
        alpha = config.get('qrf_alpha', np.ones((M, N)))
        UN_qrf, QN_qrf = qrf_noblo_mmi_linear(MAPs, N, rt, alpha)

    elif method in ('qrf.bas', 'qrf.rsrd'):
        # Delegate to existing LP bounds methods
        if method == 'qrf.bas':
            from ...mapqn.qr_bounds_bas import mapqn_qr_bounds_bas
            from ...mapqn.parameters import QRBoundsBasParameters
            qp = config.get('qrf_params', {})
            # NOTE: the dataclass fields are _M/_N (M/N are read-only
            # properties), so these must be passed positionally-named as _M/_N.
            params = QRBoundsBasParameters(
                _M=M, _N=N, MR=qp.get('MR', 1),
                f=qp.get('f', 1),
                K=K_phases,
                F=np.array(qp.get('F', [N] * M)),
                MM=np.array(qp.get('MM', np.zeros((1, 2)))),
                MM1=np.array(qp.get('MM1', np.zeros((1, M)))),
                ZZ=np.array(qp.get('ZZ', [0])),
                BB=np.array(qp.get('BB', np.zeros((1, M)))),
                mu=[mu[i, :K_phases[i], :K_phases[i]] for i in range(M)],
                v=[v[i, :K_phases[i], :K_phases[i]] for i in range(M)],
                r=rt,
            )
            sol = mapqn_qr_bounds_bas(params, 1, 'max')
            UN_qrf = np.array([sol.get_utilization(i + 1) for i in range(M)])
            QN_qrf = _derive_qn_from_bounds(UN_qrf, M, N, S, MAPs, sn)
        else:
            from ...mapqn.qr_bounds_rsrd import mapqn_qr_bounds_rsrd
            from ...mapqn.parameters import QRBoundsRsrdParameters
            alpha_config = config.get('qrf_alpha', np.ones((M, N)))
            params = QRBoundsRsrdParameters(
                _M=M, _N=N,
                K=K_phases,
                F=np.full(M, N, dtype=int),
                mu=[mu[i, :K_phases[i], :K_phases[i]] for i in range(M)],
                v=[v[i, :K_phases[i], :K_phases[i]] for i in range(M)],
                alpha=alpha_config,
                r=rt,
            )
            sol = mapqn_qr_bounds_rsrd(params, 1, 'max')
            UN_qrf = np.array([sol.get_utilization(i + 1) for i in range(M)])
            QN_qrf = _derive_qn_from_bounds(UN_qrf, M, N, S, MAPs, sn)
    else:
        raise ValueError(f"Unknown QRF method: {method}")

    # Normalize QN to population constraint
    if np.sum(QN_qrf) > 0:
        QN_qrf = QN_qrf / np.sum(QN_qrf) * N

    # Map 1D QRF results to M x K matrices (K=1)
    UN = np.zeros((M, K))
    QN = np.zeros((M, K))
    TN = np.zeros((M, K))
    RN = np.zeros((M, K))
    XN = np.zeros((1, K))
    CN = np.zeros((1, K))

    QN[:, 0] = QN_qrf

    # see _kb/06-solver-catalog.md (BA: QRF numerical robustness notes) for
    # the UN_qrf/throughput derivation via Little's law and visit ratios
    V = _get_visit_ratios(sn)

    refstat = int(np.asarray(sn.refstat).flatten()[0])

    UN_qrf = np.asarray(UN_qrf, dtype=float).flatten()

    stimes = np.array([_map_mean_from_ph(PH[i]) for i in range(M)])
    for i in range(M):
        if not np.isinf(S[i]) and stimes[i] > 0 and V[i, 0] > 0 \
                and UN_qrf[i] > 0:
            XN[0, 0] = UN_qrf[i] * S[i] / (V[i, 0] * stimes[i])
            break
    else:
        # No finite-server station carries load: fall back to the reference
        # station, where QN = X * V * stime holds exactly for a delay.
        stime_ref = stimes[refstat]
        if stime_ref > 0 and V[refstat, 0] > 0:
            XN[0, 0] = QN[refstat, 0] / (V[refstat, 0] * stime_ref)

    # Per-station metrics
    for i in range(M):
        stime = stimes[i]
        if stime > 0:
            TN[i, 0] = XN[0, 0] * V[i, 0]
            if np.isinf(S[i]):
                UN[i, 0] = QN[i, 0]
            else:
                UN[i, 0] = UN_qrf[i]
            if TN[i, 0] > 0:
                RN[i, 0] = QN[i, 0] / TN[i, 0]

    # Cycle time
    if XN[0, 0] > 0:
        CN[0, 0] = N / XN[0, 0]

    # Clean NaN
    for arr in [QN, UN, RN, TN, XN, CN]:
        arr[np.isnan(arr)] = 0.0

    runtime = time.time() - t_start
    return QN, UN, RN, TN, CN, XN, runtime


def _get_visit_ratios(sn):
    """Get visit ratios from sn struct."""
    if hasattr(sn, 'visits') and sn.visits is not None:
        # Sum visit matrices across chains
        V = None
        if isinstance(sn.visits, dict):
            for chain_id, v_mat in sn.visits.items():
                v_arr = np.asarray(v_mat)
                if V is None:
                    V = v_arr.copy()
                else:
                    V += v_arr
        elif isinstance(sn.visits, (list, tuple)):
            for v_mat in sn.visits:
                if v_mat is not None:
                    v_arr = np.asarray(v_mat)
                    if V is None:
                        V = v_arr.copy()
                    else:
                        V += v_arr
        if V is not None:
            return V

    # Fallback: compute from routing matrix
    try:
        from ...sn.transforms import sn_refresh_visits
        sn_refresh_visits(sn)
        return _get_visit_ratios(sn)
    except Exception:
        # Final fallback: ones
        return np.ones((sn.nstations, sn.nclasses))


def _map_mean_from_ph(ph_entry):
    """Compute mean service time from PH representation.

    Returns 0 if the PH entry is None or invalid.
    """
    if ph_entry is None:
        return 0.0
    try:
        mp = _proc_entry_to_map(ph_entry)
        if mp is None:
            return 0.0
        D0 = np.atleast_2d(np.asarray(mp[0], dtype=float))
        D1 = np.atleast_2d(np.asarray(mp[1], dtype=float))

        # map_mean = 1/map_lambda = 1/(pi @ D1 @ e)
        # where pi = stationary distribution of D0+D1
        Q = D0 + D1
        n = Q.shape[0]
        if n == 1:
            rate = D1[0, 0]
            return 1.0 / rate if rate > 0 else 0.0

        # Solve pi*(D0+D1) = 0, pi*e = 1
        A = Q.T.copy()
        A[-1, :] = 1.0
        b = np.zeros(n)
        b[-1] = 1.0
        try:
            pi = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return 0.0

        lam = pi @ np.sum(D1, axis=1)
        return 1.0 / lam if lam > 0 else 0.0
    except Exception:
        return 0.0


def _derive_qn_from_bounds(U_bounds, M, N, S, MAPs, sn):
    """Derive QN from utilization bounds using visit ratios and Little's law."""
    V = _get_visit_ratios(sn)

    # Find XN from first finite-server station with nonzero utilization
    XN_est = 0.0
    for i in range(M):
        if not np.isinf(S[i]) and U_bounds[i] > 0 and V[i, 0] > 0:
            stime = _map_mean_from_ph(sn.proc[i])
            if stime > 0:
                XN_est = U_bounds[i] * S[i] / (V[i, 0] * stime)
                break

    QN_qrf = np.zeros(M)
    for i in range(M):
        stime = _map_mean_from_ph(sn.proc[i])
        if stime > 0:
            TN_i = XN_est * V[i, 0]
            if np.isinf(S[i]):
                QN_qrf[i] = TN_i * stime
            else:
                if U_bounds[i] < 1:
                    QN_qrf[i] = TN_i * stime / (1 - U_bounds[i])
                else:
                    QN_qrf[i] = N

    if np.sum(QN_qrf) > 0:
        QN_qrf = QN_qrf / np.sum(QN_qrf) * N

    return QN_qrf
