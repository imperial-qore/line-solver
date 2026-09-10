"""
Hardware-aware, profiling-calibrated memory guard for the CTMC solver.

This replaces the historical hard-coded state-space size threshold with a
model that (a) probes the amount of memory actually available on the host and
(b) calibrates the per-state cost of a sparse steady-state solve by profiling
the local linear-algebra stack once and caching the result per machine.

The public entry point is ``ctmc_memory_gate`` which, given the worst-case
log state-space size, decides whether solving is safe. It never depends on the
JVM (native-Python constraint); the JAR and MATLAB codebases carry mirrored
implementations that share the same model and cache format.
"""

import json
import math
import os
import platform
import tempfile
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# model constants, kept identical across codebases (see _kb/06-solver-catalog.md, CTMC)
BYTES_PER_NZ = 16.0  # 8 (double value) + 8 (index/pointer amortized), conservative
# Fraction of available memory the solver is allowed to target.
DEFAULT_SAFETY_FRACTION = 0.6
# Grid sizes (n = g*g states) used to fit the fill power law.
CALIB_G1 = 40
CALIB_G2 = 80
# Fallback power-law coefficients if calibration cannot run.
FALLBACK_ALPHA = BYTES_PER_NZ * 8.0
FALLBACK_BETA = 1.3
# Conservative available-memory default (bytes) when the probe fails.
FALLBACK_AVAIL_BYTES = 1.0 * 1024 ** 3
# Dense n-by-n arrays the NATIVE-PYTHON generator builder holds at once, before
# any factorization: the infinitesimal generator, the immediate-transition block
# and at least one event filtration (``handler.py`` allocates one filtration per
# sync action, so 3 is a strict lower bound). MATLAB, the JAR and C++ keep these
# sparse, which the calibrated LU-factor law already covers; this floor exists
# because ``np.zeros((n, n))`` makes the python peak quadratic in the state count
# and the LU law (n^~1.3) cannot see it.
DENSE_BUILDER_ARRAYS = 3.0
BYTES_PER_DENSE_ENTRY = 8.0


def get_available_memory_bytes():
    """Portable available-physical-memory probe.

    Uses psutil when present (works on Windows/macOS/Linux). Falls back to
    POSIX sysconf, then to a conservative constant. Returns bytes, or the
    fallback constant when nothing is available. Never raises.
    """
    try:
        import psutil  # optional dependency; not a JVM dependency
        avail = psutil.virtual_memory().available
        if avail and avail > 0:
            return float(avail)
    except Exception:
        pass
    try:
        # POSIX only; absent on Windows.
        pages = os.sysconf('SC_AVPHYS_PAGES')
        page_size = os.sysconf('SC_PAGE_SIZE')
        if pages > 0 and page_size > 0:
            return float(pages) * float(page_size)
    except Exception:
        pass
    return FALLBACK_AVAIL_BYTES


def _machine_signature():
    return "{}|{}|{}".format(platform.machine(),
                             os.cpu_count() or 1,
                             platform.system())


def _cache_path():
    return os.path.join(tempfile.gettempdir(), 'line_ctmc_calib_python.json')


def _build_lattice_generator(g):
    """Build an (n-1)x(n-1) nonsingular block of a 2D nearest-neighbour CTMC
    generator on a g-by-g lattice (n = g*g). The 5-point stencil mimics the
    quasi-birth-death coupling of multi-station queueing generators, whose
    sparse-LU fill is representative and, crucially, cheap and safe to
    factorize (unlike a random Erdos-Renyi matrix, whose factor is near dense).
    """
    n = g * g
    idx = np.arange(n).reshape(g, g)
    rows = []
    cols = []
    # right/left/up/down neighbours with unit rates
    for shift, axis in ((1, 1), (-1, 1), (1, 0), (-1, 0)):
        src = idx
        dst = np.roll(idx, shift, axis=axis)
        # drop wrap-around edges to keep it a genuine lattice
        if axis == 1:
            valid = np.ones_like(idx, dtype=bool)
            if shift == 1:
                valid[:, -1] = False
            else:
                valid[:, 0] = False
        else:
            valid = np.ones_like(idx, dtype=bool)
            if shift == 1:
                valid[-1, :] = False
            else:
                valid[0, :] = False
        rows.append(src[valid].ravel())
        cols.append(dst[valid].ravel())
    r = np.concatenate(rows)
    c = np.concatenate(cols)
    data = np.ones(r.size)
    Q = sp.csc_matrix((data, (r, c)), shape=(n, n))
    # make an infinitesimal generator: diagonal = -rowsum
    rowsum = np.asarray(Q.sum(axis=1)).ravel()
    Q = Q - sp.diags(rowsum)
    # drop last state to obtain a nonsingular system (as ctmc_solve does)
    A = Q[:n - 1, :n - 1].tocsc()
    return A


def _profile_point(g):
    """Factorize the lattice block of side g; return (n, factor_bytes, secs)."""
    A = _build_lattice_generator(g)
    n = A.shape[0]
    t0 = time.perf_counter()
    lu = spla.splu(A)
    secs = time.perf_counter() - t0
    factnnz = lu.L.nnz + lu.U.nnz
    return n, BYTES_PER_NZ * float(factnnz), secs


def _fit_power_law(x1, y1, x2, y2):
    """Fit y = alpha * x^beta through two positive points."""
    if x1 <= 0 or x2 <= 0 or y1 <= 0 or y2 <= 0 or x1 == x2:
        return None
    beta = math.log(y2 / y1) / math.log(x2 / x1)
    alpha = y1 / (x1 ** beta)
    return alpha, beta


def _load_cache(sig):
    try:
        with open(_cache_path(), 'r') as fh:
            data = json.load(fh)
        if data.get('sig') == sig:
            return data
    except Exception:
        pass
    return None


def _store_cache(data):
    try:
        with open(_cache_path(), 'w') as fh:
            json.dump(data, fh)
    except Exception:
        pass


def get_calibration(verbose=False):
    """Return calibrated {alpha_mem,beta_mem,alpha_t,beta_t}, using a per-machine
    on-disk cache. Profiles the local sparse solver once on a cache miss.
    Falls back to conservative coefficients if profiling fails.
    """
    sig = _machine_signature()
    cached = _load_cache(sig)
    if cached is not None:
        return cached
    try:
        n1, b1, t1 = _profile_point(CALIB_G1)
        n2, b2, t2 = _profile_point(CALIB_G2)
        mem = _fit_power_law(n1, b1, n2, b2)
        tim = _fit_power_law(n1, max(t1, 1e-9), n2, max(t2, 1e-9))
        if mem is None:
            raise ValueError("degenerate memory fit")
        data = {
            'sig': sig,
            'alpha_mem': mem[0], 'beta_mem': mem[1],
            'alpha_t': (tim[0] if tim else 0.0),
            'beta_t': (tim[1] if tim else 1.0),
            'timestamp': time.time(),
        }
        _store_cache(data)
        if verbose:
            print("CTMC calibration: bytes ~ {:.3g}*N^{:.3f}".format(
                mem[0], mem[1]))
        return data
    except Exception as exc:
        if verbose:
            print("CTMC calibration failed ({}); using fallback model".format(exc))
        return {
            'sig': sig,
            'alpha_mem': FALLBACK_ALPHA, 'beta_mem': FALLBACK_BETA,
            'alpha_t': 0.0, 'beta_t': 1.0,
            'timestamp': time.time(),
        }


def predict_log_bytes(log_nstates, calib):
    """Log of the predicted peak solve bytes, as the larger of two footprints:
    the calibrated sparse-LU factor and the dense n-by-n arrays the native
    generator builder materializes before the factorization runs."""
    log_factor = math.log(calib['alpha_mem']) + calib['beta_mem'] * log_nstates
    log_dense = (math.log(DENSE_BUILDER_ARRAYS * BYTES_PER_DENSE_ENTRY)
                 + 2.0 * log_nstates)
    return max(log_factor, log_dense)


def predict_bytes(log_nstates, calib):
    """Predicted peak solve bytes given the natural-log worst-case state count."""
    return math.exp(min(predict_log_bytes(log_nstates, calib), 700.0))


# Largest (m_1..m_K) grid the exact ordered-buffer sum will enumerate. Beyond
# this the estimator falls back to a looser geometric bound rather than spend
# unbounded time inside a memory GATE.
ORDER_GRID_MAX = 1.0e6


def _log_sum_exp(vals):
    finite = [v for v in vals if v > -float('inf')]
    if not finite:
        return -float('inf')
    top = max(finite)
    return top + math.log(sum(math.exp(v - top) for v in finite))


def _log_binom(n, k):
    if k < 0 or k > n:
        return -float('inf')
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _station_terms(caps, cap_total):
    """Allowed occupancy vectors at ONE order-preserving station, with the log of
    the ordering multiplicity of each. `caps[k]` bounds class k; `cap_total`
    bounds the SUM, which is what a finite station capacity actually constrains
    -- it is a buffer of that many slots, not that many PER CLASS. Ignoring it
    charged a capacity-3 station for buffers of length 15 across 5 classes."""
    import itertools
    out = []
    for m in itertools.product(*[range(c + 1) for c in caps]):
        t = sum(m)
        if cap_total is not None and t > cap_total:
            continue
        out.append((m, math.lgamma(t + 1) - sum(math.lgamma(x + 1) for x in m)))
    return out


def _log_ordered_joint(caps_per_station, cap_total_per_station, njobs, m_rem):
    """Log count of (placement, ordering) configurations over ALL
    order-preserving stations at once, with the leftovers spread over `m_rem`
    share stations.

    POPULATION IS CONSERVED. A closed class has N jobs to share out; an OPEN
    class truncated at `cutoff` has at most `cutoff` jobs IN THE NETWORK, which
    is exactly how the plain stars-and-bars term has always treated it. Either
    way a job placed at one station is not available to another. Charging every
    ordered station the full population independently priced a two-station PAS
    model at 1957 * 1957 = 3829849 against a true 5040, and treating an open
    class as per-station independent priced prio_hol_open (3 ordered stations,
    cutoff 1) at 32768 and refused a model that solves in a fraction of a second.
    """
    K = len(njobs)
    state = {tuple(int(n) for n in njobs): 0.0}
    for caps, cap_tot in zip(caps_per_station, cap_total_per_station):
        nxt = {}
        for rem, acc in state.items():
            avail = [min(int(caps[k]), rem[k]) for k in range(K)]
            for m, lmult in _station_terms(avail, cap_tot):
                key = tuple(rem[k] - m[k] for k in range(K))
                v = acc + lmult
                prev = nxt.get(key)
                nxt[key] = v if prev is None else _log_sum_exp([prev, v])
        state = nxt
        if not state:
            return -float('inf')
    vals = []
    for rem, acc in state.items():
        v = acc
        if m_rem >= 1:
            for r in rem:
                v += _log_binom(r + m_rem - 1, r)
        elif any(r > 0 for r in rem):
            continue
        vals.append(v)
    return _log_sum_exp(vals)


def state_space_log_size(sn, options=None):
    """Worst-case log-size of the CTMC state space induced by ``sn``.

    The estimate multiplies four factors, summed in log space: stars-and-bars
    job placements per class (open classes truncated at the cutoff), the
    class-sequence multiplicity of every order-preserving buffer, the
    service-phase multiplicity at each station, and one routing pointer per
    (node, class) doing RROBIN or WRROBIN. Mirrors MATLAB
    ctmc_state_space_logsize.m and JAR MemoryGuard.stateSpaceLogSize.

    Args:
        sn: the network structure, after PH conversion.
        options: solver options carrying the cutoff.

    Returns:
        The natural log of the worst-case number of states.
    """
    from ...sn.network_struct import SchedStrategy, RoutingStrategy

    M = sn.nstations
    K = sn.nclasses
    NK = np.asarray(sn.njobs).flatten() if sn.njobs is not None else np.ones(K)

    cutoff = getattr(options, 'cutoff', None) if options is not None else None
    cutoff_arr = np.asarray(cutoff, dtype=float).flatten() if cutoff is not None else np.array([])
    finite_cutoff = cutoff_arr[np.isfinite(cutoff_arr)] if cutoff_arr.size else cutoff_arr
    if finite_cutoff.size:
        cutoff_scalar = float(np.max(finite_cutoff))
    else:
        # Same default the CTMC analyzer installs for open/mixed models.
        cutoff_scalar = math.ceil(6000 ** (1.0 / (M * K)))

    share_sched = (SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                   SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO,
                   SchedStrategy.GPSPRIO, SchedStrategy.LPS)

    # ORDERED BUFFERS. A station outside the share family keeps the class
    # SEQUENCE of the jobs it holds, so one occupancy vector is as many states as
    # it has sequences. Omitting the factor let gallery_mmap1_multiclass be priced
    # at 726 states against 225840 enumerated -- a 311x UNDER-estimate that neither
    # ctmc_max_states nor AUTO's cap can catch (both sit at 3e6, ABOVE the harmful
    # count), and the solve then attempts one dense (225840,225840) array = 380 GiB.
    #
    # It is computed EXACTLY, not bounded. A geometric bound (sum_t Kb^t) times a
    # placement count over the remaining stations DOUBLE COUNTS -- the sequence
    # term already places the jobs it orders -- and that overlap priced
    # sdroute_twoclasses_closed at 38880 against 1980 true, refusing a working
    # model. The exact joint sum below is
    #     sum over m (jobs at the ordered station, m_k <= n_k) of
    #         multinomial(sum m; m) * prod_k C(n_k - m_k + Mrem - 1, n_k - m_k)
    # which counts placement and ordering together and cannot overlap. Measured:
    # sdroute 9504 (passes, and the residual 4.8x over truth is the phase and
    # round-robin conservatism, a separate open item), mmap1 4.2e6 (refuses).
    #
    # A G-network signal is consumed on arrival and never occupies a buffer slot,
    # so it takes no position in a sequence and keeps its ordinary placement term.
    issignal = getattr(sn, 'issignal', None)
    if issignal is not None and np.size(issignal):
        issignal = np.asarray(issignal).astype(bool).flatten()
    else:
        issignal = np.zeros(K, dtype=bool)
    buffered = [k for k in range(K) if not (k < issignal.size and issignal[k])]

    ord_idx = []
    if len(buffered) > 1:
        for i in range(M):
            sched = sn.sched.get(i) if isinstance(sn.sched, dict) else None
            if sched == SchedStrategy.EXT or sched in share_sched:
                continue
            ord_idx.append(i)
    n_ord = len(ord_idx)

    # A finite station capacity bounds the buffer TOTAL, and a per-class capacity
    # bounds each class; both are read because the ordering term enumerates
    # sequences and is therefore sensitive to them in a way the plain placement
    # count never was.
    cap_v = np.asarray(sn.cap).flatten() if getattr(sn, 'cap', None) is not None else np.array([])
    ccap = np.atleast_2d(sn.classcap) if getattr(sn, 'classcap', None) is not None and np.size(sn.classcap) else None

    def _admitting(k, stations):
        """Stations of `stations` that can hold a job of class k at all.

        A ZERO per-class capacity means the class is disabled there, so it never
        occupies a slot and the placement term must not spread it over that
        station. ld_whittle_bandwidth disables each of its three PS routes for
        the other two classes; counting all M=4 stations priced it at
        C(9,6)^3 = 592704 states, 7852 GB under the quadratic byte model, and the
        gate refused a model whose true space is 7^3 = 343 and solves instantly.
        The ordered branch already reads classcap through caps_per; only the
        placement terms were blind to it.
        """
        n = 0
        for i in stations:
            if ccap is not None and i < ccap.shape[0] and k < ccap.shape[1] and ccap[i, k] == 0:
                continue
            n += 1
        return n

    log_nstates = 0.0
    nk_eff = np.zeros(K)
    for k in range(K):
        nk_eff[k] = cutoff_scalar if not np.isfinite(NK[k]) else float(NK[k])

    if n_ord == 0:
        for k in range(K):
            mk = max(_admitting(k, range(M)), 1)
            log_nstates += (math.lgamma(1 + nk_eff[k] + mk - 1)
                            - math.lgamma(1 + mk - 1) - math.lgamma(1 + nk_eff[k]))
    else:
        rem_idx = [i for i in range(M) if i not in ord_idx]
        for k in range(K):
            if k not in buffered:
                mk = max(_admitting(k, range(M)), 1)
                log_nstates += (math.lgamma(1 + nk_eff[k] + mk - 1)
                                - math.lgamma(1 + mk - 1) - math.lgamma(1 + nk_eff[k]))
        njobs_b = [int(nk_eff[k]) for k in buffered]
        m_rem = len(rem_idx)
        caps_per, captot_per = [], []
        for i in ord_idx:
            per = []
            for bi, k in enumerate(buffered):
                c = njobs_b[bi]
                if ccap is not None and i < ccap.shape[0] and k < ccap.shape[1] and np.isfinite(ccap[i, k]):
                    c = min(c, int(ccap[i, k]))
                per.append(c)
            caps_per.append(per)
            ct = None
            if i < cap_v.size and np.isfinite(cap_v[i]) and cap_v[i] >= 0:
                ct = int(cap_v[i])
            captot_per.append(ct)
        # The DP below carries ONE entry per remaining-population vector and steps
        # it once per ordered station, so the box it walks is prod_k (n_k+1) --
        # the same proxy MATLAB, the JAR and C++ use. Multiplying the per-station
        # box over EVERY station instead is exponential in the station count and
        # contradicts this constant's own name: on mqn_multiserver_fcfs it read
        # 16^5 = 1048576, tripped the 1e6 guard, and sent a DP costing some 1280
        # steps to the geometric fallback, which priced the model at 2.4e13 GB.
        grid = 1.0
        for c in njobs_b:
            grid *= (c + 1)
        if grid <= ORDER_GRID_MAX:
            log_nstates += _log_ordered_joint(caps_per, captot_per, njobs_b, m_rem)
        else:
            total = float(sum(njobs_b))
            log_k = math.log(len(njobs_b))
            log_nstates += n_ord * ((total + 1.0) * log_k - math.log(len(njobs_b) - 1.0)
                                    + math.log1p(-math.exp(-(total + 1.0) * log_k)))
            for k in buffered:
                mk = _admitting(k, rem_idx)
                if mk >= 1:
                    log_nstates += (math.lgamma(1 + nk_eff[k] + mk - 1)
                                    - math.lgamma(1 + mk - 1)
                                    - math.lgamma(1 + nk_eff[k]))

    phasessz = getattr(sn, 'phasessz', None)
    if phasessz is not None and np.size(phasessz):
        phasessz = np.atleast_2d(phasessz)
        nservers = np.asarray(sn.nservers).flatten() if sn.nservers is not None else np.array([])
        for i in range(min(M, phasessz.shape[0])):
            sched = sn.sched.get(i) if isinstance(sn.sched, dict) else None
            for k in range(min(K, phasessz.shape[1])):
                p = float(phasessz[i, k])
                if not np.isfinite(p) or p <= 1:
                    continue
                if sched == SchedStrategy.EXT:
                    m = 1.0
                elif sched in share_sched:
                    m = nk_eff[k]
                else:
                    c = float(nservers[i]) if i < nservers.size else 1.0
                    m = min(nk_eff[k], c)
                if not np.isfinite(m):
                    m = nk_eff[k]
                log_nstates += (math.lgamma(1 + m + p - 1) - math.lgamma(1 + p - 1)
                                - math.lgamma(1 + m))

    routing = getattr(sn, 'routing', None)
    connmatrix = getattr(sn, 'connmatrix', None)
    if routing is not None and np.size(routing) and connmatrix is not None and np.size(connmatrix):
        routing = np.atleast_2d(routing)
        connmatrix = np.atleast_2d(connmatrix)
        for ind in range(min(routing.shape[0], connmatrix.shape[0])):
            nout = int(np.count_nonzero(connmatrix[ind, :]))
            if nout <= 1:
                continue
            nrr = int(np.count_nonzero((routing[ind, :] == RoutingStrategy.RROBIN)
                                       | (routing[ind, :] == RoutingStrategy.WRROBIN)))
            if nrr > 0:
                log_nstates += nrr * math.log(nout)

    return log_nstates


def ctmc_memory_gate(log_nstates, force=False, verbose=False,
                     safety_fraction=DEFAULT_SAFETY_FRACTION):
    """Hardware-aware, calibrated pre-gate for the CTMC solver.

    Args:
        log_nstates: natural log of the worst-case state-space size.
        force: bypass the hard stop (still emits a warning).
        verbose: print the estimate.
        safety_fraction: fraction of available memory allowed as the budget.

    Returns:
        (ok, message). ok is False only when the predicted footprint exceeds
        the budget and force is not set; the caller should then abort.
    """
    avail = get_available_memory_bytes()
    budget = safety_fraction * avail
    calib = get_calibration(verbose=verbose)

    log_pred = predict_log_bytes(log_nstates, calib)
    log_budget = math.log(max(budget, 1.0))
    pred_gb = math.exp(min(log_pred, 700.0)) / 1024 ** 3
    budget_gb = budget / 1024 ** 3

    if log_pred > log_budget:
        msg = ("CTMC predicted peak memory ~{:.2f} GB exceeds the safe budget "
               "~{:.2f} GB ({:.0f}% of {:.2f} GB available). Reduce the state "
               "space (e.g. lower 'cutoff'), use another solver (MVA/NC/FLD), "
               "or set force=true to override.").format(
                   pred_gb, budget_gb, 100 * safety_fraction, avail / 1024 ** 3)
        if not force:
            return False, msg
        if verbose:
            print("Warning (forced): " + msg)
        return True, msg

    if verbose and log_pred > math.log(max(0.5 * budget, 1.0)):
        print("CTMC predicted peak memory ~{:.2f} GB (budget ~{:.2f} GB).".format(
            pred_gb, budget_gb))
    return True, ""
