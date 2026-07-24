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


def predict_bytes(log_nstates, calib):
    """Predicted peak solve bytes given the natural-log worst-case state count."""
    return math.exp(math.log(calib['alpha_mem']) + calib['beta_mem'] * log_nstates)


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

    log_pred = math.log(calib['alpha_mem']) + calib['beta_mem'] * log_nstates
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
