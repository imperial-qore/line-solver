"""
lang="java" dispatch for native Python solvers.

This module lets a native solver delegate a single solve to the canonical
``jline.jar`` (the same JAR used by MATLAB and the Python wrapper) instead of
running the pure-Python numerics. It mirrors the MATLAB mechanism in
``@NetworkSolver/getAvgTable.m`` where ``options.lang='java'`` hands the solve
to JLINE.

Transport is subprocess + JSON, never JPype/in-process JVM:

    native model --(linemodel_io.save_model)--> model.json
                 --(java -jar jline.jar -i json -s <solver> -a avg -o json)-->
    results JSON --(parsed here)--> pandas AvgTable

CRITICAL (JVM-free guarantee): nothing in this module runs, and this module is
not even imported, unless the user explicitly sets ``lang='java'``. The Java
binary and ``jline.jar`` are resolved lazily at call time, never at import or
package install. A user who never sets ``lang='java'`` needs no JVM.
"""

import json
import math as _math
import os
import shutil
import subprocess
import tempfile

import numpy as np


# Map a native solver's getName() to the CLI solver token accepted by
# jline.jar's LineCLI (-s option). Tokens must match validateSolver().
_SOLVER_TOKENS = {
    'MVA': 'mva',
    'NC': 'nc',
    'CTMC': 'ctmc',
    'MAM': 'mam',
    'FLUID': 'fluid',
    'FLD': 'fluid',
    'SSA': 'ssa',
    'JMT': 'jmt',
    'LDES': 'ldes',
    'AUTO': 'auto',
    'LN': 'ln',        # refined by _ln_token: 'ln' is NC layers, 'ln.mva' is MVA layers
    'LQNS': 'lqns',
}

# Native SolverOptions defaults (see python/line_solver/solvers.py SolverOptions).
# These are the generic, solver-type-agnostic defaults; the JAR applies its own
# per-solver-type defaults (e.g. CTMC iter_max=100, LQN iter_max=200). We forward
# a numeric option only when it differs from the generic native default, so an
# untouched option lets the JAR use its per-solver default (matching the wrapper),
# while a user override is honored.
_NATIVE_DEFAULT_TOL = 1e-4
_NATIVE_DEFAULT_ITER_TOL = 1e-4
_NATIVE_DEFAULT_ITER_MAX = 1000

def find_java():
    """
    Locate a Java launcher lazily. Resolution order:
      1. $LINE_JAVA (explicit override)
      2. $JAVA_HOME/bin/java
      3. `java` on PATH

    Raises RuntimeError with actionable guidance if none is found.
    """
    cand = os.environ.get('LINE_JAVA')
    if cand:
        # Explicit override is authoritative: fail loudly rather than silently
        # falling back, so a misconfiguration is not masked.
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
        raise RuntimeError("$LINE_JAVA is set to '%s' but it is not an executable file." % cand)

    java_home = os.environ.get('JAVA_HOME')
    if java_home:
        cand = os.path.join(java_home, 'bin', 'java')
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand

    cand = shutil.which('java')
    if cand:
        return cand

    raise RuntimeError(
        "lang='java' requires a Java runtime, but no 'java' launcher was found. "
        "Install a JRE (>= 8) and put it on PATH, or set $LINE_JAVA / $JAVA_HOME. "
        "Native solving needs no JVM: leave lang at its default ('python')."
    )


def find_jar():
    """
    Locate jline.jar lazily. Resolution order:
      1. $LINE_JLINE_JAR (explicit override)
      2. a jar bundled inside the installed package (line_solver/**/jline.jar)
      3. the source-tree copy at <repo>/common/jline.jar, found by walking up
      4. Auto-download from SourceForge (https://line-solver.sourceforge.net/latest/jline.jar)

    Raises RuntimeError with actionable guidance if none is found and download fails.
    """
    cand = os.environ.get('LINE_JLINE_JAR')
    if cand:
        # Explicit override is authoritative: fail loudly rather than silently
        # falling back, so a misconfiguration is not masked.
        if os.path.isfile(cand):
            return cand
        raise RuntimeError("$LINE_JLINE_JAR is set to '%s' but no such file exists." % cand)

    here = os.path.dirname(os.path.abspath(__file__))

    # Bundled inside the package (e.g. installed via a [java] extra).
    pkg_root = os.path.dirname(here)  # .../line_solver
    for sub in ('jline.jar', os.path.join('jar', 'jline.jar'),
                os.path.join('common', 'jline.jar')):
        cand = os.path.join(pkg_root, sub)
        if os.path.isfile(cand):
            return cand

    # Source checkout: walk up looking for common/jline.jar.
    d = here
    for _ in range(8):
        cand = os.path.join(d, 'common', 'jline.jar')
        if os.path.isfile(cand):
            return cand
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent

    # Try to auto-download from SourceForge
    try:
        return _download_jline_jar()
    except Exception as e:
        raise RuntimeError(
            "lang='java' requires jline.jar, but it was not found and download failed.\n"
            f"Error: {e}\n"
            "Options:\n"
            "  1. Set $LINE_JLINE_JAR to the path of jline.jar\n"
            "  2. Build it: cd jar && mvn clean package -P b\n"
            "  3. Download manually: https://line-solver.sourceforge.net/latest/jline.jar\n"
            "  4. Use native solving (lang='python') which needs no JAR"
        )


def _download_jline_jar() -> str:
    """Download versioned jline.jar from SourceForge and rename to jline.jar.

    Works on all OS (Windows, macOS, Linux) using platform-agnostic operations.
    """
    import urllib.request
    import shutil

    # Determine target location (prefer common/ in repo, else package dir)
    here = os.path.dirname(os.path.abspath(__file__))

    # Try to find/create common/ directory
    d = here
    common_dir = None
    for _ in range(8):
        cand = os.path.join(d, 'common')
        if os.path.isdir(cand) or (not os.path.exists(cand)):
            common_dir = cand
            break
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent

    if not common_dir:
        # Fall back to package directory
        pkg_root = os.path.dirname(here)
        common_dir = os.path.join(pkg_root, 'common')

    os.makedirs(common_dir, exist_ok=True)

    # Get version from GlobalConstants (with fallback to generic jline.jar)
    version = None
    try:
        from ...constants import GlobalConstants
        version = GlobalConstants.Version
    except Exception:
        pass

    jline_path = os.path.join(common_dir, 'jline.jar')

    if version:
        versioned_jar_filename = f'jline-{version}.jar'
        versioned_jar_path = os.path.join(common_dir, versioned_jar_filename)

        # If versioned jar already exists, just copy to jline.jar
        if os.path.isfile(versioned_jar_path):
            if not os.path.isfile(jline_path):
                shutil.copy2(versioned_jar_path, jline_path)
            return jline_path

        # Download versioned jar from SourceForge
        jline_url = f'https://line-solver.sourceforge.net/latest/{versioned_jar_filename}'
        print(f"Downloading {versioned_jar_filename} from SourceForge...")

        try:
            urllib.request.urlretrieve(jline_url, versioned_jar_path)
            if not os.path.isfile(versioned_jar_path):
                raise RuntimeError(f"Download completed but file not found: {versioned_jar_path}")

            file_size = os.path.getsize(versioned_jar_path) / (1024*1024)
            print(f"Downloaded {versioned_jar_filename} ({file_size:.1f} MB)")

            # Copy to jline.jar (works on all OS)
            shutil.copy2(versioned_jar_path, jline_path)
            print(f"Copied: {versioned_jar_filename} -> jline.jar")

            return jline_path
        except Exception as e:
            # Clean up partial files if download failed
            for path in (versioned_jar_path, jline_path):
                if os.path.exists(path):
                    try:
                        os.remove(path)
                    except Exception:
                        pass
            raise RuntimeError(f"Failed to download {versioned_jar_filename}: {e}")
    else:
        # Fallback: download generic jline.jar if version detection fails
        jline_url = 'https://line-solver.sourceforge.net/latest/jline.jar'
        print(f"Downloading jline.jar from SourceForge (version detection unavailable)...")

        try:
            urllib.request.urlretrieve(jline_url, jline_path)
            if not os.path.isfile(jline_path):
                raise RuntimeError(f"Download completed but file not found: {jline_path}")

            file_size = os.path.getsize(jline_path) / (1024*1024)
            print(f"Downloaded jline.jar ({file_size:.1f} MB)")
            return jline_path
        except Exception as e:
            if os.path.exists(jline_path):
                try:
                    os.remove(jline_path)
                except Exception:
                    pass
            raise RuntimeError(f"Failed to download jline.jar: {e}")


def _get_model(solver):
    """Resolve the native Network model from a solver, tolerating the two
    attribute conventions in the codebase (.model and .network)."""
    for attr in ('model', 'network'):
        m = getattr(solver, attr, None)
        if m is not None:
            return m
    raise RuntimeError(
        "lang='java' could not locate the model on solver '%s' "
        "(neither .model nor .network)." % type(solver).__name__
    )


def _scalar_cutoff(cutoff):
    """
    Return a finite scalar cutoff to forward to the JAR CLI, or None.

    The native CTMC/SSA ``options.cutoff`` is usually a scalar (int/float) but
    may be a per-(station,class) array. Only a finite uniform scalar can be
    forwarded via the CLI's ``--cutoff`` flag; anything else (None, infinite,
    NaN, or a non-uniform array) is left to the JAR's own default so behaviour
    is never silently misrepresented.
    """
    if cutoff is None:
        return None
    try:
        import numpy as _np
        arr = _np.asarray(cutoff, dtype=float)
    except Exception:
        return None
    if arr.size == 0:
        return None
    vals = arr.reshape(-1)
    first = vals[0]
    if not _math.isfinite(first):
        return None
    for v in vals:
        if not _math.isfinite(v) or v != first:
            return None
    return float(first)


def _forward_if_overridden(value, native_default, rtol=1e-12):
    """
    Return ``float(value)`` when it is a finite scalar that differs from the
    generic native default, else None.

    A ``None`` (unset) or a value equal to the native default is left for the
    JAR to fill from its own per-solver-type default, so an untouched option
    does not clobber the JAR's tuned default (matching wrapper behaviour); only
    a genuine user override is forwarded.
    """
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if not _math.isfinite(v):
        return None
    if abs(v - float(native_default)) <= abs(float(native_default)) * rtol:
        return None
    return v


def _ln_token(solver):
    """
    Resolve the CLI token for a layered solve from the layer solver actually in
    use. LineCLI's 'ln' token builds NC layer solvers and 'ln.mva' builds MVA
    ones, whereas the native SolverLN defaults to MVA and only uses NC when the
    caller passes a factory that produces it. Sending 'ln' unconditionally
    silently substituted the layer solver, which changes results: an RR cache
    layer, for instance, is exact under MVA's 'fpi' but approximate under NC's
    'spm'.
    """
    layer_name = None
    layers = getattr(solver, 'solvers', None)
    if layers:
        probe = layers[0]
        if hasattr(probe, 'getName'):
            try:
                layer_name = probe.getName()
            except Exception:
                layer_name = None
        if layer_name is None:
            layer_name = type(probe).__name__.replace('Solver', '')
    if layer_name is None:
        # No layers built yet: fall back to the native default layer solver.
        layer_name = 'MVA'
    return 'ln' if str(layer_name).upper() == 'NC' else 'ln.mva'


def _solver_token(solver):
    """Resolve the CLI solver token from a native solver instance."""
    name = None
    if hasattr(solver, 'getName'):
        try:
            name = solver.getName()
        except Exception:
            name = None
    if name is None:
        name = type(solver).__name__.replace('Solver', '')
    if str(name).upper() == 'LN':
        return _ln_token(solver)
    token = _SOLVER_TOKENS.get(str(name).upper())
    if token is None:
        raise RuntimeError(
            "lang='java' is not available for solver '%s'. Supported: %s."
            % (name, ', '.join(sorted(set(_SOLVER_TOKENS.values()))))
        )
    return token


def _extract_json_object(text):
    """
    Extract the first complete JSON object from CLI stdout, ignoring any
    leading/trailing non-JSON text (e.g. residual solver messages).
    """
    start = text.find('{')
    if start < 0:
        raise ValueError("no JSON object found in solver output")
    obj, _ = json.JSONDecoder().raw_decode(text[start:])
    return obj


def solve_via_jar(solver, analysis='avg', timeout=None, node=None, jclass=None,
                  percentiles=None, events=None, state=None):
    """
    Serialize the solver's model, run jline.jar on it, and return the parsed
    results dict (keyed by analysis type, e.g. {'avg': {...}}).

    Args:
        solver:      a native NetworkSolver instance exposing .model and getName().
        analysis:    CLI analysis type (default 'avg').
        timeout:     optional subprocess timeout in seconds.
        node:        0-based node index for prob/sample analyses (CLI -n).
        jclass:      0-based class index for prob-marg (CLI -c).
        percentiles: list of percentile levels 0..100 for perct-respt.
        events:      number of events for sampling analyses (CLI --events).
        state:       list of state values for prob analyses (CLI --state).
    """
    from ..io.linemodel_io import save_model

    java = find_java()
    jar = find_jar()
    token = _solver_token(solver)
    model = _get_model(solver)

    seed = None
    cutoff = None
    samples = None
    warmupfrac = None
    method = None
    tol = None
    iter_tol = None
    iter_max = None
    opts = getattr(solver, 'options', None)
    if opts is not None:
        seed = getattr(opts, 'seed', None)
        cutoff = _scalar_cutoff(getattr(opts, 'cutoff', None))
        # Simulation / Monte Carlo sample budget matters for every simulation-
        # based token (SSA, LDES, JMT) plus NC's sampling methods and LQNS
        # (lqsim). Forward it so the user's samples option is honored under
        # lang="java" instead of the JAR default. For the DES solvers (SSA,
        # LDES) options.events (the event budget) overrides samples, matching
        # the in-process solvers.
        if token in ('ssa', 'ldes', 'jmt', 'nc', 'lqns'):
            s = getattr(opts, 'samples', None)
            if token in ('ssa', 'ldes'):
                ev = getattr(opts, 'events', None)
                if ev:
                    s = ev
            if s is not None:
                try:
                    samples = int(s)
                except (TypeError, ValueError):
                    samples = None
        # SSA warmup discard fraction: forward so the JAR applies the same
        # mean-estimate and CI transient discard as the native engine.
        if token == 'ssa':
            w = getattr(opts, 'warmupfrac', None)
            try:
                if w is not None and float(w) > 0:
                    warmupfrac = float(w)
            except (TypeError, ValueError):
                warmupfrac = None
        # Solution method (e.g. 'exact', 'amva', 'comom', 'lin'): forward whenever
        # the user picked something other than the generic default, so lang="java"
        # honors the requested algorithm instead of the JAR's default method.
        m = getattr(opts, 'method', None)
        if m is not None:
            m = str(m)
            if m and m != 'default':
                method = m
        # Numeric convergence controls: forward only when they differ from the
        # generic native default, so an untouched option lets the JAR apply its
        # own per-solver-type default (mirroring the wrapper).
        tol = _forward_if_overridden(getattr(opts, 'tol', None), _NATIVE_DEFAULT_TOL)
        iter_tol = _forward_if_overridden(getattr(opts, 'iter_tol', None), _NATIVE_DEFAULT_ITER_TOL)
        im = _forward_if_overridden(getattr(opts, 'iter_max', None), _NATIVE_DEFAULT_ITER_MAX)
        if im is not None:
            try:
                iter_max = int(im)
            except (TypeError, ValueError):
                iter_max = None

    tmpdir = tempfile.mkdtemp(prefix='line_jar_')
    model_path = os.path.join(tmpdir, 'model.json')
    try:
        save_model(model, model_path)

        cmd = [java, '-jar', jar,
               '-f', model_path,
               '-i', 'json',
               '-s', token,
               '-a', analysis,
               '-o', 'json',
               '-v', 'silent']
        if seed is not None:
            cmd += ['-d', str(int(seed))]
        if cutoff is not None:
            cmd += ['--cutoff', repr(float(cutoff))]
        if samples is not None and samples > 0:
            cmd += ['--samples', str(samples)]
        if warmupfrac is not None:
            cmd += ['--warmupfrac', repr(warmupfrac)]
        if method is not None:
            cmd += ['--method', method]
        if tol is not None:
            cmd += ['--tol', repr(float(tol))]
        if iter_tol is not None:
            cmd += ['--iter_tol', repr(float(iter_tol))]
        if iter_max is not None and iter_max > 0:
            cmd += ['--iter_max', str(iter_max)]
        if node is not None:
            cmd += ['-n', str(int(node))]
        if jclass is not None:
            cmd += ['-c', str(int(jclass))]
        if percentiles is not None:
            cmd += ['--percentiles', ','.join(str(float(p)) for p in percentiles)]
        if events is not None:
            cmd += ['--events', str(int(events))]
        if state is not None:
            cmd += ['--state', ','.join(str(int(s)) for s in state)]

        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              timeout=timeout)
        stdout = proc.stdout.decode('utf-8', errors='replace')
        stderr = proc.stderr.decode('utf-8', errors='replace')
        if proc.returncode != 0:
            raise RuntimeError(
                "jline.jar exited with code %d for solver '%s'.\nstderr:\n%s"
                % (proc.returncode, token, stderr.strip())
            )
        try:
            return _extract_json_object(stdout)
        except Exception as e:
            raise RuntimeError(
                "could not parse jline.jar JSON output for solver '%s': %s\n"
                "stdout:\n%s\nstderr:\n%s"
                % (token, e, stdout.strip(), stderr.strip())
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


class _JavaAvgResult:
    """
    Result container populated from JAR avg output. Exposes the station x class
    metric matrices under every field-name convention the native solvers use:
    item access (``r['QN']``, MVA) and attribute access (``.Q/.U/.R/.T``, used by
    NC/CTMC/SSA; ``.QN/.UN/.RN/.TN``, used by MAM/FLD). This lets a single
    JAR-populated result drive every downstream getter (tables, matrices, scalar
    metrics, chain/node aggregations) unchanged.
    """

    def __init__(self, QN, UN, RN, TN, AN, WN, XN, CN):
        self.QN = self.Q = QN
        self.UN = self.U = UN
        self.RN = self.R = RN
        self.TN = self.T = TN
        self.AN = self.A = AN
        self.WN = self.W = WN
        self.XN = self.X = XN   # system throughput per class
        self.CN = self.C = CN   # system response (cycle) time per class
        self._d = {'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'AN': AN, 'WN': WN,
                   'XN': XN, 'CN': CN,
                   'Q': QN, 'U': UN, 'R': RN, 'T': TN, 'A': AN, 'W': WN,
                   'X': XN, 'C': CN}

    # Two- and one-letter aliases that name the same metric matrix. Writing one
    # (e.g. res['UN'] = ... from SolverMVA._cap_unstable_open_util) must update
    # both dict keys and both attributes so every downstream getter sees it.
    _ALIAS = {'QN': ('QN', 'Q'), 'Q': ('QN', 'Q'),
              'UN': ('UN', 'U'), 'U': ('UN', 'U'),
              'RN': ('RN', 'R'), 'R': ('RN', 'R'),
              'TN': ('TN', 'T'), 'T': ('TN', 'T'),
              'AN': ('AN', 'A'), 'A': ('AN', 'A'),
              'WN': ('WN', 'W'), 'W': ('WN', 'W'),
              'XN': ('XN', 'X'), 'X': ('XN', 'X'),
              'CN': ('CN', 'C'), 'C': ('CN', 'C')}

    def __getitem__(self, k):
        return self._d[k]

    def __setitem__(self, k, v):
        for alias in _JavaAvgResult._ALIAS.get(k, (k,)):
            self._d[alias] = v
            setattr(self, alias, v)

    def get(self, k, default=None):
        return self._d.get(k, default)

    def __contains__(self, k):
        return k in self._d

    def keys(self):
        return self._d.keys()

    # Solver-internal structures the native getters read off the result object
    # (CTMC generator/state space, fluid state vector). Fetched lazily from the
    # JAR on first access so the common avg path stays a single JVM call.
    _LAZY_ATTRS = ('infgen', 'space', 'space_aggr', 'spaceAggr', 'nodeSpace',
                   'eventFilt', 'pi', 'xvec')

    def __getattr__(self, name):
        if name in _JavaAvgResult._LAZY_ATTRS:
            solver = self.__dict__.get('_solver')
            if solver is not None:
                cache = self.__dict__.setdefault('_lazy', {})
                if name not in cache:
                    cache.update(_fetch_internal_structures(solver, name))
                if name in cache:
                    return cache[name]
        raise AttributeError(name)


def _fetch_internal_structures(solver, name):
    """Fetch the group of solver-internal structures that includes ``name`` from
    the JAR: the CTMC generator/state space (generator token) or the fluid state
    vector (statevec token). Returns a dict of native-named attributes."""
    if name == 'xvec':
        return statevec_via_jar(solver)
    return generator_via_jar(solver)


def generator_via_jar(solver, timeout=None):
    """Delegate the CTMC infinitesimal generator, state space, stationary
    distribution and event filters to jline.jar. Returns native-named
    attributes: ``infgen`` (Q), ``space`` and ``spaceAggr``/``space_aggr``
    (state spaces), ``pi`` (stationary distribution), ``eventFilt`` (list of
    matrices)."""
    import numpy as np

    res = solve_via_jar(solver, analysis='generator', timeout=timeout)
    g = res.get('generator', res)
    if not isinstance(g, dict) or 'infgen' not in g:
        raise RuntimeError(
            "jline.jar did not return a generator structure (keys: %s)."
            % (list(g.keys()) if isinstance(g, dict) else type(g)))

    def mat(x):
        return np.atleast_2d(np.asarray(x, dtype=float)) if x is not None else None

    # pi is indexed by the generator and state spaces the analyzer solved on: on
    # a reducible chain that is the component supporting pi, not the full chain.
    # The native result keeps all of them aligned, so adopt the work triple
    # (generator, space, aggregated space) whenever the JAR reports it -- ALL
    # THREE or none, never a mix. A full-space generator beside a restricted pi
    # is worse than either alone: pi @ Q cannot even be formed, and any consumer
    # that indexes Q by a state-space row reads the wrong state. On an
    # irreducible chain the work fields are the full ones, so this is a no-op.
    pi = mat(g.get('pi'))
    if pi is not None and 1 in pi.shape:
        pi = pi.reshape(-1)
    work = (g.get('infgen_work'), g.get('space_work'), g.get('space_aggr_work'))
    if all(w is not None for w in work):
        infgen, space, space_aggr = (mat(w) for w in work)
    else:
        infgen = mat(g.get('infgen'))
        space = mat(g.get('space'))
        space_aggr = mat(g.get('space_aggr'))
    return {
        'infgen': infgen,
        'space': space,
        'space_aggr': space_aggr,
        'spaceAggr': space_aggr,
        'nodeSpace': None,
        'eventFilt': [mat(m) for m in (g.get('eventFilt') or [])],
        'pi': pi,
    }


def statevec_via_jar(solver, timeout=None):
    """Delegate the fluid solver's ODE steady-state vector to jline.jar.
    Returns the native-named attribute ``xvec``."""
    import numpy as np

    res = solve_via_jar(solver, analysis='statevec', timeout=timeout)
    g = res.get('statevec', res)
    x = g.get('xvec') if isinstance(g, dict) else None
    xv = np.asarray(x, dtype=float) if x is not None else None
    if xv is not None and xv.ndim == 2 and 1 in xv.shape:
        xv = xv.reshape(-1)
    return {'xvec': xv}


# JAR avg-table column -> native station-matrix field.
_COLMAP = {'QN': 'QLen', 'UN': 'Util', 'RN': 'RespT',
           'TN': 'Tput', 'AN': 'ArvR', 'WN': 'ResidT'}


def _resolve_struct(solver):
    """Get the network struct the solver's getters read, building it if needed."""
    sn = getattr(solver, '_sn', None)
    if sn is None:
        sn = getattr(solver, 'sn', None)
    if sn is None:
        sn = _get_model(solver).getStruct()
    return sn


def station_matrices_via_jar(solver, sn, timeout=None):
    """
    Run the JAR avg analysis and reduce its node-level table to station x class
    metric matrices (QN, UN, RN, TN, AN, WN), indexed exactly as the native
    result containers are.
    """
    import numpy as np

    results = solve_via_jar(solver, analysis='avg', timeout=timeout)
    table = results.get('avg', results)
    if not isinstance(table, dict) or 'QLen' not in table:
        raise RuntimeError(
            "jline.jar did not return a structured AvgTable (got keys: %s)."
            % (list(table.keys()) if isinstance(table, dict) else type(table))
        )

    M = int(sn.nstations)
    K = int(sn.nclasses)
    nodenames = [str(x) for x in list(sn.nodenames)]
    classnames = [str(x) for x in list(sn.classnames)]
    node_to_station = np.asarray(sn.nodeToStation).flatten()

    name_to_station = {}
    for node_idx, nm in enumerate(nodenames):
        st = int(node_to_station[node_idx]) if node_idx < len(node_to_station) else -1
        if st >= 0:
            name_to_station[nm] = st
    class_to_idx = {nm: i for i, nm in enumerate(classnames)}

    mats = {k: np.zeros((M, K)) for k in _COLMAP}
    stations = table['Station']
    jobclasses = table['JobClass']
    for i in range(len(stations)):
        st = name_to_station.get(str(stations[i]))
        cl = class_to_idx.get(str(jobclasses[i]))
        if st is None or cl is None:
            continue
        for field, col in _COLMAP.items():
            v = table[col][i]
            mats[field][st, cl] = 0.0 if v is None else float(v)

    # System throughput per class XN[r] = throughput at class r's reference
    # station; system cycle time CN[r] = sum of station response times.
    refstat = np.asarray(sn.refstat).flatten().astype(int)
    XN = np.zeros(K)
    for r in range(K):
        st = int(refstat[r]) if r < len(refstat) else -1
        if 0 <= st < M:
            XN[r] = mats['TN'][st, r]
    mats['XN'] = XN
    mats['CN'] = mats['RN'].sum(axis=0)
    # The executed method (e.g. 'nrm' vs 'serial') is emitted by the CLI as a
    # sibling of 'avg' for SSA; carry it so populate_java_result can set
    # result.method (native solvers expose the method actually run).
    mats['_method'] = results.get('method') if isinstance(results, dict) else None
    return mats


def populate_java_result(solver, timeout=None):
    """
    Populate a native solver's result container (and struct) from the JAR, so
    every downstream getter returns JAR-derived values. Called from each
    solver's runAnalyzer when options.lang == 'java'.
    """
    import numpy as np

    sn = _resolve_struct(solver)
    for attr in ('_sn', 'sn'):
        try:
            setattr(solver, attr, sn)
        except Exception:
            pass

    # Name attributes some getters read (native runAnalyzer sets these; the
    # java path skips it, so provide them from the struct).
    try:
        classnames = [str(x) for x in list(sn.classnames)]
        node_to_station = np.asarray(sn.nodeToStation).flatten()
        nodenames = [str(x) for x in list(sn.nodenames)]
        station_names = [''] * int(sn.nstations)
        for node_idx, nm in enumerate(nodenames):
            st = int(node_to_station[node_idx]) if node_idx < len(node_to_station) else -1
            if 0 <= st < len(station_names):
                station_names[st] = nm
        solver.class_names = classnames
        solver.station_names = station_names
    except Exception:
        pass

    mats = station_matrices_via_jar(solver, sn, timeout=timeout)
    container = _JavaAvgResult(mats['QN'], mats['UN'], mats['RN'], mats['TN'],
                               mats['AN'], mats['WN'], mats['XN'], mats['CN'])
    # Method actually executed by the JAR (SSA reports 'nrm' vs 'serial'); native
    # solvers expose result.method, and tests assert it to catch a silent serial
    # fallback. None for solvers/CLIs that do not emit it.
    container.__dict__['method'] = mats.get('_method')
    # Reference kept so solver-internal structures (CTMC generator/state space,
    # fluid state vector) can be fetched lazily from the JAR on first access.
    container.__dict__['_solver'] = solver
    # Cache nodes: the hit/miss class-switch throughput lives only on the Cache
    # node's output edge (never at a station), so the station-level avg table the
    # JAR returns cannot carry it. Seed each Cache node's actual hit/miss
    # probabilities from the JAR cache table so the native getAvgNode() derives
    # the hit/miss node throughputs (arvTput*prob) exactly as in lang='python'.
    _inject_cache_hitprob_via_jar(solver, sn, timeout=timeout)
    for attr in ('_result', 'result'):
        try:
            setattr(solver, attr, container)
        except Exception:
            pass
    return container


def _inject_cache_hitprob_via_jar(solver, sn, timeout=None):
    """
    Seed each Cache node's actualhitprob/actualmissprob (per read class) from the
    JAR cache table, so the native getAvgNode()/sn_get_node_tput_from_tput derives
    the Cache-node hit/miss class throughputs (arvTput * prob) under lang='java'
    exactly as under lang='python'. The JAR station-level avg table cannot carry
    them because the hit/miss class-switch throughput exists only on the Cache
    node's output edge. No-op when the model has no Cache node.
    """
    import numpy as np
    from ..api.sn.network_struct import NodeType

    nodetype = getattr(sn, 'nodetype', None)
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodetype is None or nodeparam is None:
        return
    cache_inds = [ind for ind in range(int(sn.nnodes))
                  if ind < len(nodetype) and nodetype[ind] == NodeType.CACHE]
    if not cache_inds:
        return

    try:
        df = cache_table_via_jar(solver, timeout=timeout)
    except Exception:
        return
    if df is None or not hasattr(df, 'iterrows'):
        return

    classnames = [str(x) for x in list(sn.classnames)]
    nodenames = [str(x) for x in list(sn.nodenames)]
    class_to_idx = {nm: i for i, nm in enumerate(classnames)}
    name_to_node = {nm: i for i, nm in enumerate(nodenames)}
    R = int(sn.nclasses)

    # (node_ind, read_class_idx) -> (hitprob, missprob) from the JAR cache table.
    for _, row in df.iterrows():
        ind = name_to_node.get(str(row.get('Node')))
        r = class_to_idx.get(str(row.get('JobClass')))
        if ind is None or r is None or ind not in cache_inds:
            continue
        try:
            hp = float(row.get('HitProb'))
            mp = float(row.get('MissProb'))
        except (TypeError, ValueError):
            continue
        if not np.isfinite(hp):
            continue
        # nodeparam is a dict keyed by node index (cache/router/etc.), not a
        # dense per-node list, so index by key rather than position.
        if isinstance(nodeparam, dict):
            cp = nodeparam.get(ind)
        else:
            cp = nodeparam[ind] if ind < len(nodeparam) else None
        if cp is None:
            continue
        ahp = getattr(cp, 'actualhitprob', None)
        if ahp is None or not hasattr(ahp, '__len__') or len(ahp) < R:
            ahp = np.zeros(R)
        else:
            ahp = np.asarray(ahp, dtype=float).copy()
        amp = getattr(cp, 'actualmissprob', None)
        if amp is None or not hasattr(amp, '__len__') or len(amp) < R:
            amp = np.zeros(R)
        else:
            amp = np.asarray(amp, dtype=float).copy()
        if r < R:
            ahp[r] = hp
            amp[r] = mp if np.isfinite(mp) else (1.0 - hp)
        cp.actualhitprob = ahp
        cp.actualmissprob = amp


# =========================================================================
# Specialized analyses (distribution / probability / reward), delegated to the
# JAR CLI via the same subprocess+JSON transport. Each returns the native
# method's contract so the calling getter can return it unchanged.
# =========================================================================

def _parse_distribution(results, token):
    """Parse a JAR DistributionResult into the native List[Dict] contract:
    one dict per (station, class) with a non-empty CDF, keys 1-based
    'station'/'class' and numpy 't' (times) and 'p' (CDF values)."""
    import numpy as np

    dr = results.get(token, results)
    if not isinstance(dr, dict) or 'cdf' not in dr:
        raise RuntimeError(
            "jline.jar did not return a structured DistributionResult (keys: %s)."
            % (list(dr.keys()) if isinstance(dr, dict) else type(dr)))
    out = []
    for i, row in enumerate(dr['cdf']):
        for j, cell in enumerate(row):
            t = np.asarray(cell.get('t', []), dtype=float)
            p = np.asarray(cell.get('p', []), dtype=float)
            if t.size == 0:
                continue
            out.append({'station': i + 1, 'class': j + 1, 't': t, 'p': p})
    return out


def cdf_respt_via_jar(solver, timeout=None):
    """Delegate getCdfRespT() to jline.jar (response-time CDF)."""
    return _parse_distribution(
        solve_via_jar(solver, analysis='cdf-respt', timeout=timeout), 'cdf-respt')


def cdf_passt_via_jar(solver, timeout=None):
    """Delegate getCdfPassT() to jline.jar (passage-time CDF)."""
    return _parse_distribution(
        solve_via_jar(solver, analysis='cdf-passt', timeout=timeout), 'cdf-passt')


def _station_to_stateful(solver, station0):
    """Map a 0-based station index to the 0-based stateful-node index the JAR
    CLI's -n option expects. Falls back to the station index when a direct
    mapping is unavailable (true when every stateful node is a station)."""
    import numpy as np
    sn = getattr(solver, '_sn', None)
    if sn is None:
        sn = getattr(solver, 'sn', None)
    if sn is None or station0 is None:
        return station0
    try:
        node = int(np.asarray(sn.stationToNode).flatten()[station0])
        nts = np.asarray(sn.nodeToStateful).flatten()
        sf = int(nts[node])
        if sf >= 0:
            return sf
    except Exception:
        pass
    return station0


def _node_to_stateful(solver, node0):
    """Map a 0-based node index to the 0-based stateful-node index the JAR
    CLI's -n option expects. Falls back to the node index when a direct mapping
    is unavailable (true when every stateful node is also the leading node)."""
    import numpy as np
    sn = getattr(solver, '_sn', None)
    if sn is None:
        sn = getattr(solver, 'sn', None)
    if sn is None or node0 is None:
        return node0
    try:
        sf = int(np.asarray(sn.nodeToStateful).flatten()[node0])
        if sf >= 0:
            return sf
    except Exception:
        pass
    return node0


def prob_via_jar(solver, token, ist=None, jclass=None, kind='scalar',
                 onebased=False, raw_station=False, timeout=None):
    """
    Delegate a probability query to jline.jar and return it in the native
    contract. `kind` selects the shape: 'scalar' -> float(probability),
    'logtuple' -> (logProbability, probability), 'vector' -> probabilityVector.
    `ist`/`jclass` are the native station/class indices (converted using the
    solver's 0/1-based convention via `onebased`). `raw_station=True` passes the
    0-based station index directly (used by prob-marg, whose JAR path treats -n
    as a station index rather than a stateful-node index).
    """
    import numpy as np

    node = None
    if ist is not None:
        station0 = int(ist) - (1 if onebased else 0)
        node = station0 if raw_station else _station_to_stateful(solver, station0)
    jc = None
    if jclass is not None:
        jc = int(jclass) - (1 if onebased else 0)

    results = solve_via_jar(solver, analysis=token, node=node, jclass=jc,
                            timeout=timeout)
    pr = results.get(token, results)
    if not isinstance(pr, dict):
        raise RuntimeError(
            "jline.jar did not return a ProbabilityResult for '%s' (got %s)."
            % (token, type(pr)))
    scalar = float(pr.get('scalar', 0.0))
    if kind == 'scalar':
        return scalar
    if kind == 'logtuple':
        logp = float(np.log(scalar)) if scalar > 0 else float('-inf')
        return (logp, scalar)
    if kind == 'vector':
        return np.asarray(pr.get('probability', []), dtype=float).ravel()
    raise ValueError("unknown prob kind: %s" % kind)


_CACHE_COLS = ['Node', 'JobClass', 'List', 'ListCap', 'Items', 'HitProb',
               'DelayedHitProb', 'MissProb', 'HitRate', 'DelayedHitRate',
               'MissRate', 'ArvR', 'ResidT']
_ITEM_COLS = ['Node', 'Item', 'List', 'ListCap', 'Prob']
# Columns of the cache/item tables that carry names rather than metrics.
_TABLE_STRING_COLS = ('Node', 'JobClass')


def _table_via_jar(solver, token, cols, timeout=None):
    """Run a JAR table analysis and rebuild a DataFrame with the given columns."""
    import pandas as pd
    results = solve_via_jar(solver, analysis=token, timeout=timeout)
    t = results.get(token, results)
    if not isinstance(t, dict) or cols[0] not in t:
        raise RuntimeError(
            "jline.jar did not return a structured '%s' table (keys: %s)."
            % (token, list(t.keys()) if isinstance(t, dict) else type(t)))
    n = len(t.get(cols[0], []))
    rows = [{c: t[c][i] for c in cols} for i in range(n)]
    df = pd.DataFrame(rows, columns=cols)
    # Coerce the metric columns to float so a JSON null becomes NaN, as in the
    # native tables. Left alone, a null arrives as None and pandas types the whole
    # column as object, which silently breaks arithmetic on the result.
    for c in cols:
        if c not in _TABLE_STRING_COLS:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


_LN_COLS = ['Node', 'NodeType', 'QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']


def ln_avg_table_via_jar(solver, timeout=None):
    """
    Delegate a layered-network (LQN) getAvgTable() to jline.jar. Returns a
    DataFrame with one row per LQN element (host/task/entry/activity) and the
    native layered columns (Node, NodeType, QLen, Util, RespT, ResidT, ArvR,
    Tput). Used by the native SolverLN/SolverLQNS lang="java" dispatch: the
    LayeredNetwork is serialized to JSON and solved by the canonical JAR.
    """
    import pandas as pd

    results = solve_via_jar(solver, analysis='avg', timeout=timeout)
    table = results.get('avg', results)
    if not isinstance(table, dict) or 'Node' not in table:
        raise RuntimeError(
            "jline.jar did not return a structured layered AvgTable (keys: %s)."
            % (list(table.keys()) if isinstance(table, dict) else type(table)))
    # The layered fixed point may stop on iter_max without converging, so mirror
    # the JAR's flag onto the native solver: delegating must not silently report
    # convergence that the underlying solve never reached.
    if isinstance(results, dict) and 'hasconverged' in results:
        solver.hasconverged = bool(results['hasconverged'])
    n = len(table.get('Node', []))
    rows = [{c: table[c][i] for c in _LN_COLS if c in table} for i in range(n)]
    df = pd.DataFrame(rows, columns=_LN_COLS)
    # Coerce numeric metric columns to float so JSON nulls become NaN (matching
    # the native layered table), leaving the string Node/NodeType columns intact.
    for c in ('QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput'):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def ln_ensemble_avg_via_jar(solver, timeout=None):
    """
    Delegate a layered-network (LQN) getEnsembleAvg() to jline.jar, returning
    (QN, UN, RN, TN, AN, WN) as (nidx+1,) arrays indexed by LQN absolute index,
    matching the native SolverLN.get_ensemble_avg contract.

    This solves the whole LayeredNetwork in ONE JAR subprocess. Without it, a
    lang="java" ensemble solve falls through to the native LN fixed-point loop
    whose per-layer sub-solves each dispatch to their own `java -jar` process
    (hundreds of JVM startups per LN solve).
    """
    df = ln_avg_table_via_jar(solver, timeout=timeout)
    lqn = solver.lqn
    nidx = lqn.nidx
    QN, UN, RN, TN, AN, WN = (np.full(nidx + 1, np.nan) for _ in range(6))
    # The JAR's Node column carries the same hashnames the native layered table
    # builds, so index by name rather than relying on row order.
    name2idx = {solver._get_hashname(idx): idx for idx in range(1, nidx + 1)}
    for row in df.itertuples(index=False):
        idx = name2idx.get(row.Node)
        if idx is None:
            continue
        QN[idx] = row.QLen
        UN[idx] = row.Util
        RN[idx] = row.RespT
        TN[idx] = row.Tput
        AN[idx] = row.ArvR
        WN[idx] = row.ResidT
    return QN, UN, RN, TN, AN, WN


def cache_table_via_jar(solver, timeout=None):
    """Delegate getAvgCacheTable() to jline.jar."""
    return _table_via_jar(solver, 'cache', _CACHE_COLS, timeout=timeout)


def item_table_via_jar(solver, timeout=None):
    """Delegate getAvgItemTable() to jline.jar."""
    return _table_via_jar(solver, 'item', _ITEM_COLS, timeout=timeout)


def sample_via_jar(solver, node, num_events=None, timeout=None):
    """
    Delegate SSA sample() to jline.jar. Returns a SamplePath with the ARV/DEP
    event trace (and t/state), matching the native contract. Node is the native
    1-based node argument; the JAR CLI -n is a 0-based stateful-node index.
    """
    import numpy as np
    from .solver_ssa.solver_ssa import SamplePath, SampleEvent

    # Resolve a 0-based *node* index from the native node argument (an int is
    # the 1-based node index; a node object exposes get_index0()/get_index()).
    node0 = None
    if node is not None:
        if isinstance(node, (int, np.integer)):
            node0 = int(node) - 1
        else:
            v = getattr(node, 'get_index0', getattr(node, 'getIndex0', None))
            if v is not None:
                node0 = int(v() if callable(v) else v)
            else:
                for attr in ('get_index', 'getIndex', 'index'):
                    v = getattr(node, attr, None)
                    if v is not None:
                        node0 = int(v() if callable(v) else v) - 1
                        break

    # The JAR CLI -n indexes model.getStatefulNodes() (0-based stateful index).
    node_idx = _node_to_stateful(solver, node0)

    ev = int(num_events) if num_events else 10000
    results = solve_via_jar(solver, analysis='sample', node=node_idx,
                            events=ev, timeout=timeout)
    sns = results.get('sample', results)
    sp = SamplePath()
    sp.handle = node
    if isinstance(sns, dict):
        sp.t = np.asarray(sns.get('t', []), dtype=float)
        sp.state = np.asarray(sns.get('state', []), dtype=float)
        sp.isaggregate = bool(sns.get('isaggregate', False))
        events = []
        for e in sns.get('event', []) or []:
            et = e.get('event')
            if et in ('ARV', 'DEP'):
                events.append(SampleEvent(t=float(e['t']),
                                          node=int(e['node']) + 1,
                                          class_idx=int(e['jobclass']) + 1,
                                          event=et))
        sp.event = events
    return sp


def tran_avg_via_jar(solver, timeout=None):
    """
    Delegate getTranAvg() to jline.jar. Returns (QNt, UNt, TNt), each a nested
    [M][K] list of TranResult(t, values), matching the native fluid/CTMC
    transient contract.
    """
    import numpy as np
    from ..constants import TranResult

    results = solve_via_jar(solver, analysis='tran-avg', timeout=timeout)
    tr = results.get('tran-avg', results)
    if not isinstance(tr, dict) or 't' not in tr:
        raise RuntimeError(
            "jline.jar did not return a transient result (keys: %s)."
            % (list(tr.keys()) if isinstance(tr, dict) else type(tr)))

    t = np.asarray(tr.get('t', []), dtype=float).ravel()

    def _grid(key):
        g = tr.get(key, []) or []
        out = []
        for row in g:
            cells = []
            for cell in row:
                a = np.asarray(cell, dtype=float)
                if a.ndim == 2 and a.shape[1] >= 1:
                    vals = a[:, 0]           # column 0 is the metric value
                elif a.ndim == 2:
                    vals = a.ravel()
                else:
                    vals = a
                cells.append(TranResult(t, vals))
            out.append(cells)
        return out

    return _grid('QNt'), _grid('UNt'), _grid('TNt')


def perct_via_jar(solver, percentiles=None, timeout=None):
    """
    Delegate getPerctRespT() to jline.jar (fork-join response-time percentiles).
    Returns (PercRT_list, DataFrame) matching the native contract: a list of
    per-fork-join dicts and a [Percentile, RespT] table.
    """
    import numpy as np
    import pandas as pd

    if percentiles is None:
        percentiles = [50, 75, 90, 95, 99]
    results = solve_via_jar(solver, analysis='perct-respt',
                            percentiles=percentiles, timeout=timeout)
    payload = results.get('perct-respt', results)
    items = payload if isinstance(payload, list) else [payload]

    perc_rt = []
    rows = []
    for it in items:
        if not isinstance(it, dict):
            continue
        ps = np.asarray(it.get('percentiles', []), dtype=float)
        rt = np.asarray(it.get('RTp', []), dtype=float)
        perc_rt.append({'station': 'ForkJoin', 'class': 1,
                        'percentiles': ps, 'values': rt})
        for p, v in zip(ps, rt):
            rows.append({'Percentile': p, 'RespT': v})
    return perc_rt, pd.DataFrame(rows, columns=['Percentile', 'RespT'])


def avg_reward_via_jar(solver, timeout=None):
    """
    Delegate getAvgReward() to jline.jar (CTMC). Returns (values, names) where
    values is a numpy array of steady-state reward values and names the reward
    names, matching the native getAvgReward contract.
    """
    import numpy as np

    results = solve_via_jar(solver, analysis='reward-steady', timeout=timeout)
    rewards = results.get('reward-steady', results)
    if not isinstance(rewards, dict):
        raise RuntimeError(
            "jline.jar did not return a reward map (got %s)." % type(rewards))
    names = list(rewards.keys())
    values = np.array([float(rewards[n]) for n in names])
    return values, names
