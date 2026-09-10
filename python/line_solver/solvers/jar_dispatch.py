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
import platform
import shutil
import re
import subprocess
import tempfile

import numpy as np


# Map a native solver's getName() to the CLI solver method name accepted by
# jline.jar's LineCLI (-s option). Method names must match validateSolver().
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
    # `ag` and `ba` became LineCLI method names on 2026-08-27; SolverAG (RCAT/INAP) and
    # SolverBA (bounds) had existed in the JAR all along, and only the CLI's
    # validSolvers list stood between them and this bridge. Until then the
    # refusal here was correct, and the five `ag_*` parity rows skipped by name
    # ("lang='java' is not available for solver 'SolverAG'").
    'AG': 'ag',
    'BA': 'ba',
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
# Native CTMC transient output step (solver_ctmc.py: timestep or 0.1). It is NOT
# a generic default: MATLAB and the JAR leave the grid to the adaptive ODE
# solver, so lang="java" must ask for the uniform grid explicitly or the curve
# comes back on ODE nodes that no caller-chosen query point lands on.
_NATIVE_DEFAULT_CTMC_TIMESTEP = 0.1

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
        # java.exe on Windows: a bare "java" is neither a file nor executable
        # there, so this branch used to be skipped with a valid JDK installed.
        exe_name = 'java.exe' if platform.system() == 'Windows' else 'java'
        cand = os.path.join(java_home, 'bin', exe_name)
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


def _cutoff_arg(cutoff):
    """
    Return the ``--cutoff`` argument to forward to the JAR CLI, or None.

    The native CTMC/SSA ``options.cutoff`` is a scalar (int/float) or a
    per-(station,class) array; both are forwarded, a scalar as one number and an
    array as the CLI's ``';'``-separated rows of ``','``-separated cells. None,
    an empty array and any non-finite entry are left to the JAR's own default,
    which is the only case where it has no cutoff of the caller's to honour.

    A NON-UNIFORM ARRAY USED TO BE DROPPED HERE, and dropping it is not
    neutral: the JAR then truncated at its own default 10 and reported that
    answer under the caller's name. oqn_cs_routing pins [1,1,0;3,3,0;0,0,3] and
    the two truncations disagree by 9% on Tput and 32% on QLen -- a difference
    in the QUESTION, which no tolerance should absorb.
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
    for v in vals:
        if not _math.isfinite(v):
            return None
    if arr.ndim < 2 or all(v == vals[0] for v in vals):
        return repr(float(vals[0]))
    return ';'.join(','.join(repr(float(v)) for v in row) for row in _np.atleast_2d(arr))


def _finite_timespan(timespan):
    """
    Return the (T0, T1) pair to forward to the JAR CLI, or None.

    A transient analysis is the caller's time span, and the CLI has no other way
    to learn it: left unset the JAR falls back to its [Inf, Inf] default and
    ``getTranAvg`` substitutes ``30/minRate``, which is a different analysis.
    Only a finite pair is forwarded; an unset or infinite span is left to the
    JAR so its own substitution still applies.
    """
    if timespan is None:
        return None
    try:
        t0 = float(timespan[0])
        t1 = float(timespan[1])
    except (TypeError, ValueError, IndexError):
        return None
    if not _math.isfinite(t0) or not _math.isfinite(t1):
        return None
    return (t0, t1)


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
    Resolve the CLI method name for a layered solve from the layer solver actually in
    use. LineCLI spells the layer solver in the method name: 'ln.nc' (and its alias
    'ln.comom') builds NC layers, 'ln.mva' builds MVA ones, and BARE 'ln' IS
    MVA, not NC -- `buildLayeredSolver` sends it to SolverLN's own default. The
    native SolverLN likewise defaults to MVA and uses NC only when the caller's
    factory produces it. Substituting the layer solver changes results, not just
    the label: an RR cache layer is exact under MVA's 'fpi' and approximate
    under NC's 'spm' (lcq_threehosts, hit 0.5 against 0.48331), so the token
    must name the layer solver rather than lean on a default.
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
        # No layers built yet -- the normal case here, since a delegated solve
        # never enters iterate(). Ask the solver to name what its FACTORY
        # produces before falling back to the native default: reading MVA off an
        # empty list turned every `LN(model, lambda m: NC(m, opts))` into an
        # MVA-layered solve (lcq_threehosts, cache hit 0.5 against 0.48331).
        probe = getattr(solver, 'probe_layer_solver_name', None)
        if probe is not None:
            layer_name = probe()
    if layer_name is None:
        layer_name = 'MVA'
    return 'ln.nc' if str(layer_name).upper() in ('NC', 'COMOM') else 'ln.mva'


def _solver_token(solver):
    """Resolve the CLI solver method name from a native solver instance.

    `getName()` IS NOT UNIFORM across the native solvers: SolverMVA answers
    'MVA' while SolverAG and SolverBA answer 'SolverAG' and 'SolverBA'. The
    fallback below strips the prefix and the getName path did not, so a solver
    with the second spelling could never match this map however it was keyed --
    which read as "lang='java' is not available for solver 'SolverAG'", a
    message about the bridge that named the inconsistency instead. Normalised
    here, exactly as `cpp_dispatch._cpp_solver_name` already does, so the map
    stays keyed by the solver rather than by that quirk.
    """
    name = None
    if hasattr(solver, 'getName'):
        try:
            name = solver.getName()
        except Exception:
            name = None
    if name is None:
        name = type(solver).__name__
    name = str(name)
    if name.upper().startswith('SOLVER') and len(name) > 6:
        name = name[6:]
    if str(name).upper() == 'LN':
        return _ln_token(solver)
    token = _SOLVER_TOKENS.get(str(name).upper())
    if token is None:
        raise RuntimeError(
            "lang='java' is not available for solver '%s'. Supported: %s."
            % (name, ', '.join(sorted(set(_SOLVER_TOKENS.values()))))
        )
    return token


# The JAR logs through java.util.logging with a SimpleLineFormatter, i.e. lines of
# the form "WARNING: [caller] message" on STDERR. solve_via_jar keeps only the JSON
# object from stdout, so without this relay every warning the JAR raises is dropped
# and a lang=java solve looks clean where the native solve warns. That matters
# because some of these warnings report a CORRECTNESS limitation (e.g. a
# matrix-exponential service answered by a phase-type approximation), so silently
# discarding them misrepresents the result rather than merely losing a log line.
_JAR_LOG_RE = re.compile(r'^\s*(WARNING|SEVERE):\s*\[([^\]]+)\]\s*(.*)$')


def _relay_jar_logs(stdout, stderr):
    """Re-emit the JAR's warnings through the native logger."""
    from ..api.io.logging import line_warning_always
    for stream in (stderr, stdout):
        if not stream:
            continue
        for line in stream.splitlines():
            m = _JAR_LOG_RE.match(line)
            if m is None:
                continue
            # '%s' with the text as an argument, never the text as the format
            # string: a stray % in a JAR message would otherwise raise here.
            line_warning_always(m.group(2), '%s', m.group(3))


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
    timespan = None
    timestep = None
    samples = None
    warmupfrac = None
    method = None
    tol = None
    iter_tol = None
    iter_max = None
    multiserver = None
    map_env = None
    map_env_method = None
    pstar = None
    opts = getattr(solver, 'options', None)
    if opts is not None:
        seed = getattr(opts, 'seed', None)
        cutoff = _cutoff_arg(getattr(opts, 'cutoff', None))
        timespan = _finite_timespan(getattr(opts, 'timespan', None))
        if token == 'ctmc':
            ts = getattr(opts, 'timestep', None)
            timestep = float(ts) if ts else _NATIVE_DEFAULT_CTMC_TIMESTEP
        # Simulation / Monte Carlo sample budget matters for every simulation-
        # based method name (SSA, LDES, JMT) plus NC's sampling methods and LQNS
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
        # AMVA multiserver rule: it decides WHICH algorithm serves a multiserver
        # model (Seidmann transform, solver_amvald, conway, krzesinski), so a
        # delegated solve that never receives it answers under the default rule
        # while the caller believes its own was applied -- an agreement between the
        # langs that is an artifact of the option being dropped.
        _cfg = getattr(opts, 'config', None)
        _ms = None
        if isinstance(_cfg, dict):
            _ms = _cfg.get('multiserver')
        elif _cfg is not None:
            _ms = getattr(_cfg, 'multiserver', None)
        if _ms is not None and str(_ms) and str(_ms) != 'default':
            multiserver = str(_ms)
        # Non-renewal random-environment gate. map_env='off' asks for the model to
        # be REJECTED rather than approximated; a delegated solve that never
        # receives it runs the JAR's own environment fallback and answers, so the
        # same request refuses under lang='python' and returns numbers under
        # lang='java'.
        _me = _cfg.get('map_env') if isinstance(_cfg, dict) else getattr(_cfg, 'map_env', None)
        if _me is not None and str(_me) and str(_me) != 'auto':
            map_env = str(_me)
        _mem = (_cfg.get('map_env_method') if isinstance(_cfg, dict)
                else getattr(_cfg, 'map_env_method', None))
        if _mem is not None and str(_mem) and str(_mem) != 'auto':
            map_env_method = str(_mem)
        # Fluid p-norm smoothing exponent: it selects the DRIFT the matrix method
        # integrates, so a delegated solve that never receives it returns the
        # unsmoothed mean-field fixed point (lambda for an M/M/1) where the caller
        # asked for the smoothed one.
        _ps = getattr(opts, 'pstar', None)
        if _ps is None:
            _ps = (_cfg.get('pstar') if isinstance(_cfg, dict)
                   else getattr(_cfg, 'pstar', None))
        if _ps is not None:
            try:
                if np.isscalar(_ps):
                    pstar = [float(_ps)]
                else:
                    pstar = [float(p) for p in np.asarray(_ps).ravel()]
                if not pstar:
                    pstar = None
            except (TypeError, ValueError):
                pstar = None
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
               # NOT 'silent': that maps to Level.OFF in the JAR's configureLogger,
               # which discards its warnings outright and leaves a lang=java solve
               # looking clean where the native one warns. Any non-silent level sets
               # the handler to Level.WARNING, so this emits warnings and errors ONLY
               # -- no INFO chatter -- and java.util.logging writes them to stderr,
               # so the JSON on stdout stays clean. _relay_jar_logs re-emits them.
               '-v', 'normal']
        if seed is not None:
            cmd += ['-d', str(int(seed))]
        if cutoff is not None:
            cmd += ['--cutoff', cutoff]
        if timespan is not None:
            cmd += ['--timespan', '%r,%r' % (timespan[0], timespan[1])]
        if timestep is not None and timestep > 0:
            cmd += ['--timestep', repr(float(timestep))]
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
        if multiserver is not None:
            cmd += ['--multiserver', multiserver]
        if map_env is not None:
            cmd += ['--map-env', map_env]
        if map_env_method is not None:
            cmd += ['--map-env-method', map_env_method]
        if pstar is not None:
            cmd += ['--pstar', ','.join(repr(float(p)) for p in pstar)]
        if node is not None:
            cmd += ['-n', str(int(node))]
        if jclass is not None:
            cmd += ['-c', str(int(jclass))]
        if percentiles is not None:
            cmd += ['--percentiles', ','.join(str(float(p)) for p in percentiles)]
        if events is not None:
            cmd += ['--events', str(int(events))]
        if state is not None:
            # A system state is a LIST OF ROWS, one per station; a node state is
            # one flat row. The CLI spells the first with ';' between rows.
            rows = state if state and isinstance(state[0], (list, tuple)) else [state]
            cmd += ['--state',
                    ';'.join(','.join(str(int(v)) for v in row) for row in rows)]

        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              timeout=timeout)
        stdout = proc.stdout.decode('utf-8', errors='replace')
        stderr = proc.stderr.decode('utf-8', errors='replace')
        if proc.returncode != 0:
            raise RuntimeError(
                "jline.jar exited with code %d for solver '%s'.\nstderr:\n%s"
                % (proc.returncode, token, stderr.strip())
            )
        _relay_jar_logs(stdout, stderr)
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
                   'eventFilt', 'pi', 'xvec', 'station_col_ranges')

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
    """Fetch the group of solver-internal structures that includes ``name``: the
    CTMC generator/state space (generator token) or the fluid state vector
    (statevec token). Returns a dict of native-named attributes.

    IT DISPATCHES ON `options.lang`, and must. This container is shared by both
    delegating backends, so a lang='cpp' solver whose lazy attribute was fetched
    here unconditionally returned JAR-derived numbers: `SolverCTMC(m,
    lang='cpp').getInfGen()` answered with jline.jar's generator, under a lang
    that is an assertion about which engine produced it. The fluid state vector
    goes to `-a statevec` for the same reason -- it used to be refused as having
    no C++ arm, which stopped being true when that arm landed, and the stale
    refusal was what kept getCdfRespT and the percentile getters unreachable
    under lang='cpp'."""
    lang = getattr(getattr(solver, 'options', None), 'lang', 'python')
    if str(lang) == 'cpp':
        from .cpp_dispatch import generator_via_cpp, statevec_via_cpp
        if name == 'xvec':
            return statevec_via_cpp(solver)
        return generator_via_cpp(solver)
    if name == 'xvec':
        return statevec_via_jar(solver)
    return generator_via_jar(solver)


def env_avg_via_jar(solver, timeout=None):
    """
    Delegate a random-environment solve to `jline.cli.LineCLI -s env`, returning
    (QN, UN, TN) as (nstations x nclasses) arrays in the STAGE model's indexing.

    THE WHOLE ENVIRONMENT CROSSES, not one stage at a time, and that is the
    point. Running the python coupling with JAR stage solves sends each stage's
    ENTRY MARGINAL over `model.json`, which carries an `initialState` for a
    Place, a PAS queue and a Cache node only -- so every stage started from the
    default state and the loop converged to the UNCOUPLED answer: on
    renv_node_breakdown, Server QLen 0.39634 at throughput 0.71574, the latter
    refuting itself against a source admitting 0.8. Delegating the coupling
    itself leaves nothing to marshal per iteration. Twin of
    `cpp_dispatch.env_avg_via_cpp`, which does the same for `line-cli`.

    RN and AN are not returned: `SolverENV.getEnsembleAvg` defines them as NaN
    because the environment analyzer computes no response time, and the CLI
    emits nan in those columns for the same reason.
    """
    from ..io.linemodel_io import save_model

    java = find_java()
    jar = find_jar()
    env = getattr(solver, 'env_model', None)
    if env is None:
        raise RuntimeError("lang='java' needs the SolverENV's Environment; this solver carries none.")
    stage = None
    for m in list(getattr(solver, 'ensemble', []) or []):
        if m is not None:
            stage = m
            break
    if stage is None:
        raise RuntimeError("lang='java' needs a stage model to index the returned matrices by; "
                           "this Environment has no non-empty stage.")

    opts = getattr(solver, 'options', None) or {}
    method = opts.get('method') if isinstance(opts, dict) else getattr(opts, 'method', None)
    timespan = solver.stage_timespan()

    tmpdir = tempfile.mkdtemp(prefix='line_jar_env_')
    model_path = os.path.join(tmpdir, 'model.json')
    try:
        save_model(env, model_path)
        cmd = [java, '-cp', jar, 'jline.cli.LineCLI',
               '-f', model_path, '-i', 'json', '-s', 'env', '-a', 'avg', '-o', 'json']
        if method is not None and str(method) and str(method) != 'default':
            cmd += ['--method', str(method)]
        # The STAGE horizon, where the caller states a finite one: SolverENV
        # couples transients, so a stage left to pick its own horizon is a
        # different computation (LineCLI defaults it to [0,100] for exactly
        # that reason). It is read off the STAGE SOLVER, which is where an
        # example states it; see SolverENV.stage_timespan.
        ts = timespan
        if ts is not None:
            cmd += ['--timespan', '%r,%r' % (float(ts[0]), float(ts[1]))]
        # WHICH SOLVER RUNS EACH STAGE, which model.json does not carry: the
        # envelope holds the stage Networks and the transition process, and the
        # ensemble's solver choice lives on the SolverENV. The JAR engine already
        # handled a non-fluid stage (`roundMarginalForDiscreteSolver` runs for
        # exactly that case); its CLI hardwired the fluid factory until
        # `--stage-solver` landed beside it.
        from .cpp_dispatch import _env_stage_solver
        st = _env_stage_solver(solver)
        if st is not None:
            cmd += ['--stage-solver', st]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              timeout=timeout)
        stdout = proc.stdout.decode('utf-8', errors='replace')
        stderr = proc.stderr.decode('utf-8', errors='replace')
        if proc.returncode != 0:
            raise RuntimeError(
                "jline.jar exited with code %d for solver 'env'.\nstderr:\n%s"
                % (proc.returncode, stderr.strip()))
        _relay_jar_logs(stdout, stderr)
        results = _extract_json_object(stdout)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    sn = stage.getStruct() if hasattr(stage, 'getStruct') else stage.get_struct()
    mats = avg_matrices_from_results(results, sn, source='jline.jar')
    return mats['QN'], mats['UN'], mats['TN']


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
        'station_col_ranges': station_col_ranges_from_widths(solver, g.get('spaceWidths')),
    }


def station_col_ranges_from_widths(solver, widths):
    """Per-station (start, end) column ranges of a state row, from the widths of
    the per-stateful-node blocks a delegated engine reports.

    A state row is the stateful nodes' blocks laid end to end, and only the
    widths say where one ends. Without them the native getters fall back to one
    column per (station, class), which is right only when every block happens to
    be one column wide: a 2-server FCFS station carries a buffer column AND a
    service column for a single class, and the fallback returned the buffer
    alone, so tut03_repairmen reported the station state as [1] where it is
    [1 waiting, 2 in service]. Returns None when the engine sent no widths.
    """
    import numpy as np

    if widths is None:
        return None
    sn = _resolve_struct(solver)
    if sn is None:
        return None
    M = int(sn.nstations)
    ranges = [(0, 0)] * M
    col = 0
    node_to_stateful = np.asarray(sn.nodeToStateful).reshape(-1)
    node_to_station = np.asarray(sn.nodeToStation).reshape(-1)
    for isf, width in enumerate(widths):
        width = int(width)
        node = next((n for n in range(int(sn.nnodes))
                     if int(node_to_stateful[n]) == isf), -1)
        ist = int(node_to_station[node]) if node >= 0 else -1
        if 0 <= ist < M:
            ranges[ist] = (col, col + width)
        col += width
    return ranges


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


def moments_via_jar(solver, timeout=None):
    """Delegate the moment-closure second-order report to jline.jar, in the
    field layout `SolverFLD.getMoments` returns natively.

    Returns None when the delegated solve ran a first-order method, which
    computes no second moment at all -- the same answer the native getter gives.
    The JAR flattens the per-(station,class) index blocks as ``i*K+k``; they are
    folded back to the nested native shape here, because a caller reading
    ``classBlock[i][k]`` must not have to know the transport's layout.
    """
    import numpy as np

    res = solve_via_jar(solver, analysis='moments', timeout=timeout)
    m = res.get('moments', res)
    if not isinstance(m, dict) or m.get('Sigma') is None:
        return None

    def arr(x, ravel=False):
        if x is None:
            return None
        a = np.asarray(x, dtype=float)
        return a.reshape(-1) if ravel else np.atleast_2d(a)

    sigma = arr(m.get('Sigma'))
    qvar = arr(m.get('QVar'))
    out = {
        'Sigma': sigma,
        'QVar': qvar,
        'QStd': None if qvar is None else np.sqrt(np.maximum(0.0, qvar)),
        'sigma2': arr(m.get('sigma2'), ravel=True),
        # The variance the drift was closed at. getCdfRespT re-integrates the
        # passage-time ODE on that same drift, so dropping it here answers with
        # the first-order min() and reports the empty-station service law.
        'sigma2Drift': arr(m.get('sigma2Drift'), ravel=True),
        'outerIters': int(m.get('outerIters') or 0),
    }

    blocks = [np.asarray(b, dtype=np.int64).reshape(-1)
              for b in (m.get('stationBlock') or [])]
    out['stationBlock'] = blocks
    flat = [np.asarray(b, dtype=np.int64).reshape(-1)
            for b in (m.get('classBlock') or [])]
    M = len(blocks)
    K = (len(flat) // M) if M else 0
    out['classBlock'] = [[flat[i * K + k] for k in range(K)] for i in range(M)]

    caches = m.get('cache')
    if caches:
        out['cache'] = [{
            'node': int(c.get('node', -1)),
            'pi0': arr(c.get('pi0'), ravel=True),
            'Sigma': arr(c.get('Sigma')),
            'pi0Var': arr(c.get('pi0Var'), ravel=True),
            'missProbVar': arr(c.get('missProbVar'), ravel=True),
        } for c in caches]
    return out


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
    results = solve_via_jar(solver, analysis='avg', timeout=timeout)
    return avg_matrices_from_results(results, sn, source='jline.jar')


def avg_matrices_from_results(results, sn, source='jline.jar'):
    """
    Reduce a CLI avg payload to station x class metric matrices.

    TRANSPORT-AGNOSTIC ON PURPOSE. `line-cli -a avg -o json` emits the same
    column-oriented object the JAR does, key for key, so both delegating
    backends (lang='java' and lang='cpp') share this reduction instead of
    carrying a parser each. Two parsers for one wire format is two places for a
    column rename to be missed, and the second one would be found by a wrong
    number rather than by an error.

    Args:
        results: the parsed CLI JSON, {'avg': {...}} or the table itself.
        sn:      the network struct whose indexing the matrices must match.
        source:  binary name, used only in the error message.
    """
    import numpy as np

    table = results.get('avg', results)
    if not isinstance(table, dict) or 'QLen' not in table:
        raise RuntimeError(
            "%s did not return a structured AvgTable (got keys: %s)."
            % (source, list(table.keys()) if isinstance(table, dict) else type(table))
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
    # Iteration count and convergence flag of the delegated solve, emitted by the
    # CLI as siblings of 'avg'. They are what lets the delegating solver raise the
    # same non-convergence warning its native path raises; without them the same
    # model warns under lang='python' and is silent under lang='java', and the
    # silence reads as convergence. `converged` is absent when the handler reports
    # no flag (then the count is the sound signal), which is not the same as False.
    mats['_lG'] = results.get('lG') if isinstance(results, dict) else None
    mats['_iter'] = results.get('iter') if isinstance(results, dict) else None
    mats['_converged'] = results.get('converged') if isinstance(results, dict) else None
    return mats


def populate_java_result(solver, timeout=None):
    """
    Populate a native solver's result container (and struct) from the JAR, so
    every downstream getter returns JAR-derived values. Called from each
    solver's runAnalyzer when options.lang == 'java'.
    """
    import time

    import numpy as np

    from .base import print_delegated_banner

    _t0 = time.time()
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
    # The scalar metadata is written through __setitem__, which sets the dict key
    # AND the attribute: a native result is a dict-like that answers res['method']
    # as readily as res.method, and a container that carries only the attribute
    # raises KeyError on the other spelling (the map-env fallback reports itself
    # through result['method'], and read that way it looked absent under
    # lang='java' while the JAR had emitted it).
    #
    # Method actually executed by the JAR (SSA reports 'nrm' vs 'serial'); native
    # solvers expose result.method, and tests assert it to catch a silent serial
    # fallback. None for solvers/CLIs that do not emit it.
    container['method'] = mats.get('_method')
    # Log normalizing constant of a normalizing-constant solve; None for the
    # solvers that compute none, which is what the native result reports too.
    _lg = mats.get('_lG')
    container['lG'] = float(_lg) if _lg is not None else None
    # Iteration count and convergence flag of the delegated solve; native results
    # expose result.iter, and the delegating runAnalyzer reads both to decide
    # whether to raise the non-convergence warning. None when the CLI omits them.
    _it = mats.get('_iter')
    container['iter'] = int(_it) if _it is not None else None
    # SolverNC's native result names the same count `it`; both spellings are read
    # by callers, so a container that carries one of them answers half of them.
    container['it'] = container.__dict__['iter']
    _cv = mats.get('_converged')
    container['converged'] = bool(_cv) if _cv is not None else None
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
    print_delegated_banner(solver, 'java', mats.get('_method'), time.time() - _t0)
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
    model = _get_model(solver)
    model_nodes = model.get_nodes() if hasattr(model, 'get_nodes') else []

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
        # `Cache.get_hit_ratio()` reads the NODE, not the struct, so a solve that
        # seeds only `nodeparam` leaves every hit/miss ratio the model reports
        # None -- cache_compare_replc prints nan for all nine of its ratios.
        node = model_nodes[ind] if ind < len(model_nodes) else None
        if node is not None and hasattr(node, 'set_result_hit_prob'):
            node.set_result_hit_prob(ahp)
            node.set_result_miss_prob(amp)

    _refresh_cache_visits(sn)


def _refresh_cache_visits(sn):
    """Re-route each Cache read class at its actual hit/miss split and refresh
    the visits, as the native cacheqn analyzer does after it converges.

    `sn.visits` carries the NOMINAL split the routing matrix was linked with --
    hit and miss equally likely -- until a solve replaces it. Native SolverMVA
    rewrites `sn.rtnodes` at the converged split, recomputes `sn.rt` by
    stochastic complement and calls `sn_refresh_visits`; a delegated solve never
    runs that pass, so every getter that multiplies by a visit ratio answers at
    the nominal split even though the delegated engine solved the right one.
    Measured on cache_replc_routing: ResidT at (Delay1, HitClass) read 0.025
    (visit 0.5*0.5) where MATLAB, native Python and the JAR's own node table
    read 0.0197 (visit 0.394*0.5).

    THE DEFECT IS THIS BRIDGE'S, NOT THE JAR'S, and that is what localises it:
    drive the JAR CLI on the SAME exported model.json and it prints 0.0197
    itself. So the delegated engine had the right answer and the host discarded
    it. `lang='cpp'` reaches the same pass, so the body lives in api/sn.
    """
    from ..api.sn.transforms import sn_refresh_cacheqn_visits
    sn_refresh_cacheqn_visits(sn)


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
                 onebased=False, raw_station=False, timeout=None, state=None):
    """
    Delegate a probability query to jline.jar and return it in the native
    contract. `kind` selects the shape: 'scalar' -> float(probability),
    'logtuple' -> (logProbability, probability), 'vector' -> probabilityVector.
    `ist`/`jclass` are the native station/class indices (converted using the
    solver's 0/1-based convention via `onebased`). `raw_station=True` passes the
    0-based station index directly (used by prob-marg, whose JAR path treats -n
    as a station index rather than a stateful-node index). `state` names the
    per-class job counts being asked about; it must be given whenever the model's
    own initial state matters, because the model interchange carries no state and
    the JAR would otherwise answer at its default initialization. It is a flat
    list of per-class counts for a node-level query, and a list of one such row
    per station for the system-level ones.
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
                            state=state, timeout=timeout)
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
_ITEM_COLS = ['Node', 'Item', 'List', 'ListCap', 'Prob',
              'DelayedHitQLen', 'DelayedHitQLenFull']
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
    (QN, UN, RN, TN, AN, WN) as (nidx,) arrays indexed by LQN absolute index,
    matching the native SolverLN.get_ensemble_avg contract.

    This solves the whole LayeredNetwork in ONE JAR subprocess. Without it, a
    lang="java" ensemble solve falls through to the native LN fixed-point loop
    whose per-layer sub-solves each dispatch to their own `java -jar` process
    (hundreds of JVM startups per LN solve).
    """
    df = ln_avg_table_via_jar(solver, timeout=timeout)
    lqn = solver.lqn
    nidx = lqn.nidx
    QN, UN, RN, TN, AN, WN = (np.full(nidx, np.nan) for _ in range(6))
    # The JAR's Node column carries the same hashnames the native layered table
    # builds, so index by name rather than relying on row order. Element indices
    # run 0..nidx-1, the space the native ensemble arrays use.
    name2idx = {solver._get_hashname(idx): idx for idx in range(nidx)}
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
    """Delegate getAvgCacheTable() to jline.jar.

    THE DELEGATED TABLE MUST HAVE THE NATIVE TABLE'S SHAPE, not merely its
    numbers. `build_cache_avg_table` emits a `ListCost` column and INTEGER
    `List` and `Items`; the JAR's `NetworkAvgCacheTable` carries neither -- it
    has no ListCost key and reports both indices as doubles. Two things follow,
    and the second is the one that bit:

    - a caller reading `get_avg_cache_table()['ListCost']` gets a column under
      `lang='python'` and a KeyError under `lang='java'`;
    - `to_string` then prints the WHOLE FRAME differently. pandas picks its float
      format from the frame, not from each column alone, so dropping one all-NaN
      column or widening two integer columns moves every metric from five
      SIGNIFICANT digits to five DECIMAL places: ArvR 0.024273 became 0.02427.
      The values are bit-identical -- 0.02427335858096497 in both -- but the
      parity harness reads the PRINTED table, so retrieval_ps [P2J] failed on a
      cell where the two engines agree to sixteen digits.

    So the columns are added rather than the golden loosened.
    """
    df = _table_via_jar(solver, 'cache', _CACHE_COLS, timeout=timeout)
    import numpy as np
    import pandas as pd
    for c in ('List', 'Items'):
        if c in df.columns:
            v = pd.to_numeric(df[c], errors='coerce')
            # Only where every entry is a whole number: a fractional index is
            # not an index, and truncating one would hide that.
            if v.notna().all() and np.all(np.equal(np.mod(v.to_numpy(), 1), 0)):
                df[c] = v.astype('int64')
    if 'ListCost' not in df.columns:
        # The JAR reports no per-list storage cost; the native table carries the
        # column as NaN when the model declares no item sizes, which is the same
        # statement.
        df['ListCost'] = np.nan
    # AND THE SAME TYPE, which is what actually decides the printed digits:
    # `IndexedTable.to_string` renders 5 SIGNIFICANT digits where pandas renders
    # `display.precision` DECIMAL places, and its own docstring names this cell.
    # A bare DataFrame here is a table that prints differently from the native
    # one and cannot be filtered by Station/JobClass objects either.
    from ..indexed_table import IndexedTable
    df = IndexedTable(df)
    # The native build_cache_avg_table() solves through getAvgNode(), which seeds
    # every Cache node's result fields as a side effect; delegating must reproduce
    # that side effect or a caller reading node.get_delayed_hit_qlen() after
    # getAvgCacheTable() sees None under lang='java' and numbers under 'python'.
    _inject_cache_delayed_hit_qlen_via_jar(solver, timeout=timeout)
    return df


def item_table_via_jar(solver, timeout=None):
    """Delegate getAvgItemTable() to jline.jar.

    Wrapped in `IndexedTable` for the reason `cache_table_via_jar` is: the native
    builder returns one, and the wrapper is what prints 5 significant digits.
    """
    df = _table_via_jar(solver, 'item', _ITEM_COLS, timeout=timeout)
    _seed_cache_delayed_hit_qlen(solver, df)
    from ..indexed_table import IndexedTable
    return IndexedTable(df)


def _seed_cache_delayed_hit_qlen(solver, df):
    """
    Copy the per-item delayed-hit queue length out of a JAR item table onto the
    native Cache nodes and their nodeparam, matching what the native
    SolverCTMC._compute_cache_delayed_hit_qlen() sets. The item table carries one
    row per (item, list), so the per-item value repeats across lists and the first
    row for an item is taken.
    """
    import numpy as np

    if df is None or not hasattr(df, 'iterrows'):
        return
    if 'DelayedHitQLen' not in df.columns or 'DelayedHitQLenFull' not in df.columns:
        return
    model = getattr(solver, 'model', None)
    sn = getattr(solver, '_sn', None)
    if sn is None and model is not None:
        sn = model.get_struct()
    if sn is None or model is None:
        return
    nodes = model.get_nodes() if hasattr(model, 'get_nodes') else getattr(model, '_nodes', None)
    if not nodes:
        return
    nodenames = [str(x) for x in list(sn.nodenames)]
    name_to_node = {nm: i for i, nm in enumerate(nodenames)}

    # node index -> {item index (0-based): (d1, dfull)}
    per_node = {}
    for _, row in df.iterrows():
        ind = name_to_node.get(str(row.get('Node')))
        if ind is None:
            continue
        try:
            i = int(row.get('Item')) - 1
            d1 = float(row.get('DelayedHitQLen'))
            dfull = float(row.get('DelayedHitQLenFull'))
        except (TypeError, ValueError):
            continue
        if i < 0 or not np.isfinite(d1) or not np.isfinite(dfull):
            continue
        per_node.setdefault(ind, {}).setdefault(i, (d1, dfull))

    nodeparam = getattr(sn, 'nodeparam', None)
    for ind, items in per_node.items():
        n = max(items) + 1
        d1 = np.zeros(n)
        dfull = np.zeros(n)
        for i, (a, b) in items.items():
            d1[i] = a
            dfull[i] = b
        if ind < len(nodes) and hasattr(nodes[ind], 'set_result_delayed_hit_qlen'):
            nodes[ind].set_result_delayed_hit_qlen(d1, dfull)
        cp = nodeparam.get(ind) if isinstance(nodeparam, dict) else None
        if cp is not None:
            cp.delayedhitqlen = d1
            cp.delayedhitqlenfull = dfull
            cp.delayedhitprobitem = dfull - d1


def _inject_cache_delayed_hit_qlen_via_jar(solver, timeout=None):
    """
    Seed the native Cache nodes' delayed-hit queue length from the JAR item table.
    Only a retrieval-system cache has one, so the extra JAR call is gated on
    retrieval_system_capacity > 0 and never fires for an ordinary cache model.
    """
    from ..api.sn.network_struct import NodeType

    model = getattr(solver, 'model', None)
    sn = getattr(solver, '_sn', None)
    if sn is None and model is not None:
        sn = model.get_struct()
    nodetype = getattr(sn, 'nodetype', None)
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodetype is None or not isinstance(nodeparam, dict):
        return
    has_retrieval = False
    for ind in range(int(sn.nnodes)):
        if ind >= len(nodetype) or nodetype[ind] != NodeType.CACHE:
            continue
        cp = nodeparam.get(ind)
        if cp is not None and int(getattr(cp, 'retrieval_system_capacity', 0) or 0) > 0:
            has_retrieval = True
            break
    if not has_retrieval:
        return
    _seed_cache_delayed_hit_qlen(solver, _table_via_jar(solver, 'item', _ITEM_COLS, timeout=timeout))


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
