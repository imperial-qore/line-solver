"""
Delegation to the C++ multiprecision solver (`line-cli`), for `options.lang='cpp'`.

Transport is subprocess + JSON, the SAME transport `jar_dispatch` uses for
lang='java':

    native model --(linemodel_io.save_model)--> model.json
        --> line-cli -f model.json -i json -s <solver> -a avg -o json
        --> {"avg": {...}} --(avg_matrices_from_results)--> result container

WHY A SUBPROCESS AND NOT AN EXTENSION MODULE. The C++ port's whole host boundary
is one function, `line::reg::api_invoke`, plus the CLI's model-solving path; a
nanobind module would call exactly the same code behind exactly the same JSON
conversion policy (`cpp/include/line/reg/api_json.h`). Starting with the
subprocess costs one process per solve and buys: no build-time dependency on the
Python headers, no ABI coupling, no wheel to publish per interpreter, and an
identical code path to the already-proven lang='java' route. An in-process
binding is a drop-in replacement for `_run_line_cli` when the numbers are worth
the packaging -- nothing above that function knows how the JSON arrived.

WHAT FALLS BACK AND WHAT DOES NOT. `lang='cpp'` is an assertion about what
produced the numbers, so this module degrades to native Python for exactly one
reason: the binary is absent or unrunnable on this platform (no wheel for this
arch, a stale build, a missing shared library). Every other failure -- a
construct the C++ analyzer refuses, a non-zero exit, unparseable output -- is
raised. Silently answering with a different engine than the one the caller named
would hide precisely the class of defect that cross-solver comparison exists to
find.
"""

import json
import os
import shutil
import subprocess
import tempfile

import numpy as np

from .jar_dispatch import (_JavaAvgResult, _extract_json_object, _forward_if_overridden,
                           _get_model, _resolve_struct, avg_matrices_from_results,
                           _LN_COLS, _NATIVE_DEFAULT_ITER_MAX, _NATIVE_DEFAULT_ITER_TOL,
                           _NATIVE_DEFAULT_TOL, station_col_ranges_from_widths)


class LineCliNotAvailable(RuntimeError):
    """The `line-cli` binary could not be located or executed on this platform."""


# The C++ model-solving path's own method name set (`solve_model_dispatch`), which is
# NOT the JAR's: it has `ba`, which the JAR CLI's validateSolver lacks.
# Resolved here rather than through `jar_dispatch._solver_token` so an
# unsupported solver is refused with the reason that applies to THIS backend.
#
# `jmt` USED TO BE MISSING FROM THIS TABLE, with a comment saying the port did
# not carry it. It does: `line-cli -s jmt` drives JSIM through the same JMT jar
# this wrapper does, and its own `-s jmt -a prob` arm reads the simulation logs.
# The table is the authority for what lang='cpp' will accept, so a token absent
# from it makes the CLI's arm unreachable no matter what the CLI grows -- which
# is the decay this whole file's refusals are prone to. Checked against the
# CLI's own list, which is: mva, nc, ctmc, mam, ba, ssa, fluid, ldes, jmt, uq,
# env and qns. `ldes` is deliberately still absent: the Python LDES wrapper is
# ALREADY a JSON subprocess client of that same engine, so lang='cpp' there
# would name the backend it is on. `lqns` the CLI genuinely does not serve.
_CPP_SOLVER_TOKENS = {
    'MVA': 'mva',
    'NC': 'nc',
    'CTMC': 'ctmc',
    'MAM': 'mam',
    'FLUID': 'fluid',
    'FLD': 'fluid',
    'SSA': 'ssa',
    'BA': 'ba',
    'AG': 'ag',
    'AUTO': 'auto',
    'JMT': 'jmt',
}


def _cpp_solver_name(solver):
    """
    The solver's short name, with any `Solver` prefix stripped.

    getName() is NOT uniform across the native solvers: SolverMVA answers 'MVA'
    and SolverBA answers 'SolverBA'. Normalising here is what keeps the method name map
    keyed by the solver and not by that inconsistency.
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
    return name[6:] if name.upper().startswith('SOLVER') else name


def _cpp_token(solver):
    """Resolve the `line-cli -s` method name from a native solver instance."""
    name = _cpp_solver_name(solver)
    token = _CPP_SOLVER_TOKENS.get(str(name).upper())
    if token is None:
        raise RuntimeError(
            "lang='cpp' is not available for solver '%s': the C++ model-solving path serves "
            "%s here. LQNS the CLI does not carry; the LDES wrapper is already a subprocess "
            "client of the same engine, so lang='cpp' would name the backend it is on; and a "
            "LayeredNetwork goes through the layered path instead. Use lang='python' or "
            "lang='java'." % (name, ', '.join(sorted(set(_CPP_SOLVER_TOKENS.values())))))
    return token


def cpp_unsupported(solver, method, reason):
    """
    Raise for a getter this bridge cannot serve, naming the getter and the reason.

    THE REASON IS THE POINT. It used to be one reason for every refusal --
    `line-cli` printed `-o json` for `-a avg` only -- and that is no longer true:
    every analysis the CLI implements now emits a payload, so what is left
    unreachable is unreachable for a SPECIFIC, per-getter cause, and collapsing
    those into "avg only" would tell a caller to stop asking for something the
    port may in fact be one arm away from answering. The three causes that
    actually arise:

    * the C++ port has no such function (the CTMC per-state `getProb` family,
      `getTranAvg`'s transient means);
    * the CLI arm exists but not for this solver (`-a prob` under MAM or SSA:
      only mva and nc have it);
    * the quantity would have to be computed here, from another engine's numbers,
      to be returned at all -- which is the one thing `lang='cpp'` may not do.
    """
    name = _cpp_solver_name(solver)
    raise RuntimeError(
        "lang='cpp' cannot serve %s%s: %s. Use lang='python' or lang='java' for it."
        % (('Solver%s.' % name) if name else '', method, reason))


# Kept as the old name so a caller that still imports it gets the same refusal
# with the same message shape; new call sites name their own reason.
def cpp_only_avg(solver, method, analysis=None):
    """Refuse a getter whose analysis the C++ CLI serves for other solvers only."""
    return cpp_unsupported(
        solver, method if analysis is None else '%s (-a %s)' % (method, analysis),
        "the C++ solver does not implement this analysis for this solver")


def _payload(results, analysis):
    """
    The payload of the analysis that was ASKED FOR, or a refusal naming what came.

    Every `line-cli` answer is keyed by its own `-a`, so this is also the check
    that the process answered the question the caller put: a payload under
    another key is a wrong answer, and reading it positionally would turn that
    into wrong numbers instead of an error.
    """
    if not isinstance(results, dict) or analysis not in results:
        raise RuntimeError(
            "line-cli did not answer -a %s (keys: %s)."
            % (analysis, list(results.keys()) if isinstance(results, dict) else type(results)))
    p = results[analysis]
    if not isinstance(p, dict):
        raise RuntimeError("line-cli's -a %s payload is not an object (got %s)."
                           % (analysis, type(p)))
    return p


def _assert_default_state(solver, what):
    """
    Refuse a state-dependent analysis whose initial state the wire cannot carry.

    THE SINGLE STATE NOW TRAVELS. `linemodel_save` emits the (stateSpace,
    statePrior) pair for every stateful node that carries a state, including the
    one-row space with prior [1] that `setState` and `initFromMarginal` leave,
    and the C++ `network_reader` stores it; `analyzer_detail::default_init_state`
    takes that row in preference to the default marking. So `-a prob`, the
    transient arms and a sampled trajectory all start where the caller put the
    jobs, and this used to refuse them all.

    WHAT IS STILL REFUSED IS A GENUINE DISTRIBUTION over several rows, AND ONLY
    ON THE ARMS THAT SEED ONE STATE. `analyzer_detail::init_state_index`, which
    `-a sample` walks from and `-a prob` reports about, takes the DEFAULT
    marking, so a prior over k > 1 states would be answered as one state under
    the name of a question about a mixture -- the silent substitution this guard
    exists to prevent.

    IT DOES NOT APPLY TO THE TRANSIENT ARMS. `solver_ctmc_transient_analyzer`
    seeds `analyzer_detail::init_state_distribution`, the product of the
    declared per-node priors over the enumerated space, and integrates the
    forward equation from it; every quantity it reports is linear in pi(t) and
    pi(t) in pi(0), so that single integration IS the reference's weighted sum
    over the prior's support. `-a tran`, `-a tranprob` and `-a tranreward`
    therefore take a mixture and must not be gated here.
    """
    model = _get_model(solver)
    named = []
    for node in list(model.get_nodes()):
        prior = getattr(node, '_state_prior', None)
        space = getattr(node, '_state_space', None)
        rows = 0
        if space is not None:
            space = np.atleast_2d(np.asarray(space))
            rows = space.shape[0] if space.size else 0
        if (prior is not None and np.asarray(prior).size > 1) or rows > 1:
            named.append(str(getattr(node, 'name', '?')))
    if named:
        raise RuntimeError(
            "lang='cpp' cannot serve %s on a model whose initial state is a distribution over "
            "several states (%s): the C++ reader takes a one-row state space as the initial "
            "state and rebuilds the default marking for any other, so it would answer for one "
            "state under a question about a mixture. The reference weights the analysis over "
            "the prior's support; use lang='python' or lang='java' for it."
            % (what, ', '.join(named)))


def analysis_via_cpp(solver, analysis, timeout=None, flags=None):
    """Run one non-average `line-cli` analysis and return its payload."""
    return _payload(solve_via_cpp(solver, analysis=analysis, timeout=timeout, flags=flags),
                    analysis)


def find_line_cli():
    """
    Locate the `line-cli` binary.

    Search order, most explicit first: the LINE_CLI_BINARY environment variable,
    the checkout's common/ directory, the in-tree cpp/build directories, then
    PATH. This is the order of MATLAB's CPPLINE.findLineCli, and common/ comes
    before cpp/build because `cpp/make.sh` installs there: a stale cpp/build
    left over from an earlier manual build must not shadow the fresh artifact.
    Raises LineCliNotAvailable rather than returning None, so the one recoverable
    failure has one type and the caller's `except` cannot accidentally swallow a
    solver refusal alongside it.
    """
    env = os.environ.get('LINE_CLI_BINARY')
    if env:
        if os.path.isfile(env) and os.access(env, os.X_OK):
            return env
        raise LineCliNotAvailable(
            "LINE_CLI_BINARY is set to '%s', which is not an executable file." % env)

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(here, '..', '..', '..'))
    for rel in ('common/line-cli', 'cpp/build/line-cli', 'cpp/build-gmp/line-cli'):
        cand = os.path.join(root, rel)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand

    found = shutil.which('line-cli')
    if found:
        return found

    raise LineCliNotAvailable(
        "the C++ solver binary 'line-cli' was not found. Set LINE_CLI_BINARY to its "
        "path, put it on PATH, or build it with "
        "`cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build`.")


def _run_line_cli(binary, cmd, timeout=None):
    """
    Run one `line-cli` invocation and return its parsed JSON object.

    THE SEAM AN IN-PROCESS BINDING WOULD REPLACE. Everything above this call
    works on the parsed object and does not know a process was involved.
    """
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              timeout=timeout)
    except OSError as e:
        # An unrunnable binary is the platform failure lang='cpp' may degrade
        # from; a solver refusal below is not.
        raise LineCliNotAvailable("could not execute '%s': %s" % (binary, e))

    stdout = proc.stdout.decode('utf-8', errors='replace')
    stderr = proc.stderr.decode('utf-8', errors='replace')
    if proc.returncode != 0:
        # A STALE BINARY IS THE FIRST THING TO RULE OUT, and it announces itself
        # with one specific refusal: `-o json` on the model-solving path did not
        # exist before the wrapper needed it, so a binary predating it says so.
        # Reported as its own message because the raw refusal reads like a
        # missing feature rather than an out-of-date build, and this harness hit
        # it the first time it ran against the tree's `common/line-cli`.
        if 'prints only -o readable' in stderr:
            raise RuntimeError(
                "'%s' is a STALE line-cli: it predates the -o json output path this "
                "wrapper needs. Rebuild it (`cmake --build <builddir> --target line-cli`) "
                "and refresh the copy this wrapper resolves to, or point LINE_CLI_BINARY "
                "at the fresh one." % binary)
        # exit 2 is line-cli's `line::Error` channel: a refusal by name, with the
        # construct in the message. It is relayed verbatim rather than
        # translated, because the message IS the port's contract with a caller.
        raise RuntimeError(
            "line-cli exited with code %d.\nstderr:\n%s" % (proc.returncode, stderr.strip()))
    if stderr.strip():
        import warnings
        warnings.warn("line-cli: %s" % stderr.strip())
    try:
        return _extract_json_object(stdout)
    except Exception as e:
        # THE SECOND STALE-BINARY SIGNATURE, and it is silent where the one above
        # is loud: `-s ssa` and `-s fluid` used to ACCEPT `-o json` and print the
        # readable table anyway (their arms hand-rolled a printer that never
        # consulted the output mode). A binary predating that fix exits 0 and
        # emits a table with no brace in it, so the failure surfaces here as an
        # unparseable answer rather than as the out-of-date build it is.
        if 'Station' in stdout and 'JobClass' in stdout and '{' not in stdout:
            raise RuntimeError(
                "'%s' answered -o json with the READABLE table, which means it predates the "
                "shared JSON emitter for the -s ssa and -s fluid arms. Rebuild it "
                "(`cpp/make.sh`, or `cmake --build <builddir> --target line-cli`) and point "
                "LINE_CLI_BINARY at the fresh one." % binary)
        raise RuntimeError(
            "could not parse line-cli JSON output: %s\nstdout:\n%s\nstderr:\n%s"
            % (e, stdout.strip(), stderr.strip()))


def _qrf_params_json(qp):
    """options.config.qrf_params as the --qrf-params document.

    ZM is not sent: the CLI derives it from ZZ, as every codebase now does.
    The tables go out as arrays of rows, so a one-row table stays a table.
    """
    get = qp.get if isinstance(qp, dict) else lambda k, d=None: getattr(qp, k, d)
    missing = [k for k in ('f', 'MR', 'BB', 'MM', 'MM1', 'ZZ') if get(k) is None]
    if missing:
        raise ValueError(
            'options.config.qrf_params must contain field(s) %s for the QRF blocking bounds'
            % ', '.join(missing))
    doc = {
        'f': int(get('f')),
        'MR': int(get('MR')),
        'BB': np.atleast_2d(np.asarray(get('BB'), dtype=int)).tolist(),
        'MM': np.atleast_2d(np.asarray(get('MM'), dtype=int)).tolist(),
        'MM1': np.atleast_2d(np.asarray(get('MM1'), dtype=int)).tolist(),
        'ZZ': np.asarray(get('ZZ'), dtype=int).ravel().tolist(),
    }
    F = get('F')
    if F is not None and len(np.asarray(F).ravel()) > 0:
        doc['F'] = np.asarray(F, dtype=int).ravel().tolist()
    return json.dumps(doc)


def solve_via_cpp(solver, analysis='avg', timeout=None, arith=None, flags=None):
    """
    Serialize the solver's model, run `line-cli` on it, return the parsed results.

    Only the knobs the C++ CLI actually honours are forwarded, and only when the
    caller overrode the native default -- an untouched option lets the port apply
    its own SolverOptions default rather than having the wrapper's default
    imposed on it. That is the same rule `solve_via_jar` follows.

    `flags` are the analysis's OWN argv entries, assembled by the caller that
    knows the getter: `--tspan` for a transient query, `--node` for a per-node
    one, `--notation` for the ODE export. They are passed through rather than
    derived from options because they are arguments of the CALL and not settings
    of the solver -- `getTranProbSys(t)` names its horizon in the call, and there
    is no `options.tspan` to read it from.

    KNOBS ARE ALSO GATED BY TOKEN, which the JAR route does not have to do.
    line-cli REFUSES an option the chosen solver does not have (`--tol` on `-s
    ba`: "a bound is a closed form, with nothing to converge"; `--samples` on
    anything but `ssa`; `--cutoff` outside `ctmc`). That discipline is what makes
    its answers attributable, but it means a solver whose native default merely
    differs from the generic one would be refused for an option the caller never
    set. So each flag is sent only to the tokens that have it.
    """
    from ..io.linemodel_io import save_model

    binary = find_line_cli()
    token = _cpp_token(solver)
    model = _get_model(solver)

    method = None
    tol = iter_tol = iter_max = samples = seed = cutoff = None
    opts = getattr(solver, 'options', None)
    if opts is not None:
        m = getattr(opts, 'method', None)
        if m is not None and str(m) and str(m) != 'default':
            method = str(m)
        # NOTHING IS RESOLVED HERE FOR THE FLUID SOLVER ANY MORE. It used to be:
        # the two ports resolved 'default' differently, because the reference (and
        # the native solver, and the JAR) prefer the second-order closure
        # 'minnormal' wherever it applies while the C++ resolved to 'matrix', so a
        # closed PS network came back QLen [2, 1] under lang='cpp' against
        # [1.575, 1.425] under lang='python' -- both correct for the method that
        # ran. The C++ port now carries fluid_resolve_default_method and
        # fluid_minnormal_applicable, so `line-cli` resolves 'default' itself and
        # forwarding a name computed on this side would OVERRIDE that resolution
        # with one taken from a different model struct.
        # ba/ssa/ctmc have no convergence loop to control; mam and ag carry tol
        # and iter_max but not iter_tol.
        if token not in ('ba', 'ssa', 'ctmc'):
            tol = _forward_if_overridden(getattr(opts, 'tol', None), _NATIVE_DEFAULT_TOL)
            # SolverAGOptions and SolverMAMOptions spell the sweep budget
            # `max_iter`; every other options class spells it `iter_max`.
            # Reading one name only dropped the other's budget SILENTLY -- a
            # caller who asked for 500 sweeps got the port's own default and was
            # told nothing. Both ports default to 100 there, so the alias
            # changes no existing answer, only whether an override arrives.
            im = _forward_if_overridden(
                getattr(opts, 'iter_max', None) if getattr(opts, 'iter_max', None) is not None
                else getattr(opts, 'max_iter', None),
                _NATIVE_DEFAULT_ITER_MAX)
            if im is not None:
                try:
                    iter_max = int(im)
                except (TypeError, ValueError):
                    iter_max = None
            if token not in ('mam', 'ag'):
                iter_tol = _forward_if_overridden(getattr(opts, 'iter_tol', None),
                                                 _NATIVE_DEFAULT_ITER_TOL)
        # A sample path has a length and a stream; nothing else does. events
        # overrides samples for the DES solvers, matching solve_via_jar. `-s ctmc
        # -a sample` takes both too, but they arrive through `flags`: its length
        # is an argument of sampleSys(numEvents) and not an option, and the CLI
        # refuses --samples on every other CTMC analysis.
        if token == 'ssa':
            s = getattr(opts, 'events', None) or getattr(opts, 'samples', None)
            if s is not None:
                try:
                    samples = int(s)
                except (TypeError, ValueError):
                    samples = None
            sd = getattr(opts, 'seed', None)
            if sd is not None:
                try:
                    sd = int(sd)
                except (TypeError, ValueError):
                    sd = None
                # SolverSSAOptions.seed defaults to 0, which means "unset" here and
                # not "seed zero": line-cli requires a positive integer and refuses
                # 0 outright, so forwarding the default would fail every SSA solve
                # that passed no seed at all.
                seed = sd if sd and sd > 0 else None
        # The cutoff bounds an open population inside a state space. It is a
        # SolverCTMC option and reaches every CTMC analysis, unlike the sample
        # length above.
        if token == 'ctmc':
            c = getattr(opts, 'cutoff', None)
            if c is not None:
                try:
                    c = float(np.max(np.asarray(c, dtype=float)))
                except Exception:
                    c = None
                if c is not None and c > 0:
                    cutoff = c
        if arith is None:
            arith = getattr(opts, 'arith', None)

    tmpdir = tempfile.mkdtemp(prefix='line_cpp_')
    model_path = os.path.join(tmpdir, 'model.json')
    try:
        save_model(model, model_path)
        cmd = [binary, '-f', model_path, '-i', 'json', '-s', token, '-a', analysis, '-o', 'json']
        # --arith is the one flag with no JAR counterpart: it is what the C++
        # port exists for. Absent, the port runs in double, which is the only
        # setting comparable with lang='python' or lang='java'.
        if arith:
            cmd += ['--arith', str(arith)]
        if method is not None:
            cmd += ['--method', method]
        if tol is not None:
            cmd += ['--tol', repr(float(tol))]
        if iter_tol is not None:
            cmd += ['--iter_tol', repr(float(iter_tol))]
        if iter_max is not None and iter_max > 0:
            cmd += ['--iter_max', str(iter_max)]
        # AMVA multiserver rule: same reason as the JAR transport -- a rule that
        # does not reach the solver leaves it answering under 'default' while the
        # caller believes its own choice was applied.
        _cfg = getattr(opts, 'config', None) if opts is not None else None
        _ms = None
        if isinstance(_cfg, dict):
            _ms = _cfg.get('multiserver')
        elif _cfg is not None:
            _ms = getattr(_cfg, 'multiserver', None)
        if _ms is not None and str(_ms) and str(_ms) != 'default':
            cmd += ['--multiserver', str(_ms)]
        # The QRF reduction bounds of SolverBA read their blocking tables and
        # their load-dependent scaling from options.config. Without them the
        # port refuses qrf.bas and qrf.rsrd, so a caller that set them and saw
        # them dropped here would get a refusal it did not cause.
        if token == 'ba':
            _lvl = getattr(opts, 'level', None) if opts is not None else None
            if _lvl is not None and int(_lvl) > 0:
                cmd += ['--level', str(int(_lvl))]
            _qp = _cfg.get('qrf_params') if isinstance(_cfg, dict) else getattr(
                _cfg, 'qrf_params', None) if _cfg is not None else None
            if _qp is not None:
                cmd += ['--qrf-params', _qrf_params_json(_qp)]
            _qa = _cfg.get('qrf_alpha') if isinstance(_cfg, dict) else getattr(
                _cfg, 'qrf_alpha', None) if _cfg is not None else None
            if _qa is not None and len(_qa) > 0:
                cmd += ['--qrf-alpha', json.dumps(np.atleast_2d(np.asarray(
                    _qa, dtype=float)).tolist())]
        # SolverAG's own knob. maxStates is a TRUNCATION LEVEL and so part of
        # the answer, not a budget: a run truncated at the port's 100 when the
        # caller asked for 500 is a different number reported as theirs. The
        # execution backends (config 'exec', 'nworkers', 'endpoints') do not
        # travel -- line-cli exposes no flag for them and the port always sweeps
        # serially -- so a non-serial one is refused rather than dropped.
        if token == 'ag':
            _mx = _cfg.get('maxStates') if isinstance(_cfg, dict) else getattr(
                _cfg, 'maxStates', None) if _cfg is not None else None
            if _mx is not None and int(_mx) > 0 and int(_mx) != 100:
                cmd += ['--max-states', str(int(_mx))]
            _exec = _cfg.get('exec') if isinstance(_cfg, dict) else getattr(
                _cfg, 'exec', None) if _cfg is not None else None
            if _exec is not None and str(_exec) and str(_exec).lower() != 'serial':
                raise RuntimeError(
                    "lang='cpp' cannot serve SolverAG with config['exec']='%s': line-cli exposes "
                    "no execution-backend flag, so the port always sweeps the agents serially. "
                    "The backends produce the SAME iterates, so this refuses only the placement "
                    "of the work, not an answer: use lang='python' or lang='java' when the pool "
                    "or the remote workers are the point." % str(_exec))
        if samples is not None and samples > 0:
            cmd += ['--samples', str(samples)]
        if seed is not None:
            cmd += ['--seed', str(seed)]
        if cutoff is not None:
            cmd += ['--cutoff', repr(float(cutoff))]
        if flags:
            cmd += [str(f) for f in flags]
        return _run_line_cli(binary, cmd, timeout=timeout)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def station_matrices_via_cpp(solver, sn, timeout=None):
    """Run the C++ avg analysis and reduce it to station x class matrices."""
    results = solve_via_cpp(solver, analysis='avg', timeout=timeout)
    mats = avg_matrices_from_results(results, sn, source='line-cli')
    mats['_method'] = results.get('method') if isinstance(results, dict) else None
    mats['_arith'] = results.get('arith') if isinstance(results, dict) else None
    # Per-Cache hit/miss/delayed-hit vectors ride with the avg body; None when the
    # solver computed none, which the consumer reads as "clear the node".
    _avg = results.get('avg') if isinstance(results, dict) else None
    mats['_cache'] = _avg.get('Cache') if isinstance(_avg, dict) else None
    return mats


_ENV_STAGE_TOKENS = {'SolverFluid': 'fluid', 'SolverFLD': 'fluid', 'Fluid': 'fluid',
                     'FLD': 'fluid', 'SolverCTMC': 'ctmc', 'CTMC': 'ctmc'}


def _env_stage_solver(solver):
    """
    The `--stage-solver` method name of this ensemble's stage solvers, or None.

    ONE TOKEN FOR EVERY STAGE, and a mixed ensemble returns None rather than the
    first stage's: the engine runs one stage solver for the whole coupling, so
    sending one stage's name would answer for the others under it. The caller's
    own gate (`_assert_stages_delegable`) is what turns that None into a named
    refusal; returning None here keeps this function free of that decision.
    """
    names = set()
    for s in list(getattr(solver, '_solvers', []) or []):
        if s is None:
            continue
        tok = _ENV_STAGE_TOKENS.get(type(s).__name__)
        if tok is None:
            return None
        names.add(tok)
    if len(names) != 1:
        return None
    return names.pop()


def _env_stage_cutoff(solver):
    """The stage solvers' shared `options.cutoff`, or None when they disagree."""
    cuts = set()
    for s in list(getattr(solver, '_solvers', []) or []):
        c = getattr(getattr(s, 'options', None), 'cutoff', None) if s is not None else None
        if c is None or not np.isfinite(c):
            continue
        cuts.add(float(c))
    return cuts.pop() if len(cuts) == 1 else None


def env_avg_via_cpp(solver, timeout=None, arith=None):
    """
    Delegate a random-environment solve to `line-cli -s env`, returning
    (QN, UN, TN) as (nstations x nclasses) arrays in the STAGE model's indexing.

    RN and AN are not returned: `SolverENV.getEnsembleAvg` defines them as NaN
    because the environment analyzer computes no response time, and the CLI
    emits nan in those columns for the same reason.

    The Environment is serialized whole -- `save_model` emits the envelope the
    C++ `environment_reader` consumes -- rather than one model.json per stage,
    because the stage-transition process is what makes the answer an environment
    answer and it lives outside any single stage model.
    """
    from ..io.linemodel_io import save_model

    binary = find_line_cli()
    env = getattr(solver, 'env_model', None)
    if env is None:
        raise RuntimeError("lang='cpp' needs the SolverENV's Environment; this solver carries none.")
    stage = None
    for m in list(getattr(solver, 'ensemble', []) or []):
        if m is not None:
            stage = m
            break
    if stage is None:
        raise RuntimeError("lang='cpp' needs a stage model to index the returned matrices by; "
                           "this Environment has no non-empty stage.")

    opts = getattr(solver, 'options', None) or {}
    method = opts.get('method') if isinstance(opts, dict) else getattr(opts, 'method', None)
    if arith is None:
        arith = opts.get('arith') if isinstance(opts, dict) else getattr(opts, 'arith', None)

    tmpdir = tempfile.mkdtemp(prefix='line_cpp_env_')
    model_path = os.path.join(tmpdir, 'model.json')
    try:
        save_model(env, model_path)
        cmd = [binary, '-f', model_path, '-i', 'json', '-s', 'env', '-a', 'avg', '-o', 'json']
        if method is not None and str(method) and str(method) != 'default':
            cmd += ['--method', str(method)]
        # The STAGE horizon, where the caller states a finite one. It is read off
        # the STAGE SOLVER, which is where an example states it (see
        # SolverENV.stage_timespan); left unstated the engine integrates to its
        # own default 100 and the exit averages are a different quadrature.
        ts = solver.stage_timespan()
        if ts is not None:
            cmd += ['--tspan', '%r:%r' % (float(ts[0]), float(ts[1]))]
        # WHICH SOLVER RUNS EACH STAGE, which model.json does not carry: the
        # envelope holds the stage NETWORKS and the transition process, and the
        # ensemble's solver choice lives on the SolverENV rather than in the
        # model. The engine defaults to the fluid transient, so an ensemble built
        # on SolverCTMC stages had to be refused outright; naming it here makes
        # the C++ solve the ensemble the caller actually built. A cutoff travels
        # with it, since an open stage's chain needs one.
        st = _env_stage_solver(solver)
        if st is not None:
            cmd += ['--stage-solver', st]
            cut = _env_stage_cutoff(solver)
            if cut is not None:
                cmd += ['--cutoff', repr(cut)]
        if arith:
            cmd += ['--arith', str(arith)]
        results = _run_line_cli(binary, cmd, timeout=timeout)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    sn = stage.getStruct() if hasattr(stage, 'getStruct') else stage.get_struct()
    mats = avg_matrices_from_results(results, sn, source='line-cli')
    return mats['QN'], mats['UN'], mats['TN']


_CACHE_PROB_KEYS = (('actualhitprob', 'HitProb'),
                    ('actualmissprob', 'MissProb'),
                    ('actualdelayedhitprob', 'DelayedHitProb'))


def _inject_cache_hitprob_via_cpp(solver, sn, cacheblock):
    """Seed each Cache node's actualhitprob/actualmissprob from the avg payload.

    The hit/miss class-switch throughput exists only on the Cache node's output
    edge, so the station x class AvgTable cannot carry it; the native
    getAvgNode()/`sn_get_node_tput_from_tput` derives it as arvTput * prob once
    the probabilities sit on the node. Twin of `_inject_cache_hitprob_via_jar`.

    Read from `avg`'s own `Cache` array rather than a second `-a cache` run: the
    engine already computed it during THIS solve, that array is keyed by node
    NAME and carries one vector per class, and the cache TABLE would have to be
    filtered by its List column first -- a multi-list cache emits per-list rows
    after the total one, whose miss column is nan by construction.

    ABSENT CLEARS. A solver whose solution type carries no cache block (ssa,
    fluid) emits no vector, and leaving the node's copy alone would report the
    PREVIOUS engine's split under this one's banner -- the defect this exists to
    end. A key present for one class and nan for another means that class does
    not read the cache, which is what the getters already test for.
    """
    import numpy as np
    from ..api.sn import NodeType, sn_refresh_cacheqn_visits

    nodeparam = getattr(sn, 'nodeparam', None)
    if nodeparam is None:
        return
    nodetype = list(getattr(sn, 'nodetype', []) or [])
    cache_inds = [ind for ind in range(int(sn.nnodes))
                  if ind < len(nodetype) and nodetype[ind] == NodeType.CACHE]
    if not cache_inds:
        return

    nodenames = [str(x) for x in list(sn.nodenames)]
    name_to_node = {nm: i for i, nm in enumerate(nodenames)}
    R = int(sn.nclasses)
    # `Cache.get_hit_ratio()` reads the NODE, not the struct, so a solve that
    # seeds only `nodeparam` leaves every ratio the model reports None.
    model = getattr(solver, 'model', None) or getattr(solver, 'network', None)
    model_nodes = list(model.get_nodes()) if model is not None and hasattr(
        model, 'get_nodes') else []

    seeded = {}
    for entry in (cacheblock or []):
        if not isinstance(entry, dict):
            continue
        ind = name_to_node.get(str(entry.get('name')))
        if ind is None or ind not in cache_inds:
            continue
        seeded[ind] = entry

    for ind in cache_inds:
        cp = (nodeparam.get(ind) if isinstance(nodeparam, dict)
              else (nodeparam[ind] if ind < len(nodeparam) else None))
        if cp is None:
            continue
        entry = seeded.get(ind)
        written = {}
        for attr, key in _CACHE_PROB_KEYS:
            vec = entry.get(key) if entry is not None else None
            if vec is None:
                setattr(cp, attr, None)
                written[attr] = None
                continue
            arr = np.full(R, np.nan)
            vals = np.asarray(vec, dtype=float).flatten()
            arr[:min(R, vals.size)] = vals[:min(R, vals.size)]
            setattr(cp, attr, arr)
            written[attr] = arr
        node = model_nodes[ind] if ind < len(model_nodes) else None
        if node is not None and hasattr(node, 'set_result_hit_prob'):
            node.set_result_hit_prob(written['actualhitprob'])
            node.set_result_miss_prob(written['actualmissprob'])
            if hasattr(node, 'set_result_delayed_hit_prob'):
                node.set_result_delayed_hit_prob(written['actualdelayedhitprob'])

    # The visits are DERIVED from the split, so seeding it is only half the job.
    sn_refresh_cacheqn_visits(sn)


def populate_cpp_result(solver, timeout=None):
    """
    Populate a native solver's result container from `line-cli`, so every
    downstream getter returns C++-derived values. Called from a solver's
    runAnalyzer when options.lang == 'cpp'.

    Mirrors `populate_java_result`, including the cache hit-probability
    injection, which the avg payload's own `Cache` array carries. Without it
    `sn_get_node_tput_from_tput` finds no actualhitprob on the node and falls
    back to the station throughput, and a node left holding a PREVIOUS native
    solve's split reports that engine's hit rate under this one's `lang`.

    The two cache TABLES do not go through the node at all: they delegate to
    `-a cache` / `-a item` via `cache_table_via_cpp` and `item_table_via_cpp`.
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

    mats = station_matrices_via_cpp(solver, sn, timeout=timeout)
    container = _JavaAvgResult(mats['QN'], mats['UN'], mats['RN'], mats['TN'],
                               mats['AN'], mats['WN'], mats['XN'], mats['CN'])
    container.__dict__['method'] = mats.get('_method')
    container.__dict__['arith'] = mats.get('_arith')
    # Same contract as the JAR delegation: the count and the convergence flag ride
    # the shared avg payload, and stay None when the CLI does not emit them.
    _it = mats.get('_iter')
    container.__dict__['iter'] = int(_it) if _it is not None else None
    _cv = mats.get('_converged')
    container.__dict__['converged'] = bool(_cv) if _cv is not None else None
    container.__dict__['_solver'] = solver
    for attr in ('_result', 'result'):
        try:
            setattr(solver, attr, container)
        except Exception:
            pass
    # Cache nodes: the hit/miss class-switch throughput lives only on the Cache
    # node's output edge, so the station-level results cannot carry it.
    _inject_cache_hitprob_via_cpp(solver, sn, mats.get('_cache'))
    print_delegated_banner(solver, 'cpp', mats.get('_method'), time.time() - _t0)
    return container


# --- the analyses that are not the average table ---------------------------
#
# Each function below mirrors ONE @Solver getter, and returns what that getter's
# native implementation returns -- shape for shape, index base for index base --
# so a solver's `lang='cpp'` branch is a delegation and not a translation its
# caller has to know about. The reshaping that does happen here is transport
# undoing: `line-cli` sends a generator as triplets because a dense generator is
# quadratic on the wire, and sends 0-based indices because the consumer is code.
#
# WHAT IS *NOT* HERE IS AS DELIBERATE AS WHAT IS. There is no adapter for the
# CTMC per-state `getProb` family or for `getTranAvg`: the port has no such
# function, and the only way to return a number would be to compute it natively
# and label it C++. Those getters refuse, by name, through `cpp_unsupported`.


def _dense_from_triplets(payload, n, key_from='From', key_to='To', key_rate='Rate'):
    """Rebuild a dense (n x n) matrix from the wire's sparse triplets."""
    M = np.zeros((int(n), int(n)))
    rows = payload.get(key_from) or []
    cols = payload.get(key_to) or []
    vals = payload.get(key_rate) or []
    if not (len(rows) == len(cols) == len(vals)):
        raise RuntimeError(
            "line-cli sent %d row indices, %d column indices and %d rates for one sparse "
            "matrix; the triplets must be of one length." % (len(rows), len(cols), len(vals)))
    for k in range(len(vals)):
        M[int(rows[k]), int(cols[k])] = float(vals[k])
    return M


def generator_via_cpp(solver, timeout=None):
    """
    The CTMC generator, state spaces and stationary law, as `generator_via_jar`
    returns them: `infgen`, `space`, `space_aggr`/`spaceAggr`, `nodeSpace`,
    `eventFilt`, `pi`.

    TWO INVOCATIONS, ONE PER GETTER. `-a gen` is getInfGen and `-a states` is
    getStateSpace; the CLI keeps them apart because they are different getters,
    so the bridge asks for both and CHECKS THAT THEY AGREE on the state count
    before pairing Q with pi. An unchecked pairing is the one way this can go
    wrong silently: Q from one enumeration indexed by another's states is not a
    chain, and every consumer that reads a row of the space would be reading the
    wrong state rather than getting an error.
    """
    gen = analysis_via_cpp(solver, 'gen', timeout=timeout)
    st = analysis_via_cpp(solver, 'states', timeout=timeout)
    n = int(gen.get('size', 0))
    space = np.atleast_2d(np.asarray(st.get('space') or [], dtype=float))
    if space.size == 0:
        space = np.zeros((0, 0))
    if space.shape[0] != n:
        raise RuntimeError(
            "line-cli returned a %d-state generator and a %d-row state space for the same "
            "model; they cannot be paired." % (n, space.shape[0]))
    space_aggr = np.atleast_2d(np.asarray(st.get('spaceAggr') or [], dtype=float))
    pi = np.asarray(st.get('pi') or [], dtype=float).reshape(-1)
    filt = []
    for e in gen.get('sync') or []:
        filt.append(_dense_from_triplets(e, n))
    return {
        'infgen': _dense_from_triplets(gen, n),
        'space': space,
        'space_aggr': space_aggr,
        'spaceAggr': space_aggr,
        # The C++ keeps a state as its per-node blocks and reports their widths;
        # it does not send the split matrices, and slicing `space` by the widths
        # here would be this module deciding what a node's block is. The native
        # getter slices it from `station_col_ranges` instead, so None is the
        # honest answer and is what the JAR route also returns.
        'nodeSpace': None,
        'eventFilt': filt,
        'pi': pi,
        # `NodeWidths` is where one stateful node's block of a state row ends and
        # the next begins; the native getters cut the row with it.
        'station_col_ranges': station_col_ranges_from_widths(solver, st.get('NodeWidths')),
    }


def prob_aggr_via_cpp(solver, timeout=None):
    """
    `getProbAggr` per station and `getProbSysAggr`, from `-a prob`.

    Available under mva, nc and ctmc, and the three DISAGREE BY CONSTRUCTION:
    SolverMVA fits a binomial to its own means (Schmidt 1997), SolverNC takes a
    ratio of normalizing constants (the product-form model's own probability),
    and SolverCTMC reports the stationary law of the chain itself, exact for any
    model the chain represents. The bridge does not pick between them -- it asks
    the solver the caller named.

    THE DECLARED STATE TRAVELS WITH THE MODEL. `save_model` emits the one-row
    (stateSpace, statePrior) pair that `setState` and `initFromMarginal` leave on
    a node, so all three arms answer about the state the caller set rather than
    rebuilding the default marking. It used not to: measured on the 3-state
    Delay->PS chain, lang='python' answered 0.4 / 0.2 / 0.4 for Q1 at [1,1] /
    [2,0] / [0,2] while lang='cpp' answered 0.4 for all three, because the state
    stopped at the writer. What `_assert_default_state` still refuses is a prior
    over SEVERAL rows, which is a mixture rather than a state.
    """
    _assert_default_state(solver, 'getProbAggr/getProbSys/getProbSysAggr')
    p = analysis_via_cpp(solver, 'prob', timeout=timeout)
    return {
        'stations': [str(x) for x in (p.get('Station') or [])],
        'prob': np.asarray(p.get('Prob') or [], dtype=float),
        'probAggr': np.asarray(p.get('ProbAggr') or [], dtype=float),
        'probSys': float(p['ProbSys']) if p.get('ProbSys') is not None else float('nan'),
        # `-s mva` reports no joint: the binomial fit is a per-station law and
        # the arm sends JSON null rather than a number. nan says "not reported"
        # where float(None) would raise from inside the bridge and read as a
        # transport failure.
        'probSysAggr': (float(p['ProbSysAggr'])
                        if p.get('ProbSysAggr') is not None else float('nan')),
    }


def jmt_prob_aggr_via_cpp(solver, timeout=None):
    """
    `getProbAggr` / `getProbSysAggr` under `-s jmt -a prob`.

    A SIMULATED probability, not an analytic one: the C++ JSIM wrapper runs one
    logged simulation and dwell-weights the fraction of time the target state
    holds, which is what the native path here does with its own sample. The
    `seen` flags come back beside the values because a state the run never
    entered reports 0, and 0 from "never observed" is a different statement from
    0 from "impossible" -- the reference warns rather than letting the two look
    alike.
    """
    _assert_default_state(solver, 'getProbAggr/getProbSysAggr')
    p = analysis_via_cpp(solver, 'prob', timeout=timeout, flags=None)
    return {
        'stations': [str(x) for x in (p.get('Station') or [])],
        'probAggr': np.asarray(p.get('ProbAggr') or [], dtype=float),
        'probSysAggr': (float(p['ProbSysAggr'])
                        if p.get('ProbSysAggr') is not None else float('nan')),
        'stateSeen': [bool(x) for x in (p.get('StateSeen') or [])],
        'sysStateSeen': bool(p.get('SysStateSeen', True)),
    }


def prob_sys_via_cpp(solver, timeout=None):
    """`getProbSys`: the joint probability of the declared state, phases and all."""
    return prob_aggr_via_cpp(solver, timeout=timeout)['probSys']


def prob_marg_via_cpp(solver, ist, jobclass=None, states=None, timeout=None):
    """
    `getProbMarg`'s curve for station `ist` (1-based), from `-a marg`.

    THE TWO SOLVERS ASK DIFFERENT QUESTIONS UNDER ONE NAME. SolverMVA's
    getProbMarg is P(n jobs OF CLASS r) and takes the class and the job counts
    to report; SolverNC's is P(n jobs IN TOTAL at the station), keyed by station
    alone. The payload key differs with them (`marginal` against `curves`), and
    both are read here so neither is quoted under the other's name.
    """
    sn = solver.sn if getattr(solver, 'sn', None) is not None else solver.model.getStruct()
    flags = ['--node', str(int(np.ravel(sn.stationToNode)[int(ist) - 1]) + 1)]
    if jobclass is not None:
        flags += ['--class', str(int(jobclass))]
    if states is not None and len(states):
        flags += ['--marg-states', ','.join(str(int(v)) for v in states)]
    p = analysis_via_cpp(solver, 'marg', timeout=timeout, flags=flags)
    for e in (p.get('marginal') or p.get('curves') or []):
        if int(e['station']) + 1 != int(ist):
            continue
        return (np.asarray(e['P'], dtype=float), np.asarray(e['logP'], dtype=float))
    raise RuntimeError("line-cli's -a marg answer carries no curve for station %d" % int(ist))


def mam_cdf_respt_via_cpp(solver, timeout=None):
    """
    `getCdfRespT` under `-s mam -a cdf`, as the native List[Dict].

    THE WIRE ENTRIES CARRY A CLASS AND NO STATION, and that is not an omission:
    solver_mam_passage_time covers exactly two stations, a Source and one queue,
    so the station is fixed by the model rather than by the class. It is resolved
    the same way `CPPLINE.mamCdfRespT` resolves it -- the station scheduling
    FCFS, HOL, FCFSPRPRIO or PS -- and a model that is not of that shape is an
    error rather than a guessed row.
    """
    sn = solver.sn if getattr(solver, 'sn', None) is not None else solver.model.getStruct()
    wanted = ('FCFS', 'HOL', 'FCFSPRPRIO', 'PS')
    rows = [i for i in range(int(sn.nstations))
            if str(getattr(sn.sched[i], 'name', sn.sched[i])).upper() in wanted]
    if len(rows) != 1:
        raise RuntimeError(
            "the C++ MAM response-time law is defined for a Source and ONE queue, and this "
            "model has %d stations scheduling FCFS, HOL, FCFSPRPRIO or PS, so the curve "
            "line-cli sent could not be placed" % len(rows))
    p = analysis_via_cpp(solver, 'cdf', timeout=timeout)
    out = []
    for e in (p.get('respt') or []):
        out.append({'station': rows[0] + 1,
                    'class': int(e['jobclass']) + 1,
                    't': np.asarray(e['t'], dtype=float),
                    'p': np.asarray(e['F'], dtype=float)})
    return out


def mam_prob_via_cpp(solver, node, timeout=None):
    """
    `getProb` and `getProbMarg` under `-s mam -a prob`, for the 1-based `node`.

    BOTH VIEWS OF ONE QBD LAW come back together: the (level x phase) joint table
    and the per-class marginal. The multi-queue refusal is the REFERENCE'S own --
    QBD is a single-queue method -- and is raised by the CLI rather than
    reproduced here, so the two backends decline the same models with the same
    stated reason.
    """
    flags = ['--node', str(int(node))]
    cutoff = getattr(solver.options, 'cutoff', None)
    if cutoff is not None and np.isfinite(np.max(cutoff)) and np.max(cutoff) > 0:
        flags += ['--cutoff', repr(float(np.max(cutoff)))]
    p = analysis_via_cpp(solver, 'prob', timeout=timeout, flags=flags)
    joint = np.atleast_2d(np.asarray(p.get('joint') or [], dtype=float))
    marginal = {}
    for e in (p.get('marginal') or []):
        marginal[int(e['jobclass'])] = np.asarray(e['P'], dtype=float)
    return {'joint': joint, 'marginal': marginal}


def norm_const_aggr_via_cpp(solver, timeout=None):
    """`getProbNormConstAggr`: log G, from `-a normconst`."""
    p = analysis_via_cpp(solver, 'normconst', timeout=timeout)
    return float(p['logNormConstAggr'])


def statevec_via_cpp(solver, timeout=None):
    """
    The fluid solver's solved ODE state vector (`result.xvec`), from
    `-a statevec`, as the native-named attribute the lazy container caches.

    THE COORDINATE ORDER IS THE CONTRACT, and it is the one every consumer
    assumes: station-major, then class, then phase, which is what the C++ emits
    and what it labels the coordinates with when the layout accounts for the
    whole vector. `labelled: false` says the solver's state is WIDER than that
    layout (an augmented arm), and the passage-time getters index by it, so a
    silently unlabelled vector is refused rather than read positionally.
    """
    p = analysis_via_cpp(solver, 'statevec', timeout=timeout)
    if p.get('labelled') is False:
        raise RuntimeError(
            "line-cli returned an ODE state vector whose coordinates it could not label as "
            "(station, class, phase); the passage-time getters index by that layout, so the "
            "vector cannot be read positionally. Use lang='python' or lang='java' for it.")
    xv = np.asarray(p.get('xvec') or [], dtype=float)
    if xv.ndim == 2 and 1 in xv.shape:
        xv = xv.reshape(-1)
    return {'xvec': xv}


def cdf_respt_via_cpp(solver, timeout=None):
    """
    `getCdfRespT`, from `-a cdf`, in the native List[Dict] contract:
    1-based 'station'/'class' and numpy 't'/'p', as `cdf_respt_via_jar` returns.

    THIS IS THE EXACT LAW, not an exponential of the mean. `-a cdf` builds the
    tagged chain per (station, class) and integrates it, which is what MATLAB's
    `@SolverCTMC/getCdfRespT` does; the native Python getter fits an exponential
    to the mean response time instead, so the two differ in KIND and not only in
    the last digits -- an exponential law and the true one agree on the mean by
    construction and nowhere else. That gap is the native implementation's, and
    delegating leaves it visible rather than reproducing it.

    A pair the chain never visits carries NO curve on the wire (there is no
    arrival event to condition on) and is simply absent from the list, which is
    also how the native and JAR contracts report it.
    """
    p = analysis_via_cpp(solver, 'cdf', timeout=timeout)
    out = []
    for e in p.get('respt') or []:
        out.append({'station': int(e['station']) + 1,
                    'class': int(e['jobclass']) + 1,
                    't': np.asarray(e['t'], dtype=float),
                    'p': np.asarray(e['F'], dtype=float)})
    return out


def cdf_sys_respt_via_cpp(solver, timeout=None):
    """
    `getCdfSysRespT`, from the `sysrespt` block of `-a cdf`, as List[Dict] with
    1-based 'chain' and numpy 't'/'p'.

    ONE LAW PER CHAIN, which is what the native getter and MATLAB's
    `@SolverCTMC/getCdfSysRespT` also return (RD = cell(1, sn.nchains)). The two
    implementations tag the same chain and integrate the same absorption law, so
    this is a delegation and not a substitution.

    A chain with no law on the wire (no tagged arrival to condition on) is
    ABSENT from the list rather than carried as a flat zero curve, matching the
    `respt` block above and the native contract.
    """
    p = analysis_via_cpp(solver, 'cdf', timeout=timeout)
    out = []
    for e in p.get('sysrespt') or []:
        out.append({'chain': int(e['chain']) + 1,
                    't': np.asarray(e['t'], dtype=float),
                    'p': np.asarray(e['F'], dtype=float)})
    return out


class _UnbridgeableReward(Exception):
    """Raised when a declared reward has no serializable form; the caller
    evaluates it locally against the C++ stationary law (see
    `avg_reward_via_cpp`). Carries the reward names, for the message a caller
    that cannot do that must print."""

    def __init__(self, names):
        self.names = list(names)
        Exception.__init__(
            self,
            "the reward(s) %s are defined by a bare callable (or Reward.custom), which has no "
            "serializable form" % ', '.join("'%s'" % n for n in self.names))


class _CppLaw(object):
    """The stationary law `_compute_avg_reward` reads, filled from line-cli.

    It carries exactly the three attributes that getter touches, so a solver
    whose `_result` is swapped for one of these evaluates the caller's reward
    functions against the C++ chain and nothing else changes."""

    def __init__(self, pi, space, space_aggr):
        self.pi = pi
        self.space = space
        self.space_aggr = space_aggr


def _avg_reward_locally(solver, timeout=None):
    """Evaluate the model's reward functions against the C++ stationary law.

    Used when a reward has no serializable form. `generator_via_cpp` supplies
    `pi` and the aggregated state space out of line-cli, and the solver's own
    `_compute_avg_reward` applies the caller's Python functions to them, so the
    numbers are the C++ chain's and only the reward map is local.
    """
    law = generator_via_cpp(solver, timeout=timeout)
    pi = np.asarray(law.get('pi'), dtype=float).reshape(-1)
    space_aggr = law.get('space_aggr')
    if pi.size == 0 or space_aggr is None or space_aggr.size == 0:
        raise RuntimeError(
            "line-cli returned no stationary law for this model, so the reward(s) defined by a "
            "bare callable cannot be evaluated against it.")
    if space_aggr.shape[0] != pi.size:
        raise RuntimeError(
            "line-cli returned a %d-entry stationary law over a %d-row aggregated state space; "
            "a reward cannot be paired with it." % (pi.size, space_aggr.shape[0]))
    saved = getattr(solver, '_result', None)
    try:
        solver._result = _CppLaw(pi, law.get('space'), space_aggr)
        return solver._compute_avg_reward()
    finally:
        solver._result = saved


def avg_reward_via_cpp(solver, timeout=None):
    """
    `getAvgReward`: the steady-state expectation of every DECLARED reward, and
    their names, from `-a reward`.

    Only a reward built from a `Reward.*` template survives the wire: a bare
    lambda has no serializable form, and `linemodel_save` refuses to emit one
    rather than write a reward that would be wrong on reload. So a model whose
    rewards are all custom functions is refused by the writer, before the C++ is
    reached, with that reward named.
    """
    # A REWARD THE WRITER DROPS MUST STOP THE SOLVE, NOT SHORTEN THE ANSWER.
    # `_rewards_to_json` omits a reward built from a bare callable, because a
    # Python closure has no serializable form, and warns. If the solve continues
    # anyway the C++ answers for the SURVIVING rewards only, and the caller pairs
    # that shorter vector with its own full declaration list: on
    # rewardModel_multiclass that raised IndexError out of the example, and on
    # rewardModel_aggregation -- where every reward is a closure -- line-cli was
    # handed a model with no reward at all and exited 2. Refuse by name instead,
    # which is what every other unbridgeable construct does.
    if _unbridgeable_rewards(solver):
        # A CLOSURE CANNOT CROSS THE WIRE, and there is nothing to serialize it
        # into -- but the reward is a function OF THE STATIONARY LAW, and that
        # law is the C++ one: `pi` and `space_aggr` reach the result container
        # through `generator_via_cpp`. Evaluating the caller's own Python
        # function against the C++ chain therefore answers under the caller's
        # `lang` -- the engine that computed the answer is still line-cli, and
        # this module only applies a function the caller supplied and no wire
        # format can express. Refusing instead left rewardModel_multiclass and
        # rewardModel_aggregation, whose every reward is a closure, unsolvable.
        #
        # EVERY reward is evaluated locally here, not just the unbridgeable
        # ones: `-a reward` answers for the serializable subset only, and
        # splicing two vectors computed over the same law adds nothing but a
        # pairing that can rotate. The law is the C++ one either way.
        return _avg_reward_locally(solver, timeout=timeout)

    p = analysis_via_cpp(solver, 'reward', timeout=timeout)
    vals = np.asarray(p.get('E') or [], dtype=float)
    names = [str(x) for x in (p.get('Reward') or [])]

    # THE WIRE SORTS THE REWARDS; THE MODEL DOES NOT. `_rewards_to_json` emits
    # them in NAME order so the three writers byte-match, so the C++ answers in
    # that order while the native path answers in DECLARATION order. A caller
    # that pairs the values with its own declaration list -- which is what the
    # rewardModel examples do -- then reads every reward under the wrong name:
    # on rewardModel_templates that printed QueueLength 0.014699 (the blocking
    # probability) and BlockingProb 0.738976 (the utilization), a rotation, not a
    # numeric error. Restore the model's order here so lang='cpp' is
    # indistinguishable from the native getter.
    declared = getattr(solver.model, '_rewards', None)
    if declared and names and len(names) == vals.size:
        pos = dict((nm, i) for i, nm in enumerate(names))
        order = [pos[nm] for nm in declared.keys() if nm in pos]
        if len(order) == len(names):
            vals = vals[order]
            names = [names[i] for i in order]
    return vals, names


def _unbridgeable_rewards(solver):
    """The declared rewards no wire format can carry, by name.

    A reward survives `_rewards_to_json` only as a `Reward.*` template: a bare
    callable, or a descriptor of kind Custom, has no serializable form and the
    writer omits it with a warning. The list is what decides between the two
    reward paths below, and is shared by the steady-state and transient arms so
    they cannot disagree about which model is bridgeable.
    """
    declared = getattr(solver.model, '_rewards', None) or {}
    try:
        from ..lang.reward import RewardDescriptor
    except ImportError:
        return []
    return [nm for nm, fn in declared.items()
            if not isinstance(fn, RewardDescriptor) or fn.kind == 'Custom' or fn.node is None]


def tran_reward_at_via_cpp(solver, t, reward_vector, timeout=None):
    """
    `getTranReward(t, reward_vector)`: the scalar E[r(X(t))] of a reward given as
    a vector over the DETAILED state space.

    A vector indexed by state cannot cross the wire under a state ordering the
    caller has not seen, so the enumeration is fetched first (`-a states`, via
    `generator_via_cpp`) and the default reward -- the total jobs in each state
    -- is built from it. The occupancy itself is line-cli's pi(t); this function
    only pairs the two, which is the same division of labour as
    `_avg_reward_locally`.
    """
    law = generator_via_cpp(solver, timeout=timeout)
    space = law.get('space')
    saved = getattr(solver, '_result', None)
    try:
        solver._result = _CppLaw(law.get('pi'), space, law.get('space_aggr'))
        return solver._compute_tran_reward(t, reward_vector)
    finally:
        solver._result = saved


def tran_reward_via_cpp(solver, rewards_dict, name=None, timeout=None):
    """
    `getTranReward`/`get_tran_reward`: E[r(X(t))] over `options.timespan`.

    TWO ROUTES, ONE LAW, chosen exactly as `avg_reward_via_cpp` chooses. Where
    every reward is a `Reward.*` template the C++ carries the reward map itself
    and `-a tranreward` returns the trajectories; where any reward is a bare
    callable nothing can serialize it, so `-a tranprob` supplies the transient
    law and its AGGREGATED labels and the caller's own functions are applied to
    them here. Either way the chain, the integration and the initial
    distribution are line-cli's -- only the reward map is local, and only when
    no wire format can express it.

    THE HORIZON IS THE CALLER'S, as it is for `-a tran`: line-cli refuses an
    unstated one rather than inventing a bound, so it is checked here where
    `options.timespan` can be named.

    Returns the native getter's triple `(Rt, t, names)`, or `(Rt[i], t,
    names[i])` when `name` selects one reward.
    """
    ts = getattr(getattr(solver, 'options', None), 'timespan', None)
    ts = np.asarray(ts, dtype=float).ravel() if ts is not None else np.asarray([])
    if ts.size != 2 or not np.all(np.isfinite(ts)) or not (ts[0] < ts[1]):
        raise RuntimeError(
            "getTranReward under lang='cpp' needs a finite horizon: set options.timespan = "
            "[t0, t1] with t0 < t1 (for example CTMC(model, timespan=[0, 5])).")
    span = ['--tspan', '%r:%r' % (float(ts[0]), float(ts[1]))]

    names = list(rewards_dict.keys())
    if _unbridgeable_rewards(solver):
        p = analysis_via_cpp(solver, 'tranprob', timeout=timeout, flags=span)
        t = np.asarray(p.get('t') or [], dtype=float).ravel()
        pit = np.atleast_2d(np.asarray(p.get('pitAggr') or [], dtype=float))
        space = np.atleast_2d(np.asarray(p.get('labelsAggr') or [], dtype=float))
        if t.size == 0 or pit.size == 0 or space.size == 0:
            raise RuntimeError(
                "line-cli returned no transient law for this model, so the reward(s) defined by "
                "a bare callable cannot be evaluated against it.")
        if pit.shape[1] != space.shape[0]:
            raise RuntimeError(
                "line-cli returned a %d-column transient law over a %d-row aggregated state "
                "space; a reward cannot be paired with it." % (pit.shape[1], space.shape[0]))
        names, Rmat = solver.reward_matrix_over(space, rewards_dict, space.shape[0])
        metrics = pit @ Rmat.T                      # (ntimes x nrewards)
        Rt = [{'t': t, 'metric': np.asarray(metrics[:, ri]).ravel(), 'name': names[ri]}
              for ri in range(len(names))]
    else:
        p = analysis_via_cpp(solver, 'tranreward', timeout=timeout, flags=span)
        t = np.asarray(p.get('t') or [], dtype=float).ravel()
        E = [np.asarray(row, dtype=float).ravel() for row in (p.get('E') or [])]
        wire = [str(x) for x in (p.get('Reward') or [])]
        if len(wire) != len(E):
            raise RuntimeError(
                "line-cli answered -a tranreward with %d names and %d trajectories."
                % (len(wire), len(E)))
        # THE WIRE SORTS THE REWARDS; THE MODEL DOES NOT -- the same rotation
        # `avg_reward_via_cpp` undoes, and for the same reason: the examples pair
        # the answer with their own declaration list.
        pos = dict((nm, i) for i, nm in enumerate(wire))
        order = [pos[nm] for nm in names if nm in pos]
        if len(order) != len(wire):
            names, order = wire, list(range(len(wire)))
        Rt = [{'t': t, 'metric': E[order[ri]], 'name': names[ri]} for ri in range(len(names))]

    if name is not None:
        if name not in names:
            raise ValueError('Reward "%s" not found. Available rewards: %s'
                             % (name, ', '.join(names)))
        idx = names.index(name)
        return Rt[idx], t, names[idx]
    return Rt, t, names


def tran_prob_via_cpp(solver, t, node=None, timeout=None):
    """
    `getTranProb`/`getTranProbAggr`/`getTranProbSys`/`getTranProbSysAggr`, from
    `-a tranprob` integrated over [0, t].

    ONE INTEGRATION SERVES ALL FOUR: the CLI takes both label sets off the same
    forward solve, so the detailed and the aggregate view of one horizon cost one
    process and cannot come from two different integrations.

    The native getters return pi AT t, while the CLI returns the whole trajectory
    over the horizon; the LAST row is the answer, and the time it was taken at is
    returned with it so a caller can see what it actually got rather than trust
    that the integrator landed on t.
    """
    t = float(t)
    if not (t > 0.0):
        raise RuntimeError(
            "getTranProb* under lang='cpp' needs a positive horizon (got %r): pi(t) on an "
            "unstated or empty span is not a quantity, and line-cli refuses --tspan 0:0." % t)
    flags = ['--tspan', '0:%s' % repr(t)]
    if node is not None:
        flags += ['--node', str(int(node) + 1)]
    p = analysis_via_cpp(solver, 'tranprob', timeout=timeout, flags=flags)
    times = np.asarray(p.get('t') or [], dtype=float)
    return {
        't': times,
        'tlast': float(times[-1]) if times.size else float('nan'),
        'pit': np.atleast_2d(np.asarray(p.get('pit') or [], dtype=float)),
        'pitAggr': np.atleast_2d(np.asarray(p.get('pitAggr') or [], dtype=float)),
        'labels': np.atleast_2d(np.asarray(p.get('labels') or [], dtype=float)),
        'labelsAggr': np.atleast_2d(np.asarray(p.get('labelsAggr') or [], dtype=float)),
    }


def tran_avg_via_cpp(solver, timeout=None):
    """
    `getTranAvg`, from `-a tran` over `options.timespan`.

    Returns (QNt, UNt, TNt), each an (nstations x nclasses) nested list of
    TranResult(t, metric) -- the same objects the native getter returns, so a
    caller reads `QNt[i][r].metric[-1]` without knowing which engine answered.

    THE HORIZON IS THE CALLER'S. line-cli refuses `-a tran` without `--tspan`
    rather than inventing one, so an unbounded or unstated `options.timespan` is
    refused HERE, where the option can be named, instead of arriving as an argv
    error. This is the twin of MATLAB's CPPLINE.tranAvg, down to the payload
    reading, and both leave the metricVal wrapping to their own callers.

    A (station, class) pair the engine sends no curve for stays None, exactly as
    the native getter leaves a disabled handle: an absent curve is not a zero
    one, and filling it here would invent the only series in the table.
    """
    from ..constants import TranResult

    ts = getattr(getattr(solver, 'options', None), 'timespan', None)
    ts = np.asarray(ts, dtype=float).ravel() if ts is not None else np.asarray([])
    if ts.size != 2 or not np.all(np.isfinite(ts)) or not (ts[0] < ts[1]):
        raise RuntimeError(
            "getTranAvg under lang='cpp' needs a finite horizon: set options.timespan = "
            "[t0, t1] with t0 < t1 (for example Fluid(model, timespan=[0, 50])).")
    p = analysis_via_cpp(solver, 'tran', timeout=timeout,
                         flags=['--tspan', '%r:%r' % (float(ts[0]), float(ts[1]))])

    sn = _resolve_struct(solver)
    M, K = int(sn.nstations), int(sn.nclasses)
    QNt = [[None] * K for _ in range(M)]
    UNt = [[None] * K for _ in range(M)]
    TNt = [[None] * K for _ in range(M)]
    for e in (p.get('curves') or []):
        # indexBase is 0 on the wire, as everywhere else in this protocol.
        i, r = int(e['station']), int(e['jobclass'])
        t = np.asarray(e.get('t') or [], dtype=float).ravel()
        QNt[i][r] = TranResult(t, np.asarray(e.get('QLen') or [], dtype=float).ravel())
        UNt[i][r] = TranResult(t, np.asarray(e.get('Util') or [], dtype=float).ravel())
        TNt[i][r] = TranResult(t, np.asarray(e.get('Tput') or [], dtype=float).ravel())
    return QNt, UNt, TNt


def sample_path_via_cpp(solver, numEvents, node=None, seed=None, timeout=None):
    """
    `sampleSys`/`sampleSysAggr`, and with `node` also `sample`/`sampleAggr`, from
    `-a sample`.

    THE SEED IS PART OF THE ANSWER. A trajectory is a measurement only with the
    stream that produced it, so the seed actually used comes back in the result
    rather than being left implicit in the CLI's default.
    """
    _assert_default_state(solver, 'sample/sampleSys')
    flags = ['--samples', str(int(numEvents))]
    if seed is not None:
        flags += ['--seed', str(int(seed))]
    if node is not None:
        flags += ['--node', str(int(node) + 1)]
    p = analysis_via_cpp(solver, 'sample', timeout=timeout, flags=flags)
    ev = [None if x is None else int(x) for x in (p.get('event') or [])]
    idx = np.asarray(p.get('state') or [], dtype=int)
    space = np.atleast_2d(np.asarray(p.get('space') or [], dtype=float))
    out = {
        't': np.asarray(p.get('t') or [], dtype=float),
        'state': idx,
        # The states THEMSELVES, from the space the payload carries: a trajectory
        # of indices into a chain the caller cannot see is not a trajectory.
        'stateRows': space[idx, :] if space.size and idx.size else np.zeros((0, 0)),
        'space': space,
        'event': ev,
        'sysAggr': np.atleast_2d(np.asarray(p.get('sysAggr') or [], dtype=float)),
        'seed': int(p.get('seed', 0)),
        'events': int(p.get('events', 0)),
    }
    if node is not None:
        out['nodeState'] = np.atleast_2d(np.asarray(p.get('nodeState') or [], dtype=float))
        out['nodeAggr'] = np.atleast_2d(np.asarray(p.get('nodeAggr') or [], dtype=float))
    return out


_CACHE_COLS = ['Node', 'JobClass', 'List', 'ListCap', 'Items', 'HitProb',
               'DelayedHitProb', 'MissProb', 'HitRate', 'DelayedHitRate',
               'MissRate', 'ArvR', 'ResidT', 'ListCost']
_ITEM_COLS = ['Node', 'Item', 'List', 'ListCap', 'Size', 'Prob', 'Cost',
              'DelayedHitQLen', 'DelayedHitQLenFull']
# Columns of the two tables that carry names rather than metrics.
_TABLE_STRING_COLS = ('Node', 'JobClass')


def _table_via_cpp(solver, analysis, cols, timeout=None):
    """Rebuild one of the cache tables as the native getters' IndexedTable."""
    import pandas as pd

    from ..indexed_table import IndexedTable

    p = analysis_via_cpp(solver, analysis, timeout=timeout)
    if cols[0] not in p:
        raise RuntimeError(
            "line-cli did not return a structured '%s' table (keys: %s)."
            % (analysis, sorted(p.keys())))
    n = len(p.get(cols[0]) or [])
    df = pd.DataFrame([{c: (p.get(c) or [None] * n)[i] for c in cols} for i in range(n)],
                      columns=cols)
    # A JSON null must become NaN, as in the native tables: left alone it arrives
    # as None, pandas types the whole column as object, and arithmetic on the
    # result silently breaks.
    for c in cols:
        if c not in _TABLE_STRING_COLS:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    table = IndexedTable(df)
    if not getattr(solver, '_table_silent', False):
        print(table)
    return table


def cache_table_via_cpp(solver, timeout=None):
    """
    `getAvgCacheTable`, from `-a cache`.

    DELEGATED RATHER THAN REBUILT, because the native builder reads the hit,
    delayed-hit and miss probabilities off the Cache NODE OBJECTS, which only a
    native solve writes. Under `lang='cpp'` nothing writes them, so building the
    table here reported a plain cache -- hit and miss summing to one and the
    delayed-hit column flat zero -- for a model that has a retrieval system.
    That is the failure mode this bridge exists to prevent: numbers from one
    engine under another engine's `lang`.
    """
    return _table_via_cpp(solver, 'cache', _CACHE_COLS, timeout=timeout)


def item_table_via_cpp(solver, timeout=None):
    """`getAvgItemTable`, from `-a item`. Same reasoning as the cache table."""
    return _table_via_cpp(solver, 'item', _ITEM_COLS, timeout=timeout)


def export_odes_via_cpp(solver, notation='scalar', timeout=None):
    """
    `exportODEs`: the LaTeX document of the fluid drift, from `-a odes`.

    `--notation` REACHES THE EXPORTER, so 'matrix' returns the matrix document
    and not the scalar one under another name. A method whose drift has no single
    exported form -- `rmf`, where a Cache model's routing is rewritten between
    sweeps -- is refused by the C++ by name, as is a method whose symbolic export
    is unported (`minnormal`, which is what `default` resolves to on many models,
    so an explicit `method=` is often needed to reach this at all).
    """
    p = analysis_via_cpp(solver, 'odes', timeout=timeout,
                         flags=['--notation', str(notation)])
    tex = p.get('latex')
    if not isinstance(tex, str) or not tex.strip():
        raise RuntimeError("line-cli's -a odes payload carries no LaTeX document.")
    return tex


# --- layered networks (SolverLN) -------------------------------------------
#
# THE LAYERED PATH HAS TWO WIRE FORMATS AND PICKS BY WHAT THE MODEL CARRIES.
# `.lqnx` is the default: it is the LQNS file format, so the same bytes are what
# an external tool would read, and line-cli's `lqn_reader` covers all of it. It
# is however LOSSIER than the JSON interchange -- the schema accepts think-time
# on a REFERENCE task only, and `writeXML` warns and drops a non-reference
# task's -- so a model carrying one is serialized with `save_model` instead and
# read back by line-cli's `lqn_json_reader`, which takes think time on any task.
# The CLI decides by content (`is_layered_json`), so no `-i` is passed with it.
#
# The choice is per-model rather than always-JSON because the two readers do not
# cover the same models: `lqn_json_reader` refuses `fanIn`/`fanOut` replication,
# which the .lqnx reader carries. Defaulting to .lqnx keeps those solvable and
# moves only the models the XML cannot express.
#
# The C++ LnOptions defaults (iter_max=200, iter_tol=5e-3, interlocking, relax
# 'fixed'/0.5, srvn layering) are the same as the native SolverLN defaults, so an
# untouched option is not forwarded and the two engines start from one setting.

_LN_DEFAULT_ITER_MAX = 200
_LN_DEFAULT_ITER_TOL = 5e-3


def _ln_layer_solver(solver):
    """
    Resolve the C++ `--layer-solver` value from the layer solvers the native
    SolverLN actually built, refusing anything the port has no layer engine for.

    The C++ layered path runs its layers under SolverMVA, SolverNC, SolverFluid
    or SolverSSA (`solve_layer_*` in cpp/include/line/solvers/ln/solver_ln.h),
    whereas the native factory may also hand back CTMC or LDES. Those are
    REFUSED BY NAME rather than served with MVA layers: the layer solver is what
    the fixed point is a fixed point of, so substituting it answers a different
    question. NC was refused here long after the C++ engine gained
    `--layer-solver nc`, and the stale refusal, not a missing engine, is why
    lqn_twotasks and lqn_ofbiz reported no table at all on the P2C row.
    """
    names = []
    cls = getattr(solver, '_layer_solver_cls', None)
    if cls is not None:
        # THE FACTORY IS THE DECLARATION; the built list is the declaration AFTER
        # per-layer substitution. When the caller named a class, that class alone
        # decides the ensemble knob.
        return _ln_resolve_one(getattr(cls, '__name__', str(cls)))
    for probe in (getattr(solver, 'solvers', None) or []):
        if probe is None:
            continue
        nm = None
        getname = getattr(probe, 'getName', None)
        if getname is not None:
            try:
                nm = getname()
            except Exception:
                nm = None
        names.append(nm if nm else type(probe).__name__)
    # No layers built yet -- the normal case here, since a delegated solve never
    # enters iterate(). Ask the solver to name what its FACTORY produces before
    # falling back to the native default: a lambda factory records no
    # `_layer_solver_cls`, so defaulting to MVA ran `--layer-solver mva` for
    # every `LN(model, lambda m: NC(m, opts))`.
    if not names:
        probe = getattr(solver, 'probe_layer_solver_name', None)
        probed = probe() if probe is not None else None
        names.append(probed if probed else 'MVA')

    resolved = set()
    for nm in names:
        n = str(nm).replace('Solver', '').upper()
        if n in ('MAM', 'CTMC'):
            # NOT A CHOICE, A SUBSTITUTION, and one the C++ makes for itself:
            # `solve_layer` sends a layer carrying a finite capacity region to
            # CTMC and a layer carrying a SetupTask's setup / delay-off times to
            # MAM, exactly as the native SolverLN does before it reaches the
            # user's factory. Reading it back off the BUILT layer list as if the
            # caller had asked for it refused lqn_setup, whose factory returns
            # MVA and whose host layer the substitution moved to MAM.
            continue
        if n == 'MVA':
            resolved.add('mva')
        elif n in ('FLD', 'FLUID'):
            resolved.add('fluid')
        elif n in ('NC', 'COMOM'):
            resolved.add('nc')
        elif n == 'SSA':
            resolved.add('ssa')
        else:
            raise RuntimeError(
                "lang='cpp' runs the LQN layers under SolverMVA, SolverNC, SolverFluid "
                "or SolverSSA; this SolverLN builds a '%s' layer solver, which the C++ "
                "port has no layered engine for. Solve it with lang='python' (or "
                "lang='java'), or build the SolverLN with one of those layer "
                "factories." % nm)
    if len(resolved) > 1:
        raise RuntimeError(
            "lang='cpp' takes ONE layer solver for the whole ensemble (the C++ "
            "--layer-solver knob is per-ensemble, like the reference's factory argument), "
            "but this SolverLN mixes %s across layers." % sorted(resolved))
    # Every layer was an automatic substitution, so the ensemble knob keeps its
    # default; the C++ re-derives the same per-layer choice itself.
    return resolved.pop() if resolved else 'mva'


def _ln_resolve_one(nm):
    """The `--layer-solver` value for ONE named layer-solver class."""
    n = str(nm).replace('Solver', '').upper()
    if n == 'MVA':
        return 'mva'
    if n in ('FLD', 'FLUID'):
        return 'fluid'
    if n in ('NC', 'COMOM'):
        return 'nc'
    if n == 'SSA':
        return 'ssa'
    raise RuntimeError(
        "lang='cpp' runs the LQN layers under SolverMVA, SolverNC, SolverFluid or "
        "SolverSSA; this SolverLN builds a '%s' layer solver, which the C++ port has no "
        "layered engine for. Solve it with lang='python' (or lang='java'), or build the "
        "SolverLN with one of those layer factories." % nm)


def _lqnx_lossy_tasks(solver):
    """
    Name the tasks whose think time an `.lqnx` serialization of this model would
    drop, i.e. the non-reference ones that carry one.

    The lqnx schema accepts think-time on a REFERENCE task only, and `writeXML`
    reports the loss and writes the file anyway (lqns would reject it
    otherwise). LINE does give a non-reference task's think time to its callers,
    so a model carrying one solves to DIFFERENT numbers through that transport:
    on gallery_lqn_basic, T3's think time caps its throughput at
    multiplicity/think = 25/4 and drops it from 66.4 to 6.22. A non-empty
    return therefore does not refuse the solve -- it routes it through
    `save_model` instead, whose reader takes think time on any task.
    """
    lqn = getattr(solver, 'lqn', None)
    if lqn is None:
        return []
    think = getattr(lqn, 'think', None)
    if think is None:
        return []
    items = think.items() if isinstance(think, dict) else enumerate(think)
    offenders = []
    for idx, mean_or_dist in items:
        if mean_or_dist is None:
            continue
        idx = int(idx)
        # A non-dict `think` is 0-indexed with position 0 unused, matching how
        # _construct() reads it.
        if not isinstance(think, dict):
            idx += 1
        if idx <= int(getattr(lqn, 'tshift', 0)):
            continue
        if isinstance(mean_or_dist, (int, float)):
            mean = float(mean_or_dist)
        else:
            getter = getattr(mean_or_dist, 'getMean', getattr(mean_or_dist, 'get_mean', None))
            try:
                mean = float(getter()) if getter is not None else 0.0
            except Exception:
                mean = 0.0
        if mean > 0 and not solver._is_ref_task(idx):
            offenders.append(solver._get_hashname(idx))
    return offenders


def _ln_knobs(solver):
    """
    Translate the native SolverLN options into layered-path CLI flags, refusing
    every setting the flag set cannot carry.

    A knob with no CLI counterpart is an error and not a silent drop: `relax`,
    `relax_factor` and `layering` all change the fixed point the native solver
    converges to, so running without them would answer a different question
    under the caller's `lang`.
    """
    args = []
    opts = getattr(solver, 'options', None)
    if opts is None:
        return args

    method = getattr(opts, 'method', None)
    method = str(method).lower() if method is not None else 'default'
    if method not in ('', 'default'):
        # --method is refused by line-cli on the layered path (it names an
        # algorithm inside a Network solver), and the native LN methods that are
        # not 'default' -- 'nc', the 'mwba.*' bounds -- are separate engines
        # rather than options of this one.
        raise RuntimeError(
            "lang='cpp' serves the layered solver's default method only; options.method='%s' "
            "selects a different layered engine (the C++ layered path takes no --method). Use "
            "lang='python' for it." % method)

    iter_tol = _forward_if_overridden(getattr(opts, 'iter_tol', None), _LN_DEFAULT_ITER_TOL)
    if iter_tol is not None:
        args += ['--iter_tol', repr(float(iter_tol))]
    im = _forward_if_overridden(getattr(opts, 'iter_max', None), _LN_DEFAULT_ITER_MAX)
    if im is not None:
        try:
            im = int(im)
        except (TypeError, ValueError):
            im = None
        if im is not None and im > 0:
            args += ['--iter_max', str(im)]

    cfg = getattr(opts, 'config', None)
    if cfg is not None:
        def _cfg(key, default):
            try:
                v = cfg.get(key, default)
            except AttributeError:
                v = getattr(cfg, key, default)
            return default if v is None else v

        if not bool(_cfg('interlocking', True)):
            args.append('--no-interlocking')
        layering = str(_cfg('layering', 'srvn')).lower()
        if layering != 'srvn':
            raise RuntimeError(
                "lang='cpp' implements the 'srvn' layering (one submodel per server); "
                "config['layering']='%s' builds a different layer decomposition, which the "
                "C++ port does not carry. Use lang='python' for it." % layering)
        relax = str(_cfg('relax', 'fixed')).lower()
        if relax != 'fixed':
            raise RuntimeError(
                "lang='cpp' runs the C++ layered solver's fixed under-relaxation; "
                "config['relax']='%s' has no CLI counterpart and changes the iterate, so it "
                "cannot be honoured here. Use lang='python' for it." % relax)
        try:
            relax_factor = float(_cfg('relax_factor', 0.5))
        except (TypeError, ValueError):
            relax_factor = 0.5
        if abs(relax_factor - 0.5) > 1e-12:
            raise RuntimeError(
                "lang='cpp' cannot forward config['relax_factor']=%g (the layered CLI exposes "
                "no relaxation flag), and the C++ solver would relax by 0.5 instead. Use "
                "lang='python' for it." % relax_factor)
    return args


def solve_lqn_via_cpp(solver, timeout=None, arith=None):
    """
    Serialize the solver's LayeredNetwork, run `line-cli` on it and return the
    parsed JSON object ({model, arith, layers, iterations, converged, seconds,
    rows}).

    The wire format is `.lqnx` unless that would drop a non-reference task's
    think time, in which case the model.json interchange carries it instead --
    see `_lqnx_lossy_tasks`. Caches, item entries and setup / delay-off times
    ride the `.lqnx` dialect itself (see LayeredNetwork.writeXML).

    ONE subprocess solves the whole ensemble, as in `ln_ensemble_avg_via_jar`:
    letting the native fixed point run and dispatching each layer separately
    would spawn one process per layer per iteration.
    """
    binary = find_line_cli()
    model = _get_model(solver)
    lossy = _lqnx_lossy_tasks(solver)
    layer_solver = _ln_layer_solver(solver)

    opts = getattr(solver, 'options', None)
    if arith is None and opts is not None:
        arith = getattr(opts, 'arith', None)
    if arith and layer_solver == 'fluid' and str(arith) != 'double':
        # line-cli refuses this itself; naming it here says which of the two
        # settings has to change instead of surfacing a bare CLI refusal.
        raise RuntimeError(
            "the C++ fluid layer solver is double-precision only; arith='%s' cannot be "
            "combined with Fluid layers." % arith)

    tmpdir = tempfile.mkdtemp(prefix='line_cpp_lqn_')
    try:
        if lossy:
            # No `-i`: line-cli takes the layered path off the document's own
            # `type`, and naming `-i json` would send it to the Network reader.
            from ..io.linemodel_io import save_model
            model_path = os.path.join(tmpdir, 'model.json')
            save_model(model, model_path)
            fmt_args = []
        else:
            model_path = os.path.join(tmpdir, 'model.lqnx')
            model.writeXML(model_path)
            fmt_args = ['-i', 'lqnx']
        # -s ln and -s ln.mva are the same engine in line-cli (the layer solver
        # comes from --layer-solver, not from the token), so the token is left at
        # 'ln' and the layer engine is stated once, explicitly.
        cmd = [binary, '-f', model_path] + fmt_args + ['-s', 'ln', '-a', 'avg', '-o', 'json',
               '--layer-solver', layer_solver]
        if arith:
            cmd += ['--arith', str(arith)]
        cmd += _ln_knobs(solver)
        return _run_line_cli(binary, cmd, timeout=timeout)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def ln_avg_table_via_cpp(solver, timeout=None):
    """
    Delegate a layered getAvgTable() to `line-cli`, returning a DataFrame with
    the native layered columns (Node, NodeType, QLen, Util, RespT, ResidT, ArvR,
    Tput), one row per LQN element in LQN index order.

    ROWS ARE REORDERED INTO NATIVE LQN INDEX ORDER, keyed by element name. The
    two orders are NOT the same: line-cli indexes the elements as the `.lqnx`
    file declares them, which groups tasks under their processor, while the
    native struct numbers tasks in model-declaration order. On a model whose
    tasks are declared in one order and hosted in another (T1 on P1, T2 on P3,
    T3 on P2) the two tables carry the same numbers against different rows, and a
    positional comparison reads that as a 6x utilization error.

    ArvR is all-NaN because the C++ layered result has no arrival-rate vector --
    and so does the native table, whose AN is never filled ("not available yet"),
    so the two agree by construction rather than by omission.
    """
    import numpy as np
    import pandas as pd

    results = solve_lqn_via_cpp(solver, timeout=timeout)
    rows = results.get('rows') if isinstance(results, dict) else None
    if not isinstance(rows, list):
        raise RuntimeError(
            "line-cli did not return a layered AvgTable (keys: %s)."
            % (list(results.keys()) if isinstance(results, dict) else type(results)))
    # The layered fixed point may stop on iter_max; mirror the flag onto the
    # native solver so a delegated solve cannot report convergence it never
    # reached, and record the iteration count the same way iterate() does.
    if 'converged' in results:
        solver.hasconverged = bool(results['converged'])
    if 'iterations' in results:
        try:
            solver.it = int(results['iterations'])
        except (TypeError, ValueError):
            pass

    by_name = {}
    for r in rows:
        by_name[str(r.get('node'))] = r
    lqn = solver.lqn
    out = []
    missing = []
    # LQN element indices run 0..nidx-1, the same space the native table walks
    for idx in range(int(lqn.nidx)):
        name = solver._get_hashname(idx)
        r = by_name.get(name)
        if r is None:
            missing.append(name)
            continue
        kind = str(r.get('type'))
        native_kind = solver._get_type_name(idx)
        if kind != native_kind:
            # The two ports disagree on what this element IS, which is a
            # structural difference in how the model was read, not a numeric one.
            raise RuntimeError(
                "line-cli reports element '%s' as a %s where the native struct has a %s; "
                "the two ports read the .lqnx differently." % (name, kind, native_kind))
        out.append({'Node': name, 'NodeType': native_kind,
                    'QLen': r.get('QLen'), 'Util': r.get('Util'), 'RespT': r.get('RespT'),
                    'ResidT': r.get('ResidT'), 'ArvR': np.nan, 'Tput': r.get('Tput')})
    if missing:
        raise RuntimeError(
            "line-cli returned no row for LQN element(s) %s; the .lqnx interchange did not "
            "carry the whole model." % ', '.join(missing))
    df = pd.DataFrame(out, columns=_LN_COLS)
    # A JSON null is an undefined metric and must read as NaN, as in the native
    # table; left alone it arrives as None and types the column as object.
    for c in ('QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def ln_ensemble_avg_via_cpp(solver, timeout=None):
    """
    Delegate a layered getEnsembleAvg() to `line-cli`, returning
    (QN, UN, RN, TN, AN, WN) as (nidx,) arrays indexed by LQN absolute index,
    matching the native SolverLN.get_ensemble_avg contract.
    """
    import numpy as np

    df = ln_avg_table_via_cpp(solver, timeout=timeout)
    lqn = solver.lqn
    nidx = int(lqn.nidx)
    QN, UN, RN, TN, AN, WN = (np.full(nidx, np.nan) for _ in range(6))
    # Index by element name, not by row order: line-cli emits the elements in
    # LQN index order, and asserting that here rather than assuming it turns a
    # future reordering into missing rows instead of silently shifted metrics.
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
