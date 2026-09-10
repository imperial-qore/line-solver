"""Backend dispatch for the JMT command line.

One analysis of jmt.commandline.Jmt is run on a model file, and the result file
is left where the JMT CLI itself would leave it, that is at
``<model>-result.jsim`` for mode "sim" and ``<model>-result.jmva`` for mode
"mva". Every backend satisfies that contract, so the result parsers do not know
or care which one ran.

Backend selection, in order:

1. ``options.rest_url`` non-empty: POST to a JMT REST server (the
   imperialqore/jmt-rest container). Nothing is executed locally.
2. a local JVM plus common/JMT.jar: the default, unchanged.
3. no local JVM, but Docker is usable: ask once per session whether to pull and
   use the JMT image, and dispatch through it.

see _kb/06-solver-catalog.md (Wrappers: three ways to reach an external binary)
"""

import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from typing import Optional, Tuple

from ...io.viewers import get_java_exe

DEFAULT_IMAGES = ('imperialqore/jmt-rest:latest', 'imperialqore/jmt-rest')

RESULT_SUFFIX = {'sim': '-result.jsim', 'mva': '-result.jmva'}

# Session state: the JVM probe and the Docker consent are both asked once.
_JAVA_EXE = False  # False = not probed yet, None = no JVM, str = launcher path
# How long `java -version` gets to answer. Generous because it is a liveness
# check on a shared host, not a benchmark: a cold JVM behind NFS on a box
# running a parallel suite takes far longer than an idle one. Breaching it is
# not a verdict -- see java_exe.
_JAVA_PROBE_TIMEOUT = 60
_DOCKER_DECISION = None


class JMTBackendError(RuntimeError):
    """Raised when no backend can run JMT, or the chosen backend failed."""


def result_path_for(model_path: str, mode: str) -> str:
    """Return the path JMT writes its result to for this model and mode."""
    if mode not in RESULT_SUFFIX:
        raise ValueError(f"Unknown JMT analysis mode: {mode}")
    return model_path + RESULT_SUFFIX[mode]


def java_exe() -> Optional[str]:
    """Return a runnable Java launcher, or None when no JVM is reachable.

    Resolution is $LINE_JAVA, then $JAVA_HOME/bin/java, then the PATH, so a
    host whose PATH carries no java (the usual case on Windows) is still
    served. Probed once per session: the check costs a process launch, and a
    JVM does not appear or vanish mid-session.

    A SLOW PROBE IS NOT AN ABSENT JVM, and reading it as one is how host load
    became a wrong answer -- the same mistake `run_jmt` names below about its own
    budget. `get_java_exe` returns a launcher that RESOLVED and exists, so the
    only thing `-version` adds is whether it runs; a `TimeoutExpired` says the
    machine is busy, not that the launcher is missing. It used to be swallowed
    with FileNotFoundError and cached as None for the life of the process, and
    every later JMT solve in that process then died on "No JVM was found on the
    path" with java sitting on the PATH all along. Measured 2026-09-07 on
    picard06, the python parity block under `pytest -n 4`: five rows failed that
    way across two runs (fj_deep_nesting, fj_route_overlap, cqn_twoclass_erl,
    oqn_cs_routing, sdroute_twoclasses_closed -- every one of them a golden with
    a JMT row) while 66 other JMT rows in the same run passed, which is what a
    per-process cached false negative looks like. A JVM that is genuinely broken
    still fails, but at the real solve and with java's own message, which is the
    better diagnosis anyway.
    """
    global _JAVA_EXE
    if _JAVA_EXE is False:
        cand = get_java_exe()
        if cand is None:
            _JAVA_EXE = None
        else:
            try:
                proc = subprocess.run([cand, '-version'], capture_output=True,
                                      timeout=_JAVA_PROBE_TIMEOUT)
                _JAVA_EXE = cand if proc.returncode == 0 else None
            except subprocess.TimeoutExpired:
                _JAVA_EXE = cand
            except (FileNotFoundError, OSError):
                _JAVA_EXE = None
    return _JAVA_EXE


def has_java() -> bool:
    """Report whether a JVM is reachable."""
    return java_exe() is not None


def _docker_available() -> bool:
    """Report whether the Docker daemon is reachable."""
    if sys.platform.startswith('win'):
        # Bind-mount dispatch is supported on unix hosts only, as for LQNS.
        return False
    try:
        proc = subprocess.run(['docker', 'info'], capture_output=True, timeout=20)
        return proc.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


def _image_present(image: str) -> bool:
    """Report whether an image is already on the host."""
    try:
        proc = subprocess.run(['docker', 'images', '-q', image],
                              capture_output=True, timeout=20)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False
    return proc.returncode == 0 and bool(proc.stdout.strip())


def _ask_docker_consent(image: str) -> bool:
    """Ask whether the JMT Docker image may be pulled and used.

    LINE_JMT_DOCKER answers for the user in unattended runs: "1" consents, "0"
    refuses. Without it, a non-interactive session refuses rather than blocking:
    a solver must never hang on a prompt nobody can answer, and must never
    download without being asked.
    """
    env = os.environ.get('LINE_JMT_DOCKER', '').strip().lower()
    if env:
        return env in ('1', 'true', 'yes', 'y')

    if not sys.stdin or not sys.stdin.isatty():
        return False

    print("\nSolverJMT needs Java, which was not found on this host.")
    print(f"Docker is available and can run JMT from the image {image} instead.")
    try:
        answer = input("Pull and use that image? [y/N]: ")
    except EOFError:
        return False
    return answer.strip().lower() in ('y', 'yes')


def resolve_docker_image(options=None) -> Optional[str]:
    """Resolve the JMT Docker image to dispatch through, with user consent.

    Returns None when Docker is unusable, when no image can be obtained, or
    when the user declines. The decision is remembered for the session: a sweep
    over many models must not ask once per model.
    """
    global _DOCKER_DECISION

    if not _docker_available():
        return None

    requested = getattr(options, 'container', None) if options is not None else None
    if not requested:
        requested = os.environ.get('LINE_JMT_IMAGE', '')
    candidates = (requested,) if requested else DEFAULT_IMAGES

    # An image already on the host needs no pull and no question.
    for image in candidates:
        if _image_present(image):
            return image

    if _DOCKER_DECISION == 'no':
        return None

    target = candidates[0]
    if _DOCKER_DECISION != 'yes' and not _ask_docker_consent(target):
        _DOCKER_DECISION = 'no'
        return None
    _DOCKER_DECISION = 'yes'

    # Storage check before pulling (shared guard used by the LQNS/QNS wrappers).
    from ...io import docker_util
    if not docker_util.has_storage_for(target):
        _DOCKER_DECISION = 'no'
        print(f"[LINE] Skipping docker pull of {target}: insufficient free space at the "
              "Docker storage location. Install Java, or set options.rest_url to a running "
              "JMT REST server.", file=sys.stderr)
        return None

    print(f"Pulling {target} (this happens once)...")
    try:
        proc = subprocess.run(['docker', 'pull', target], timeout=1800)
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError) as exc:
        raise JMTBackendError(f"Could not pull {target}: {exc}") from None
    if proc.returncode != 0:
        _DOCKER_DECISION = 'no'
        raise JMTBackendError(
            f"Could not pull {target}. Install Java, or set options.rest_url to a "
            "running JMT REST server.")
    return target


def run_jmt(mode: str, model_path: str, seed: int, jmt_jar: Optional[str] = None,
            options=None, timeout: Optional[float] = None,
            verbose: bool = False) -> Tuple[int, str, str]:
    """Run one JMT analysis and leave the result next to the model.

    Args:
        mode: "sim" (JSIM) or "mva" (JMVA)
        model_path: path of the model file to analyse
        seed: simulation seed, ignored by the JMVA engine
        jmt_jar: path of JMT.jar, required by the local backend
        options: solver options, read for rest_url and container
        timeout: wall-clock budget in seconds, or None for the default
        verbose: print the command that is run

    Returns:
        (returncode, stdout, stderr) of the backend that ran.

    Raises:
        JMTBackendError: if no backend can run JMT, or the remote one failed.
    """
    if mode not in RESULT_SUFFIX:
        raise ValueError(f"Unknown JMT analysis mode: {mode}")

    rest_url = getattr(options, 'rest_url', None) if options is not None else None
    if rest_url:
        return _run_rest(rest_url, mode, model_path, seed, timeout, verbose)

    # NO BUDGET MEANS NO BUDGET. This used to fall back to 600 s, which is not a
    # safety net: a jsim run that is merely SLOW -- a big model, or a loaded host --
    # was killed at ten minutes and its caller handed an empty table, so host load
    # turned into a wrong answer. MATLAB runs the local JVM arm untimed and the JAR
    # defaults simulationTimeoutSeconds to 0 ("No timeout by default"); only a
    # FINITE options.timeout, which the caller asked for, is a budget here.
    sub_timeout = timeout if (timeout is not None and math.isfinite(timeout)
                              and timeout > 0) else None

    java = java_exe()
    if java is not None:
        if not jmt_jar:
            raise JMTBackendError("A local JVM was found but JMT.jar was not supplied")
        cmd = [java, '-cp', jmt_jar, 'jmt.commandline.Jmt', mode, model_path,
               '-seed', str(int(seed))]
        if verbose:
            print("SolverJMT command: " + ' '.join(cmd))
        proc = subprocess.run(cmd, capture_output=True,
                              cwd=os.path.dirname(model_path) or None,
                              timeout=sub_timeout)
        return (proc.returncode,
                proc.stdout.decode('utf-8', errors='ignore'),
                proc.stderr.decode('utf-8', errors='ignore'))

    image = resolve_docker_image(options)
    if not image:
        raise JMTBackendError(
            "SolverJMT requires a Java runtime and JMT.jar. No JVM was found on the "
            "path, and Docker is not usable either. Install Java, or start a JMT REST "
            "server and set options.rest_url.")

    return _run_docker(image, mode, model_path, seed, sub_timeout, verbose)


def _run_docker(image: str, mode: str, model_path: str, seed: int,
                timeout: float, verbose: bool) -> Tuple[int, str, str]:
    """Run the analysis in the JMT container and copy the result back.

    The model is staged under HOME rather than run in place: snap-confined
    Docker cannot bind-mount the system temp dir where the JMT model normally
    lives. In mva mode JMT rewrites the model file itself, which is why the
    staged copy, not the original, is what the container touches.
    """
    stage_root = os.path.join(os.path.expanduser('~'), '.line', 'line_workspace', 'jmt-docker')
    os.makedirs(stage_root, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix='tmp_', dir=stage_root)
    try:
        base = os.path.basename(model_path)
        staged = os.path.join(workdir, base)
        shutil.copyfile(model_path, staged)
        os.chmod(staged, 0o644)

        cmd = ['docker', 'run', '--rm',
               '--user', f"{os.getuid()}:{os.getgid()}",
               '-v', f"{workdir}:{workdir}", '-w', workdir,
               image, mode, base, '-seed', str(int(seed))]
        if verbose:
            print("SolverJMT command: " + ' '.join(cmd))
        proc = subprocess.run(cmd, capture_output=True, timeout=timeout)

        staged_result = result_path_for(staged, mode)
        if os.path.isfile(staged_result):
            shutil.copyfile(staged_result, result_path_for(model_path, mode))

        return (proc.returncode,
                proc.stdout.decode('utf-8', errors='ignore'),
                proc.stderr.decode('utf-8', errors='ignore'))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _run_rest(rest_url: str, mode: str, model_path: str, seed: int,
              timeout: Optional[float], verbose: bool) -> Tuple[int, str, str]:
    """Solve through a JMT REST server and write its result next to the model.

    The request carries the same JSIM or JMVA document the CLI reads, and the
    response carries the same result document the CLI writes, so a fixed seed
    gives the same numbers as the local backend.
    """
    url = rest_url.rstrip('/')
    if not re.search(r'/api/v\d+/solve/(sim|mva)$', url):
        url = url + '/api/v1/solve/' + mode

    with open(model_path, 'r') as f:
        model_text = f.read()
    payload = {'model': {'content': model_text, 'base64': False}}
    # JMVA takes its algorithm and tolerance from the model document, so the
    # seed is only meaningful for the simulation route; sending it to
    # /solve/mva would be rejected by the server's option allow-list.
    if mode == 'sim':
        payload['options'] = {'seed': int(seed)}

    if verbose:
        print(f"SolverJMT REST: POST {url}")

    req_timeout = (timeout + 30.0) if (timeout is not None and math.isfinite(timeout)
                                       and timeout > 0) else 3600
    req = urllib.request.Request(url, data=json.dumps(payload).encode('utf-8'),
                                 method='POST',
                                 headers={'Content-Type': 'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=req_timeout) as resp:
            response = json.loads(resp.read().decode('utf-8'))
    except urllib.error.HTTPError as e:
        detail = e.read().decode('utf-8', errors='ignore')
        raise JMTBackendError(f"JMT REST solve failed (HTTP {e.code}): {detail}") from None
    except urllib.error.URLError as e:
        raise JMTBackendError(f"JMT REST request to {url} failed: {e.reason}") from None

    if response.get('status') != 'completed':
        raise JMTBackendError("JMT REST solve failed: %s"
                              % (response.get('error') or 'unspecified error'))

    raw = response.get('raw_output') or {}
    result_xml = raw.get('result_xml')
    if not result_xml:
        raise JMTBackendError(
            "JMT REST response carries no result document. The server was asked to "
            "include the raw output; check that include_raw_output is not disabled.")

    with open(result_path_for(model_path, mode), 'w') as f:
        f.write(result_xml)

    return 0, raw.get('stdout') or '', raw.get('stderr') or ''
