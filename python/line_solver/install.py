"""
Environment check for the native Python LINE package.

Mirrors the MATLAB ``lineInstall`` script: it verifies that the optional
external dependencies used by some solvers are reachable and warns, without
failing, when one is missing -- a Java runtime, the LQNS binaries behind the
LQNS wrapper (lqns, lqsim and qnsolver), ``JMT.jar`` behind the JMT solver, and the SageMath
symbolic backend (the ``imperialqore/line-sage-rest`` Docker image, required by
the symbolic methods of SolverCTMC/SolverFluid).

The JMT probe asks where JMT already is rather than calling the resolver that
fetches it, so that a check reports a missing JMT instead of spending 50MB
acquiring one. An environment check must not change the environment it is
reporting on.

Run it as ``line-install`` on the command line or call ``line_install()``.
"""

import shutil
import subprocess
import warnings

from line_solver.api.sym.sage_client import DOCKER_IMAGES, _find_image
from line_solver.api.solvers.jmt import jmt_find_path


def _run_ok(cmd):
    """True if the command runs and exits zero, False on any failure."""
    try:
        out = subprocess.run(cmd, capture_output=True, timeout=30)
        return out.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def line_install():
    """Check optional dependencies and warn about missing ones.

    Returns True when everything is in place, False when at least one warning
    was issued. Never raises: a missing dependency degrades a subset of solvers
    but leaves the rest of LINE fully usable.
    """
    has_warnings = False

    print("Checking Java...")
    if shutil.which("java") is None:
        warnings.warn(
            "the Java Runtime Environment (JRE) is not on PATH, this is required "
            "by lang='java' and the JMT/LDES solvers.")
        has_warnings = True

    print("Checking LQNS...")
    if shutil.which("lqns") is None:
        warnings.warn(
            "LQNS is not on PATH, so SolverLQNS cannot run. It needs the "
            "lqns, lqsim and qnsolver "
            "commands. Download them at: https://github.com/layeredqueuing/V6")
        has_warnings = True

    print("Checking JMT...")
    jmt_path = jmt_find_path()
    if jmt_path is None:
        warnings.warn(
            "JMT.jar was not found, so the JMT simulation solver cannot run "
            "yet. It is about 50MB and is downloaded automatically on the "
            "first call to the solver; set LINE_JMT_JAR to use an existing "
            "copy instead, or LINE_JMT_DOWNLOAD=0 to refuse the download.")
        has_warnings = True
    else:
        print("  %s" % jmt_path)

    print("Checking C++ engine (line-cli)...")
    # _find_existing_line_cli, not find_line_cli: asking whether lang='cpp' can
    # run must not be what installs it, exactly as jmt_find_path is used above.
    from .solvers.cpp_dispatch import _find_existing_line_cli, line_cli_platform_tag
    cli_path = _find_existing_line_cli()
    if cli_path is not None:
        print("  %s" % cli_path)
    elif line_cli_platform_tag() is None:
        warnings.warn(
            "no 'line-cli' build is published for this platform, so lang='cpp' "
            "cannot run. Build one with cpp/make.sh -O and set LINE_CLI_BINARY "
            "to it, or leave lang at its default ('python').")
        has_warnings = True
    else:
        warnings.warn(
            "the line-cli binary was not found, so lang='cpp' cannot run yet. "
            "It is about 36MB compressed (94MB once unpacked) and is downloaded "
            "automatically on the first lang='cpp' call; set LINE_CLI_BINARY to "
            "use an existing copy instead, or LINE_CLI_DOWNLOAD=0 to refuse the "
            "download.")
        has_warnings = True

    print("Checking symbolic backend (line-sage-rest)...")
    if shutil.which("docker") is None or not _run_ok(["docker", "info"]):
        warnings.warn(
            "Docker is not available, so the SageMath symbolic backend cannot "
            "start. It is required by the symbolic methods of SolverCTMC/"
            "SolverFluid (set_backend('sage')). Install Docker, then run: "
            "docker pull %s" % DOCKER_IMAGES[0])
        has_warnings = True
    elif _find_image() is None:
        warnings.warn(
            "the line-sage-rest image is not present locally, this may be "
            "required by some LINE methods. Pull it with: "
            "docker pull %s" % DOCKER_IMAGES[0])
        has_warnings = True

    if has_warnings:
        print("Completed. LINE has warnings.")
    else:
        print("Success. LINE is ready to use.")
    return not has_warnings


def main():
    line_install()


if __name__ == "__main__":
    main()
