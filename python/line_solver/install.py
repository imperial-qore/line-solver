"""
Environment check for the native Python LINE package.

Mirrors the MATLAB ``lineInstall`` script: it verifies that the optional
external dependencies used by some solvers are reachable and warns, without
failing, when one is missing. In particular it reports whether the SageMath
symbolic backend (the ``imperialqore/line-sage-rest`` Docker image, required by
the symbolic methods of SolverCTMC/SolverFluid) can be started on first use.

Run it as ``line-install`` on the command line or call ``line_install()``.
"""

import shutil
import subprocess
import warnings

from line_solver.api.sym.sage_client import DOCKER_IMAGES, _find_image


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

    print("Checking symbolic backend (line-sage-rest)...")
    if shutil.which("docker") is None or not _run_ok(["docker", "info"]):
        warnings.warn(
            "Docker is not available, so the SageMath symbolic backend cannot "
            "start. It is required by the symbolic methods of SolverCTMC/"
            "SolverFluid (set_backend('sage')). Install Docker, then run: "
            "docker run -d -p 8080:8080 %s" % DOCKER_IMAGES[0])
        has_warnings = True
    elif _find_image() is None:
        warnings.warn(
            "the line-sage-rest image is not present locally, this may be "
            "required by some LINE methods. Pull it with: "
            "docker run -d -p 8080:8080 %s" % DOCKER_IMAGES[0])
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
