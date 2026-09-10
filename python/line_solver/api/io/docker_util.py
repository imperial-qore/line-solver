"""
Docker primitives shared by the backends that legitimately ship an image: the
JMT backend (``imperialqore/jmt-rest``) and the Sage symbolic engine
(``imperialqore/line-sage-rest``). Each owns its image name; this module only
answers whether the daemon is up, whether an image is present, whether there is
room to pull it, and performs the pull. Mirrors ``jline.io.DockerImage``.

LQNS, lqsim and qnsolver are deliberately absent: their licence is an evaluation
agreement that forbids redistribution, so LINE runs them only from a binary the
user installed. ``run-tests.sh --lqns-docker`` puts a shim on the PATH when a
containerised build is what should be exercised.

Every function is a no-op / returns empty on Windows.
"""

import os
import platform
import re
import subprocess
import sys

# Conservative free-space floor required before a pull (2 GiB).
_DEFAULT_MIN_FREE_BYTES = 2 * 1024 * 1024 * 1024


def _is_windows():
    return platform.system() == 'Windows'


def daemon_available():
    """True if the Docker daemon is reachable (unix only)."""
    if _is_windows():
        return False
    try:
        return subprocess.run(['docker', 'info'], stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL, timeout=15).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


_command_cache = {}


def has_local_image(image):
    """True if the named image is already present in the local Docker store."""
    if _is_windows() or not image:
        return False
    try:
        out = subprocess.run(['docker', 'images', '-q', image],
                             capture_output=True, timeout=15)
        return out.returncode == 0 and out.stdout.decode().strip() != ''
    except (OSError, subprocess.SubprocessError):
        return False


def pull(image):
    """Pull an image, streaming Docker's progress to the console. No timeout."""
    if _is_windows() or not image:
        return False
    try:
        return subprocess.run(['docker', 'pull', image]).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def has_storage_for(image):
    """True if the Docker storage location has room for ``image``."""
    free = _free_bytes_at_docker_root()
    if free < 0:
        return True  # could not determine; do not block the pull
    return free >= _required_bytes(image)


def _required_bytes(image):
    override = _override_min_free_bytes()
    if override > 0:
        return override
    return max(_DEFAULT_MIN_FREE_BYTES, _estimate_image_bytes(image))


def _override_min_free_bytes():
    v = os.environ.get('LINE_DOCKER_MIN_FREE_BYTES')
    if v:
        try:
            return int(v.strip())
        except ValueError:
            pass
    return -1


def _free_bytes_at_docker_root():
    """Usable bytes on the filesystem backing the Docker root dir; -1 if unknown."""
    root = None
    try:
        out = subprocess.run(['docker', 'info', '--format', '{{.DockerRootDir}}'],
                             capture_output=True, timeout=15)
        if out.returncode == 0:
            root = out.stdout.decode().strip()
    except (OSError, subprocess.SubprocessError):
        root = None
    if not root:
        root = '/var/lib/docker'
    # The root dir may not exist for / be readable by this user: walk up to an
    # existing ancestor so statvfs reports a real filesystem.
    while root and not os.path.exists(root):
        parent = os.path.dirname(root)
        if parent == root:
            break
        root = parent
    if not root or not os.path.exists(root):
        root = '/'
    try:
        st = os.statvfs(root)
        free = st.f_bavail * st.f_frsize
        return free if free > 0 else -1
    except OSError:
        return -1


def _estimate_image_bytes(image):
    """On-disk estimate (compressed layer sizes x3), or 0 if it cannot be determined."""
    try:
        out = subprocess.run(['docker', 'manifest', 'inspect', image],
                             capture_output=True, timeout=30)
        if out.returncode != 0:
            return 0
        json_txt = out.stdout.decode()
    except (OSError, subprocess.SubprocessError):
        return 0
    total = 0
    for m in re.findall(r'"size"\s*:\s*(\d+)', json_txt):
        try:
            total += int(m)
        except ValueError:
            pass
    return total * 3 if total > 0 else 0


def _capture(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, timeout=10)
        if out.returncode == 0:
            return out.stdout.decode().strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return ''
