"""
JMT viewer functions.

This module provides functions to open JMT model files in JMT GUI applications:
- jsimg_view: Open model in JSIMgraph (graphical editor)
- jsimw_view: Open model in JSIMwiz (wizard interface)

Port from:
    - matlab/src/io/jsimgView.m
    - matlab/src/io/jsimwView.m

Note: These functions require JMT.jar to be installed and Java to be available
on the system PATH.
"""

import os
import subprocess
import platform
from typing import Optional

from .logging import line_warning, line_error, get_logger, VerboseLevel


def get_jmt_path() -> Optional[str]:
    """
    Get the path to JMT.jar.

    Searches for JMT.jar in common locations:
    - LINE common/ directory
    - JMT installation directory
    - Current directory

    Returns:
        Path to directory containing JMT.jar, or None if not found
    """
    # Common locations to check
    search_paths = []

    # Check LINE common directory (relative to this module)
    module_dir = os.path.dirname(os.path.abspath(__file__))
    line_root = os.path.abspath(os.path.join(module_dir, '..', '..', '..', '..'))
    search_paths.append(os.path.join(line_root, 'common'))

    # Check JMT_HOME environment variable
    jmt_home = os.environ.get('JMT_HOME')
    if jmt_home:
        search_paths.append(jmt_home)

    # Check common installation paths
    if platform.system() == 'Windows':
        search_paths.extend([
            os.path.join(os.environ.get('PROGRAMFILES', 'C:\\Program Files'), 'JMT'),
            os.path.join(os.environ.get('LOCALAPPDATA', ''), 'JMT'),
        ])
    elif platform.system() == 'Darwin':  # macOS
        search_paths.extend([
            '/Applications/JMT',
            os.path.expanduser('~/Applications/JMT'),
        ])
    else:  # Linux
        search_paths.extend([
            '/usr/share/jmt',
            '/opt/jmt',
            os.path.expanduser('~/.jmt'),
        ])

    # Current directory
    search_paths.append(os.getcwd())

    # Search for JMT.jar
    for path in search_paths:
        if path and os.path.exists(path):
            jmt_jar = os.path.join(path, 'JMT.jar')
            if os.path.exists(jmt_jar):
                return path

    return None


def jsimg_view(filename: str, jmt_path: Optional[str] = None,
               suppress_output: bool = True) -> bool:
    """
    Open a JSIM model file in JSIMgraph (JMT graphical editor).

    JSIMgraph provides a graphical interface for viewing and editing
    queueing network models in JMT format.

    Args:
        filename: Path to the JSIM/JSIMG file to open
        jmt_path: Optional path to JMT.jar directory. If None, searches
                  common locations.
        suppress_output: Whether to suppress JMT output (default True)

    Returns:
        True if JMT was launched successfully, False otherwise

    Example:
        >>> # Export model to JSIMG and view it
        >>> from line_solver.api.io import qn2jsimg, jsimg_view
        >>> jsimg_file = qn2jsimg(model)
        >>> jsimg_view(jsimg_file)

    References:
        MATLAB: matlab/src/io/jsimgView.m
    """
    return _launch_jmt_viewer(filename, 'jsimg', jmt_path, suppress_output)


def jsimw_view(filename: str, jmt_path: Optional[str] = None,
               suppress_output: bool = True) -> bool:
    """
    Open a JSIM model file in JSIMwiz (JMT wizard interface).

    JSIMwiz provides a wizard-style interface for configuring and running
    simulations of queueing network models.

    Args:
        filename: Path to the JSIM/JSIMG file to open
        jmt_path: Optional path to JMT.jar directory. If None, searches
                  common locations.
        suppress_output: Whether to suppress JMT output (default True)

    Returns:
        True if JMT was launched successfully, False otherwise

    Example:
        >>> # Export model to JSIMG and view it in wizard
        >>> from line_solver.api.io import qn2jsimg, jsimw_view
        >>> jsimg_file = qn2jsimg(model)
        >>> jsimw_view(jsimg_file)

    References:
        MATLAB: matlab/src/io/jsimwView.m
    """
    return _launch_jmt_viewer(filename, 'jsimw', jmt_path, suppress_output)


def jmva_view(filename: str, jmt_path: Optional[str] = None,
              suppress_output: bool = True) -> bool:
    """
    Open a JMVA model file in JMVA (JMT MVA solver).

    JMVA provides an interface for Mean Value Analysis of queueing networks.

    Args:
        filename: Path to the JMVA file to open
        jmt_path: Optional path to JMT.jar directory. If None, searches
                  common locations.
        suppress_output: Whether to suppress JMT output (default True)

    Returns:
        True if JMT was launched successfully, False otherwise
    """
    return _launch_jmt_viewer(filename, 'jmva', jmt_path, suppress_output)


def _launch_jmt_viewer(filename: str, viewer_type: str,
                       jmt_path: Optional[str] = None,
                       suppress_output: bool = True) -> bool:
    """
    Launch a JMT viewer for the specified file.

    Args:
        filename: Path to the model file
        viewer_type: Type of viewer ('jsimg', 'jsimw', or 'jmva')
        jmt_path: Path to JMT.jar directory
        suppress_output: Whether to suppress JMT output

    Returns:
        True if launched successfully, False otherwise
    """
    logger = get_logger()

    # Resolve absolute path for filename
    if not os.path.isabs(filename):
        filename = os.path.abspath(filename)

    # Check if file exists
    if not os.path.exists(filename):
        line_error('jmt_viewer', f'File not found: {filename}')
        return False

    # Find JMT.jar
    if jmt_path is None:
        jmt_path = get_jmt_path()

    if jmt_path is None:
        line_error('jmt_viewer', 'JMT.jar not found. Please set JMT_HOME environment variable or provide jmt_path.')
        return False

    jmt_jar = os.path.join(jmt_path, 'JMT.jar')
    if not os.path.exists(jmt_jar):
        line_error('jmt_viewer', f'JMT.jar not found at: {jmt_jar}')
        return False

    # Build command
    java_cmd = 'java'
    classpath = jmt_jar
    main_class = 'jmt.commandline.Jmt'

    # Check verbosity level
    if hasattr(logger, 'level'):
        suppress_output = logger.level != VerboseLevel.DEBUG

    # Build the command list
    cmd = [java_cmd, '-cp', classpath, main_class, viewer_type, filename]

    # Determine output redirection
    if suppress_output:
        if platform.system() == 'Windows':
            stdout = subprocess.DEVNULL
            stderr = subprocess.DEVNULL
        else:
            stdout = subprocess.DEVNULL
            stderr = subprocess.DEVNULL
    else:
        stdout = None
        stderr = None

    try:
        # Try to launch JMT
        process = subprocess.Popen(
            cmd,
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,  # Detach from parent process
        )

        # Check if process started (give it a moment)
        try:
            # Poll to see if it immediately failed
            retcode = process.poll()
            if retcode is not None and retcode != 0:
                # Process failed immediately, try with --illegal-access=permit
                cmd_retry = [java_cmd, '--illegal-access=permit', '-cp', classpath,
                           main_class, viewer_type, filename]
                process = subprocess.Popen(
                    cmd_retry,
                    stdout=stdout,
                    stderr=stderr,
                    start_new_session=True,
                )
                retcode = process.poll()
                if retcode is not None and retcode != 0:
                    line_warning('jmt_viewer', f'JMT process exited with code {retcode}')
                    return False
        except Exception:
            pass  # Process is still running, which is good

        return True

    except FileNotFoundError:
        line_error('jmt_viewer', 'Java not found. Please ensure Java is installed and on the PATH.')
        return False
    except Exception as e:
        line_error('jmt_viewer', f'Failed to launch JMT: {e}')
        return False


__all__ = [
    'get_jmt_path',
    'jsimg_view',
    'jsimw_view',
    'jmva_view',
]
