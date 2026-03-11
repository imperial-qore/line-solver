"""
Lazy JVM initialization for LINE.

This module provides on-demand JVM loading so that the LINE package
can be imported without starting the JVM. The JVM is started only
when Java functionality is actually needed.

Usage:
    from line_solver._jvm import ensure_jvm, is_jvm_started

    # Check if JVM is running
    if not is_jvm_started():
        print("JVM not started yet")

    # Start JVM if not already running
    ensure_jvm()
"""

import os
import sys
from pathlib import Path
from typing import Optional

# Global state
_jvm_started = False
_jline_package = None


def get_jar_path() -> str:
    """
    Get the path to the jline.jar file.

    Returns:
        Path to jline.jar

    Raises:
        FileNotFoundError: If jline.jar cannot be found
    """
    # Try relative to this package first
    package_dir = Path(__file__).parent
    jar_candidates = [
        package_dir / '..' / '..' / 'common' / 'jline.jar',  # From python/line_solver/ to common/
        package_dir.parent.parent / 'common' / 'jline.jar',  # Explicit path
        Path(os.environ.get('JLINE_JAR', '')) if os.environ.get('JLINE_JAR') else None,
    ]

    for jar_path in jar_candidates:
        if jar_path and jar_path.exists():
            return str(jar_path.resolve())

    # Try finding it via importlib.resources for installed packages
    try:
        import importlib.resources as pkg_resources
        with pkg_resources.path('line_solver', 'jline.jar') as jar:
            if jar.exists():
                return str(jar)
    except (ImportError, FileNotFoundError, TypeError):
        pass

    raise FileNotFoundError(
        "Could not find jline.jar. Set JLINE_JAR environment variable "
        "or ensure it's in the common/ directory."
    )


def get_ldes_jar_path() -> Optional[str]:
    """
    Get the path to the ldes.jar file, downloading it if necessary.

    Returns:
        Path to ldes.jar, or None if unavailable
    """
    package_dir = Path(__file__).parent
    common_dir = (package_dir / '..' / '..' / 'common').resolve()
    ldes_jar = common_dir / 'ldes.jar'

    if ldes_jar.exists():
        return str(ldes_jar)

    # Try alternative path
    ldes_jar_alt = package_dir.parent.parent / 'common' / 'ldes.jar'
    if ldes_jar_alt.exists():
        return str(ldes_jar_alt.resolve())

    # Auto-download from SourceForge
    print("ldes.jar not found in", str(common_dir))
    print("Attempting to download ldes.jar...")
    try:
        from urllib.request import urlretrieve
        ldes_url = 'https://line-solver.sourceforge.net/latest/ldes.jar'
        common_dir.mkdir(parents=True, exist_ok=True)
        urlretrieve(ldes_url, str(ldes_jar))
        print("Successfully downloaded ldes.jar to", str(common_dir))
        return str(ldes_jar)
    except Exception as e:
        print(f"Warning: Failed to download ldes.jar: {e}")
        print("LDES solver will use the version bundled in jline.jar")
        return None


def ensure_jvm(jar_path: Optional[str] = None, jvm_args: Optional[list] = None):
    """
    Ensure the JVM is started with the LINE JAR.

    This function is idempotent - calling it multiple times has no effect
    if the JVM is already started.

    Args:
        jar_path: Optional explicit path to jline.jar
        jvm_args: Optional JVM arguments (e.g., ['-Xmx2g'])

    Returns:
        The jline package for accessing Java classes
    """
    global _jvm_started, _jline_package

    if _jvm_started:
        return _jline_package

    import jpype
    import jpype.imports

    if jar_path is None:
        jar_path = get_jar_path()

    if jvm_args is None:
        jvm_args = []

    # Add default JVM args if not specified
    if not any(arg.startswith('-Xmx') for arg in jvm_args):
        jvm_args.append('-Xmx4g')

    # Build classpath: ldes.jar first (overrides bundled LDES in jline.jar)
    classpath = []
    ldes_jar = get_ldes_jar_path()
    if ldes_jar:
        classpath.append(ldes_jar)
    classpath.append(jar_path)

    if not jpype.isJVMStarted():
        jpype.startJVM(*jvm_args, classpath=classpath)
    else:
        # JVM already started, just add our classpath
        for cp in classpath:
            jpype.addClassPath(cp)

    # Import the jline package
    _jline_package = jpype.JPackage('jline')
    _jvm_started = True

    # Set library attribution shown flag
    try:
        _jline_package.lang.GlobalConstants.setLibraryAttributionShown(True)
    except Exception:
        pass  # Ignore if not available

    return _jline_package


def is_jvm_started() -> bool:
    """
    Check if the JVM has been started.

    Returns:
        True if JVM is running, False otherwise
    """
    return _jvm_started


def get_jline_package():
    """
    Get the jline Java package.

    Returns:
        The jline package for accessing Java classes

    Raises:
        RuntimeError: If JVM has not been started
    """
    if not _jvm_started:
        raise RuntimeError("JVM not started. Call ensure_jvm() first.")
    return _jline_package


def shutdown_jvm():
    """
    Shutdown the JVM.

    Note: After calling this, the JVM cannot be restarted in the same process.
    """
    global _jvm_started, _jline_package

    if _jvm_started:
        import jpype
        try:
            jpype.shutdownJVM()
        except Exception:
            pass
        _jvm_started = False
        _jline_package = None
