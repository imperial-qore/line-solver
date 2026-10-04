"""Filesystem helpers shared by the subprocess-backed solvers.

Port from:
    - matlab/lineTempName.m
    - jar/src/main/java/jline/io/SysUtils.java (lineTempName)
"""

import os
import tempfile


def line_temp_name(solvername: str = '', dockermountable: bool = False) -> str:
    """Create and return a private staging directory for one solver run.

    Mirrors MATLAB ``lineTempName`` and the JAR's ``SysUtils.lineTempName``:
    the directory is ``<root>/line_workspace/<solvername>/tmp_<random>``, and
    the root is, in order of precedence,

      1. ``$LINE_WORKSPACE_ROOT`` when set. run-tests.sh sets it when a solver is
         wrapped in a container, so the staging path is one the daemon can
         bind-mount; no solver needs to know a container is involved.
      2. ``$HOME/.line`` when ``dockermountable``, because snap-confined Docker
         cannot bind-mount the system temp dir.
      3. the system temp dir.

    Honouring ``$LINE_WORKSPACE_ROOT`` is not cosmetic: a solver that stages
    under /tmp regardless hands the container a path it cannot see, and the run
    fails with the model file "missing" rather than with anything naming the
    cause.
    """
    baseroot = tempfile.gettempdir()
    if dockermountable:
        home = os.path.expanduser('~')
        if home:
            baseroot = os.path.join(home, '.line')

    envroot = os.environ.get('LINE_WORKSPACE_ROOT', '').strip()
    if envroot:
        baseroot = envroot

    basedir = os.path.join(baseroot, 'line_workspace', solvername or '')
    os.makedirs(basedir, exist_ok=True)
    return tempfile.mkdtemp(prefix='tmp_', dir=basedir)
