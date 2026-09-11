#!/usr/bin/env python3
"""Launcher for the LINE Solver MCP server, kept at its historical path.

The server itself lives in the package, at ``line_solver/mcp_server.py``, so
that ``pip install line-solver`` ships it and the ``line-mcp`` console script
can start it. This file remains so that existing client configurations naming
``python/mcp_server.py`` keep working, and so that the server can be started
from a source checkout in which the package has not been installed: the
checkout is put ahead of any installed copy, as it was before the move.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from line_solver.mcp_server import main

if __name__ == "__main__":
    main()
