"""
Wrapper for rmf_tool - Refined Mean Field Approximation library.

Vendored from https://github.com/ngast/rmf_tool (MIT License)
Author: Nicolas Gast (nicolas.gast@inria.fr)

Provides:
- Mean field approximation for population processes
- Refined mean field approximation (1/N corrections)
- Transient and steady-state analysis
"""
import sys
import os

# Add the rmf_tool/src directory to the path so rmftool package is importable
_rmf_src = os.path.join(os.path.dirname(__file__), 'rmf_tool', 'src')
if _rmf_src not in sys.path:
    sys.path.insert(0, _rmf_src)

from rmftool import *
