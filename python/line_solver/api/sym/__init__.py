"""
Computer algebra backend for the native Python implementation.

sympy covers the symbolic work by default. This package adds SageMath behind
the line-sage-rest service as an alternative engine, shared with the JAR
(``jline.api.sym``) and MATLAB (``SAGE.m``): it is exact and fast over the
rational function field, and it gives the three codebases one normal form, so
a symbolic result computed in one is comparable with another.

It is opt in. Set ``LINE_SYMBOLIC_BACKEND=sage``, or call
``set_backend("http://host:port")``, and the symbolic branch of ``ctmc_solve``
routes through the service; otherwise nothing changes.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from .sage_client import (DOCKER_IMAGES, MODE_ENV, URL_ENV, SageError,
                          SageRestEngine, backend_mode, prefers_sage, require,
                          resolve, set_backend, stop_container)

__all__ = ["SageRestEngine", "SageError", "resolve", "require", "set_backend",
           "backend_mode", "prefers_sage", "stop_container", "DOCKER_IMAGES",
           "URL_ENV", "MODE_ENV"]
