# Lazy imports to avoid loading all modules at once
# Only import essential modules that have been updated to use local butools
try:
    from .check import CheckPHRepresentation
except ImportError:
    CheckPHRepresentation = None

# Optional imports - only load if external butools is available
try:
    from .baseph import *
    from .misc import *
    from .monocyclic import *
    from .orders import *
    from .canonical import *
    from .appie import *
    from .phfromme import *
except ImportError:
    pass
