# Lazy imports to avoid loading all modules at once
# Only import essential modules that have been updated to use local butools
try:
    from .check import CheckMAPRepresentation, CheckMMAPRepresentation
except ImportError:
    CheckMAPRepresentation = None
    CheckMMAPRepresentation = None

# Optional imports - only load if external butools is available
try:
    from .basemap import *
    from .misc import *
    from .matching import *
    from .minimalrep import *
    from .mapfromrap import *
except ImportError:
    pass
