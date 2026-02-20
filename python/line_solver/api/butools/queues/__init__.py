# Lazy imports to avoid loading all modules at once
# Only import MMAPPH1FCFS which has been updated to use local butools
try:
    from .mtfcfs import MMAPPH1FCFS
except ImportError:
    MMAPPH1FCFS = None

# Optional imports - only load if external butools is available
try:
    from .dfcfs import *
    from .cfcfs import *
    from .prprio import *
    from .npprio import *
    from .flufluqueue import *
except ImportError:
    pass
