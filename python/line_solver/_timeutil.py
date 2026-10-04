"""Wall-clock time-budget helpers for native solver cooperative checkpoints.

A solver captures ``t0 = time.time()`` at launch and calls ``timeout_exceeded``
at loop boundaries. When the elapsed wall-clock time exceeds ``timeout`` (in
seconds) the solver stops and returns an interim solution if available, else an
empty result with a warning. A non-finite or non-positive ``timeout`` means no
budget.
"""

import math
import time


def timeout_exceeded(t0, timeout):
    """Return True if ``time.time() - t0`` exceeds ``timeout`` seconds."""
    if timeout is None:
        return False
    try:
        budget = float(timeout)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(budget) or budget <= 0:
        return False
    return (time.time() - t0) > budget
