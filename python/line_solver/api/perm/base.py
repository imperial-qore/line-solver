"""
Common base for the matrix permanent solvers.

Port of jline.lib.perm.PermSolver: every permanent algorithm, exact or
approximate, exposes the same triple (value, time, memory) so that the
algorithms can be benchmarked against each other on the same matrix.
"""

import time
import tracemalloc
from dataclasses import dataclass

import numpy as np


@dataclass
class PermResult:
    """Result of permanent computation."""
    value: float
    time_ms: float
    n: int


class PermSolver:
    """
    Abstract base class for permanent computation solvers.

    Concrete solvers implement compute() and store the outcome in self.value.
    solve() wraps compute() with wall-clock and peak-allocation measurement,
    mirroring the JAR PermSolver.solve() instrumentation.

    Attributes:
        matrix: Input square matrix as a float ndarray
        n: Matrix order
        value: Permanent (or its approximation) after compute()
        time_ms: Wall-clock duration of the last solve() in milliseconds
        memory_bytes: Peak allocation of the last solve() in bytes
    """

    def __init__(self, matrix: np.ndarray):
        self.matrix = np.asarray(matrix, dtype=float)
        self.n = self.matrix.shape[0]
        self.value = 0.0
        self.time_ms = 0.0
        self.memory_bytes = 0

    def get_matrix(self) -> np.ndarray:
        """Return the input matrix."""
        return self.matrix

    def get_n(self) -> int:
        """Return the matrix order."""
        return self.n

    def get_value(self) -> float:
        """Return the computed permanent value."""
        return self.value

    def get_time(self) -> float:
        """Return the duration of the last solve() in milliseconds."""
        return self.time_ms

    def get_memory(self) -> int:
        """Return the peak allocation of the last solve() in bytes."""
        return self.memory_bytes

    def compute(self):
        """Compute the permanent or its approximation; sets self.value."""
        raise NotImplementedError

    def solve(self):
        """Run compute() while measuring time and peak memory."""
        tracing = tracemalloc.is_tracing()
        if not tracing:
            tracemalloc.start()
        baseline = tracemalloc.get_traced_memory()[0]
        start_time = time.perf_counter()

        self.compute()

        self.time_ms = (time.perf_counter() - start_time) * 1000.0
        self.memory_bytes = max(tracemalloc.get_traced_memory()[1] - baseline, 0)
        if not tracing:
            tracemalloc.stop()

    def get_result(self) -> PermResult:
        """Get the computation result."""
        return PermResult(value=self.value, time_ms=self.time_ms, n=self.n)


def require_full_support(matrix, caller: str) -> None:
    """
    Refuse a matrix the permanent approximations cannot take.

    The four approximate engines -- BethePermanent, HeuristicPermanent,
    HuberLawSampler and AdaPartSampler -- all rest, directly or through the
    Sinkhorn scaling they share, on a strictly positive matrix.

    Why this is refused rather than floored. Each of these used to replace a
    zero by a small constant eps before doing anything else, and that
    substitution is not invertible: every permutation picks exactly one entry
    from each row, so a matrix with an identically zero row has permanent 0
    while the floored matrix has permanent n! * eps times the permanent of the
    rest. n! outruns eps quickly -- with eps = 2.22e-16 the fabricated value
    passes 1% at n = 17 and 1 at n = 18, and at n = 20 the floor alone
    manufactures a permanent of about 540 where the truth is exactly zero. The
    order of the replicated demand matrix in pfqn_jointmarg is sum(N), so a
    closed model with 18 jobs and one structurally zero demand is already
    there.

    Positivity is sufficient but not necessary: the sharp precondition of the
    Sinkhorn scaling is TOTAL SUPPORT (every positive entry lies on a positive
    permutation), which a matrix with a strictly positive permanent can still
    fail -- [[J3, 0], [J3, J3]] has permanent 36, no total support, and the
    scaling exits on its tolerance rather than on convergence. Positivity is
    used here because it is O(n^2) and is the contract the C++ header states.

    Args:
        matrix: Matrix to check
        caller: Name of the engine, used in the message

    Raises:
        ValueError: If any entry is not strictly positive
    """
    m = np.asarray(matrix, dtype=float)
    bad = np.argwhere(~(m > 0.0))
    if bad.size > 0:
        i, j = int(bad[0][0]), int(bad[0][1])
        raise ValueError(
            "The '" + caller + "' permanent approximation requires a strictly "
            "positive matrix: entry (" + str(i) + ", " + str(j) + ") is "
            + str(m[i, j]) + ", so the matrix has no full support. Flooring it "
            "would change the permanent by n!*eps, which is O(1) by n=18. Use "
            "the exact engine (permanent / perm).")
