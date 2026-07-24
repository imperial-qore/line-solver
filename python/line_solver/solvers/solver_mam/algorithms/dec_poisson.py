"""
Decomposition with Poisson arrivals (dec.poisson) algorithm.

Simplified variant of dec.source that uses Poisson (single-state exponential)
approximations for all arrival processes by setting space_max=1.

This is faster than dec.source but less accurate, as it eliminates the
correlated structure in multi-class arrival processes.

Algorithm:
Same as dec.source but with MMAP order limited to 1 (exponential only).
"""

from typing import Tuple, Optional

from . import MAMAlgorithm, MAMResult
from .dec_source import DecSourceAlgorithm, DecSourceOptions


class DecPoissonAlgorithm(MAMAlgorithm):
    """Decomposition with Poisson arrivals (simplified, faster version)."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if network can be solved by dec.poisson.

        dec.poisson supports same networks as dec.source.
        It's a simplified variant, not a restricted variant.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        return DecSourceAlgorithm.supports_network(sn)

    def solve(self, sn, options=None) -> MAMResult:
        """Solve using dec.source with Poisson (space_max=1) approximation.

        This is a wrapper around dec.source with space_max=1 to force
        all arrival processes to be single-state exponential (Poisson).

        Args:
            sn: NetworkStruct
            options: Solver options (or DecSourceOptions)

        Returns:
            MAMResult
        """
        # Create options with space_max=1 for Poisson approximation
        if options is None:
            opts = DecSourceOptions(space_max=1)
        else:
            # Make a copy and override space_max
            if hasattr(options, 'space_max'):
                options.space_max = 1
            else:
                opts = DecSourceOptions(
                    tol=getattr(options, 'tol', 1e-6),
                    max_iter=getattr(options, 'max_iter', 100),
                    space_max=1,
                    verbose=getattr(options, 'verbose', False)
                )
                options = opts

        # Use dec.source solver with space_max=1
        solver = DecSourceAlgorithm()
        result = solver.solve(sn, options)

        # Override method name
        result.method = "dec.poisson"

        return result


from .dec_source import DecSourceOptions
