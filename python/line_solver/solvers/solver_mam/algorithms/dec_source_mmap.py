"""
MMAP fork-join decomposition (dec.source.mmap) algorithm.

Refines the per-node departure processes of the network by a parametric
decomposition, superposing and splitting them as marked MAPs and synchronizing
the branches of a fork-join construct with mmap_max. Unlike dec.source, the
arrival process seen by each queue keeps its marked multi-class structure, and
the join contributes a synchronization delay derived from the parallel branch
response times.

References:
    MATLAB: matlab/src/solvers/MAM/solver_mam_basic_mmap.m and the
            solver_mam_basic_mmap_inner / _closed / solver_mam_traffic_mmap
            family it dispatches to.
"""

import time
from typing import Optional, Tuple

import numpy as np

from . import MAMAlgorithm, MAMResult
from .dec_source import DecSourceAlgorithm
from ....api.solvers.mam.handler import SolverMAMOptions as HandlerOptions
from ....api.solvers.mam.mmap_fj import solver_mam_basic_mmap


class DecSourceMMAPAlgorithm(MAMAlgorithm):
    """MMAP fork-join decomposition with mmap_max join synchronization."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if the network can be solved by dec.source.mmap.

        The method is the fork-join-aware generalization of dec.source and
        accepts the same networks; a fork-join topology is what it is for, not a
        restriction.
        """
        return DecSourceAlgorithm.supports_network(sn)

    def solve(self, sn, options=None) -> MAMResult:
        """Solve the network with the MMAP fork-join decomposition.

        Args:
            sn: NetworkStruct
            options: SolverMAMOptions (solver level) or None

        Returns:
            MAMResult with the QN, UN, RN, TN metrics
        """
        start_time = time.time()

        # see _kb/06-solver-catalog.md (MAM: "dec.source.mmap config defaults
        # mirror MATLAB") for why these handler defaults are not overridden
        handler_opts = HandlerOptions()
        if options is not None:
            if hasattr(options, 'tol'):
                # MATLAB drives the departure-process fixed point on iter_tol,
                # which SolverMAM.defaultOptions sets equal to tol.
                handler_opts.tol = options.tol
                handler_opts.iter_tol = options.tol
            if hasattr(options, 'max_iter'):
                handler_opts.iter_max = options.max_iter
            if hasattr(options, 'verbose'):
                handler_opts.verbose = options.verbose

        ret = solver_mam_basic_mmap(sn, handler_opts)

        TN = ret.T if ret.T is not None else np.zeros_like(ret.Q)

        return MAMResult(
            QN=ret.Q,
            UN=ret.U,
            RN=ret.R,
            TN=TN,
            CN=ret.C,
            XN=ret.X,
            totiter=ret.it,
            method="dec.source.mmap",
            runtime=time.time() - start_time,
        )
