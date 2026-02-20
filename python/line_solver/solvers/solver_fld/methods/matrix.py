"""
Matrix method wrapper for FLD solver.

Wraps the existing api/solvers/fld/handler.py implementation.
"""

import numpy as np
from typing import Optional
from ...sn import NetworkStruct
from ..options import SolverFLDOptions, FLDResult
from ..utils import extract_metrics_from_handler_result


class MatrixMethod:
    """Matrix method for fluid analysis using p-norm smoothing.

    Default method that uses the Ruuskanen et al. 2021 matrix formulation
    with p-norm smoothing for capacity constraints.

    Wraps: api/solvers/fld/handler.py
    """

    def __init__(self, sn: NetworkStruct, options: SolverFLDOptions):
        """Initialize matrix method.

        Args:
            sn: Compiled NetworkStruct
            options: SolverFLDOptions
        """
        self.sn = sn
        self.options = options

    def solve(self) -> FLDResult:
        """Solve using matrix method.

        Returns:
            FLDResult with performance metrics
        """
        from ....api.solvers.fld.handler import solver_fld, SolverFLDOptions as HandlerFLDOptions

        # Convert options to handler format
        handler_opts = HandlerFLDOptions(
            method=self.options.method if self.options.method != 'default' else 'matrix',
            tol=self.options.tol,
            verbose=self.options.verbose,
            stiff=self.options.stiff,
            iter_max=self.options.iter_max,
            timespan=self.options.timespan,
            pstar=[self.options.pstar] * self.sn.nstations,
            num_cdf_pts=200
        )

        # Solve using handler
        handler_result = solver_fld(self.sn, handler_opts)

        # Extract metrics
        QN, UN, RN, TN, CN, XN = extract_metrics_from_handler_result(handler_result, self.sn)

        # Extract transient data from handler result
        QNt = {}
        UNt = {}
        TNt = {}
        if hasattr(handler_result, 'Qt') and handler_result.Qt is not None:
            M = len(handler_result.Qt)
            K = len(handler_result.Qt[0]) if M > 0 else 0
            for i in range(M):
                for r in range(K):
                    if handler_result.Qt[i][r] is not None:
                        QNt[(i, r)] = np.array(handler_result.Qt[i][r])
                    if hasattr(handler_result, 'Ut') and handler_result.Ut is not None:
                        UNt[(i, r)] = np.array(handler_result.Ut[i][r])
                    if hasattr(handler_result, 'Tt') and handler_result.Tt is not None:
                        TNt[(i, r)] = np.array(handler_result.Tt[i][r])

        # Build result
        result = FLDResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            CN=CN,
            XN=XN,
            t=handler_result.t,
            QNt=QNt,
            UNt=UNt,
            TNt=TNt,
            xvec=handler_result.odeStateVec,
            iterations=handler_result.it,
            runtime=0.0,  # Runtime set by caller
            method='matrix'
        )

        return result


def solve_matrix(sn: NetworkStruct, options: Optional[SolverFLDOptions] = None) -> FLDResult:
    """Convenience function to solve using matrix method.

    Args:
        sn: Compiled NetworkStruct
        options: SolverFLDOptions (uses defaults if None)

    Returns:
        FLDResult
    """
    if options is None:
        options = SolverFLDOptions(method='matrix')

    method = MatrixMethod(sn, options)
    return method.solve()
