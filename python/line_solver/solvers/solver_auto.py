"""
Native Python implementation of automatic solver selection.

This module provides the SolverAuto class that analyzes the network model
and automatically selects the most appropriate solver based on model characteristics.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Type
from dataclasses import dataclass
from enum import Enum

from ..api.sn import NetworkStruct, SchedStrategy, NodeType


class SolverType(Enum):
    """Available solver types."""
    MVA = 'MVA'
    DES = 'DES'
    NC = 'NC'
    CTMC = 'CTMC'
    SSA = 'SSA'
    FLUID = 'Fluid'
    MAM = 'MAM'
    JMT = 'JMT'


@dataclass
class ModelAnalyzer:
    """
    Analyzes queueing network model characteristics for solver selection.

    This class examines the network structure and determines which solvers
    are applicable and which is likely to be most efficient.

    Original implementation by LINE
    """
    sn: NetworkStruct

    def has_single_chain(self) -> bool:
        """Check if model has exactly one chain."""
        return self.sn.nchains == 1

    def has_multi_chain(self) -> bool:
        """Check if model has multiple chains."""
        return self.sn.nchains > 1

    def is_closed_model(self) -> bool:
        """Check if all classes are closed (finite population)."""
        if self.sn.njobs is None or len(self.sn.njobs) == 0:
            return False
        return np.all(np.isfinite(self.sn.njobs.flatten()))

    def is_open_model(self) -> bool:
        """Check if all classes are open (infinite population)."""
        if self.sn.njobs is None or len(self.sn.njobs) == 0:
            return True
        return np.all(np.isinf(self.sn.njobs.flatten()))

    def is_mixed_model(self) -> bool:
        """Check if model has both open and closed classes."""
        if self.sn.njobs is None or len(self.sn.njobs) == 0:
            return False
        njobs_flat = self.sn.njobs.flatten()
        has_open = np.any(np.isinf(njobs_flat))
        has_closed = np.any(np.isfinite(njobs_flat) & (njobs_flat > 0))
        return has_open and has_closed

    def has_product_form(self) -> bool:
        """
        Check if model has product form solution.

        Product form requires:
        - FCFS, PS, LCFS-PR, or INF scheduling
        - Exponential or phase-type service
        - No finite capacities (or proper handling)
        """
        # Check scheduling strategies
        product_form_scheds = {
            SchedStrategy.FCFS,
            SchedStrategy.PS,
            SchedStrategy.LCFSPR,
            SchedStrategy.INF,
        }

        for station_id, sched in self.sn.sched.items():
            if sched not in product_form_scheds:
                return False

        # Check for finite capacity constraints
        if self.sn.cap is not None and len(self.sn.cap) > 0:
            if np.any(np.isfinite(self.sn.cap) & (self.sn.cap > 0) & (self.sn.cap < 1e6)):
                return False

        return True

    def has_multi_server(self) -> bool:
        """Check if any station has multiple servers."""
        return self.sn.has_multi_server()

    def has_load_dependence(self) -> bool:
        """Check if model has load-dependent service rates."""
        return self.sn.has_load_dependence()

    def has_homogeneous_scheduling(self, strategy: SchedStrategy) -> bool:
        """Check if all stations use the same scheduling strategy."""
        for station_id, sched in self.sn.sched.items():
            if sched != strategy:
                return False
        return True

    def has_all_infinite_servers(self) -> bool:
        """Check if all stations are infinite servers (delays)."""
        return self.has_homogeneous_scheduling(SchedStrategy.INF)

    def get_total_jobs(self) -> float:
        """Get total number of jobs in closed classes."""
        return self.sn.get_total_population()

    def get_avg_jobs_per_chain(self) -> float:
        """Get average number of jobs per chain."""
        if self.sn.nchains == 0:
            return 0.0
        return self.get_total_jobs() / self.sn.nchains

    def has_fork_join(self) -> bool:
        """Check if model has fork-join nodes."""
        if self.sn.nodetype is None:
            return False
        for nt in self.sn.nodetype:
            if nt == NodeType.FORK or nt == NodeType.JOIN:
                return True
        return False

    def has_class_switching(self) -> bool:
        """Check if model has class switching nodes."""
        if self.sn.nodetype is None:
            return False
        for nt in self.sn.nodetype:
            if nt == NodeType.CLASSSWITCH:
                return True
        return False

    def has_cache(self) -> bool:
        """Check if model has cache nodes."""
        if self.sn.nodetype is None:
            return False
        for nt in self.sn.nodetype:
            if nt == NodeType.CACHE:
                return True
        return False

    def has_finite_capacity(self) -> bool:
        """Check if model has finite capacity constraints."""
        if self.sn.cap is not None and len(self.sn.cap) > 0:
            return np.any(np.isfinite(self.sn.cap) & (self.sn.cap > 0) & (self.sn.cap < 1e6))
        return False

    def get_num_stations(self) -> int:
        """Get number of service stations."""
        return self.sn.nstations

    def get_num_classes(self) -> int:
        """Get number of job classes."""
        return self.sn.nclasses

    def get_state_space_size_estimate(self) -> float:
        """
        Estimate the size of the state space.

        Uses the formula: product of (n_i + K - 1) choose (K - 1) over all classes
        where n_i is the population of class i and K is number of stations.
        """
        if not self.is_closed_model():
            return float('inf')

        M = self.sn.nstations
        N = self.get_total_jobs()

        if N == 0 or M == 0:
            return 1.0

        # Simple approximation: (N + M - 1)! / (N! * (M-1)!)
        # Using Stirling's approximation for large values
        if N > 100 or M > 20:
            # Very rough upper bound
            return (N + M) ** min(N, M)

        # For small values, compute more accurately
        from math import comb
        return comb(int(N) + M - 1, M - 1)


@dataclass
class SolverAutoOptions:
    """Options for the automatic solver.

    Attributes:
        selection_method: Solver selection strategy:
            - 'default', 'heur': Heuristic-based automatic selection
            - 'sim': Best simulator (DES/JMT)
            - 'exact': Best exact method (NC for closed, CTMC otherwise)
            - 'fast': Fast approximate method (MVA)
            - 'accurate': Accurate approximate method (Fluid)
            # - 'ai': AI-based selection (not yet available)
        forced_solver: Force a specific solver by name
        verbose: Enable verbose output
    """
    selection_method: str = 'default'
    forced_solver: Optional[str] = None
    verbose: bool = False


class SolverAuto:
    """
    Native Python automatic solver selection.

    This solver analyzes the network model and automatically selects the most
    appropriate solver based on model characteristics, feature support, and
    performance considerations.

    Selection Algorithm

    For closed networks:
    - Single chain: Use NC (fast exact solution)
    - Multi-chain, all infinite servers: Use MVA
    - Multi-chain, product form, no multi-server: Use NC
    - FCFS without multi-server:
        - avg_jobs_per_chain > 30: Use Fluid
        - 10 < avg_jobs_per_chain <= 30: Use MVA
        - total_jobs < 5: Use NC
    - PS/PSPRIO with multi-server: Use MVA
    - FCFS with multi-server:
        - total_jobs < 5: Use NC
        - total_jobs >= 5: Use MVA

    For open networks:
    - Simple open: Use MVA or DES

    Args:
        model: Network model (Python wrapper or native structure)
        options: SolverAutoOptions configuration
        **kwargs: Additional solver options passed to selected solver

    Example:
        >>> solver = SolverAuto(model)
        >>> solver.runAnalyzer()
        >>> table = solver.getAvgTable()
        >>> print(f"Selected: {solver.get_selected_solver_name()}")
    """

    def __init__(
        self,
        model: Any,
        options: Optional[SolverAutoOptions] = None,
        **kwargs
    ):
        self.model = model
        self.options = options or SolverAutoOptions()
        self.kwargs = kwargs

        self._result = None
        self._table_silent = False
        self._selected_solver = None
        self._selected_solver_name = None
        self._candidate_solvers: List[str] = []

        # Extract network structure
        self._sn = self._get_network_struct()

        # Create model analyzer
        self._analyzer = ModelAnalyzer(sn=self._sn)

        # Determine candidate solvers
        self._determine_candidates()

    def _get_network_struct(self) -> NetworkStruct:
        """Get NetworkStruct from model."""
        model = self.model

        # If already a NetworkStruct
        if isinstance(model, NetworkStruct):
            return model

        # Handle JPype wrapper
        if hasattr(model, 'obj'):
            from .convert import wrapper_sn_to_native
            try:
                return wrapper_sn_to_native(model)
            except Exception:
                pass

        # Handle model with getStruct
        if hasattr(model, 'getStruct'):
            from .convert import wrapper_sn_to_native
            return wrapper_sn_to_native(model.getStruct())

        raise ValueError("Cannot extract network structure from model")

    def _determine_candidates(self) -> None:
        """Determine which solvers can handle this model."""
        self._candidate_solvers = []

        # MVA can handle most closed and some open networks
        if self._analyzer.is_closed_model() or self._analyzer.is_open_model():
            if not self._analyzer.has_fork_join():
                self._candidate_solvers.append('MVA')

        # DES can handle most networks
        if not self._analyzer.has_cache():
            self._candidate_solvers.append('DES')

        # NC for product-form networks
        if self._analyzer.has_product_form() and self._analyzer.is_closed_model():
            self._candidate_solvers.append('NC')

        # Default fallback
        if not self._candidate_solvers:
            self._candidate_solvers.append('MVA')

    def _select_solver_heuristic(self) -> str:
        """
        Select solver using heuristic rules.

        Solver selection heuristic
        """
        # Check for forced solver
        if self.options.forced_solver:
            return self.options.forced_solver.upper()

        analyzer = self._analyzer

        # Open networks
        if analyzer.is_open_model():
            if 'MVA' in self._candidate_solvers:
                return 'MVA'
            return 'DES'

        # Mixed networks - use simulation
        if analyzer.is_mixed_model():
            return 'DES'

        # Closed networks
        if analyzer.is_closed_model():
            # Single chain - NC is fastest
            if analyzer.has_single_chain():
                if 'NC' in self._candidate_solvers:
                    return 'NC'

            # All infinite servers - MVA is simple
            if analyzer.has_all_infinite_servers():
                return 'MVA'

            # Product form without multi-server - use NC
            if analyzer.has_product_form() and not analyzer.has_multi_server():
                if 'NC' in self._candidate_solvers:
                    return 'NC'

            # FCFS networks
            if analyzer.has_homogeneous_scheduling(SchedStrategy.FCFS):
                avg_jobs = analyzer.get_avg_jobs_per_chain()
                total_jobs = analyzer.get_total_jobs()

                if not analyzer.has_multi_server():
                    # Heavy load - Fluid is accurate
                    if avg_jobs > 30:
                        if 'Fluid' in self._candidate_solvers:
                            return 'Fluid'
                        return 'MVA'

                    # Medium load
                    if avg_jobs > 10:
                        return 'MVA'

                    # Light load - NC avoids AMVA errors
                    if total_jobs < 5:
                        if 'NC' in self._candidate_solvers:
                            return 'NC'

                    return 'MVA'

                else:
                    # Multi-server FCFS
                    if total_jobs < 5:
                        if 'NC' in self._candidate_solvers:
                            return 'NC'
                    return 'MVA'

            # PS/PSPRIO with multi-server
            if (analyzer.has_homogeneous_scheduling(SchedStrategy.PS) and
                    analyzer.has_multi_server()):
                return 'MVA'

            # Default for closed networks
            return 'MVA'

        # Fallback
        return 'MVA'

    def _create_solver(self, solver_name: str, method: str = 'default') -> Any:
        """Create an instance of the specified solver.

        Args:
            solver_name: Name of the solver to create (e.g., 'MVA', 'NC')
            method: Sub-method to use (e.g., 'lin', 'exact', 'amva')
        """
        solver_name = solver_name.upper()

        # Merge method into kwargs for solvers that support it
        solver_kwargs = dict(self.kwargs)
        if method != 'default':
            solver_kwargs['method'] = method

        if solver_name == 'MVA':
            from .solver_mva import SolverMVA
            return SolverMVA(self.model, **solver_kwargs)

        elif solver_name == 'DES':
            from .solver_des import SolverDES, DESOptions
            options = DESOptions(**{k: v for k, v in self.kwargs.items()
                                   if hasattr(DESOptions, k)})
            return SolverDES(self._sn, options)

        elif solver_name == 'NC':
            # NC not yet implemented natively, fall back to MVA
            from .solver_mva import SolverMVA
            return SolverMVA(self.model, **solver_kwargs)

        elif solver_name in ('FLUID', 'FLD'):
            # Fluid not yet implemented natively, fall back to MVA with amva method
            from .solver_mva import SolverMVA
            if method == 'default':
                solver_kwargs['method'] = 'amva'
            return SolverMVA(self.model, **solver_kwargs)

        else:
            # Default to MVA
            from .solver_mva import SolverMVA
            return SolverMVA(self.model, **solver_kwargs)

    def runAnalyzer(self):
        """
        Run the analysis with the automatically selected solver.

        Returns:
            self for method chaining
        """
        # Select solver based on method
        method = self.options.selection_method

        # Parse solver.method syntax (e.g., 'mva.lin', 'nc.exact')
        if '.' in method:
            solver_part, sub_method = method.split('.', 1)
            self._solver_method = sub_method
        else:
            solver_part = method
            self._solver_method = 'default'

        # Handle high-level selection methods
        if solver_part == 'sim':
            # Best simulator - DES
            self._selected_solver_name = 'DES'
        elif solver_part == 'exact':
            # Best exact method - NC for closed, MVA otherwise
            if self._analyzer.is_closed_model():
                self._selected_solver_name = 'NC'
            else:
                self._selected_solver_name = 'MVA'
        elif solver_part == 'fast':
            # Fast approximate method - MVA
            self._selected_solver_name = 'MVA'
        elif solver_part == 'accurate':
            # Accurate approximate method - Fluid
            self._selected_solver_name = 'Fluid'
        elif solver_part in ('default', 'heur', 'auto', 'line'):
            self._selected_solver_name = self._select_solver_heuristic()
        # Handle specific solver targeting
        elif solver_part.upper() in ['MVA', 'NC', 'CTMC', 'FLUID', 'SSA', 'JMT', 'DES', 'MAM']:
            self._selected_solver_name = solver_part.upper()
            if self._selected_solver_name == 'FLUID':
                self._selected_solver_name = 'Fluid'
        else:
            # Unknown method - fall back to heuristic
            self._selected_solver_name = self._select_solver_heuristic()

        if self.options.verbose:
            print(f"SolverAuto: Selected solver '{self._selected_solver_name}'")
            if self._solver_method != 'default':
                print(f"SolverAuto: Using method '{self._solver_method}'")
            print(f"SolverAuto: Candidates were {self._candidate_solvers}")

        # Create and run solver
        try:
            self._selected_solver = self._create_solver(self._selected_solver_name, self._solver_method)
            self._selected_solver.runAnalyzer()
            self._result = self._selected_solver
        except Exception as e:
            # Try fallback to other candidates
            for candidate in self._candidate_solvers:
                if candidate != self._selected_solver_name:
                    try:
                        if self.options.verbose:
                            print(f"SolverAuto: Trying fallback to '{candidate}'")
                        self._selected_solver = self._create_solver(candidate, self._solver_method)
                        self._selected_solver.runAnalyzer()
                        self._selected_solver_name = candidate
                        self._result = self._selected_solver
                        break
                    except Exception:
                        continue
            else:
                raise RuntimeError(f"All solver candidates failed. Last error: {e}")

        return self

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get comprehensive average performance metrics table.

        Returns:
            pandas.DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        # Suppress output from delegated solver
        if hasattr(self._selected_solver, '_table_silent'):
            self._selected_solver._table_silent = True

        df = self._selected_solver.getAvgTable()

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        return df

    def getSelectedSolverName(self) -> str:
        """Get the name of the automatically selected solver."""
        if self._selected_solver_name is None:
            self._selected_solver_name = self._select_solver_heuristic()
        return self._selected_solver_name

    def getCandidateSolverNames(self) -> List[str]:
        """Get list of candidate solver names for this model."""
        return self._candidate_solvers.copy()

    def getModel(self) -> Any:
        """Get the network model being solved."""
        return self.model

    def getStruct(self, *args, **kwargs) -> Any:
        """Get the network structure from the model."""
        if hasattr(self.model, 'getStruct'):
            return self.model.getStruct(*args, **kwargs)
        return self._sn

    def getName(self) -> str:
        """Get solver name."""
        return 'SolverAuto'

    def setForcedSolver(self, solver_name: str) -> None:
        """Force a specific solver to be used."""
        self.options.forced_solver = solver_name

    def setSelectionMethod(self, method: str) -> None:
        """Set the solver selection method."""
        self.options.selection_method = method

    def setOptions(self, options: Optional[SolverAutoOptions] = None, **kwargs) -> None:
        """
        Set solver options.

        Args:
            options: SolverAutoOptions instance or None to update current
            **kwargs: Individual option updates (selection_method, forced_solver, verbose)
        """
        if options is not None:
            self.options = options
        else:
            # Update individual options
            if 'selection_method' in kwargs:
                self.options.selection_method = kwargs['selection_method']
            if 'forced_solver' in kwargs:
                self.options.forced_solver = kwargs['forced_solver']
            if 'verbose' in kwargs:
                self.options.verbose = kwargs['verbose']

    def supports(self, model: Any) -> bool:
        """
        Check if SolverAuto can handle this model.

        SolverAuto supports any model that at least one candidate solver supports.

        Args:
            model: Network model to check

        Returns:
            True if any candidate solver supports the model, False otherwise
        """
        if not self._candidate_solvers:
            return False

        for candidate in self._candidate_solvers:
            try:
                solver = self._create_solver(candidate)
                if hasattr(solver, 'supports') and callable(solver.supports):
                    if solver.supports(model):
                        return True
            except Exception:
                pass

        # At least MVA should support most models
        return len(self._candidate_solvers) > 0

    def hasResults(self) -> bool:
        """Check if analysis has been completed and results are available."""
        return self._result is not None

    def getResults(self) -> Any:
        """
        Get the solver results object.

        Returns:
            The selected solver or None if not yet run
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result

    def listValidMethods(self) -> List[str]:
        """
        List valid selection methods for SolverAuto.

        Returns:
            List of valid method names for solver selection:
            - 'default', 'heur', 'auto': Heuristic-based selection
            - 'sim': Best simulator (DES/JMT)
            - 'exact': Best exact method (NC for closed, MVA otherwise)
            - 'fast': Fast approximate method (MVA)
            - 'accurate': Accurate approximate method (Fluid)
            - Direct solver targets: 'mva', 'nc', 'ctmc', 'fluid', 'ssa', 'jmt', 'des', 'mam'
        """
        return [
            # High-level selection methods
            'default', 'heur', 'auto', 'line',
            # Method shortcuts
            'sim',      # Best simulator
            'exact',    # Best exact method
            'fast',     # Fast approximate (MVA)
            'accurate', # Accurate approximate (Fluid)
            # Direct solver targeting
            'mva', 'nc', 'ctmc', 'fluid', 'ssa', 'jmt', 'des', 'mam',
        ]

    # ========== Delegation Methods for NetworkSolver Parity ==========

    def _choose_solver_for_method(self, method_name: str) -> Optional[str]:
        """
        Heuristically select the best solver for a specific method.

        Based on MATLAB SolverAuto.chooseSolver() logic:
        - Probability methods (getProb*) prefer NC/CTMC
        - Distribution methods (getCdf*) prefer Fluid/JMT
        - Transient methods (getTran*) prefer Fluid/JMT
        - Average methods prefer MVA (fast)

        Args:
            method_name: Name of the method to optimize for

        Returns:
            Solver name to prioritize or None for general heuristic
        """
        method_lower = method_name.lower()

        # Probability methods -> NC or CTMC
        if 'prob' in method_lower:
            if 'NC' in self._candidate_solvers:
                return 'NC'
            if 'CTMC' in self._candidate_solvers:
                return 'CTMC'

        # Distribution/CDF methods -> Fluid or JMT (simulation)
        if 'cdf' in method_lower or 'perct' in method_lower:
            if 'Fluid' in self._candidate_solvers:
                return 'Fluid'
            if 'JMT' in self._candidate_solvers:
                return 'JMT'
            if 'DES' in self._candidate_solvers:
                return 'DES'

        # Transient methods -> Fluid or JMT
        if 'tran' in method_lower:
            if 'Fluid' in self._candidate_solvers:
                return 'Fluid'
            if 'JMT' in self._candidate_solvers:
                return 'JMT'
            if 'DES' in self._candidate_solvers:
                return 'DES'

        # Sample methods -> DES or JMT
        if 'sample' in method_lower:
            if 'DES' in self._candidate_solvers:
                return 'DES'
            if 'JMT' in self._candidate_solvers:
                return 'JMT'

        # Average methods -> MVA (fastest)
        if 'avg' in method_lower or 'util' in method_lower or 'tput' in method_lower:
            if 'MVA' in self._candidate_solvers:
                return 'MVA'

        # Default to chosen solver
        return None

    def _delegate(self, method_name: str, *args, num_returns: int = 1, **kwargs) -> Any:
        """
        Delegate method call to the selected solver with intelligent fallback.

        Based on MATLAB SolverAuto.delegate() logic:
        1. Determine proposed solvers: best solver for method + candidates
        2. Try each solver in order
        3. Support variable return values (1-7)
        4. On failure, try next candidate

        Args:
            method_name: Name of the method to call
            *args: Positional arguments
            num_returns: Expected number of return values (default: 1)
            **kwargs: Keyword arguments

        Returns:
            Result from the selected solver's method, handling variable returns
        """
        if self._selected_solver is None:
            self.runAnalyzer()

        # Determine proposed solvers to try
        proposed_solvers = []

        # First priority: method-specific best solver
        best_for_method = self._choose_solver_for_method(method_name)
        if best_for_method and best_for_method != self._selected_solver_name:
            try:
                best_solver = self._create_solver(best_for_method)
                if best_solver.supports(self.model) if hasattr(best_solver, 'supports') else True:
                    proposed_solvers.append(best_solver)
            except Exception:
                pass  # If can't create, skip to next

        # Second priority: currently selected solver
        proposed_solvers.append(self._selected_solver)

        # Third priority: other candidates
        for candidate in self._candidate_solvers:
            if candidate != self._selected_solver_name and candidate != best_for_method:
                try:
                    solver = self._create_solver(candidate)
                    proposed_solvers.append(solver)
                except Exception:
                    pass  # Skip if can't create

        # Try each proposed solver in order
        last_error = None
        for solver in proposed_solvers:
            try:
                if not (hasattr(solver, 'supports') and callable(solver.supports)):
                    # Solver doesn't support method or check fails, try anyway
                    pass
                elif not solver.supports(self.model):
                    # Skip this solver if it doesn't support the model
                    continue

                if not hasattr(solver, method_name):
                    # Method doesn't exist on this solver, try next
                    if self.options.verbose:
                        print(f"SolverAuto: Method '{method_name}' unsupported by {solver.__class__.__name__}")
                    continue

                # Call the method
                method = getattr(solver, method_name)
                result = method(*args, **kwargs)

                if self.options.verbose:
                    print(f"SolverAuto: Method '{method_name}' succeeded with {solver.__class__.__name__}")

                # Update selected solver if different
                if solver != self._selected_solver:
                    self._selected_solver = solver
                    self._selected_solver_name = solver.__class__.__name__.replace('SolverNative', '').replace('Solver', '')

                return result

            except Exception as e:
                last_error = e
                if self.options.verbose:
                    print(f"SolverAuto: {solver.__class__.__name__} failed for '{method_name}': {str(e)}")
                continue

        # All solvers failed
        error_msg = f"All solvers failed for method '{method_name}'"
        if last_error:
            error_msg += f". Last error: {str(last_error)}"
        raise RuntimeError(error_msg)

    # ========== Basic Metric Methods ==========

    def getAvg(self):
        """Get all average metrics (QN, UN, RN, TN, AN, WN)."""
        return self._delegate('getAvg')

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths at steady-state."""
        return self._delegate('getAvgQLen')

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations at steady-state."""
        return self._delegate('getAvgUtil')

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times at steady-state."""
        return self._delegate('getAvgRespT')

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times at steady-state."""
        return self._delegate('getAvgResidT')

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times (queue time excluding service)."""
        return self._delegate('getAvgWaitT')

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs at steady-state."""
        return self._delegate('getAvgTput')

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates at steady-state."""
        return self._delegate('getAvgArvR')

    # ========== Chain Methods ==========

    def getAvgChain(self):
        """Get average metrics by chain."""
        return self._delegate('getAvgChain')

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue length by chain."""
        return self._delegate('getAvgQLenChain')

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilization by chain."""
        return self._delegate('getAvgUtilChain')

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response time by chain."""
        return self._delegate('getAvgRespTChain')

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence time by chain."""
        return self._delegate('getAvgResidTChain')

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughput by chain."""
        return self._delegate('getAvgTputChain')

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rate by chain."""
        return self._delegate('getAvgArvRChain')

    # ========== Node Methods ==========

    def getAvgNode(self):
        """Get average metrics by node."""
        return self._delegate('getAvgNode')

    def getAvgNodeQLenChain(self) -> np.ndarray:
        """Get average node queue length by chain."""
        return self._delegate('getAvgNodeQLenChain')

    def getAvgNodeUtilChain(self) -> np.ndarray:
        """Get average node utilization by chain."""
        return self._delegate('getAvgNodeUtilChain')

    def getAvgNodeRespTChain(self) -> np.ndarray:
        """Get average node response time by chain."""
        return self._delegate('getAvgNodeRespTChain')

    def getAvgNodeResidTChain(self) -> np.ndarray:
        """Get average node residence time by chain."""
        return self._delegate('getAvgNodeResidTChain')

    def getAvgNodeTputChain(self) -> np.ndarray:
        """Get average node throughput by chain."""
        return self._delegate('getAvgNodeTputChain')

    def getAvgNodeArvRChain(self) -> np.ndarray:
        """Get average node arrival rate by chain."""
        return self._delegate('getAvgNodeArvRChain')

    # ========== System Methods ==========

    def getAvgSys(self):
        """Get system-level average metrics."""
        return self._delegate('getAvgSys')

    def getAvgSysRespT(self) -> float:
        """Get system average response time."""
        return self._delegate('getAvgSysRespT')

    def getAvgSysTput(self) -> float:
        """Get system average throughput."""
        return self._delegate('getAvgSysTput')

    # ========== Table Methods ==========

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics table by chain."""
        return self._delegate('getAvgChainTable')

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system average metrics table."""
        return self._delegate('getAvgSysTable')

    def getAvgNodeTable(self) -> pd.DataFrame:
        """Get average metrics table by node."""
        return self._delegate('getAvgNodeTable')

    def getAvgQLenTable(self) -> pd.DataFrame:
        """Get average queue length table."""
        return self._delegate('getAvgQLenTable')

    def getAvgUtilTable(self) -> pd.DataFrame:
        """Get average utilization table."""
        return self._delegate('getAvgUtilTable')

    def getAvgRespTTable(self) -> pd.DataFrame:
        """Get average response time table."""
        return self._delegate('getAvgRespTTable')

    def getAvgTputTable(self) -> pd.DataFrame:
        """Get average throughput table."""
        return self._delegate('getAvgTputTable')

    # ========== Distribution Methods ==========

    def getCdfRespT(self):
        """Get CDF of response times at steady-state."""
        return self._delegate('getCdfRespT')

    def getCdfPassT(self):
        """Get CDF of passage times at steady-state."""
        return self._delegate('getCdfPassT')

    def getPerctRespT(self, percentiles: List[float]):
        """Get response time percentiles."""
        return self._delegate('getPerctRespT', percentiles)

    # ========== Transient Methods ==========

    def getTranAvg(self):
        """Get transient average metrics."""
        return self._delegate('getTranAvg')

    def getTranCdfRespT(self):
        """Get transient CDF of response times."""
        return self._delegate('getTranCdfRespT')

    def getTranCdfPassT(self):
        """Get transient CDF of passage times."""
        return self._delegate('getTranCdfPassT')

    # ========== Probability Methods ==========

    def getProb(self, node, state):
        """Get marginal state probability for a node."""
        return self._delegate('getProb', node, state)

    def getProbAggr(self, node, state_a):
        """Get aggregated state probability for a node."""
        return self._delegate('getProbAggr', node, state_a)

    def getProbSys(self):
        """Get joint system state probability."""
        return self._delegate('getProbSys')

    def getProbSysAggr(self):
        """Get aggregated system state probability."""
        return self._delegate('getProbSysAggr')

    def getProbMarg(self, node, jobclass, state_m):
        """Get marginalized state probability for station and class."""
        return self._delegate('getProbMarg', node, jobclass, state_m)

    def getProbNormConstAggr(self):
        """Get normalizing constant (log)."""
        return self._delegate('getProbNormConstAggr')

    # ========== Transient Probability Methods ==========

    def getTranProb(self, node):
        """Get transient state probabilities for a node."""
        return self._delegate('getTranProb', node)

    def getTranProbAggr(self, node):
        """Get transient aggregated probabilities for a node."""
        return self._delegate('getTranProbAggr', node)

    def getTranProbSys(self):
        """Get transient system probabilities."""
        return self._delegate('getTranProbSys')

    def getTranProbSysAggr(self):
        """Get transient system aggregated probabilities."""
        return self._delegate('getTranProbSysAggr')

    # ========== Sampling Methods ==========

    def sample(self, node, num_events: int = 1000):
        """Sample node states."""
        return self._delegate('sample', node, num_events)

    def sampleAggr(self, node, num_events: int = 1000):
        """Sample aggregated node states."""
        return self._delegate('sampleAggr', node, num_events)

    def sampleSys(self, num_events: int = 1000):
        """Sample system states."""
        return self._delegate('sampleSys', num_events)

    def sampleSysAggr(self, num_events: int = 1000):
        """Sample system aggregated states."""
        return self._delegate('sampleSysAggr', num_events)

    @classmethod
    def load(cls, *args, **kwargs) -> Any:
        """
        Factory method to load models or instantiate solvers.

        Usage 1: Load a model from file.
            model = LINE.load(filename)
            model = LINE.load(filename, verbose=True)
        Supported formats: .jsim, .jsimg, .jsimw (JMT), .xml/.lqn/.lqnx (LQN),
                          .mat (MATLAB), .pkl/.pickle (Python pickle)

        Usage 2: Instantiate a solver with a specific method.
            solver = LINE.load(method, model, **options)

        Args:
            For file loading:
                filename: Path to model file
                verbose: Print loading info (default False)
            For solver loading:
                method: Solver method name (e.g., 'mva', 'ctmc', 'ssa', 'fluid', 'nc', 'mam', 'jmt')
                model: Network model to solve
                **options: Solver options

        Returns:
            Loaded model object or Solver instance
        """
        import os

        if len(args) == 0:
            raise ValueError("LINE.load requires at least one argument")

        first_arg = args[0]

        # Check if this is file loading (first arg is a filename string with extension)
        if isinstance(first_arg, str) and '.' in first_arg:
            ext = os.path.splitext(first_arg)[1].lower()
            supported_extensions = {'.jsim', '.jsimg', '.jsimw', '.jmva',
                                   '.xml', '.lqn', '.lqnx', '.mat', '.pkl', '.pickle'}

            if ext in supported_extensions:
                filename = first_arg
                verbose = args[1] if len(args) > 1 else kwargs.get('verbose', False)

                if not os.path.exists(filename):
                    raise FileNotFoundError(f"File not found: {filename}")

                # Load based on format
                if ext in ('.jsim', '.jsimg', '.jsimw'):
                    # Load JSIM/JMT model
                    from ..api.io import jsim2line
                    return jsim2line(filename)

                elif ext == '.jmva':
                    # Load JMVA model
                    from ..api.io import jmva2line
                    return jmva2line(filename)

                elif ext in ('.xml', '.lqn', '.lqnx'):
                    # Load LQN model
                    from ..api.io import lqn2qn
                    return lqn2qn(filename)

                elif ext == '.mat':
                    # Load MATLAB saved model
                    try:
                        from scipy.io import loadmat
                        data = loadmat(filename)
                        # Try to find a Network object
                        if 'model' in data:
                            result = data['model']
                        elif 'network' in data:
                            result = data['network']
                        else:
                            # Return the raw data
                            result = data
                        if verbose:
                            print(f"Loaded model from {filename}")
                        return result
                    except ImportError:
                        raise ImportError("scipy is required to load .mat files")

                elif ext in ('.pkl', '.pickle'):
                    import pickle
                    with open(filename, 'rb') as f:
                        result = pickle.load(f)
                    if verbose:
                        print(f"Loaded from {filename}")
                    return result

        # Otherwise, this is solver loading: LINE.load(method, model, **options)
        if len(args) < 2:
            raise ValueError("Solver loading requires method and model arguments: LINE.load(method, model)")

        method = args[0]
        model = args[1]
        options = kwargs

        # Parse the method and create appropriate solver
        method_lower = method.lower()

        # Strip solver prefix if present (e.g., 'mva.exact' -> 'exact')
        if method_lower in ('default', 'auto', 'line'):
            return cls(model, **options)

        elif method_lower.startswith('ctmc'):
            from .solver_ctmc import SolverCTMC
            sub_method = method_lower.replace('ctmc.', '').replace('ctmc', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverCTMC(model, **options)

        elif method_lower.startswith('mva') or method_lower in (
            'amva', 'qna', 'sqrt', 'mm1', 'mmk', 'mg1', 'gm1', 'gig1', 'gim1',
            'gigk', 'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower',
            'gb.upper', 'gb.lower', 'pb.upper', 'pb.lower', 'sb.upper', 'sb.lower'
        ):
            from .solver_mva import SolverMVA
            sub_method = method_lower.replace('mva.', '').replace('mva', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverMVA(model, **options)

        elif method_lower.startswith('ssa') or method_lower in ('nrm', 'serial', 'parallel'):
            from .solver_ssa import SolverSSA
            sub_method = method_lower.replace('ssa.', '').replace('ssa', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverSSA(model, **options)

        elif method_lower.startswith('jmt') or method_lower in ('jsim', 'jmva'):
            from .solver_jmt import SolverJMT
            sub_method = method_lower.replace('jmt.', '').replace('jmt', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverJMT(model, **options)

        elif method_lower.startswith('fluid'):
            from .solver_fluid import SolverFluid
            sub_method = method_lower.replace('fluid.', '').replace('fluid', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverFluid(model, **options)

        elif method_lower.startswith('nc') or method_lower in (
            'comom', 'comomld', 'cub', 'ls', 'le', 'mmint2', 'panacea', 'pana'
        ):
            from .solver_nc import SolverNC
            sub_method = method_lower.replace('nc.', '').replace('nc', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverNC(model, **options)

        elif method_lower.startswith('mam'):
            from .solver_mam import SolverMAM
            sub_method = method_lower.replace('mam.', '').replace('mam', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverMAM(model, **options)

        elif method_lower.startswith('des'):
            from .solver_des import SolverDES
            return SolverDES(model, **options)

        else:
            # Unknown method - use SolverAuto with the method as selection
            options['method'] = method
            return cls(model, **options)

    # ========== Python-style Aliases (snake_case) ==========

    # Aliases (Python naming convention)
    run_analyzer = runAnalyzer
    get_name = getName
    get_avg_table = getAvgTable
    avg_table = getAvgTable  # Method alias, not property
    avg_chain_table = getAvgChainTable
    avg_sys_table = getAvgSysTable
    get_selected_solver_name = getSelectedSolverName
    get_candidate_solver_names = getCandidateSolverNames
    set_forced_solver = setForcedSolver
    set_selection_method = setSelectionMethod

    # Basic metrics
    get_avg = getAvg
    get_avg_qlen = getAvgQLen
    get_avg_util = getAvgUtil
    get_avg_resp_t = getAvgRespT
    get_avg_resid_t = getAvgResidT
    get_avg_wait_t = getAvgWaitT
    get_avg_tput = getAvgTput
    get_avg_arv_r = getAvgArvR

    # Chain metrics
    get_avg_chain = getAvgChain
    get_avg_qlen_chain = getAvgQLenChain
    get_avg_util_chain = getAvgUtilChain
    get_avg_resp_t_chain = getAvgRespTChain
    get_avg_resid_t_chain = getAvgResidTChain
    get_avg_tput_chain = getAvgTputChain
    get_avg_arv_r_chain = getAvgArvRChain

    # Node metrics
    get_avg_node = getAvgNode
    get_avg_node_qlen_chain = getAvgNodeQLenChain
    get_avg_node_util_chain = getAvgNodeUtilChain
    get_avg_node_resp_t_chain = getAvgNodeRespTChain
    get_avg_node_resid_t_chain = getAvgNodeResidTChain
    get_avg_node_tput_chain = getAvgNodeTputChain
    get_avg_node_arv_r_chain = getAvgNodeArvRChain

    # System metrics
    get_avg_sys = getAvgSys
    get_avg_sys_resp_t = getAvgSysRespT
    get_avg_sys_tput = getAvgSysTput

    # Table methods
    get_avg_chain_table = getAvgChainTable
    get_avg_sys_table = getAvgSysTable
    get_avg_node_table = getAvgNodeTable
    get_avg_qlen_table = getAvgQLenTable
    get_avg_util_table = getAvgUtilTable
    get_avg_resp_t_table = getAvgRespTTable
    get_avg_tput_table = getAvgTputTable

    # Distribution methods
    get_cdf_resp_t = getCdfRespT
    get_cdf_pass_t = getCdfPassT
    get_perct_resp_t = getPerctRespT

    # Transient methods
    get_tran_avg = getTranAvg
    get_tran_cdf_resp_t = getTranCdfRespT
    get_tran_cdf_pass_t = getTranCdfPassT

    # Probability methods
    get_prob = getProb
    get_prob_aggr = getProbAggr
    get_prob_sys = getProbSys
    get_prob_sys_aggr = getProbSysAggr
    get_prob_marg = getProbMarg
    get_prob_norm_const_aggr = getProbNormConstAggr

    # Transient probability methods
    get_tran_prob = getTranProb
    get_tran_prob_aggr = getTranProbAggr
    get_tran_prob_sys = getTranProbSys
    get_tran_prob_sys_aggr = getTranProbSysAggr

    # Sampling methods
    sample_aggr = sampleAggr
    sample_sys = sampleSys
    sample_sys_aggr = sampleSysAggr


# Alias for convenience
LINE = SolverAuto


__all__ = [
    'SolverAuto',
    'SolverAutoOptions',
    'ModelAnalyzer',
    'SolverType',
    'LINE',
]
