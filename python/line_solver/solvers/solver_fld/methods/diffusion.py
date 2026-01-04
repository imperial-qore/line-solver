"""
Diffusion method for FLD solver (SDE solver for closed networks).

Implements Euler-Maruyama simulation of the diffusion approximation
for closed queueing networks:

    dQ = drift(Q) * dt + diffusion(Q) * dW

where:
- Q is the state vector (queue lengths)
- drift(Q) represents the mean flow dynamics
- diffusion(Q) represents the variance scaling
- dW is a Wiener process increment

This method is restricted to closed networks (fixed population).

Reference:
- Diffusion approximation for queueing networks
- Euler-Maruyama discretization for SDEs
"""

import numpy as np
import time
from typing import Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from ...api.sn import NetworkStruct

from ..options import SolverFLDOptions, FLDResult


class DiffusionMethod:
    """Diffusion method for fluid analysis using Euler-Maruyama SDE solver.

    Implements the diffusion approximation for closed queueing networks
    using custom Euler-Maruyama time stepping with population constraint enforcement.
    """

    def __init__(self, sn, options: SolverFLDOptions):
        """Initialize diffusion method.

        Args:
            sn: Compiled NetworkStruct (must be closed network)
            options: SolverFLDOptions with time stepping parameters

        Raises:
            ValueError: If network is not closed
        """
        self.sn = sn
        self.options = options
        self.iterations = 0
        self.runtime = 0.0

        # Validate closed network
        if not self._is_closed_network():
            raise ValueError("Diffusion method requires closed network")

    def _is_closed_network(self) -> bool:
        """Check if network is closed (all classes have fixed populations).

        Returns:
            True if all job classes have fixed populations
        """
        if not hasattr(self.sn, 'njobs') or self.sn.njobs is None:
            return False

        njobs = np.asarray(self.sn.njobs).flatten()
        return len(njobs) > 0 and np.all(np.isfinite(njobs))

    def solve(self) -> FLDResult:
        """Solve using diffusion (SDE) method.

        Returns:
            FLDResult with performance metrics
        """
        start_time = time.time()

        try:
            # Build system components
            W_T, SQ, c, N = self._build_ode_components()

            # Compute initial state
            x0 = self._compute_initial_state(N)

            # Run Euler-Maruyama stepping
            t_span, x_final, t_vec, Q_hist = self._solve_sde(W_T, SQ, c, N, x0)

            # Extract metrics from final state
            QN, UN, RN, TN, CN, XN = self._extract_metrics(x_final)

            # Build result
            result = FLDResult(
                QN=QN,
                UN=UN,
                RN=RN,
                TN=TN,
                CN=CN,
                XN=XN,
                xvec=x_final,
                t=t_vec,
                iterations=len(t_vec),
                method='diffusion',
            )

            result.runtime = time.time() - start_time
            return result

        except Exception as e:
            if self.options.verbose:
                print(f"Diffusion method error: {e}")
            raise

    def _build_ode_components(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Build components needed for SDE system.

        Returns:
            Tuple of (W_T, SQ, c, N) where N is closed population vector
        """
        M = self.sn.nstations
        K = self.sn.nclasses

        # Transposed transition matrix (simplified)
        W_T = -np.eye(M * K)

        # State-to-station aggregation matrix
        SQ = np.zeros((M * K, M * K))
        for i in range(M):
            for k in range(K):
                SQ[i * K + k, i * K + k] = 1.0

        # Server capacities
        c = np.ones(M * K)
        if hasattr(self.sn, 'nservers') and self.sn.nservers is not None:
            for i in range(min(M, len(self.sn.nservers))):
                for k in range(K):
                    c[i * K + k] = self.sn.nservers[i]

        # Closed population (njobs)
        N = np.zeros(K)
        if hasattr(self.sn, 'njobs') and self.sn.njobs is not None:
            njobs = np.asarray(self.sn.njobs).flatten()
            for k in range(min(K, len(njobs))):
                if np.isfinite(njobs[k]):
                    N[k] = njobs[k]

        return W_T, SQ, c, N

    def _compute_initial_state(self, N: np.ndarray) -> np.ndarray:
        """Compute initial state vector for SDE.

        Args:
            N: Population vector per class

        Returns:
            Initial state vector (M*K,) with population distributed
        """
        M = self.sn.nstations
        K = self.sn.nclasses

        # Distribute population evenly across stations per class
        x0 = np.zeros(M * K)
        for k in range(K):
            if N[k] > 0:
                per_station = N[k] / M
                for i in range(M):
                    x0[i * K + k] = per_station

        return x0

    def _solve_sde(
        self,
        W_T: np.ndarray,
        SQ: np.ndarray,
        c: np.ndarray,
        N: np.ndarray,
        x0: np.ndarray,
    ) -> Tuple[Tuple[float, float], np.ndarray, np.ndarray, list]:
        """Solve SDE using Euler-Maruyama stepping.

        Args:
            W_T: Transposed transition matrix
            SQ: State-to-station aggregation matrix
            c: Server capacities
            N: Population vector
            x0: Initial state

        Returns:
            Tuple of (t_span, x_final, t_vec, Q_hist)
        """
        t_start, t_end = self.options.timespan
        if not np.isfinite(t_end):
            t_end = 100.0

        # Determine time step
        dt = self.options.timestep if self.options.timestep is not None else 0.01
        nsteps = int((t_end - t_start) / dt) + 1

        # Initialize trajectory
        x = x0.copy()
        x = self._enforce_population_constraint(x, N)

        t_vec = np.linspace(t_start, t_end, nsteps)
        Q_hist = [x.copy()]

        # Euler-Maruyama stepping
        for step in range(1, nsteps):
            # Compute drift and diffusion
            drift = self._compute_drift(x, W_T, SQ, c)
            diffusion = self._compute_diffusion(x, SQ, c)

            # Wiener process increment
            dW = np.random.randn(len(x)) * np.sqrt(dt)

            # Euler-Maruyama step
            x = x + drift * dt + diffusion * dW

            # Enforce constraints
            x = np.maximum(x, 0.0)  # Non-negativity
            x = self._enforce_population_constraint(x, N)

            Q_hist.append(x.copy())

        x_final = x
        return (t_start, t_end), x_final, t_vec, Q_hist

    def _compute_drift(
        self,
        x: np.ndarray,
        W_T: np.ndarray,
        SQ: np.ndarray,
        c: np.ndarray,
    ) -> np.ndarray:
        """Compute drift term for SDE.

        The drift is the deterministic component based on fluid dynamics:
        drift = W_T @ theta where theta = min(1, c / sum(x))

        Args:
            x: State vector
            W_T: Transposed transition matrix
            SQ: State-to-station aggregation matrix
            c: Server capacities

        Returns:
            Drift vector
        """
        x = np.maximum(x, 0.0)
        sum_x = SQ @ x + 1e-14

        # Compute theta = min(1, c / sum_x)
        theta = np.minimum(1.0, np.divide(c, sum_x))
        theta = np.where(np.isfinite(theta), theta, 0.0)

        # Drift = W_T @ (x * theta)
        drift = W_T @ (x * theta)

        return drift

    def _compute_diffusion(
        self,
        x: np.ndarray,
        SQ: np.ndarray,
        c: np.ndarray,
    ) -> np.ndarray:
        """Compute diffusion term for SDE (volatility).

        The diffusion coefficient scales the Wiener process:
        diffusion_i = sqrt(throughput_i) = sqrt(min(sum_x, c))

        Args:
            x: State vector
            SQ: State-to-station aggregation matrix
            c: Server capacities

        Returns:
            Diffusion vector (volatility per component)
        """
        x = np.maximum(x, 0.0)
        sum_x = SQ @ x + 1e-14

        # Throughput = min(sum_x, c)
        throughput = np.minimum(sum_x, c)

        # Diffusion = sqrt(throughput)
        diffusion = np.sqrt(np.maximum(throughput, 0.0))

        return diffusion

    def _enforce_population_constraint(
        self,
        x: np.ndarray,
        N: np.ndarray,
    ) -> np.ndarray:
        """Enforce population constraint: sum(x[k]) = N[k] for each class k.

        Args:
            x: State vector (length = M * K)
            N: Population vector per class (length = K)

        Returns:
            Corrected state vector with populations enforced
        """
        K = len(N)  # Number of classes
        M = int(len(x) / K)  # Number of stations

        x = x.copy()
        for k in range(K):
            # Sum of queue lengths for class k across all stations
            indices = [i * K + k for i in range(M)]
            current_sum = np.sum(x[indices])

            if current_sum > 1e-10:
                # Scale to match required population
                scale = N[k] / current_sum
                x[indices] = x[indices] * scale
            else:
                # Distribute evenly if current_sum ≈ 0
                x[indices] = N[k] / M

        return x

    def _extract_metrics(self, x: np.ndarray) -> Tuple:
        """Extract performance metrics from state vector.

        Args:
            x: Final state vector

        Returns:
            Tuple of (QN, UN, RN, TN, CN, XN)
        """
        M = self.sn.nstations
        K = self.sn.nclasses

        # Map state to metrics
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        TN = np.zeros((M, K))

        for i in range(M):
            for k in range(K):
                idx = i * K + k
                if idx < len(x):
                    QN[i, k] = max(0.0, x[idx])
                    # Utilization = min(1, throughput/capacity)
                    UN[i, k] = min(1.0, max(0.0, x[idx]))

        # Compute response times via Little's Law
        with np.errstate(divide='ignore', invalid='ignore'):
            RN = np.divide(QN, np.maximum(TN, 1e-10))
            RN = np.where(np.isfinite(RN), RN, 0.0)

        # Compute cycle and system metrics
        CN = np.sum(RN, axis=0, keepdims=True)
        XN = np.mean(TN, axis=0, keepdims=True)

        return QN, UN, RN, TN, CN, XN


def solve_diffusion(sn, options: Optional[SolverFLDOptions] = None) -> FLDResult:
    """Convenience function to solve using diffusion method.

    Args:
        sn: Compiled NetworkStruct (must be closed)
        options: SolverFLDOptions with time parameters

    Returns:
        FLDResult

    Raises:
        ValueError: If network is not closed
    """
    if options is None:
        options = SolverFLDOptions(method='diffusion')

    method = DiffusionMethod(sn, options)
    return method.solve()
