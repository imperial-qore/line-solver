"""
FCFS approximation utilities for closing method.

Implements moment matching and Coxian distribution fitting for
iterative refinement of non-exponential service processes.

The closing method uses FCFS queue approximations to represent
service processes:

1. Extract throughput (T) and utilization (U) from ODE solution
2. Compute first and second moments using throughput data
3. Fit Coxian distributions to matched moments
4. Update service process representations
5. Check for convergence
"""

import numpy as np
from typing import Tuple, Dict, Optional


def compute_service_moments(
    throughput: np.ndarray,
    utilization: np.ndarray,
    service_rate: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute mean and variance of service times from flow metrics.

    Using Little's Law and utilization:
    - L = lambda * W (queue length = arrival rate * wait time)
    - U = lambda * S (utilization = arrival rate * service time)
    - Therefore: S = U / lambda

    Args:
        throughput: Throughput per station-class (M, K)
        utilization: Utilization per station-class (M, K)
        service_rate: Base service rate per station-class (M, K)

    Returns:
        Tuple of (mean_service_time, variance_service_time) as (M, K) arrays
    """
    M, K = throughput.shape

    # Mean service time: S = U / lambda
    # (where lambda is the effective arrival rate = throughput)
    with np.errstate(divide='ignore', invalid='ignore'):
        mean_S = np.divide(utilization, throughput)
        mean_S = np.where(np.isfinite(mean_S), mean_S, 1.0 / service_rate)
        mean_S = np.maximum(mean_S, 0.0)

    # Variance: Use Poisson assumption (variance = mean for exponential)
    # For general distributions, variance ≥ mean²
    # Start with exponential variance as baseline
    var_S = mean_S ** 2  # Exponential: Var[S] = (E[S])²

    return mean_S, var_S


def fit_coxian_2(mean: float, variance: float) -> Tuple[float, float, float]:
    """
    Fit a 2-phase Coxian distribution to mean and variance.

    A 2-phase Coxian distribution has 3 parameters:
    - mu_1: rate of first phase
    - mu_2: rate of second phase
    - p: probability of going to phase 2 (vs exiting after phase 1)

    Given mean and variance, we solve for these parameters.

    Args:
        mean: Target mean
        variance: Target variance

    Returns:
        Tuple of (mu_1, mu_2, p) parameters
    """
    if mean <= 0:
        return 1.0, 1.0, 0.0  # Default exponential

    if variance < mean ** 2:
        # Variance is too low (less than exponential)
        # Use exponential approximation
        mu = 1.0 / mean
        return mu, mu, 0.0

    # Compute SCV (squared coefficient of variation)
    scv = variance / (mean ** 2)

    if scv <= 1.0:
        # Use exponential
        mu = 1.0 / mean
        return mu, mu, 0.0
    else:
        # Hyperexponential-like fit
        # For 2-phase Coxian: E[S] = 1/mu_1 + p/mu_2
        #                    Var[S] = 1/mu_1² + p*(1+p)/mu_2²

        # Simple approximation: use two exponential phases with different rates
        # Let mu_1 > mu_2 (first phase faster)
        # Set p = (scv - 1) / (scv + 1) as approximation

        p = min(0.9, (scv - 1.0) / (scv + 1.0))  # Cap p at 0.9

        # Solve for mu_1, mu_2 from mean constraint
        # E[S] = 1/mu_1 + p/mu_2 = mean
        # Approximate: mu_1 ≈ 2/mean, mu_2 ≈ 1/(p*mean)

        mu_1 = 2.0 / mean
        mu_2 = max(0.1, 1.0 / (max(p, 0.1) * mean))

        return mu_1, mu_2, p


def update_service_process(
    mean_service_time: float,
    variance_service_time: float,
    current_rate: float,
) -> Dict:
    """
    Update service process representation based on fitted moments.

    Returns a service process dictionary containing distribution type
    and parameters for use in subsequent ODE iterations.

    Args:
        mean_service_time: Mean of service time distribution
        variance_service_time: Variance of service time distribution
        current_rate: Current service rate (for normalization)

    Returns:
        Dictionary with keys: 'type', 'mu', 'variance', 'rate'
    """
    # Fit Coxian distribution
    mu_1, mu_2, p = fit_coxian_2(mean_service_time, variance_service_time)

    # Compute SCV for diagnostics
    scv = variance_service_time / (mean_service_time ** 2) if mean_service_time > 0 else 1.0

    # Choose representation based on SCV
    if scv <= 1.05:
        # Approximately exponential
        process_type = 'Exp'
        params = {'mu': 1.0 / mean_service_time}
    elif scv <= 2.0:
        # Use 2-phase Coxian
        process_type = 'Coxian2'
        params = {'mu_1': mu_1, 'mu_2': mu_2, 'p': p}
    else:
        # Highly variable - use hyperexponential
        process_type = 'HyperExp'
        params = {'mu_1': mu_1, 'mu_2': mu_2, 'p': p}

    return {
        'type': process_type,
        'params': params,
        'mean': mean_service_time,
        'variance': variance_service_time,
        'scv': scv,
        'rate': 1.0 / mean_service_time if mean_service_time > 0 else 1.0,
    }


class FCFSApproximator:
    """Outer iteration loop for closing method FCFS approximation."""

    def __init__(
        self,
        sn,
        options,
        solver_method,
    ):
        """Initialize FCFS approximator.

        Args:
            sn: NetworkStruct
            options: SolverFLDOptions
            solver_method: ODE solver method to iterate
        """
        self.sn = sn
        self.options = options
        self.solver_method = solver_method
        self.M = sn.nstations if hasattr(sn, 'nstations') else 1
        self.K = sn.nclasses if hasattr(sn, 'nclasses') else 1

    def iterate(self, initial_service_rates: Optional[np.ndarray] = None) -> Tuple:
        """
        Run iterative FCFS approximation until convergence.

        Args:
            initial_service_rates: Starting service rates (M, K)

        Returns:
            Tuple of (result, iterations_completed)
        """
        if initial_service_rates is None:
            initial_service_rates = np.ones((self.M, self.K))

        # Initialize service process representations
        eta = initial_service_rates.copy()
        eta_history = [eta.copy()]

        # Main iteration loop
        for iteration in range(self.options.iter_max):
            # Solve ODE with current service processes
            result = self.solver_method.solve()

            if result is None or result.TN is None:
                break

            # Extract metrics
            T = np.mean(result.TN, axis=1, keepdims=True)  # Throughput per station
            U = np.mean(result.UN, axis=1, keepdims=True)  # Utilization per station

            # Compute new service times via moment matching
            new_eta = np.zeros_like(eta)
            for i in range(self.M):
                for r in range(self.K):
                    if T[i, r] > 1e-10:
                        # Service time: S = U / lambda
                        service_time = U[i, r] / T[i, r]
                        # Update rate: eta = 1 / S
                        new_eta[i, r] = 1.0 / max(service_time, 1e-6)
                    else:
                        new_eta[i, r] = eta[i, r]

            # Check convergence
            if self._has_converged(eta, new_eta):
                if self.options.verbose:
                    print(f"FCFS approximation converged at iteration {iteration + 1}")
                return result, iteration + 1

            eta = new_eta
            eta_history.append(eta.copy())

        if self.options.verbose:
            print(f"FCFS approximation completed {self.options.iter_max} iterations")

        return result, self.options.iter_max

    def _has_converged(self, eta_prev: np.ndarray, eta_curr: np.ndarray) -> bool:
        """Check if iteration has converged.

        Args:
            eta_prev: Previous iteration service rates
            eta_curr: Current iteration service rates

        Returns:
            True if relative change < tolerance
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_change = np.abs(1.0 - eta_curr / (eta_prev + 1e-14))
            rel_change = np.where(np.isfinite(rel_change), rel_change, 1e10)

        return np.max(rel_change) < self.options.iter_tol
