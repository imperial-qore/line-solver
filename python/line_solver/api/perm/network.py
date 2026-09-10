"""
Permanent-based marginal probabilities of closed queueing networks.

Twin of jline.lib.perm.QueueingNetwork, jline.lib.perm.NetworkNoThink and
jline.lib.perm.NetworkThink. The marginal probability of a queue-length vector
in a closed product-form network is a permanent of the service demand matrix
replicated according to the state, so any of the permanent solvers of this
package can be plugged in, exact or approximate.
"""

from dataclasses import dataclass
from itertools import product

import numpy as np

from .approx import BethePermanent
from .base import PermSolver, require_full_support
from .exact import NaivePermanent, RyzerPermanent
from .sampling import AdaPartSampler, HuberLawSampler

MIN_VALUE = 2.220446049250314e-16


def preprocessing_ds(m: np.ndarray, tolerance: float = 0.001,
                     max_iterations: int = 10000):
    """
    Make a matrix doubly stochastic with the Sinkhorn algorithm.

    Args:
        m: Square nonnegative matrix
        tolerance: Row and column sum error below which the iteration stops
        max_iterations: Sweep budget; exceeding it raises rather than returning

    Returns:
        Tuple of (doubly stochastic matrix, rescaling factor) with DS = X A Y,
        so that perm(A) = perm(DS) / factor

    Raises:
        ValueError: If the scaling has not converged within max_iterations.
            Sinkhorn converges to a doubly stochastic limit if and only if the
            matrix has total support; a matrix with a strictly positive
            permanent can still fail that, and then the margins stall on the
            tolerance instead of converging.
    """
    m = np.asarray(m, dtype=float)
    n = m.shape[0]
    # A zero used to be floored to MIN_VALUE here. That is the same
    # non-invertible substitution the four approximations used to make: it
    # changes the permanent by n!*eps, and it also hides the real precondition,
    # since flooring manufactures total support and the scaling then converges
    # on a matrix that never had it. Refuse instead.
    if n > 0:
        require_full_support(m, 'preprocessing_ds')
    result = m.copy()

    x = np.eye(n)
    y = np.eye(n)

    max_row_error = np.inf
    max_col_error = np.inf

    sweeps = 0
    while max_row_error > tolerance or max_col_error > tolerance:
        sweeps += 1
        if sweeps > max_iterations:
            raise ValueError(
                "preprocessing_ds did not converge in " + str(max_iterations) +
                " sweeps (row error " + str(max_row_error) + ", column error " +
                str(max_col_error) + " against a tolerance of " + str(tolerance) +
                "). The usual cause is a matrix without total support.")
        col_sums = result.sum(axis=0)
        for j in range(n):
            if col_sums[j] > 0:
                result[:, j] /= col_sums[j]
                y[:, j] /= col_sums[j]

        row_sums = result.sum(axis=1)
        for i in range(n):
            if row_sums[i] > 0:
                result[i, :] /= row_sums[i]
                x[i, :] /= row_sums[i]

        max_col_error = float(np.max(np.abs(result.sum(axis=0) - 1.0)))
        max_row_error = float(np.max(np.abs(result.sum(axis=1) - 1.0)))

    rescaling_factor = float(np.prod(np.diagonal(x) * np.diagonal(y)))
    return result, rescaling_factor


@dataclass
class MarginalResult:
    """Marginal probability of a state with the cost of computing it."""
    probability: float
    time_ms: float
    memory_bytes: int


class NetworkNoThink:
    """
    Closed queueing network without think time.

    Port of jline.lib.perm.NetworkNoThink. joint() evaluates the unnormalized
    product-form joint probability of a per-queue per-class state, marginal()
    evaluates the unnormalized marginal probability of a per-queue state as a
    permanent of the replicated demand matrix.
    """

    def __init__(self, number_of_queues: int, number_of_classes: int,
                 number_per_class, mean_service_demand, progress: bool = True):
        """
        Initialize the network.

        Args:
            number_of_queues: Number of queueing stations
            number_of_classes: Number of job classes
            number_per_class: Population of each class
            mean_service_demand: Mean service demand per queue and per queue
            progress: If True, generate_marginal prints the visited states
        """
        self.number_of_queues = number_of_queues
        self.number_of_classes = number_of_classes
        self.number_per_class = np.asarray(number_per_class, dtype=int)
        self.mean_service_demand = np.asarray(mean_service_demand, dtype=float)
        self.progress = progress

    def joint(self, state) -> float:
        """
        Unnormalized joint probability of a queue-by-class state.

        Args:
            state: (number_of_queues, number_of_classes) job counts

        Returns:
            Unnormalized probability
        """
        state = np.asarray(state, dtype=int)
        term_msd = float(np.prod(self.mean_service_demand ** state))
        inv_factorial = float(np.prod([_factorial(v) for v in state.ravel()]))
        factorial_term = float(np.prod([_factorial(int(row.sum())) for row in state]))
        return term_msd * factorial_term / inv_factorial

    def marginal(self, solver: PermSolver, state, preprocessing: bool = False) -> MarginalResult:
        """
        Unnormalized marginal probability of a per-queue state.

        Args:
            solver: Permanent solver instance selecting the algorithm to use
            state: Job count of each queue
            preprocessing: If True, Sinkhorn scale the matrix before solving

        Returns:
            MarginalResult with the probability and the solver cost
        """
        matrix = _replicate_demands(self.mean_service_demand, state, self.number_of_queues)
        return _marginal_from_matrix(solver, matrix, state, preprocessing)

    def generate_marginal(self, solver: PermSolver, preprocessing: bool = False):
        """
        Marginal probability of every reachable per-queue state.

        Args:
            solver: Permanent solver instance selecting the algorithm to use
            preprocessing: If True, Sinkhorn scale the matrix before solving

        Returns:
            Dict mapping each state tuple to its MarginalResult
        """
        results = {}
        for state in _all_marginal_states(self.number_of_queues, int(self.number_per_class.sum())):
            results[state] = self.marginal(solver, state, preprocessing)
            if self.progress:
                print("Processed state: " + str(list(state)))
        return results


class NetworkThink:
    """
    Closed queueing network with a think time per class.

    Port of jline.lib.perm.NetworkThink. Identical to NetworkNoThink except
    that the jobs not in a queue are thinking, which contributes think-time
    columns to the permanent matrix and a think-time factor to the joint.
    """

    def __init__(self, number_of_queues: int, number_of_classes: int, number_per_class,
                 think_time, mean_service_demand, progress: bool = True):
        """
        Initialize the network.

        Args:
            number_of_queues: Number of queueing stations
            number_of_classes: Number of job classes
            number_per_class: Population of each class
            think_time: Think time of each class
            mean_service_demand: Mean service demand per queue and per queue
            progress: If True, generate_marginal prints the visited states
        """
        self.number_of_queues = number_of_queues
        self.number_of_classes = number_of_classes
        self.number_per_class = np.asarray(number_per_class, dtype=int)
        self.think_time = np.asarray(think_time, dtype=float)
        self.mean_service_demand = np.asarray(mean_service_demand, dtype=float)
        self.progress = progress

    def joint(self, state) -> float:
        """
        Unnormalized joint probability of a queue-by-class state.

        Args:
            state: (number_of_queues, number_of_classes) job counts

        Returns:
            Unnormalized probability
        """
        state = np.asarray(state, dtype=int)
        term_msd = float(np.prod(self.mean_service_demand ** state))
        inv_factorial = float(np.prod([_factorial(v) for v in state.ravel()]))
        factorial_term = float(np.prod([_factorial(int(row.sum())) for row in state]))
        thinking = self.number_per_class - state.sum(axis=0)
        think_time_term = float(np.prod(self.think_time ** thinking))
        return term_msd * factorial_term * think_time_term / inv_factorial

    def marginal(self, solver: PermSolver, state, preprocessing: bool = False) -> MarginalResult:
        """
        Unnormalized marginal probability of a per-queue state.

        Args:
            solver: Permanent solver instance selecting the algorithm to use
            state: Job count of each queue
            preprocessing: If True, Sinkhorn scale the matrix before solving

        Returns:
            MarginalResult with the probability and the solver cost
        """
        state = np.asarray(state, dtype=int)
        total_jobs = int(state.sum())
        matrix = np.zeros((total_jobs, total_jobs))

        row_index = 0
        for i in range(self.number_of_queues):
            for _ in range(int(state[i])):
                col_index = 0
                for j in range(self.number_of_queues):
                    for _ in range(int(state[j])):
                        matrix[row_index, col_index] = self.mean_service_demand[i, j]
                        col_index += 1
                for j in range(self.number_of_classes):
                    thinking_jobs = int(self.number_per_class[j]) - total_jobs
                    for _ in range(thinking_jobs):
                        matrix[row_index, col_index] = self.think_time[j]
                        col_index += 1
                row_index += 1

        return _marginal_from_matrix(solver, matrix, state, preprocessing)

    def generate_marginal(self, solver: PermSolver, preprocessing: bool = False):
        """
        Marginal probability of every reachable per-queue state.

        Args:
            solver: Permanent solver instance selecting the algorithm to use
            preprocessing: If True, Sinkhorn scale the matrix before solving

        Returns:
            Dict mapping each state tuple to its MarginalResult
        """
        results = {}
        for state in _all_marginal_states(self.number_of_queues, int(self.number_per_class.sum())):
            results[state] = self.marginal(solver, state, preprocessing)
            if self.progress:
                print("Processed state: " + str(list(state)))
        return results


def _replicate_demands(mean_service_demand: np.ndarray, state, number_of_queues: int) -> np.ndarray:
    """Replicate the demand matrix, one row and column per job in the state."""
    state = np.asarray(state, dtype=int)
    total = int(state.sum())
    matrix = np.zeros((total, total))
    row_index = 0
    for i in range(number_of_queues):
        for _ in range(int(state[i])):
            col_index = 0
            for j in range(number_of_queues):
                for _ in range(int(state[j])):
                    matrix[row_index, col_index] = mean_service_demand[i, j]
                    col_index += 1
            row_index += 1
    return matrix


def _marginal_from_matrix(solver: PermSolver, matrix: np.ndarray, state,
                          preprocessing: bool) -> MarginalResult:
    """Solve the replicated matrix with a fresh instance of the solver's class."""
    if preprocessing:
        processed_matrix, rescaling_factor = preprocessing_ds(matrix)
    else:
        processed_matrix, rescaling_factor = matrix, 1.0

    if isinstance(solver, BethePermanent):
        solver_instance = BethePermanent(processed_matrix, 1e-3, 1000, True)
    elif isinstance(solver, NaivePermanent):
        solver_instance = NaivePermanent(processed_matrix, True)
    elif isinstance(solver, RyzerPermanent):
        solver_instance = RyzerPermanent(processed_matrix, "default", True)
    elif isinstance(solver, AdaPartSampler):
        solver_instance = AdaPartSampler(processed_matrix)
        solver_instance.solve()
    elif isinstance(solver, HuberLawSampler):
        solver_instance = HuberLawSampler(processed_matrix)
        solver_instance.solve()
    else:
        raise ValueError("Unsupported solver type")

    factorial_normalization = float(np.prod([_factorial(int(v)) for v in np.asarray(state).ravel()]))
    # preprocessing_ds returns f = prod(diag(X)*diag(Y)) with DS = X A Y, so
    # perm(A) = perm(DS) / f. This used to MULTIPLY, which is wrong by f^2: on
    # the two-station one-class three-job model of the unit test it returned
    # 0.03029 against a truth of 0.464. Both the test that exercises
    # preprocessing_ds and _kb/03-api-layer.md already stated the direction.
    probability = solver_instance.value / rescaling_factor / factorial_normalization
    return MarginalResult(probability, solver_instance.time_ms, solver_instance.memory_bytes)


def _all_marginal_states(number_of_queues: int, total_jobs: int):
    """Enumerate the per-queue job counts summing to total_jobs."""
    states = []
    for prefix in product(*[range(total_jobs + 1)] * (number_of_queues - 1)):
        remaining = total_jobs - sum(prefix)
        if remaining >= 0:
            states.append(tuple(list(prefix) + [remaining]))
    return states


def _factorial(n: int) -> float:
    """Compute n! in double precision."""
    if n <= 1:
        return 1.0
    result = 1.0
    for i in range(2, int(n) + 1):
        result *= i
    return result
