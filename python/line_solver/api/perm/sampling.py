"""
Randomized approximations of the matrix permanent.

AdaPartSampler is the adaptive partitioning rejection sampler built on the
Soules column bound, the twin of jline.lib.perm.AdaPartSampler. HuberLawSampler
is the Huber-Law acceptance-rejection importance sampler on the doubly
stochastic rescaling of the matrix, the twin of jline.lib.perm.HuberLawSampler.

Both samplers take a seed and draw from a numpy Generator, so a run is
reproducible; the JAR twins use an unseeded java.util.Random, hence the two
codebases agree in distribution but not sample by sample.
"""

import time

import numpy as np

from scipy.optimize import linear_sum_assignment

from .base import PermSolver, require_full_support


class AdaPartSampler(PermSolver):
    """
    Adaptive partitioning (AdaPart) sampler for the permanent.

    Port of jline.lib.perm.AdaPartSampler. The state space of the permutations
    is recursively partitioned, each part is bounded by the Soules bound, and a
    part is drawn with probability proportional to its bound; the acceptance
    ratio of the resulting rejection sampler scales the root bound into an
    unbiased estimate of the permanent.

    Modes: 'classic' stops after maximum_accepted_samples acceptances, 'time'
    after maximum_time milliseconds, 'sample' after maximum_samples draws.
    """

    def __init__(self, matrix: np.ndarray, maximum_accepted_samples: int = 100,
                 maximum_time: float = 30000.0, maximum_samples: int = 450,
                 mode: str = "classic", solve: bool = False, seed=None,
                 maximum_draws: int = 1000000):
        """
        Initialize the AdaPart sampler.

        Args:
            matrix: Nonnegative matrix for which to approximate the permanent
            maximum_accepted_samples: Acceptance budget of the 'classic' mode
            maximum_time: Time budget in milliseconds of the 'time' mode
            maximum_samples: Draw budget of the 'sample' mode
            mode: One of 'classic', 'time', 'sample'
            solve: If True, compute immediately
            seed: Seed of the numpy Generator used for the draws
        """
        super().__init__(matrix)
        self.maximum_accepted_samples = maximum_accepted_samples
        self.maximum_time = maximum_time
        self.maximum_samples = maximum_samples
        self.maximum_draws = maximum_draws
        self.mode = mode
        self.random = np.random.default_rng(seed)
        if self.n > 0:
            require_full_support(self.matrix, 'adapart')

        self.sample_accepted = []
        self.sample_time = []
        self.perm_step = []

        if solve:
            self.solve()

    def compute(self):
        """Run the sampler in the configured mode."""
        if self.mode == "time":
            self.value = self._sampling_time()
        elif self.mode == "sample":
            self.value = self._sampling_sample()
        else:
            self.value = self._sampling_classic()

    def _sampling_classic(self) -> float:
        """Sample until maximum_accepted_samples draws have been accepted."""
        z_ub = self._soules_bound(self.matrix)
        accepted = 0
        total = 0
        start_time = time.perf_counter()
        accepted_list = []
        time_list = []
        # Bounded independently of the scaling: with perm(A) = 0 the acceptance
        # probability is 0 and this loop would never terminate. A cap that
        # RETURNS a number would be a workaround, so it raises.
        while accepted < self.maximum_accepted_samples:
            if total >= self.maximum_draws:
                raise ValueError(
                    "Only " + str(accepted) + " of the " +
                    str(self.maximum_accepted_samples) + " required acceptances "
                    "were obtained in " + str(total) + " draws. Raise "
                    "maximum_draws or use the exact engine.")
            sample = self._sample()
            accepted += sample
            total += 1
            accepted_list.append(sample)
            time_list.append((time.perf_counter() - start_time) * 1000.0)
        self.sample_accepted = accepted_list
        self.sample_time = time_list
        self._compute_perm_step(z_ub)
        return z_ub * accepted / float(total)

    def _sampling_time(self) -> float:
        """Sample until the time budget is exhausted."""
        z_ub = self._soules_bound(self.matrix)
        accepted = 0
        total = 0
        start_time = time.perf_counter()
        accepted_list = []
        time_list = []
        while (time.perf_counter() - start_time) * 1000.0 < self.maximum_time:
            sample = self._sample()
            accepted += sample
            total += 1
            accepted_list.append(sample)
            time_list.append((time.perf_counter() - start_time) * 1000.0)
        self.sample_accepted = accepted_list
        self.sample_time = time_list
        self._compute_perm_step(z_ub)
        return z_ub * accepted / float(total) if total > 0 else 0.0

    def _sampling_sample(self) -> float:
        """Sample until maximum_samples draws have been taken."""
        z_ub = self._soules_bound(self.matrix)
        accepted = 0
        total = 0
        start_time = time.perf_counter()
        accepted_list = []
        time_list = []
        while len(accepted_list) < self.maximum_samples:
            sample = self._sample()
            accepted += sample
            total += 1
            accepted_list.append(sample)
            time_list.append((time.perf_counter() - start_time) * 1000.0)
        self.sample_accepted = accepted_list
        self.sample_time = time_list
        self._compute_perm_step(z_ub)
        return z_ub * accepted / float(total)

    def _sample(self) -> int:
        """Draw one partition path, returning 1 if accepted and 0 if rejected."""
        n = self.n
        s_set = {tuple([n] * n): None}
        start_time = time.perf_counter()
        while self._any_contains_n(s_set) and self._within_time(start_time):
            s_init = next(iter(s_set))
            s_matrix = self._modify_matrix(self.matrix, s_init)
            ub = self._soules_bound(s_matrix)
            zub_s = ub
            init = True
            while (ub >= zub_s or init) and self._within_time(start_time):
                init = False
                # Only elements with a free position can be refined; expanding a
                # complete assignment yields no children and stalls the sampler.
                s_list = [t for t in s_set if n in t]
                if not s_list:
                    break
                s_sub = s_list[int(self.random.integers(len(s_list)))]
                del s_set[s_sub]
                sub_matrix = self._modify_matrix(self.matrix, s_sub)
                # Discount the bound of the element actually removed, not the
                # root bound; see _kb/03-api-layer.md.
                sub_ub = self._soules_bound(sub_matrix)
                new_ub, j = self._select_column(sub_matrix, sub_ub, ub, s_sub)
                for i in range(n):
                    s_add = list(s_sub)
                    if i not in s_sub:
                        s_add[j] = i
                        s_set[tuple(s_add)] = None
                no_progress = new_ub >= ub
                ub = new_ub
                # The Soules bound is tight on matrices with equal entries, so
                # refinement cannot improve it and "refine until improved" would
                # never exit. Stop on the first non-improving expansion instead;
                # a tight bound means the draw is accepted with probability 1.
                if no_progress:
                    break
            c = self._compute_probabilities(s_set, zub_s)
            if c == len(s_set):
                return 0
            s_set = self._subset(s_set, c)
        return 1

    def _within_time(self, start_time: float) -> bool:
        """Time guard, active in 'time' mode only, mirroring the JAR."""
        if self.mode != "time":
            return True
        return (time.perf_counter() - start_time) * 1000.0 < self.maximum_time

    def _any_contains_n(self, s_set) -> bool:
        """True while some partition element still has an unassigned column."""
        for s in s_set:
            if self.n in s:
                return True
        return False

    def _select_column(self, s_matrix: np.ndarray, removed_ub: float, ub: float,
                       s_sub):
        """Pick the column whose expansion minimizes the summed Soules bound.

        Only columns still unassigned in s_sub are candidates. Scoring an
        already-assigned column just re-derives its own constraint, which always
        looks cheapest, so the sampler would re-split the same column forever
        and never complete an assignment.

        Args:
            s_matrix: Sub-problem matrix being expanded
            removed_ub: Soules bound of that sub-problem, leaving the partition
            ub: Running Soules bound summed over the whole partition
            s_sub: Partial assignment being expanded, self.n marking free columns

        Returns:
            Tuple of (updated partition bound, chosen column)
        """
        n = self.n
        ubi = np.full(n, np.inf)
        for i in range(n):
            if s_sub[i] != n:
                continue
            total = 0.0
            for j in range(n):
                assignment = [n] * n
                assignment[i] = j
                f_b = self._modify_matrix(s_matrix, assignment)
                total += self._soules_bound(f_b)
            ubi[i] = total
        j = int(np.argmin(ubi))
        return ub - removed_ub + float(ubi[j]), j

    def _compute_probabilities(self, s_set, zub_s: float) -> int:
        """Draw a partition element, or the slack index that means rejection."""
        s_list = list(s_set)
        p = [self._soules_bound(self._modify_matrix(self.matrix, s)) for s in s_list]
        if sum(p) > 0:
            p = [v / zub_s for v in p]
        slack = 1.0 - sum(p)
        p.append(slack)
        total_prob = sum(p)
        if total_prob > 0:
            p = [abs(v) / total_prob for v in p]
        rand = float(self.random.random())
        cum_sum = 0.0
        for i, v in enumerate(p):
            cum_sum += v
            if rand <= cum_sum:
                return i
        return len(p) - 1

    def _subset(self, s_set, c: int):
        """Keep the drawn element, completing it when one column is left."""
        n = self.n
        s_list = list(s_set)
        s_inter = list(s_list[c])
        if s_inter.count(n) == 1:
            existing_values = set(v for v in s_inter if v != n)
            missing_value = -1
            for v in range(n):
                if v not in existing_values:
                    missing_value = v
                    break
            s_inter[s_inter.index(n)] = missing_value
        return {tuple(s_inter): None}

    def _modify_matrix(self, m: np.ndarray, t) -> np.ndarray:
        """Zero out the entries excluded by the partial assignment t."""
        n = self.n
        mask = np.zeros((n, n), dtype=bool)
        for j in range(len(t)):
            if t[j] != n:
                mask[t[j], j] = True
        for i in range(n):
            if i not in t:
                for j in range(len(t)):
                    if t[j] == n:
                        mask[i, j] = True
        return np.where(mask, m, 0.0)

    def _soules_bound(self, m: np.ndarray) -> float:
        """Soules upper bound of the permanent, a product of column bounds."""
        n = self.n
        gamma = np.zeros(n + 1)
        factorial = 1.0
        for k in range(1, n + 1):
            factorial *= k
            gamma[k] = factorial ** (1.0 / k)
        delta = np.array([gamma[n - i] - gamma[n - i - 1] for i in range(n)])
        sorted_columns = np.sort(m, axis=0)
        m_sum = delta @ sorted_columns
        return float(np.prod(m_sum))

    def _compute_perm_step(self, z_ub: float):
        """Running estimate of the permanent after each draw."""
        steps = []
        cum_accepted = 0
        for i, accepted in enumerate(self.sample_accepted):
            cum_accepted += accepted
            steps.append(z_ub * cum_accepted / (i + 1))
        self.perm_step = steps


class HuberLawSampler(PermSolver):
    """
    Huber-Law acceptance-rejection sampler for the permanent.

    Port of jline.lib.perm.HuberLawSampler. The matrix is rescaled to be doubly
    stochastic, a permutation is drawn column by column under the Huber-Law
    upper bound on the remaining permanent, and the acceptance ratio times the
    rescaling constant estimates the permanent.

    Modes: 'classic' stops once K = 14 delta^-2 log(2/epsilon) draws have been
    accepted, 'time' after maximum_time milliseconds, 'sample' after
    number_of_samples draws.
    """

    def __init__(self, matrix: np.ndarray, delta: float = 0.1, alpha2: float = 0.000001,
                 epsilon: float = 0.1, mode: str = "classic", number_of_samples: int = 1000,
                 maximum_time: float = 30000.0, solve: bool = False, seed=None,
                 maximum_draws: int = 1000000, maximum_sinkhorn: int = 10000):
        """
        Initialize the Huber-Law sampler.

        Args:
            matrix: Nonnegative matrix for which to approximate the permanent
            delta: Relative accuracy target, sets the acceptance budget K
            alpha2: Convergence threshold of the doubly stochastic rescaling
            epsilon: Failure probability target, sets the acceptance budget K
            mode: One of 'classic', 'time', 'sample'
            number_of_samples: Draw budget of the 'sample' mode
            maximum_time: Time budget in milliseconds of the 'time' mode
            solve: If True, compute immediately
            seed: Seed of the numpy Generator used for the draws
        """
        super().__init__(matrix)
        self.delta = delta
        self.alpha2 = alpha2
        self.epsilon = epsilon
        self.mode = mode
        self.number_of_samples = number_of_samples
        self.maximum_time = maximum_time
        self.maximum_draws = maximum_draws
        self.maximum_sinkhorn = maximum_sinkhorn
        self.random = np.random.default_rng(seed)
        if self.n > 0:
            require_full_support(self.matrix, 'huberlaw')

        self.c_matrix = None
        self.rescaling_constant = 1.0
        self.sample_accepted = []
        self.sample_time = []
        self.perm_step = []

        if solve:
            self.solve()

    def get_sample_accepted(self):
        """Return the per-draw acceptance indicators of the last run."""
        return self.sample_accepted

    def get_sample_time(self):
        """Return the per-draw elapsed times in milliseconds of the last run."""
        return self.sample_time

    def get_perm_step(self):
        """Return the running permanent estimate after each draw."""
        return self.perm_step

    def compute(self):
        """Run the sampler in the configured mode."""
        if self.mode == "time":
            self.value = self._sampling_time()
        elif self.mode == "sample":
            self.value = self._sampling_sample()
        else:
            self.value = self._sampling_classic()

    def _sampling_classic(self) -> float:
        """Sample until K draws have been accepted."""
        self._rescale()
        k = int(14.0 * self.delta ** (-2) * np.log(2.0 / self.epsilon))
        start_time = time.perf_counter()
        accepted_list = []
        time_list = []
        accepted_count = 0
        # Bounded independently of the scaling, as in AdaPartSampler.
        while accepted_count < k:
            if len(accepted_list) >= self.maximum_draws:
                raise ValueError(
                    "Only " + str(accepted_count) + " of the " + str(k) +
                    " required acceptances were obtained in " +
                    str(len(accepted_list)) + " draws. The acceptance probability "
                    "is too low for this budget; raise maximum_draws, relax "
                    "delta, or use the exact engine.")
            sigma = self._sample()
            is_accepted = 0 if np.any(sigma == self.n) else 1
            accepted_list.append(is_accepted)
            time_list.append((time.perf_counter() - start_time) * 1000.0)
            if is_accepted == 1:
                accepted_count += 1
        self.sample_accepted = accepted_list
        self.sample_time = time_list
        self._compute_perm_step()
        return sum(self.sample_accepted) / float(len(self.sample_accepted)) * self.rescaling_constant

    def _sampling_time(self) -> float:
        """Sample until the time budget is exhausted."""
        self._rescale()
        start_time = time.perf_counter()
        accepted_list = []
        time_list = []
        while (time.perf_counter() - start_time) * 1000.0 < self.maximum_time:
            sigma = self._sample()
            is_accepted = 0 if np.any(sigma == self.n) else 1
            accepted_list.append(is_accepted)
            time_list.append((time.perf_counter() - start_time) * 1000.0)
        self.sample_accepted = accepted_list
        self.sample_time = time_list
        self._compute_perm_step()
        if not self.sample_accepted:
            return 0.0
        return sum(self.sample_accepted) / float(len(self.sample_accepted)) * self.rescaling_constant

    def _sampling_sample(self) -> float:
        """Sample until number_of_samples draws have been taken."""
        self._rescale()
        start_time = time.perf_counter()
        accepted_list = []
        time_list = []
        while len(accepted_list) < self.number_of_samples:
            sigma = self._sample()
            is_accepted = 0 if np.any(sigma == self.n) else 1
            accepted_list.append(is_accepted)
            time_list.append((time.perf_counter() - start_time) * 1000.0)
        self.sample_accepted = accepted_list
        self.sample_time = time_list
        self._compute_perm_step()
        return sum(self.sample_accepted) / float(len(self.sample_accepted)) * self.rescaling_constant

    def _sample(self) -> np.ndarray:
        """Draw one permutation, returning a vector of n entries when rejected."""
        n = self.n
        m = self.c_matrix.copy()
        sigma = np.zeros(n, dtype=int)
        for j in range(n):
            row_sums = m.sum(axis=1)
            ub = float(np.prod([self._h(r) for r in row_sums])) / np.exp(n)
            p = self._precomputing(m, j)
            normalized_p = p / ub
            prob = np.zeros(n + 1)
            prob[:n] = normalized_p
            prob[n] = 1.0 - float(normalized_p.sum())
            if prob[n] < 0:
                pos_sum = float(prob[:n].sum())
                prob[:n] = prob[:n] / pos_sum
                prob[n] = 0.0
            rand_val = float(self.random.random())
            cum_sum = 0.0
            selected_i = n
            for i in range(n + 1):
                cum_sum += prob[i]
                if rand_val <= cum_sum:
                    selected_i = i
                    break
            if selected_i == n:
                return np.full(n, n, dtype=int)
            sigma[j] = selected_i
            new_matrix = np.zeros((n, n))
            keep = np.ones((n, n), dtype=bool)
            keep[selected_i, :] = False
            keep[:, j] = False
            new_matrix[keep] = m[keep]
            new_matrix[selected_i, j] = m[selected_i, j]
            m = new_matrix
        return sigma

    @staticmethod
    def _h(r: float) -> float:
        """Huber-Law bound factor of a row with remaining mass r."""
        if r >= 1.0:
            return r + 0.5 * np.log(max(r, 1.0)) + np.e - 1
        return 1 + (np.e - 1) * r

    def _precomputing(self, m: np.ndarray, j: int) -> np.ndarray:
        """Unnormalized selection weights of each row for column j."""
        n = self.n
        c = m[:, j].copy()
        r = m.sum(axis=1) - c
        hr = np.array([self._h(v) for v in r])
        hr_product = float(np.prod(hr))
        exp_factor = np.exp(n - 1.0)
        return hr_product / hr * c / exp_factor

    def _rescale(self):
        """Rescale the matrix to doubly stochastic and set the scaling constant."""
        n = self.n
        # The matrix is strictly positive here (require_full_support in the
        # constructor), so no log floor is needed, and the assignment is the
        # real maximum-weight one rather than a row-by-row greedy.
        log_matrix = np.log(self.matrix)
        assignment = self._hungarian_assignment(log_matrix)

        max_element = float(self.matrix.max())
        m_scaled = self.matrix / max_element

        # alpha3 is a permanent lower bound of the scaled matrix, the one floored below
        alpha3 = 1.0
        for i in range(n):
            alpha3 *= m_scaled[i, assignment[i]]

        alpha1 = alpha3 * self.delta / 3 / _factorial(n)
        m_scaled = np.maximum(m_scaled, alpha1)

        doubly_stochastic, x, y = self._make_doubly_stochastic(m_scaled)

        z = np.eye(n)
        for i in range(n):
            z[i, i] = 1.0 / float(doubly_stochastic[i, :].max())
        self.c_matrix = z @ doubly_stochastic

        row_sums = self.c_matrix.sum(axis=1)
        h_product = 1.0
        for i in range(n):
            h_product *= self._h(row_sums[i]) / np.e

        diagonal_product = 1.0
        for i in range(n):
            diagonal_product *= x[i, i] * y[i, i] * z[i, i]
        self.rescaling_constant = h_product / diagonal_product * max_element ** n

    def _hungarian_assignment(self, cost_matrix: np.ndarray) -> np.ndarray:
        """
        Maximum-weight perfect assignment of a square cost matrix.

        This used to be a row-by-row greedy despite the name, and the greedy
        returns a zero-weight assignment on inputs that admit a positive one:
        on [[1, 2], [0, 3]] row 0 takes the larger entry in column 1, leaving
        row 1 with the zero in column 0. The weight is alpha3 in _rescale, and
        alpha3 sets the flooring level alpha1 = alpha3*delta/(3*n!) of the
        Huber-Law sampler, so a suboptimal assignment silently weakens the
        method's own guarantee. A correct assignment makes alpha3 > 0 whenever
        the matrix has a positive permanent.

        The MATLAB and JAR twins solve the same problem and agree on the
        optimal VALUE, which is all alpha3 depends on; they need not return the
        same permutation when the optimum is degenerate.
        """
        rows, cols = linear_sum_assignment(cost_matrix, maximize=True)
        assignment = np.full(self.n, -1, dtype=int)
        assignment[rows] = cols
        return assignment

    def _make_doubly_stochastic(self, m: np.ndarray):
        """Alternate column and row normalization, returning (B, X, Y)."""
        n = self.n
        result = m.copy()
        x = np.eye(n)
        y = np.eye(n)
        max_row_error = np.inf
        max_col_error = np.inf
        sweeps = 0
        # Capped: a row that sums to zero leaves max_row_error at 1 forever and
        # the guarded normalization below skips it, so this loop used to spin
        # without terminating. A cap that RETURNS is a workaround; this raises.
        while max_row_error > self.alpha2 or max_col_error > self.alpha2:
            sweeps += 1
            if sweeps > self.maximum_sinkhorn:
                raise ValueError(
                    "The doubly stochastic rescaling did not converge in " +
                    str(self.maximum_sinkhorn) + " sweeps (row error " +
                    str(max_row_error) + ", column error " + str(max_col_error) +
                    " against a tolerance of " + str(self.alpha2) + "). The "
                    "usual cause is a matrix without total support.")
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
        return result, x, y

    def _compute_perm_step(self):
        """Running estimate of the permanent after each draw."""
        steps = []
        cum_accepted = 0
        for i, accepted in enumerate(self.sample_accepted):
            cum_accepted += accepted
            steps.append(self.rescaling_constant * cum_accepted / (i + 1))
        self.perm_step = steps


def _factorial(n: int) -> float:
    """Compute n! in double precision by the plain product, as in the JAR."""
    result = 1.0
    for i in range(1, n + 1):
        result *= i
    return result
