"""
Deterministic approximations of the matrix permanent.

BethePermanent is the sum-product (belief propagation) approximation of the
permanent, the twin of jline.lib.perm.BethePermanent. HeuristicPermanent is the
Sinkhorn scaling plus mean-field / van der Waerden and Gurvits capacity bound
heuristic, the twin of MATLAB perm_heur.m and of jline.lib.perm.HeuristicPermanent.
SaddlePointPermanent is the Laplace expansion of the coefficient integral, the
homogeneous variant of cache_spm and the twin of MATLAB perm_spm.m.
"""

import math

import numpy as np

from .base import PermSolver, require_full_support

# MIN_VALUE guards the LOGARITHMS of message products in _bethe against
# underflow. It is deliberately NOT applied to the input matrix: flooring the
# input is what fabricates a permanent of n!*eps where the truth is zero, so a
# non-positive entry is refused outright instead. See require_full_support.
MIN_VALUE = 2.220446049250314e-16


class BethePermanent(PermSolver):
    """
    Sum Product Algorithm (SPA) approximation of the Bethe permanent.

    Port of jline.lib.perm.BethePermanent. Two families of messages, the right
    going messages r and the left going messages l, are iterated to a fixed
    point on the square root of the matrix, and the Bethe free energy of that
    fixed point is exponentiated to give the estimate. The Bethe permanent is a
    lower bound of the permanent for nonnegative matrices.
    """

    def __init__(self, matrix: np.ndarray, epsilon: float = 0.001,
                 max_iteration: int = 200000, solve: bool = False):
        """
        Initialize the Bethe permanent solver.

        Args:
            matrix: Matrix for which to approximate the permanent
            epsilon: Squared message change below which the iteration stops
            max_iteration: Maximum number of message passing iterations
            solve: If True, compute immediately
        """
        super().__init__(matrix)
        self.epsilon = epsilon
        self.max_iteration = max_iteration
        if self.n > 0:
            require_full_support(self.matrix, 'bethe')
        self.matrix_sqrt = np.sqrt(self.matrix)

        if solve:
            self.solve()

    def compute(self):
        """Run the sum product algorithm and evaluate the Bethe free energy."""
        self.value = self._spa()

    def _spa(self) -> float:
        """Iterate the messages to a fixed point and return the Bethe estimate."""
        n = self.n
        r_past = np.ones((n, n))
        l_past = np.ones((n, n))

        r, l = self._update(l_past)

        iteration = 0
        while self._convergence(r_past, l_past, r, l) > self.epsilon and iteration < self.max_iteration:
            iteration += 1
            r_past = r.copy()
            l_past = l.copy()
            r, l = self._update(l_past)

        return self._bethe(l, r)

    def _update(self, l: np.ndarray):
        """Update the right going message and the left going message."""
        s = self.matrix_sqrt
        with np.errstate(divide='ignore', invalid='ignore'):
            row_terms = s * l
            denom_r = row_terms.sum(axis=1) - np.diagonal(row_terms)
            r1 = s / denom_r[:, np.newaxis]

            col_terms = s * r1
            denom_l = col_terms.sum(axis=0) - np.diagonal(col_terms)
            l1 = s / denom_l[np.newaxis, :]
        return r1, l1

    @staticmethod
    def _convergence(r0: np.ndarray, l0: np.ndarray, r1: np.ndarray, l1: np.ndarray) -> float:
        """Squared difference between past and present messages."""
        return float(((r0 - r1) ** 2).sum() + ((l0 - l1) ** 2).sum())

    def _bethe(self, l: np.ndarray, r: np.ndarray) -> float:
        """Compute the Bethe permanent from the left and right going messages."""
        s = self.matrix_sqrt
        term1 = np.maximum((s * l).sum(axis=1), MIN_VALUE)
        term2 = np.maximum((s * r).sum(axis=0), MIN_VALUE)
        term3 = np.maximum(r * l + 1.0, MIN_VALUE)

        with np.errstate(divide='ignore', invalid='ignore'):
            log_value = np.log(term1).sum() + np.log(term2).sum() - np.log(term3).sum()
            result = float(np.exp(log_value))
        if not np.isfinite(result):
            return 0.0
        return result


class HeuristicPermanent(PermSolver):
    """
    Heuristic approximation to the permanent of a nonnegative matrix.

    Port of MATLAB perm_heur.m and of jline.lib.perm.HeuristicPermanent. The
    algorithm Sinkhorn scales the matrix to be approximately doubly stochastic,
    averages a mean-field van der Waerden estimate with a Gurvits capacity
    bound on the scaled matrix, and undoes the scaling.
    """

    def __init__(self, matrix: np.ndarray, tolerance: float = 1e-10,
                 max_iterations: int = 1000, solve: bool = False):
        """
        Initialize the heuristic permanent solver.

        Args:
            matrix: Nonnegative matrix for which to approximate the permanent
            tolerance: Convergence threshold of the Sinkhorn scaling
            max_iterations: Maximum number of Sinkhorn iterations
            solve: If True, compute immediately

        Raises:
            ValueError: If the matrix has a negative element
        """
        super().__init__(matrix)
        self.tolerance = tolerance
        self.max_iterations = max_iterations

        negatives = np.argwhere(self.matrix < 0.0)
        if negatives.size > 0:
            i, j = negatives[0]
            raise ValueError("Matrix must be non-negative. Found negative element at ("
                             + str(i) + ", " + str(j) + "): " + str(self.matrix[i, j]))
        if self.n > 0:
            require_full_support(self.matrix, 'heur')

        if solve:
            self.solve()

    def compute(self):
        """Compute the heuristic permanent approximation."""
        self.value = self._compute_heuristic_permanent()

    def _compute_heuristic_permanent(self) -> float:
        """Sinkhorn scaling followed by the mean-field and capacity estimates."""
        n = self.n
        # A zero used to be replaced by 1e-15 here. That is not invertible: it
        # changes the permanent by n!*eps, which is O(1) by n=18. The
        # constructor refuses such a matrix instead.
        working_matrix = self.matrix.copy()

        b, r, c = self._sinkhorn_scaling(working_matrix)

        row_sums = b.sum(axis=1)
        row_prod = float(np.prod(row_sums))
        p_meanfield = _factorial(n) * (row_prod / float(n) ** n)

        cap = np.exp(np.log(row_sums).sum() / n)
        p_gurvits = _factorial(n) * (cap / n) ** n

        p_est = 0.5 * (p_meanfield + p_gurvits)

        scale_factor = float(np.prod(1.0 / r) * np.prod(1.0 / c))
        return float(p_est * scale_factor)

    def _sinkhorn_scaling(self, input_matrix: np.ndarray):
        """Scale a matrix towards double stochasticity, returning (B, r, c)."""
        n = self.n
        b = input_matrix
        r = np.ones(n)
        c = np.ones(n)

        converged = False
        for _ in range(self.max_iterations):
            r = 1.0 / (b @ c)
            c = 1.0 / (b.T @ r)
            if np.max(np.abs(r * (b @ c) - 1.0)) < self.tolerance:
                converged = True
                break
        if not converged:
            raise ValueError(
                "The Sinkhorn scaling did not converge to a doubly stochastic "
                "matrix in " + str(self.max_iterations) + " sweeps (margin error "
                + str(float(np.max(np.abs(r * (b @ c) - 1.0)))) + " against a "
                "tolerance of " + str(self.tolerance) + "). The estimate below "
                "assumes convergence, so no value is returned. The usual cause "
                "is a matrix without total support.")

        scaled_b = r[:, np.newaxis] * b * c[np.newaxis, :]
        return scaled_b, r, c


class SaddlePointPermanent(PermSolver):
    """
    Saddle-point (SPM) approximation of the permanent of a positive matrix.

    Port of MATLAB perm_spm.m and of jline.lib.perm.SaddlePointPermanent. This
    is the HOMOGENEOUS variant of cache_spm: both evaluate the same Cauchy
    integral by Laplace's method and differ only in the generating function
    whose coefficient they extract,

        cache_spm  E(m) = prod_l m_l! [prod_l z_l^m_l] prod_k (1 + sum_l g_kl z_l)
        perm_spm   P    = prod_l m_l! [prod_l z_l^m_l] prod_k (    sum_l A_kl z_l)

    The cache factor carries a "+1" because an item may stay out of the cache,
    so the coefficient it extracts is a rectangular permanent over n items and
    sum(m) < n slots. Dropping the "+1" forces every row to be matched, which
    is exactly the permanent and requires sum(m) == n. That is the one case
    cache_spm cannot serve: at n == sum(m) its multipliers diverge and it falls
    back on cache_erec. Here the integrand is homogeneous, and the saddle point
    is interior in the h-1 directions that survive.

    Method. With z_l = xi_l exp(i th_l) the saddle point in xi solves

        sum_k A_kl xi_l / (sum_j A_kj xi_j) = m_l,   l = 1..h,

    that is, P_kl = A_kl xi_l / S_k with S = A @ xi is the diagonal scaling of A
    to row sums 1 and column sums m (Sinkhorn scaling; doubly stochastic when m
    is all ones). There phi = sum_k log S_k - sum_l m_l log xi_l is the log of
    the Gurvits capacity, an upper bound on the log permanent. The Gaussian
    correction uses H = diag(m) - P.T @ P, a weighted graph Laplacian on the
    columns: H @ ones = 0, which is the invariance of the integrand under
    th -> th + c*ones that homogeneity creates. That direction is a full period
    rather than a Gaussian, so it contributes 2*pi and leaves an (h-1)
    dimensional Laplace integral. Any principal (h-1) submatrix serves, since
    all cofactors of a Laplacian are equal, and

        log P = sum_l log(m_l!) - (h-1)/2 log(2 pi) + phi - 1/2 log det(H_red).

    Accuracy. Exact for h == 1, where the permanent is n! prod_k A[k,0]. It is a
    genuine asymptotic expansion as min(m) grows with h fixed, the ratio to the
    exact permanent falling from 1.11 at m = (2,2,2) to 1.02 at m = (3,3). At
    m = ones the dimension of the integral grows with the expansion parameter
    and the leading term keeps a systematic bias: on the n x n matrix of ones it
    returns (2 pi)^(-(n-1)/2) n^(n+1/2) against the exact n!, a ratio tending to
    (e/sqrt(2 pi))^n = 1.084^n, and random positive matrices track that closely
    (1.31 at n = 4, 1.87 at n = 8). So at m = ones it OVERESTIMATES, with a
    spread across matrices far tighter than the bias itself, and it is not a
    bound in either direction. The raw capacity is off by 352x on the same n = 8
    instances, and BethePermanent is a genuine lower bound.

    Attributes:
        value: the approximate permanent
        log_value: its logarithm, correct even when value overflows
        xi: the saddle point, unit geometric mean, zero on a dropped column
        log_capacity: log Gurvits capacity, an upper bound on the log permanent
    """

    def __init__(self, matrix: np.ndarray, m=None, tolerance: float = 1e-11,
                 max_iterations: int = 10000, solve: bool = False):
        """
        Initialize the saddle-point permanent solver.

        Args:
            matrix: n x h strictly positive matrix
            m: column multiplicities, non-negative integers summing to n
               (default: ones, which requires a square matrix)
            tolerance: margin on the column sums at which the scaling stops
            max_iterations: maximum number of scaling sweeps
            solve: If True, compute immediately

        Raises:
            ValueError: if the matrix is negative, is not strictly positive, or
                the multiplicities do not sum to the number of rows
        """
        super().__init__(matrix)
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.log_value = 0.0
        self.log_capacity = 0.0

        if self.matrix.size == 0:
            self.h = 0
            self.m = np.zeros(0)
            self.xi = np.zeros(0)
            self.value = 1.0    # the permanent of the empty matrix is 1
            return

        if self.matrix.ndim != 2:
            raise ValueError("perm_spm requires a two-dimensional matrix.")
        n, h = self.matrix.shape
        self.h = h
        negatives = np.argwhere(self.matrix < 0.0)
        if negatives.size > 0:
            i, j = negatives[0]
            raise ValueError("Matrix must be non-negative. Found negative element at ("
                             + str(i) + ", " + str(j) + "): " + str(self.matrix[i, j]))
        if m is None:
            if h != n:
                raise ValueError("Without column multiplicities the matrix must be square; "
                                 "it is " + str(n) + "x" + str(h) + ".")
            m = np.ones(n)
        m = np.asarray(m, dtype=float).ravel()
        if m.size != h:
            raise ValueError("The multiplicity vector has " + str(m.size) + " entries against "
                             + str(h) + " columns.")
        if np.any(m < 0) or np.any(np.abs(m - np.round(m)) > 0):
            raise ValueError("Column multiplicities must be non-negative integers.")
        if int(round(float(np.sum(m)))) != n:
            raise ValueError("The column multiplicities must sum to the number of rows (sum(m) = "
                             + str(float(np.sum(m))) + " against " + str(n) + " rows). The "
                             "integrand is homogeneous of degree " + str(n) + ", so every other "
                             "coefficient of it is exactly zero.")
        require_full_support(self.matrix, 'spm')
        self.m = m
        self.xi = np.zeros(h)

        if solve:
            self.solve()

    def compute(self):
        """Locate the saddle point and evaluate the Laplace expansion there."""
        if self.matrix.size == 0:
            self.value = 1.0
            self.log_value = 0.0
            return
        self.log_value, self.log_capacity = self._expand()
        self.value = float(np.exp(self.log_value))

    def _expand(self):
        """Return (log estimate, log capacity) at the saddle point."""
        # A column repeated zero times leaves the permanent unchanged, and its
        # xi is a boundary of the Laplace integral rather than a direction of
        # it, so it must leave the expansion. Dropping it is exact: setting
        # z_l = 0 removes the column, and prod_l m_l! is unchanged since 0! = 1.
        keep = np.flatnonzero(self.m > 0)
        a = self.matrix[:, keep]
        mk = self.m[keep]
        hk = keep.size

        xik = self._scale(a, mk)
        self.xi = np.zeros(self.h)
        self.xi[keep] = xik

        s = a @ xik
        p = (a * xik[np.newaxis, :]) / s[:, np.newaxis]
        log_capacity = float(np.sum(np.log(s)) - np.dot(mk, np.log(xik)))

        # H = diag(mk) - P'P is a Laplacian, so it is singular along ones and
        # all of its principal cofactors are equal; the last index is dropped
        # only because one has to be. Strict positivity of A makes the column
        # graph complete, hence H_red positive definite and Cholesky the right
        # factor.
        if hk > 1:
            h_full = np.diag(mk) - p.T @ p
            try:
                chol = np.linalg.cholesky(h_full[:hk - 1, :hk - 1])
            except np.linalg.LinAlgError:
                raise ValueError("The reduced Hessian is not positive definite, so the saddle "
                                 "point is degenerate and the Gaussian factor does not exist.")
            log_det = 2.0 * float(np.sum(np.log(np.diag(chol))))
        else:
            # no direction survives the homogeneity, and det of the empty matrix is 1
            log_det = 0.0

        log_fact = float(sum(math.lgamma(v + 1.0) for v in mk))
        log_value = (log_fact - 0.5 * (hk - 1) * math.log(2.0 * math.pi)
                     + log_capacity - 0.5 * log_det)
        return log_value, log_capacity

    def _scale(self, a: np.ndarray, mk: np.ndarray) -> np.ndarray:
        """Scale a to row sums 1 and column sums mk, returning the multipliers."""
        xik = np.ones(mk.size)
        margin = np.inf
        for _ in range(self.max_iterations):
            s = a @ xik
            colsum = xik * (a.T @ (1.0 / s))
            margin = float(np.max(np.abs(colsum - mk)))
            if margin < self.tolerance:
                return xik
            xik = xik * (mk / colsum)
            xik = xik / np.exp(np.mean(np.log(xik)))    # the saddle is a ray; pin its scale
        raise ValueError(
            "The scaling to row sums 1 and column sums m did not converge in "
            + str(self.max_iterations) + " sweeps (margin error " + str(margin)
            + " against a tolerance of " + str(self.tolerance) + "). The expansion assumes the "
            "saddle point, so no value is returned. The usual cause is a matrix without total "
            "support.")


def perm_heur(matrix: np.ndarray) -> float:
    """
    Heuristic approximation to the permanent of a nonnegative matrix.

    Twin of MATLAB perm_heur.m, using the default Sinkhorn tolerance 1e-10 and
    1000 iterations.

    Args:
        matrix: Nonnegative matrix

    Returns:
        Approximate permanent value
    """
    solver = HeuristicPermanent(matrix, solve=True)
    return solver.value


def perm_bethe(matrix: np.ndarray, epsilon: float = 0.001, max_iteration: int = 200000) -> float:
    """
    Bethe (sum product) approximation of the permanent.

    Args:
        matrix: Nonnegative matrix
        epsilon: Squared message change below which the iteration stops
        max_iteration: Maximum number of message passing iterations

    Returns:
        Approximate permanent value
    """
    solver = BethePermanent(matrix, epsilon, max_iteration, solve=True)
    return solver.value


def perm_spm(matrix: np.ndarray, m=None, tolerance: float = 1e-11,
             max_iterations: int = 10000) -> float:
    """
    Saddle-point (SPM) approximation of the permanent of a positive matrix.

    Twin of MATLAB perm_spm.m. Approximates the permanent of the matrix built
    from the n x h matrix by repeating column l exactly m[l] times, so sum(m)
    must equal n; with m omitted the matrix must be square and the permanent of
    the matrix itself is approximated. See SaddlePointPermanent for the method
    and for what the estimate does and does not guarantee.

    Args:
        matrix: n x h strictly positive matrix
        m: column multiplicities summing to n (default: ones, requires square)
        tolerance: margin on the column sums at which the scaling stops
        max_iterations: maximum number of scaling sweeps

    Returns:
        Approximate permanent value
    """
    solver = SaddlePointPermanent(matrix, m, tolerance, max_iterations, solve=True)
    return solver.value


def _factorial(n: int) -> float:
    """
    Compute n! in double precision.

    Exact for n <= 170, the largest factorial representable as a double;
    beyond that Stirling's approximation avoids the overflow.
    """
    if n <= 1:
        return 1.0
    if n <= 170:
        result = 1.0
        for i in range(2, n + 1):
            result *= i
        return result
    return float(np.sqrt(2.0 * np.pi * n) * (n / np.e) ** n)
